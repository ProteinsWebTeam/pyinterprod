import os
from typing import Tuple

import cx_Oracle

from .. import logger, orautils


def _condense(matches: dict):
    for entry_acc in matches:
        fragments = []
        start = end = None
        for s, e in sorted(matches[entry_acc]):
            if start is None:
                # Leftmost fragment
                start = s
                end = e
            elif s > end:
                """
                      end
                    ----] [----
                          s
                -> new fragment
                """
                fragments.append((start, end))
                start = s
                end = e
            elif e > end:
                """
                        end
                    ----]
                      ------]
                            e
                -> extend
                """
                end = e

        fragments.append((start, end))
        matches[entry_acc] = fragments


def build_xref_condensed(user: str, dsn: str):
    logger.info("building XREF_CONDENSED")
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()

    orautils.drop_table(cur, "INTERPRO", "XREF_CONDENSED")
    cur.execute(
        """
        CREATE TABLE INTERPRO.XREF_CONDENSED
        (
            PROTEIN_AC VARCHAR2(15) NOT NULL,
            ENTRY_AC VARCHAR2(9) NOT NULL,
            ENTRY_TYPE CHAR(1) NOT NULL,
            ENTRY_NAME VARCHAR2(100) NOT NULL,
            POS_FROM NUMBER(5) NOT NULL,
            POS_TO NUMBER(5) NOT NULL
        )
        PARTITION BY LIST (ENTRY_TYPE) (
          PARTITION PART_A VALUES ('A'),
          PARTITION PART_B VALUES ('B'),
          PARTITION PART_C VALUES ('C'),
          PARTITION PART_D VALUES ('D'),
          PARTITION PART_F VALUES ('F'),
          PARTITION PART_H VALUES ('H'),
          PARTITION PART_P VALUES ('P'),
          PARTITION PART_R VALUES ('R'),
          PARTITION PART_U VALUES ('U')
        ) COMPRESS NOLOGGING
        """
    )

    cur.execute(
        """
        SELECT E.ENTRY_AC, E.NAME, E.ENTRY_TYPE, EM.METHOD_AC
        FROM INTERPRO.ENTRY E
        INNER JOIN INTERPRO.ENTRY2METHOD EM
          ON E.ENTRY_AC = EM.ENTRY_AC
        WHERE E.CHECKED = 'Y'
        """
    )

    entries = {}
    signatures = {}
    for entry_acc, name, entry_type, method_acc in cur:
        signatures[method_acc] = entry_acc
        entries[entry_acc] = (entry_type, name)

    cur.execute(
        """
        SELECT PROTEIN_AC, METHOD_AC, POS_FROM, POS_TO, FRAGMENTS
        FROM INTERPRO.MATCH
        ORDER BY PROTEIN_AC
        """
    )

    matches = {}
    _protein_acc = None
    query = """
      INSERT /*+ APPEND */ INTO INTERPRO.XREF_CONDENSED
      VALUES (:1, :2, :3, :4, :5, :6)
    """
    table = orautils.TablePopulator(con, query, autocommit=True)

    for protein_acc, method_acc, pos_from, pos_to, fragments in cur:
        if protein_acc != _protein_acc:
            if matches:
                _condense(matches)
                for entry_acc, frags in matches.items():
                    entry_type, name = entries[entry_acc]
                    for start, end in frags:
                        table.insert((_protein_acc, entry_acc, entry_type,
                                      name, start, end))

            _protein_acc = protein_acc
            matches = {}

        try:
            entry_acc = signatures[method_acc]
        except KeyError:
            # Not integrated or integrated in an unchecked entry
            continue

        if entry_acc in matches:
            entry = matches[entry_acc]
        else:
            entry = matches[entry_acc] = []

        if fragments is not None:
            for frag in fragments.split(','):
                start, end, _ = frag.split('-')
                start = int(start)
                end = int(end)

                if start <= end:
                    entry.append((start, end))

        if not entry:
            entry.append((pos_from, pos_to))

    if matches:
        _condense(matches)
        for entry_acc, frags in matches.items():
            entry_type, name = entries[entry_acc]
            for start, end in frags:
                table.insert((_protein_acc, entry_acc, entry_type,
                              name, start, end))

    table.close()

    logger.info("gathering statistics")
    orautils.gather_stats(cur, "INTERPRO", "XREF_CONDENSED")

    logger.info("creating indexes")
    for col in ("PROTEIN_AC", "ENTRY_AC"):
        cur.execute(
            """
            CREATE INDEX I_XREF_CONDENSED${0}
            ON INTERPRO.XREF_CONDENSED ({0}) NOLOGGING
            """.format(col)
        )

    orautils.grant(cur, "INTERPRO", "XREF_CONDENSED", "SELECT", "KRAKEN")
    cur.close()
    con.close()


def build_xref_summary(user: str, dsn: str):
    logger.info("building XREF_SUMMARY")
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    orautils.drop_table(cur, "INTERPRO", "XREF_SUMMARY")
    cur.execute(
        """
        CREATE TABLE INTERPRO.XREF_SUMMARY
        (
            DBCODE CHAR(1) NOT NULL,
            PROTEIN_AC VARCHAR2(15) NOT NULL,
            ENTRY_AC VARCHAR2(9) NOT NULL,
            SHORT_NAME VARCHAR2(30),
            METHOD_AC VARCHAR2(25) NOT NULL,
            METHOD_NAME VARCHAR2(100),
            POS_FROM NUMBER(5) NOT NULL,
            POS_TO NUMBER(5) NOT NULL,
            MATCH_STATUS CHAR(1) NOT NULL,
            SCORE FLOAT(*),
            FRAGMENTS VARCHAR2(400)
        )
        PARTITION BY LIST (DBCODE) (
          PARTITION PART_B VALUES ('B'),
          PARTITION PART_D VALUES ('D'),
          PARTITION PART_F VALUES ('F'),
          PARTITION PART_H VALUES ('H'),
          PARTITION PART_J VALUES ('J'),
          PARTITION PART_M VALUES ('M'),
          PARTITION PART_N VALUES ('N'),
          PARTITION PART_P VALUES ('P'),
          PARTITION PART_Q VALUES ('Q'),
          PARTITION PART_R VALUES ('R'),
          PARTITION PART_U VALUES ('U'),
          PARTITION PART_V VALUES ('V'),
          PARTITION PART_X VALUES ('X'),
          PARTITION PART_Y VALUES ('Y')
        ) COMPRESS NOLOGGING
        """
    )

    cur.execute(
        """
        INSERT /*+ APPEND */ INTO INTERPRO.XREF_SUMMARY
        SELECT
          MA.DBCODE,
          MA.PROTEIN_AC,
          E.ENTRY_AC,
          E.SHORT_NAME,
          MA.METHOD_AC,
          ME.NAME,
          MA.POS_FROM,
          MA.POS_TO,
          MA.STATUS,
          MA.SCORE,
          MA.FRAGMENTS
        FROM INTERPRO.MATCH MA
        INNER JOIN INTERPRO.METHOD ME ON MA.METHOD_AC = ME.METHOD_AC
        INNER JOIN INTERPRO.ENTRY2METHOD EM ON ME.METHOD_AC = EM.METHOD_AC
        INNER JOIN INTERPRO.ENTRY E ON EM.ENTRY_AC = E.ENTRY_AC
        WHERE EM.EVIDENCE != 'EXC'
        AND E.CHECKED = 'Y'
        """
    )
    con.commit()

    logger.info("gathering statistics")
    orautils.gather_stats(cur, "INTERPRO", "XREF_SUMMARY")

    logger.info("creating indexes")
    for col in ("PROTEIN_AC", "ENTRY_AC", "METHOD_AC"):
        cur.execute(
            """
            CREATE INDEX I_XREF_SUMMARY${0}
            ON INTERPRO.XREF_SUMMARY ({0}) NOLOGGING
            """.format(col)
        )

    orautils.grant(cur, "INTERPRO", "XREF_SUMMARY", "SELECT", "KRAKEN")
    cur.close()
    con.close()


def export_databases(user: str, dsn: str, dst: str):
    logger.info("exporting dat/tab files")
    os.makedirs(dst, exist_ok=True)
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    cur.execute(
        """
        SELECT
          PROTEIN_AC,
          METHOD_AC,
          MIN(METHOD_NAME),
          MIN(DBCODE),
          MIN(ENTRY_AC),
          MIN(SHORT_NAME),
          COUNT(*)
        FROM INTERPRO.XREF_SUMMARY
        GROUP BY PROTEIN_AC, METHOD_AC
        ORDER BY PROTEIN_AC
        """
    )

    dbcodes = {
        'B': ("SFLD", "SFLD"),
        'D': ("ProDom", "PD"),
        'F': ("PRINTS", "PP"),
        'H': ("Pfam", "PF"),
        'J': ("CDD", "CDD"),
        'M': ("PROSITE", "PR"),
        'N': ("TIGRFAMs", "TF"),
        'P': ("PROSITE", "PR"),
        'Q': ("HAMAP", "HP"),
        'R': ("SMART", "SM"),
        'U': ("PIRSF", "PI"),
        'V': ("PANTHER", "PTHR"),
        'X': ("Gene3D", "G3D"),
        'Y': ("SUPFAM", "SF"),
    }

    files = []
    handlers = {}
    for dbcode in dbcodes:
        if dbcode != 'P':
            # P (PROSITE patterns) -> same file than M (PROSITE profiles)
            dbname, dbkey = dbcodes[dbcode]

            files.append(os.path.join(dst, dbname.lower() + ".tab"))
            fh1 = open(files[-1], "wt")

            files.append(os.path.join(dst, dbname + ".dat"))
            fh2 = open(files[-1], "wt")

            handlers[dbcode] = (fh1, fh2)

    handlers['P'] = handlers['M']

    files.append(os.path.join(dst, "interpro.tab"))
    ifh1 = open(files[-1], "wt")
    files.append(os.path.join(dst, "InterPro.dat"))
    ifh2 = open(files[-1], "wt")

    entries = {}
    _protein = None
    for row in cur:
        protein_acc = row[0]
        signature_acc = row[1]
        signature_name = row[2]
        dbcode = row[3]
        entry_acc = row[4]
        entry_name = row[5]
        num_matches = int(row[6])

        dbname, dbkey = dbcodes[dbcode]
        fh1, fh2 = handlers[dbcode]

        row = (protein_acc, dbkey, signature_acc, signature_name, num_matches)
        fh1.write("{}\t{}\t{}\t{}\t{}\n".format(*row))

        if dbcode == 'X':
            # Gene3D: transform accession and do not print signature name
            fh2.write("{}   DR   {}; {}; -; {}.\n".format(*_convert_gene3d(row)))
        elif dbcode == 'F':
            # PRINTS: do not print match count
            fh2.write("{}   DR   {}; {}; {}.\n".format(*row))
        else:
            fh2.write("{}   DR   {}; {}; {}; {}.\n".format(*row))

        if protein_acc != _protein:
            for _entry in sorted(entries):
                row = (_protein, _entry, entries[_entry])
                ifh1.write("{}\t{}\t{}\n".format(*row))
                ifh2.write("{}    DR   InterPro; {}; {}.\n".format(*row))

            entries = {}
            _protein = protein_acc

        if entry_acc not in entries:
            entries[entry_acc] = entry_name

    for _entry in sorted(entries):
        row = (_protein, _entry, entries[_entry])
        ifh1.write("{}\t{}\t{}\n".format(*row))
        ifh2.write("{}    DR   InterPro; {}; {}.\n".format(*row))

    ifh1.close()
    ifh2.close()

    for fh1, fh2 in handlers.values():
        fh1.close()
        fh2.close()

    for path in files:
        os.chmod(path, 0o775)

    logger.info("complete")


def _convert_gene3d(row: Tuple[str, str, str, str, int]) -> Tuple[str, str,
                                                                  str, int]:
    # G3DSA:3.90.1580.10 -> 3.90.1580.10
    return row[0], row[1], row[2][6:], row[4]


def build_aa_iprscan(user: str, dsn: str):
    logger.info("creating AA_IPRSCAN")
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    orautils.drop_mview(cur, "IPRSCAN", "AA_IPRSCAN")
    orautils.drop_table(cur, "IPRSCAN", "AA_IPRSCAN")

    select_stmt = """
        SELECT
          UPI,
          ANALYSIS_ID AS LIBRARY_ID,
          METHOD_AC AS SIGNATURE,
          SEQ_START,
          SEQ_END,
          SEQ_FEATURE
    """

    cur.execute(
        """
        CREATE TABLE IPRSCAN.AA_IPRSCAN COMPRESS NOLOGGING
        AS {} FROM IPRSCAN.MV_IPRSCAN WHERE 1 = 0
        """.format(select_stmt)
    )

    partitions = ("TMHMM", "SIGNALP_EUK", "SIGNALP_GRAM_POSITIVE",
                  "SIGNALP_GRAM_NEGATIVE", "COILS", "PROSITE_PATTERNS",
                  "PROSITE_PROFILES", "MOBIDBLITE")
    for p in partitions:
        logger.debug("inserting partition {}".format(p))
        cur.execute(
            """
            INSERT /*+ APPEND */ INTO IPRSCAN.AA_IPRSCAN
            {}
            FROM IPRSCAN.MV_IPRSCAN PARTITION({})
            """.format(select_stmt, p)
        )
        con.commit()

    logger.debug("inserting partition PHOBIUS")
    cur.execute(
        """
        INSERT /*+ APPEND */ INTO IPRSCAN.AA_IPRSCAN
        {}
        FROM IPRSCAN.MV_IPRSCAN PARTITION(PHOBIUS)
        WHERE METHOD_AC IN ('SIGNAL_PEPTIDE','TRANSMEMBRANE')
        """.format(select_stmt)
    )
    con.commit()

    orautils.gather_stats(cur, "IPRSCAN", "AA_IPRSCAN")
    orautils.grant(cur, "IPRSCAN", "AA_IPRSCAN", "SELECT", "KRAKEN")

    logger.info("creating indices")
    for col in ("UPI", "SIGNATURE"):
        logger.debug("index on {}".format(col))
        cur.execute(
            """
            CREATE INDEX I_AA_IPRSCAN${0}
            ON IPRSCAN.AA_IPRSCAN ({0}) NOLOGGING
            """.format(col)
        )

    cur.close()
    con.close()
    logger.info("AA_IPRSCAN is ready")
