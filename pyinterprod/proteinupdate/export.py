import heapq
import pickle
import os
from tempfile import mkdtemp
from typing import Generator, Optional, Tuple

import cx_Oracle

from .. import logger, orautils


def _condense(matches: dict):
    for entry_acc in matches:
        entry_matches = []
        start = end = None
        for s, e in sorted(matches[entry_acc]):
            if start is None:
                # Leftmost match
                start = s
                end = e
            elif s > end:
                """
                      end
                    ----] [----
                          s
                -> new match
                """
                entry_matches.append((start, end))
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

        entry_matches.append((start, end))
        matches[entry_acc] = entry_matches


def _iter_file(filepath: str) -> Generator[tuple, None, None]:
    with open(filepath, "rb") as fh:
        while True:
            try:
                item = pickle.load(fh)
            except EOFError:
                break
            else:
                yield item


def _iter_matches(cur: cx_Oracle.Cursor, tmpdir: str, chunk_size: int) -> Generator[tuple, None, None]:
    # Export matches to files
    cur.execute(
        """
        SELECT M.PROTEIN_AC, EM.ENTRY_AC, M.POS_FROM, M.POS_TO
        FROM INTERPRO.MATCH M
        INNER JOIN INTERPRO.ENTRY2METHOD EM
          ON M.METHOD_AC = EM.METHOD_AC
         WHERE EM.ENTRY_AC IN (
           SELECT ENTRY_AC 
           FROM INTERPRO.ENTRY 
           WHERE CHECKED = 'Y'
        )
        """
    )

    chunk = []
    files = []
    for row in cur:
        chunk.append(row)

        if len(chunk) == chunk_size:
            chunk.sort()
            filepath = os.path.join(tmpdir, str(len(files)))
            files.append(filepath)
            with open(filepath, "wb") as fh:
                for item in chunk:
                    pickle.dump(item, fh)

            chunk = []

    if chunk:
        chunk.sort()
        filepath = os.path.join(tmpdir, str(len(files)))
        files.append(filepath)
        with open(filepath, "wb") as fh:
            for item in chunk:
                pickle.dump(item, fh)

    # Merge/sort
    for row in heapq.merge(*[_iter_file(filepath) for filepath in files]):
        yield row

    # Clean temporary files
    total_size = 0
    for filepath in files:
        total_size += os.path.getsize(filepath)
        os.remove(filepath)

    logger.info("temporary files: {:.0f} MB".format(total_size/1024/1024))


def build_xref_condensed(user: str, dsn: str, dir: Optional[str]=None,
                         chunk_size: int=1000000):
    tmpdir = mkdtemp(dir=dir)

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

    query = """
      INSERT /*+ APPEND */ INTO INTERPRO.XREF_CONDENSED
      VALUES (:1, :2, :3, :4, :5, :6)
    """
    table = orautils.TablePopulator(con, query, autocommit=True)

    cur.execute(
        """
        SELECT ENTRY_AC, ENTRY_TYPE, NAME
        FROM INTERPRO.ENTRY
        WHERE CHECKED = 'Y'
        """
    )
    entries = {row[0]: (row[1], row[2]) for row in cur}

    _protein_acc = None
    protein_matches = []
    for protein_acc, entry_acc, pos_from, pos_to in _iter_matches(cur, tmpdir, chunk_size):
        if protein_acc != _protein_acc:
            if protein_matches:
                # Condense in-place
                _condense(protein_matches)
                for _entry_acc, entry_matches in protein_matches.items():
                    entry_type, entry_name = entries[_entry_acc]

                    for _pos_from, _pos_end in entry_matches:
                        table.insert((_protein_acc, _entry_acc, entry_type,
                                      entry_name, _pos_from, _pos_end))

            protein_matches = {}
            _protein_acc = protein_acc

        if entry_acc in protein_matches:
            protein_matches[entry_acc].append((pos_from, pos_to))
        else:
            protein_matches[entry_acc] = [(pos_from, pos_to)]

    # Last protein
    _condense(protein_matches)
    for _entry_acc, entry_matches in protein_matches.items():
        entry_type, entry_name = entries[_entry_acc]

        for _pos_from, _pos_end in entry_matches:
            table.insert((_protein_acc, _entry_acc, entry_type,
                          entry_name, _pos_from, _pos_end))

    table.close()

    # Removing temporary directory (now empty)
    os.rmdir(tmpdir)

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
