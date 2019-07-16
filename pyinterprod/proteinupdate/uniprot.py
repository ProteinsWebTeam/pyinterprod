# -*- coding: utf-8 -*-

import os
import re
from typing import Generator, Tuple

import cx_Oracle

from .sendmail import send_mail
from .. import logger, orautils


class Entry(object):
    def __init__(self):
        self.identifier = None
        self.is_reviewed = None
        self.length = None
        self.accession = None
        self.secondary_accessions = []
        self.date = None
        self.is_fragment = None
        self.taxon_id = None
        self.crc64 = None

        self.prog_id = re.compile("^ID\s+(\w+)\s+([a-zA-Z]+);\s+(\d+)",
                                  flags=re.M)
        self.prog_ac = re.compile("^AC\s+(.+;$)", flags=re.M)
        self.prog_ac2 = re.compile("\w+", flags=re.M)
        self.prog_dt = re.compile(
            "DT\s+(\d+-[A-Z]{3}-\d{4}), sequence version \d+\.",
            flags=re.M
        )
        self.prog_de = re.compile("^DE\s+Flags:\s*Fragment;", flags=re.M)
        self.prog_ox = re.compile("^OX\s+NCBI_TaxID=(\d+)", flags=re.M)
        self.prog_sq = re.compile("^SQ.+?(\w+)\s*CRC64;", flags=re.M)

    def update(self, buffer: str) -> tuple:
        self.identifier, status, length = self.prog_id.search(buffer).groups()
        self.is_reviewed = status == "Reviewed"
        self.length = int(length)

        accessions = []
        for m in self.prog_ac.finditer(buffer):
            for acc in self.prog_ac2.findall(m.group(1)):
                accessions.append(acc)
        self.accession, *self.secondary_accessions = accessions

        # date_string = self.prog_dt.search(buffer).group(1)
        # self.date = datetime.strptime(date_string, "%d-%b-%Y")

        self.is_fragment = self.prog_de.search(buffer) is not None
        self.taxon_id = int(self.prog_ox.search(buffer).group(1))

        self.crc64 = self.prog_sq.search(buffer).group(1)

        return (
            self.accession,
            self.identifier,
            1 if self.is_reviewed else 0,
            self.crc64,
            self.length,
            1 if self.is_fragment else 0,
            self.taxon_id
        )


def read_flat_file(path: str) -> Generator:
    with open(path, "rt") as fh:
        buffer = ""
        e = Entry()
        for line in fh:
            if line[:2] == "//":
                yield e.update(buffer)
                buffer = ""
            else:
                buffer += line

        if buffer:  # in case the last line of the file is not //
            yield e.update(buffer)


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
                   [----] [----]
                          s
                -> new match
                """
                entry_matches.append((start, end))
                start = s
                end = e
            elif e > end:
                """
                        end
                   [----]
                     [------]
                            e
                -> extend
                """
                end = e

        entry_matches.append((start, end))
        matches[entry_acc] = entry_matches


def build_xref_condensed(user: str, dsn: str):
    logger.info("building XREF_CONDENSED")
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()

    orautils.drop_table(cur, "INTERPRO", "XREF_CONDENSED", purge=True)
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

    table = orautils.TablePopulator(con,
                                    query="INSERT /*+ APPEND */ "
                                          "INTO INTERPRO.XREF_CONDENSED "
                                          "VALUES (:1, :2, :3, :4, :5, :6)",
                                    autocommit=True)

    cur.execute(
        """
        SELECT PROTEIN_AC, METHOD_AC, POS_FROM, POS_TO, FRAGMENTS
        FROM INTERPRO.MATCH
        ORDER BY PROTEIN_AC
        """
    )
    _protein_acc = None
    protein_matches = []
    for protein_acc, method_acc, pos_from, pos_to, fragments in cur:
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

        try:
            entry_acc = signatures[method_acc]
        except KeyError:
            # Signature not integrated or integrated in an unchecked entry
            continue

        if entry_acc in protein_matches:
            e = protein_matches[entry_acc]
        else:
            e = protein_matches[entry_acc] = []

        """
        As of May 2019, UniProt does not use discontinuous domains
        because their collaborators need to be able to
        distinguish between repeated matches and fragmented matches
        """
        # if fragments is not None:
        #     for frag in fragments.split(','):
        #         start, end, _ = frag.split('-')
        #         e.append((int(start), int(end)))
        # else:
        #     e.append((pos_from, pos_to))
        e.append((pos_from, pos_to))

    # Last protein
    _condense(protein_matches)
    for _entry_acc, entry_matches in protein_matches.items():
        entry_type, entry_name = entries[_entry_acc]

        for _pos_from, _pos_end in entry_matches:
            table.insert((_protein_acc, _entry_acc, entry_type,
                          entry_name, _pos_from, _pos_end))

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
    orautils.drop_table(cur, "INTERPRO", "XREF_SUMMARY", purge=True)
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
    cur.execute("SELECT VERSION FROM INTERPRO.DB_VERSION WHERE DBCODE = 'u'")
    release = cur.fetchone()[0]
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
        # 'D': ("ProDom", "PD"),  # ProDom removed from InterPro
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
            fh2.write(
                "{}   DR   {}; {}; -; {}.\n".format(*_convert_gene3d(row)))
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

    send_mail(
        to_addrs=["uniprot-database@ebi.ac.uk"],
        cc_addrs=["uniprot-prod@ebi.ac.uk"],
        subject="InterPro XREF files are ready",
        content="""\
Dear UniProt team,

The InterPro cross-references files for {} are available in the following directory:
  {}

Best regards,
The InterPro Production Team
        """.format(release, dst)
    )

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
    orautils.drop_table(cur, "IPRSCAN", "AA_IPRSCAN", purge=True)

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


def ask_to_snapshot(user: str, dsn: str):
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    cur.execute("SELECT VERSION FROM INTERPRO.DB_VERSION WHERE DBCODE = 'u'")
    release = cur.fetchone()[0]
    cur.execute("SELECT MAX(UPI) FROM UNIPARC.PROTEIN@UAREAD")
    max_upi = cur.fetchone()[0]
    cur.execute(
        """
        SELECT I.IPRSCAN_SIG_LIB_REL_ID, D.DBNAME, V.VERSION
        FROM INTERPRO.IPRSCAN2DBCODE I
        INNER JOIN INTERPRO.CV_DATABASE D ON I.DBCODE = D.DBCODE
        LEFT OUTER JOIN INTERPRO.DB_VERSION V ON I.DBCODE = V.DBCODE
        WHERE I.IPRSCAN_SIG_LIB_REL_ID != I.PREV_ID
        ORDER BY D.DBNAME
        """
    )
    updates = []
    for analysis_id, name, version in cur:
        if version:
            name += " (" + version + ")"
        updates.append((name, analysis_id))

    content = """\
INTERPRO.XREF_SUMMARY and INTERPRO.XREF_CONDENSED are ready.

Please, take a snapshot of IPPRO, and inform UniProt they can refresh IPREADU.

Recipients
----------
  To: uniprot-database@ebi.ac.uk
  Cc: automated_annotation@ebi.ac.uk, uniprot-prod@ebi.ac.uk

Subject
-------
Protein update {} completed: please refresh IPREADU

Body (change text between square brackets)
------------------------------------------
Dear UniProt team,

You may refresh IPREADU with the snapshot of IPPRO from [date/time].

The max UPI we processed up to in this update is {}.
""".format(release, max_upi)

    if updates:
        content += """
We have updated the following databases since the previous protein update:
"""
        for db in updates:
            content += "  - {:<30}new analysis ID: {}\n".format(*db)
    else:
        content += """
There are no changes to the analysis IDs for this protein update.
"""

    content += """
Best regards,
The InterPro Production Team
"""

    try:
        send_mail(
            to_addrs=["interpro-team@ebi.ac.uk"],
            subject="Protein update {}: please snapshot IPPRO".format(release),
            content=content
        )

        # Update PREV_ID only if mail sent
        cur.execute(
            """
            UPDATE INTERPRO.IPRSCAN2DBCODE
            SET PREV_ID = IPRSCAN_SIG_LIB_REL_ID
            """
        )
        con.commit()
    finally:
        cur.close()
        con.close()


def report_integration_changes(user: str, dsn: str):
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    cur.execute(
        """
        SELECT DISTINCT METHOD_AC, ENTRY_AC
        FROM INTERPRO.XREF_SUMMARY
        """
    )
    previous_signatures = dict(cur.fetchall())
    cur.execute(
        """
        SELECT M.METHOD_AC, EM.ENTRY_AC
        FROM INTERPRO.METHOD M
        LEFT OUTER JOIN INTERPRO.ENTRY2METHOD EM
          ON M.METHOD_AC = EM.METHOD_AC
        """
    )
    current_signatures = dict(cur.fetchall())
    cur.execute(
        """
        SELECT ENTRY_AC, CHECKED
        FROM INTERPRO.ENTRY
        """
    )
    current_entries = {row[0]: row[1] == 'Y' for row in cur}
    cur.execute("SELECT VERSION FROM INTERPRO.DB_VERSION WHERE DBCODE = 'u'")
    release = cur.fetchone()[0]
    cur.close()
    con.close()

    changes = {}
    for s_acc, e_acc in previous_signatures.items():
        if s_acc in current_signatures:
            # Signature still exists
            ne_acc = current_signatures.pop(s_acc)
            if ne_acc:
                # Signature still integrated
                if e_acc != ne_acc:
                    # In another entry
                    if current_entries[ne_acc]:
                        # The other entry is checked
                        changes[s_acc] = ("moved", e_acc, ne_acc)
                    else:
                        # The other entry is unchecked
                        changes[s_acc] = ("moved to unchecked", e_acc, ne_acc)
                elif current_entries[e_acc]:
                    pass  # No changes
                else:
                    # In the same (now unchecked) entry
                    changes[s_acc] = ("entry unchecked", e_acc, '')
            else:
                # Signature unintegrated
                changes[s_acc] = ("unintegrated", e_acc, '')
        else:
            # Signature does not exist any more (member database update)
            changes[s_acc] = ("deleted", e_acc, '')

    for s_acc, ne_acc in current_signatures.items():
        if ne_acc and current_entries[ne_acc]:
            changes[s_acc] = ("integrated", '', ne_acc)

    content = """\
Dear UniProt team,

Please find below the list of recent integration changes.

"""
    content += "{:<30}{:<30}{:<20}{:<20}\n".format("Signature", "Comment",
                                                   "Previous entry",
                                                   "Current entry")
    content += '-' * 100 + '\n'

    for acc in sorted(changes, key=lambda k: k.lower()):
        content += "{:<30}{:<30}{:<20}{:<20}\n".format(acc, *changes[acc])

    content += "\nBest regards,\nThe InterPro Production Team\n"

    send_mail(
        to_addrs=["aa_dev@ebi.ac.uk"],
        cc_addrs=["unirule@ebi.ac.uk"],
        bcc_addrs=["interpro-team@ebi.ac.uk"],
        subject="Protein update {}: integration changes".format(release),
        content=content
    )
