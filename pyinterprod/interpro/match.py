# -*- coding: utf-8 -*-

import json
import os
from typing import Dict

import cx_Oracle

from pyinterprod import logger
from pyinterprod.utils import oracle


_ENTRIES = "entries_proteins.json"
_DATABASES = "databases_matches.json"
ENTRIES = "entries_changes.tsv"
DATABASES = "databases_changes.tsv"


def update_matches(url: str, outdir: str):
    os.makedirs(outdir, exist_ok=True)

    con = cx_Oracle.connect(url)
    prepare_matches(con)
    check_matches(con, outdir)
    insert_matches(con)
    track_count_changes(con, outdir)
    con.close()


def prepare_matches(con: cx_Oracle.Connection):
    cur = con.cursor()

    logger.info("populating MATCH_NEW")
    oracle.drop_table(cur, "INTERPRO.MATCH_NEW", purge=True)

    cur.execute(
        """
        CREATE TABLE INTERPRO.MATCH_NEW NOLOGGING
        AS
        SELECT * 
        FROM INTERPRO.MATCH WHERE 1 = 0
        """
    )

    cur.execute(
        """
        INSERT /*+ APPEND */ INTO INTERPRO.MATCH_NEW
        SELECT
          P.PROTEIN_AC, M.METHOD_AC, M.SEQ_START, M.SEQ_END, 'T',
          D.DBCODE, D.EVIDENCE,
          SYSDATE, SYSDATE, SYSDATE, 'INTERPRO',
          M.EVALUE, M.MODEL_AC, M.FRAGMENTS
        FROM INTERPRO.PROTEIN_TO_SCAN P
        INNER JOIN IPRSCAN.MV_IPRSCAN M
          ON P.UPI = M.UPI
        INNER JOIN INTERPRO.IPRSCAN2DBCODE D
          ON M.ANALYSIS_ID = D.IPRSCAN_SIG_LIB_REL_ID
        WHERE D.DBCODE NOT IN ('g', 'j', 'n', 'q', 's', 'v', 'x')
        AND M.SEQ_START != M.SEQ_END
        """
    )
    con.commit()

    logger.info("indexing MATCH_NEW")
    for col in ("DBCODE", "PROTEIN_AC"):
        cur.execute(
            f"""
            CREATE INDEX I_MATCH_NEW${col}
            ON INTERPRO.MATCH_NEW ({col}) NOLOGGING
            """
        )

    logger.info("gathering statistics")
    oracle.gather_stats(cur, "INTERPRO", "MATCH_NEW")

    logger.info("deleting SUPERFAMILY duplicated matches")
    cur.execute(
        """
        DELETE FROM INTERPRO.MATCH_NEW M1
        WHERE EXISTS(
          SELECT 1
          FROM INTERPRO.MATCH_NEW M2
          WHERE M2.DBCODE = 'Y'
          AND M1.PROTEIN_AC = M2.PROTEIN_AC
          AND M1.METHOD_AC = M2.METHOD_AC
          AND M1.POS_FROM = M2.POS_FROM
          AND M1.POS_TO = M2.POS_TO
          AND M1.SCORE > M2.SCORE
        )
        """
    )
    logger.info(f"{cur.rowcount} SUPERFAMILY matches deleted")
    con.commit()
    cur.close()


def check_matches(con: cx_Oracle.Connection, outdir: str):
    cur = con.cursor()

    # Matches outside of the protein
    logger.info("checking out-of-bound matches")
    cur.execute(
        """
        SELECT M.PROTEIN_AC, M.METHOD_AC, M.POS_TO, P.LEN
        FROM INTERPRO.MATCH_NEW M
        INNER JOIN INTERPRO.PROTEIN P
          ON M.PROTEIN_AC = P.PROTEIN_AC
        WHERE M.POS_TO > P.LEN
        """
    )
    n = 0
    for row in cur:
        logger.critical("out-of-bound: {}\t{}\t{}\t{}".format(*row))
        n += 1

    if n:
        cur.close()
        con.close()
        raise RuntimeError("{} out-of-bound matches".format(n))

    # Matches with invalid start/end positions
    logger.info("checking invalid matches")
    cur.execute(
        """
        SELECT PROTEIN_AC, METHOD_AC, POS_FROM, POS_TO
        FROM INTERPRO.MATCH_NEW
        WHERE POS_FROM < 1 OR POS_FROM > POS_TO
        """
    )
    n = 0
    for row in cur:
        logger.critical("invalid: {}\t{}\t{}\t{}".format(*row))
        n += 1

    if n:
        cur.close()
        con.close()
        raise RuntimeError("{} invalid matches".format(n))

    with open(os.path.join(outdir, _ENTRIES), "wt") as fh:
        json.dump(get_entries_proteins_count(cur), fh)

    with open(os.path.join(outdir, _DATABASES), "wt") as fh:
        json.dump(get_databases_matches_count(cur), fh)

    cur.close()


def insert_matches(con: cx_Oracle.Connection):
    cur = con.cursor()

    logger.info("updating MATCH")
    cur.execute(
        """
        DELETE FROM INTERPRO.MATCH
        WHERE PROTEIN_AC IN (
          SELECT PROTEIN_AC
          FROM INTERPRO.PROTEIN_TO_SCAN
        )
        """
    )
    logger.info(f"{cur.rowcount} rows deleted")

    cur.execute(
        """
        INSERT INTO INTERPRO.MATCH
        SELECT * FROM INTERPRO.MATCH_NEW
        """
    )
    logger.info("{cur.rowcount} rows inserted")
    con.commit()

    oracle.drop_table(cur, "INTERPRO.MATCH_NEW", purge=True)
    cur.close()


def track_count_changes(con: cx_Oracle.Connection, outdir: str):
    logger.info("tracking changes")
    cur = con.cursor()

    with open(os.path.join(outdir, _ENTRIES), "rt") as fh:
        previous = json.load(fh)

    changes = []
    for accession, count in get_entries_proteins_count(cur).items():
        try:
            prev_count = previous[accession]
        except KeyError:
            # This entry did not have matches: skip
            continue

        change = (count - prev_count) / prev_count
        if abs(change) >= 0.5:
            changes.append((accession, prev_count, count, change))

    changes.sort(key=lambda x: abs(x[3]))
    with open(os.path.join(outdir, ENTRIES), "wt") as fh:
        fh.write("# Accession\tPrevious protein count\t"
                 "New protein count\tChange (%)\n")

        for accession, prev_count, count, change in changes:
            fh.write(f"{accession}\t{prev_count}\t{count}\t{change*100:.0f}\n")

    with open(os.path.join(outdir, _DATABASES), "rt") as fh:
        previous = json.load(fh)

    for code, count in get_databases_matches_count(cur).items():
        prev_count = previous.pop(code, 0)
        changes.append((code, prev_count, count))

    # Only if a database does not have any matches any more
    for code, prev_count in previous.items():
        changes.append((code, prev_count, 0))

    with open(os.path.join(outdir, DATABASES), "wt") as fh:
        fh.write("# Code\tPrevious match count\tNew match count\n")

        for code, prev_count, count in sorted(changes):
            fh.write(f"{code}\t{prev_count}\t{count}\n")

    cur.close()


def get_entries_proteins_count(cur: cx_Oracle.Cursor) -> Dict[str, int]:
    cur.execute(
        """
        SELECT E.ENTRY_AC, COUNT(DISTINCT M.PROTEIN_AC)
        FROM INTERPRO.ENTRY2METHOD E
        INNER JOIN INTERPRO.MATCH M
          ON E.METHOD_AC = M.METHOD_AC
        GROUP BY E.ENTRY_AC
        """
    )
    return dict(cur.fetchall())


def get_databases_matches_count(cur: cx_Oracle.Cursor) -> Dict[str, int]:
    cur.execute(
        """
        SELECT DBCODE, COUNT(*)
        FROM INTERPRO.MATCH
        GROUP BY DBCODE
        """
    )
    return dict(cur.fetchall())
