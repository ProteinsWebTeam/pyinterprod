# -*- coding: utf-8 -*-

import os
import pickle
from typing import Dict

import cx_Oracle

from pyinterprod import logger
from pyinterprod.utils import oracle


_ENTRIES_DAT = "entries_proteins.dat"
_DATABASES_DAT = "databases_matches.dat"
ENTRIES_TSV = "entries_changes.tsv"
DATABASES_TSV = "databases_changes.tsv"


def update_matches(url: str, outdir: str):
    os.makedirs(outdir, exist_ok=True)

    con = cx_Oracle.connect(url)
    _prepare_matches(con)
    _check_matches(con, outdir)
    _insert_matches(con)
    _track_count_changes(con, outdir)
    con.close()


def update_feature_matches(url: str):
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    logger.info("updating FEATURE_MATCH")
    cur.execute(
        """
        DELETE FROM INTERPRO.FEATURE_MATCH
        WHERE PROTEIN_AC IN (
          SELECT PROTEIN_AC
          FROM INTERPRO.PROTEIN_TO_SCAN
        )
        """
    )
    logger.info(f"{cur.rowcount} rows deleted")

    cur.execute(
        """
        INSERT INTO INTERPRO.FEATURE_MATCH
        SELECT
          P.PROTEIN_AC, M.METHOD_AC, M.SEQ_FEATURE, M.SEQ_START, M.SEQ_END,
          D.DBCODE
        FROM INTERPRO.PROTEIN_TO_SCAN P
        INNER JOIN IPRSCAN.MV_IPRSCAN M
          ON P.UPI = M.UPI
        INNER JOIN INTERPRO.IPRSCAN2DBCODE D
          ON M.ANALYSIS_ID = D.IPRSCAN_SIG_LIB_REL_ID
        -- Include MobiDB-Lite, Phobius, SignalP (Euk, Gram+, Gram-), TMHMM, COILS
        WHERE D.DBCODE IN ('g', 'j', 'n', 's', 'v', 'q', 'x')
        """
    )
    logger.info(f"{cur.rowcount} rows inserted")
    con.commit()

    cur.close()
    con.close()


def update_variant_matches(url: str):
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    logger.info("updating VARSPLIC_MASTER")
    oracle.truncate_table(cur, "INTERPRO.VARSPLIC_MASTER", reuse_storage=True)
    cur.execute(
        """
        INSERT INTO INTERPRO.VARSPLIC_MASTER
        SELECT 
          SUBSTR(X.AC, 1, INSTR(X.AC, '-') - 1),
          SUBSTR(X.AC, INSTR(X.AC, '-') + 1),
          P.CRC64,
          P.LEN
        FROM UNIPARC.XREF X
        INNER JOIN UNIPARC.PROTEIN P ON X.UPI = P.UPI
        WHERE X.DBID IN (24, 25) -- SWISSPROT_VARSPLIC, TREMBL_VARSPLIC
        AND X.DELETED = 'N'
        """
    )
    logger.info(f"{cur.rowcount} rows inserted")

    logger.info("updating VARSPLIC_MATCH")
    oracle.truncate_table(cur, "INTERPRO.VARSPLIC_MATCH", reuse_storage=True)
    cur.execute(
        """
        INSERT INTO INTERPRO.VARSPLIC_MATCH
        SELECT 
          X.AC, MV.METHOD_AC, MV.SEQ_START, MV.SEQ_END, 'T' AS STATUS, 
          I2D.DBCODE, I2D.EVIDENCE, SYSDATE, SYSDATE, SYSDATE, 
          'INTERPRO', MV.EVALUE, MV.MODEL_AC, MV.FRAGMENTS
        FROM UNIPARC.XREF X
        INNER JOIN IPRSCAN.MV_IPRSCAN MV
          ON X.UPI = MV.UPI
        INNER JOIN INTERPRO.IPRSCAN2DBCODE I2D 
          ON MV.ANALYSIS_ID = I2D.IPRSCAN_SIG_LIB_REL_ID
        INNER JOIN INTERPRO.METHOD M
          ON MV.METHOD_AC = M.METHOD_AC
        WHERE X.DBID IN (24, 25)  -- SWISSPROT_VARSPLIC, TREMBL_VARSPLIC
        AND X.DELETED = 'N'
        -- Exclude MobiDB-Lite, Phobius, SignalP (Euk, Gram+, Gram-), TMHMM, COILS
        AND I2D.DBCODE NOT IN ('g', 'j', 'n', 'q', 's', 'v', 'x')
        AND M.SKIP_FLAG = 'N'
        """
    )
    logger.info(f"{cur.rowcount} rows inserted")
    con.commit()

    cur.close()
    con.close()


def update_site_matches(url: str):
    con = cx_Oracle.connect(url)
    cur = con.cursor()

    logger.info("populating SITE_MATCH_NEW")
    oracle.drop_table(cur, "INTERPRO.SITE_MATCH_NEW", purge=True)

    cur.execute(
        """
        CREATE TABLE INTERPRO.SITE_MATCH_NEW
        AS
        SELECT
            P.PROTEIN_AC, S.METHOD_AC, S.LOC_START, S.LOC_END, S.DESCRIPTION,
            S.RESIDUE, S.RESIDUE_START, S.RESIDUE_END, S.NUM_SITES, D.DBCODE
        FROM INTERPRO.PROTEIN_TO_SCAN P
        INNER JOIN IPRSCAN.SITE S
          ON P.UPI = S.UPI
        INNER JOIN INTERPRO.IPRSCAN2DBCODE D
          ON S.ANALYSIS_ID = D.IPRSCAN_SIG_LIB_REL_ID
        """
    )

    logger.info("indexing")
    cur.execute(
        """
        CREATE INDEX I_SITE_NEW
        ON INTERPRO.SITE_MATCH_NEW (PROTEIN_AC, METHOD_AC, LOC_START, LOC_END)
        """
    )

    logger.info("gathering statistics")
    oracle.gather_stats(cur, "INTERPRO", "SITE_MATCH_NEW")

    logger.info("checking")
    cur.execute(
        """
        SELECT COUNT(*)
        FROM (
            SELECT DISTINCT PROTEIN_AC, METHOD_AC, LOC_START, LOC_END
            FROM INTERPRO.SITE_MATCH_NEW
            MINUS (
              SELECT DISTINCT PROTEIN_AC, METHOD_AC, POS_FROM, POS_TO
              FROM INTERPRO.MATCH PARTITION (MATCH_DBCODE_J)
              UNION ALL
              SELECT DISTINCT PROTEIN_AC, METHOD_AC, POS_FROM, POS_TO
              FROM INTERPRO.MATCH PARTITION (MATCH_DBCODE_B)
            )
        )
        """
    )

    cnt, = cur.fetchone()
    if cnt:
        cur.close()
        con.close()
        raise RuntimeError(f"{cnt} matches in SITE_MATCH_NEW "
                           f"that are not in MATCH")

    logger.info("updating SITE_MATCH")
    cur.execute(
        """
        DELETE FROM INTERPRO.SITE_MATCH
        WHERE PROTEIN_AC IN (
          SELECT PROTEIN_AC
          FROM INTERPRO.PROTEIN_TO_SCAN
        )
        """
    )
    logger.info(f"{cur.rowcount} rows deleted")

    cur.execute(
        """
        INSERT INTO INTERPRO.SITE_MATCH
        SELECT * FROM INTERPRO.SITE_MATCH_NEW
        """
    )
    logger.info(f"{cur.rowcount} rows inserted")
    con.commit()

    oracle.drop_table(cur, "INTERPRO.SITE_MATCH_NEW", purge=True)

    cur.close()
    con.close()


def _prepare_matches(con: cx_Oracle.Connection):
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
        -- Exclude MobiDB-Lite, Phobius, SignalP (Euk, Gram+, Gram-), TMHMM, COILS
        WHERE D.DBCODE NOT IN ('g', 'j', 'n', 's', 'v', 'q', 'x')
        AND M.SEQ_START != M.SEQ_END
        """
    )
    con.commit()

    logger.info("indexing")
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


def _check_matches(con: cx_Oracle.Connection, outdir: str):
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
    cnt = 0
    for row in cur:
        logger.critical("out-of-bound: {}\t{}\t{}\t{}".format(*row))
        cnt += 1

    if cnt:
        cur.close()
        con.close()
        raise RuntimeError(f"{cnt} out-of-bound matches")

    # Matches with invalid start/end positions
    logger.info("checking invalid matches")
    cur.execute(
        """
        SELECT PROTEIN_AC, METHOD_AC, POS_FROM, POS_TO
        FROM INTERPRO.MATCH_NEW
        WHERE POS_FROM < 1 OR POS_FROM > POS_TO
        """
    )
    cnt = 0
    for row in cur:
        logger.critical("invalid: {}\t{}\t{}\t{}".format(*row))
        cnt += 1

    if cnt:
        cur.close()
        con.close()
        raise RuntimeError(f"{cnt} invalid matches")

    with open(os.path.join(outdir, _ENTRIES_DAT), "wb") as fh:
        pickle.dump(_get_entries_proteins_count(cur), fh)

    with open(os.path.join(outdir, _DATABASES_DAT), "wb") as fh:
        pickle.dump(_get_databases_matches_count(cur), fh)

    cur.close()


def _insert_matches(con: cx_Oracle.Connection):
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
    logger.info(f"{cur.rowcount} rows inserted")
    con.commit()

    oracle.drop_table(cur, "INTERPRO.MATCH_NEW", purge=True)
    cur.close()


def _track_count_changes(con: cx_Oracle.Connection, outdir: str):
    logger.info("tracking changes")
    cur = con.cursor()

    with open(os.path.join(outdir, _ENTRIES_DAT), "rb") as fh:
        previous = pickle.load(fh)

    changes = []
    for accession, count in _get_entries_proteins_count(cur).items():
        try:
            prev_count = previous[accession]
        except KeyError:
            # This entry did not have matches: skip
            continue

        change = (count - prev_count) / prev_count
        if abs(change) >= 0.5:
            changes.append((accession, prev_count, count, change))

    changes.sort(key=lambda x: abs(x[3]))
    with open(os.path.join(outdir, ENTRIES_TSV), "wt") as fh:
        fh.write("# Accession\tPrevious protein count\t"
                 "New protein count\tChange (%)\n")

        for accession, prev_count, count, change in changes:
            fh.write(f"{accession}\t{prev_count}\t{count}\t{change*100:.0f}\n")

    with open(os.path.join(outdir, _ENTRIES_DAT), "rb") as fh:
        previous = pickle.load(fh)

    changes = []
    for code, count in _get_databases_matches_count(cur).items():
        prev_count = previous.pop(code, 0)
        changes.append((code, prev_count, count))

    # Only if a database does not have any matches any more
    for code, prev_count in previous.items():
        changes.append((code, prev_count, 0))

    with open(os.path.join(outdir, DATABASES_TSV), "wt") as fh:
        fh.write("# Code\tPrevious match count\tNew match count\n")

        for code, prev_count, count in sorted(changes):
            fh.write(f"{code}\t{prev_count}\t{count}\n")

    cur.close()


def _get_entries_proteins_count(cur: cx_Oracle.Cursor) -> Dict[str, int]:
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


def _get_databases_matches_count(cur: cx_Oracle.Cursor) -> Dict[str, int]:
    cur.execute(
        """
        SELECT DBCODE, COUNT(*)
        FROM INTERPRO.MATCH
        GROUP BY DBCODE
        """
    )
    return dict(cur.fetchall())
