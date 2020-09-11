# -*- coding: utf-8 -*-

import os
import pickle
from typing import Dict

import cx_Oracle

from pyinterprod import logger
from pyinterprod.utils import oracle
from .database import Database


_ENTRIES_TMP_FILE = "entries_proteins.dat"
_DATABASES_TMP_FILE = "databases_matches.dat"
ENTRIES_CHANGES_FILE = "entries_changes.dat"
DATABASES_CHANGES_FILE = "databases_changes.dat"


def update_database_matches(url: str, database: Database, partition: str):
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    oracle.drop_table(cur, "INTERPRO.MATCH_NEW", purge=True)
    cur.execute(
        """
        CREATE TABLE INTERPRO.MATCH_NEW NOLOGGING
        AS SELECT * FROM INTERPRO.MATCH WHERE 1 = 0
        """
    )
    cur.execute(
        """
        INSERT /*+ APPEND */ INTO INTERPRO.MATCH_NEW
        SELECT
          X.AC, M.METHOD_AC, M.SEQ_START, M.SEQ_END, 'T',
          D.DBCODE, D.EVIDENCE,
          SYSDATE, SYSDATE, SYSDATE, 'INTERPRO',
          M.EVALUE, M.MODEL_AC, M.FRAGMENTS
        FROM IPRSCAN.MV_IPRSCAN M
        INNER JOIN UNIPARC.XREF X 
          ON M.UPI = X.UPI
        INNER JOIN INTERPRO.IPRSCAN2DBCODE D
          ON M.ANALYSIS_ID = D.IPRSCAN_SIG_LIB_REL_ID
        WHERE X.DBID IN (2, 3)  -- Swiss-Prot or TrEMBL
        AND X.DELETED = 'N'
        AND M.ANALYSIS_ID = :1 
        AND M.SEQ_START != M.SEQ_END
        """, (database.analysis_id,)
    )
    con.commit()

    # Add constraints/indexes to be able to exchange partition
    logger.debug("creating indexes and constraints")
    cur.execute(
        """
        CREATE INDEX MATCH_NEW$D 
        ON INTERPRO.MATCH_NEW (DBCODE) 
        NOLOGGING
        """
    )
    cur.execute(
        """
        CREATE INDEX MATCH_NEW$E 
        ON INTERPRO.MATCH_NEW (EVIDENCE) 
        NOLOGGING
        """
    )
    cur.execute(
        """
        CREATE INDEX MATCH_NEW$S 
        ON INTERPRO.MATCH_NEW (STATUS) 
        NOLOGGING
        """
    )
    cur.execute(
        """
        CREATE INDEX MATCH_NEW$M 
        ON INTERPRO.MATCH_NEW (METHOD_AC) 
        NOLOGGING
        """
    )
    cur.execute(
        """
        ALTER TABLE INTERPRO.MATCH_NEW
        ADD CONSTRAINT MATCH_NEW$CK1 
        CHECK ( POS_FROM >= 1 ) 
        """
    )
    cur.execute(
        """
        ALTER TABLE INTERPRO.MATCH_NEW
        ADD CONSTRAINT MATCH_NEW$CK2
        CHECK ( POS_TO - POS_FROM > 0 ) 
        """
    )
    cur.execute(
        """
        ALTER TABLE INTERPRO.MATCH_NEW
        ADD CONSTRAINT MATCH_NEW$CK3
        CHECK ( STATUS != 'N' OR (STATUS = 'N' AND DBCODE IN ('P', 'M', 'Q')) )
        """
    )
    cur.execute(
        """
        ALTER TABLE INTERPRO.MATCH_NEW
        ADD CONSTRAINT MATCH_NEW$PK
        PRIMARY KEY (PROTEIN_AC, METHOD_AC, POS_FROM, POS_TO)
        """
    )
    cur.execute(
        f"""
        ALTER TABLE INTERPRO.MATCH_NEW
        ADD CONSTRAINT MATCH_NEW$FK1 FOREIGN KEY (DBCODE) 
        REFERENCES INTERPRO.CV_DATABASE (DBCODE)
        """
    )
    cur.execute(
        """
        ALTER TABLE INTERPRO.MATCH_NEW
        ADD CONSTRAINT MATCH_NEW$FK2 FOREIGN KEY (EVIDENCE) 
        REFERENCES INTERPRO.CV_EVIDENCE (CODE)
        """
    )
    cur.execute(
        f"""
        ALTER TABLE INTERPRO.MATCH_NEW
        ADD CONSTRAINT MATCH_NEW$FK3 FOREIGN KEY (METHOD_AC)
        REFERENCES INTERPRO.METHOD (METHOD_AC)
        """
    )
    cur.execute(
        f"""
        ALTER TABLE INTERPRO.MATCH_NEW
        ADD CONSTRAINT MATCH_NEW$FK4 FOREIGN KEY (PROTEIN_AC) 
        REFERENCES INTERPRO.PROTEIN (PROTEIN_AC)
        """
    )
    cur.execute(
        f"""
        ALTER TABLE INTERPRO.MATCH_NEW
        ADD CONSTRAINT MATCH_NEW$FK5 FOREIGN KEY (STATUS) 
        REFERENCES INTERPRO.CV_STATUS (CODE)
        """
    )
    cur.execute(
        f"""
        ALTER TABLE INTERPRO.MATCH_NEW
        ADD CONSTRAINT MATCH_NEW$CK4 
        CHECK (PROTEIN_AC IS NOT NULL )
        """
    )
    cur.execute(
        f"""
        ALTER TABLE INTERPRO.MATCH_NEW
        ADD CONSTRAINT MATCH_NEW$CK5 
        CHECK (METHOD_AC IS NOT NULL )
        """
    )

    cur.execute("SELECT COUNT(*) FROM INTERPRO.MATCH_NEW")
    cnt, = cur.fetchone()
    if not cnt:
        raise RuntimeError(f"no rows inserted "
                           f"for analysis ID {database.analysis_id}")

    logger.debug(f"exchanging partition")
    cur.execute(
        f"""
        ALTER TABLE INTERPRO.MATCH
        EXCHANGE PARTITION ({partition}) 
        WITH TABLE INTERPRO.MATCH_NEW
        WITHOUT VALIDATION
        """
    )
    oracle.drop_table(cur, "INTERPRO.MATCH_NEW", purge=True)

    cur.close()
    con.close()


def update_matches(url: str, outdir: str):
    """
    Add protein matches for recently added/modified sequences

    :param url: Oracle connection string
    :param outdir: output directory for data files
    """
    os.makedirs(outdir, exist_ok=True)

    con = cx_Oracle.connect(url)
    _prepare_matches(con)
    _check_matches(con, outdir)
    _insert_matches(con)
    _track_count_changes(con, outdir)
    con.close()


def update_feature_matches(url: str):
    """
    Add protein feature matches for recently added/modified sequences

    :param url: Oracle connection string
    """
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
    """
    Recreate splice-variants table with the most recent data

    :param url: Oracle connection string
    """
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
    """
    Add protein site matches for recently added/modified sequences

    :param url: Oracle connection string
    """
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
    """
    Import protein matches in a staging table

    :param con: Oracle connection object
    """
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
    """
    Check there are not errors in imported matches

    :param con: Oracle connection object
    :param outdir: output directory for data files
    """
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

    with open(os.path.join(outdir, _ENTRIES_TMP_FILE), "wb") as fh:
        pickle.dump(_get_entries_proteins_count(cur), fh)

    with open(os.path.join(outdir, _DATABASES_TMP_FILE), "wb") as fh:
        pickle.dump(_get_databases_matches_count(cur), fh)

    cur.close()


def _insert_matches(con: cx_Oracle.Connection):
    """
    Update the MATCH table with data from the staging table

    :param con: Oracle connection object
    """
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
    """
    Evaluate how the protein counts changed between two UniProt versions

    :param con: Oracle connection object
    :param outdir: directory for data files
    """
    logger.info("tracking changes")
    cur = con.cursor()

    # Find entries with significant protein count changes
    with open(os.path.join(outdir, _ENTRIES_TMP_FILE), "rb") as fh:
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

    with open(os.path.join(outdir, ENTRIES_CHANGES_FILE), "wb") as fh:
        pickle.dump(changes, fh)

    # List the previous/new match counts for all databases
    with open(os.path.join(outdir, _DATABASES_TMP_FILE), "rb") as fh:
        previous = pickle.load(fh)

    changes = []
    for code, count in _get_databases_matches_count(cur).items():
        prev_count = previous.pop(code, 0)
        changes.append((code, prev_count, count))

    # Only if a database does not have any matches any more
    for code, prev_count in previous.items():
        changes.append((code, prev_count, 0))

    with open(os.path.join(outdir, DATABASES_CHANGES_FILE), "wb") as fh:
        pickle.dump(changes, fh)

    cur.close()

    os.remove(os.path.join(outdir, _DATABASES_TMP_FILE))
    os.remove(os.path.join(outdir, _ENTRIES_TMP_FILE))


def _get_entries_proteins_count(cur: cx_Oracle.Cursor) -> Dict[str, int]:
    """
    Return the number of protein matched by each InterPro entry

    :param cur: Oracle cursor object
    :return: dictionary
    """
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
    """
    Return the number of matches per member database
    :param cur: Oracle cursor object
    :return: dictionary
    """
    cur.execute(
        """
        SELECT DBCODE, COUNT(*)
        FROM INTERPRO.MATCH
        GROUP BY DBCODE
        """
    )
    return dict(cur.fetchall())
