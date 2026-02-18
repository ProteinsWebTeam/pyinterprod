import os
import pickle

from collections import defaultdict

import oracledb
import psycopg

from pyinterprod import logger
from pyinterprod.utils import oracle
from pyinterprod.utils.pg import url2dict
from .contrib import toad
from .database import Database


FILE_ENTRY_PROT_COUNTS = "entries.prot.counts.pickle"


def export_entries_protein_counts(cur: oracledb.Cursor, pg_url: str, data_dir: str):
    with open(os.path.join(data_dir, FILE_ENTRY_PROT_COUNTS), "wb") as fh:
        pickle.dump(_get_entries_protein_counts(cur, pg_url), fh)


def update_database_matches(uri: str, databases: list[Database]):
    """

    :param uri:
    :param databases: list of Database objects
    :return:
    """
    con = oracledb.connect(uri)
    cur = con.cursor()

    partitions = {}
    for p in oracle.get_partitions(cur, "INTERPRO", "MATCH"):
        dbcode = p["value"][1:-1]  # 'X' -> X
        partitions[dbcode] = p["name"]

    for database in databases:
        logger.info(f"{database.name}")
        oracle.drop_table(cur, "INTERPRO.MATCH_NEW", purge=True)
        logger.info(f"\tpopulating MATCH_NEW "
                    f"(ANALYSIS_ID: {database.analysis_id})")
        cur.execute(
            """
            CREATE TABLE INTERPRO.MATCH_NEW NOLOGGING
            AS SELECT * FROM INTERPRO.MATCH WHERE 1 = 0
            """
        )

        if database.identifier == "V":
            # PANTHER: import annotation node ID
            feature = "M.SEQ_FEATURE"
        else:
            feature = "NULL"

        cur.execute(
            f"""
            INSERT /*+ APPEND */ INTO INTERPRO.MATCH_NEW
            SELECT
              X.AC, M.METHOD_AC, M.SEQ_START, M.SEQ_END, 'T',
              D.DBCODE, D.EVIDENCE,
              SYSDATE, SYSDATE, SYSDATE, 'INTERPRO',
              M.EVALUE, M.MODEL_AC, M.FRAGMENTS, {feature}
            FROM IPRSCAN.MV_IPRSCAN M
            INNER JOIN UNIPARC.XREF X
              ON M.UPI = X.UPI
            INNER JOIN INTERPRO.IPRSCAN2DBCODE D
              ON M.ANALYSIS_ID = D.IPRSCAN_SIG_LIB_REL_ID
            WHERE M.ANALYSIS_ID = :1
            AND M.SEQ_START != M.SEQ_END
            AND X.DBID IN (2, 3)  -- Swiss-Prot or TrEMBL
            AND X.DELETED = 'N'
            """,
            [database.analysis_id]
        )
        con.commit()

        # Add constraints/indexes to be able to exchange partition
        logger.info("\tcreating indexes and constraints")
        for col in ("PROTEIN_AC", "METHOD_AC", "STATUS", "DBCODE", "EVIDENCE"):
            cur.execute(
                f"""
                CREATE INDEX MATCH_NEW${col[0]}
                ON INTERPRO.MATCH_NEW ({col})
                TABLESPACE INTERPRO_IND
                NOLOGGING
                """
            )

        cur.execute(
            """
            ALTER TABLE INTERPRO.MATCH_NEW
            ADD CONSTRAINT CK_MATCH_NEW$FROM
            CHECK (POS_FROM >= 1)
            """
        )
        cur.execute(
            """
            ALTER TABLE INTERPRO.MATCH_NEW
            ADD CONSTRAINT CK_MATCH_NEW$NEG
            CHECK (POS_TO - POS_FROM > 0)
            """
        )
        cur.execute(
            """
            ALTER TABLE INTERPRO.MATCH_NEW
            ADD CONSTRAINT CK_MATCH_NEW$STATUS
            CHECK (STATUS != 'N' OR (STATUS = 'N' AND DBCODE IN ('P','M','Q')))
            """
        )
        cur.execute(
            """
            ALTER TABLE INTERPRO.MATCH_NEW
            ADD CONSTRAINT PK_MATCH_NEW
            PRIMARY KEY (PROTEIN_AC, METHOD_AC, POS_FROM, POS_TO)
            """
        )
        cur.execute(
            """
            ALTER TABLE INTERPRO.MATCH_NEW
            ADD CONSTRAINT FK_MATCH_NEW$DBCODE
            FOREIGN KEY (DBCODE) REFERENCES INTERPRO.CV_DATABASE (DBCODE)
            """
        )
        cur.execute(
            """
            ALTER TABLE INTERPRO.MATCH_NEW
            ADD CONSTRAINT FK_MATCH_NEW$EVI
            FOREIGN KEY (EVIDENCE) REFERENCES INTERPRO.CV_EVIDENCE (CODE)
            """
        )
        cur.execute(
            """
            ALTER TABLE INTERPRO.MATCH_NEW
            ADD CONSTRAINT FK_MATCH_NEW$METHOD
            FOREIGN KEY (METHOD_AC) REFERENCES INTERPRO.METHOD (METHOD_AC)
            """
        )
        cur.execute(
            """
            ALTER TABLE INTERPRO.MATCH_NEW
            ADD CONSTRAINT FK_MATCH_NEW$PROTEIN
            FOREIGN KEY (PROTEIN_AC) REFERENCES INTERPRO.PROTEIN (PROTEIN_AC)
            """
        )
        cur.execute(
            """
            ALTER TABLE INTERPRO.MATCH_NEW
            ADD CONSTRAINT FK_MATCH_NEW$STATUS
            FOREIGN KEY (STATUS) REFERENCES INTERPRO.CV_STATUS (CODE)
            """
        )
        cur.execute(
            """
            ALTER TABLE INTERPRO.MATCH_NEW
            ADD CONSTRAINT CK_MATCH_NEW$PROTEIN
            CHECK (PROTEIN_AC IS NOT NULL)
            """
        )
        cur.execute(
            """
            ALTER TABLE INTERPRO.MATCH_NEW
            ADD CONSTRAINT CK_MATCH_NEW$METHOD
            CHECK (METHOD_AC IS NOT NULL)
            """
        )

        cur.execute("SELECT COUNT(*) FROM INTERPRO.MATCH_NEW")
        cnt, = cur.fetchone()
        if not cnt:
            raise RuntimeError(f"no rows inserted "
                               f"for analysis ID {database.analysis_id}")

        logger.info(f"\texchanging partition")
        partition = partitions[database.identifier]
        cur.execute(
            f"""
            ALTER TABLE INTERPRO.MATCH
            EXCHANGE PARTITION ({partition})
            WITH TABLE INTERPRO.MATCH_NEW
            """
        )
        oracle.drop_table(cur, "INTERPRO.MATCH_NEW", purge=True)

        # logger.info("\tgathering statistics")
        # oracle.gather_stats(cur, "INTERPRO", "MATCH", partition)

    cur.close()
    con.close()

    logger.info("complete")


def update_database_feature_matches(uri: str, databases: list[Database]):
    """

    :param uri:
    :param databases: list of Database objects
    :return:
    """
    con = oracledb.connect(uri)
    cur = con.cursor()

    partitions = {}
    for p in oracle.get_partitions(cur, "INTERPRO", "FEATURE_MATCH"):
        dbcode = p["value"][1:-1]  # 'X' -> X
        partitions[dbcode] = p["name"]

    for database in databases:
        logger.info(f"{database.name}")
        oracle.drop_table(cur, "INTERPRO.FEATURE_MATCH_NEW", purge=True)
        cur.execute(
            """
            CREATE TABLE INTERPRO.FEATURE_MATCH_NEW NOLOGGING
            AS SELECT * FROM INTERPRO.FEATURE_MATCH WHERE 1 = 0
            """
        )

        if database.identifier == "d":
            # Pfam-N matches updated in update-features task
            cur.execute(
                """
                INSERT /*+ APPEND */ INTO INTERPRO.FEATURE_MATCH_NEW
                SELECT *
                FROM (
                    SELECT M.PROTEIN_ID,
                           M.METHOD_AC,
                           NULL,
                           M.POS_FROM,
                           CASE WHEN M.POS_TO > P.LEN 
                                THEN P.LEN ELSE M.POS_TO END POS_TO,
                           :1
                    FROM INTERPRO.PFAMN_MATCH M
                    INNER JOIN INTERPRO.PROTEIN P 
                        ON M.PROTEIN_ID = P.PROTEIN_AC
                )
                WHERE POS_TO >= POS_FROM
                """,
                [database.identifier]
            )
            con.commit()
        elif database.identifier == "l":
            # ELM matches updated in update-features task
            cur.execute(
                """
                INSERT /*+ APPEND */ INTO INTERPRO.FEATURE_MATCH_NEW
                SELECT M.PROTEIN_ID,M.METHOD_AC,NULL,M.POS_FROM,M.POS_TO,:1
                FROM INTERPRO.ELM_MATCH M
                WHERE EXISTS (
                    SELECT 1
                    FROM INTERPRO.PROTEIN P
                    WHERE M.PROTEIN_ID = P.PROTEIN_AC
                )                
                """,
                [database.identifier]
            )
            con.commit()
        else:
            logger.info(f"\tpopulating FEATURE_MATCH_NEW "
                        f"(ANALYSIS_ID: {database.analysis_id})")
            cur.execute(
                """
                INSERT /*+ APPEND */ INTO INTERPRO.FEATURE_MATCH_NEW
                SELECT
                  X.AC, M.METHOD_AC, M.SEQ_FEATURE, M.SEQ_START, M.SEQ_END,
                  D.DBCODE
                FROM IPRSCAN.MV_IPRSCAN M
                INNER JOIN UNIPARC.XREF X
                  ON M.UPI = X.UPI
                INNER JOIN INTERPRO.IPRSCAN2DBCODE D
                  ON M.ANALYSIS_ID = D.IPRSCAN_SIG_LIB_REL_ID
                WHERE M.ANALYSIS_ID = :1
                AND X.DBID IN (2, 3)  -- Swiss-Prot or TrEMBL
                AND X.DELETED = 'N'
                """,
                [database.analysis_id]
            )
            con.commit()

        # Add indexes to be able to exchange partition
        logger.info("\tcreating constraints")
        cur.execute(
            """
            ALTER TABLE INTERPRO.FEATURE_MATCH_NEW
            ADD CONSTRAINT CK_FMATCH_NEW$FROM
            CHECK ( POS_FROM >= 1 )
            """
        )
        cur.execute(
            """
            ALTER TABLE INTERPRO.FEATURE_MATCH_NEW
            ADD CONSTRAINT CK_FMATCH_NEW$NEG
            CHECK ( POS_TO >= POS_FROM )
            """
        )
        cur.execute(
            """
            ALTER TABLE INTERPRO.FEATURE_MATCH_NEW
            ADD CONSTRAINT PK_FMATCH_NEW
            PRIMARY KEY (PROTEIN_AC, METHOD_AC, POS_FROM, POS_TO, DBCODE)
            """
        )
        cur.execute(
            f"""
            ALTER TABLE INTERPRO.FEATURE_MATCH_NEW
            ADD CONSTRAINT FK_FMATCH_NEW$M$D
            FOREIGN KEY (METHOD_AC, DBCODE)
                REFERENCES INTERPRO.FEATURE_METHOD (METHOD_AC, DBCODE)
            """
        )
        cur.execute(
            """
            ALTER TABLE INTERPRO.FEATURE_MATCH_NEW
            ADD CONSTRAINT FK_FMATCH_NEW$P
            FOREIGN KEY (PROTEIN_AC) 
                REFERENCES INTERPRO.PROTEIN (PROTEIN_AC)
            """
        )

        cur.execute("SELECT COUNT(*) FROM INTERPRO.FEATURE_MATCH_NEW")
        cnt, = cur.fetchone()
        if not cnt:
            raise RuntimeError(f"no rows inserted "
                               f"for analysis ID {database.analysis_id}")

        logger.info(f"\texchanging partition")
        partition = partitions[database.identifier]
        cur.execute(
            f"""
            ALTER TABLE INTERPRO.FEATURE_MATCH
            EXCHANGE PARTITION ({partition})
            WITH TABLE INTERPRO.FEATURE_MATCH_NEW
            """
        )
        oracle.drop_table(cur, "INTERPRO.FEATURE_MATCH_NEW", purge=True)

        logger.info("\tgathering statistics")
        oracle.gather_stats(cur, "INTERPRO", "FEATURE_MATCH", partition)

    cur.close()
    con.close()

    logger.info("complete")


def update_database_site_matches(uri: str, databases: list[Database]):
    """

    :param uri:
    :param databases: list of Database objects
    :return:
    """
    con = oracledb.connect(uri)
    cur = con.cursor()

    match_partitions = {}
    for p in oracle.get_partitions(cur, "INTERPRO", "MATCH"):
        dbcode = p["value"][1:-1]  # 'X' -> X
        match_partitions[dbcode] = p["name"]

    site_partitions = {}
    for p in oracle.get_partitions(cur, "INTERPRO", "SITE_MATCH"):
        dbcode = p["value"][1:-1]
        site_partitions[dbcode] = p["name"]

    for database in databases:
        logger.info(database.name)

        oracle.drop_table(cur, "INTERPRO.SITE_MATCH_NEW", purge=True)
        cur.execute(
            """
            CREATE TABLE INTERPRO.SITE_MATCH_NEW NOLOGGING
            AS SELECT * FROM INTERPRO.SITE_MATCH WHERE 1 = 0
            """
        )

        logger.info(f"\tinserting site matches")
        cur.execute(
            f"""
            INSERT /*+ APPEND */ INTO INTERPRO.SITE_MATCH_NEW
            SELECT
                X.AC, S.METHOD_AC, D.DBCODE, S.LOC_START, S.LOC_END, 
                S.DESCRIPTION, S.RESIDUE, S.RESIDUE_START, S.RESIDUE_END, 
                S.NUM_SITES
            FROM IPRSCAN.SITE S
            INNER JOIN UNIPARC.XREF X
              ON S.UPI = X.UPI
            INNER JOIN INTERPRO.IPRSCAN2DBCODE D
              ON S.ANALYSIS_ID = D.IPRSCAN_SIG_LIB_REL_ID
            WHERE S.ANALYSIS_ID = :1
            AND X.DBID IN (2, 3)  -- Swiss-Prot or TrEMBL
            AND X.DELETED = 'N'
            """,
            [database.analysis_id]
        )
        con.commit()

        logger.info(f"\tindexing")
        cur.execute(
            """
            CREATE INDEX I_SITE_MATCH_NEW
            ON INTERPRO.SITE_MATCH_NEW (
                PROTEIN_AC, METHOD_AC, LOC_START, LOC_END
            )
            TABLESPACE INTERPRO_IND
            NOLOGGING
            """
        )

        if database.identifier in match_partitions:
            logger.info(f"\tchecking matches")
            partition = match_partitions[database.identifier]
            cur.execute(
                f"""
                SELECT COUNT(*)
                FROM (
                    SELECT DISTINCT PROTEIN_AC, METHOD_AC, LOC_START, LOC_END
                    FROM INTERPRO.SITE_MATCH_NEW
                    MINUS
                    SELECT DISTINCT PROTEIN_AC, METHOD_AC, POS_FROM, POS_TO
                    FROM INTERPRO.MATCH PARTITION ({partition})
                )
                """
            )

            cnt, = cur.fetchone()
            if cnt:
                cur.close()
                con.close()
                raise RuntimeError(f"{database.name}: {cnt} matches "
                                   f"in SITE_MATCH_NEW that are not in MATCH")

        logger.info(f"\tadding constraint")
        cur.execute(
            """
            ALTER TABLE INTERPRO.SITE_MATCH_NEW
            ADD CONSTRAINT FK_SITE_MATCH_NEW$P
            FOREIGN KEY (PROTEIN_AC) REFERENCES INTERPRO.PROTEIN (PROTEIN_AC)
            """
        )
        cur.execute(
            """
            ALTER TABLE INTERPRO.SITE_MATCH_NEW
            ADD CONSTRAINT FK_SITE_MATCH_NEW$D
            FOREIGN KEY (DBCODE) REFERENCES INTERPRO.CV_DATABASE (DBCODE)
            """
        )

        logger.info(f"\texchanging partition")
        partition = site_partitions[database.identifier]
        cur.execute(
            f"""
            ALTER TABLE INTERPRO.SITE_MATCH
            EXCHANGE PARTITION ({partition})
            WITH TABLE INTERPRO.SITE_MATCH_NEW
            """
        )

        oracle.drop_table(cur, "INTERPRO.SITE_MATCH_NEW", purge=True)

    cur.close()
    con.close()

    logger.info("complete")


def rebuild_indexes(uri: str, table: str):
    con = oracledb.connect(uri)
    cur = con.cursor()

    for index in oracle.get_indexes(cur, "INTERPRO", table):
        if index["is_unusable"]:
            logger.info(f"rebuilding index {index['name']}")
            oracle.rebuild_index(cur, index["name"])

    for index in oracle.get_partitioned_indexes(cur, "INTERPRO", table):
        if index["is_unusable"]:
            logger.info(f"rebuilding index {index['name']}, "
                        f"partition {index['partition']}")
            oracle.rebuild_index(cur, index["name"],
                                 partition=index["partition"])

    cur.close()
    con.close()

    logger.info("complete")


def update_matches(uri: str):
    """
    Add protein matches for recently added/modified sequences

    :param uri: Oracle connection string
    """
    con = oracledb.connect(uri)
    _prepare_matches(con)
    _check_matches(con)
    _insert_matches(con)
    con.close()


def update_feature_matches(uri: str):
    """
    Add protein feature matches for recently added/modified sequences

    :param uri: Oracle connection string
    """
    con = oracledb.connect(uri)
    cur = con.cursor()
    logger.info("deleting matches")
    cur.execute(
        """
        DELETE FROM INTERPRO.FEATURE_MATCH
        WHERE PROTEIN_AC IN (
          SELECT PROTEIN_AC
          FROM INTERPRO.PROTEIN_TO_SCAN
        )
        """
    )
    con.commit()
    logger.info(f"  {cur.rowcount} rows deleted")

    for p in oracle.get_partitions(cur, "INTERPRO", "FEATURE_MATCH"):
        dbcode = p["value"][1:-1]  # 'X' -> X
        partition = p["name"]

        logger.info(f"inserting matches ({partition})")
        cur.execute(
            """
            INSERT INTO INTERPRO.FEATURE_MATCH
            SELECT P.PROTEIN_AC, M.METHOD_AC, M.SEQ_FEATURE, M.SEQ_START, 
                   M.SEQ_END, D.DBCODE
            FROM INTERPRO.PROTEIN_TO_SCAN P
            INNER JOIN IPRSCAN.MV_IPRSCAN M
              ON P.UPI = M.UPI
            INNER JOIN INTERPRO.IPRSCAN2DBCODE D
              ON M.ANALYSIS_ID = D.IPRSCAN_SIG_LIB_REL_ID
            WHERE D.DBCODE = :1
            """,
            [dbcode]
        )
        con.commit()
        logger.info(f"  {cur.rowcount} rows inserted")

    cur.close()
    con.close()
    logger.info("complete")


def update_variant_matches(uri: str):
    """
    Recreate splice-variants table with the most recent data
    from SWISSPROT_VARSPLIC. TREMBL_VARSPLIC (DBID=25) is obsolete and
    only contains deleted cross-references.

    :param uri: Oracle connection string
    """
    con = oracledb.connect(uri)
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
        WHERE X.DBID = 24
        AND X.DELETED = 'N'
        """
    )
    logger.info(f"{cur.rowcount} rows inserted")

    logger.info("updating VARSPLIC_MATCH")
    oracle.truncate_table(cur, "INTERPRO.VARSPLIC_MATCH", reuse_storage=True)

    dbcodes = []
    for p in oracle.get_partitions(cur, "INTERPRO", "MATCH"):
        dbcode = p["value"][1:-1]  # 'X' -> X
        dbcodes.append(dbcode)

    params = ",".join([f":{i+1}" for i in range(len(dbcodes))])
    cur.execute(
        f"""
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
        WHERE X.DBID = 24
          AND X.DELETED = 'N'
          AND I2D.DBCODE IN ({params})
          AND M.SKIP_FLAG = 'N'
        """,
        dbcodes
    )
    logger.info(f"{cur.rowcount} rows inserted")
    con.commit()

    cur.close()
    con.close()


def update_site_matches(uri: str):
    """
    Add protein site matches for recently added/modified sequences

    :param uri: Oracle connection string
    """
    con = oracledb.connect(uri)
    cur = con.cursor()

    logger.info("populating SITE_MATCH_NEW")
    oracle.drop_table(cur, "INTERPRO.SITE_MATCH_NEW", purge=True)

    cur.execute(
        """
        CREATE TABLE INTERPRO.SITE_MATCH_NEW NOLOGGING
        AS
        SELECT
            P.PROTEIN_AC, S.METHOD_AC, D.DBCODE, S.LOC_START, S.LOC_END, 
            S.DESCRIPTION, S.RESIDUE, S.RESIDUE_START, S.RESIDUE_END, 
            S.NUM_SITES
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
        TABLESPACE INTERPRO_IND
        NOLOGGING
        """
    )

    logger.info("gathering statistics")
    oracle.gather_stats(cur, "INTERPRO", "SITE_MATCH_NEW")

    logger.info("checking")
    match_partitions = {}
    for p in oracle.get_partitions(cur, "INTERPRO", "MATCH"):
        dbcode = p["value"][1:-1]  # 'X' -> X
        match_partitions[dbcode] = p["name"]

    queries = []
    params = []
    for p in oracle.get_partitions(cur, "INTERPRO", "SITE_MATCH"):
        dbcode = p["value"][1:-1]

        if dbcode in match_partitions:
            partition = match_partitions[dbcode]
            queries.append(
                f"""
                SELECT DISTINCT PROTEIN_AC, METHOD_AC, POS_FROM, POS_TO
                FROM INTERPRO.MATCH PARTITION ({partition})
                """
            )
            params.append(dbcode)

    if queries:
        in_cond = [f":{i+1}" for i in range(len(params))]
        cur.execute(
            f"""
            SELECT COUNT(*)
            FROM (
                SELECT DISTINCT PROTEIN_AC, METHOD_AC, LOC_START, LOC_END
                FROM INTERPRO.SITE_MATCH_NEW
                WHERE DBCODE IN ({', '.join(in_cond)})
                MINUS (
                  {' UNION ALL '.join(queries)}
                )
            )
            """,
            params
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


def _prepare_matches(con: oracledb.Connection):
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

    binds = []
    params = []
    for p in oracle.get_partitions(cur, "INTERPRO", "MATCH"):
        dbcode = p["value"][1:-1]  # 'X' -> X
        binds.append(":1")
        params.append(dbcode)

    cur.execute(
        f"""
        INSERT /*+ APPEND */ INTO INTERPRO.MATCH_NEW
        SELECT
          P.PROTEIN_AC, M.METHOD_AC, M.SEQ_START, M.SEQ_END, 'T',
          D.DBCODE, D.EVIDENCE,
          SYSDATE, SYSDATE, SYSDATE, 'INTERPRO',
          M.EVALUE, M.MODEL_AC, M.FRAGMENTS, 
          CASE WHEN D.DBCODE = 'V' THEN M.SEQ_FEATURE ELSE NULL END
        FROM INTERPRO.PROTEIN_TO_SCAN P
        INNER JOIN IPRSCAN.MV_IPRSCAN M
          ON P.UPI = M.UPI
        INNER JOIN INTERPRO.IPRSCAN2DBCODE D
          ON M.ANALYSIS_ID = D.IPRSCAN_SIG_LIB_REL_ID
        WHERE D.DBCODE IN ({','.join(binds)})
        AND M.SEQ_START != M.SEQ_END
        """,
        params
    )
    con.commit()

    logger.info("indexing")
    for col in ("DBCODE", "PROTEIN_AC"):
        cur.execute(
            f"""
            CREATE INDEX I_MATCH_NEW${col}
            ON INTERPRO.MATCH_NEW ({col})
            TABLESPACE INTERPRO_IND
            NOLOGGING
            """
        )

    # logger.info("gathering statistics")
    # oracle.gather_stats(cur, "INTERPRO", "MATCH_NEW")

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


def _check_matches(con: oracledb.Connection):
    """
    Check there are not errors in imported matches

    :param con: Oracle connection object
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

    cur.close()


def _insert_matches(con: oracledb.Connection):
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
    con.commit()

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


def track_entry_changes(
    cur: oracledb.Cursor,
    pg_uri: str,
    data_dir: str,
    threshold: float
) -> list:
    """
    Find entries with significant protein count changes

    :param cur: Oracle cursor object
    :param pg_uri: PostgreSQL connection string
    :param data_dir: directory containing the file for protein counts
                     before the update
    :param threshold: report entries with at least (threshold/100)% change
                      in protein counts
    :return: list of entries with a significant changes in proteins count
    """

    with open(os.path.join(data_dir, FILE_ENTRY_PROT_COUNTS), "rb") as fh:
        old_counts = pickle.load(fh)

    new_counts = _get_entries_protein_counts(cur, pg_uri)

    changes = []
    for acc in sorted(old_counts):
        entry_old_counts = old_counts[acc]
        entry_new_counts = new_counts.pop(acc, {"total": {}, "swissprot": 0, "pdb": 0})

        # Total number of proteins matched
        entry_old_total = sum(entry_old_counts["total"].values())
        entry_new_total = sum(entry_new_counts["total"].values())
        change_total = (entry_new_total - entry_old_total) / entry_old_total if entry_old_total else 0

        # Total number of PDB-mapped proteins matched
        entry_old_pdb = entry_old_counts["pdb"]
        entry_new_pdb = entry_new_counts["pdb"]
        change_pdb = (entry_new_pdb - entry_old_pdb) / entry_old_pdb if entry_old_pdb else 0

        # Total number of swissprot proteins
        entry_old_swiss = entry_old_counts["swissprot"]
        entry_new_swiss = entry_old_counts["swissprot"]
        change_swiss = (entry_new_swiss - entry_old_swiss) / entry_old_swiss if entry_old_swiss else 0

        # If the entry does not have any matches anymore,
        # we want to report it
        if entry_new_total != 0 and abs(change_total) < threshold:
            continue

        # If there were not matches before and no matches now
        # we don't want to report it
        if not entry_old_total and not entry_new_total:
            continue

        entry_superkingdoms = {}
        for superkingdom, old_cnt in entry_old_counts["total"].items():
            new_cnt = entry_new_counts["total"].pop(superkingdom, 0)
            entry_superkingdoms[superkingdom] = (old_cnt, new_cnt)

        # superkingdoms with proteins only matched in new UniProt release
        for superkingdom, new_cnt in entry_new_counts["total"].items():
            entry_superkingdoms[superkingdom] = (0, new_cnt)

        changes.append((
            acc,
            entry_old_total,
            entry_new_total,
            change_total,
            entry_old_pdb,
            entry_new_pdb,
            change_pdb,
            entry_old_swiss,
            entry_new_swiss,
            change_swiss,
            entry_superkingdoms
        ))

    return changes


def _get_taxon2superkingdom(cur: oracledb.Cursor) -> dict[int, str]:
    # Load all taxa
    cur.execute(
        """
        SELECT TAX_ID, SCIENTIFIC_NAME, RANK, PARENT_ID
        FROM INTERPRO.ETAXI
        """
    )
    taxa = {}
    for tax_id, name, rank, parent_id in cur:
        if tax_id in (1, 131567):
            """
            Skip root and cellular organisms (131567) which contains:
                * Bacteria (2)
                * Archaea (2157)
                * Eukaryota (2759)
            """
            continue
        elif parent_id in (1, 131567):
            rank = "domain"
            parent_id = None

        taxa[tax_id] = (name, rank, parent_id)

    # For each taxon, find its root (i.e. superkingdom)
    taxon2superkingdom = {}
    for tax_id in taxa:
        name, rank, parent_id = taxa[tax_id]

        while parent_id is not None:
            name, rank, parent_id = taxa[parent_id]
            if rank == "domain":
                taxon2superkingdom[tax_id] = name
                break

    return taxon2superkingdom


def _get_entries_protein_counts(
        cur: oracledb.Cursor,
        pg_url: str
) -> dict[str, dict[str, dict[str, int]]]:
    """
    Return a dict with the 'total' number of proteins matched
    by each InterPro entry per superkingdom, the number of
    `swissprot` proteins that match the entry, and the number of
    proteins with at least one 'pdb' structure.

    Only complete sequences are considered.

    :param cur: Oracle cursor object
    :param pg_url: PostgreSQL connection string
    :return: dictionary
    """
    taxon2superkingdom = _get_taxon2superkingdom(cur)

    # Get total protein counts per entry per superkingdom
    cur.execute(
        """
        SELECT EM.ENTRY_AC, P.TAX_ID, COUNT(DISTINCT P.PROTEIN_AC)
        FROM INTERPRO.PROTEIN P
        INNER JOIN INTERPRO.MATCH M ON P.PROTEIN_AC = M.PROTEIN_AC
        INNER JOIN INTERPRO.ENTRY2METHOD EM ON EM.METHOD_AC = M.METHOD_AC
        WHERE P.FRAGMENT = 'N'
        GROUP BY EM.ENTRY_AC, P.TAX_ID
        """
    )
    counts = {}
    for entry_acc, tax_id, n_proteins in cur:
        try:
            e = counts[entry_acc]
        except KeyError:
            e = counts[entry_acc] = {"total": {}, "swissprot": 0, "pdb": 0}
        superkingdom = taxon2superkingdom[tax_id]
        try:
            e["total"][superkingdom] += n_proteins
        except KeyError:
            e["total"][superkingdom] = n_proteins

    # Get number of assoiated UniProt entries with at least one PDB
    integrated = defaultdict(list)
    cur.execute("SELECT ENTRY_AC, METHOD_AC FROM INTERPRO.ENTRY2METHOD")
    for entry_acc, method_acc in cur:
        integrated[entry_acc].append(method_acc)

    pg_con = psycopg.connect(**url2dict(pg_url))
    pg_cur = pg_con.cursor()

    for entry_acc in integrated:
        try:
            e = counts[entry_acc]
        except KeyError:
            e = counts[entry_acc] = {"total": {}, "swissprot": 0, "pdb": 0}
        pg_cur.execute(
            """
            SELECT COUNT(DISTINCT S.PROTEIN_ACC)
            FROM SIGNATURE2STRUCTURE S
            WHERE S.SIGNATURE_ACC = ANY(%s)
            """,
            (integrated[entry_acc],)
        )
        pdb_count, = pg_cur.fetchone()
        e["pdb"] = pdb_count

    # Get the number of specifically Swissprot proteins
    for entry_acc in integrated:
        pg_cur.execute(
            """
            SELECT COUNT(DISTINCT SP.PROTEIN_ACC)
            FROM SIGNATURE2PROTEIN SP
            WHERE SP.SIGNATURE_ACC = ANY(%s) AND IS_REVIEWED
            """,
            (list(integrated[entry_acc]),)
        )
        pdb_count, = pg_cur.fetchone()
        e = counts[entry_acc]
        e["swissprot"] = pdb_count

    pg_cur.close()
    pg_con.close()

    return counts


def get_sig_protein_counts(cur: oracledb.Cursor, pg_url: str,
                           dbcode: str) -> dict[str, dict[str, int]]:
    """
    Return the number of protein matches and the number of proteins with
    at least one 'pdb' structure by each member database signature.
    Only complete sequences are considered

    :param cur: Oracle cursor object
    :param pg_url: PostgreSQL connection string
    :param dbcode: member database code
    :return: dictionary
    """
    taxon2superkingdom = _get_taxon2superkingdom(cur)
    cur.execute(
        f"""
        SELECT M.METHOD_AC, P.TAX_ID, COUNT(DISTINCT P.PROTEIN_AC)
        FROM INTERPRO.MATCH M
        INNER JOIN INTERPRO.PROTEIN P
            ON P.PROTEIN_AC = M.PROTEIN_AC
        WHERE M.DBCODE = :1 
          AND P.FRAGMENT = 'N'
        GROUP BY M.METHOD_AC, P.TAX_ID
        """,
        [dbcode]
    )
    counts = {}
    for sig_acc, tax_id, n_proteins in cur:
        try:
            sig = counts[sig_acc]
        except KeyError:
            sig = counts[sig_acc] = {"total": {}, "pdb": 0}

        superkingdom = taxon2superkingdom[tax_id]
        try:
            sig["total"][superkingdom] += n_proteins
        except KeyError:
            sig["total"][superkingdom] = n_proteins

    pg_con = psycopg.connect(**url2dict(pg_url))
    pg_cur = pg_con.cursor()

    signatures = list(counts.keys())
    if signatures:
        pg_cur.execute(
            """
            SELECT SIGNATURE_ACC, COUNT(DISTINCT PROTEIN_ACC)
            FROM SIGNATURE2STRUCTURE
            WHERE SIGNATURE_ACC = ANY(%s)
            GROUP BY SIGNATURE_ACC
            """,
            (signatures,)
        )

    for sig_acc, pdb_count in pg_cur:
        try:
            sig = counts[sig_acc]
        except KeyError:
            sig = counts[sig_acc] = {"total": {}, "pdb": 0}

        sig["pdb"] = pdb_count

    pg_cur.close()
    pg_con.close()

    return counts


def insert_toad_matches(uri: str,
                        databases: list[Database],
                        files: dict[str, str],
                        tmpdir: str | None = None):
    """
    Update the TOAD matches table with matches extracted from TAR archives
    :param uri: Oracle connection string
    :param databases: list of Database objects of member databases to update
    :param files: Dictionary of database name -> tar file path
    :param tmpdir: Path to directory for temporary files
    """
    _databases = {}
    for db in databases:
        _databases[db.identifier] = files[db.identifier]

    toad.load_matches(uri, _databases, tmpdir=tmpdir)
    rebuild_indexes(uri, "TOAD_MATCH")
    logger.info("done")


def update_toad_matches(uri: str):
    """
    Delete TOAD matches for proteins recently changed
    :param uri: Oracle connection string
    """
    con = oracledb.connect(uri)
    cur = con.cursor()
    logger.info("deleting matches")
    cur.execute(
        """
        DELETE FROM INTERPRO.TOAD_MATCH
        WHERE PROTEIN_AC IN (
          SELECT PROTEIN_AC
          FROM INTERPRO.PROTEIN_TO_SCAN
        )
        """
    )
    con.commit()
    logger.info(f"  {cur.rowcount} rows deleted")

    cur.close()
    con.close()
    logger.info("done")
