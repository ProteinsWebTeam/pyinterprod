import os
import pickle

import oracledb

from pyinterprod import logger
from pyinterprod.utils import oracle


FILE_ENTRY_PROT_COUNTS = "entries.prot.counts.pickle"

"""
To add a partition:
SQL> ALTER TABLE <TABLE> 
     ADD PARTITION <NAME> VALUES ('<DBCODE>');
"""
MATCH_PARTITIONS = {
    "B": "MATCH_DBCODE_B",  # SFLD
    "F": "MATCH_DBCODE_F",  # PRINTS
    "H": "MATCH_DBCODE_H",  # Pfam
    "J": "MATCH_DBCODE_J",  # CDD
    "M": "MATCH_DBCODE_M",  # PROSITE profiles
    "N": "MATCH_DBCODE_N",  # NCBIfam
    "P": "MATCH_DBCODE_P",  # PROSITE patterns
    "Q": "MATCH_DBCODE_Q",  # HAMAP
    "R": "MATCH_DBCODE_R",  # SMART
    "U": "MATCH_DBCODE_U",  # PIRSF
    "V": "MATCH_DBCODE_V",  # PANTHER
    "X": "MATCH_DBCODE_X",  # CATH-Gene3D
    "Y": "MATCH_DBCODE_Y",  # SUPERFAMILY
}
FEATURE_MATCH_PARTITIONS = {
    "a": "ANTIFAM",
    "d": "PFAM_N",
    "f": "FUNFAM",
    "g": "MOBIDBLITE",
    "j": "PHOBIUS",
    "l": "ELM",
    "n": "SIGNALP_E",
    "q": "TMHMM",
    "s": "SIGNALP_GP",
    "v": "SIGNALP_GN",
    "x": "COILS"
}
SITE_PARTITIONS = {
    # DB identifier -> (str:partition, bool:check against MATCH table)
    "B": ("SFLD", True),
    "J": ("CDD", True),
    "z": ("PIRSR", False)
}


def export_entries_protein_counts(cur: oracledb.Cursor, data_dir: str):
    with open(os.path.join(data_dir, FILE_ENTRY_PROT_COUNTS), "wb") as fh:
        pickle.dump(_get_entries_protein_counts(cur), fh)


def update_database_matches(uri: str, databases: list):
    """

    :param uri:
    :param databases: list of Database objects
    :return:
    """
    con = oracledb.connect(uri)
    cur = con.cursor()

    for database in databases:
        logger.info(f"{database.name}")
        oracle.drop_table(cur, "INTERPRO.MATCH_NEW", purge=True)
        logger.info(f"\tpopulating MATCH_MEW "
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
        partition = MATCH_PARTITIONS[database.identifier]
        cur.execute(
            f"""
            ALTER TABLE INTERPRO.MATCH
            EXCHANGE PARTITION ({partition})
            WITH TABLE INTERPRO.MATCH_NEW
            """
        )
        oracle.drop_table(cur, "INTERPRO.MATCH_NEW", purge=True)

        logger.info("\tgathering statistics")
        oracle.gather_stats(cur, "INTERPRO", "MATCH", partition)

    for index in oracle.get_indexes(cur, "INTERPRO", "MATCH"):
        if index["is_unusable"]:
            logger.info(f"rebuilding index {index['name']}")
            oracle.rebuild_index(cur, index["name"])

    for index in oracle.get_partitioned_indexes(cur, "INTERPRO", "MATCH"):
        if index["is_unusable"]:
            logger.info(f"rebuilding index {index['name']}, "
                        f"partition {index['partition']}")
            oracle.rebuild_index(cur, index["name"],
                                 partition=index["partition"])

    cur.close()
    con.close()

    logger.info("complete")


def update_database_feature_matches(uri: str, databases: list):
    """

    :param uri:
    :param databases: list of Database objects
    :return:
    """
    con = oracledb.connect(uri)
    cur = con.cursor()

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
        partition = FEATURE_MATCH_PARTITIONS[database.identifier]
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

    for index in oracle.get_indexes(cur, "INTERPRO", "FEATURE_MATCH"):
        if index["is_unusable"]:
            logger.info(f"rebuilding index {index['name']}")
            oracle.rebuild_index(cur, index["name"])

    for index in oracle.get_partitioned_indexes(cur, "INTERPRO",
                                                "FEATURE_MATCH"):
        if index["is_unusable"]:
            logger.info(f"rebuilding index {index['name']}, "
                        f"partition {index['partition']}")
            oracle.rebuild_index(cur, index["name"],
                                 partition=index["partition"])

    cur.close()
    con.close()

    logger.info("complete")


def update_database_site_matches(uri: str, databases: list):
    """

    :param uri:
    :param databases: list of Database objects
    :return:
    """
    con = oracledb.connect(uri)
    cur = con.cursor()

    for database in databases:
        logger.info(database.name)

        """
        Same partition names in:
            - IPRSCAN.SITE
            - INTERPRO.SITE_MATCH
        """
        site_partition, ck_matches = SITE_PARTITIONS[database.identifier]

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
            FROM IPRSCAN.SITE PARTITION ({site_partition}) S
            INNER JOIN UNIPARC.XREF X
              ON S.UPI = X.UPI
            INNER JOIN INTERPRO.IPRSCAN2DBCODE D
              ON S.ANALYSIS_ID = D.IPRSCAN_SIG_LIB_REL_ID
            WHERE S.ANALYSIS_ID = :1
            AND X.DBID IN (2, 3)  -- Swiss-Prot or TrEMBL
            AND X.DELETED = 'N'
            """, (database.analysis_id,)
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

        if ck_matches:
            logger.info(f"\tchecking matches")
            match_partition = MATCH_PARTITIONS[database.identifier]
            cur.execute(
                f"""
                SELECT COUNT(*)
                FROM (
                    SELECT DISTINCT PROTEIN_AC, METHOD_AC, LOC_START, LOC_END
                    FROM INTERPRO.SITE_MATCH_NEW
                    MINUS
                    SELECT DISTINCT PROTEIN_AC, METHOD_AC, POS_FROM, POS_TO
                    FROM INTERPRO.MATCH PARTITION ({match_partition})
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
        cur.execute(
            f"""
            ALTER TABLE INTERPRO.SITE_MATCH
            EXCHANGE PARTITION ({site_partition})
            WITH TABLE INTERPRO.SITE_MATCH_NEW
            """
        )

        oracle.drop_table(cur, "INTERPRO.SITE_MATCH_NEW", purge=True)

    for index in oracle.get_indexes(cur, "INTERPRO", "SITE_MATCH"):
        if index["is_unusable"]:
            logger.info(f"rebuilding index {index['name']}")
            oracle.rebuild_index(cur, index["name"])

    for index in oracle.get_partitioned_indexes(cur, "INTERPRO", "SITE_MATCH"):
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

    for dbcode, partition in FEATURE_MATCH_PARTITIONS.items():
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
        logger.info(f"  {cur.rowcount} rows inserted")

    con.commit()
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

    dbcodes = list(MATCH_PARTITIONS.keys())
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
    params = []
    queries = []
    for identifier, (_, ck_matches) in SITE_PARTITIONS.items():
        if ck_matches:
            partition = MATCH_PARTITIONS[identifier]
            params.append(identifier)
            queries.append(
                f"""
                SELECT DISTINCT PROTEIN_AC, METHOD_AC, POS_FROM, POS_TO
                FROM INTERPRO.MATCH PARTITION ({partition})
                """
            )

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
            """, params
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

    params = ",".join([":" + str(i + 1)
                       for i in range(len(MATCH_PARTITIONS))])
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
        WHERE D.DBCODE IN ({params})
        AND M.SEQ_START != M.SEQ_END
        """,
        list(MATCH_PARTITIONS.keys())
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


def track_entry_changes(cur: oracledb.Cursor, data_dir: str,
                        threshold: float) -> list:
    """
    Find entries with significant protein count changes

    :param cur: Oracle cursor object
    :param data_dir: directory containing the file for protein counts
                     before the update
    :param threshold: report entries with at least (threshold/100)% change
                      in protein counts
    :return: list of entries with a significant changes in proteins count
    """

    with open(os.path.join(data_dir, FILE_ENTRY_PROT_COUNTS), "rb") as fh:
        old_counts = pickle.load(fh)

    new_counts = _get_entries_protein_counts(cur)
    changes = []
    for acc in sorted(old_counts):
        entry_old_counts = old_counts[acc]
        entry_new_counts = new_counts.pop(acc, {})

        # Total number of proteins matched
        entry_old_total = sum(entry_old_counts.values())
        entry_new_total = sum(entry_new_counts.values())

        change = (entry_new_total - entry_old_total) / entry_old_total

        # If the entry does not have any matches anymore,
        # we want to report it
        if entry_new_total != 0 and abs(change) < threshold:
            continue

        entry_superkingdoms = {}
        for superkingdom, old_cnt in entry_old_counts.items():
            new_cnt = entry_new_counts.pop(superkingdom, 0)
            entry_superkingdoms[superkingdom] = (old_cnt, new_cnt)

        # superkingdoms with proteins only matched in new UniProt release
        for superkingdom, new_cnt in entry_new_counts.items():
            entry_superkingdoms[superkingdom] = (0, new_cnt)

        changes.append((
            acc,
            entry_old_total,
            entry_new_total,
            change,
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
            Skip root and meta-superkingdom (131567) which contains:
                * Bacteria (2)
                * Archaea (2157)
                * Eukaryota (2759)
            """
            continue
        elif parent_id in (1, 131567):
            rank = "superkingdom"
            parent_id = None

        taxa[tax_id] = (name, rank, parent_id)

    # For each taxon, find its root (i.e. superkingdom)
    taxon2superkingdom = {}
    for tax_id in taxa:
        name, rank, parent_id = taxa[tax_id]

        while parent_id is not None:
            name, rank, parent_id = taxa[parent_id]
            if rank == "superkingdom":
                taxon2superkingdom[tax_id] = name
                break

    return taxon2superkingdom


def _get_entries_protein_counts(
        cur: oracledb.Cursor
) -> dict[str, dict[str, int]]:
    """
    Return the number of protein matched by each InterPro entry.
    Only complete sequences are considered.

    :param cur: Oracle cursor object
    :return: dictionary
    """
    taxon2superkingdom = _get_taxon2superkingdom(cur)
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
            e = counts[entry_acc] = {}

        superkingdom = taxon2superkingdom[tax_id]
        try:
            e[superkingdom] += n_proteins
        except KeyError:
            e[superkingdom] = n_proteins

    return counts


def get_sig_protein_counts(cur: oracledb.Cursor,
                           dbid: str) -> dict[str, dict[str, int]]:
    """
    Return the number of protein matches by each member database signature.
    Only complete sequences are considered

    :param cur: Oracle cursor object
    :param dbid: member database identifier
    :return: dictionary
    """
    partition = MATCH_PARTITIONS[dbid]
    taxon2superkingdom = _get_taxon2superkingdom(cur)
    cur.execute(
        f"""
        SELECT M.METHOD_AC, P.TAX_ID, COUNT(DISTINCT P.PROTEIN_AC)
        FROM INTERPRO.MATCH PARTITION ({partition}) M
        INNER JOIN INTERPRO.PROTEIN P
            ON P.PROTEIN_AC = M.PROTEIN_AC
        WHERE P.FRAGMENT = 'N'
        GROUP BY M.METHOD_AC, P.TAX_ID
        """
    )
    counts = {}
    for sig_acc, tax_id, n_proteins in cur:
        try:
            sig = counts[sig_acc]
        except KeyError:
            sig = counts[sig_acc] = {}

        superkingdom = taxon2superkingdom[tax_id]
        try:
            sig[superkingdom] += n_proteins
        except KeyError:
            sig[superkingdom] = n_proteins

    return counts


# def _get_databases_matches_count(cur: oracledb.Cursor) -> dict[str, int]:
#     """
#     Return the number of matches per member database
#     :param cur: Oracle cursor object
#     :return: dictionary
#     """
#     cur.execute(
#         """
#         SELECT DBCODE, COUNT(*)
#         FROM INTERPRO.MATCH
#         GROUP BY DBCODE
#         """
#     )
#     return dict(cur.fetchall())


# def add_site_subpartitions(uri: str, owner: str, table: str, partition: str,
#                            stop: str, prefix: str = ""):
#     if len(stop) != 8 or stop[:3] != "UPI" or not stop[3:].isalnum():
#         raise ValueError(f"Invalid range stop: {stop}. "
#                          f"Expected format: UPIxxxxx, with x being digits")
#
#     con = oracledb.connect(uri)
#     cur = con.cursor()
#
#     subpartitions = set()
#     for subpart in oracle.get_subpartitions(cur, owner, table, partition):
#         subpartitions.add(subpart["name"])
#
#     new_subpartitions = {}
#     # range_upi yields start/stop, but since step = 1, start == stop
#     for value, _ in range_upi("UPI00000", stop, 1):
#         name = prefix + value
#         if name not in subpartitions:
#             new_subpartitions[name] = value
#
#     for name, value in new_subpartitions.items():
#         cur.execute(
#             f"""
#             ALTER TABLE {table} MODIFY PARTITION {partition}
#             ADD SUBPARTITION {name} VALUES ('{value}')
#             """
#         )
#
#     cur.close()
#     con.close()
