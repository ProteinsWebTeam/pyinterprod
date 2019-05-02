# -*- coding: utf-8 -*-

import cx_Oracle

from .. import logger, orautils


def update_taxonomy(user: str, dsn: str):
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()

    logger.info("creating TAXONOMY_LOAD")
    orautils.drop_table(cur, "INTERPRO", "TAXONOMY_LOAD")
    cur.execute(
        """
        CREATE TABLE INTERPRO.TAXONOMY_LOAD
        NOLOGGING
        AS
        SELECT 
          TAX_ID, 
          PARENT_ID, 
          SPTR_SCIENTIFIC AS SCIENTIFIC_NAME,
          RANK,
          NVL(N.SPTR_COMMON, N.NCBI_COMMON) AS COMMON_NAME
        FROM TAXONOMY.V_PUBLIC_NODE@SWPREAD 
        """
    )
    orautils.gather_stats(cur, "INTERPRO", "TAXONOMY_LOAD")

    logger.info("creating ETAXI")
    orautils.drop_table(cur, "INTERPRO", "ETAXI")
    cur.execute(
        """
        CREATE TABLE INTERPRO.ETAXI 
        AS
          SELECT 
            N.TAX_ID, N.PARENT_ID, N.SCIENTIFIC_NAME, 
            'X' COMPLETE_GENOME_FLAG, N.RANK, 0 HIDDEN, 
            LR.TREE_LEFT LEFT_NUMBER, LR.TREE_RIGHT RIGHT_NUMBER, 
            'X' ANNOTATION_SOURCE, 
            N.SCIENTIFIC_NAME || CASE WHEN N.COMMON_NAME IS NULL THEN '' ELSE ' (' || N.COMMON_NAME || ')' END FULL_NAME
        FROM INTERPRO.TAXONOMY_LOAD N
        INNER JOIN (
            SELECT 
              TAX_ID, 
              MIN(TREE_NUMBER) TREE_LEFT, 
              MAX(TREE_NUMBER) TREE_RIGHT
            FROM (
                SELECT PARENT_ID AS TAX_ID, ROWNUM AS TREE_NUMBER
                FROM (
                    SELECT TAX_ID, PARENT_ID
                    FROM (
                        SELECT TAX_ID, PARENT_ID 
                        FROM INTERPRO.TAXONOMY_LOAD
                        UNION ALL
                        SELECT 9999999 AS TAX_ID, TAX_ID AS PARENT_ID 
                        FROM INTERPRO.TAXONOMY_LOAD
                        UNION ALL
                        SELECT 0 AS TAX_ID, TAX_ID AS PARENT_ID 
                        FROM INTERPRO.TAXONOMY_LOAD
                    )
                    START WITH TAX_ID = 1
                    CONNECT BY PRIOR TAX_ID=PARENT_ID
                    ORDER SIBLINGS BY TAX_ID
                )
                WHERE TAX_ID IN (9999999, 0)
            )
            GROUP BY TAX_ID
        ) LR
        ON LR.TAX_ID = N.TAX_ID
        """
    )
    orautils.gather_stats(cur, "INTERPRO", "ETAXI")
    cur.execute(
        """
        CREATE INDEX I_ETAXI$L$R$T 
        ON INTERPRO.ETAXI (LEFT_NUMBER, RIGHT_NUMBER, TAX_ID)
        """
    )
    cur.execute(
        """
        CREATE INDEX I_ETAXI$P$T$R
        ON INTERPRO.ETAXI (PARENT_ID, TAX_ID, RANK)
        """
    )
    cur.execute(
        """
        CREATE INDEX I_ETAXI$T$P$R
        ON INTERPRO.ETAXI (TAX_ID, PARENT_ID, RANK)
        """
    )

    # Dropping temporary table
    orautils.drop_table(cur, "INTERPRO", "TAXONOMY_LOAD")

    orautils.grant(cur, "INTERPRO", "ETAXI", "SELECT", "PUBLIC")
    cur.close()
    con.close()
