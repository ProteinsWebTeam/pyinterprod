# -*- coding: utf-8 -*-

import cx_Oracle

from pyinterprod import logger
from pyinterprod.utils import oracle, Table


def refresh_taxonomy(ipr_url: str, swp_url: str):
    ipr_con = cx_Oracle.connect(ipr_url)

    logger.info("creating TAXONOMY_LOAD")
    ipr_cur = ipr_con.cursor()
    oracle.drop_table(ipr_cur, "INTERPRO.TAXONOMY_LOAD", purge=True)
    ipr_cur.execute(
        """
        CREATE TABLE INTERPRO.TAXONOMY_LOAD
        (
            TAX_ID NUMBER(10) NOT NULL,
            PARENT_ID NUMBER(10),
            SCIENTIFIC_NAME VARCHAR2(255),
            RANK VARCHAR2(50) NOT NULL,
            COMMON_NAME VARCHAR2(255)
        ) NOLOGGING          
        """
    )
    ipr_cur.close()

    req = """
        INSERT /*+ APPEND */ 
        INTO INTERPRO.TAXONOMY_LOAD 
        VALUES (:1, :2, :3, :4, :5)
    """
    with Table(ipr_con, req, autocommit=True) as taxonomy_load:
        swp_con = cx_Oracle.connect(swp_url)
        swp_cur = swp_con.cursor()
        swp_cur.execute(
            """
            SELECT 
                TAX_ID, PARENT_ID, SPTR_SCIENTIFIC, RANK, 
                NVL(SPTR_COMMON, NCBI_COMMON)
            FROM TAXONOMY.V_PUBLIC_NODE
            """
        )

        for rec in swp_cur:
            taxonomy_load.insert(rec)

        swp_cur.close()
        swp_con.close()

    ipr_cur = ipr_con.cursor()
    oracle.gather_stats(ipr_cur, "INTERPRO", "TAXONOMY_LOAD")

    logger.info("populating ETAXI")
    oracle.truncate_table(ipr_cur, "INTERPRO.ETAXI", reuse_storage=True)
    ipr_cur.execute(
        """
        INSERT INTO INTERPRO.ETAXI
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
    ipr_cur.commit()
    oracle.gather_stats(ipr_cur, "INTERPRO", "ETAXI")

    # Dropping temporary table
    oracle.drop_table(ipr_cur, "INTERPRO.TAXONOMY_LOAD", purge=True)

    ipr_cur.close()
    ipr_con.close()
    logger.info("complete")
