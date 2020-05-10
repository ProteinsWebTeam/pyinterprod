# -*- coding: utf-8 -*-

import cx_Oracle

from pyinterprod import logger
from pyinterprod.utils import oracle


def refresh_taxonomy(url: str):
    con = cx_Oracle.connect(url)
    cur = con.cursor()

    logger.info("creating TAXONOMY_LOAD")
    oracle.drop_table(cur, "INTERPRO.TAXONOMY_LOAD", purge=True)
    cur.execute(
        """
        CREATE TABLE INTERPRO.TAXONOMY_LOAD NOLOGGING 
        AS SELECT TAX_ID, PARENT_ID, SPTR_SCIENTIFIC AS SCIENTIFIC_NAME, RANK, 
                  NVL(SPTR_COMMON, NCBI_COMMON) AS COMMON_NAME
        FROM TAXONOMY.V_PUBLIC_NODE@SWPREAD
        """
    )
    oracle.gather_stats(cur, "INTERPRO", "TAXONOMY_LOAD")

    logger.info("populating ETAXI")
    oracle.truncate_table(cur, "INTERPRO.ETAXI", reuse_storage=True)
    cur.execute(
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
    con.commit()
    oracle.gather_stats(cur, "INTERPRO", "ETAXI")

    # Dropping temporary table
    oracle.drop_table(cur, "INTERPRO.TAXONOMY_LOAD", purge=True)

    # TODO: stop refreshing the tables below when we stop supporting InterPro6

    logger.info("populating UNIPROT_TAXONOMY")
    oracle.truncate_table(cur, "INTERPRO.UNIPROT_TAXONOMY", reuse_storage=True)
    cur.execute(
        """
        INSERT INTO INTERPRO.UNIPROT_TAXONOMY
        SELECT P.PROTEIN_AC, P.TAX_ID, NVL(ET.LEFT_NUMBER, 0) LEFT_NUMBER, 
               NVL(ET.RIGHT_NUMBER, 0) RIGHT_NUMBER
        FROM INTERPRO.PROTEIN P
        LEFT OUTER JOIN INTERPRO.ETAXI ET ON P.TAX_ID = ET.TAX_ID
        """
    )
    con.commit()
    oracle.gather_stats(cur, "INTERPRO", "UNIPROT_TAXONOMY")

    # logger.info("populating MV_TAX_ENTRY_COUNT")
    # oracle.truncate_table(cur, "INTERPRO.MV_TAX_ENTRY_COUNT", reuse_storage=True)
    # cur.execute(
    #     """
    #     INSERT INTO  INTERPRO.MV_TAX_ENTRY_COUNT
    #     WITH QUERY1 AS (
    #       SELECT ENTRY_AC, ANC.PARENT AS TAX_ID, COUNT(1) AS COUNT
    #       FROM INTERPRO.UNIPROT_TAXONOMY UT
    #       JOIN INTERPRO.MV_ENTRY2PROTEIN_TRUE MVEP
    #         ON UT.PROTEIN_AC=MVEP.PROTEIN_AC
    #       JOIN (
    #         SELECT NVL(
    #                  SUBSTR(
    #                    SYS_CONNECT_BY_PATH(TAX_ID, '.'),
    #                    2,
    #                    INSTR(SYS_CONNECT_BY_PATH (TAX_ID,'.'),'.',2) - 2
    #                  ),
    #                  TAX_ID
    #                ) AS CHILD,
    #                TAX_ID AS PARENT
    #         FROM INTERPRO.ETAXI ET
    #         CONNECT BY PRIOR PARENT_ID=TAX_ID
    #       ) ANC ON ANC.CHILD=UT.TAX_ID
    #       GROUP BY ENTRY_AC, ANC.PARENT
    #     ),
    #     QUERY2 AS (
    #       SELECT ENTRY_AC, TAX_ID, COUNT(1) AS COUNT
    #       FROM INTERPRO.UNIPROT_TAXONOMY UT
    #       JOIN INTERPRO.MV_ENTRY2PROTEIN_TRUE MVEP
    #       ON UT.PROTEIN_AC=MVEP.PROTEIN_AC
    #       GROUP BY ENTRY_AC, TAX_ID
    #     )
    #     SELECT QUERY1.ENTRY_AC, QUERY1.TAX_ID, QUERY1.COUNT AS COUNT,
    #            QUERY2.COUNT AS COUNT_SPECIFIED_TAX_ID
    #     FROM QUERY1
    #     LEFT OUTER JOIN QUERY2
    #       ON QUERY1.ENTRY_AC = QUERY2.ENTRY_AC
    #       AND QUERY1.TAX_ID = QUERY2.TAX_ID
    #     """
    # )
    # con.commit()
    # oracle.gather_stats(cur, "INTERPRO", "MV_TAX_ENTRY_COUNT")

    cur.close()
    con.close()
    logger.info("complete")
