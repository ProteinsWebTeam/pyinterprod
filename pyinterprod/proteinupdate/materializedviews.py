import cx_Oracle

from .. import logger, orautils


def refresh(url: str):
    con = cx_Oracle.connect(url)
    cur = con.cursor()

    logger.info("creating MV_METHOD2PROTEIN")
    orautils.drop_table(cur, "INTERPRO", "MV_METHOD2PROTEIN")
    cur.execute(
        """
        CREATE TABLE INTERPRO.MV_METHOD2PROTEIN (
            METHOD_AC NOT NULL VARCHAR2(25),
            PROTEIN_AC NOT NULL VARCHAR2(15),
            MATCH_COUNT NOT NULL NUMBER(7)
        ) NOLOGGING
        """
    )
    cur.execute(
        """
        INSERT /*+ APPEND */ INTO INTERPRO.MV_METHOD2PROTEIN
        SELECT METHOD_AC, PROTEIN_AC, COUNT(*)
        FROM INTERPRO.MATCH
        GROUP BY METHOD_AC, PROTEIN_AC
        """
    )
    con.commit()

    logger.info("indexing MV_METHOD2PROTEIN")
    cur.execute(
        """
        CREATE INDEX I_MV_METHOD2PROTEIN
        ON INTERPRO.MV_METHOD2PROTEIN (METHOD_AC)
        """
    )

    logger.info("creating MV_METHOD_MATCH")
    orautils.drop_table(cur, "INTERPRO", "MV_METHOD_MATCH")
    cur.execute(
        """
        CREATE TABLE INTERPRO.MV_METHOD_MATCH (
            METHOD_AC NOT NULL VARCHAR2(25),
            PROTEIN_COUNT NOT NULL NUMBER(8),
            MATCH_COUNT NOT NULL NUMBER(8)
        ) NOLOGGING
        """
    )
    cur.execute(
        """
        INSERT /*+ APPEND */ INTO INTERPRO.MV_METHOD_MATCH
        SELECT METHOD_AC, COUNT(*), SUM(MATCH_COUNT)
        FROM INTERPRO.MV_METHOD2PROTEIN
        GROUP BY METHOD_AC
        """
    )
    con.commit()

    logger.info("indexing MV_METHOD_MATCH")
    cur.execute(
        """
        CREATE UNIQUE INDEX UI_MV_METHOD_MATCH
        ON INTERPRO.MV_METHOD_MATCH (METHOD_AC)
        """
    )

    logger.info("creating MV_ENTRY2PROTEIN")
    orautils.drop_table(cur, "INTERPRO", "MV_ENTRY2PROTEIN")
    cur.execute(
        """
        CREATE TABLE INTERPRO.MV_ENTRY2PROTEIN (
            ENTRY_AC NOT NULL VARCHAR2(9),
            PROTEIN_AC NOT NULL VARCHAR2(15),
            MATCH_COUNT NOT NULL NUMBER(7)
        ) NOLOGGING
        """
    )
    cur.execute(
        """
        INSERT /*+ APPEND */ INTO INTERPRO.MV_ENTRY2PROTEIN
        SELECT E.ENTRY_AC, M.PROTEIN_AC, COUNT(*)
        FROM INTERPRO.ENTRY2MATCH E
        INNER JOIN INTERPRO.MV_METHOD2PROTEIN M
            ON E.METHOD_AC = M.METHOD_AC
        GROUP BY E.ENTRY_AC, M.PROTEIN_AC
        """
    )
    con.commit()

    logger.info("indexing MV_ENTRY2PROTEIN")
    cur.execute(
        """
        CREATE INDEX I_MV_ENTRY2PROTEIN
        ON INTERPRO.MV_ENTRY2PROTEIN (ENTRY_AC)
        """
    )

    cur.close()
    con.close()
