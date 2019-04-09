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
            METHOD_AC VARCHAR2(25) NOT NULL ,
            PROTEIN_AC VARCHAR2(15) NOT NULL,
            MATCH_COUNT NUMBER(7) NOT NULL
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
            METHOD_AC VARCHAR2(25) NOT NULL,
            PROTEIN_COUNT NUMBER(8) NOT NULL,
            MATCH_COUNT NUMBER(8) NOT NULL
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

    logger.info("creating MV_ENTRY2PROTEIN_TRUE")
    orautils.drop_table(cur, "INTERPRO", "MV_ENTRY2PROTEIN_TRUE")
    cur.execute(
        """
        CREATE TABLE INTERPRO.MV_ENTRY2PROTEIN_TRUE (
            ENTRY_AC VARCHAR2(9) NOT NULL,
            PROTEIN_AC VARCHAR2(15) NOT NULL,
            MATCH_COUNT NUMBER(7) NOT NULL
        ) NOLOGGING
        """
    )
    cur.execute(
        """
        INSERT /*+ APPEND */ INTO INTERPRO.MV_ENTRY2PROTEIN_TRUE
        SELECT E.ENTRY_AC, M.PROTEIN_AC, COUNT(*)
        FROM INTERPRO.ENTRY2METHOD E
        INNER JOIN INTERPRO.MV_METHOD2PROTEIN M
            ON E.METHOD_AC = M.METHOD_AC
        GROUP BY E.ENTRY_AC, M.PROTEIN_AC
        """
    )
    con.commit()

    logger.info("indexing MV_ENTRY2PROTEIN_TRUE")
    cur.execute(
        """
        CREATE INDEX I_MV_ENTRY2PROTEIN_TRUE
        ON INTERPRO.MV_ENTRY2PROTEIN (ENTRY_AC)
        """
    )

    cur.close()
    con.close()
