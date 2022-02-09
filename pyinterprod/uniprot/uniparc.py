import cx_Oracle

from pyinterprod import logger
from pyinterprod.utils import oracle, Table


def update(ipr_url: str, unp_url: str):
    update_databases(ipr_url, unp_url)
    update_proteins(ipr_url, unp_url)
    update_xrefs(ipr_url, unp_url)


def update_databases(ipr_url: str, unp_url: str):
    logger.info("creating table CV_DATABASE")
    ipr_con = cx_Oracle.connect(ipr_url)
    ipr_cur = ipr_con.cursor()
    oracle.drop_table(ipr_cur, "UNIPARC.CV_DATABASE", purge=True)
    ipr_cur.execute(
        """
        CREATE TABLE UNIPARC.CV_DATABASE
        (
            ID NUMBER(3) NOT NULL,
            TIMESTAMP DATE NOT NULL,
            USERSTAMP VARCHAR2(30) NOT NULL,
            DESCR VARCHAR2(32) NOT NULL,
            CURRENT_RELEASE NUMBER(6),
            FULL_DESCR VARCHAR2(512),
            ALIVE VARCHAR2(1) NOT NULL,
            FOR_RELEASE CHAR,
            DISPLAY_NAME VARCHAR2(40),
            UNIREF_UNLOAD VARCHAR2(1) NOT NULL,
            UNIPROT VARCHAR2(1) NOT NULL,
            MULTIPLE_TAXID VARCHAR2(1),
            STOP_DAILY_LOAD VARCHAR2(1)
        ) NOLOGGING
        """
    )
    ipr_cur.close()

    req = """
        INSERT /*+ APPEND */ 
        INTO UNIPARC.CV_DATABASE 
        VALUES (:1, :2, :3, :4, :5, :6, :7, :8, :9, :10, :11, :12, :13)
    """
    with Table(ipr_con, req, autocommit=True) as cv_database:
        unp_con = cx_Oracle.connect(unp_url)
        unp_cur = unp_con.cursor()
        unp_cur.execute("SELECT * FROM UNIPARC.CV_DATABASE")
        for rec in unp_cur:
            cv_database.insert(rec)

        unp_cur.close()
        unp_con.close()

    ipr_cur = ipr_con.cursor()
    ipr_cur.execute("GRANT SELECT ON UNIPARC.CV_DATABASE TO PUBLIC")
    ipr_cur.execute(
        """
        CREATE UNIQUE INDEX PK_CV_DATABASE
        ON UNIPARC.CV_DATABASE (ID)
        TABLESPACE UNIPARC_IND
        NOLOGGING
        """
    )
    ipr_cur.execute(
        """
        CREATE UNIQUE INDEX UQ_CV_DATABASE$DESCR
        ON UNIPARC.CV_DATABASE (DESCR)
        TABLESPACE UNIPARC_IND
        NOLOGGING
        """
    )
    ipr_cur.execute(
        """
        ALTER TABLE UNIPARC.CV_DATABASE
        ADD (
          CONSTRAINT PK_CV_DATABASE PRIMARY KEY(ID)
            USING INDEX PK_CV_DATABASE,
          CONSTRAINT UQ_CV_DATABASE$DESCR UNIQUE (DESCR)
            USING INDEX UQ_CV_DATABASE$DESCR
        )
        """
    )
    oracle.gather_stats(ipr_cur, "UNIPARC", "CV_DATABASE")
    ipr_cur.close()
    ipr_con.close()
    logger.info("complete")


def update_proteins(ipr_url: str, unp_url: str):
    logger.info("creating table PROTEIN")
    ipr_con = cx_Oracle.connect(ipr_url)
    ipr_cur = ipr_con.cursor()
    oracle.drop_table(ipr_cur, "UNIPARC.PROTEIN", purge=True)
    ipr_cur.execute(
        """
        CREATE TABLE UNIPARC.PROTEIN
        (
            UPI CHAR(13) NOT NULL,
            TIMESTAMP DATE NOT NULL,
            USERSTAMP VARCHAR2(30) NOT NULL,
            CRC64 CHAR(16) NOT NULL,
            LEN NUMBER(6) NOT NULL,
            SEQ_SHORT VARCHAR2(4000),
            SEQ_LONG CLOB,
            MD5 VARCHAR2(32) NOT NULL
        ) NOLOGGING
        """
    )
    ipr_cur.close()

    req = """
        INSERT /*+ APPEND */ 
        INTO UNIPARC.PROTEIN
        VALUES (:1, :2, :3, :4, :5, :6, :7, :8)
    """
    with Table(ipr_con, req, autocommit=True, buffer_size=1000) as protein:
        unp_con = cx_Oracle.connect(unp_url)
        unp_cur = unp_con.cursor()
        unp_cur.outputtypehandler = oracle.clob_as_str
        unp_cur.execute(
            """
            SELECT UPI, TIMESTAMP, USERSTAMP, CRC64, LEN, SEQ_SHORT, 
                   SEQ_LONG, MD5
            FROM UNIPARC.PROTEIN          
            """
        )
        for rec in unp_cur:
            protein.insert(rec)

        unp_cur.close()
        unp_con.close()

    ipr_cur = ipr_con.cursor()
    ipr_cur.execute("GRANT SELECT ON UNIPARC.PROTEIN TO PUBLIC")

    logger.info("creating index PK_PROTEIN")
    ipr_cur.execute(
        """
        CREATE UNIQUE INDEX PK_PROTEIN
        ON UNIPARC.PROTEIN (UPI)
        TABLESPACE UNIPARC_IND
        """
    )

    logger.info("gathering statistics on table PROTEIN")
    oracle.gather_stats(ipr_cur, "UNIPARC", "PROTEIN")
    ipr_cur.close()
    ipr_con.close()
    logger.info("complete")


def update_xrefs(ipr_url: str, unp_url: str):
    logger.info("creating table XREF")
    ipr_con = cx_Oracle.connect(ipr_url)
    ipr_cur = ipr_con.cursor()
    oracle.drop_table(ipr_cur, "UNIPARC.XREF", purge=True)
    ipr_cur.execute(
        """
        CREATE TABLE UNIPARC.XREF
        (
            UPI CHAR(13) NOT NULL,
            AC VARCHAR2(70) NOT NULL,
            DBID NUMBER(3) NOT NULL,
            DELETED CHAR NOT NULL,
            VERSION NUMBER(6)
        ) NOLOGGING
        """
    )
    ipr_cur.close()

    req = """
        INSERT /*+ APPEND */ 
        INTO UNIPARC.XREF
        VALUES (:1, :2, :3, :4, :5)
    """
    with Table(ipr_con, req, autocommit=True) as xref:
        unp_con = cx_Oracle.connect(unp_url)
        unp_cur = unp_con.cursor()
        unp_cur.execute(
            """
            SELECT UPI, AC, DBID, DELETED, VERSION
            FROM UNIPARC.XREF          
            """
        )
        for rec in unp_cur:
            xref.insert(rec)

        unp_cur.close()
        unp_con.close()

    ipr_cur = ipr_con.cursor()
    ipr_cur.execute("GRANT SELECT ON UNIPARC.XREF TO PUBLIC")

    for col in ["UPI", "AC", "DBID"]:
        logger.info(f"creating index I_XREF${col}")
        ipr_cur.execute(
            f"""
            CREATE INDEX I_XREF${col}
            ON UNIPARC.XREF ({col})
            TABLESPACE UNIPARC_IND
            """
        )

    logger.info("gathering statistics on table XREF")
    oracle.gather_stats(ipr_cur, "UNIPARC", "XREF")
    ipr_cur.close()
    ipr_con.close()
    logger.info("complete")
