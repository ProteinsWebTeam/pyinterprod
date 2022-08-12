from typing import Optional

import cx_Oracle

from pyinterprod import logger
from pyinterprod.utils import oracle


def update(ipr_uri: str, unp_uri: str):
    update_databases(ipr_uri, unp_uri)
    update_proteins(ipr_uri, unp_uri)
    update_xrefs(ipr_uri, unp_uri)


def update_databases(ipr_uri: str, unp_uri: str):
    logger.info("creating table CV_DATABASE")
    con = cx_Oracle.connect(unp_uri)
    cur = con.cursor()
    cur.execute("SELECT * FROM UNIPARC.CV_DATABASE")
    records = cur.fetchall()
    cur.close()
    con.close()

    con = cx_Oracle.connect(ipr_uri)
    cur = con.cursor()
    oracle.drop_table(cur, "UNIPARC.CV_DATABASE", purge=True)
    cur.execute(
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

    cur.executemany(
        """
            INSERT /*+ APPEND */ 
            INTO UNIPARC.CV_DATABASE 
            VALUES (:1, :2, :3, :4, :5, :6, :7, :8, :9, :10, :11, :12, :13)
        """,
        records
    )

    con.commit()

    cur.execute("GRANT SELECT ON UNIPARC.CV_DATABASE TO PUBLIC")
    cur.execute(
        """
        CREATE UNIQUE INDEX PK_CV_DATABASE
        ON UNIPARC.CV_DATABASE (ID)
        TABLESPACE UNIPARC_IND
        NOLOGGING
        """
    )
    cur.execute(
        """
        CREATE UNIQUE INDEX UQ_CV_DATABASE$DESCR
        ON UNIPARC.CV_DATABASE (DESCR)
        TABLESPACE UNIPARC_IND
        NOLOGGING
        """
    )
    cur.execute(
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
    oracle.gather_stats(cur, "UNIPARC", "CV_DATABASE")
    cur.close()
    con.close()
    logger.info("complete")


def update_proteins(ipr_uri: str, unp_uri: str, top_up: bool = False):
    logger.info("creating table PROTEIN")
    con = cx_Oracle.connect(ipr_uri)
    cur = con.cursor()

    if top_up:
        cur.execute("SELECT MAX(UPI) FROM UNIPARC.PROTEIN")
        max_upi, = cur.fetchone()
    else:
        max_upi = None
        oracle.drop_table(cur, "UNIPARC.PROTEIN", purge=True)
        cur.execute(
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

    cnt = 0
    records = []
    req = """
        INSERT /*+ APPEND */ 
        INTO UNIPARC.PROTEIN
        VALUES (:1, :2, :3, :4, :5, :6, :7, :8)
    """

    for rec in iter_proteins(unp_uri, greather_than=max_upi):
        # First element of record is ID, which we do not need
        records.append(rec[1:])
        cnt += 1

        if len(records) == 1000:
            cur.executemany(req, records)
            con.commit()
            records.clear()

    if records:
        cur.executemany(req, records)
        con.commit()
        records.clear()

    if not top_up:
        cur.execute("GRANT SELECT ON UNIPARC.PROTEIN TO PUBLIC")
        cur.execute(
            """
            CREATE UNIQUE INDEX PK_PROTEIN
            ON UNIPARC.PROTEIN (UPI)
            TABLESPACE UNIPARC_IND
            """
        )

    logger.info("gathering statistics on table PROTEIN")
    oracle.gather_stats(cur, "UNIPARC", "PROTEIN")
    cur.close()
    con.close()
    logger.info("complete")


def update_xrefs(ipr_uri: str, unp_uri: str):
    logger.info("creating table XREF")
    ipr_con = cx_Oracle.connect(ipr_uri)
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
    records = []

    unp_con = cx_Oracle.connect(unp_uri)
    unp_cur = unp_con.cursor()
    unp_cur.execute(
        """
        SELECT UPI, AC, DBID, DELETED, VERSION
        FROM UNIPARC.XREF          
        """
    )
    for rec in unp_cur:
        records.append(rec)

        if len(records) == 1000:
            ipr_cur.executemany(req, records)
            ipr_con.commit()
            records.clear()

    unp_cur.close()
    unp_con.close()

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


def iter_proteins(uri: str, greather_than: Optional[str] = None):
    con = cx_Oracle.connect(uri)
    cur = con.cursor()
    cur.outputtypehandler = oracle.clob_as_str

    sql = """
        SELECT ID, UPI, TIMESTAMP, USERSTAMP, CRC64, LEN, SEQ_SHORT, 
               SEQ_LONG, MD5
        FROM UNIPARC.PROTEIN      
    """

    if greather_than is not None:
        sql += "WHERE UPI > :1"
        params = [greather_than]
    else:
        params = []

    cur.execute(sql, params)
    yield from cur
    cur.close()
    con.close()
