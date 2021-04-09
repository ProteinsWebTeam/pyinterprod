# -*- coding: utf-8 -*-

import cx_Oracle

from pyinterprod import logger
from pyinterprod.utils import oracle


def update(url: str):
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    update_databases(cur)
    update_proteins(cur)
    update_xrefs(cur)
    cur.close()
    con.close()
    logger.info("complete")


def update_databases(cur: cx_Oracle.Cursor):
    logger.info("creating table CV_DATABASE")
    oracle.drop_table(cur, "UNIPARC.CV_DATABASE", purge=True)
    cur.execute(
        """
        CREATE TABLE UNIPARC.CV_DATABASE
        TABLESPACE UNIPARC_TAB
        NOLOGGING
        AS
        SELECT *
        FROM UNIPARC.CV_DATABASE@UAREAD
        """
    )
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


def update_proteins(cur: cx_Oracle.Cursor):
    logger.info("creating table PROTEIN")
    oracle.drop_mview(cur, "UNIPARC.PROTEIN")
    oracle.drop_table(cur, "UNIPARC.PROTEIN", purge=True)
    cur.execute(
        """
        CREATE TABLE UNIPARC.PROTEIN
        TABLESPACE UNIPARC_TAB
        NOLOGGING
        AS
        SELECT UPI, TIMESTAMP, USERSTAMP, CRC64, LEN, SEQ_SHORT, SEQ_LONG, MD5
        FROM UNIPARC.PROTEIN@UAREAD
        """
    )
    cur.execute("GRANT SELECT ON UNIPARC.PROTEIN TO PUBLIC")

    logger.info("creating index PK_PROTEIN")
    cur.execute(
        """
        CREATE UNIQUE INDEX PK_PROTEIN
        ON UNIPARC.PROTEIN (UPI)
        TABLESPACE UNIPARC_IND
        """
    )

    logger.info("gathering statistics on table PROTEIN")
    oracle.gather_stats(cur, "UNIPARC", "PROTEIN")


def update_xrefs(cur: cx_Oracle.Cursor):
    logger.info("creating table XREF")
    oracle.drop_table(cur, "UNIPARC.XREF", purge=True)
    oracle.drop_table(cur, "UNIPARC.XREF_OLD", purge=True)
    cur.execute(
        """
        CREATE TABLE UNIPARC.XREF
        TABLESPACE UNIPARC_TAB
        NOLOGGING
        AS
        SELECT UPI, AC, DBID, DELETED, VERSION
        FROM UNIPARC.XREF@UAREAD
        """
    )
    cur.execute("GRANT SELECT ON UNIPARC.XREF TO PUBLIC")

    for col in ["UPI", "AC", "DBID"]:
        logger.info(f"creating index I_XREF${col}")
        cur.execute(
            f"""
            CREATE INDEX I_XREF${col}
            ON UNIPARC.XREF ({col})
            TABLESPACE UNIPARC_IND
            """
        )

    logger.info("gathering statistics on table XREF")
    oracle.gather_stats(cur, "UNIPARC", "XREF")
