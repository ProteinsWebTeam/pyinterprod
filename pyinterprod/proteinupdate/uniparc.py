# -*- coding: utf-8 -*-

import cx_Oracle

from .. import logger, orautils


def update(user: str, dsn: str):
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    _update_xref(cur)
    _update_database(cur)
    cur.close()
    con.close()


def _update_database(cur: cx_Oracle.Cursor):
    logger.info("creating UNIPARC.CV_DATABASE")

    orautils.drop_table(cur, "UNIPARC", "CV_DATABASE", purge=True)
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
    orautils.grant(cur, "UNIPARC", "CV_DATABASE", "SELECT", "PUBLIC")

    logger.info("analyzing UNIPARC.CV_DATABASE")
    orautils.gather_stats(cur, "UNIPARC", "CV_DATABASE")

    logger.info("indexing UNIPARC.CV_DATABASE")
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

    logger.info("UNIPARC.CV_DATABASE ready")


def _update_protein(user: str, dsn: str):
    logger.info("creating UNIPARC.PROTEIN")

    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    orautils.drop_table(cur, "UNIPARC", "PROTEIN", purge=True)
    orautils.drop_mview(cur, "UNIPARC", "PROTEIN")
    cur.execute(
        """
        CREATE TABLE UNIPARC.PROTEIN
        TABLESPACE UNIPARC_TAB
        NOLOGGING
        AS
        SELECT
            UPI, TIMESTAMP, USERSTAMP, CRC64, LEN, SEQ_SHORT, SEQ_LONG, MD5
        FROM UNIPARC.PROTEIN@UAPRO
        """
    )
    orautils.grant(cur, "UNIPARC", "PROTEIN", "SELECT", "PUBLIC")

    logger.info("analyzing UNIPARC.PROTEIN")
    orautils.gather_stats(cur, "UNIPARC", "PROTEIN")

    logger.info("indexing UNIPARC.PROTEIN")
    cur.execute(
        """
        CREATE UNIQUE INDEX PK_PROTEIN
        ON UNIPARC.PROTEIN (UPI)
        TABLESPACE UNIPARC_IND
        """
    )

    cur.close()
    con.close()

    logger.info("UNIPARC.PROTEIN ready")


def _update_xref(cur: cx_Oracle.Cursor):
    logger.info("creating UNIPARC.XREF")

    # XREF_OLD: legacy table, just to be sure it does not exist any more
    orautils.drop_table(cur, "UNIPARC", "XREF", purge=True)
    orautils.drop_table(cur, "UNIPARC", "XREF_OLD", purge=True)
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
    orautils.grant(cur, "UNIPARC", "XREF", "SELECT", "PUBLIC")

    logger.info("analyzing UNIPARC.XREF")
    orautils.gather_stats(cur, "UNIPARC", "XREF")

    for col in ("UPI", "AC", "DBID"):
        logger.info("creating index INDEX I_XREF${}".format(col))
        cur.execute(
            """
            CREATE INDEX I_XREF${0}
            ON UNIPARC.XREF({0})
            TABLESPACE UNIPARC_IND
            NOLOGGING
            """.format(col)
        )

    logger.info("UNIPARC.XREF ready")
