import cx_Oracle

from .. import logger, orautils


def update(user: str, dsn: str):
    logger.info("creating UNIPARC.XREF")
    
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    
    # XREF_OLD: legacy table, just to be sure it does not exist any more
    orautils.drop_table(cur, "UNIPARC", "XREF_OLD"
    orautils.drop_table(cur, "UNIPARC", "XREF_OLD")

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
                        
    cur.execute(
        """
        CREATE INDEX I_XREF$UPI
        ON UNIPARC.XREF(UPI)
        TABLESPACE UNIPARC_IND
        NOLOGGING
        """
    )
    cur.execute(
        """
        CREATE INDEX I_XREF$AC
        ON UNIPARC.XREF(AC)
        TABLESPACE UNIPARC_IND
        NOLOGGING
        """
    )
    cur.execute(
        """
        CREATE INDEX I_XREF$DBID
        ON UNIPARC.XREF(DBID)
        TABLESPACE UNIPARC_IND
        NOLOGGING
        """
    )

    logger.info("creating UNIPARC.CV_DATABASE")
    try:
        cur.execute("DROP TABLE UNIPARC.CV_DATABASE")
    except cx_Oracle.DatabaseError as e:
        if e.args[0].code != 942:
            # ORA-00942: table or view does not exist
            raise e
    cur.execute(
        """
        CREATE TABLE UNIPARC.CV_DATABASE TABLESPACE UNIPARC_TAB AS
        SELECT *
        FROM UNIPARC.CV_DATABASE@UAREAD
        """
    )
    orautils.grant(cur, "UNIPARC", "CV_DATABASE", "SELECT", "PUBLIC")
    cur.execute(
        """
        CREATE UNIQUE INDEX PK_CV_DATABASE
        ON UNIPARC.CV_DATABASE (ID)
        TABLESPACE UNIPARC_IND
        """
    )
    cur.execute(
        """
        CREATE UNIQUE INDEX UQ_CV_DATABASE$DESCR
        ON UNIPARC.CV_DATABASE (DESCR)
        TABLESPACE UNIPARC_IND
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
    cur.close()
    con.close()
