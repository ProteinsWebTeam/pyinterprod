import cx_Oracle


def update(url: str):
    con = cx_Oracle.connect(url)
    cur = con.cursor()

    cur.execute("DROP TABLE UNIPARC.XREF")
    cur.execute(
        """
        CREATE TABLE UNIPARC.XREF TABLESPACE UNIPARC_TAB AS
        SELECT UPI, AC, DBID, DELETED, VERSION
        FROM UNIPARC.XREF@UAREAD
        """
    )
    cur.execute(
        """
        CREATE INDEX I_XREF$AC
        ON XREF(AC)
        TABLESPACE UNIPARC_IND
        """
    )
    cur.execute(
        """
        CREATE INDEX I_XREF$UPI
        ON XREF(UPI)
        TABLESPACE UNIPARC_IND
        """
    )
    cur.execute("GRANT SELECT ON UNIPARC.XREF TO PUBLIC")

    cur.execute("DROP TABLE UNIPARC.CV_DATABASE")
    cur.execute(
        """
        CREATE TABLE UNIPARC.CV_DATABASE TABLESPACE UNIPARC_TAB AS
        SELECT *
        FROM UNIPARC.CV_DATABASE@UAREAD
        """
    )
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
    cur.execute("GRANT SELECT ON UNIPARC.CV_DATABASE TO PUBLIC")

    cur.close()
    con.close()
