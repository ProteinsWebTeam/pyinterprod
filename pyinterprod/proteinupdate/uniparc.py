import cx_Oracle

from .. import logger


def update(url: str):
    con = cx_Oracle.connect(url)
    cur = con.cursor()

    logger.info("creating UNIPARC.XREF")
    # XREF_OLD: legacy table, just to be sure it does not exist any more
    try:
        cur.execute("DROP TABLE UNIPARC.XREF_OLD")
    except cx_Oracle.DatabaseError as e:
        if e.args[0].code != 942:
            # ORA-00942: table or view does not exist
            raise e

    try:
        cur.execute("DROP TABLE UNIPARC.XREF")
    except cx_Oracle.DatabaseError as e:
        if e.args[0].code != 942:
            # ORA-00942: table or view does not exist
            raise e

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

    cur.execute(
        """
        SELECT COUNT(*)
        FROM UNIPARC.XREF_OLD UX
        INNER JOIN INTERPRO.PROTEIN IP
          ON UX.AC = IP.PROTEIN_AC
        INNER JOIN UNIPARC.PROTEIN UP
          ON UX.UPI = UP.UPI
        WHERE UX.DELETED = 'N' 
        AND IP.CRC64 != UP.CRC64 
        """
    )
    n_rows = cur.fetchone()[0]
    logger.info("{} proteins with mismatched CRC64".format(n_rows))

    if n_rows:
        cur.execute(
            """
            DELETE FROM INTERPRO.XREF_SUMMARY
            WHERE PROTEIN_AC IN (
                SELECT UX.AC
                FROM UNIPARC.XREF_OLD UX
                INNER JOIN INTERPRO.PROTEIN IP
                  ON UX.AC = IP.PROTEIN_AC
                INNER JOIN UNIPARC.PROTEIN UP
                  ON UX.UPI = UP.UPI
                WHERE UX.DELETED = 'N' 
                AND IP.CRC64 != UP.CRC64 
            )
            """
        )

        logger.info("{} rows deleted".format(cur.rowcount))
        con.commit()

    cur.close()
    con.close()
