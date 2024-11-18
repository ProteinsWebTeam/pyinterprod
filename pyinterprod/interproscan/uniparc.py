import oracledb

from pyinterprod import logger
from pyinterprod.uniprot.uniparc import iter_proteins
from pyinterprod.utils import oracle


def import_sequences(ispro_uri: str, uniparc_uri: str, top_up: bool = False,
                     max_upi: str | None = None):
    logger.info("importing sequences from UniParc")
    con = oracledb.connect(ispro_uri)
    cur = con.cursor()

    if top_up:
        cur.execute("SELECT MAX(UPI) FROM UNIPARC.PROTEIN")
        current_max_upi, = cur.fetchone()
    else:
        current_max_upi = None
        oracle.drop_table(cur, "UNIPARC.PROTEIN", purge=True)
        cur.execute(
            """
                CREATE TABLE UNIPARC.PROTEIN
                (
                    ID NUMBER(15) NOT NULL,
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

    logger.info(f"\thighest UPI: {current_max_upi or 'N/A'}")

    cnt = 0
    records = []
    req = """
        INSERT /*+ APPEND */ 
        INTO UNIPARC.PROTEIN
        VALUES (:1, :2, :3, :4, :5, :6, :7, :8, :9)
    """

    for rec in iter_proteins(uniparc_uri, gt=current_max_upi, le=max_upi):
        records.append(rec)
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
        cur.execute("CREATE UNIQUE INDEX PK_PROTEIN ON UNIPARC.PROTEIN (UPI)")

    cur.execute("SELECT MAX(UPI) FROM UNIPARC.PROTEIN")
    current_max_upi, = cur.fetchone()
    logger.info(f"\tnew highest UPI: {current_max_upi or 'N/A'}")

    cur.close()
    con.close()

    logger.info(f"\t{cnt:,} sequences imported")
