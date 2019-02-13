import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Optional

import cx_Oracle

from .. import logger, orautils


def get_max_upi(cur: cx_Oracle.Cursor, analysis_id: int) -> str:
    upis = []

    cur.execute(
        """
        SELECT MIN(JOB_START)
        FROM RUNNING_JOBS@ISPRO
        WHERE ANALYSIS_ID = :1
        """, (analysis_id,)
    )
    row = cur.fetchone()
    if row:
        upis.append(row[0])

    cur.execute(
        """
        SELECT MIN(JOB_START)
        FROM COMPLETED_JOBS@ISPRO
        WHERE ANALYSIS_ID = :1
        AND PERSISTED = 0
        """, (analysis_id,)
    )
    row = cur.fetchone()
    if row:
        upis.append(row[0])

    cur.execute(
        """
        SELECT MAX(JOB_END)
        FROM COMPLETED_JOBS@ISPRO
        WHERE ANALYSIS_ID = :1
        AND PERSISTED = 1
        """, (analysis_id,)
    )
    row = cur.fetchone()
    if row:
        upis.append(row[0])

    cur.execute(
        """
        SELECT MAX(JOB_END)
        FROM COMPLETED_PERSIST_JOBS@ISPRO
        WHERE ANALYSIS_ID = :1
        AND PERSISTED = 1
        """, (analysis_id,)
    )
    row = cur.fetchone()
    if row:
        upis.append(row[0])

    upis = [v for v in upis if v is not None]
    return max(upis) if upis else None


def get_analysis_max_upi(url: str, table: str) -> str:
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.execute(
        """
        SELECT MAX(UPI)
        FROM IPRSCAN.{}@ISPRO
        """.format(table)
    )
    row = cur.fetchone()
    cur.close()
    con.close()
    return row[0]


def check_ispro(url: str, max_upi_read: str) -> Optional[dict]:
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.execute(
        """
        SELECT ANALYSIS_ID, ANALYSIS_NAME, MATCH_TABLE
        FROM (
          SELECT
            ANALYSIS_ID, 
            ANALYSIS_NAME, 
            MATCH_TABLE, 
            ROW_NUMBER() OVER (
              PARTITION 
              BY MATCH_TABLE 
              ORDER BY ANALYSIS_ID DESC
            ) CNT
          FROM IPM_ANALYSIS@ISPRO
            INNER JOIN IPM_ANALYSIS_MATCH_TABLE@ISPRO
              ON ANALYSIS_MATCH_TABLE_ID = ID
          WHERE ACTIVE = 1
        )
        WHERE CNT = 1

        """
    )

    analyses = []
    for row in cur:
        analyses.append({
            "id": row[0],
            "name": row[1],
            "table": row[2]
        })

    cur.close()
    con.close()

    with ThreadPoolExecutor() as executor:
        fs = {}
        for e in analyses:
            f = executor.submit(get_analysis_max_upi, url, e["table"])
            fs[f] = e

        for f in as_completed(fs):
            e = fs[f]
            upi = None

            try:
                upi = f.result()
            except Exception as exc:
                logger.error("{} ({}) exited: "
                             "{}".format(e["name"], e["table"], exc))
            finally:
                fs[f]["upi"] = upi

    tables = {}
    for e in analyses:
        if e["upi"] is None or e["upi"] < max_upi_read:
            return None
        elif e["table"] in tables:
            tables[e["table"]].append(e["id"])
        else:
            tables[e["table"]] = [e["id"]]

    return tables


def import_ispro(url):
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.execute("SELECT MAX(UPI) FROM UNIPARC.PROTEIN@UAREAD")
    max_upi_read = cur.fetchone()[0]
    cur.close()
    con.close()

    while True:
        tables = check_ispro(url, max_upi_read)
        if tables:
            break
        time.sleep(60 * 15)

    with ThreadPoolExecutor() as executor:
        fs = {}
        for table in tables:
            mview = "MV_" + table
            f = executor.submit(orautils.refresh_mview, url, mview)
            fs[f] = mview

        success = True
        for f in as_completed(fs):
            mview = fs[f]
            e = f.exception()
            if e:
                success = False
                logger.error("{}: {}".format(mview, e))
            else:
                logger.info("{}: refreshed".format(mview))

    if not success:
        raise RuntimeError()

    # con = cx_Oracle.connect(url)
    # cur = con.cursor()
    #
    # for mview in fs.values():
    #     cur.execute(
    #         """
    #         ALTER TABLE IPRSCAN.MV_IPRSCAN
    #         EXCHANGE PARTITION {} WITH TABLE {}
    #         INCLUDING INDEXES
    #         WITHOUT VALIDATION
    #         """.format()
    #     )
    #
    # cur.close()
    # con.close()
    #
    #


def import_mv_iprscan(url_src, url_dst):
    obj = orautils.parse_url(url_src)

    con = cx_Oracle.connect(url_dst)
    cur = con.cursor()
    logger.info("DROP TABLE")
    cur.execute("DROP TABLE IPRSCAN.MV_IPRSCAN")

    orautils.create_db_link(cur, "IPPRO", obj["username"], obj["password"],
                            "{}:{}/{}".format(obj["host"], obj["port"],
                                              obj["service"])
                            )

    logger.info("CREATE TABLE")
    cur.execute(
        """
        CREATE TABLE IPRSCAN.MV_IPRSCAN 
        TABLESPACE IPRSCAN_TAB 
        AS
        SELECT * 
        FROM IPRSCAN.MV_IPRSCAN@IPPRO
        """
    )

    logger.info("MV_IPRSCAN_ANALYSIS_ID_MAJORX")
    cur.execute(
        """
        CREATE INDEX IPRSCAN.MV_IPRSCAN_ANALYSIS_ID_MAJORX
        ON IPRSCAN.MV_IPRSCAN (ANALYSIS_ID, RELNO_MAJOR)
        TABLESPACE IPRSCAN_IND
        """
    )

    logger.info("MV_IPRSCAN_ANALYSIS_ID_UPIX")
    cur.execute(
        """
        CREATE INDEX IPRSCAN.MV_IPRSCAN_ANALYSIS_ID_UPIX
        ON IPRSCAN.MV_IPRSCAN (ANALYSIS_ID, UPI)
        TABLESPACE IPRSCAN_IND
        """
    )

    logger.info("MV_IPRSCAN_UPI_METHOD_ACX")
    cur.execute(
        """
        CREATE INDEX IPRSCAN.MV_IPRSCAN_UPI_METHOD_ACX
        ON IPRSCAN.MV_IPRSCAN (UPI, METHOD_AC)
        TABLESPACE IPRSCAN_IND
        """
    )

    logger.info("MV_IPRSCAN_UPI_METHOD_AN_IDX")
    cur.execute(
        """
        CREATE INDEX IPRSCAN.MV_IPRSCAN_UPI_METHOD_AN_IDX
        ON IPRSCAN.MV_IPRSCAN (UPI, METHOD_AC, ANALYSIS_ID)
        TABLESPACE IPRSCAN_IND
        """
    )

    cur.close()
    con.close()


def prepare_matches(url: str):
    con = cx_Oracle.connect(url)
    cur = con.cursor()

    logger.info("adding new matches")
    cur.execute("TRUNCATE TABLE INTERPRO.MATCH_NEW")
    cur.execute(
        """
        INSERT /*+ APPEND */ INTO INTERPRO.MATCH_NEW (
          PROTEIN_AC, METHOD_AC, POS_FROM, POS_TO, STATUS, 
          DBCODE, EVIDENCE,
          SEQ_DATE, MATCH_DATE, TIMESTAMP, USERSTAMP, 
          SCORE, MODEL_AC, FRAGMENTS
        ) 
        SELECT 
          P.PROTEIN_AC, M.METHOD_AC, M.SEQ_START, M.SEQ_END, 'T', 
          D.DBCODE, D.EVIDENCE,
          SYSDATE, SYSDATE, SYSDATE, 'INTERPRO',
          M.EVALUE, M.MODEL_AC, M.FRAGMENTS
        FROM INTERPRO.PROTEIN_TO_SCAN P
        INNER JOIN IPRSCAN.MV_IPRSCAN M 
          ON P.UPI = M.UPI
        INNER JOIN INTERPRO.IPRSCAN2DBCODE D 
          ON M.ANALYSIS_ID = D.IPRSCAN_SIG_LIB_REL_ID
        WHERE D.DBCODE NOT IN ('g', 'j', 'n', 'q', 's', 'v', 'x')
        AND M.SEQ_START != M.SEQ_END
        """
    )
    con.commit()

    logger.info("deleting duplicated SUPERFAMILY matches")
    cur.execute(
        """
        DELETE FROM INTERPRO.MATCH_NEW M1
        WHERE EXISTS(
          SELECT 1
          FROM INTERPRO.MATCH_NEW M2
          WHERE M2.DBCODE = 'Y'
          AND M1.PROTEIN_AC = M2.PROTEIN_AC
          AND M1.METHOD_AC = M2.METHOD_AC
          AND M1.POS_FROM = M2.POS_FROM
          AND M1.POS_TO = M2.POS_TO
          AND M1.SCORE > M2.SCORE
        )
        """
    )
    con.commit()

    for idx in orautils.get_indexes(cur, "INTERPRO", "MATCH_NEW"):
        logger.debug("rebuilding {}".format(idx))
        cur.execute("ALTER INDEX INTERPRO.{} REBUILD".format(idx))
        cur.execute("ALTER INDEX INTERPRO.{} COMPUTE STATISTICS".format(idx))

    cur.close()
    con.close()


def check_matches(url: str):
    con = cx_Oracle.connect(url)
    cur = con.cursor()

    # Match outside of the protein
    cur.execute(
        """
        SELECT COUNT(*)
        FROM INTERPRO.MATCH_NEW M
        INNER JOIN INTERPRO.PROTEIN P
          ON M.PROTEIN_AC = P.PROTEIN_AC
        WHERE M.POS_TO > P.LEN
        """
    )
    n = cur.fetchone()[0]
    if n:
        cur.close()
        con.close()
        raise RuntimeError("{} out-of-bound matches".format(n))

    # Match with invalid start/end positions
    cur.execute(
        """
        SELECT COUNT(*)
        FROM INTERPRO.MATCH_NEW
        WHERE POS_FROM < 1 OR POS_FROM > POS_TO
        """
    )
    if n:
        cur.close()
        con.close()
        raise RuntimeError("{} invalid matches".format(n))

    # Signature matching proteins for the 1st time
    cur.execute(
        """
        SELECT METHOD_AC, COUNT(DISTINCT PROTEIN_AC) CNT
        FROM INTERPRO.MATCH_NEW
        WHERE METHOD_AC IN (
          SELECT METHOD_AC FROM INTERPRO.METHOD
          MINUS 
          SELECT DISTINCT METHOD_AC FROM INTERPRO.MATCH
        )
        GROUP BY METHOD_AC
        """
    )
    signatures = cur.fetchall()

    cur.close()
    con.close()
