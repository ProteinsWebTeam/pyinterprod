import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Optional

import cx_Oracle

from .. import logger
from ..orautils import create_db_link, refresh_mview, parse_url


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
            f = executor.submit(refresh_mview, url, mview)
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
    obj = parse_url(url_src)

    con = cx_Oracle.connect(url_dst)
    cur = con.cursor()
    cur.execute("TRUNCATE TABLE IPRSCAN.MV_IPRSCAN")

    create_db_link(cur, "IPPRO", obj["username"], obj["password"],
                   "{}:{}/{}".format(obj["host"], obj["port"], obj["service"])
                   )

    cur.execute(
        """
        INSERT /*+ APPEND NOLOGGING */ INTO IPRSCAN.MV_IPRSCAN
        SELECT * FROM IPRSCAN.MV_IPRSCAN@IPPRO
        """
    )

    con.commit()
    cur.close()
    con.close()
