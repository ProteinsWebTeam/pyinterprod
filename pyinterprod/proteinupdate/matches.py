from concurrent.futures import ThreadPoolExecutor, as_completed

import cx_Oracle

from .. import logger
from ..orautils import create_db_link


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


def check_ispro(url: str, iprscan_user: str, ispro_url: str,
                uniparc_user: str, uaread_url: str):

    con = cx_Oracle.connect(url)
    cur = con.cursor()
    create_db_link(cur, "ISPRO", *iprscan_user.split('/'), ispro_url)
    create_db_link(cur, "UAREAD", *uniparc_user.split('/'), uaread_url)

    cur.execute(
        """
        SELECT 
          ANALYSIS_ID, ANALYSIS_NAME, MATCH_TABLE
        FROM IPM_ANALYSIS@ISPRO
        INNER JOIN IPM_ANALYSIS_MATCH_TABLE@ISPRO
          ON ANALYSIS_MATCH_TABLE_ID = ID
        WHERE ACTIVE = 1
        """
    )

    analyses = []
    for row in cur:
        analyses.append({
            "id": row[0],
            "name": row[1],
            "table": row[2]
        })

    cur.execute("SELECT MAX(UPI) FROM UNIPARC.PROTEIN@UAREAD")
    max_upi_read = cur.fetchone()[0]

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

    print(max_upi_read)

    for e in analyses:
        if e["upi"] < max_upi_read:
            print(e["name"], e["upi"])

    return all([e["upi"] and e["upi"] >= max_upi_read for e in analyses])
