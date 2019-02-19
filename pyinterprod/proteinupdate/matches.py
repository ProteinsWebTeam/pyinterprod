import json
import os
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Optional

import cx_Oracle

from .. import logger, orautils


def get_max_upi(url: str, analysis_id: int) -> Optional[str]:
    upis = []
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.execute(
        """
        SELECT MAX(HWM_SUBMITTED)
        FROM IPRSCAN.HWM@ISPRO
        WHERE ANALYSIS_ID = :1
        """, (analysis_id,)
    )
    row = cur.fetchone()
    if row:
        upis.append(row[0])

    cur.execute(
        """
        SELECT MAX(JOB_END)
        FROM IPRSCAN.IPM_RUNNING_JOBS@ISPRO
        WHERE ANALYSIS_ID = :1
        """, (analysis_id,)
    )
    row = cur.fetchone()
    if row:
        upis.append(row[0])

    cur.execute(
        """
        SELECT MIN(JOB_END)
        FROM (
            SELECT MIN(JOB_END) JOB_END
            FROM IPRSCAN.IPM_COMPLETED_JOBS@ISPRO
            WHERE ANALYSIS_ID = :analysisid
            AND PERSISTED = 0
            UNION ALL
            SELECT MAX(JOB_END) JOB_END
            FROM IPRSCAN.IPM_COMPLETED_JOBS@ISPRO
            WHERE ANALYSIS_ID = :analysisid
            AND PERSISTED = 1
        )
        """, dict(analysisid=analysis_id)
    )
    row = cur.fetchone()
    if row:
        upis.append(row[0])

    cur.execute(
        """
        SELECT MIN(JOB_END)
        FROM (
            SELECT MIN(JOB_END) JOB_END
            FROM IPRSCAN.IPM_PERSISTED_JOBS@ISPRO
            WHERE ANALYSIS_ID = :analysisid
            AND PERSISTED = 0
            UNION ALL
            SELECT MAX(JOB_END) JOB_END
            FROM IPRSCAN.IPM_PERSISTED_JOBS@ISPRO
            WHERE ANALYSIS_ID = :analysisid
            AND PERSISTED = 1
        )
        """, dict(analysisid=analysis_id)
    )
    row = cur.fetchone()
    if row:
        upis.append(row[0])

    try:
        return max([upi for upi in upis if upi is not None])
    except ValueError:
        return None


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
            f = executor.submit(get_max_upi, url, e["id"])
            fs[f] = e

        for f in as_completed(fs):
            e = fs[f]
            upi = None

            try:
                upi = f.result()
            except Exception as exc:
                logger.error("{} ({}) exited: "
                             "{}".format(e["name"], e["id"], exc))
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

    logger.info("SUPERFAMILY: deleting duplicated matches")
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
    logger.info("SUPERFAMILY: {} matches deleted".format(cur.rowcount))
    con.commit()

    for idx in orautils.get_indexes(cur, "INTERPRO", "MATCH_NEW"):
        logger.debug("rebuilding {}".format(idx))
        cur.execute("ALTER INDEX INTERPRO.{} REBUILD".format(idx))
        cur.execute("ALTER INDEX INTERPRO.{} COMPUTE STATISTICS".format(idx))

    cur.close()
    con.close()


def get_match_counts(cur: cx_Oracle.Cursor) -> tuple:
    logger.info("counting entries protein counts")
    cur.execute(
        """
        SELECT E.ENTRY_AC, COUNT(DISTINCT M.PROTEIN_AC)
        FROM INTERPRO.ENTRY2METHOD E
        INNER JOIN INTERPRO.MATCH M
          ON E.METHOD_AC = M.METHOD_AC
        GROUP BY E.ENTRY_AC
        """
    )
    entries = dict(cur.fetchall())

    logger.info("counting databases match counts")
    cur.execute(
        """
        SELECT DBCODE, COUNT(*)
        FROM INTERPRO.MATCH
        GROUP BY DBCODE
        """
    )
    databases = dict(cur.fetchall())

    return entries, databases


def check_matches(url: str, outdir: str):
    os.makedirs(outdir, exist_ok=True)

    con = cx_Oracle.connect(url)
    cur = con.cursor()

    # Match outside of the protein
    logger.info("checking out-of-bound matches")
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
    logger.info("checking invalid matches")
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

    entries, databases = get_match_counts(cur)
    cur.close()
    con.close()

    with open(os.path.join(outdir, "counts.json"), "wt") as fh:
        json.dump(dict(entries=entries, databases=databases), fh)


def update_matches(url: str):
    logger.info("updating matches")
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.execute(
        """
        DELETE FROM INTERPRO.MATCH
        WHERE PROTEIN_AC IN (
          SELECT PROTEIN_AC
          FROM INTERPRO.PROTEIN_TO_SCAN
        )
        """
    )
    con.commit()
    logger.info("{} matches deleted".format(cur.rowcount))

    cur.execute(
        """
        INSERT INTO INTERPRO.MATCH
        SELECT * FROM INTERPRO.MATCH_NEW
        """
    )
    con.commit()
    logger.info("{} matches inserted".format(cur.rowcount))

    cur.close()
    con.close()


def update_feature_matches(url: str):
    logger.info("updating feature matches")
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.execute(
        """
        DELETE FROM INTERPRO.FEATURE_MATCH
        WHERE PROTEIN_AC IN (
          SELECT PROTEIN_AC
          FROM INTERPRO.PROTEIN_TO_SCAN
        )
        """
    )
    con.commit()
    logger.info("{} feature matches deleted".format(cur.rowcount))

    cur.execute(
        """
        INSERT INTO INTERPRO.FEATURE_MATCH (
          PROTEIN_AC, METHOD_AC, SEQ_FEATURE, POS_FROM, POS_TO,
          DBCODE, SEQ_DATE, MATCH_DATE, TIMESTAMP, USERSTAMP
        )
        SELECT
          P.PROTEIN_AC, M.METHOD_AC, M.SEQ_FEATURE, M.SEQ_START, M.SEQ_END,
          D.DBCODE, SYSDATE, SYSDATE, SYSDATE, 'INTERPRO'
        FROM INTERPRO.PROTEIN_TO_SCAN P
        INNER JOIN IPRSCAN.MV_IPRSCAN M
          ON P.UPI = M.UPI
        INNER JOIN INTERPRO.IPRSCAN2DBCODE D
          ON M.ANALYSIS_ID = D.IPRSCAN_SIG_LIB_REL_ID
        WHERE D.DBCODE IN ('g', 'j', 'n', 'q', 's', 'v', 'x')
        """
    )
    con.commit()
    logger.info("{} feature matches inserted".format(cur.rowcount))

    cur.close()
    con.close()


def post_matches_update(url: str, outdir: str):
    with open(os.path.join(outdir, "counts.json"), "rt") as fh:
        prev = json.load(fh)

    con = cx_Oracle.connect(url)
    cur = con.cursor()
    entries, databases = get_match_counts(cur)
    cur.close()
    con.close()

    changes = []
    for entry_acc, new_count in entries.items():
        prev_count = prev["entries"].pop(entry_acc, 0)

        try:
            change = (new_count - prev_count) / prev_count
        except ZeroDivisionError:
            changes.append((entry_acc, prev_count, new_count, "N/A"))
        else:
            if abs(change) >= 0.5:
                changes.append((entry_acc, prev_count, new_count, change))

    with open(os.path.join(outdir, "entries_changes.tsv"), "wt") as fh:
        def _sort_entries(e):
            return 1 if e[3] == "N/A" else 0, e[3], e[0]

        fh.write("# Accession\tPrevious protein count\t"
                 "New protein count\tChange (%)\n")

        for ac, pc, nc, c in sorted(changes, key=_sort_entries):
            if c != "N/A":
                c = round(c * 100, 2)

            fh.write("{}\t{}\t{}\t{}\n".format(ac, pc, nc, c))

    changes = []
    for dbcode, new_count in databases.items():
        prev_count = prev["databases"].pop(dbcode, 0)
        changes.append((dbcode, prev_count, new_count))

    for dbcode, prev_count in prev["databases"].items():
        changes.append((dbcode, prev_count, 0))

    with open(os.path.join(outdir, "databases_changes.tsv"), "wt") as fh:
        fh.write("# Code\tPrevious match count\tNew match count\n")

        for dbcode, pc, nc in sorted(changes):
            fh.write("{}\t{}\t{}\n".format(dbcode, pc, nc))
