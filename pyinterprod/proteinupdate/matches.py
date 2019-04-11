import json
import os
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Collection, List, Optional

import cx_Oracle

from . import materializedviews as mviews
from .. import logger, orautils


def _get_max_upi(url: str, analysis_id: int) -> Optional[str]:
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


def _get_ispro_upis(url: str) -> List[dict]:
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    # TODO: select the SITE_TABLE as well
    cur.execute(
        """
        SELECT *
        FROM (
          SELECT
            A.ANALYSIS_ID,
            A.ANALYSIS_NAME,
            A.ANALYSIS_TYPE,
            B.MATCH_TABLE,
            ROW_NUMBER() OVER (
              PARTITION
              BY A.ANALYSIS_MATCH_TABLE_ID
              ORDER BY A.ANALYSIS_ID DESC
            ) CNT
          FROM IPM_ANALYSIS@ISPRO A
            INNER JOIN IPM_ANALYSIS_MATCH_TABLE@ISPRO B
              ON A.ANALYSIS_MATCH_TABLE_ID = B.ID
          WHERE A.ACTIVE = 1
        )
        ORDER BY CNT
        """
    )

    tables = set()
    analyses = []
    for row in cur:
        """
        SignalP has EUK/Gram+/Gram- matches in the same tables so we need
        several times the same time for different analyses.
        
        But we should accept any other duplicates 
        (ACTIVE might not always be up-to-date)
        """
        table = row[3]
        if table in tables and table != "ipm_signalp_match":
            continue

        tables.add(table)
        analyses.append({
            "id": row[0],
            "full_name": row[1],
            "name": row[2],
            "table": table,
            "upi": None
        })

    cur.close()
    con.close()

    with ThreadPoolExecutor() as executor:
        fs = {}
        for e in analyses:
            f = executor.submit(_get_max_upi, url, e["id"])
            fs[f] = e

        for f in as_completed(fs):
            e = fs[f]

            try:
                upi = f.result()
            except Exception as exc:
                logger.error("{} ({}) exited: "
                             "{}".format(e["name"], e["id"], exc))
            else:
                e["upi"] = upi

    return analyses


def _check_ispro(url: str, max_attempts: int=1, secs: int=3600,
                 exclude: Optional[Collection[str]]=None) -> List[dict]:
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.execute("SELECT MAX(UPI) FROM UNIPARC.PROTEIN@UAREAD")
    max_upi_read = cur.fetchone()[0]
    cur.close()
    con.close()

    num_attempts = 0
    while True:
        analyses = _get_ispro_upis(url)
        not_ready = 0
        for e in analyses:
            if exclude and e["name"] in exclude:
                status = "ready (excluded)"
            elif e["upi"] and e["upi"] >= max_upi_read:
                status = "ready"
            else:
                status = "not ready"
                not_ready += 1

            logger.debug("{id:<5}{full_name:<30}{table:<30}{upi:<20}"
                         "{status}".format(**e, status=status))

        num_attempts += 1
        if not_ready == 0:
            break
        elif num_attempts == max_attempts:
            raise RuntimeError("ISPRO tables not ready "
                               "after {} attempts".format(num_attempts))
        else:
            logger.info("{} analyses not ready".format(not_ready))
            time.sleep(secs)

    return analyses


def _import_table(url: str, owner: str, name: str):
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    orautils.drop_mview(cur, owner, name)
    orautils.drop_table(cur, owner, name)
    cur.execute(
        """
        CREATE TABLE {0}.{1} NOLOGGING
        AS
        SELECT *
        FROM {0}.{1}@ISPRO
        """.format(owner, name)
    )
    orautils.gather_stats(cur, owner, name)
    cur.execute(
        """
        CREATE INDEX {0}.I_{1} 
        ON {0}.{1}(ANALYSIS_ID) 
        NOLOGGING
        """.format(owner, name)
    )
    cur.close()
    con.close()


def import_ispro(user: str, dsn: str, **kwargs):
    url = orautils.make_connect_string(user, dsn)
    analyses = _check_ispro(url, **kwargs)

    logger.info("building tables with data from ISPRO")
    with ThreadPoolExecutor() as executor:
        fs = {}
        for analysis in analyses:
            table_name = analysis["table"].upper()
            f = executor.submit(_import_table, url, "IPRSCAN", table_name)
            fs[f] = table_name

        num_errors = 0
        for f in as_completed(fs):
            table_name = fs[f]
            exc = f.exception()
            if exc is None:
                logger.info("{:<30}: refreshed".format(table_name))
            else:
                logger.error("{:<30}: {}".format(table_name, exc))
                num_errors += 1

        if num_errors:
            raise RuntimeError("{} tables "
                               "were not imported".format(num_errors))

        functions = {
            "ipm_cdd_match": mviews.update_cdd,
            "ipm_coils_match": mviews.update_coils,
            "ipm_gene3d_match": mviews.update_gene3d,
            "ipm_hamap_match": mviews.update_hamap,
            "ipm_mobidblite_match": mviews.update_mobidblite,
            "ipm_panther_match": mviews.update_panther,
            "ipm_pfam_match": mviews.update_pfam,
            "ipm_phobius_match": mviews.update_phobius,
            "ipm_pirsf_match": mviews.update_pirsf,
            "ipm_prints_match": mviews.update_prints,
            "ipm_prodom_match": mviews.update_prodom,
            "ipm_prosite_patterns_match": mviews.update_prosite_patterns,
            "ipm_prosite_profiles_match": mviews.update_prosite_profiles,
            "ipm_smart_match": mviews.update_smart,
            "ipm_superfamily_match": mviews.update_superfamily,
            "ipm_tigrfam_match": mviews.update_tigrfam,
            "ipm_tmhmm_match": mviews.update_tmhmm
        }

        signalp = {
            "signalp_euk": mviews.update_signalp_euk,
            "signalp_gram_negative": mviews.update_signalp_gram_neg,
            "signalp_gram_positive": mviews.update_signalp_gram_pos
        }

        fs = {}
        for analysis in analyses:
            _id = analysis["id"]
            name = analysis["name"]
            table_name = analysis["table"]

            if table_name == "ipm_signalp_match":
                fn = signalp[name]
            else:
                fn = functions[table_name]

            f = executor.submit(fn, url, _id)
            fs[f] = analysis["full_name"]

        num_errors = 0
        for f in as_completed(fs):
            full_name = fs[f]
            exc = f.exception()
            if exc is None:
                logger.info("{:<30}: refreshed".format(full_name))
            else:
                logger.error("{:<30}: {}".format(full_name, exc))
                num_errors += 1

        if num_errors:
            raise RuntimeError("{} partitions "
                               "were not exchanged".format(num_errors))


def _import_mv_iprscan(url: str, url_src: str):
    # Get table partitions (source database)
    con = cx_Oracle.connect(url_src)
    cur = con.cursor()
    partitions = orautils.get_partitions(cur, "IPRSCAN", "MV_IPRSCAN")
    column = orautils.get_partitioning_key(cur, "IPRSCAN", "MV_IPRSCAN")
    cur.close()
    con.close()

    # Open connection (target database)
    con = cx_Oracle.connect(url)
    cur = con.cursor()

    # Create link to source database
    link = orautils.parse_url(url_src)
    orautils.create_db_link(cur,
                            link="IPPRO",
                            user=link["username"],
                            passwd=link["password"],
                            dsn="{host}:{port}/{service}".format(**link)
                            )

    # Building TABLE
    logger.info("creating table")
    query = "CREATE TABLE IPRSCAN.MV_IPRSCAN"
    if partitions:
        query += " PARTITION BY LIST ({}) ({})".format(
            column,
            ','.join([
                "PARTITION {name} VALUES ({value})".format(**p)
                for p in partitions
            ])
        )

    query += " AS SELECT * FROM IPRSCAN.MV_IPRSCAN@IPPRO"

    orautils.drop_table(cur, "IPRSCAN", "MV_IPRSCAN")
    cur.execute(query)

    # Building indices
    logger.info("creating indices")
    logger.debug("\tMV_IPRSCAN_ANALYSIS_ID_MAJORX")
    cur.execute(
        """
        CREATE INDEX IPRSCAN.MV_IPRSCAN_ANALYSIS_ID_MAJORX
        ON IPRSCAN.MV_IPRSCAN (ANALYSIS_ID, RELNO_MAJOR)
        TABLESPACE IPRSCAN_IND
        """
    )

    logger.debug("\tMV_IPRSCAN_ANALYSIS_ID_UPIX")
    cur.execute(
        """
        CREATE INDEX IPRSCAN.MV_IPRSCAN_ANALYSIS_ID_UPIX
        ON IPRSCAN.MV_IPRSCAN (ANALYSIS_ID, UPI)
        TABLESPACE IPRSCAN_IND
        """
    )

    logger.debug("\tMV_IPRSCAN_UPI_METHOD_ACX")
    cur.execute(
        """
        CREATE INDEX IPRSCAN.MV_IPRSCAN_UPI_METHOD_ACX
        ON IPRSCAN.MV_IPRSCAN (UPI, METHOD_AC)
        TABLESPACE IPRSCAN_IND
        """
    )

    logger.debug("\tMV_IPRSCAN_UPI_METHOD_AN_IDX")
    cur.execute(
        """
        CREATE INDEX IPRSCAN.MV_IPRSCAN_UPI_METHOD_AN_IDX
        ON IPRSCAN.MV_IPRSCAN (UPI, METHOD_AC, ANALYSIS_ID)
        TABLESPACE IPRSCAN_IND
        """
    )

    cur.close()
    con.close()


def prepare_matches(user: str, dsn: str):
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
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


def _get_match_counts(cur: cx_Oracle.Cursor) -> tuple:
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


def check_matches(user: str, dsn: str, outdir: str):
    os.makedirs(outdir, exist_ok=True)

    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
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

    entries, databases = _get_match_counts(cur)
    cur.close()
    con.close()

    with open(os.path.join(outdir, "counts.json"), "wt") as fh:
        json.dump(dict(entries=entries, databases=databases), fh)


def update_matches(user: str, dsn: str):
    logger.info("updating matches")
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
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


def update_feature_matches(user: str, dsn: str):
    logger.info("updating feature matches")
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()

    # Make sure legacy tables are dropped
    for t in ("FEATURE_MATCH_NEW", "FEATURE_MATCH_NEW_STG"):
        orautils.drop_table(cur, "INTERPRO", t)

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


def track_count_changes(user: str, dsn: str, outdir: str):
    with open(os.path.join(outdir, "counts.json"), "rt") as fh:
        prev = json.load(fh)

    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    entries, databases = _get_match_counts(cur)
    cur.close()
    con.close()

    changes = []
    for entry_acc, new_count in entries.items():
        prev_count = prev["entries"].pop(entry_acc, 0)

        try:
            change = (new_count - prev_count) / prev_count * 100
        except ZeroDivisionError:
            changes.append((entry_acc, prev_count, new_count, "N/A"))
        else:
            if abs(change) >= 50:
                changes.append((entry_acc, prev_count, new_count, change))

    with open(os.path.join(outdir, "entries_changes.tsv"), "wt") as fh:
        def _sort_entries(e):
            return 1 if isinstance(e[3], str) else 0, e[3], e[0]

        fh.write("# Accession\tPrevious protein count\t"
                 "New protein count\tChange (%)\n")

        for ac, pc, nc, c in sorted(changes, key=_sort_entries):
            if isinstance(c, str):
                fh.write("{}\t{}\t{}\t{}\n".format(ac, pc, nc, c))
            else:
                fh.write("{}\t{}\t{}\t{:.0f}\n".format(ac, pc, nc, c))

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


def update_alt_splicing_matches(user: str, dsn: str):
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    logger.info("dropping table")
    for t in ("VARSPLIC_MASTER", "VARSPLIC_MATCH", "VARSPLIC_NEW"):
        orautils.drop_table(cur, "INTERPRO", t)

    logger.info("building talbe")
    cur.execute(
        """
        CREATE TABLE INTERPRO.VARSPLIC (
            PROTEIN_AC VARCHAR2(15) NOT NULL,
            VARIANT NUMBER(3) NOT NULL,
            LEN NUMBER(5) NOT NULL,
            METHOD_AC VARCHAR2(25) NOT NULL,
            DBCODE CHAR(1),
            POS_FROM NUMBER NOT NULL,
            POS_TO NUMBER NOT NULL,
            FRAGMENTS VARCHAR2(400),
            MODEL_AC VARCHAR2(255)
        ) NOLOGGING 
        """
    )

    cur.execute(
        """
        INSERT /*+ APPEND */ INTO INTERPRO.VARSPLIC
        SELECT 
            SUBSTR(X.AC, 1, INSTR(X.AC, '-') - 1),
            SUBSTR(X.AC, INSTR(X.AC, '-') + 1),
            P.LEN,
            M.METHOD_AC,
            I2D.DBCODE,
            M.SEQ_START,
            M.SEQ_END,
            M.FRAGMENTS,
            M.MODEL_AC
        FROM UNIPARC.XREF X
        INNER JOIN UNIPARC.PROTEIN P
          ON X.UPI = P.UPI
        INNER JOIN IPRSCAN.MV_IPRSCAN M
          ON X.UPI = M.UPI
        INNER JOIN INTERPRO.IPRSCAN2DBCODE I2D
          ON M.ANALYSIS_ID = I2D.IPRSCAN_SIG_LIB_REL_ID
        WHERE X.DBID IN (24, 25)
        AND X.DELETED = 'N'
        AND I2D.DBCODE NOT IN ('g', 'j', 'n', 'q', 's', 'v', 'x')
        """
    )
    con.commit()

    logger.info("indexing table")
    cur.execute(
        """
        CREATE INDEX I_VARSPLIC$P
        ON INTERPRO.VARSPLIC (PROTEIN_AC) NOLOGGING
        """
    )

    cur.execute(
        """
        CREATE INDEX I_VARSPLIC$M
        ON INTERPRO.VARSPLIC (METHOD_AC) NOLOGGING
        """
    )

    cur.close()
    con.close()
