# -*- coding: utf-8 -*-

import os
import time
from concurrent.futures import ThreadPoolExecutor
from typing import Callable, List, Optional, Tuple

import cx_Oracle

from .. import logger, orautils


PREFIX = "MV_"
SUFFIX = "_TMP"


def import_matches(user: str, dsn: str, max_workers: int=0,
                   checkpoint: Optional[str]=None, force: bool=False):
    if checkpoint:
        open(checkpoint, "w").close()
        os.chmod(checkpoint, 0o775)
        while os.path.isfile(checkpoint):
            time.sleep(600)

    url = orautils.make_connect_string(user, dsn)
    analyses = _get_analyses(url, datatype="matches", force=force)

    if not isinstance(max_workers, int) or max_workers < 1:
        max_workers = len(analyses)

    functions = {
        "ipm_cdd_match": (_insert_cdd_matches, "CDD"),
        "ipm_coils_match": (_insert_coils, "COILS"),
        "ipm_gene3d_match": (_insert_gene3d, "GENE3D"),
        "ipm_hamap_match": (_insert_hamap, "HAMAP"),
        "ipm_mobidblite_match": (_insert_mobidblite, "MOBIDBLITE"),
        "ipm_panther_match": (_insert_panther, "PANTHER"),
        "ipm_pfam_match": (_insert_pfam, "PFAM"),
        "ipm_phobius_match": (_insert_phobius, "PHOBIUS"),
        "ipm_pirsf_match": (_insert_pirsf, "PIRSF"),
        "ipm_prints_match": (_insert_prints, "PRINTS"),
        # "ipm_prodom_match": (_insert_prodom, "PRODOM"),
        "ipm_prosite_patterns_match": (_insert_prosite_patterns, "PROSITE_PATTERNS"),
        "ipm_prosite_profiles_match": (_insert_prosite_profiles, "PROSITE_PROFILES"),
        "ipm_sfld_match": (_insert_sfld_matches, "SFLD"),
        "ipm_smart_match": (_insert_smart, "SMART"),
        "ipm_superfamily_match": (_insert_superfamily, "SUPERFAMILY"),
        "ipm_tigrfam_match": (_insert_tigrfam, "TIGRFAM"),
        "ipm_tmhmm_match": (_insert_tmhmm, "TMHMM"),
    }
    signalp_partitions = {
        "signalp_euk": "SIGNALP_EUK",
        "signalp_gram_negative": "SIGNALP_GRAM_NEGATIVE",
        "signalp_gram_positive": "SIGNALP_GRAM_POSITIVE"
    }

    table2name = {}
    for analysis in analyses:
        table = analysis["match_table"].lower()
        if table == "ipm_signalp_match":
            table2name[table] = "SIGNALP"
        else:
            table2name[table] = analysis["full_name"]

    logger.info("updating MV_IPRSCAN")
    num_errors = 0
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        fs = {}
        pending = set()
        running = set()
        done = set()
        ignored = set()
        timestamp = time.time()

        while True:
            signalp_actions = []
            for analysis in analyses:
                _id = analysis["id"]
                table = analysis["match_table"].lower()

                if table in running or table in done or table in ignored:
                    continue
                elif table not in functions and table != "ipm_signalp_match":
                    logger.warning(f"ignored analysis: {analysis['full_name']}")
                    ignored.add(table)
                    continue
                elif analysis["ready"]:
                    try:
                        pending.remove(table)
                    except KeyError:
                        pass

                    if table == "ipm_signalp_match":
                        # SignalP has one source table, but three analyses
                        partition = signalp_partitions[analysis["name"]]
                        signalp_actions.append((_id, partition))
                    else:
                        fn, partition = functions[table]
                        f = executor.submit(_import_member_db, url, "IPRSCAN",
                                            table, "MV_IPRSCAN", partition,
                                            _id, fn)
                        fs[f] = table
                        running.add(table)
                        logger.info(f"  {analysis['full_name']:<40}ready")
                else:
                    pending.add(table)

            if len(signalp_actions) == len(signalp_partitions):
                f = executor.submit(_import_signalp, url, "IPRSCAN",
                                    "ipm_signalp_match", "MV_IPRSCAN",
                                    signalp_actions)
                fs[f] = "ipm_signalp_match"
                running.add("ipm_signalp_match")
                logger.info(f"  {'SIGNALP':<40}ready")

            _fs = {}
            for f in fs:
                table = fs[f]
                if f.done():
                    running.remove(table)
                    done.add(table)
                    exc = f.exception()
                    name = table2name[table]
                    if exc is None:
                        logger.info(f"  {name:<40}imported")
                    else:
                        exc_name = exc.__class__.__name__
                        logger.error(f"  {name:<40}failed ({exc_name}: {exc})")
                        num_errors += 1

                    timestamp = time.time()
                else:
                    _fs[f] = table

            fs = _fs

            if not pending and not running:
                break

            time.sleep(600)
            analyses = _get_analyses(url, datatype="matches", force=force)

            # if time.time() - timestamp > 6 * 3600:
            #     timestamp = time.time()
            #
            #     if pending:
            #         logger.info(f"  {len(pending)} pending: "
            #                     f"{', '.join(pending)}")
            #
            #     if running:
            #         logger.info(f"  {len(running)} running: "
            #                     f"{', '.join(running)}")

    if num_errors:
        raise RuntimeError("{} analyses failed".format(num_errors))

    logger.info("indexing MV_IPRSCAN")
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()

    for idx in orautils.get_indices(cur, "IPRSCAN", "MV_IPRSCAN"):
        # cur.execute("ALTER INDEX {} REBUILD NOLOGGING".format(idx))
        orautils.drop_index(cur, "IPRSCAN", idx["name"])

    orautils.catch_temp_error(cur, """
        CREATE INDEX I_MV_IPRSCAN$UPI
        ON IPRSCAN.MV_IPRSCAN (UPI) NOLOGGING
        """)

    orautils.gather_stats(cur, "IPRSCAN", "MV_IPRSCAN")
    cur.close()
    con.close()

    logger.info("MV_IPRSCAN is ready")


def import_sites(user: str, dsn: str, max_workers: int=0,
                 checkpoint: Optional[str]=None, force: bool=False):
    if checkpoint:
        open(checkpoint, "w").close()
        os.chmod(checkpoint, 0o775)
        while os.path.isfile(checkpoint):
            time.sleep(600)

    url = orautils.make_connect_string(user, dsn)
    analyses = _get_analyses(url, datatype="sites", force=force)

    if not isinstance(max_workers, int) or max_workers < 1:
        max_workers = len(analyses)

    partitions = {
        "ipm_cdd_site": "CDD",
        "ipm_sfld_site": "SFLD"
    }

    logger.info("updating SITE")
    num_errors = 0
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        fs = {}
        pending = set()
        running = set()
        done = set()
        while True:
            for analysis in analyses:
                _id = analysis["id"]
                name = analysis["name"]
                full_name = analysis["full_name"]
                table = analysis["site_table"]

                if not table:
                    continue
                elif full_name in running or full_name in done:
                    continue
                elif not analysis["ready"]:
                    pending.add(full_name)
                    continue
                else:
                    try:
                        pending.remove(full_name)
                    except KeyError:
                        pass

                    partition = partitions[table.lower()]
                    f = executor.submit(_import_member_db, url, "IPRSCAN",
                                        table, "SITE", partition, _id,
                                        _insert_sites)
                    fs[f] = full_name
                    running.add(full_name)
                    logger.info(f"  {full_name:<40}ready")

            _fs = {}
            for f in fs:
                full_name = fs[f]
                if f.done():
                    running.remove(full_name)
                    done.add(full_name)
                    exc = f.exception()
                    if exc is None:
                        logger.info(f"  {full_name:<40}imported")
                    else:
                        exc_name = exc.__class__.__name__
                        logger.error(f"  {full_name:<40}failed ({exc_name}: {exc})")
                        num_errors += 1
                else:
                    _fs[f] = full_name
            fs = _fs

            if not pending and not running:
                break

            time.sleep(600)
            analyses = _get_analyses(url, datatype="sites", force=force)

    if num_errors:
        raise RuntimeError("{} analyses failed".format(num_errors))

    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()

    for idx in orautils.get_indices(cur, "IPRSCAN", "SITE"):
        # cur.execute("ALTER INDEX {} REBUILD NOLOGGING".format(idx))
        orautils.drop_index(cur, "IPRSCAN", idx["name"])

    logger.info("indexing SITE")
    cur.execute(
        """
        CREATE INDEX I_SITE$UPI
        ON IPRSCAN.SITE (UPI) NOLOGGING
        """
    )
    orautils.gather_stats(cur, "IPRSCAN", "SITE")
    cur.close()
    con.close()

    logger.info("SITE is ready")


def _get_max_persisted_job(cur: cx_Oracle.Cursor, analysis_id: int,
                           max_upi: str, persist_flag: int) -> Optional[str]:
    cur.execute(
        """
        SELECT SUM(CNT)
        FROM (
            SELECT COUNT(*) AS CNT
            FROM IPRSCAN.IPM_RUNNING_JOBS@ISPRO
            WHERE ANALYSIS_ID = :analysisid
              AND JOB_START <= :maxupi
            UNION ALL
            SELECT COUNT(*) AS CNT
            FROM IPRSCAN.IPM_COMPLETED_JOBS@ISPRO
            WHERE ANALYSIS_ID = :analysisid
              AND JOB_START <= :maxupi AND PERSISTED < :persisted
            UNION ALL
            SELECT COUNT(*) AS CNT
            FROM IPRSCAN.IPM_PERSISTED_JOBS@ISPRO
            WHERE ANALYSIS_ID = :analysisid
              AND JOB_START <= :maxupi AND PERSISTED < :persisted
        )
        """,
        dict(analysisid=analysis_id, persisted=persist_flag, maxupi=max_upi)
    )
    if cur.fetchone()[0]:
        return None

    cur.execute(
        """
        SELECT MAX(JOB_END)
        FROM IPRSCAN.IPM_PERSISTED_JOBS@ISPRO
        WHERE ANALYSIS_ID = :1
        """, (analysis_id, )
    )
    row = cur.fetchone()
    return row[0] if row else None


def _get_analyses(url: str, datatype: str="matches", force: bool=False) -> List[dict]:
    if datatype not in ("matches", "sites"):
        raise ValueError("invalid value for 'datatype' argument "
                         "(expects 'matches' or 'sites')")

    con = cx_Oracle.connect(url)
    cur = con.cursor()
    max_upi = _get_max_upi(cur, "UNIPARC", "PROTEIN", "UAREAD")
    cur.execute(
        """
        SELECT *
        FROM (
          SELECT
            A.ANALYSIS_ID,
            A.ANALYSIS_NAME,
            A.ANALYSIS_TYPE,
            B.MATCH_TABLE,
            B.SITE_TABLE,
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

    match_tables = set()
    analyses = []
    for row in cur:
        """
        SignalP has EUK/Gram+/Gram- matches in the same tables so we need
        several times the same time for different analyses.

        But we shouldn't accept any other duplicates
        (ACTIVE might not always be up-to-date)
        """
        match_table = row[3]
        if match_table in match_tables and match_table != "ipm_signalp_match":
            continue

        match_tables.add(match_table)
        analyses.append({
            "id": row[0],
            "full_name": row[1],
            "name": row[2],
            "match_table": match_table,
            "site_table": row[4],
            "upi": None,
            "ready": False
        })

    for e in analyses:
        """
        CDD/SFLD
            - PERSISTED=1 when matches are ready
            - PERSISTED=2 when matches and sites are ready

        Others:
            - PERSISTED=2 when matches are ready
        """
        if datatype == "matches":
            key = "match_table"
            flag = 1 if e["name"] in ("cdd", "sfld") else 2
        else:
            key = "site_table"
            flag = 2

        if e[key] is not None:
            e["upi"] = _get_max_persisted_job(cur, e["id"], max_upi, flag)
            e["ready"] = e["upi"] and e["upi"] >= max_upi

        if force:
            e["ready"] = True

    cur.close()
    con.close()

    return analyses


def check_ispro(url: str, datatype: str):
    analyses = _get_analyses(url, datatype)
    for e in analyses:
        e["status"] = "ready" if e["ready"] else "not ready"

        if not e["site_table"]:
            e["site_table"] = ""

        if not e["upi"]:
            e["upi"] = ""

        print("{id:<5}{full_name:<30}{match_table:<30}"
              "{site_table:<30}{upi:<20}{status}".format(**e))


def _get_max_upi(cur: cx_Oracle.Cursor, owner: str, table: str,
                 dblink: Optional[str]=None) -> Optional[str]:
    try:
        cur.execute(
            """
            SELECT MAX(UPI)
            FROM {}.{}{}
            """.format(owner, table, '@' + dblink if dblink else '')
        )
    except cx_Oracle.DatabaseError as exc:
        error, = exc.args
        if error.code == 942:
            # ORA-00942 (table or view does not exist)
            return None
        else:
            raise exc
    else:
        return cur.fetchone()[0]


def _import_signalp(url: str, owner: str, table_src: str, table_dst: str,
                    actions: List[Tuple[int, str]]):
    table_stg = PREFIX + table_src
    table_tmp = table_src + SUFFIX

    con = cx_Oracle.connect(url)
    cur = con.cursor()
    upi = _get_max_upi(cur, owner, table_stg)
    max_upi = _get_max_upi(cur, "UNIPARC", "PROTEIN", "UAREAD")
    if not upi or upi < max_upi:
        # Not the same UPI: import table
        orautils.drop_mview(cur, owner, table_stg)
        orautils.drop_table(cur, owner, table_stg, purge=True)
        cur.execute(
            """
            CREATE TABLE {0}.{1} NOLOGGING
            AS
            SELECT *
            FROM {0}.{2}@ISPRO
            """.format(owner, table_stg, table_src)
        )
        cur.execute(
            """
            CREATE INDEX {0}.{1}$ID ON {0}.{2} (ANALYSIS_ID) NOLOGGING
            """.format(owner, table_src, table_stg)
        )
        cur.execute(
            """
            CREATE INDEX {0}.{1}$UPI ON {0}.{2} (UPI) NOLOGGING
            """.format(owner, table_src, table_stg)
        )

    # Create temporary table for partition exchange
    orautils.drop_table(cur, owner, table_tmp, purge=True)
    cur.execute(
        """
        CREATE TABLE {0}.{1} NOLOGGING
        AS
        SELECT *
        FROM {0}.{2}
        WHERE 1=0
        """.format(owner, table_tmp, table_dst)
    )
    for analysis_id, partition in actions:
        orautils.truncate_table(cur, owner, table_tmp)
        cur.execute(
            """
            INSERT /*+ APPEND */ INTO {0}.{1}
            SELECT
              ANALYSIS_ID, UPI, METHOD_AC, RELNO_MAJOR, RELNO_MINOR, SEQ_START,
              SEQ_END, 0, 0, 0, NULL, SEQSCORE, SEQSCORE, 0, 0, 0, 0, MODEL_AC,
              NULL, FRAGMENTS
            FROM {0}.{2}
            WHERE ANALYSIS_ID = :1
            """.format(owner, table_tmp, table_stg),
            (analysis_id,)
        )
        con.commit()
        orautils.exchange_partition(cur, owner, table_tmp, table_dst, partition, stats=False)

    orautils.drop_table(cur, owner, table_tmp, purge=True)
    cur.close()
    con.close()


def _import_member_db(url: str, owner: str, table_src: str, table_dst: str,
                      partition: str, analysis_id: int, func: Callable):
    table_stg = PREFIX + table_src
    table_tmp = table_src + SUFFIX

    con = cx_Oracle.connect(url)
    cur = con.cursor()
    upi = _get_max_upi(cur, owner, "{} partition ({})".format(table_dst, partition))
    max_upi = _get_max_upi(cur, "UNIPARC", "PROTEIN", "UAREAD")
    if upi >= max_upi:
        # Nothing to do here
        cur.close()
        con.close()
        return

    upi = _get_max_upi(cur, owner, table_stg)
    if not upi or upi < max_upi:
        # Not the same UPI: import table
        orautils.drop_mview(cur, owner, table_stg)
        orautils.drop_table(cur, owner, table_stg, purge=True)
        cur.execute(
            """
            CREATE TABLE {0}.{1} NOLOGGING
            AS
            SELECT *
            FROM {0}.{2}@ISPRO
            """.format(owner, table_stg, table_src)
        )
        cur.execute(
            """
            CREATE INDEX {0}.{1}$ID ON {0}.{2} (ANALYSIS_ID) NOLOGGING
            """.format(owner, table_src, table_stg)
        )
        cur.execute(
            """
            CREATE INDEX {0}.{1}$UPI ON {0}.{2} (UPI) NOLOGGING
            """.format(owner, table_src, table_stg)
        )

    # Create temporary table for partition exchange
    orautils.drop_table(cur, owner, table_tmp, purge=True)
    cur.execute(
        """
        CREATE TABLE {0}.{1} NOLOGGING
        AS
        SELECT *
        FROM {0}.{2}
        WHERE 1=0
        """.format(owner, table_tmp, table_dst)
    )

    func(cur, owner, table_stg, table_tmp, analysis_id)
    con.commit()  # important to commit here as func() does not

    # Exchange partition and drop temporary table
    orautils.exchange_partition(cur, owner, table_tmp, table_dst, partition, stats=False)
    orautils.drop_table(cur, owner, table_tmp, purge=True)
    cur.close()
    con.close()


def _insert_cdd_matches(cur: cx_Oracle.Cursor, owner: str, table_src: str,
                        table_dst: str, analysis_id: int):
    cur.execute(
        """
        INSERT /*+ APPEND */ INTO {0}.{1}
        SELECT
          ANALYSIS_ID, UPI, METHOD_AC, RELNO_MAJOR, RELNO_MINOR, SEQ_START,
          SEQ_END, 0, 0, 0, NULL, SEQSCORE, SEQSCORE, SEQEVALUE, SEQEVALUE,
          0, 0, MODEL_AC, NULL, FRAGMENTS
        FROM {0}.{2}
        WHERE ANALYSIS_ID = :1
        """.format(owner, table_dst, table_src),
        (analysis_id,)
    )


def _insert_sites(cur: cx_Oracle.Cursor, owner: str, table_src: str,
                  table_dst: str, analysis_id: int):
    cur.execute(
        """
        INSERT /*+ APPEND */ INTO {0}.{1}
        SELECT
            UPI, ANALYSIS_ID, METHOD_AC, LOC_START, LOC_END, NUM_SITES,
            RESIDUE, RES_START, RES_END, DESCRIPTION
        FROM {0}.{2}
        WHERE ANALYSIS_ID = :1
        """.format(owner, table_dst, table_src),
        (analysis_id,)
    )


def _insert_coils(cur: cx_Oracle.Cursor, owner: str, table_src: str,
                  table_dst: str, analysis_id: int):
    cur.execute(
        """
        INSERT /*+ APPEND */ INTO {0}.{1}
        SELECT
          ANALYSIS_ID, UPI, METHOD_AC, RELNO_MAJOR, RELNO_MINOR, SEQ_START,
          SEQ_END, 0, 0, 0, NULL, 0, 0, 0, 0, 0, 0, MODEL_AC, NULL, FRAGMENTS
        FROM {0}.{2}
        WHERE ANALYSIS_ID = :1
        """.format(owner, table_dst, table_src),
        (analysis_id,)
    )


def _insert_gene3d(cur: cx_Oracle.Cursor, owner: str, table_src: str,
                   table_dst: str, analysis_id: int):
    cur.execute(
        """
        INSERT /*+ APPEND */ INTO {0}.{1}
        SELECT
          ANALYSIS_ID, UPI, METHOD_AC, RELNO_MAJOR, RELNO_MINOR, SEQ_START,
          SEQ_END, HMM_START, HMM_END, HMM_LENGTH, HMM_BOUNDS, SCORE,
          SEQSCORE, EVALUE, SEQEVALUE, ENV_START, ENV_END, MODEL_AC,
          NULL, FRAGMENTS
        FROM {0}.{2}
        WHERE ANALYSIS_ID = :1
        """.format(owner, table_dst, table_src),
        (analysis_id,)
    )


def _insert_hamap(cur: cx_Oracle.Cursor, owner: str, table_src: str,
                  table_dst: str, analysis_id: int):
    cur.execute(
        """
        INSERT /*+ APPEND */ INTO {0}.{1}
        SELECT
          ANALYSIS_ID, UPI, METHOD_AC, SUBSTR(RELNO_MAJOR, 1, 4),
          SUBSTR(RELNO_MAJOR, 6, 7), SEQ_START, SEQ_END, 0, 0, 0, NULL , 0,
          SEQSCORE, 0, 0, 0, 0, MODEL_AC, NULL, FRAGMENTS
        FROM {0}.{2}
        WHERE ANALYSIS_ID = :1
        """.format(owner, table_dst, table_src),
        (analysis_id,)
    )


def _insert_mobidblite(cur: cx_Oracle.Cursor, owner: str, table_src: str,
                       table_dst: str, analysis_id: int):
    cur.execute(
        """
        INSERT /*+ APPEND */ INTO {0}.{1}
        SELECT
          ANALYSIS_ID, UPI, METHOD_AC, RELNO_MAJOR, RELNO_MINOR, SEQ_START,
          SEQ_END, 0, 0, 0, NULL , 0, 0, 0, 0, 0, 0, MODEL_AC, SEQ_FEATURE,
          FRAGMENTS
        FROM {0}.{2}
        WHERE ANALYSIS_ID = :1
        """.format(owner, table_dst, table_src),
        (analysis_id,)
    )


def _insert_panther(cur: cx_Oracle.Cursor, owner: str, table_src: str,
                    table_dst: str, analysis_id: int):
    cur.execute(
        """
        INSERT /*+ APPEND */ INTO {0}.{1}
        SELECT
          ANALYSIS_ID, UPI, METHOD_AC, RELNO_MAJOR, RELNO_MINOR, SEQ_START,
          SEQ_END, HMM_START, HMM_END, HMM_LENGTH, HMM_BOUNDS, SEQSCORE,
          SEQSCORE, SEQEVALUE, SEQEVALUE, ENV_START, ENV_END, MODEL_AC,
          NULL, FRAGMENTS
        FROM {0}.{2}
        WHERE ANALYSIS_ID = :1
        """.format(owner, table_dst, table_src),
        (analysis_id,)
    )


def _insert_pfam(cur: cx_Oracle.Cursor, owner: str, table_src: str,
                 table_dst: str, analysis_id: int):
    cur.execute(
        """
        INSERT /*+ APPEND */ INTO {0}.{1}
        SELECT
          ANALYSIS_ID, UPI, METHOD_AC, RELNO_MAJOR, RELNO_MINOR, SEQ_START,
          SEQ_END, HMM_START, HMM_END, HMM_LENGTH, HMM_BOUNDS, SCORE,
          SEQSCORE, EVALUE, SEQEVALUE, ENV_START, ENV_END, MODEL_AC,
          NULL, FRAGMENTS
        FROM {0}.{2}
        WHERE ANALYSIS_ID = :1
        """.format(owner, table_dst, table_src),
        (analysis_id,)
    )


def _insert_phobius(cur: cx_Oracle.Cursor, owner: str, table_src: str,
                    table_dst: str, analysis_id: int):
    cur.execute(
        """
        INSERT /*+ APPEND */ INTO {0}.{1}
        SELECT
          ANALYSIS_ID, UPI, METHOD_AC, RELNO_MAJOR, RELNO_MINOR, SEQ_START,
          SEQ_END, 0, 0, 0, NULL, 0, 0, 0, 0, 0, 0, MODEL_AC, NULL, FRAGMENTS
        FROM {0}.{2}
        WHERE ANALYSIS_ID = :1
        """.format(owner, table_dst, table_src),
        (analysis_id,)
    )


def _insert_pirsf(cur: cx_Oracle.Cursor, owner: str, table_src: str,
                  table_dst: str, analysis_id: int):
    cur.execute(
        """
        INSERT /*+ APPEND */ INTO {0}.{1}
        SELECT
          ANALYSIS_ID, UPI, METHOD_AC, RELNO_MAJOR, RELNO_MINOR, SEQ_START,
          SEQ_END, HMM_START, HMM_END, HMM_LENGTH, HMM_BOUNDS, SCORE,
          SEQSCORE, EVALUE, SEQEVALUE, ENV_START, ENV_END, MODEL_AC,
          NULL, FRAGMENTS
        FROM {0}.{2}
        WHERE ANALYSIS_ID = :1
        """.format(owner, table_dst, table_src),
        (analysis_id,)
    )


def _insert_prints(cur: cx_Oracle.Cursor, owner: str, table_src: str,
                   table_dst: str, analysis_id: int):
    cur.execute(
        """
        INSERT /*+ APPEND */ INTO {0}.{1}
        SELECT
          ANALYSIS_ID, UPI, METHOD_AC, RELNO_MAJOR, RELNO_MINOR, SEQ_START,
          SEQ_END, 0, 0, MOTIF_NUMBER, NULL, 0, SEQSCORE, PVALUE, SEQEVALUE,
          0, 0, MODEL_AC, GRAPHSCAN, FRAGMENTS
        FROM {0}.{2}
        WHERE ANALYSIS_ID = :1
        """.format(owner, table_dst, table_src),
        (analysis_id,)
    )


def _insert_prodom(cur: cx_Oracle.Cursor, owner: str, table_src: str,
                   table_dst: str, analysis_id: int):
    cur.execute(
        """
        INSERT /*+ APPEND */ INTO {0}.{1}
        SELECT
          ANALYSIS_ID, UPI, METHOD_AC, RELNO_MAJOR, RELNO_MINOR, SEQ_START,
          SEQ_END, 0, 0, 0, NULL, SEQSCORE, SEQSCORE, SEQEVALUE, SEQEVALUE,
          0, 0, MODEL_AC, NULL, FRAGMENTS
        FROM {0}.{2}
        WHERE ANALYSIS_ID = :1
        """.format(owner, table_dst, table_src),
        (analysis_id,)
    )


def _insert_prosite_patterns(cur: cx_Oracle.Cursor, owner: str, table_src: str,
                             table_dst: str, analysis_id: int):
    cur.execute(
        """
        INSERT /*+ APPEND */ INTO {0}.{1}
        SELECT
          ANALYSIS_ID, UPI, METHOD_AC, SUBSTR(RELNO_MAJOR, 1, 4),
          SUBSTR(RELNO_MAJOR, 6, 7), SEQ_START,
          SEQ_END, 0, 0, 0, LOCATION_LEVEL, 0, 0, 0, 0, 0, 0, MODEL_AC,
          ALIGNMENT, FRAGMENTS
        FROM {0}.{2}
        WHERE ANALYSIS_ID = :1
        """.format(owner, table_dst, table_src),
        (analysis_id,)
    )


def _insert_prosite_profiles(cur: cx_Oracle.Cursor, owner: str, table_src: str,
                             table_dst: str, analysis_id: int):
    cur.execute(
        """
        INSERT /*+ APPEND */ INTO {0}.{1}
        SELECT
          ANALYSIS_ID, UPI, METHOD_AC, SUBSTR(RELNO_MAJOR, 1, 4),
          SUBSTR(RELNO_MAJOR, 6, 7), SEQ_START,
          SEQ_END, 0, 0, 0, NULL, 0, SEQSCORE, 0, 0, 0, 0, MODEL_AC,
          ALIGNMENT, FRAGMENTS
        FROM {0}.{2}
        WHERE ANALYSIS_ID = :1
        """.format(owner, table_dst, table_src),
        (analysis_id,)
    )


def _insert_sfld_matches(cur: cx_Oracle.Cursor, owner: str, table_src: str,
                         table_dst: str, analysis_id: int):
    cur.execute(
        """
        INSERT /*+ APPEND */ INTO {0}.{1}
        SELECT
          ANALYSIS_ID, UPI, METHOD_AC, RELNO_MAJOR, RELNO_MINOR, SEQ_START,
          SEQ_END, HMM_START, HMM_END, HMM_LENGTH, HMM_BOUNDS, SCORE,
          SEQSCORE, EVALUE, SEQEVALUE, ENV_START, ENV_END, MODEL_AC, NULL,
          FRAGMENTS
        FROM {0}.{2}
        WHERE ANALYSIS_ID = :1
        """.format(owner, table_dst, table_src),
        (analysis_id,)
    )


def _insert_smart(cur: cx_Oracle.Cursor, owner: str, table_src: str,
                  table_dst: str, analysis_id: int):
    cur.execute(
        """
        INSERT /*+ APPEND */ INTO {0}.{1}
        SELECT
          ANALYSIS_ID, UPI, METHOD_AC, RELNO_MAJOR, RELNO_MINOR, SEQ_START,
          SEQ_END, HMM_START, HMM_END, HMM_LENGTH, HMM_BOUNDS, SCORE,
          SEQSCORE, EVALUE, SEQEVALUE, 0, 0, MODEL_AC, NULL, FRAGMENTS
        FROM {0}.{2}
        WHERE ANALYSIS_ID = :1
        """.format(owner, table_dst, table_src),
        (analysis_id,)
    )


def _insert_superfamily(cur: cx_Oracle.Cursor, owner: str, table_src: str,
                        table_dst: str, analysis_id: int):
    cur.execute(
        """
        INSERT /*+ APPEND */ INTO {0}.{1}
        SELECT
          ANALYSIS_ID, UPI, METHOD_AC, RELNO_MAJOR, RELNO_MINOR, SEQ_START,
          SEQ_END, 0, 0, HMM_LENGTH, NULL, 0, 0, SEQEVALUE, SEQEVALUE, 0, 0,
          MODEL_AC, NULL, FRAGMENTS
        FROM {0}.{2}
        WHERE ANALYSIS_ID = :1
        """.format(owner, table_dst, table_src),
        (analysis_id,)
    )


def _insert_tigrfam(cur: cx_Oracle.Cursor, owner: str, table_src: str,
                    table_dst: str, analysis_id: int):
    cur.execute(
        """
        INSERT /*+ APPEND */ INTO {0}.{1}
        SELECT
          ANALYSIS_ID, UPI, METHOD_AC, RELNO_MAJOR, RELNO_MINOR, SEQ_START,
          SEQ_END, HMM_START, HMM_END, HMM_LENGTH, HMM_BOUNDS, SCORE,
          SEQSCORE, EVALUE, SEQEVALUE, ENV_START, ENV_END, MODEL_AC,
          NULL, FRAGMENTS
        FROM {0}.{2}
        WHERE ANALYSIS_ID = :1
        """.format(owner, table_dst, table_src),
        (analysis_id,)
    )


def _insert_tmhmm(cur: cx_Oracle.Cursor, owner: str, table_src: str,
                  table_dst: str, analysis_id: int):
    cur.execute(
        """
        INSERT /*+ APPEND */ INTO {0}.{1}
        SELECT
          ANALYSIS_ID, UPI, METHOD_AC, RELNO_MAJOR, RELNO_MINOR, SEQ_START,
          SEQ_END, 0, 0, 0, NULL, SEQSCORE, SEQSCORE, 0, 0, 0, 0, MODEL_AC,
          NULL, FRAGMENTS
        FROM {0}.{2}
        WHERE ANALYSIS_ID = :1
        """.format(owner, table_dst, table_src),
        (analysis_id,)
    )
