# -*- coding: utf-8 -*-

import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Callable, Collection, List, Optional, Tuple

import cx_Oracle

from .. import logger, orautils


def import_matches(user: str, dsn: str, **kwargs):
    max_workers = kwargs.pop("max_workers", 0)

    url = orautils.make_connect_string(user, dsn)
    analyses = _check_ispro(url, **kwargs)

    if not isinstance(max_workers, int) or max_workers < 1:
        max_workers = len(analyses)

    logger.info("updating MV_IPRSCAN")
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        fs = {}
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
        signalp_actions = []

        for analysis in analyses:
            _id = analysis["id"]
            name = analysis["name"]
            full_name = analysis["full_name"]
            table = analysis["match_table"].lower()

            if table == "ipm_signalp_match":
                # SignalP has one source table, but three dbcodes
                partition = signalp_partitions[name]
                signalp_actions.append((_id, partition))
            elif table in functions:
                fn, partition = functions[table]
                f = executor.submit(_import_member_db, url, "IPRSCAN", table,
                                    "MV_IPRSCAN", partition, _id, fn)
                fs[f] = full_name

        if signalp_actions:
            f = executor.submit(_import_signalp, url, "IPRSCAN",
                                "ipm_signalp_match", "MV_IPRSCAN",
                                signalp_actions)
            fs[f] = "SIGNALP"

        num_errors = 0
        for f in as_completed(fs):
            name = fs[f]
            exc = f.exception()
            if exc is None:
                logger.info("\t{} is ready".format(name))
            else:
                exc_name = exc.__class__.__name__
                logger.error("\t{} failed ({}: {})".format(name, exc_name, exc))
                num_errors += 1

    if num_errors:
        raise RuntimeError("{} analyses failed".format(num_errors))

    logger.info("indexing MV_IPRSCAN")
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()

    for idx in orautils.get_indices(cur, "IPRSCAN", "MV_IPRSCAN"):
        # cur.execute("ALTER INDEX {} REBUILD NOLOGGING".format(idx))
        orautils.drop_index(cur, "IPRSCAN", idx["name"])

    cur.execute(
        """
        CREATE INDEX I_MV_IPRSCAN$UPI
        ON IPRSCAN.MV_IPRSCAN (UPI) NOLOGGING
        """
    )

    cur.close()
    con.close()

    logger.info("MV_IPRSCAN is ready")


def import_sites(user: str, dsn: str, **kwargs):
    max_workers = kwargs.pop("max_workers", 0)

    url = orautils.make_connect_string(user, dsn)
    analyses = _check_ispro(url, **kwargs)

    if not isinstance(max_workers, int) or max_workers < 1:
        max_workers = len(analyses)

    logger.info("updating SITE")
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        fs = {}
        partitions = {
            "ipm_cdd_site": "CDD",
            "ipm_sfld_site": "SFLD"
        }
        for analysis in analyses:
            _id = analysis["id"]
            full_name = analysis["full_name"]
            table = analysis["site_table"]

            if not table:
                continue

            partition = partitions[table.lower()]
            f = executor.submit(_import_member_db, url, "IPRSCAN", table,
                                "SITE", partition, _id, _insert_sites)
            fs[f] = full_name

        num_errors = 0
        for f in as_completed(fs):
            name = fs[f]
            exc = f.exception()
            if exc is None:
                logger.info("\t{} is ready".format(name))
            else:
                exc_name = exc.__class__.__name__
                logger.error("\t{} failed ({}: {})".format(name, exc_name, exc))
                num_errors += 1

    if num_errors:
        raise RuntimeError("{} analyses failed".format(num_errors))

    logger.info("indexing SITE")
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()

    for idx in orautils.get_indices(cur, "IPRSCAN", "SITE"):
        # cur.execute("ALTER INDEX {} REBUILD NOLOGGING".format(idx))
        orautils.drop_index(cur, "IPRSCAN", idx["name"])

    cur.execute(
        """
        CREATE INDEX I_SITE$UPI
        ON IPRSCAN.SITE (UPI) NOLOGGING
        """
    )

    cur.close()
    con.close()

    logger.info("SITE is ready")


def _get_max_persisted_job(cur: cx_Oracle.Cursor, analysis_id: int) -> Optional[str]:
    cur.execute(
        """
        SELECT COUNT(*)
        FROM IPRSCAN.IPM_RUNNING_JOBS@ISPRO
        WHERE ANALYSIS_ID = :1
        """, (analysis_id,)
    )
    if cur.fetchone()[0]:
        return None

    cur.execute(
        """
        SELECT MAX(JOB_END)
        FROM IPRSCAN.IPM_PERSISTED_JOBS@ISPRO
        WHERE ANALYSIS_ID = :1
        AND PERSISTED > 0
        """, (analysis_id, )
    )
    row = cur.fetchone()
    return row[0] if row else None

    # upis = []
    # cur.execute(
    #     """
    #     SELECT MAX(HWM_SUBMITTED)
    #     FROM IPRSCAN.HWM@ISPRO
    #     WHERE ANALYSIS_ID = :1
    #     """, (analysis_id,)
    # )
    # row = cur.fetchone()
    # if row:
    #     upis.append(row[0])
    #
    # cur.execute(
    #     """
    #     SELECT MAX(JOB_END)
    #     FROM IPRSCAN.IPM_RUNNING_JOBS@ISPRO
    #     WHERE ANALYSIS_ID = :1
    #     """, (analysis_id,)
    # )
    # row = cur.fetchone()
    # if row:
    #     upis.append(row[0])
    #
    # cur.execute(
    #     """
    #     SELECT MIN(JOB_END)
    #     FROM (
    #         SELECT MIN(JOB_END) JOB_END
    #         FROM IPRSCAN.IPM_COMPLETED_JOBS@ISPRO
    #         WHERE ANALYSIS_ID = :analysisid
    #         AND PERSISTED = 0
    #         UNION ALL
    #         SELECT MAX(JOB_END) JOB_END
    #         FROM IPRSCAN.IPM_COMPLETED_JOBS@ISPRO
    #         WHERE ANALYSIS_ID = :analysisid
    #         AND PERSISTED = 1
    #     )
    #     """, dict(analysisid=analysis_id)
    # )
    # row = cur.fetchone()
    # if row:
    #     upis.append(row[0])
    #
    # cur.execute(
    #     """
    #     SELECT MIN(JOB_END)
    #     FROM (
    #         SELECT MIN(JOB_END) JOB_END
    #         FROM IPRSCAN.IPM_PERSISTED_JOBS@ISPRO
    #         WHERE ANALYSIS_ID = :analysisid
    #         AND PERSISTED = 0
    #         UNION ALL
    #         SELECT MAX(JOB_END) JOB_END
    #         FROM IPRSCAN.IPM_PERSISTED_JOBS@ISPRO
    #         WHERE ANALYSIS_ID = :analysisid
    #         AND PERSISTED = 1
    #     )
    #     """, dict(analysisid=analysis_id)
    # )
    # row = cur.fetchone()
    # if row:
    #     upis.append(row[0])
    #
    # try:
    #     return max([upi for upi in upis if upi is not None])
    # except ValueError:
    #     return None


def _get_analyses(url: str) -> List[dict]:
    con = cx_Oracle.connect(url)
    cur = con.cursor()
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
            "upi": None
        })

    for e in analyses:
        e["upi"] = _get_max_persisted_job(cur, e["id"])

    cur.close()
    con.close()

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
        analyses = _get_analyses(url)
        not_ready = 0
        for e in analyses:
            if exclude and e["name"] in exclude:
                status = "ready (excluded)"
            elif e["upi"] and e["upi"] >= max_upi_read:
                status = "ready"
            else:
                status = "not ready"
                not_ready += 1

            if e["site_table"]:
                _site_table = e["site_table"]
            else:
                _site_table = ""

            logger.debug("{id:<5}{full_name:<30}{match_table:<30}"
                         "{_site_table:<30}{upi:<20}"
                         "{status}".format(**e, status=status,
                                           _site_table=_site_table))

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
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    dst_upi = _get_max_upi(cur, owner, table_src)
    if not dst_upi or dst_upi != _get_max_upi(cur, owner, table_src, "ISPRO"):
        # Not the same UPI: import table
        orautils.drop_mview(cur, owner, table_src)
        orautils.drop_table(cur, owner, table_src)
        cur.execute(
            """
            CREATE TABLE {0}.{1} NOLOGGING
            AS
            SELECT *
            FROM {0}.{1}@ISPRO
            """.format(owner, table_src)
        )
        orautils.gather_stats(cur, owner, table_src)
        cur.execute(
            """
            CREATE INDEX {0}.I_{1}
            ON {0}.{1}(ANALYSIS_ID)
            NOLOGGING
            """.format(owner, table_src)
        )

    table_tmp = table_src + "_TMP"
    orautils.drop_table(cur, owner, table_tmp)
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
            """.format(owner, table_tmp, table_src),
            (analysis_id,)
        )
        con.commit()
        orautils.exchange_partition(cur, owner, table_tmp, table_dst, partition)

    orautils.drop_table(cur, owner, table_tmp)
    cur.close()
    con.close()


def _import_member_db(url: str, owner: str, table_src: str, table_dst: str,
                      partition: str, analysis_id: int, func: Callable):
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    dst_upi = _get_max_upi(cur, owner, table_src)
    if not dst_upi or dst_upi != _get_max_upi(cur, owner, table_src, "ISPRO"):
        # Not the same UPI: import table
        orautils.drop_mview(cur, owner, table_src)
        orautils.drop_table(cur, owner, table_src)
        cur.execute(
            """
            CREATE TABLE {0}.{1} NOLOGGING
            AS
            SELECT *
            FROM {0}.{1}@ISPRO
            """.format(owner, table_src)
        )
        orautils.gather_stats(cur, owner, table_src)
        cur.execute(
            """
            CREATE INDEX {0}.I_{1}
            ON {0}.{1}(ANALYSIS_ID)
            NOLOGGING
            """.format(owner,  table_src)
        )

    # Create temporary table for partition exchange
    table_tmp = table_src + "_TMP"
    orautils.drop_table(cur, owner, table_tmp)
    cur.execute(
        """
        CREATE TABLE {0}.{1} NOLOGGING
        AS
        SELECT *
        FROM {0}.{2}
        WHERE 1=0
        """.format(owner, table_tmp, table_dst)
    )

    func(cur, owner, table_src, table_tmp, analysis_id)
    con.commit()  # important to commit here as func() does not

    # Exchange partition and drop temporary table
    orautils.exchange_partition(cur, owner, table_tmp, table_dst, partition)
    orautils.drop_table(cur, owner, table_tmp)
    cur.close()
    con.close()


def _import_table(url: str, owner: str, name: str):
    index = "I_" + name
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.execute(
        """
        SELECT MAX(UPI)
        FROM {}.{}@ISPRO
        """.format(owner, name)
    )
    src_upi = cur.fetchone()[0]

    try:
        cur.execute(
            """
            SELECT MAX(UPI)
            FROM {}.{}
            """.format(owner, name)
        )
    except cx_Oracle.DatabaseError as exc:
        error, = exc.args
        if error.code == 942:
            # ORA-00942 (table or view does not exist)
            dst_upi = None
        else:
            raise exc
    else:
        dst_upi = cur.fetchone()[0]

    if src_upi != dst_upi:
        orautils.drop_mview(cur, owner, name)
        orautils.drop_table(cur, owner, name)
        orautils.drop_index(cur, owner, index)
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
            CREATE INDEX {0}.{1}
            ON {0}.{2}(ANALYSIS_ID)
            NOLOGGING
            """.format(owner, index, name)
        )
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
