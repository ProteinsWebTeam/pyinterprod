# -*- coding: utf-8 -*-

import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from typing import List, Sequence, Tuple

import cx_Oracle
from cx_Oracle import Cursor

from pyinterprod import logger
from pyinterprod.interpro.database import Database
from pyinterprod.utils import oracle

PREFIX = "MV_"

# Columns to select when inserting matches in MV_IPRSCAN
# Keys are partitions names
MATCH_SELECT = {
    "cdd": [
        'ANALYSIS_ID', 'UPI', 'METHOD_AC', 'RELNO_MAJOR', 'RELNO_MINOR',
        'SEQ_START', 'SEQ_END', '0', '0', '0',
        'NULL', 'SEQSCORE',  'SEQSCORE', 'SEQEVALUE', 'SEQEVALUE',
        '0', '0', 'MODEL_AC', 'NULL', 'FRAGMENTS'
    ],
    "coils": [
        'ANALYSIS_ID', 'UPI', 'METHOD_AC', 'RELNO_MAJOR', 'RELNO_MINOR',
        'SEQ_START', 'SEQ_END', '0', '0', '0',
        'NULL', '0', '0', '0', '0',
        '0', '0', 'MODEL_AC', 'NULL', 'FRAGMENTS'
    ],
    "gene3d": [
        'ANALYSIS_ID', 'UPI', 'METHOD_AC', 'RELNO_MAJOR', 'RELNO_MINOR',
        'SEQ_START', 'SEQ_END', 'HMM_START', 'HMM_END', 'HMM_LENGTH',
        'HMM_BOUNDS', 'SCORE', 'SEQSCORE', 'EVALUE', 'SEQEVALUE',
        'ENV_START', 'ENV_END', 'MODEL_AC', 'NULL', 'FRAGMENTS'
    ],
    "hamap": [
        'ANALYSIS_ID', 'UPI', 'METHOD_AC', 'SUBSTR(RELNO_MAJOR', '1', '4)', 'SUBSTR(RELNO_MAJOR', '6', '7)',
        'SEQ_START', 'SEQ_END', '0', '0', '0',
        'NULL', '0', 'SEQSCORE', '0', '0',
        '0', '0', 'MODEL_AC', 'NULL', 'FRAGMENTS'
    ],
    "mobidblite": [
        'ANALYSIS_ID', 'UPI', 'METHOD_AC', 'RELNO_MAJOR', 'RELNO_MINOR',
        'SEQ_START', 'SEQ_END', '0', '0', '0',
        'NULL', '0', '0', '0', '0',
        '0', '0', 'MODEL_AC', 'SEQ_FEATURE', 'FRAGMENTS'
    ],
    "panther": [
        'ANALYSIS_ID', 'UPI', 'METHOD_AC', 'RELNO_MAJOR', 'RELNO_MINOR',
        'SEQ_START', 'SEQ_END', 'HMM_START', 'HMM_END', 'HMM_LENGTH',
        'HMM_BOUNDS', 'SEQSCORE', 'SEQSCORE', 'SEQEVALUE', 'SEQEVALUE',
        'ENV_START', 'ENV_END', 'MODEL_AC', 'NULL', 'FRAGMENTS'
    ],
    "pfam": [
        'ANALYSIS_ID', 'UPI', 'METHOD_AC', 'RELNO_MAJOR', 'RELNO_MINOR',
        'SEQ_START', 'SEQ_END', 'HMM_START', 'HMM_END', 'HMM_LENGTH',
        'HMM_BOUNDS', 'SCORE', 'SEQSCORE', 'EVALUE', 'SEQEVALUE',
        'ENV_START', 'ENV_END', 'MODEL_AC', 'NULL', 'FRAGMENTS'
    ],
    "phobius": [
        'ANALYSIS_ID', 'UPI', 'METHOD_AC', 'RELNO_MAJOR', 'RELNO_MINOR',
        'SEQ_START', 'SEQ_END', '0', '0', '0',
        'NULL', '0', '0', '0', '0',
        '0', '0', 'MODEL_AC', 'NULL', 'FRAGMENTS'
    ],
    "pirsf": [
        'ANALYSIS_ID', 'UPI', 'METHOD_AC', 'RELNO_MAJOR', 'RELNO_MINOR',
        'SEQ_START', 'SEQ_END', 'HMM_START', 'HMM_END', 'HMM_LENGTH',
        'HMM_BOUNDS', 'SCORE', 'SEQSCORE', 'EVALUE', 'SEQEVALUE',
        'ENV_START', 'ENV_END', 'MODEL_AC', 'NULL', 'FRAGMENTS'
    ],
    "prints": [
        'ANALYSIS_ID', 'UPI', 'METHOD_AC', 'RELNO_MAJOR', 'RELNO_MINOR',
        'SEQ_START', 'SEQ_END', '0', '0', 'MOTIF_NUMBER',
        'NULL', '0', 'SEQSCORE', 'PVALUE', 'SEQEVALUE',
        '0', '0', 'MODEL_AC', 'GRAPHSCAN', 'FRAGMENTS'
    ],
    # "prodom": [
    #     'ANALYSIS_ID', 'UPI', 'METHOD_AC', 'RELNO_MAJOR', 'RELNO_MINOR',
    #     'SEQ_START', 'SEQ_END', '0', '0', '0',
    #     'NULL', 'SEQSCORE', 'SEQSCORE', 'SEQEVALUE', 'SEQEVALUE',
    #     '0', '0', 'MODEL_AC', 'NULL', 'FRAGMENTS'
    # ],
    "prosite_patterns": [
        'ANALYSIS_ID', 'UPI', 'METHOD_AC', 'SUBSTR(RELNO_MAJOR', '1', '4)', 'SUBSTR(RELNO_MAJOR', '6', '7)',
        'SEQ_START', 'SEQ_END', '0', '0', '0',
        'LOCATION_LEVEL', '0', '0', '0', '0',
        '0', '0', 'MODEL_AC', 'ALIGNMENT', 'FRAGMENTS'
    ],
    "prosite_profiles": [
        'ANALYSIS_ID', 'UPI', 'METHOD_AC', 'SUBSTR(RELNO_MAJOR', '1', '4)', 'SUBSTR(RELNO_MAJOR', '6', '7)',
        'SEQ_START', 'SEQ_END', '0', '0', '0',
        'NULL', '0', 'SEQSCORE', '0', '0',
        '0', '0', 'MODEL_AC', 'ALIGNMENT', 'FRAGMENTS'
    ],
    "sfld": [
        'ANALYSIS_ID', 'UPI', 'METHOD_AC', 'RELNO_MAJOR', 'RELNO_MINOR',
        'SEQ_START', 'SEQ_END', 'HMM_START', 'HMM_END', 'HMM_LENGTH',
        'HMM_BOUNDS', 'SCORE', 'SEQSCORE', 'EVALUE', 'SEQEVALUE',
        'ENV_START', 'ENV_END', 'MODEL_AC', 'NULL', 'FRAGMENTS'
    ],
    "signalp_euk": [
        'ANALYSIS_ID', 'UPI', 'METHOD_AC', 'RELNO_MAJOR', 'RELNO_MINOR',
        'SEQ_START', 'SEQ_END', '0', '0', '0',
        'NULL', 'SEQSCORE', 'SEQSCORE', '0', '0',
        '0', '0', 'MODEL_AC', 'NULL', 'FRAGMENTS'
    ],
    "signalp_gram_positive": [
        'ANALYSIS_ID', 'UPI', 'METHOD_AC', 'RELNO_MAJOR', 'RELNO_MINOR',
        'SEQ_START', 'SEQ_END', '0', '0', '0',
        'NULL', 'SEQSCORE', 'SEQSCORE', '0', '0',
        '0', '0', 'MODEL_AC', 'NULL', 'FRAGMENTS'
    ],
    "signalp_gram_negative": [
        'ANALYSIS_ID', 'UPI', 'METHOD_AC', 'RELNO_MAJOR', 'RELNO_MINOR',
        'SEQ_START', 'SEQ_END', '0', '0', '0',
        'NULL', 'SEQSCORE', 'SEQSCORE', '0', '0',
        '0', '0', 'MODEL_AC', 'NULL', 'FRAGMENTS'
    ],
    "smart": [
        'ANALYSIS_ID', 'UPI', 'METHOD_AC', 'RELNO_MAJOR', 'RELNO_MINOR',
        'SEQ_START', 'SEQ_END', 'HMM_START', 'HMM_END', 'HMM_LENGTH',
        'HMM_BOUNDS', 'SCORE', 'SEQSCORE', 'EVALUE', 'SEQEVALUE',
        '0', '0', 'MODEL_AC', 'NULL', 'FRAGMENTS'
    ],
    "superfamily": [
        'ANALYSIS_ID', 'UPI', 'METHOD_AC', 'RELNO_MAJOR', 'RELNO_MINOR',
        'SEQ_START', 'SEQ_END', '0', '0', 'HMM_LENGTH',
        'NULL', '0', '0', 'SEQEVALUE', 'SEQEVALUE',
        '0', '0', 'MODEL_AC', 'NULL', 'FRAGMENTS'
    ],
    "tigrfam": [
        'ANALYSIS_ID', 'UPI', 'METHOD_AC', 'RELNO_MAJOR', 'RELNO_MINOR',
        'SEQ_START', 'SEQ_END', 'HMM_START', 'HMM_END', 'HMM_LENGTH',
        'HMM_BOUNDS', 'SCORE', 'SEQSCORE', 'EVALUE', 'SEQEVALUE',
        'ENV_START', 'ENV_END', 'MODEL_AC', 'NULL', 'FRAGMENTS'
    ],
    "tmhmm": [
        'ANALYSIS_ID', 'UPI', 'METHOD_AC', 'RELNO_MAJOR', 'RELNO_MINOR',
        'SEQ_START', 'SEQ_END', '0', '0', '0',
        'NULL', 'SEQSCORE', 'SEQSCORE', '0', '0',
        '0', '0', 'MODEL_AC', 'NULL', 'FRAGMENTS'
    ],
}

# Partition in SITE
SITE_PARITIONS = {
    "cdd": "CDD",
    "sfld": "SFLD"
}

# Columns to select when inserting site matches in SITE
SITE_SELECT = ['UPI', 'ANALYSIS_ID', 'METHOD_AC', 'LOC_START', 'LOC_END',
               'NUM_SITES', 'RESIDUE', 'RES_START', 'RES_END', 'DESCRIPTION']


@dataclass
class Analysis:
    id: int
    type: str
    name: str
    is_active: bool
    table: str
    persisted: int

    def is_ready(self, cur: Cursor, max_upi: str) -> bool:
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
            dict(analysisid=self.id, persisted=self.persisted, maxupi=max_upi)
        )
        if cur.fetchone()[0]:
            return False

        cur.execute(
            """
            SELECT MAX(JOB_END)
            FROM IPRSCAN.IPM_PERSISTED_JOBS@ISPRO
            WHERE ANALYSIS_ID = :1 AND PERSISTED >= :2
            """, (self.id, self.persisted)
        )
        row = cur.fetchone()
        return row and row[0] >= max_upi


def get_analyses(cur: Cursor, mode: str, **kwargs) -> List[Analysis]:
    ids = kwargs.get("ids", [])
    stable = kwargs.get("stable", True)

    if mode not in ("matches", "sites"):
        raise ValueError("supported modes: matches, sites")

    if ids:
        params = [':' + str(i+1) for i in range(len(ids))]
        sql_filter = f"A.ANALYSIS_ID IN ({','.join(params)})"
        params = tuple(ids)
    elif stable:
        sql_filter = ("A.ANALYSIS_ID IN (SELECT IPRSCAN_SIG_LIB_REL_ID "
                      "FROM INTERPRO.IPRSCAN2DBCODE)")
        params = tuple()
    else:
        sql_filter = "A.ACTIVE = 1"
        params = tuple()

    cur.execute(
        f"""
        SELECT
            A.ANALYSIS_ID,
            A.ANALYSIS_NAME,
            A.ANALYSIS_TYPE,
            A.ACTIVE,
            B.MATCH_TABLE,
            B.SITE_TABLE
          FROM IPM_ANALYSIS@ISPRO A
          INNER JOIN IPM_ANALYSIS_MATCH_TABLE@ISPRO B
            ON A.ANALYSIS_MATCH_TABLE_ID = B.ID
          WHERE {sql_filter}
          ORDER BY A.ANALYSIS_NAME
        """, params
    )

    analyses = []
    for a_id, a_name, a_type, a_active, match_table, site_table in cur:
        """
        CDD/SFLD
            - PERSISTED=1 when matches are ready
            - PERSISTED=2 when matches and sites are ready

        Others:
            - PERSISTED=2 when matches are ready
        """
        if mode == "matches":
            persisted = 1 if a_type in ("cdd", "sfld") else 2
            table = match_table.upper()
        elif site_table:
            persisted = 2
            table = site_table.upper()
        else:
            continue

        analyses.append(Analysis(
            id=a_id,
            type=a_type,
            name=a_name,
            is_active=a_active == 1,
            table=table,
            persisted=persisted
        ))

    return analyses


def import_from_ispro(cur: Cursor, src: str, dst: str):
    oracle.drop_mview(cur, dst)
    oracle.drop_table(cur, dst, purge=True)
    cur.execute(
        f"""
        CREATE TABLE IPRSCAN.{dst} NOLOGGING
        AS
        SELECT *
        FROM IPRSCAN.{src}@ISPRO
        """
    )

    """
    Use remote table name to have fewer characters (no prefix)
    as Oracle < 12.2 do not allow object names longer than 30 characters
    """
    cur.execute(
        f"""
        CREATE INDEX {src}$ID
        ON IPRSCAN.{dst} (ANALYSIS_ID)
        NOLOGGING
        """
    )
    cur.execute(
        f"""
        CREATE INDEX {src}$UPI
        ON IPRSCAN.{dst} (UPI)
        NOLOGGING
        """
    )


def update_analyses(url: str, remote_table: str, partitioned_table: str,
                    analyses: Sequence[Tuple[int, str, Sequence[str]]],
                    force_import: bool = False):
    """
    Update matches for member database analyses.
    :param url: Oracle connection string
    :param remote_table: Match table in ISPRO
    :param partitioned_table: Partitioned table in production database
    :param analyses: Sequence of analyses (analysis ID, partition name
                     in `partitioned_table`, columns to select)
    :param force_import: If True, import data from ISPRO regardless of the UPI
    """
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.execute("SELECT MAX(UPI) FROM UNIPARC.PROTEIN")
    max_upi, = cur.fetchone()

    part2id = {}
    up_to_date = 0
    for analysis_id, partition, columns in analyses:
        cur.execute(
            f"""
            SELECT MAX(ANALYSIS_ID)
            FROM IPRSCAN.{partitioned_table} PARTITION ({partition})
            """
        )
        row = cur.fetchone()
        part2id[partition] = (row[0] if row else None, analysis_id)

        cur.execute(
            f"""
            SELECT MAX(UPI)
            FROM IPRSCAN.{partitioned_table} PARTITION ({partition})
            WHERE ANALYSIS_ID = :1
            """, (analysis_id,)
        )
        row = cur.fetchone()
        upi = row[0] if row else None

        if upi and upi >= max_upi:
            # Data in `partitioned_table` >= UniParc: no need to refresh data
            logger.debug(f"{partition} ({analysis_id}): up-to-date")
            up_to_date += 1
        else:
            logger.debug(f"{partition} ({analysis_id}): outdated")

    if up_to_date == len(analyses):
        cur.close()
        con.close()
        return

    local_table = PREFIX + remote_table
    try:
        cur.execute(f"SELECT MAX(UPI) FROM IPRSCAN.{local_table}")
    except cx_Oracle.DatabaseError as exc:
        error, = exc.args
        if error.code == 942:
            # ORA-00942: Table or view does not exist
            upi_loc = None
        else:
            raise exc
    else:
        row = cur.fetchone()
        upi_loc = row[0] if row else None

    if not upi_loc or upi_loc < max_upi or force_import:
        """
        No matches for the highest UPI: need to import table from ISPRO
        All analyses for this table are imported
            (i.e. previous versions are not ignored)
        """
        logger.debug(f"importing {remote_table}@ISPRO -> {local_table}")
        import_from_ispro(cur, remote_table, local_table)

    # Create temporary table for the partition exchange
    tmp_table = f"IPRSCAN.{remote_table}"
    oracle.drop_table(cur, tmp_table, purge=True)
    cur.execute(
        f"""
        CREATE TABLE {tmp_table} NOLOGGING
        AS
        SELECT *
        FROM IPRSCAN.{partitioned_table}
        WHERE 1 = 0
        """
    )

    for analysis_id, partition, columns in analyses:
        logger.debug(f"{partition} ({analysis_id}): updating")

        # Truncate table (if several analyses for one table, e.g. SignalP)
        oracle.truncate_table(cur, tmp_table, reuse_storage=True)

        # Insert only one analysis ID
        cur.execute(
            f"""
            INSERT /*+ APPEND */ INTO {tmp_table}
            SELECT {', '.join(columns)}
            FROM IPRSCAN.{local_table}
            WHERE ANALYSIS_ID = :1
            """, (analysis_id,)
        )
        con.commit()

        prev_val, new_val = part2id[partition]
        if prev_val is not None and prev_val != new_val:
            """
            Different ANALYSIS_ID (database update):
            1. TRUNCATE the partition, to remove rows with the old ANALYSIS_ID
            2. Modify the partition (remove old value)
            3. Modify the partition (add new value)
            """
            logger.debug(f"{partitioned_table} ({partition}): "
                         f"{prev_val} -> {new_val}")
            cur.execute(
                f"""
                ALTER TABLE IPRSCAN.{partitioned_table} 
                TRUNCATE PARTITION {partition}
                """
            )
            cur.execute(
                f"""
                ALTER TABLE IPRSCAN.{partitioned_table}
                MODIFY PARTITION {partition}
                ADD VALUES ({new_val})
                """
            )
            cur.execute(
                f"""
                ALTER TABLE IPRSCAN.{partitioned_table}
                MODIFY PARTITION {partition}
                DROP VALUES ({prev_val})
                """
            )

        # Exchange partition with temp table
        cur.execute(
            f"""
            ALTER TABLE IPRSCAN.{partitioned_table}
            EXCHANGE PARTITION {partition}
            WITH TABLE {tmp_table}
            """
        )

    # Drop temporary table
    oracle.drop_table(cur, tmp_table, purge=True)

    cur.close()
    con.close()


def import_matches(url: str, **kwargs):
    databases = kwargs.get("databases", [])
    force_import = kwargs.get("force_import", False)
    threads = kwargs.get("threads", 1)

    if databases:  # expects a sequence of Database objects
        databases = {db.analysis_id for db in databases}

    con = cx_Oracle.connect(url)
    cur = con.cursor()

    pending = {}
    for analysis in get_analyses(cur, mode="matches"):
        if databases and analysis.id not in databases:
            continue

        try:
            columns = MATCH_SELECT[analysis.type]
        except KeyError:
            logger.warning(f"ignoring analysis {analysis.name}")
            continue

        try:
            pending[analysis.table].append((analysis, analysis.type, columns))
        except KeyError:
            pending[analysis.table] = [(analysis, analysis.type, columns)]

    cur.execute("SELECT MAX(UPI) FROM UNIPARC.PROTEIN")
    max_upi, = cur.fetchone()
    cur.close()
    con.close()

    if not pending:
        logger.info("No databases to import")
        return

    with ThreadPoolExecutor(max_workers=threads) as executor:
        running = []
        failed = 0

        while True:
            con = cx_Oracle.connect(url)
            cur = con.cursor()

            tmp = []
            for f, table, names in running:
                if not f.done():
                    tmp.append((f, table, names))
                    continue

                try:
                    f.result()
                except Exception as exc:
                    for name in names:
                        logger.info(f"{name:<35} failed ({exc})")
                    failed += 1
                else:
                    for name in names:
                        logger.info(f"{name:<35} done")

            running = tmp

            tmp = {}
            for table, analyses in pending.items():
                ready = []
                for analysis, partition, columns in analyses:
                    if analysis.is_ready(cur, max_upi):
                        ready.append((analysis.id, partition, columns))

                if len(ready) < len(analyses):
                    # Not ready
                    tmp[table] = analyses
                    continue

                names = [e[0].name for e in analyses]
                for name in names:
                    logger.info(f"{name:<35} ready")

                args = (url, table, "MV_IPRSCAN", ready, force_import)
                f = executor.submit(update_analyses, *args)

                running.append((f, table, names))

            pending = tmp

            cur.close()
            con.close()

            if pending or running:
                time.sleep(600)
            else:
                break

    if failed:
        raise RuntimeError(f"{failed} errors")

    con = cx_Oracle.connect(url)
    cur = con.cursor()

    for index in oracle.get_indexes(cur, "IPRSCAN", "MV_IPRSCAN"):
        if index["unusable"]:
            logger.info(f"rebuilding index {index['name']}")
            oracle.catch_temp_error(fn=oracle.rebuild_index,
                                    args=(cur, index["name"]))

    logger.info("gathering statistics")
    oracle.gather_stats(cur, "IPRSCAN", "MV_IPRSCAN")

    cur.close()
    con.close()

    logger.info("complete")


def import_sites(url: str, **kwargs):
    databases = kwargs.get("databases", [])
    force_import = kwargs.get("force_import", False)
    threads = kwargs.get("threads", 1)

    if databases:  # expects a sequence of Database objects
        databases = {db.analysis_id for db in databases}

    con = cx_Oracle.connect(url)
    cur = con.cursor()

    pending = {}
    for analysis in get_analyses(cur, mode="sites"):
        if databases and analysis.id not in databases:
            continue

        try:
            partition = SITE_PARITIONS[analysis.type]
        except KeyError:
            logger.warning(f"ignoring analysis {analysis.name}")
            continue

        try:
            pending[analysis.table].append((analysis, partition, SITE_SELECT))
        except KeyError:
            pending[analysis.table] = [(analysis, partition, SITE_SELECT)]

    cur.execute("SELECT MAX(UPI) FROM UNIPARC.PROTEIN")
    max_upi, = cur.fetchone()
    cur.close()
    con.close()

    if not pending:
        logger.info("No databases to import")
        return

    with ThreadPoolExecutor(max_workers=threads) as executor:
        running = []
        failed = 0

        while True:
            con = cx_Oracle.connect(url)
            cur = con.cursor()

            tmp = []
            for f, table, names in running:
                if not f.done():
                    tmp.append((f, table, names))
                    continue

                try:
                    f.result()
                except Exception as exc:
                    logger.error(f"{names}: failed ({exc})")
                    failed += 1
                else:
                    logger.info(f"{names}: done")

            running = tmp

            tmp = {}
            for table, analyses in pending.items():
                ready = []
                for analysis, partition, columns in analyses:
                    if analysis.is_ready(cur, max_upi):
                        ready.append((analysis.id, partition, columns))

                if len(ready) < len(analyses):
                    # Not ready
                    tmp[table] = analyses
                    continue

                names = ', '.join(e[0].name for e in analyses)
                logger.info(f"{names}: ready")

                args = (url, table, "SITE", ready, force_import)
                f = executor.submit(update_analyses, *args)

                running.append((f, table, names))

            pending = tmp

            cur.close()
            con.close()

            if pending or running:
                time.sleep(600)
            else:
                break

    if failed:
        raise RuntimeError(f"{failed} errors")

    con = cx_Oracle.connect(url)
    cur = con.cursor()

    for index in oracle.get_indexes(cur, "IPRSCAN", "SITE"):
        if index["unusable"]:
            logger.info(f"rebuilding index {index['name']}")
            oracle.catch_temp_error(fn=oracle.rebuild_index,
                                    args=(cur, index["name"]))

    logger.info("gathering statistics")
    oracle.gather_stats(cur, "IPRSCAN", "SITE")

    cur.close()
    con.close()

    logger.info("complete")


# # TODO: check if function can be deleted
# def update_mv_iprscan_mini(url: str, databases: Sequence[Database]):
#     """
#     Import the previous and the new analyses (for one or more member databases)
#     from ISPRO
#
#     :param url: Oracle connection string (IPRSCAN user)
#     :param databases: sequence of member databases
#     """
#     con = cx_Oracle.connect(url)
#     cur = con.cursor()
#
#     cur.execute("SELECT MAX(UPI) FROM UNIPARC.PROTEIN")
#     max_upi, = cur.fetchone()
#
#     tables = {}
#     for db in databases:
#         """
#         Get up to two most recent analyses
#         We assume that all analyses for a given member database use the same
#         table in ISPRO
#         """
#         cur.execute(
#             f"""
#             SELECT ANALYSIS_ID
#             FROM IPRSCAN.IPM_ANALYSIS@ISPRO
#             WHERE ANALYSIS_MATCH_TABLE_ID = (
#               SELECT ANALYSIS_MATCH_TABLE_ID
#               FROM IPRSCAN.IPM_ANALYSIS@ISPRO
#               WHERE ANALYSIS_ID = :1
#             )
#             ORDER BY ANALYSIS_ID
#             """, (db.analysis_id,)
#         )
#         ids = [analysis_id for analysis_id, in cur][-2:]  # up to two
#         analyses = get_analyses(cur, use_matches=True, ids=ids)
#         if len(analyses) == 2:
#             old_analyse, new_analyse = analyses
#
#             # # Update IPRSCAN2DBCODE
#             # cur.execute(
#             #     """
#             #     UPDATE INTERPRO.IPRSCAN2DBCODE
#             #     SET IPRSCAN_SIG_LIB_REL_ID = :1
#             #     WHERE DBCODE = :2
#             #     """, (new_analyse.id, dbcode)
#             # )
#
#             o_ready = old_analyse.is_ready(cur, max_upi)
#             n_ready = new_analyse.is_ready(cur, max_upi)
#             if not (o_ready and n_ready):
#                 cur.close()
#                 con.close()
#                 raise RuntimeError(f"{db.name}: not ready")
#             elif new_analyse.type in MATCH_SELECT:
#                 columns = MATCH_SELECT[new_analyse.type]
#             elif old_analyse.type in MATCH_SELECT:
#                 columns = MATCH_SELECT[old_analyse.type]
#             else:
#                 cur.close()
#                 con.close()
#                 raise RuntimeError(f"{db.name}: not supported")
#
#             tab = new_analyse.table  # old_analyse.table should be the same
#             obj = (old_analyse.id, new_analyse.id, columns)
#             try:
#                 tables[tab]["analyses"].append(obj)
#             except KeyError:
#                 tables[tab] = {
#                     "analyses": [obj],
#                     "name": db.name
#                 }
#         else:
#             cur.close()
#             con.close()
#             raise NotImplementedError("new member database")  # todo
#
#     oracle.drop_table(cur, "IPRSCAN.MV_IPRSCAN_MINI", purge=True)
#     cur.execute(
#         f"""
#         CREATE TABLE IPRSCAN.MV_IPRSCAN_MINI NOLOGGING
#         AS
#         SELECT *
#         FROM IPRSCAN.MV_IPRSCAN
#         WHERE 1 = 0
#         """
#     )
#
#     for remote_table, obj in tables.items():
#         logger.info(f"importing {obj['name']}")
#         local_table = PREFIX + remote_table
#
#         # Import data from ISPRO
#         import_from_ispro(cur, remote_table, local_table)
#
#         for old_id, new_id, columns in obj["analyses"]:
#             cur.execute(
#                 f"""
#                 INSERT /*+ APPEND */ INTO IPRSCAN.MV_IPRSCAN_MINI
#                 SELECT {', '.join(columns)}
#                 FROM IPRSCAN.{local_table}
#                 WHERE ANALYSIS_ID IN (:1, :2)
#                 """, (old_id, new_id)
#             )
#             con.commit()
#
#     logger.info("indexing")
#     cur.execute(
#         f"""
#         CREATE INDEX IPRSCAN.MV_IPRSCAN_MINI$ID
#         ON IPRSCAN.MV_IPRSCAN_MINI (ANALYSIS_ID)
#         NOLOGGING
#         """
#     )
#
#     cur.close()
#     con.close()
#     logger.info("complete")
