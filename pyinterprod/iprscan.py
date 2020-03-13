# -*- coding: utf-8 -*-

import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from typing import List, Optional, Sequence, Tuple

import cx_Oracle
from cx_Oracle import Cursor

from pyinterprod import logger
from pyinterprod.utils import oracle


# Columns to select when inserting matches in MV_IPRSCAN
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

# Partitions in MV_IPRSCAN
MATCH_PARTITIONS = {
    "cdd": "CDD",
    "coils": "COILS",
    "gene3d": "GENE3D",
    "hamap": "HAMAP",
    "mobidblite": "MOBIDBLITE",
    "panther": "PANTHER",
    "pfam": "PFAM",
    "phobius": "PHOBIUS",
    "pirsf": "PIRSF",
    "prints": "PRINTS",
    # "prodom": "PRODOM",
    "prosite_patterns": "PROSITE_PATTERNS",
    "prosite_profiles": "PROSITE_PROFILES",
    "sfld": "SFLD",
    "signalp_euk": "SIGNALP_EUK",
    "signalp_gram_positive": "SIGNALP_GRAM_POSITIVE",
    "signalp_gram_negative": "SIGNALP_GRAM_NEGATIVE",
    "smart": "SMART",
    "superfamily": "SUPERFAMILY",
    "tigrfam": "TIGRFAM",
    "tmhmm": "TMHMM"
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
    name: str
    full_name: str
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
            WHERE ANALYSIS_ID = :1 AND PERSISTED >- :2
            """, (self.id, self.persisted)
        )
        row = cur.fetchone()
        return row and row[0] >= max_upi


def get_analyses(url: str, use_matches: bool=True) -> List[Analysis]:
    con = cx_Oracle.connect(url)
    cur = con.cursor()

    cur.execute(
        """
        SELECT
            A.ANALYSIS_ID,
            A.ANALYSIS_NAME,
            A.ANALYSIS_TYPE,
            B.MATCH_TABLE,
            B.SITE_TABLE
          FROM IPM_ANALYSIS@ISPRO A
          INNER JOIN IPM_ANALYSIS_MATCH_TABLE@ISPRO B
            ON A.ANALYSIS_MATCH_TABLE_ID = B.ID
          WHERE A.ACTIVE = 1
        """
    )

    analyses = []
    for analysis_id, name, analysis_type, match_table, site_table in cur:
        """
        CDD/SFLD
            - PERSISTED=1 when matches are ready
            - PERSISTED=2 when matches and sites are ready

        Others:
            - PERSISTED=2 when matches are ready
        """
        if use_matches:
            persisted = 1 if analysis_type in ("cdd", "sfld") else 2
            table = match_table.upper()
        elif site_table:
            persisted = 2
            table = site_table.upper()
        else:
            continue

        analysis = Analysis(analysis_id, analysis_type, name, table, persisted)
        analyses.append(analysis)

    cur.close()
    con.close()
    return analyses


def get_max_upi(cur: Cursor, sql: str) -> Optional[str]:
    try:
        cur.execute(sql)
    except cx_Oracle.DatabaseError as exc:
        error, = exc.args
        if error.code in (942, 2149):
            # ORA-00942: Table or view does not exist
            # ORA-02149: Specified partition does not exist
            return None
        else:
            raise exc
    else:
        row = cur.fetchone()
        return row[0] if row else None


def update_analyses(url: str, table: str, partitioned_table: str,
                    analyses: Sequence[Tuple[int, str, Sequence[str]]]):
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.execute("SELECT MAX(UPI) FROM UNIPARC.PROTEIN")
    max_upi, = cur.fetchone()

    up_to_date = 0
    for analysis_id, partition, columns in analyses:
        sql = f"""
            SELECT MAX(UPI)
            FROM IPRSCAN.{partitioned_table} PARTITION ({partition})
        """
        upi = get_max_upi(cur, sql)

        if upi and upi >= max_upi:
            # Data in `partitioned_table` >= UniParc: no need to refresh data
            up_to_date += 1

    if up_to_date == len(analyses):
        cur.close()
        con.close()
        return

    upi_loc = get_max_upi(cur, f"SELECT MAX(UPI) FROM IPRSCAN.{table}")
    if not upi_loc or upi_loc < max_upi:
        # No matches for the highest UPI: need to import table from ISPRO
        oracle.drop_mview(cur, table)
        oracle.drop_table(cur, table, purge=True)
        cur.execute(
            f"""
            CREATE TABLE IPRSCAN.{table} NOLOGGING
            AS
            SELECT *
            FROM IPRSCAN.{table}@ISPRO
            """
        )
        cur.execute(
            f"""
            CREATE INDEX {table}$ID 
            ON IPRSCAN.{table} (ANALYSIS_ID)
            NOLOGGING
            """
        )
        cur.execute(
            f"""
            CREATE INDEX {table}$UPI 
            ON IPRSCAN.{table} (UPI)
            NOLOGGING
            """
        )

    # Create temporary table for the partition exchange
    tmp_table = f"IPRSCAN.{table}_TMP"
    oracle.drop_table(cur, tmp_table, purge=True)
    cur.execute(
        f"""
        CREATE TABLE {tmp_table}
        AS
        SELECT *
        FROM IPRSCAN.{partitioned_table}
        WHERE 1=0
        """
    )

    for analysis_id, partition, columns in analyses:
        # Truncate table (if several analyses for one table, e.g. SignalP)
        oracle.truncate_table(cur, tmp_table, reuse_storage=True)

        # Insert only the active analysis ID
        # (i.e. data from previous versions is ignored)
        cur.execute(
            f"""
            INSERT /*+ APPEND */ INTO {tmp_table}
            SELECT {', '.join(columns)}
            FROM IPRSCAN.{table}
            WHERE ANALYSIS_ID = :1
            """, (analysis_id,)
        )
        con.commit()

        # Brings new data (in the tmp table) online
        cur.execute(
            f"""
            ALTER TABLE IPRSCAN.{partitioned_table}
            EXCHANGE PARTITION {partition}
            WITH TABLE {tmp_table}
            INCLUDING INDEXES
            WITHOUT VALIDATION
            """
        )

    # Drop temporary table
    oracle.drop_table(cur, tmp_table, purge=True)

    cur.close()
    con.close()


def import_matches(url: str, threads: int=1):
    pending = {}
    for analysis in get_analyses(url, use_matches=True):
        try:
            partition = MATCH_PARTITIONS[analysis.name]
            columns = MATCH_SELECT[analysis.name]
        except KeyError:
            logger.warning(f"ignoring analysis {analysis.full_name}")
            continue

        try:
            pending[analysis.table].append((analysis, partition, columns))
        except KeyError:
            pending[analysis.table] = [(analysis, partition, columns)]

    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.execute("SELECT MAX(UPI) FROM UNIPARC.PROTEIN")
    max_upi, = cur.fetchone()
    cur.close()
    con.close()

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

                names = ', '.join(e[0].full_name for e in analyses)
                logger.info(f"{names}: ready")

                args = (url, table, "MV_IPRSCAN", ready)
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

    logger.info("rebuilding indexes")
    for idx in oracle.get_indexes(cur, "IPRSCAN", "MV_IPRSCAN"):
        oracle.catch_temp_error(fn=oracle.rebuild_index,
                                args=(cur, idx["name"]))

    logger.info("gathering statistics")
    oracle.gather_stats(cur, "IPRSCAN", "MV_IPRSCAN")

    cur.close()
    con.close()

    logger.info("complete")


def import_sites(url: str, threads: int=1):
    pending = {}
    for analysis in get_analyses(url, use_matches=False):
        try:
            partition = SITE_PARITIONS[analysis.name]
        except KeyError:
            logger.warning(f"ignoring analysis {analysis.full_name}")
            continue

        try:
            pending[analysis.table].append((analysis, partition, SITE_SELECT))
        except KeyError:
            pending[analysis.table] = [(analysis, partition, SITE_SELECT)]

    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.execute("SELECT MAX(UPI) FROM UNIPARC.PROTEIN")
    max_upi, = cur.fetchone()
    cur.close()
    con.close()

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

                names = ', '.join(e[0].full_name for e in analyses)
                logger.info(f"{names}: ready")

                args = (url, table, "SITE", ready)
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
