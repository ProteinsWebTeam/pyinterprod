import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from typing import List, Sequence, Tuple

import cx_Oracle
from cx_Oracle import Cursor

from pyinterprod import logger
from pyinterprod.utils import oracle

PREFIX = "MV_"

"""
Keys: analysis names in ISPRO
Values: columns to select when inserting matches in MV_IPRSCAN (IPPRO)
        and partitions in MV_IPRSCAN
"""
MATCH_PARTITIONS = {
    "AntiFam": {
        "columns": [
            'ANALYSIS_ID', 'UPI', 'METHOD_AC', 'RELNO_MAJOR', 'RELNO_MINOR',
            'SEQ_START', 'SEQ_END', 'HMM_START', 'HMM_END', 'HMM_LENGTH',
            'HMM_BOUNDS', 'SCORE', 'SEQSCORE', 'EVALUE', 'SEQEVALUE',
            'ENV_START', 'ENV_END', 'MODEL_AC', 'NULL', 'FRAGMENTS'
        ],
        "partition": "antifam"
    },
    "CATH-Gene3D": {
        "columns": [
            'ANALYSIS_ID', 'UPI', 'METHOD_AC', 'RELNO_MAJOR', 'RELNO_MINOR',
            'SEQ_START', 'SEQ_END', 'HMM_START', 'HMM_END', 'HMM_LENGTH',
            'HMM_BOUNDS', 'SCORE', 'SEQSCORE', 'EVALUE', 'SEQEVALUE',
            'ENV_START', 'ENV_END', 'MODEL_AC', 'NULL', 'FRAGMENTS'
        ],
        "partition": "gene3d"
    },
    "CDD": {
        "columns": [
            'ANALYSIS_ID', 'UPI', 'METHOD_AC', 'RELNO_MAJOR', 'RELNO_MINOR',
            'SEQ_START', 'SEQ_END', '0', '0', '0',
            'NULL', 'SEQSCORE',  'SEQSCORE', 'SEQEVALUE', 'SEQEVALUE',
            '0', '0', 'MODEL_AC', 'NULL', 'FRAGMENTS'
        ],
        "partition": "cdd"
    },
    "COILS": {
        "columns": [
            'ANALYSIS_ID', 'UPI', 'METHOD_AC', 'RELNO_MAJOR', 'RELNO_MINOR',
            'SEQ_START', 'SEQ_END', '0', '0', '0',
            'NULL', '0', '0', '0', '0',
            '0', '0', 'MODEL_AC', 'NULL', 'FRAGMENTS'
        ],
        "partition": "coils"
    },
    "FunFam": {
        "columns": [
            'ANALYSIS_ID', 'UPI', 'METHOD_AC', 'RELNO_MAJOR', 'RELNO_MINOR',
            'SEQ_START', 'SEQ_END', 'HMM_START', 'HMM_END', 'HMM_LENGTH',
            'HMM_BOUNDS', 'SCORE', 'SEQSCORE', 'EVALUE', 'SEQEVALUE',
            'ENV_START', 'ENV_END', 'MODEL_AC', 'NULL', 'FRAGMENTS'
        ],
        "partition": "funfam"
    },
    "HAMAP": {
        "columns": [
            'ANALYSIS_ID', 'UPI', 'METHOD_AC', 'SUBSTR(RELNO_MAJOR', '1', '4)',
            'SUBSTR(RELNO_MAJOR', '6', '7)',
            'SEQ_START', 'SEQ_END', '0', '0', '0',
            'NULL', '0', 'SEQSCORE', '0', '0',
            '0', '0', 'MODEL_AC', 'NULL', 'FRAGMENTS'
        ],
        "partition": "hamap"
    },
    "MobiDB Lite": {
        "columns": [
            'ANALYSIS_ID', 'UPI', 'METHOD_AC', 'RELNO_MAJOR', 'RELNO_MINOR',
            'SEQ_START', 'SEQ_END', '0', '0', '0',
            'NULL', '0', '0', '0', '0',
            '0', '0', 'MODEL_AC', 'SEQ_FEATURE', 'FRAGMENTS'
        ],
        "partition": "mobidblite"
    },
    "PANTHER": {
        "columns": [
            'ANALYSIS_ID', 'UPI', 'METHOD_AC', 'RELNO_MAJOR', 'RELNO_MINOR',
            'SEQ_START', 'SEQ_END', 'HMM_START', 'HMM_END', 'HMM_LENGTH',
            'HMM_BOUNDS', 'SEQSCORE', 'SEQSCORE', 'SEQEVALUE', 'SEQEVALUE',
            'ENV_START', 'ENV_END', 'MODEL_AC', 'NULL', 'FRAGMENTS'
        ],
        "partition": "panther"
    },
    "Pfam": {
        "columns": [
            'ANALYSIS_ID', 'UPI', 'METHOD_AC', 'RELNO_MAJOR', 'RELNO_MINOR',
            'SEQ_START', 'SEQ_END', 'HMM_START', 'HMM_END', 'HMM_LENGTH',
            'HMM_BOUNDS', 'SCORE', 'SEQSCORE', 'EVALUE', 'SEQEVALUE',
            'ENV_START', 'ENV_END', 'MODEL_AC', 'NULL', 'FRAGMENTS'
        ],
        "partition": "pfam"
    },
    "Phobius": {
        "columns": [
            'ANALYSIS_ID', 'UPI', 'METHOD_AC', 'RELNO_MAJOR', 'RELNO_MINOR',
            'SEQ_START', 'SEQ_END', '0', '0', '0',
            'NULL', '0', '0', '0', '0',
            '0', '0', 'MODEL_AC', 'NULL', 'FRAGMENTS'
        ],
        "partition": "phobius"
    },
    "PIRSF": {
        "columns": [
            'ANALYSIS_ID', 'UPI', 'METHOD_AC', 'RELNO_MAJOR', 'RELNO_MINOR',
            'SEQ_START', 'SEQ_END', 'HMM_START', 'HMM_END', 'HMM_LENGTH',
            'HMM_BOUNDS', 'SCORE', 'SEQSCORE', 'EVALUE', 'SEQEVALUE',
            'ENV_START', 'ENV_END', 'MODEL_AC', 'NULL', 'FRAGMENTS'
        ],
        "partition": "pirsf"
    },
    "PRINTS": {
        "columns": [
            'ANALYSIS_ID', 'UPI', 'METHOD_AC', 'RELNO_MAJOR', 'RELNO_MINOR',
            'SEQ_START', 'SEQ_END', '0', '0', 'MOTIF_NUMBER',
            'NULL', '0', 'SEQSCORE', 'PVALUE', 'SEQEVALUE',
            '0', '0', 'MODEL_AC', 'GRAPHSCAN', 'FRAGMENTS'
        ],
        "partition": "prints"
    },
    "PROSITE patterns": {
        "columns": [
            'ANALYSIS_ID', 'UPI', 'METHOD_AC', 'SUBSTR(RELNO_MAJOR', '1', '4)',
            'SUBSTR(RELNO_MAJOR', '6', '7)',
            'SEQ_START', 'SEQ_END', '0', '0', '0',
            'LOCATION_LEVEL', '0', '0', '0', '0',
            '0', '0', 'MODEL_AC', 'ALIGNMENT', 'FRAGMENTS'
        ],
        "partition": "prosite_patterns"
    },
    "PROSITE profiles": {
        "columns": [
            'ANALYSIS_ID', 'UPI', 'METHOD_AC', 'SUBSTR(RELNO_MAJOR', '1', '4)',
            'SUBSTR(RELNO_MAJOR', '6', '7)',
            'SEQ_START', 'SEQ_END', '0', '0', '0',
            'NULL', '0', 'SEQSCORE', '0', '0',
            '0', '0', 'MODEL_AC', 'ALIGNMENT', 'FRAGMENTS'
        ],
        "partition": "prosite_profiles"
    },
    "SFLD": {
        "columns": [
            'ANALYSIS_ID', 'UPI', 'METHOD_AC', 'RELNO_MAJOR', 'RELNO_MINOR',
            'SEQ_START', 'SEQ_END', 'HMM_START', 'HMM_END', 'HMM_LENGTH',
            'HMM_BOUNDS', 'SCORE', 'SEQSCORE', 'EVALUE', 'SEQEVALUE',
            'ENV_START', 'ENV_END', 'MODEL_AC', 'NULL', 'FRAGMENTS'
        ],
        "partition": "sfld"
    },
    "SignalP_Euk": {
        "columns": [
            'ANALYSIS_ID', 'UPI', 'METHOD_AC', 'RELNO_MAJOR', 'RELNO_MINOR',
            'SEQ_START', 'SEQ_END', '0', '0', '0',
            'NULL', 'SEQSCORE', 'SEQSCORE', '0', '0',
            '0', '0', 'MODEL_AC', 'NULL', 'FRAGMENTS'
        ],
        "partition": "signalp_euk"
    },
    "SignalP_Gram_positive": {
        "columns": [
            'ANALYSIS_ID', 'UPI', 'METHOD_AC', 'RELNO_MAJOR', 'RELNO_MINOR',
            'SEQ_START', 'SEQ_END', '0', '0', '0',
            'NULL', 'SEQSCORE', 'SEQSCORE', '0', '0',
            '0', '0', 'MODEL_AC', 'NULL', 'FRAGMENTS'
        ],
        "partition": "signalp_gram_positive"
    },
    "SignalP_Gram_negative": {
        "columns": [
            'ANALYSIS_ID', 'UPI', 'METHOD_AC', 'RELNO_MAJOR', 'RELNO_MINOR',
            'SEQ_START', 'SEQ_END', '0', '0', '0',
            'NULL', 'SEQSCORE', 'SEQSCORE', '0', '0',
            '0', '0', 'MODEL_AC', 'NULL', 'FRAGMENTS'
        ],
        "partition": "signalp_gram_negative"
    },
    "SMART": {
        "columns": [
            'ANALYSIS_ID', 'UPI', 'METHOD_AC', 'RELNO_MAJOR', 'RELNO_MINOR',
            'SEQ_START', 'SEQ_END', 'HMM_START', 'HMM_END', 'HMM_LENGTH',
            'HMM_BOUNDS', 'SCORE', 'SEQSCORE', 'EVALUE', 'SEQEVALUE',
            '0', '0', 'MODEL_AC', 'NULL', 'FRAGMENTS'
        ],
        "partition": "smart"
    },
    "SUPERFAMILY": {
        "columns": [
            'ANALYSIS_ID', 'UPI', 'METHOD_AC', 'RELNO_MAJOR', 'RELNO_MINOR',
            'SEQ_START', 'SEQ_END', '0', '0', 'HMM_LENGTH',
            'NULL', '0', '0', 'SEQEVALUE', 'SEQEVALUE',
            '0', '0', 'MODEL_AC', 'NULL', 'FRAGMENTS'
        ],
        "partition": "superfamily"
    },
    "TIGRFAMs": {
        "columns": [
            'ANALYSIS_ID', 'UPI', 'METHOD_AC', 'RELNO_MAJOR', 'RELNO_MINOR',
            'SEQ_START', 'SEQ_END', 'HMM_START', 'HMM_END', 'HMM_LENGTH',
            'HMM_BOUNDS', 'SCORE', 'SEQSCORE', 'EVALUE', 'SEQEVALUE',
            'ENV_START', 'ENV_END', 'MODEL_AC', 'NULL', 'FRAGMENTS'
        ],
        "partition": "tigrfam"
    },
    "TMHMM": {
        "columns": [
            'ANALYSIS_ID', 'UPI', 'METHOD_AC', 'RELNO_MAJOR', 'RELNO_MINOR',
            'SEQ_START', 'SEQ_END', '0', '0', '0',
            'NULL', 'SEQSCORE', 'SEQSCORE', '0', '0',
            '0', '0', 'MODEL_AC', 'NULL', 'FRAGMENTS'
        ],
        "partition": "tmhmm"
    },
}

# Columns to select when inserting site matches in SITE
SITE_SELECT = ["ANALYSIS_ID", "UPI_RANGE", "UPI", "METHOD_AC", "LOC_START",
               "LOC_END", "NUM_SITES", "RESIDUE", "RES_START", "RES_END",
               "DESCRIPTION"]

SITE_PARITIONS = {
    "CDD": {
        "columns": SITE_SELECT,
        "partition": "CDD"
    },
    "PIRSR": {
        "columns": SITE_SELECT,
        "partition": "PIRSR"
    },
    "SFLD": {
        "columns": SITE_SELECT,
        "partition": "SFLD"
    },
}


@dataclass
class Analysis:
    id: int
    name: str
    version: str
    is_active: bool
    table: str

    def is_ready(self, cur: Cursor, max_upi: str) -> bool:
        # Check for incomplete jobs that are required (<= max UPI)
        cur.execute(
            """
            SELECT COUNT(*) AS CNT
            FROM IPRSCAN.ANALYSIS_JOBS@ISPRO
            WHERE ANALYSIS_ID = :1
              AND UPI_FROM <= :2
              AND END_TIME IS NULL
            """,
            [self.id, max_upi]
        )
        if cur.fetchone()[0]:
            return False

        # Check that the highest completed job is >= max UPI
        cur.execute(
            """
            SELECT MAX(UPI_TO)
            FROM IPRSCAN.ANALYSIS_JOBS@ISPRO
            WHERE ANALYSIS_ID = :1
              AND END_TIME IS NOT NULL
            """,
            [self.id]
        )
        row = cur.fetchone()
        return row[0] is not None and row[0] >= max_upi


def get_analyses(cur: Cursor, **kwargs) -> List[Analysis]:
    ids = kwargs.get("ids", [])
    match_type = kwargs.get("type", "matches")
    status = kwargs.get("status", "production")

    if match_type not in ("matches", "sites"):
        raise ValueError("supported values for type: matches, sites")
    elif status not in ("production", "active", "all"):
        raise ValueError("supported values for status: production, active, all")

    if ids:
        params = [':' + str(i+1) for i in range(len(ids))]
        sql_filter = f"WHERE A.ID IN ({','.join(params)})"
        params = tuple(ids)
    elif status == "production":
        sql_filter = ("WHERE  A.ID IN (SELECT IPRSCAN_SIG_LIB_REL_ID "
                      "FROM INTERPRO.IPRSCAN2DBCODE)")
        params = tuple()
    elif status == "active":
        sql_filter = "WHERE A.ACTIVE = 'Y'"
        params = tuple()
    else:
        sql_filter = ""
        params = tuple()

    cur.execute(
        f"""
        SELECT A.ID, A.NAME, A.VERSION, A.ACTIVE, T.MATCH_TABLE, T.SITE_TABLE
        FROM IPRSCAN.ANALYSIS@ISPRO A
        INNER JOIN IPRSCAN.ANALYSIS_TABLES@ISPRO T ON A.NAME = T.NAME
        {sql_filter}
        ORDER BY A.NAME
        """, params
    )

    analyses = []
    for a_id, a_name, a_version, a_active, match_table, site_table in cur:
        if match_type == "matches":
            table = match_table.upper()
        elif site_table:
            table = site_table.upper()
        else:
            continue

        analyses.append(Analysis(
            id=a_id,
            name=a_name,
            version=a_version,
            is_active=a_active == "Y",
            table=table
        ))

    return analyses


def import_tables(uri: str, mode: str = "matches", **kwargs):
    databases = kwargs.get("databases", [])
    force = kwargs.get("force", False)
    threads = kwargs.get("threads", 1)

    if mode not in ("matches", "sites"):
        raise ValueError(f"invalid mode '{mode}'")
    elif databases:  # expects a sequence of Database objects
        databases = {db.analysis_id for db in databases}

    con = cx_Oracle.connect(uri)
    cur = con.cursor()

    pending = {}
    for analysis in get_analyses(cur, type=mode):
        if databases and analysis.id not in databases:
            continue

        try:
            pending[analysis.table].append(analysis)
        except KeyError:
            pending[analysis.table] = [analysis]

    cur.execute("SELECT MAX(UPI) FROM UNIPARC.PROTEIN")
    max_upi, = cur.fetchone()
    cur.close()
    con.close()

    if not pending:
        logger.info("No tables to import")
        return
    elif threads < 1:
        threads = len(pending)

    with ThreadPoolExecutor(max_workers=threads) as executor:
        running = []
        failed = 0

        while True:
            con = cx_Oracle.connect(uri)
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
                        logger.error(f"{name:<38} failed: {exc}")
                    failed += 1
                else:
                    for name in names:
                        logger.info(f"{name:<40} done")

            running = tmp

            tmp = {}
            for table, analyses in pending.items():
                ready = []
                for analysis in analyses:
                    if analysis.is_ready(cur, max_upi):
                        ready.append(analysis.id)

                if len(ready) < len(analyses):
                    # Not ready
                    tmp[table] = analyses
                    continue

                f = executor.submit(_import_table, uri, table, ready, force)
                running.append((f, table, [f"{e.name} {e.version}"
                                           for e in analyses]))

            pending = tmp

            cur.close()
            con.close()

            if pending or running:
                time.sleep(60)
            else:
                break

    if failed:
        raise RuntimeError(f"{failed} errors")

    logger.info("done")


def _import_table(uri: str, remote_table: str, analyses: Sequence[int],
                  force: bool = False):
    con = cx_Oracle.connect(uri)
    cur = con.cursor()

    local_table = PREFIX + remote_table
    oracle.drop_mview(cur, local_table)
    oracle.drop_table(cur, local_table, purge=True)

    cur.execute(
        f"""
        CREATE TABLE IPRSCAN.{local_table} NOLOGGING
        AS
        SELECT *
        FROM IPRSCAN.{remote_table}@ISPRO
        """
    )

    """
    Use remote table name to have fewer characters (no prefix)
    as Oracle < 12.2 do not allow object names longer than 30 characters
    """
    cur.execute(
        f"""
            CREATE INDEX {remote_table}$ID
            ON IPRSCAN.{local_table} (ANALYSIS_ID)
            NOLOGGING
            """
    )

    cur.close()
    con.close()


def update_partitions(uri: str, mode: str = "matches", **kwargs):
    databases = kwargs.get("databases", [])
    force = kwargs.get("force", False)
    threads = kwargs.get("threads", 1)

    if mode == "matches":
        partitions = MATCH_PARTITIONS
        partitioned_table = "MV_IPRSCAN"
    elif mode == "sites":
        partitions = SITE_PARITIONS
        partitioned_table = "SITE"
    else:
        raise ValueError(f"invalid mode '{mode}'")

    if databases:  # expects a sequence of Database objects
        databases = {db.analysis_id for db in databases}

    con = cx_Oracle.connect(uri)
    cur = con.cursor()

    tables = {}
    for analysis in get_analyses(cur, type=mode):
        if databases and analysis.id not in databases:
            continue

        try:
            obj = partitions[analysis.name]
        except KeyError:
            logger.warning(f"ignoring {analysis.name} {analysis.version}")
            continue

        item = (analysis, obj["partition"], obj["columns"])

        try:
            tables[analysis.table].append(item)
        except KeyError:
            tables[analysis.table] = [item]

    cur.close()
    con.close()

    if not tables:
        logger.info("No tables to update")
        return
    elif threads < 1:
        threads = len(tables)

    with ThreadPoolExecutor(max_workers=threads) as executor:
        fs = {}
        for table, analyses in tables.items():
            args = (
                uri,
                PREFIX + table,
                partitioned_table,
                [
                    (analysis.id, partition, columns)
                    for analysis, partition, columns in analyses
                ],
                force
            )

            f = executor.submit(_update_partition, *args)
            fs[f] = [f"{analysis.name} {analysis.version}"
                     for analysis, _, _ in analyses]

        failed = 0
        for f in as_completed(fs):
            names = fs[f]

            try:
                f.result()
            except Exception as exc:
                failed += 1
                for name in names:
                    logger.error(f"{name:<38} failed: {exc}")
            else:
                for name in names:
                    logger.info(f"{name:<40} done")

    if failed:
        raise RuntimeError(f"{failed} errors")

    con = cx_Oracle.connect(uri)
    cur = con.cursor()

    for index in oracle.get_indexes(cur, "IPRSCAN", partitioned_table):
        if index["is_unusable"]:
            logger.info(f"rebuilding index {index['name']}")
            oracle.catch_temp_error(fn=oracle.rebuild_index,
                                    args=(cur, index["name"]))

    logger.info("gathering statistics")
    oracle.gather_stats(cur, "IPRSCAN", partitioned_table)

    cur.close()
    con.close()

    logger.info("done")


def _update_partition(uri: str, table: str, partitioned_table: str,
                      analyses: Sequence[Tuple[int, str, Sequence[str]]],
                      force: bool = False):
    """
    Update partitioned table with matches
    :param uri: Oracle connection string
    :param table: Match table
    :param partitioned_table: Partitioned table
    :param analyses: list of analyses to update (analysis ID, partition name
                     in `partitioned_table`, columns to select from `table`)
    :param force: If True, update partition even if up-to-date
    """
    con = cx_Oracle.connect(uri)
    cur = con.cursor()

    if not force:
        # Check if the data is already up-to-date
        up_to_date = 0

        for analysis_id, partition, columns in analyses:
            cur.execute(
                f"""
                SELECT MAX(UPI)
                FROM IPRSCAN.{table}
                WHERE ANALYSIS_ID = :1
                """,
                [analysis_id]
            )
            max_upi_1, = cur.fetchone()

            cur.execute(
                f"""
                SELECT MAX(UPI)
                FROM IPRSCAN.{partitioned_table} PARTITION ({partition})
                WHERE ANALYSIS_ID = :1
                """,
                [analysis_id]
            )
            max_upi_2, = cur.fetchone()

            if max_upi_1 == max_upi_2:
                # This partition is already up-to-date
                up_to_date += 1

        if up_to_date == len(analyses):
            # All analyses are up-to-date
            cur.close()
            con.close()
            return

    tmp_table = f"IPRSCAN.{table}_TMP"

    for analysis_id, partition, columns in analyses:
        # Create temporary table for the partition exchange
        oracle.drop_table(cur, tmp_table, purge=True)

        sql = f"CREATE TABLE {tmp_table}"
        subparts = oracle.get_subpartitions(cur, schema="IPRSCAN",
                                            table=partitioned_table,
                                            partition=partition)

        if subparts:
            """
            The target table is sub-partitioned: the staging table needs
            to be partitioned
            """
            col = subparts[0]["column"]
            subparts = [
                f"PARTITION {s['name']} VALUES ({s['value']})"
                for s in subparts
            ]

            sql += f" PARTITION BY LIST ({col}) ({', '.join(subparts)})"

        cur.execute(
            f"""{sql}
            NOLOGGING
            AS
            SELECT *
            FROM IPRSCAN.{partitioned_table}
            WHERE 1 = 0
            """
        )

        # Insert only one analysis ID
        cur.execute(
            f"""
            INSERT /*+ APPEND */ INTO {tmp_table}
            SELECT {', '.join(columns)}
            FROM IPRSCAN.{table}
            WHERE ANALYSIS_ID = :1
            """,
            [analysis_id]
        )
        con.commit()

        high_value = None
        for p in oracle.get_partitions(cur, "IPRSCAN", partitioned_table):
            if p["name"].lower() == partition.lower():
                high_value = int(p["value"])
                break

        if high_value is None:
            cur.close()
            con.close()
            raise RuntimeError(f"Partition {partition} not found "
                               f"in {partitioned_table}")

        if high_value != analysis_id:
            """
            Different ANALYSIS_ID -> database update:
            1. TRUNCATE the partition, to remove rows with the old ANALYSIS_ID
            2. Modify the partition (remove old value)
            3. Modify the partition (add new value)
            """
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
                ADD VALUES ({analysis_id})
                """
            )
            cur.execute(
                f"""
                ALTER TABLE IPRSCAN.{partitioned_table}
                MODIFY PARTITION {partition}
                DROP VALUES ({high_value})
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


def check_ispro(url: str, match_type: str = "matches",
                status: str = "production"):
    con = cx_Oracle.connect(url)
    cur = con.cursor()

    analyses = {}
    for analysis in get_analyses(cur, type=match_type, status=status):
        try:
            analyses[analysis.table].append(analysis)
        except KeyError:
            analyses[analysis.table] = [analysis]

    cur.execute("SELECT MAX(UPI) FROM UNIPARC.PROTEIN")
    max_upi, = cur.fetchone()

    for table in sorted(analyses):
        for a in analyses[table]:
            status = "ready" if a.is_ready(cur, max_upi) else "pending"
            print(f"{a.id:<3} {a.name:<40} {a.version:<30} {table:<30} "
                  f"{status}")

    cur.close()
    con.close()
