from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Optional, Union

import cx_Oracle
from mundone.task import Task

from pyinterprod import logger
from pyinterprod.uniprot.uniparc import iter_proteins
from pyinterprod.utils import oracle


def get_analyses(obj: Union[str, cx_Oracle.Cursor]) -> dict:
    if isinstance(obj, str):
        con = cx_Oracle.connect(obj)
        cur = con.cursor()
    else:
        con = None
        cur = obj

    cur.execute(
        """
        SELECT A.ID, A.NAME, A.VERSION, A.MAX_UPI, A.I5_DIR,
               T.MATCH_TABLE, T.SITE_TABLE
        FROM IPRSCAN.ANALYSIS A
        INNER JOIN IPRSCAN.ANALYSIS_TABLES T
            ON LOWER(A.NAME) = LOWER(T.NAME)
        WHERE A.ACTIVE = 'Y'
        """
    )

    analyses = {}
    for row in cur:
        analyses[row[0]] = {
            "name": row[1],
            "version": row[2],
            "max_upi": row[3],
            "i5_dir": row[4],
            "tables": {
                "matches": row[5],
                "sites": row[6],
            }
        }

    if con is not None:
        cur.close()
        con.close()

    return analyses


def clean_tables(uri: str, analysis_ids: Optional[list[int]] = None):
    con = cx_Oracle.connect(uri)
    cur = con.cursor()

    analyses = get_analyses(cur)

    table2analyses = {}
    for analysis_id, analysis in analyses.items():
        table = analysis["tables"]["matches"]
        try:
            table2analyses[table].append(analysis_id)
        except KeyError:
            table2analyses[table] = [analysis_id]

        table = analysis["tables"]["sites"]
        if table:
            try:
                table2analyses[table].add(analysis_id)
            except KeyError:
                table2analyses[table] = {analysis_id}

    actions = []
    for table in sorted(table2analyses):
        table = table.upper()
        for p in oracle.get_partitions(cur, "IPRSCAN", table):
            if p["value"] == "DEFAULT":
                continue

            analysis_id = int(p["value"])

            if analysis_ids and analysis_id not in analysis_ids:
                continue

            if analysis_id not in analyses:
                # Obsolete analysis: remove data
                actions.append((
                    f"  - {p['name']:<30}: drop partition",
                    [(
                        f"ALTER TABLE {table} DROP PARTITION {p['name']}", []
                    )]
                ))
                continue

            analysis = analyses[analysis_id]
            max_upi = analysis["max_upi"]

            if analysis_id not in table2analyses[table]:
                # Obsolete analysis: remove data
                actions.append((
                    f"  - {p['name']:<30}: drop partition",
                    [(
                        f"ALTER TABLE {table} DROP PARTITION {p['name']}", []
                    )]
                ))
            elif max_upi:
                cur.execute(
                    """
                    SELECT COUNT(*) 
                    FROM IPRSCAN.ANALYSIS_JOBS 
                    WHERE ANALYSIS_ID = :1
                      AND UPI_FROM > :2
                    """,
                    [analysis_id, max_upi]
                )
                cnt, = cur.fetchone()

                if cnt > 0:
                    # Delete jobs after the max UPI
                    actions.append((
                        f"  - {p['name']:<30}: delete jobs > {max_upi}",
                        [(
                            """
                            DELETE FROM IPRSCAN.ANALYSIS_JOBS
                            WHERE ANALYSIS_ID = :1
                            AND UPI_FROM > :2
                            """,
                            [analysis_id, max_upi]
                        ), (
                            f"""
                            DELETE FROM {table} PARTITION ({p['name']})
                            WHERE UPI_FROM > :1
                            """,
                            [max_upi]
                        )]
                    ))
            else:
                # No max UPI: remove data
                actions.append((
                    f"  - {p['name']:<30}: delete jobs",
                    [(
                        f"DELETE FROM IPRSCAN.ANALYSIS_JOBS "
                        f"WHERE ANALYSIS_ID = :1",
                        [analysis_id]
                    ), (
                        f"ALTER TABLE {table} TRUNCATE PARTITION {p['name']}",
                        []
                    )]
                ))

    if actions:
        print("The following actions will be performed:")
        for descr, queries in actions:
            print(descr)

        if input("Proceed? [y/N] ").lower().strip() == "y":
            for descr, queries in actions:
                for sql, params in queries:
                    cur.execute(sql, params)

            con.commit()
        else:
            print("Canceled")

    cur.close()
    con.close()


def import_uniparc(ispro_uri: str, uniparc_uri: str, top_up: bool = False):
    logger.info("importing sequences from UniParc")
    con = cx_Oracle.connect(ispro_uri)
    cur = con.cursor()

    if top_up:
        cur.execute("SELECT MAX(UPI) FROM UNIPARC.PROTEIN")
        max_upi, = cur.fetchone()
    else:
        max_upi = None
        oracle.drop_mview(cur, "UNIPARC.PROTEIN")
        oracle.drop_table(cur, "UNIPARC.PROTEIN", purge=True)
        cur.execute(
            """
                CREATE TABLE UNIPARC.PROTEIN
                (
                    ID NUMBER(15) NOT NULL,
                    UPI CHAR(13) NOT NULL,
                    TIMESTAMP DATE NOT NULL,
                    USERSTAMP VARCHAR2(30) NOT NULL,
                    CRC64 CHAR(16) NOT NULL,
                    LEN NUMBER(6) NOT NULL,
                    SEQ_SHORT VARCHAR2(4000),
                    SEQ_LONG CLOB,
                    MD5 VARCHAR2(32) NOT NULL
                ) NOLOGGING
            """
        )

    cnt = 0
    records = []
    req = """
        INSERT /*+ APPEND */ 
        INTO UNIPARC.PROTEIN
        VALUES (:1, :2, :3, :4, :5, :6, :7, :8, :9)
    """

    for rec in iter_proteins(uniparc_uri, greather_than=max_upi):
        records.append(rec)
        cnt += 1

        if len(records) == 1000:
            cur.executemany(req, records)
            con.commit()
            records.clear()

    if records:
        cur.executemany(req, records)
        con.commit()
        records.clear()

    if not top_up:
        cur.execute("GRANT SELECT ON UNIPARC.PROTEIN TO PUBLIC")
        cur.execute("CREATE UNIQUE INDEX PK_PROTEIN ON UNIPARC.PROTEIN (UPI)")

    cur.close()
    con.close()

    logger.info(f"\t{cnt:,} sequences imported")


def get_incomplete_jobs(cur: cx_Oracle.Cursor) -> dict:
    cur.execute(
        """
        SELECT ANALYSIS_ID, UPI_FROM, UPI_TO
        FROM IPRSCAN.ANALYSIS_JOBS
        WHERE END_TIME IS NULL
        """
    )

    incomplete_jobs = {}
    for analysis_id, upi_from, upi_to in cur:
        try:
            incomplete_jobs[analysis_id].append((upi_from, upi_to))
        except KeyError:
            incomplete_jobs[analysis_id] = [(upi_from, upi_to)]

    return incomplete_jobs


def add_job(cur: Optional[cx_Oracle.Cursor], analysis_id: int, upi_from: str,
            upi_to: str, uri: Optional[str] = None):
    if uri:
        con = cx_Oracle.connect(uri)
        cur = con.cursor()
    else:
        con = None

    cur.execute(
        """
        UPDATE IPRSCAN.ANALYSIS
        SET MAX_UPI = :upi
        WHERE ID = :analysis_id AND (MAX_UPI IS NULL OR MAX_UPI < :upi)
        """,
        analysis_id=analysis_id, upi=upi_to
    )
    cur.execute(
        """
        INSERT INTO IPRSCAN.ANALYSIS_JOBS 
            (ANALYSIS_ID, UPI_FROM, UPI_TO, SUBMIT_TIME)
        VALUES (:1, :2, :3, SYSDATE)
        """,
        (analysis_id, upi_from, upi_to)
    )
    cur.connection.commit()

    if con is not None:
        cur.close()
        con.close()


def update_job(uri: str, analysis_id: int, upi_from: str, upi_to: str,
               task: Task, max_mem: int, cpu_time: int):
    con = cx_Oracle.connect(uri)
    cur = con.cursor()
    cur.execute(
        """
        UPDATE IPRSCAN.ANALYSIS_JOBS
        SET START_TIME = :1,
            END_TIME = :2,
            MAX_MEMORY = :3,
            LIM_MEMORY = :4,
            CPU_TIME = :5
        WHERE ANALYSIS_ID = :6
            AND UPI_FROM = :7
            AND UPI_TO = :8
            AND END_TIME IS NULL
        """,
        [task.start_time, task.end_time, max_mem, int(task.scheduler["mem"]),
         cpu_time, analysis_id, upi_from, upi_to]
    )
    con.commit()
    cur.close()
    con.close()


def rebuild_indexes(uri: str, analysis_ids: Optional[list[int]] = None):
    con = cx_Oracle.connect(uri)
    cur = con.cursor()

    analyses = get_analyses(cur)

    tables = set()
    for analysis_id, analysis in analyses.items():
        if analysis_ids and analysis_id not in analysis_ids:
            continue

        tables.add(analysis["tables"]["matches"])

        if analysis["tables"]["sites"]:
            tables.add(analysis["tables"]["sites"])

    to_rebuild = set()
    for table in tables:
        for index in oracle.get_indexes(cur, "IPRSCAN", table):
            if index["is_unusable"]:
                to_rebuild.add(index["name"])

    cur.close()
    con.close()

    errors = 0
    with ThreadPoolExecutor(max_workers=8) as executor:
        fs = {}

        for index in to_rebuild:
            f = executor.submit(_rebuild_index, uri, index)
            fs[f] = index

        for f in as_completed(fs):
            index = fs[f]

            try:
                f.result()
            except Exception as exc:
                logger.error(f"failed to rebuild {index}: {exc}")
                errors += 1
            else:
                logger.info(f"{index} rebuilt")

    if errors > 0:
        raise RuntimeError(f"{errors} errors occurred")


def _rebuild_index(uri: str, name: str):
    con = cx_Oracle.connect(uri)
    cur = con.cursor()
    oracle.rebuild_index(cur, name)
    cur.close()
    con.close()
