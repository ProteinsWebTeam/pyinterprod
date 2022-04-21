from typing import Optional, Union

import cx_Oracle
from mundone.task import Task

from pyinterprod import logger
from pyinterprod.utils import oracle
from .utils import int_to_upi, upi_to_int


def init_tables(ippro_uri: str, ispro_uri: str, i5_dir: str, others=None):
    if others is None:
        others = []

    con = cx_Oracle.connect(ippro_uri)
    cur = con.cursor()

    active_analyses = {}
    cur.execute(
        """
        SELECT I2D.IPRSCAN_SIG_LIB_REL_ID, D.DBNAME, V.VERSION
        FROM INTERPRO.CV_DATABASE D
        INNER JOIN INTERPRO.IPRSCAN2DBCODE I2D
            ON I2D.DBCODE = D.DBCODE
        INNER JOIN INTERPRO.DB_VERSION V
            ON D.DBCODE = V.DBCODE
        """
    )
    for analysis_id, name, version in cur:
        active_analyses[analysis_id] = (name, version, i5_dir)

    cur.close()
    con.close()

    con = cx_Oracle.connect(ispro_uri)
    cur = con.cursor()
    cur.execute(
        """
        SELECT A.ANALYSIS_ID, T.MATCH_TABLE, T.SITE_TABLE
        FROM IPRSCAN.IPM_ANALYSIS A
        INNER JOIN IPRSCAN.IPM_ANALYSIS_MATCH_TABLE T
        ON A.ANALYSIS_MATCH_TABLE_ID = T.ID
        """
    )
    tables = {}
    for analysis_id, match_table, site_table in cur:
        try:
            name, version, i5_dir = active_analyses[analysis_id]
        except KeyError:
            continue
        else:
            tables[name] = (
                match_table.upper(),
                site_table.upper() if site_table else None
            )

    for analysis_id, name, version, match_table, site_table, i5_dir in others:
        active_analyses[analysis_id] = (name, version, i5_dir)
        tables[name] = (
            match_table.upper(),
            site_table.upper() if site_table else None
        )

    tables = [(k, *v) for k, v in tables.items()]

    cur.execute(
        """
        SELECT ANALYSIS_ID, HWM_SUBMITTED
        FROM IPRSCAN.IPM_HWM
        """
    )
    hwm = dict(cur.fetchall())

    oracle.drop_table(cur, "IPRSCAN.ANALYSIS_JOBS", purge=True)
    oracle.drop_table(cur, "IPRSCAN.ANALYSIS", purge=True)
    oracle.drop_table(cur, "IPRSCAN.ANALYSIS_TABLES", purge=True)

    cur.execute(
        """
        CREATE TABLE IPRSCAN.ANALYSIS
        (
            ID NUMBER(4) NOT NULL
                CONSTRAINT ANALYSIS_PK PRIMARY KEY,
            NAME VARCHAR2(30) NOT NULL,
            VERSION VARCHAR2(20) NOT NULL,
            ACTIVE CHAR(1) DEFAULT 'N' NOT NULL
                CONSTRAINT ANALYSIS_CHK
                CHECK (ACTIVE in ('Y', 'N')),
            MAX_UPI VARCHAR2(13) DEFAULT NULL,
            I5_DIR VARCHAR2(1000) NOT NULL,
            TIMESTAMP DATE DEFAULT SYSDATE
        )
        """
    )

    cur.execute(
        """
        CREATE TABLE IPRSCAN.ANALYSIS_TABLES
        (
            NAME VARCHAR2(30) NOT NULL PRIMARY KEY,
            MATCH_TABLE VARCHAR2(30) NOT NULL,
            SITE_TABLE VARCHAR2(30)
        )
        """
    )

    cur.execute(
        """
        CREATE TABLE IPRSCAN.ANALYSIS_JOBS
        (
            ANALYSIS_ID NUMBER(4) NOT NULL
                CONSTRAINT ANALYSIS_JOBS_FK
                REFERENCES ANALYSIS (ID),
            UPI_FROM VARCHAR2(13) NOT NULL,
            UPI_TO VARCHAR2(13) NOT NULL,
            SUBMIT_TIME DATE DEFAULT SYSDATE NOT NULL,
            START_TIME DATE DEFAULT NULL,
            END_TIME DATE DEFAULT NULL,
            MAX_MEMORY NUMBER(6) DEFAULT NULL,
            LIM_MEMORY NUMBER(6) DEFAULT NULL,
            CPU_TIME NUMBER(9) DEFAULT NULL,
            CONSTRAINT ANALYSIS_JOBS_PK
                PRIMARY KEY (ANALYSIS_ID, UPI_FROM, UPI_TO)
        )
        """
    )

    for analysis_id, (name, version, i5_dir) in active_analyses.items():
        cur.execute(
            """
            INSERT INTO IPRSCAN.ANALYSIS
                (ID, NAME, VERSION, ACTIVE, MAX_UPI, I5_DIR)
            VALUES
                (:1, :2, :3, :4, :5, :6)
            """,
            (analysis_id, name, version, 'Y', hwm.get(analysis_id), i5_dir)
        )

    cur.executemany(
        """
        INSERT INTO IPRSCAN.ANALYSIS_TABLES
        VALUES (:1, :2, :3)
        """,
        tables
    )

    con.commit()
    cur.close()
    con.close()


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
            ON A.NAME = T.NAME
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


def iter_proteins(uri: str, greather_than: Optional[str] = None):
    con = cx_Oracle.connect(uri)
    cur = con.cursor()
    cur.outputtypehandler = oracle.clob_as_str

    sql = """
        SELECT ID, UPI, TIMESTAMP, USERSTAMP, CRC64, LEN, SEQ_SHORT, 
               SEQ_LONG, MD5
        FROM UNIPARC.PROTEIN      
    """

    if greather_than is not None:
        sql += "WHERE UPI > :1"
        params = [greather_than]
    else:
        params = []

    cur.execute(sql, params)
    yield from cur
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


def split_incomplete_jobs(uri: str, new_size: int):
    con = cx_Oracle.connect(uri)
    cur = con.cursor()
    jobs = get_incomplete_jobs(cur)

    params = []
    for analysis_id, analysis_jobs in jobs.items():
        for upi_from, upi_to in analysis_jobs:
            while upi_from <= upi_to:
                start = upi_to_int(upi_from)
                new_start = start + (new_size - 1)
                new_upi_to = int_to_upi(new_start)

                params.append((analysis_id, upi_from, new_upi_to))

                upi_from = int_to_upi(new_start + 1)

    cur.execute(
        """
        DELETE FROM IPRSCAN.ANALYSIS_JOBS
        WHERE END_TIME IS NULL
        """
    )

    cur.executemany(
        """
        INSERT INTO IPRSCAN.ANALYSIS_JOBS 
            (ANALYSIS_ID, UPI_FROM, UPI_TO, SUBMIT_TIME)
        VALUES (:1, :2, :3, SYSDATE)
        """,
        params
    )

    con.commit()
    cur.close()
    con.close()


def add_job(cur: cx_Oracle, analysis_id: int, upi_from: str, upi_to: str):
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


def reset_job(uri: str, analysis_id: int, upi_from: str, upi_to: str):
    con = cx_Oracle.connect(uri)
    cur = con.cursor()
    cur.execute(
        """
        UPDATE IPRSCAN.ANALYSIS_JOBS
        SET SUBMIT_TIME = SYSDATE,
            START_TIME = NULL,
            END_TIME = NULL,
            MAX_MEMORY = NULL,
            LIM_MEMORY = NULL,
            CPU_TIME = NULL
        WHERE ANALYSIS_ID = :1 
          AND UPI_FROM = :2 
          AND UPI_TO = :3
        """,
        [analysis_id, upi_from, upi_to]
    )
    con.commit()
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
        """,
        [task.start_time, task.end_time, max_mem, int(task.scheduler["mem"]),
         cpu_time, analysis_id, upi_from, upi_to]
    )
    con.commit()
    cur.close()
    con.close()
