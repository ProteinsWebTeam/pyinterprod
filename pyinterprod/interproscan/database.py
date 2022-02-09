from typing import Sequence

import cx_Oracle

from pyinterprod.utils import oracle


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

    for table in ("IPM2_JOBS", "IPM2_ANALYSIS", "IPM2_TABLES"):
        oracle.drop_table(cur, f"IPRSCAN.{table}", purge=True)

    cur.execute(
        """
        CREATE TABLE IPRSCAN.IPM2_ANALYSIS
        (
            ID NUMBER(4) NOT NULL
                CONSTRAINT IPM2_ANALYSIS_PK PRIMARY KEY,
            NAME VARCHAR2(30) NOT NULL,
            VERSION VARCHAR2(20) NOT NULL,
            ACTIVE NUMBER(1)
                CONSTRAINT IPM2_ANALYSIS_CHK
                CHECK (ACTIVE in (0, 1)),
            MAX_UPI VARCHAR2(13) DEFAULT NULL,
            I5_DIR VARCHAR2(1000) NOT NULL,
            TIMESTAMP DATE DEFAULT SYSDATE
        )
        """
    )

    cur.execute(
        """
        CREATE TABLE IPRSCAN.IPM2_TABLES
        (
            NAME VARCHAR2(30) NOT NULL PRIMARY KEY,
            MATCH_TABLE VARCHAR2(30) NOT NULL,
            SITE_TABLE VARCHAR2(30)
        )
        """
    )

    cur.execute(
        """
        CREATE TABLE IPRSCAN.IPM2_JOBS
        (
            ANALYSIS_ID NUMBER(4) NOT NULL
                CONSTRAINT IPM2_JOBS_FK
                REFERENCES IPM2_ANALYSIS (ID),
            UPI_FROM VARCHAR2(13) NOT NULL,
            UPI_TO VARCHAR2(13) NOT NULL,
            START_TIME DATE DEFAULT NULL,
            END_TIME DATE DEFAULT NULL,
            MAX_MEMORY NUMBER(6) DEFAULT NULL,
            CONSTRAINT IPM2_JOBS_PK
                PRIMARY KEY (ANALYSIS_ID, UPI_FROM, UPI_TO)
        )
        """
    )

    for analysis_id, (name, version, i5_dir) in active_analyses.items():
        cur.execute(
            """
            INSERT INTO IPRSCAN.IPM2_ANALYSIS
                (ID, NAME, VERSION, ACTIVE, MAX_UPI, I5_DIR)
            VALUES
                (:1, :2, :3, :4, :5, :6)
            """,
            (analysis_id, name, version, 1, hwm.get(analysis_id), i5_dir)
        )

    cur.executemany(
        """
        INSERT INTO IPRSCAN.IPM2_TABLES
        VALUES (:1, :2, :3)
        """,
        tables
    )

    con.commit()
    cur.close()
    con.close()


def get_analyses(cur: cx_Oracle.Cursor) -> dict:
    cur.execute(
        """
        SELECT A.ID, A.NAME, A.VERSION, A.MAX_UPI, A.I5_DIR,
               T.MATCH_TABLE, T.SITE_TABLE
        FROM IPRSCAN.IPM2_ANALYSIS A
        INNER JOIN IPRSCAN.IPM2_TABLES T
            ON A.NAME = T.NAME
        WHERE A.ACTIVE = 1
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

    return analyses


def clean_tables(uri: str):
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

    for table in table2analyses:
        for p in oracle.get_partitions(cur, "IPRSCAN", table.upper()):
            if p["value"] == "DEFAULT":
                continue

            analysis_id = int(p["value"])
            if analysis_id not in table2analyses[table]:
                # Obsolete analysis: remove data
                cur.execute(f"ALTER TABLE {table} "
                            f"DROP PARTITION {p['name']}")
            elif not analyses[analysis_id]["max_upi"]:
                # No max UPI: remove data
                cur.execute(f"ALTER TABLE {table} "
                            f"TRUNCATE PARTITION {p['name']}")

    cur.close()
    con.close()


def reset_max_upi(uri: str, analyses: Sequence[int]):
    con = cx_Oracle.connect(uri)
    cur = con.cursor()

    for analysis_id in analyses:
        cur.execute(
            """
            UPDATE IPRSCAN.IPM2_ANALYSIS
            SET MAX_UPI = NULL
            WHERE ID = :1
            """,
            (analysis_id,)
        )

    con.commit()
    cur.close()
    con.close()
