from typing import Optional, Sequence

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

    for table in ("ANALYSIS_JOBS", "ANALYSIS", "ANALYSIS_TABLES"):
        oracle.drop_table(cur, f"IPRSCAN.{table}", purge=True)

    cur.execute(
        """
        CREATE TABLE IPRSCAN.ANALYSIS
        (
            ID NUMBER(4) NOT NULL
                CONSTRAINT ANALYSIS_PK PRIMARY KEY,
            NAME VARCHAR2(30) NOT NULL,
            VERSION VARCHAR2(20) NOT NULL,
            ACTIVE NUMBER(1)
                CONSTRAINT ANALYSIS_CHK
                CHECK (ACTIVE in (0, 1)),
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
            START_TIME DATE DEFAULT NULL,
            END_TIME DATE DEFAULT NULL,
            MAX_MEMORY NUMBER(6) DEFAULT NULL,
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
            (analysis_id, name, version, 1, hwm.get(analysis_id), i5_dir)
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


def get_analyses(cur: cx_Oracle.Cursor) -> dict:
    cur.execute(
        """
        SELECT A.ID, A.NAME, A.VERSION, A.MAX_UPI, A.I5_DIR,
               T.MATCH_TABLE, T.SITE_TABLE
        FROM IPRSCAN.ANALYSIS A
        INNER JOIN IPRSCAN.ANALYSIS_TABLES T
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
            UPDATE IPRSCAN.ANALYSIS
            SET MAX_UPI = NULL
            WHERE ID = :1
            """,
            (analysis_id,)
        )

    con.commit()
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
        sql += "WHERE UPI > :upi"
        params = {"upi": greather_than}
    else:
        params = {}

    cur.execute(sql, params)
    yield from cur
    cur.close()
    con.close()


def import_uniparc(ispro_uri: str, uniparc_uri: str, top_up: bool = False):
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

    records = []
    req = """
        INSERT /*+ APPEND */ 
        INTO UNIPARC.PROTEIN
        VALUES (:1, :2, :3, :4, :5, :6, :7, :8, :9)
    """

    for rec in iter_proteins(uniparc_uri, greather_than=max_upi):
        records.append(rec)

        if len(records) == 1000:
            cur.executemany(req, records)
            con.commit()
            records.clear()

    if records:
        cur.executemany(req, records)
        con.commit()
        records.clear()

    if not top_up:
        cur.execute(f"GRANT SELECT ON UNIPARC.PROTEIN TO PUBLIC")
        cur.execute("CREATE UNIQUE INDEX PK_PROTEIN ON UNIPARC.PROTEIN (UPI)")

    cur.close()
    con.close()
