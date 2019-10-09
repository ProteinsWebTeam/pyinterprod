# -*- coding: utf-8 -*-

import re
import subprocess
from typing import List, Optional, Tuple, Union

import cx_Oracle

from . import logger


INSERT_SIZE = 100000
DATA_PUMP_DIR = "PANDA_DATA_PUMP_DIR"


def copy_tables(user_src: str, user_dst: str, dsn: str):
    schema_src = user_src.split('/')[0]
    schema_dst = user_dst.split('/')[0]
    dumpfile = schema_src.upper()

    con = cx_Oracle.connect(make_connect_string(user_src, dsn))
    cur = con.cursor()
    tables = [f"{schema_src}.{t}" for t in get_tables(cur, schema_src)]
    cur.execute(f"UPDATE {schema_src}.CV_DATABASE SET IS_READY = 'N'")
    con.commit()
    cur.close()
    con.close()

    returncode = subprocess.call(["expdp", make_connect_string(user_src, dsn),
                                  f"TABLES={','.join(tables)}",
                                  f"DIRECTORY={DATA_PUMP_DIR}",
                                  f"DUMPFILE={dumpfile}.dmp",
                                  f"LOGFILE={dumpfile}-exp.log",
                                  "REUSE_DUMPFILES=YES"])

    if returncode:
        raise RuntimeError(f"expdp exited with code {returncode}")

    con = cx_Oracle.connect(make_connect_string(user_src, dsn))
    cur = con.cursor()
    cur.execute(f"UPDATE {schema_src}.CV_DATABASE SET IS_READY = 'Y'")
    con.commit()
    cur.close()
    con.close()

    drop_all(user_dst, dsn)
    returncode = subprocess.call(["impdp", make_connect_string(user_dst, dsn),
                                  f"TABLES={','.join(tables)}",
                                  f"DIRECTORY={DATA_PUMP_DIR}",
                                  f"DUMPFILE={dumpfile}.dmp",
                                  f"LOGFILE={dumpfile}-imp.log",
                                  f"REMAP_SCHEMA={schema_src}:{schema_dst}"])

    if returncode:
        raise RuntimeError(f"impdp exited with code {returncode}")

    con = cx_Oracle.connect(make_connect_string(user_dst, dsn))
    cur = con.cursor()
    cur.execute(f"UPDATE {schema_dst}.CV_DATABASE SET IS_READY = 'Y'")
    con.commit()
    cur.close()
    con.close()


def copy_table(user_src: str, user_dst: str, dsn: str, table: dict):
    owner_src = user_src.split('/')[0]
    owner_dst = user_dst.split('/')[0]

    con = cx_Oracle.connect(make_connect_string(user_dst, dsn))
    cur = con.cursor()

    query = f"CREATE TABLE {owner_dst}.{table['name']}"
    if table["partitions"]:
        column = table["partitions"][0]["column"]
        partitions = []
        for p in table["partitions"]:
            partitions.append(f"PARTITION {p['name']} VALUES ({p['value']})")

        query += f" PARTITION BY LIST ({column}) ({','.join(partitions)})"

    query += f" NOLOGGING AS SELECT * FROM {owner_src}.{table['name']}"
    cur.execute(query)

    for idx in table["indexes"]:
        columns = ','.join([
            f"{col['name']} {col['order']}"
            for col in idx["columns"]
        ])

        query = f"""
            CREATE {'UNIQUE' if idx['is_unique'] else ''} INDEX {idx['name']} 
            ON {owner_dst}.{table['name']} ({columns})
            {'LOGGING' if idx['logging'] else 'NOLOGGING'}
        """
        cur.execute(query)

    for c in table["constraints"]:
        if c["type"] == 'C' and not c["condition"].endswith("IS NOT NULL"):
            c_type = "CHECK"
            c_check = c["condition"]
        elif c["type"] == 'P':
            c_type = "PRIMARY KEY"
            c_check = ", ".join([col["name"] for col in c["columns"]])
        elif c["type"] == 'U':
            c_type = "UNIQUE"
            c_check = ", ".join([col["name"] for col in c["columns"]])
        else:
            continue

        query = f"""
            ALTER TABLE {owner_dst}.{table['name']}
            ADD CONSTRAINT {c['name']} {c_type} ({c_check})
        """

        if c["index_name"]:
            query += f"USING INDEX {c['index_name']}"

        cur.execute(query)

    for privilege, grantee in table["grants"]:
        grant(cur, owner_dst, table["name"], privilege, grantee)

    gather_stats(cur, owner_dst, table["name"])
    cur.close()
    con.close()


def drop_all(user: str, dsn: str):
    con = cx_Oracle.connect(make_connect_string(user, dsn))
    cur = con.cursor()
    owner = user.split('/')[0]
    for table in get_tables(cur, owner):
        drop_table(cur, owner, table, purge=True)
    cur.close()
    con.close()


def get_tables(cur: cx_Oracle.Cursor, owner: str) -> List[str]:
    cur.execute(
        """
        SELECT TABLE_NAME
        FROM ALL_TABLES
        WHERE UPPER(OWNER) = UPPER(:1)
        ORDER BY TABLE_NAME
        """, (owner,)
    )
    return [row[0] for row in cur]


def get_mviews(cur: cx_Oracle.Cursor, owner: str) -> List[str]:
    cur.execute(
        """
        SELECT MVIEW_NAME
        FROM ALL_MVIEWS
        WHERE UPPER(OWNER) = UPPER(:1)
        ORDER BY MVIEW_NAME
        """, (owner,)
    )
    return [row[0] for row in cur]


def parse_url(url: str) -> dict:
    m = re.match(r"(.+)/(.+)@(.+):(\d+)/([a-z]+)", url, re.I)
    return {
        "username": m.group(1),
        "password": m.group(2),
        "host": m.group(3),
        "port": int(m.group(4)),
        "service": m.group(5)
    }


def get_child_tables(cur: cx_Oracle.Cursor, owner: str,
                     table: str) -> List[dict]:
    cur.execute(
        """
        SELECT OWNER, TABLE_NAME, CONSTRAINT_NAME, COLUMN_NAME
        FROM ALL_CONS_COLUMNS
        WHERE OWNER = :owner
          AND CONSTRAINT_NAME IN (
              SELECT CONSTRAINT_NAME
              FROM ALL_CONSTRAINTS
              WHERE CONSTRAINT_TYPE = 'R'
                AND R_CONSTRAINT_NAME IN (
                  SELECT CONSTRAINT_NAME
                  FROM ALL_CONSTRAINTS
                  WHERE CONSTRAINT_TYPE IN ('P', 'U')
                    AND OWNER = :owner
                    AND TABLE_NAME = :tabname
          )
        )
        """,
        dict(owner=owner, tabname=table)
    )

    cols = ("owner", "name", "constraint", "column")
    return [dict(zip(cols, row)) for row in cur]


def toggle_constraint(cur: cx_Oracle.Cursor, owner: str, table: str,
                      constraint: str, enable: bool) -> bool:
    query = """
        ALTER TABLE {}.{}
        {} CONSTRAINT {}
    """.format(owner, table, "ENABLE" if enable else "DISABLE", constraint)

    try:
        cur.execute(query)
    except cx_Oracle.DatabaseError as e:
        logger.error(e)
        logger.error(query)
        return False
    else:
        return True


def get_grants(cur: cx_Oracle.Cursor, owner: str, table: str) -> List[Tuple[str, str]]:
    cur.execute(
        """
        SELECT PRIVILEGE, GRANTEE
        FROM DBA_TAB_PRIVS
        WHERE OWNER = UPPER(:1)
        AND TABLE_NAME = :2
        """, (owner, table)
    )
    return cur.fetchall()


def get_partitions(cur: cx_Oracle.Cursor, owner: str, table: str) -> List[dict]:
    cur.execute(
        """
        SELECT
          P.PARTITION_NAME, P.PARTITION_POSITION, P.HIGH_VALUE,
          K.COLUMN_NAME, K.COLUMN_POSITION
        FROM ALL_TAB_PARTITIONS P
        INNER JOIN ALL_PART_KEY_COLUMNS K
          ON P.TABLE_OWNER = K.OWNER AND P.TABLE_NAME = K.NAME
        WHERE P.TABLE_OWNER = UPPER(:1) AND P.TABLE_NAME = :2
        """,
        (owner, table)
    )
    partitions = {}
    for part_name, part_pos, part_val, col_name, col_pos in cur:
        if col_pos > 1:
            raise ValueError("Multi-column partitioning keys not supported")

        partitions[part_name] = {
            "name": part_name,
            "position": part_pos,
            "value": part_val,
            "column": col_name
        }

    return sorted(partitions.values(), key=lambda c: c["position"])


def drop_index(cur: cx_Oracle.Cursor, owner: str, name: str):
    try:
        cur.execute("DROP INDEX {}.{}".format(owner, name))
    except cx_Oracle.DatabaseError as exc:
        error, = exc.args

        # ORA-01418: specified index does not exist
        if error.code != 1418:
            raise exc


def drop_table(cur: cx_Oracle.Cursor, owner: str, name: str, purge: bool=False):
    try:
        purge = " PURGE" if purge else ""
        cur.execute("DROP TABLE {}.{}{}".format(owner, name, purge))
    except cx_Oracle.DatabaseError as exc:
        error, = exc.args

        # ORA-00942: table or view does not exist
        # ORA-08103: object no longer exists
        if error.code not in (942, 8103):
            raise exc


def drop_mview(cur: cx_Oracle.Cursor, owner: str, name: str):
    try:
        cur.execute("DROP MATERIALIZED VIEW {}.{}".format(owner, name))
    except cx_Oracle.DatabaseError as exc:
        error, = exc.args

        # ORA-08103: object no longer exists
        # ORA-12003: materialized view does not exist
        if error.code not in (8103, 12003):
            raise exc


def grant(cur: cx_Oracle.Cursor, owner: str, table: str, privilege: str,
          grantee: str):
    cur.execute("GRANT {} ON {}.{} TO {}".format(privilege, owner, table,
                                                 grantee))


def gather_stats(cur: cx_Oracle.Cursor, owner: str, table: str,
                 partition: Optional[str]=None, cascade: bool=False):
    if cascade:
        cur.execute(
            """
            BEGIN
                DBMS_STATS.GATHER_TABLE_STATS(:1, :2 :3, cascade=>TRUE);
            END;
            """, (owner, table, partition)
        )
    elif partition:
        cur.callproc("DBMS_STATS.GATHER_TABLE_STATS",
                     (owner, table, partition))
    else:
        cur.callproc("DBMS_STATS.GATHER_TABLE_STATS", (owner, table))


def truncate_table(cur: cx_Oracle.Cursor, owner: str, table: str,
                   reuse: bool=False):
    query = "TRUNCATE TABLE {}.{}".format(owner, table)
    if reuse:
        query += " REUSE STORAGE"

    cur.execute(query)


def exchange_partition(cur: cx_Oracle.Cursor, owner: str, src: str, dst: str,
                       partition: str, stats: bool=True):
    cur.execute(
        """
        ALTER TABLE {0}.{2}
        EXCHANGE PARTITION {3}
        WITH TABLE {0}.{1}
        INCLUDING INDEXES
        WITHOUT VALIDATION
        """.format(owner, src, dst, partition)
    )
    if stats:
        gather_stats(cur, owner, dst, partition)


def get_constraints(cur: cx_Oracle.Cursor, owner: str, table: str) -> List[dict]:
    cur.execute(
        """
        SELECT
          C.CONSTRAINT_NAME, C.CONSTRAINT_TYPE, C.TABLE_NAME,
          C.SEARCH_CONDITION, C.DELETE_RULE, C.STATUS,
          C.R_CONSTRAINT_NAME, C.INDEX_NAME, CC.COLUMN_NAME, CC.POSITION
        FROM ALL_CONSTRAINTS C
          INNER JOIN ALL_CONS_COLUMNS CC
            ON C.OWNER = CC.OWNER
            AND C.CONSTRAINT_NAME = CC.CONSTRAINT_NAME
            AND C.TABLE_NAME = CC.TABLE_NAME
        WHERE C.OWNER = UPPER(:1)
        AND C.TABLE_NAME = :2
        ORDER BY C.CONSTRAINT_NAME, CC.POSITION
        """, (owner, table)
    )

    constraints = {}
    for row in cur:
        name = row[0]
        if name in constraints:
            c = constraints[name]
        else:
            c = constraints[name] = {
                "name": name,
                "type": row[1],
                "table_name": row[2],
                "condition": row[3],
                "cascade": row[4] == "CASCADE",
                "is_enabled": row[5] == "ENABLED",
                "reference": row[6],
                "index_name": row[7],
                "columns": []
            }

        c["columns"].append({
            "name": row[8],
            "order": row[9]
        })

    return list(constraints.values())


def get_indices(cur: cx_Oracle.Cursor, owner: str, table: str) -> List[dict]:
    cur.execute(
        """
        SELECT
          I.OWNER, I.INDEX_NAME, I.TABLE_OWNER, I.TABLE_NAME,
          I.TABLESPACE_NAME, I.UNIQUENESS,
          I.LOGGING, IC.COLUMN_NAME, IC.DESCEND
        FROM ALL_INDEXES I
        INNER JOIN ALL_IND_COLUMNS IC
          ON I.OWNER = IC.INDEX_OWNER
          AND I.INDEX_NAME = IC.INDEX_NAME
          AND I.TABLE_NAME = IC.TABLE_NAME
        WHERE I.TABLE_OWNER = UPPER(:1)
        AND I.TABLE_NAME = :2
        ORDER BY I.INDEX_NAME, IC.COLUMN_POSITION
        """, (owner, table)
    )

    indices = {}
    for row in cur:
        name = row[1]
        if name in indices:
            idx = indices[name]
        else:
            idx = indices[name] = {
                "owner": row[0],
                "name": name,
                "table_owner": row[2],
                "table_name": row[3],
                "tablespace": row[4],
                "is_unique": row[5] == "UNIQUE",
                "logging": row[6] == "YES",
                "columns": []
            }

        idx["columns"].append({
            "name": row[7],
            "order": row[8]
        })

    return list(indices.values())


def recreate_index(cur: cx_Oracle.Cursor, index: dict):
    columns = ','.join(["{name} {order}".format(**col)
                        for col in index["columns"]])

    if index["is_unique"]:
        uniqueness = "UNIQUE"
    else:
        uniqueness = ""

    if index["tablespace"]:
        tablespace = "TABLESPACE " + index["tablespace"]
    else:
        tablespace = ""

    if index["logging"]:
        loginfo = "LOGGING"
    else:
        loginfo = "NOLOGGING"

    cur.execute(
        """
        CREATE {} INDEX {}.{}
        ON {}.{}({})
        {}
        {}
        """.format(
            uniqueness, index["owner"], index["name"],
            index["table_owner"], index["table_name"], columns,
            tablespace,
            loginfo
        )
    )


def rebuild_indices(cur: cx_Oracle.Cursor, owner: str, table: str):
    for idx in get_indices(cur, owner, table):
        logger.info("rebuilding {name}".format(**idx))

        if idx["logging"]:
            cur.execute("ALTER INDEX {name} REBUILD".format(**idx))
        else:
            cur.execute("ALTER INDEX {name} REBUILD NOLOGGING".format(**idx))


def make_connect_string(user: str, dsn: str) -> str:
    return user + '@' + dsn


def refresh_mview(url: str, name: str):
    con = cx_Oracle.connect(url)
    cur = con.cursor()

    try:
        cur.execute(
            """
            BEGIN
                DBMS_MVIEW.REFRESH(:1, method => '?');
            END;
            """, (name,)
        )
    finally:
        cur.close()
        con.close()


class TablePopulator(object):
    def __init__(self, con: cx_Oracle.Connection, query: str,
                 buffer_size: int=100000, autocommit: bool=False):
        self.con = con
        self.cur = con.cursor()
        self.query = query
        self.buffer_size = buffer_size
        self.autocommit = autocommit
        self.records = []
        self.rowcount = 0

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __del__(self):
        self.close()

    def _execute(self, record: Union[dict, tuple]):
        self.records.append(record)
        self.rowcount += 1

        if len(self.records) == self.buffer_size:
            self.flush()

    def insert(self, record: Union[dict, tuple]):
        self._execute(record)

    def update(self, record: Union[dict, tuple]):
        self._execute(record)

    def delete(self, record: Union[dict, tuple]):
        self._execute(record)

    def flush(self):
        if not self.records:
            return

        self.cur.executemany(self.query, self.records)
        self.records = []

        if self.autocommit:
            self.con.commit()

    def close(self):
        if self.cur is not None:
            self.flush()
            self.cur.close()
            self.cur = None

    def open(self):
        self.close()
        self.cur = self.con.cursor()
        self.rowcount = 0
