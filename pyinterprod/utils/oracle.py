# -*- coding: utf-8 -*-

import time
from typing import Callable, List, Optional, Sequence

from cx_Oracle import Cursor, DatabaseError

from pyinterprod import logger


def drop_mview(cur: Cursor, name: str):
    try:
        cur.execute(f"DROP MATERIALIZED VIEW {name}")
    except DatabaseError as exc:
        error, = exc.args

        # ORA-08103: object no longer exists
        # ORA-12003: materialized view does not exist
        if error.code not in (8103, 12003):
            raise exc


def drop_table(cur: Cursor, name: str, purge: bool=False):
    if purge:
        sql = f"DROP TABLE {name} PURGE"
    else:
        sql = f"DROP TABLE {name}"

    try:
        cur.execute(sql)
    except DatabaseError as exc:
        error, = exc.args

        # ORA-00942: table or view does not exist
        # ORA-08103: object no longer exists
        if error.code not in (942, 8103):
            raise exc


def gather_stats(cur: Cursor, schema: str, table: str, partition: Optional[str]=None):
    if partition:
        args = (schema, table, partition)
    else:
        args = (schema, table)

    cur.callproc("DBMS_STATS.GATHER_TABLE_STATS", args)


def get_child_tables(cur: Cursor, schema: str, name: str) -> List[tuple]:
    cur.execute(
        """
        SELECT TABLE_NAME, CONSTRAINT_NAME, COLUMN_NAME
        FROM ALL_CONS_COLUMNS
        WHERE OWNER = :schema
        AND CONSTRAINT_NAME IN (
          SELECT CONSTRAINT_NAME
          FROM ALL_CONSTRAINTS
          WHERE CONSTRAINT_TYPE = 'R'             -- Referential
          AND R_CONSTRAINT_NAME IN (
            SELECT CONSTRAINT_NAME
            FROM ALL_CONSTRAINTS
            WHERE CONSTRAINT_TYPE IN ('P', 'U')   -- Primary/Unique
            AND OWNER = :schema
            AND TABLE_NAME = :name
          )
        )
        """, dict(schema=schema, name=name)
    )

    return cur.fetchall()


def get_indexes(cur: Cursor, owner: str, name: str) -> List[dict]:
    cur.execute(
        """
        SELECT I.OWNER, I.INDEX_NAME, I.UNIQUENESS, I.TABLESPACE_NAME, 
               I.LOGGING, I.STATUS, IC.COLUMN_NAME, IC.DESCEND
        FROM ALL_INDEXES I
        INNER JOIN ALL_IND_COLUMNS IC
          ON I.OWNER = IC.INDEX_OWNER
          AND I.INDEX_NAME = IC.INDEX_NAME
          AND I.TABLE_NAME = IC.TABLE_NAME
        WHERE I.TABLE_OWNER = :1
        AND I.TABLE_NAME = :2
        ORDER BY I.INDEX_NAME, IC.COLUMN_POSITION
        """, (owner, name)
    )

    indexes = {}
    for row in cur:
        name = row[1]
        try:
            index = indexes[name]
        except KeyError:
            index = indexes[name] = {
                "owner": row[0],
                "name": name,
                "is_unique": row[2] == "UNIQUE",
                "tablespace": row[3],
                "logging": row[4] == "YES",
                "unusable": row[5] == "UNUSABLE",  # can be VALID or N/A
                "columns": []
            }

        index["columns"].append({
            "name": row[6],
            "order": row[7]
        })

    return list(indexes.values())


def get_partitions(cur: Cursor, schema: str, name: str) -> List[dict]:
    cur.execute(
        """
        SELECT P.PARTITION_NAME, P.PARTITION_POSITION, P.HIGH_VALUE,
               K.COLUMN_NAME, K.COLUMN_POSITION
        FROM ALL_TAB_PARTITIONS P
        INNER JOIN ALL_PART_KEY_COLUMNS K
          ON P.TABLE_OWNER = K.OWNER 
          AND P.TABLE_NAME = K.NAME
        WHERE P.TABLE_OWNER = :1 
        AND P.TABLE_NAME = :2
        """, (schema, name)
    )

    partitions = {}
    for row in cur:
        part_name = row[0]
        if part_name in partitions:
            raise ValueError("Multi-column partitioning keys not supported")

        partitions[part_name] = {
            "name": part_name,
            "position": row[1],
            "value": row[2],
            "column": row[3]
        }

    return sorted(partitions.values(), key=lambda x: x["position"])


def rebuild_index(cur: Cursor, name: str, parallel: bool=False):
    if parallel:
        cur.execute(f"ALTER INDEX {name} REBUILD PARALLEL")

        # Prevent the index to be accessed in parallel by default
        cur.execute(f"ALTER INDEX {name} NOPARALLEL")
    else:
        cur.execute(f"ALTER INDEX {name} REBUILD")


def toggle_constraint(cur: Cursor, table: str, constraint: str, enable: bool):
    if enable:
        sql = f"ALTER TABLE {table} ENABLE CONSTRAINT {constraint}"
    else:
        sql = f"ALTER TABLE {table} DISABLE CONSTRAINT {constraint}"

    cur.execute(sql)


def truncate_table(cur: Cursor, name: str, reuse_storage: bool=False):
    if reuse_storage:
        sql = f"TRUNCATE TABLE {name} REUSE STORAGE"
    else:
        sql = f"TRUNCATE TABLE {name}"

    cur.execute(sql)


def catch_temp_error(fn: Callable, args: Sequence, max_attempts: int=3):
    num_attempts = 0
    while True:
        try:
            fn(*args)
        except DatabaseError as exc:
            error, = exc.args
            num_attempts += 1

            # ORA-01652: unable to extend temp segment by <int> in tablespace TEMP
            if error.code != 1652 or num_attempts == max_attempts:
                raise exc
            else:
                logger.warning(exc)
                time.sleep(3600)
        else:
            break