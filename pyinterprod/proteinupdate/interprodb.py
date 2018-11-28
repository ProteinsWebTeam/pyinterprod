#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
from concurrent import futures
from typing import Generator, List

import cx_Oracle

from .io import ProteinDatabase


_MAX_ITEMS = 100000


def get_proteins(url: str) -> Generator:
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.execute(
        """
        SELECT
          NAME, PROTEIN_AC, DBCODE, LEN, FRAGMENT, TAX_ID, CRC64
        FROM INTERPRO.PROTEIN
        """
    )

    for row in cur:
        yield (
            row[0],
            row[1],
            1 if row[2] == 'S' else 0,
            row[3],
            1 if row[4] == 'Y' else 0,
            row[5],
            row[6]
        )

    cur.close()
    con.close()


def update_proteins(url: str, db: ProteinDatabase) -> int:
    con = cx_Oracle.connect(url)
    cur = con.cursor()

    count = 0
    for row in db.get_changes():
        cur.execute(
            """
            UPDATE INTERPRO.PROTEIN
            SET 
              NAME = :1,
              DBCODE = :2,
              CRC64 = :3,
              LEN = :4,
              FRAGMENT = :5,
              TAX_ID = :6,
              TIMESTAMP = SYSDATE,
              USERSTAMP = USER
              WHERE PROTEIN_AC = :7
            """,
            (
                row[1],
                'S' if row[2] else 'T',
                row[3],
                row[4],
                'Y' if row[5] else 'N',
                row[6],
                row[0]
            )
        )
        count += 1

    con.commit()
    cur.close()
    con.close()
    return count


def insert_proteins(url: str, db: ProteinDatabase) -> int:
    con = cx_Oracle.connect(url)
    cur = con.cursor()

    items = []
    count = 0
    for row in db.get_new():
        items.append((
            row[0],                     # accession
            row[1],                     # identifier
            'S' if row[2] else 'T',     # dbcode
            row[3],                     # crc64
            row[4],                     # length
            'Y' if row[5] else 'N',     # sequence status (fragment)
            row[6]                      # taxon ID
        ))
        count += 1

        if len(items) == _MAX_ITEMS:
            cur.executemany(
                """
                INSERT INTO INTERPRO.PROTEIN (
                  PROTEIN_AC, NAME, DBCODE, CRC64, LEN, FRAGMENT, 
                  STRUCT_FLAG, TAX_ID
                )
                VALUES (:1, :2, :3, :4, :5, :6, 'N', :7)
                """,
                items
            )
            items = []

    con.commit()
    cur.close()
    con.close()

    return count


def delete_proteins(url: str, table: str, column: str, step: int=1000):
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.execute()
    cur.execute(
        """
        SELECT COUNT(*)
        FROM INTERPRO.{}
        WHERE {} IN (
            SELECT PROTEIN_AC
            FROM INTERPRO.PROTEIN_TO_DELETE
        )
        """.format(table, column)
    )
    stop = cur.fetchone()[0]

    for i in range(0, stop, step):
        cur.execute(
            """
            DELETE FROM INTERPRO.{}
            WHERE {} IN (
              SELECT PROTEIN_AC
              FROM INTERPRO.PROTEIN_TO_DELETE
              WHERE ID BETWEEN :1 and :2
            )
            """.format(table, column),
            (i, i+step-1)
        )
        logging.info("{}: {} / {}".format(table, min(i+step, stop), stop))
    con.commit()
    cur.close()
    con.close()


def prepare_deletion(url: str, db: ProteinDatabase) -> int:
    con = cx_Oracle.connect(url)
    cur = con.cursor()

    try:
        cur.execute("INTERPRO.PROTEIN_TO_DELETE")
    except cx_Oracle.DatabaseError:
        pass
    finally:
        cur.execute(
            """
            CREATE TABLE INTERPRO.PROTEIN_TO_DELETE
            (
                ID NUMBER NOT NULL,
                PROTEIN_AC VARCHAR2(15) NOT NULL
            ) NOLOGGING
            """
        )

    count = 0
    items = []
    for accession in db.get_deleted():
        items.append((count, accession))
        count += 1

        if len(items) == _MAX_ITEMS:
            cur.executemany(
                """
                INSERT /*+APPEND*/ INTO INTERPRO.PROTEIN_TO_DELETE
                VALUES (:1, :2)
                """,
                items
            )
            con.commit()
            items = []

    if items:
        cur.executemany(
            """
            INSERT /*+APPEND*/ INTO INTERPRO.PROTEIN_TO_DELETE
            VALUES (:1, :2)
            """,
            items
        )
        con.commit()

    cur.execute(
        """
        ALTER TABLE INTERPRO.PROTEIN_TO_DELETE
        ADD CONSTRAINT PK_PROTEIN_TO_DELETE
        PRIMARY KEY (ID)
        """
    )

    cur.exec("DBMS_STATS.GATHER_TABLE_STATS",
             ("INTERPRO", "PROTEIN_TO_DELETE"))

    cur.close()
    con.close()
    return count


def count_rows_to_delete(url: str, table: str, column: str) -> int:
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.execute(
        """
        SELECT COUNT(*)
        FROM INTERPRO.{}
        WHERE {} IN (
            SELECT PROTEIN_AC
            FROM INTERPRO.PROTEIN_TO_DELETE
        )
        """.format(table, column)
    )
    cnt = cur.fetchone()[0]
    cur.close()
    con.close()
    return cnt


def get_child_tables(url: str, owner: str, table: str) -> List[dict]:
    con = cx_Oracle.connect(url)
    cur = con.cursor()
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
    tables = [dict(zip(cols, row)) for row in cur]

    cur.close()
    con.close()
    return tables


def toggle_constraint(url: str, owner: str, table: str, constraint: str,
                      enable: bool=True):
    con = cx_Oracle.connect(url)
    cur = con.cursor()

    try:
        cur.execute(
            """
            ALTER TABLE {}.{} {} CONSTRAINT {}
            """.format(
                owner, table,
                "ENABLE" if enable else "DISABLE",
                constraint
            )
        )
    except cx_Oracle.DatabaseError as e:
        logging.critical("failed: ALTER TABLE {}.{} {} CONSTRAINT {}".format(
            owner, table,
            "ENABLE" if enable else "DISABLE",
            constraint
        ))
        raise e
    finally:
        cur.close()
        con.close()


def get_tables_with_proteins_to_delete(url: str) -> List[dict]:
    tables = get_child_tables(url, "INTERPRO", "PROTEIN")
    tables.append({
        "owner": "INTERPRO",
        "name": "PROTEIN",
        "column": "PROTEIN_AC"
    })

    n = len(tables)
    tables_to_process = []
    with futures.ThreadPoolExecutor(max_workers=n) as executor:
        future_to_idx = {}
        for i, t in enumerate(tables):
            future = executor.submit(count_rows_to_delete, url,
                                     t["name"], t["column"])
            future_to_idx[future] = i

        done, not_done = futures.wait(future_to_idx)
        if not_done:
            raise RuntimeError("error in table(s): ".format(
                ", ".join([
                    tables[future_to_idx[f]]["name"]
                    for f in not_done]
                )
            ))

        for future in done:
            if future.result():
                i = future_to_idx[future]
                tables_to_process.append(tables[i])

    return tables_to_process
