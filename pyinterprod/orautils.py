import re
from typing import List, Optional

import cx_Oracle

from . import logger


def parse_url(url: str) -> dict:
    m = re.match(r"(.+)/(.+)@(.+):(\d+)/[a-z]+", url, re.I)
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
                      constraint: str, enable: bool=True) -> bool:
    cur.execute(
        """
        SELECT STATUS
        FROM USER_CONSTRAINTS
        WHERE OWNER = :1
        AND TABLE_NAME = :2
        AND CONSTRAINT_NAME = :3
        """,
        (owner, table, constraint)
    )
    row = cur.fetchone()
    if not row:
        logger.warning("{} unknown".format(constraint))
        return False

    is_enabled = row[0] == "ENABLE"
    if is_enabled == enable:
        logger.info("skipping {}".format(constraint))
        return True  # Already with the desired status

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


def get_partitions(cur: cx_Oracle.Cursor, owner: str, table: str) -> list:
    cur.execute(
        """
        SELECT PARTITION_NAME
        FROM ALL_TAB_PARTITIONS
        WHERE TABLE_OWNER = :1
        AND TABLE_NAME = :2
        """,
        (owner, table)
    )
    return [row[0].strip('\'') for row in cur]


def delete_iter(url: str, table: str, column: str, stop: int, step: int,
                partition: Optional[str] = None):
    if partition:
        query = """
            DELETE FROM INTERPRO.{} PARTITION ({})
            WHERE {} IN (
              SELECT PROTEIN_AC
              FROM INTERPRO.PROTEIN_TO_DELETE
              WHERE ID BETWEEN :1 and :2
            )
        """.format(table, partition, column)
        table += " ({})".format(partition)
    else:
        query = """
            DELETE FROM INTERPRO.{}
            WHERE {} IN (
              SELECT PROTEIN_AC
              FROM INTERPRO.PROTEIN_TO_DELETE
              WHERE ID BETWEEN :1 and :2
            )
        """.format(table, column)

    con = cx_Oracle.connect(url)
    cur = con.cursor()
    for i in range(0, stop, step):
        cur.execute(query, (i, i+step-1))
        logger.info("{}: {} / {}".format(table, min(i + step, stop), stop))

    con.commit()
    cur.close()
    con.close()


def create_db_link(cur: cx_Oracle.Cursor, link: str, user: str, passwd: str,
                   connect_string: str):
    cur.execute(
        """
        CREATE PUBLIC DATABASE LINK {}
        CONNECT TO {} IDENTIFIED BY {}
        USING '{}'
        """.format(link, user, passwd, connect_string)
    )
