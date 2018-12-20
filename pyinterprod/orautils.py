import logging
from typing import List, Tuple

import cx_Oracle


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


def toggle_constraint(cur: cx_Oracle.Cursor, owner: str, table: str,
                      constraint: str, enable: bool=True) -> bool:
    query = """
        ALTER TABLE {}.{}
        {} CONSTRAINT {}
    """.format(owner, table, "ENABLE" if enable else "DISABLE", constraint)

    try:
        cur.execute(query)
    except cx_Oracle.DatabaseError as e:
        logging.error(e)
        logging.error(query)
        return False
    else:
        return True


def toggle_constraints(cur: cx_Oracle.Cursor, owner: str, table: str,
                       enable: bool=True) -> List[Tuple[str, bool]]:
    cur.execute(
        """
        SELECT CONSTRAINT_NAME, STATUS
        FROM USER_CONSTRAINTS
        WHERE OWNER = :1 AND TABLE_NAME = :2
        """,
        (owner, table)
    )
    constraints = dict(cur.fetchall())

    results = []
    for c, s in constraints.items():
        is_enabled = s == "ENABLED"
        if is_enabled != enable:
            success = toggle_constraint(cur, owner, table, c, enable)
        else:
            # Already has the desired status (enabled/disabled)
            success = True

        results.append(c, success)

    return results


def get_partitions(cur: cx_Oracle.Cursor, owner: str, table: str) -> list:
    cur.execute(
        """
        SELECT HIGH_VALUE
        FROM ALL_TAB_PARTITIONS
        WHERE TABLE_OWNER = :1
        AND TABLE_NAME = :2
        """,
        (owner, table)
    )
    partitions = [row[0].strip('\'') for row in cur]
    return partitions
