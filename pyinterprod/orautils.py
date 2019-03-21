import re
from typing import List, Optional

import cx_Oracle

from . import logger


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
                partition: Optional[str]=None):
    con = cx_Oracle.connect(url)
    cur = con.cursor()

    if partition:
        dst = "{} PARTITION ({})".format(table, partition)
        name = "{} ({})".format(table, partition)
    else:
        dst = "{}".format(table)
        name = table

    cur.execute(
        """
        SELECT COUNT(*)
        FROM INTERPRO.{}
        WHERE {} IN (
          SELECT PROTEIN_AC 
          FROM INTERPRO.PROTEIN_TO_DELETE
        )
        """.format(dst, column)
    )
    num_rows = cur.fetchone()[0]

    if num_rows:
        for i in range(0, stop, step):
            cur.execute(
                """
                DELETE FROM INTERPRO.{}
                WHERE {} IN (
                  SELECT PROTEIN_AC
                  FROM INTERPRO.PROTEIN_TO_DELETE
                  WHERE ID BETWEEN :1 and :2
                )
                """.format(dst, column),
                (i, i + step - 1)
            )
            logger.debug("{}: {} / {}".format(name, min(i + step, stop), stop))

        con.commit()

    cur.close()
    con.close()


def create_db_link(cur: cx_Oracle.Cursor, link: str, user: str, passwd: str,
                   connect_string: str):
    try:
        cur.execute("DROP PUBLIC DATABASE LINK {}".format(link))
    except cx_Oracle.DatabaseError as e:
        if e.args[0].code == 2024:
            """
            ORA-02024: database link not found
            that's fine, we're about to create the link
            """
            pass
        else:
            raise e

    cur.execute(
        """
        CREATE PUBLIC DATABASE LINK {}
        CONNECT TO {} IDENTIFIED BY {}
        USING '{}'
        """.format(link, user, passwd, connect_string)
    )


def get_indexes(cur: cx_Oracle.Cursor, owner: str, table: str):
    cur.execute(
        """
        SELECT INDEX_NAME
        FROM ALL_INDEXES
        WHERE TABLE_OWNER = :1
        AND TABLE_NAME = :2
        """,
        (owner, table)
    )
    return [row[0] for row in cur]


def create_db_links(url: str, urls: List[str]):
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    for obj in map(parse_url, urls):
        if obj["service"].startswith(('V', 'v')):
            link = obj["service"][1:]
        else:
            link = obj["service"]

        create_db_link(cur,
                       link=link,
                       user=obj["username"],
                       passwd=obj["password"],
                       connect_string="{}:{}/{}".format(obj["host"],
                                                        obj["port"],
                                                        obj["service"])
                       )

    cur.close()
    con.close()


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
                 buffer_size: int=10000, autocommit: bool=False):
        self.con = con
        self.cur = con.cursor()
        self.query = query
        self.buffer_size = buffer_size
        self.autocommit = autocommit
        self.records = []

    def __del__(self):
        self.close()

    def insert(self, record: tuple):
        self.records.append(record)

        if len(self.records) == self.buffer_size:
            self.flush(commit=self.autocommit)

    def flush(self, commit: bool=False):
        if not self.records:
            return

        self.cur.executemany(self.query, self.records)
        self.records = []

        if commit:
            self.commit()

    def commit(self):
        self.con.commit()

    def close(self):
        if self.cur is not None:
            self.flush(commit=True)
            self.cur.close()
            self.cur = None






