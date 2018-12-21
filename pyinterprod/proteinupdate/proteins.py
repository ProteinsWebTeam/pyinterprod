import os
import sqlite3
from concurrent import futures
from tempfile import mkstemp
from typing import Generator, List, Optional, Tuple

import cx_Oracle

from . import sprot
from .. import logger, orautils

_MAX_ITEMS = 100000


class ProteinDatabase(object):
    def __init__(self, path: Optional[str] = None, dir: Optional[str] = None):
        if path:
            self.path = path
            self.temporary = False
        else:
            fd, self.path = mkstemp(dir=dir)
            os.close(fd)
            os.remove(self.path)
            self.temporary = True

        self._create_table("protein_old")
        self._create_table("protein_new")

    def __del__(self):
        if self.temporary:
            self.drop()

    @property
    def size(self) -> int:
        return os.path.getsize(self.path)

    def drop(self):
        try:
            os.remove(self.path)
        except FileNotFoundError:
            pass

    def _create_table(self, table_name: str):
        with sqlite3.connect(self.path) as con:
            con.execute(
                """
                CREATE TABLE IF NOT EXISTS {} (
                  accession TEXT NOT NULL PRIMARY KEY,
                  identifier TEXT NOT NULL,
                  is_reviewed INTEGER NOT NULL,
                  crc64 TEXT NOT NULL,
                  length INTEGER NOT NULL,
                  is_fragment INTEGER NOT NULL,
                  taxon_id INTEGER NOT NULL
                )
                """.format(table_name)
            )

    def _insert(self, src: Generator, table_name: str, max_items: int) -> int:
        query = """
            INSERT INTO {}
            VALUES (?, ?, ?, ?, ?, ?, ?)
        """.format(table_name)

        with sqlite3.connect(self.path) as con:
            count = 0
            items = []

            for protein in src:
                items.append(protein)
                count += 1

                if len(items) == max_items:
                    con.executemany(query, items)
                    items = []

            if items:
                con.executemany(query, items)
            con.commit()

        return count

    def insert_new(self, src: Generator, max_items: int=_MAX_ITEMS) -> int:
        return self._insert(src, "protein_new", max_items)

    def insert_old(self, src: Generator, max_items: int=_MAX_ITEMS) -> int:
        return self._insert(src, "protein_old", max_items)

    def _iter(self, table_name: str) -> Generator[Tuple, None, None]:
        with sqlite3.connect(self.path) as con:
            for row in con.execute("SELECT * FROM {}".format(table_name)):
                yield row

    def iter_new(self) -> Generator[Tuple, None, None]:
        return self._iter("protein_new")

    def iter_old(self) -> Generator[Tuple, None, None]:
        return self._iter("protein_old")

    def get_deleted(self) -> Generator[str, None, None]:
        with sqlite3.connect(self.path) as con:
            cur = con.execute(
                """
                SELECT accession
                FROM protein_old
                EXCEPT
                SELECT accession
                FROM protein_new
                """
            )

            for row in cur:
                yield row[0]

    def get_new(self) -> Generator[str, None, None]:
        with sqlite3.connect(self.path) as con:
            cur = con.execute(
                """
                SELECT
                  accession,
                  identifier,
                  is_reviewed,
                  crc64,
                  length,
                  is_fragment,
                  taxon_id
                FROM protein_new
                WHERE accession NOT IN (
                  SELECT accession
                  FROM protein_old
                )
                """
            )

            for row in cur:
                yield row

    def get_changes(self) -> Generator[Tuple, None, None]:
        with sqlite3.connect(self.path) as con:
            cur = con.execute(
                """
                SELECT
                  accession, p1.identifier, p1.is_reviewed, p1.crc64,
                  p1.length, p1.is_fragment, p1.taxon_id
                FROM protein_new AS p1
                INNER JOIN protein_old AS p2
                  USING (accession)
                WHERE p1.identifier != p2.identifier
                  OR p1.is_reviewed != p2.is_reviewed
                  OR p1.crc64 != p2.crc64
                  OR p1.length != p2.length
                  OR p1.is_fragment != p2.is_fragment
                  OR p1.taxon_id != p2.taxon_id

                """
            )

            for row in cur:
                yield row


def _get_proteins(url: str) -> Generator:
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.execute(
        """
        SELECT
          PROTEIN_AC, NAME, DBCODE, CRC64, LEN, FRAGMENT, TAX_ID
        FROM INTERPRO.PROTEIN
        """
    )

    for row in cur:
        yield (
            row[0],
            row[1],
            1 if row[2] == 'S' else 0,
            row[3],
            row[4],
            1 if row[5] == 'Y' else 0,
            row[6]
        )

    cur.close()
    con.close()


def _update_proteins(url: str, db: ProteinDatabase) -> int:
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


def _insert_proteins(url: str, db: ProteinDatabase) -> int:
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
            # TIMESTAMP and USERSTAMP will be set to their DEFAULT
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

    if items:
        # TIMESTAMP and USERSTAMP will be set to their DEFAULT
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

    con.commit()
    cur.close()
    con.close()

    return count


def _delete_proteins(url: str, table: str, column: str, stop: int,
                     step: int=_MAX_ITEMS):
    con = cx_Oracle.connect(url)
    cur = con.cursor()
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
        logger.info("{}: {} / {}".format(table, min(i+step, stop), stop))
    con.commit()
    cur.close()
    con.close()


def _prepare_deletion(url: str, db: ProteinDatabase) -> int:
    con = cx_Oracle.connect(url)
    cur = con.cursor()

    try:
        cur.execute("DROP TABLE INTERPRO.PROTEIN_TO_DELETE")
    except cx_Oracle.DatabaseError:
        pass
    finally:
        cur.execute(
            """
            CREATE TABLE INTERPRO.PROTEIN_TO_DELETE
            (
                ID NUMBER NOT NULL,
                PROTEIN_AC VARCHAR2(15) NOT NULL
            ) NOlogger
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

    cur.callproc("DBMS_STATS.GATHER_TABLE_STATS",
                 ("INTERPRO", "PROTEIN_TO_DELETE"))

    cur.close()
    con.close()
    return count


def _count_rows_to_delete(url: str, table: str, column: str) -> int:
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


def _count_proteins_to_delete(cur: cx_Oracle.Cursor) -> int:
    cur.execute(
        """
        SELECT COUNT(*)
        FROM INTERPRO.PROTEIN_TO_DELETE
        """
    )
    return cur.fetchone()[0]


def _get_tables_with_proteins_to_delete(url: str) -> List[dict]:
    tables = orautils.get_child_tables(url, "INTERPRO", "PROTEIN")
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
            future = executor.submit(_count_rows_to_delete, url,
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


def _exchange_match_partition(cur: cx_Oracle.Cursor, dbcode: str):
    cur.execute(
        """
        ALTER TABLE INTERPRO.MATCH
        EXCHANGE PARTITION MATCH_DBCODE_{0}
        WITH TABLE INTERPRO.MATCH_TMP_{0}
        WITHOUT VALIDATION
        UPDATE GLOBAL INDEXES
        """.format(dbcode)
    )


def _export_match_partition(url: str, dbcode: str):
    table = "MATCH_TMP_" + dbcode
    con = cx_Oracle.connect(url)
    cur = con.cursor()

    try:
        cur.execute("DROP TABLE INTERPRO." + table)
    except cx_Oracle.DatabaseError:
        pass

    cur.execute(
        """
        CREATE TABLE INTERPRO.{} AS
        SELECT *
        FROM INTERPRO.MATCH PARTITION (MATCH_DBCODE_{})
        WHERE PROTEIN_AC NOT IN (
          SELECT PROTEIN_AC
          FROM INTERPRO.PROTEIN_TO_DELETE
        )
        """.format(table, dbcode)
    )

    con.commit()
    cur.close()
    con.close()


def track_changes(url: str, swissprot_path: str, trembl_path: str,
                  dir: Optional[str]=None):
    if dir:
        os.makedirs(dir, exist_ok=True)

    logger.info("loading proteins")

    db = ProteinDatabase(dir=dir)
    count = db.insert_old(_get_proteins(url))
    logger.info("InterPro: {} proteins".format(count))

    count = sprot.load(swissprot_path, db.path, "protein_new")
    logger.info("Swiss-Prot: {} proteins".format(count))

    count = sprot.load(trembl_path, db.path, "protein_new")
    logger.info("TrEMBL: {} proteins".format(count))

    logger.info("disk space used: {} bytes".format(db.size))

    count = _update_proteins(url, db)
    logger.info("{} proteins updated".format(count))

    count = _insert_proteins(url, db)
    logger.info("{} proteins added".format(count))

    count = _prepare_deletion(url, db)
    logger.info("{} proteins to delete".format(count))

    db.drop()


def delete(url: str, truncate: bool=False, refresh_partitions: bool=False):
    tables = _get_tables_with_proteins_to_delete(url)

    if not tables:
        return

    con = cx_Oracle.connect(url)
    cur = con.cursor()

    if truncate:
        logger.info("truncating MV/*NEW tables")
        _tables = []

        for t in tables:
            if t["name"].startswith("MV_") or t["name"].endswith("_NEW"):
                cur.execute("TRUNCATE TABLE INTERPRO.{}".format( t["name"]))
            else:
                _tables.append(t)

        tables = _tables

    logger.info("disabling referential constraints")
    table_constraints = []
    success = True
    for t in tables:
        try:
            contraint = t["constraint"]
        except KeyError:
            # The table does not have a constraint to disable
            continue
        else:
            table_constraints.append(t)
            if not orautils.toggle_constraint(cur, t["owner"], t["name"],
                                              contraint, False):
                success = False

    if not success:
        cur.close()
        con.close()
        raise RuntimeError("one or more constraints could not be disabled")

    with futures.ThreadPoolExecutor() as executor:
        fs = {}

        if refresh_partitions:
            idx = None
            for i, t in enumerate(tables):
                if t["name"] == "MATCH":
                    idx = i
                    break

            if idx is not None:
                # Disable all MATCH constraints
                r = orautils.toggle_constraints(cur, "INTERPRO", "MATCH",
                                                False)
                if not all([s for c, s in r]):
                    cur.close()
                    con.close()
                    raise RuntimeError("one or more constraints "
                                       "could not be disabled")

                tables.pop(idx)
                dbcodes = orautils.get_partitions(cur, "INTERPRO", "MATCH")
                for dbcode in dbcodes:
                    f = executor.submit(_export_match_partition, url, dbcode)
                    fs[f] = (-1, dbcode)

        count = _count_proteins_to_delete(cur)
        for i, t in enumerate(tables):
            f = executor.submit(_delete_proteins, url, t["name"],
                                t["column"], count)
            fs[f] = (i, None)

        success = True
        dbcodes = []
        for f in futures.as_completed(fs):
            i, dbcode = fs[f]

            if i == -1:
                if f.exception() is None:
                    logger.info("partition for '{}' done".format(dbcode))
                    dbcodes.append(dbcode)
                else:
                    logger.error("partition for '{}' exited".format(dbcode))
                    success = False
            else:
                t = tables[i]
                if f.exception() is None:
                    logger.info("table '{}' done".format(t["name"]))
                else:
                    logger.error("table '{}' exited".format(t["name"]))
                    success = False

        if not success:
            cur.close()
            con.close()
            raise RuntimeError("one or more tasks were not completed")

        if dbcodes:
            for dbcode in dbcodes:
                logger.info("exchanging partition for '{}'".format(dbcode))
                _exchange_match_partition(cur, dbcode)

            # TODO: uncomment, or rebuild indexes before enabling constraints
            # r = orautils.toggle_constraints(cur, "INTERPRO", "MATCH",
            #                                 True)

        success = True
        for t in table_constraints:
            if not orautils.toggle_constraint(cur, t["owner"], t["name"],
                                              t["constraint"], True):
                success = False

        cur.close()
        con.close()

        if success:
            logger.info("complete")
        else:
            raise RuntimeError("one or more constraints could not be enabled")
