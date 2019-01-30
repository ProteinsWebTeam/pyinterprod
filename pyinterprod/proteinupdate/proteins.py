import datetime
import os
import sqlite3
from concurrent import futures
from tempfile import mkstemp
from typing import Generator, Optional, Tuple

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

    def get_annotation_changes(self) -> Generator[Tuple, None, None]:
        with sqlite3.connect(self.path) as con:
            cur = con.execute(
                """
                SELECT
                  accession, p1.identifier, p1.is_reviewed, p1.crc64,
                  p1.length, p1.is_fragment, p1.taxon_id
                FROM protein_new AS p1
                INNER JOIN protein_old AS p2
                  USING (accession)
                WHERE p1.crc64 = p2.crc64
                AND (
                     p1.identifier != p2.identifier
                  OR p1.is_reviewed != p2.is_reviewed
                  OR p1.crc64 != p2.crc64
                  OR p1.length != p2.length
                  OR p1.is_fragment != p2.is_fragment
                  OR p1.taxon_id != p2.taxon_id
                )

                """
            )

            for row in cur:
                yield row

    def get_sequence_changes(self) -> Generator[Tuple, None, None]:
        with sqlite3.connect(self.path) as con:
            cur = con.execute(
                """
                SELECT
                  accession, p1.identifier, p1.is_reviewed, p1.crc64,
                  p1.length, p1.is_fragment, p1.taxon_id
                FROM protein_new AS p1
                INNER JOIN protein_old AS p2
                  USING (accession)
                WHERE p1.crc64 != p2.crc64
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


def _update_proteins(con: cx_Oracle.Connection, db: ProteinDatabase) -> int:
    cur = con.cursor()

    count = 0
    for row in db.get_annotation_changes():
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

    items = []
    for row in db.get_sequence_changes():
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

        items.append((row[0],))
        if len(items) == _MAX_ITEMS:
            cur.executemany(
                """
                INSERT INTO INTERPRO.PROTEIN_CHANGES (NEW_PROTEIN_AC)
                VALUES (:1)
                """,
                items
            )
            items = []

    if items:
        cur.executemany(
            """
            INSERT INTO INTERPRO.PROTEIN_CHANGES (NEW_PROTEIN_AC)
            VALUES (:1)
            """,
            items
        )

    con.commit()
    cur.close()
    return count


def _insert_proteins(con: cx_Oracle.Connection, db: ProteinDatabase) -> int:
    cur = con.cursor()

    items = []
    items2 = []
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

        items2.append((row[0],))
        if len(items2) == _MAX_ITEMS:
            cur.executemany(
                """
                INSERT INTO INTERPRO.PROTEIN_CHANGES (NEW_PROTEIN_AC)
                VALUES (:1)
                """,
                items2
            )
            items2 = []

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

    if items2:
        cur.executemany(
            """
            INSERT INTO INTERPRO.PROTEIN_CHANGES (NEW_PROTEIN_AC)
            VALUES (:1)
            """,
            items2
        )

    con.commit()
    cur.close()

    return count


def _prepare_deletion(con: cx_Oracle.Connection, db: ProteinDatabase) -> int:
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

    cur.callproc("DBMS_STATS.GATHER_TABLE_STATS",
                 ("INTERPRO", "PROTEIN_TO_DELETE"))

    cur.close()
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


def _init_protein_changes(con: cx_Oracle.Connection):
    cur = con.cursor()
    cur.execute("TRUNCATE TABLE INTERPRO.PROTEIN_CHANGES")
    cur.close()


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

    con = cx_Oracle.connect(url)
    _init_protein_changes(con)

    count = _update_proteins(con, db)
    logger.info("{} proteins updated".format(count))

    count = _insert_proteins(con, db)
    logger.info("{} proteins added".format(count))

    count = _prepare_deletion(con, db)
    logger.info("{} proteins to delete".format(count))

    con.close()

    db.drop()


def delete(url: str, truncate: bool=False):
    con = cx_Oracle.connect(url)
    cur = con.cursor()

    tables = orautils.get_child_tables(cur, "INTERPRO", "PROTEIN")
    tables.append({
        "owner": "INTERPRO",
        "name": "PROTEIN",
        "constraint": None,
        "column": "PROTEIN_AC"
    })

    if truncate:
        logger.info("truncating MV*/*NEW tables")
        _tables = []

        for t in tables:
            if t["name"].startswith("MV_") or t["name"].endswith("_NEW"):
                cur.execute("TRUNCATE TABLE INTERPRO.{}".format(t["name"]))
            else:
                _tables.append(t)

        tables = _tables

    logger.info("disabling referential constraints")
    ok = True
    for t in tables:
        if t["constraint"] and not orautils.toggle_constraint(
                cur, t["owner"], t["name"], t["constraint"], False
        ):
            ok = False

    if not ok:
        cur.close()
        con.close()
        raise RuntimeError("one or more constraints could not be disabled")

    count = _count_proteins_to_delete(cur)

    with futures.ThreadPoolExecutor() as executor:
        fs = {}

        for t in tables:
            if t["name"] == "MATCH":
                # for MATCH table: delete by partition
                dbcodes = orautils.get_partitions(cur, t["owner"], t["name"])
                for dbcode in dbcodes:
                    p = dict(t)  # shallow copy
                    p["partition"] = dbcode
                    f = executor.submit(orautils.delete_iter, url, t["name"],
                                        t["column"], count, _MAX_ITEMS,
                                        dbcode)
                    fs[f] = p
            else:
                f = executor.submit(orautils.delete_iter, url, t["name"],
                                    t["column"], count, _MAX_ITEMS)
                fs[f] = t

        tables = []
        errors = 0
        for f in futures.as_completed(fs):
            t = fs[f]

            if t.get("partition"):
                name = "{} ({})".format(t["name"], t["partition"])
            else:
                name = t["name"]

            if f.exception() is None:
                logger.info("table '{}' done".format(name))
                tables.append(t)
            else:
                logger.info("table '{}' exited".format(name))
                errors += 0

        if errors:
            cur.close()
            con.close()
            raise RuntimeError("some tables were not processed")
        else:
            error = False
            constraints = set()
            for t in tables:
                if not t["constraint"] or t["constraint"] in constraints:
                    continue

                constraints.add(t["constraint"])
                logger.info("enabling: {}.{}.{}".format(t["owner"],
                                                        t["name"],
                                                        t["constraint"]))

                if not orautils.toggle_constraint(cur, t["owner"], t["name"],
                                                  t["constraint"], True):
                    error = True

            cur.close()
            con.close()
            if error:
                raise RuntimeError("One or more constraints "
                                   "could not be enabled")


def update(url: str, version: str, date: datetime.datetime):
    delete(url, truncate=True)

    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.execute("SELECT COUNT(*) FROM INTERPRO.PROTEIN WHERE DBCODE = 'S'")
    n_swissprot = cur.fetchone()[0]

    cur.execute("SELECT COUNT(*) FROM INTERPRO.PROTEIN WHERE DBCODE = 'T'")
    n_trembl = cur.fetchone()[0]

    cur.execute(
        """
        UPDATE INTERPRO.DB_VERSION
        SET
          VERSION = :1,
          ENTRY_COUNT = :2,
          FILE_DATE = :3,
          LOAD_DATE = SYSDATE
          WHERE INTERPRO.DB_VERSION.DBCODE = 'S'
        """, (version, n_swissprot, date)
    )

    cur.execute(
        """
        UPDATE INTERPRO.DB_VERSION
        SET
          VERSION = :1,
          ENTRY_COUNT = :2,
          FILE_DATE = :3,
          LOAD_DATE = SYSDATE
          WHERE INTERPRO.DB_VERSION.DBCODE = 'T'
        """, (version, n_trembl, date)
    )

    cur.execute(
        """
        UPDATE INTERPRO.DB_VERSION
        SET
          VERSION = :1,
          ENTRY_COUNT = :2,
          FILE_DATE = :3,
          LOAD_DATE = SYSDATE
          WHERE INTERPRO.DB_VERSION.DBCODE = 'u'
        """, (version, n_swissprot + n_trembl, date)
    )

    con.commit()
    cur.close()
    con.close()


def refresh_sequences_to_scan(url: str):
    logger.info("refreshing PROTEIN_TO_SCAN")

    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.execute("TRUNCATE TABLE INTERPRO.PROTEIN_TO_SCAN")
    cur.execute(
        """
        INSERT /*+ APPEND */ INTO INTERPRO.PROTEIN_TO_SCAN (
            PROTEIN_AC, DBCODE, TIMESTAMP, UPI
        )
        SELECT
            IP.PROTEIN_AC, IP.DBCODE, IP.TIMESTAMP, UPI
        FROM INTERPRO.PROTEIN IP
        LEFT OUTER JOIN (
            SELECT DISTINCT
                UPI, AC, CRC64
            FROM UNIPARC.XREF
            INNER JOIN UNIPARC.PROTEIN USING (UPI)
            WHERE DBID IN (2, 3)
        ) UP ON (IP.PROTEIN_AC = UP.AC AND IP.CRC64 = UP.CRC64)
        WHERE IP.PROTEIN_AC IN (
            SELECT NEW_PROTEIN_AC
            FROM INTERPRO.PROTEIN_CHANGES
        )
        """
    )

    logger.info("{} rows inserted".format(cur.rowcount))
    con.commit()
    cur.close()
    con.close()
