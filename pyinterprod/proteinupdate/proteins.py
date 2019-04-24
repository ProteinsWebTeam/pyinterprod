import os
import sqlite3
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
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
            try:
                os.makedirs(dir, exist_ok=True)
            except TypeError:
                # dir is None
                pass

            fd, self.path = mkstemp(dir=dir)
            os.close(fd)
            os.remove(self.path)
            self.temporary = True

        self._create_table("protein_old")
        self._create_table("protein_new")

    def __del__(self):
        if self.temporary:
            self.drop()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
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

    query = """
      INSERT INTO INTERPRO.PROTEIN_CHANGES (NEW_PROTEIN_AC)
      VALUES (:1)
    """
    table = orautils.TablePopulator(con, query, buffer_size=_MAX_ITEMS)
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
        table.insert((row[0],))

    table.close()
    con.commit()
    cur.close()
    return count


def _insert_proteins(con: cx_Oracle.Connection, db: ProteinDatabase) -> int:
    # TIMESTAMP and USERSTAMP will be set to their DEFAULT
    query1 = """
        INSERT INTO INTERPRO.PROTEIN (
          PROTEIN_AC, NAME, DBCODE, CRC64, LEN, FRAGMENT,
          STRUCT_FLAG, TAX_ID
        )
        VALUES (:1, :2, :3, :4, :5, :6, 'N', :7)
    """
    table1 = orautils.TablePopulator(con, query1, buffer_size=_MAX_ITEMS)

    query2 = """
        INSERT INTO INTERPRO.PROTEIN_CHANGES (NEW_PROTEIN_AC)
        VALUES (:1)
    """
    table2 = orautils.TablePopulator(con, query2, buffer_size=_MAX_ITEMS)

    count = 0
    for row in db.get_new():
        table1.insert((
            row[0],                     # accession
            row[1],                     # identifier
            'S' if row[2] else 'T',     # dbcode
            row[3],                     # crc64
            row[4],                     # length
            'Y' if row[5] else 'N',     # sequence status (fragment)
            row[6]                      # taxon ID
        ))
        table2.insert((row[0],))
        count += 1

    table1.close()
    table2.close()
    return count


def _prepare_deletion(con: cx_Oracle.Connection, db: ProteinDatabase) -> int:
    cur = con.cursor()
    orautils.drop_table(cur, "INTERPRO", "PROTEIN_TO_DELETE")
    cur.execute(
        """
        CREATE TABLE INTERPRO.PROTEIN_TO_DELETE
        (
            ID NUMBER NOT NULL,
            PROTEIN_AC VARCHAR2(15) NOT NULL
        ) NOLOGGING
        """
    )
    cur.close()

    query = """
        INSERT INTO INTERPRO.PROTEIN_TO_DELETE
        VALUES (:1, :2)
    """
    table = orautils.TablePopulator(con, query, buffer_size=_MAX_ITEMS)

    count = 0
    for accession in db.get_deleted():
        table.insert((count, accession))
        count += 1

    cur = con.cursor()
    orautils.gather_stats(cur, "INTERPRO", "PROTEIN_TO_DELETE")
    cur.execute(
        """
        ALTER TABLE INTERPRO.PROTEIN_TO_DELETE
        ADD CONSTRAINT PK_PROTEIN_TO_DELETE
        PRIMARY KEY (ID)
        """
    )

    cur.close()
    return count


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


def insert_new(user: str, dsn: str, swissprot_path: str, trembl_path: str,
               dir: Optional[str]=None):
    logger.info("loading proteins")
    with ProteinDatabase(dir=dir) as db:
        url = orautils.make_connect_string(user, dsn)

        count = db.insert_old(_get_proteins(url))
        logger.info("UniProt proteins in InterPro: {}".format(count))

        count = sprot.load(swissprot_path, db.path, "protein_new")
        logger.info("New Swiss-Prot: {} proteins".format(count))

        count = sprot.load(trembl_path, db.path, "protein_new")
        logger.info("New TrEMBL: {} proteins".format(count))

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


def _delete_iter(url: str, table: str, column: str, stop: int, step: int,
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

        logger.debug("{}: gathering statistics".format(name))
        orautils.gather_stats(cur, "INTERPRO", table, partition)

    cur.close()
    con.close()


def delete_obsolete(user: str, dsn: str, truncate: bool=False):
    url = orautils.make_connect_string(user, dsn)
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
    num_errors = 0
    for t in tables:
        if t["constraint"]:
            to, tn, tc = t["owner"], t["name"], t["constraint"]
            if not orautils.toggle_constraint(cur, to, tn, tc, False):
                logger.error("could not disable {}.{}".format(tn, tc))
                num_errors += 1

    if num_errors:
        cur.close()
        con.close()
        raise RuntimeError("{} constraints "
                           "could not be disabled".format(num_errors))

    count = _count_proteins_to_delete(cur)
    logger.info("{} proteins to delete".format(count))

    jobs = []
    for t in tables:
        to, tn, tc = t["owner"], t["name"], t["column"]

        if tn == "MATCH":
            for p in orautils.get_partitions(cur, to, tn):
                jobs.append((tn, tc, p["name"]))
        else:
            jobs.append((tn, tc, None))

        # # TODO: fix: FEATURE_MATCH has partitions but then
        # # ORA-02149: Specified partition does not exist is raised
        # partitions = orautils.get_partitions(cur, to, tn)
        # if partitions:
        #     for p in partitions:
        #         jobs.append((tn, tc, p["name"]))
        # else:
        #     jobs.append((tn, tc, None))

    with ThreadPoolExecutor(max_workers=len(jobs)) as executor:
        fs = {}

        for tn, tc, pn in jobs:
            f = executor.submit(_delete_iter, url, tn, tc, count, _MAX_ITEMS, pn)
            fs[f] = (tn, tc, pn)

        num_errors = 0
        for f in as_completed(fs):
            tn, tc, pn = fs[f]
            name = tn if pn is None else "{} ({})".format(tn, pn)

            try:
                f.result()  # returns None
            except Exception as exc:
                logger.error("{}: exited ({})".format(name, exc))
                num_errors += 1
            else:
                logger.debug("{}: done".format(name))

        if num_errors:
            cur.close()
            con.close()
            raise RuntimeError("{} tables failed".format(num_errors))

    logger.info("enabling referential constraints")
    num_errors = 0
    constraints = set()
    for t in tables:
        to, tn, tc = t["owner"], t["name"], t["constraint"]

        if not tc or tc in constraints:
            """
            Either no constraint
            or prevent the same constrain to be enabled several times
            """
            continue

        constraints.add(tc)
        logger.debug("enabling: {}.{}.{}".format(to, tn, tc))

        if not orautils.toggle_constraint(cur, to, tn, tc, True):
            logger.error("could not enable {}.{}".format(tn, tc))
            num_errors += 1

    cur.close()
    con.close()
    if num_errors:
        raise RuntimeError("{} constraints "
                           "could not be disabled".format(num_errors))

    logger.info("complete")


def update_database_info(user: str, dsn: str, version: str, date: str):
    date = datetime.strptime(date, "%d-%b-%Y")
    url = orautils.make_connect_string(user, dsn)

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


def find_protein_to_refresh(user: str, dsn: str):
    logger.info("PROTEIN_TO_SCAN: refreshing")
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    cur.execute("TRUNCATE TABLE INTERPRO.PROTEIN_TO_SCAN")

    # Consider CRC64 have been checked and that no mismatches were found
    cur.execute(
        """
        INSERT INTO INTERPRO.PROTEIN_TO_SCAN
        SELECT PC.NEW_PROTEIN_AC, IP.DBCODE, IP.TIMESTAMP, UX.UPI
        FROM INTERPRO.PROTEIN_CHANGES PC
        INNER JOIN INTERPRO.PROTEIN IP
          ON PC.NEW_PROTEIN_AC = IP.PROTEIN_AC
        INNER JOIN UNIPARC.XREF UX
          ON PC.NEW_PROTEIN_AC = UX.AC
        WHERE UX.DBID IN (2, 3)
        AND UX.DELETED = 'N'
        """
    )
    logger.info("PROTEIN_TO_SCAN: {} rows inserted".format(cur.rowcount))
    con.commit()

    orautils.gather_stats(cur, "INTERPRO", "PROTEIN_TO_SCAN")
    orautils.rebuild_indices(cur, "INTERPRO", "PROTEIN_TO_SCAN")

    cur.close()
    con.close()


def check_crc64(user: str, dsn: str):
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    cur.execute(
        """
        SELECT IP.PROTEIN_AC, IP.CRC64, UP.CRC64
        FROM INTERPRO.PROTEIN IP
        INNER JOIN UNIPARC.XREF UX
          ON IP.PROTEIN_AC = UX.AC
        INNER JOIN UNIPARC.PROTEIN UP
          ON UX.UPI = UP.UPI
        WHERE UX.DBID IN (2, 3)
        AND UX.DELETED = 'N'
        AND IP.CRC64 != UP.CRC64
        """
    )

    num_errors = 0
    for row in cur:
        logger.debug("{}: {} / {}".format(*row))
        num_errors += 1

    cur.close()
    con.close()

    if num_errors:
        raise RuntimeError("{} CRC64 mismatches".format(num_errors))
    else:
        logger.info("no CRC64 mismatches")
