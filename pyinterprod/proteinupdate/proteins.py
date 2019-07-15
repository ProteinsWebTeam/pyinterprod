# -*- coding: utf-8 -*-

import os
import sqlite3
from concurrent import futures
from datetime import datetime
from tempfile import mkdtemp, mkstemp
from typing import Generator, Iterable, Optional, Tuple

import cx_Oracle

from . import sprot
from .. import logger, orautils


class ProteinDatabase(object):
    def __init__(self, path: Optional[str] = None, dir: Optional[str] = None):
        if dir is not None:
            os.makedirs(dir, exist_ok=True)
        self.dir = mkdtemp(dir=dir)

        os.environ["SQLITE_TMPDIR"] = self.dir
        if path:
            self.path = path
        else:
            fd, self.path = mkstemp(dir=self.dir)
            os.close(fd)
            os.remove(self.path)

        self._create_table("protein")
        # self._create_table("protein_old")

    def __del__(self):
        self._clean()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._clean()

    def __iter__(self):
        with sqlite3.connect(self.path) as con:
            for row in con.execute("SELECT * FROM protein"):
                yield row

    @property
    def size(self) -> int:
        size = 0
        for f in os.listdir(self.dir):
            try:
                s = os.path.getsize(os.path.join(self.dir, f))
            except FileNotFoundError:
                pass
            else:
                size += s

        return size

    def _clean(self):
        if not self.dir:
            return

        for f in os.listdir(self.dir):
            os.remove(os.path.join(self.dir, f))

        os.rmdir(self.dir)
        self.dir = None

    def _create_table(self, table: str):
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
                """.format(table)
            )

    # def insert_old(self, it: Iterable) -> int:
    #     with sqlite3.connect(self.path) as con:
    #         query = "INSERT INTO protein_old VALUES (?, ?, ?, ?, ?, ?, ?)"
    #         with orautils.TablePopulator(con, query) as table:
    #             for row in it:
    #                 table.insert(row)
    #
    #         con.commit()
    #
    #     return table.rowcount
    #
    # def get_deleted(self) -> Generator[str, None, None]:
    #     with sqlite3.connect(self.path) as con:
    #         cur = con.execute(
    #             """
    #             SELECT accession
    #             FROM protein_old
    #             EXCEPT
    #             SELECT accession
    #             FROM protein
    #             """
    #         )
    #
    #         for row in cur:
    #             yield row[0]
    #
    # def get_new(self) -> Generator[str, None, None]:
    #     with sqlite3.connect(self.path) as con:
    #         cur = con.execute(
    #             """
    #             SELECT
    #               accession,
    #               identifier,
    #               is_reviewed,
    #               crc64,
    #               length,
    #               is_fragment,
    #               taxon_id
    #             FROM protein
    #             WHERE accession NOT IN (
    #               SELECT accession
    #               FROM protein_old
    #             )
    #             """
    #         )
    #
    #         for row in cur:
    #             yield row
    #
    # def get_annotation_changes(self) -> Generator[Tuple, None, None]:
    #     with sqlite3.connect(self.path) as con:
    #         cur = con.execute(
    #             """
    #             SELECT
    #               accession, p1.identifier, p1.is_reviewed, p1.crc64,
    #               p1.length, p1.is_fragment, p1.taxon_id
    #             FROM protein AS p1
    #             INNER JOIN protein_old AS p2
    #               USING (accession)
    #             WHERE p1.crc64 = p2.crc64
    #             AND (
    #                  p1.identifier != p2.identifier
    #               OR p1.is_reviewed != p2.is_reviewed
    #               OR p1.crc64 != p2.crc64
    #               OR p1.length != p2.length
    #               OR p1.is_fragment != p2.is_fragment
    #               OR p1.taxon_id != p2.taxon_id
    #             )
    #             """
    #         )
    #
    #         for row in cur:
    #             yield row
    #
    # def get_sequence_changes(self) -> Generator[Tuple, None, None]:
    #     with sqlite3.connect(self.path) as con:
    #         cur = con.execute(
    #             """
    #             SELECT
    #               accession, p1.identifier, p1.is_reviewed, p1.crc64,
    #               p1.length, p1.is_fragment, p1.taxon_id
    #             FROM protein AS p1
    #             INNER JOIN protein_old AS p2
    #               USING (accession)
    #             WHERE p1.crc64 != p2.crc64
    #             """
    #         )
    #
    #         for row in cur:
    #             yield row

    def drop(self):
        try:
            os.remove(self.path)
        except FileNotFoundError:
            pass


def load(user: str, dsn: str, swissprot_path: str, trembl_path: str,
         dir: Optional[str]=None):
    url = orautils.make_connect_string(user, dsn)

    logger.info("loading proteins")
    with ProteinDatabase(dir=dir) as db:
        count = sprot.load(swissprot_path, db.path, "protein")
        logger.info("Swiss-Prot: {} proteins".format(count))

        count = sprot.load(trembl_path, db.path, "protein")
        logger.info("TrEMBL: {} proteins".format(count))

        logger.info("disk space used: {:.0f} MB".format(db.size/1024**2))

        con = cx_Oracle.connect(url)
        logger.info("{} proteins inserted".format(_insert_all(con, db)))

    # Update annotations
    logger.info("{} annotations updated".format(_update_annotations(con)))

    # Update sequences (CRC64 hashes)
    logger.info("{} sequences updated".format(_update_sequences(con)))

    # Track deleted proteins
    logger.info("{} proteins to delete".format(_track_deleted(con)))

    # Insert new proteins
    logger.info("{} new proteins".format(_insert_new(con)))
    con.commit()

    # Delete staging table
    cur = con.cursor()
    orautils.drop_table(cur, "INTERPRO", "PROTEIN_STG", purge=True)
    cur.close()
    con.close()


def load2(user: str, dsn: str, swissp_src: str, trembl_src: str,
          dir: Optional[str]=None):
    url = orautils.make_connect_string(user, dsn)
    database_old = database_new = None
    with futures.ProcessPoolExecutor(max_workers=2) as executor:
        fs = {
            executor.submit(_export_proteins, url, dir): "old",
            executor.submit(_load_flat_files, swissp_src, trembl_src, dir): "new"
        }

        for f in futures.as_completed(fs):
            try:
                database, swissp_cnt, trembl_cnt = f.result()
            except Exception as exc:
                logger.debug(f"{fs[f]}: exited ({exc})")
            else:
                if fs[f] == "old":
                    database_old = database
                else:
                    database_new = database

                logger.debug(f"{fs[f]} counts: {swissp_cnt} (Swiss-Prot), {trembl_cnt} (TrEMBL)")

    print(database_new)
    print(database_old)


def _export_proteins(url: str, dir: Optional[str]=None) -> Tuple[str, int, int]:
    fd, database = mkstemp(dir=dir)
    os.close(fd)
    os.remove(database)

    swissp_cnt = trembl_cnt = 0
    with sqlite3.connect(database) as con1:
        con1.execute(
            """
            CREATE TABLE protein (
              accession TEXT NOT NULL PRIMARY KEY,
              identifier TEXT NOT NULL,
              is_reviewed INTEGER NOT NULL,
              crc64 TEXT NOT NULL,
              length INTEGER NOT NULL,
              is_fragment INTEGER NOT NULL,
              taxon_id INTEGER NOT NULL
            )
            """
        )

        query = "INSERT INTO protein VALUES (?, ?, ?, ?, ?, ?, ?)"
        with orautils.TablePopulator(con1, query) as table:
            con2 = cx_Oracle.connect(url)
            cur2 = con2.cursor()
            cur2.execute(
                """
                SELECT
                  PROTEIN_AC, NAME, DBCODE, CRC64, LEN, FRAGMENT, TAX_ID
                FROM INTERPRO.PROTEIN
                """
            )

            rows = []
            for row in cur2:
                if row[2] == 'S':
                    swissp_cnt += 1
                    is_reviewed = 1
                else:
                    trembl_cnt += 1
                    is_reviewed = 0

                table.insert((
                    row[0], row[1], is_reviewed, row[3], row[4],
                    1 if row[5] == 'Y' else 0, row[6]
                ))

            cur2.close()
            con2.close()

        con1.commit()

    return database, swissp_cnt, trembl_cnt


def _load_flat_files(swissp_src: str, trembl_src: str, dir: Optional[str]=None) -> Tuple[str, int, int]:
    fd, database = mkstemp(dir=dir)
    os.close(fd)
    os.remove(database)

    with sqlite3.connect(database) as con:
        con.execute(
            """
            CREATE TABLE protein (
              accession TEXT NOT NULL PRIMARY KEY,
              identifier TEXT NOT NULL,
              is_reviewed INTEGER NOT NULL,
              crc64 TEXT NOT NULL,
              length INTEGER NOT NULL,
              is_fragment INTEGER NOT NULL,
              taxon_id INTEGER NOT NULL
            )
            """
        )

    swissp_cnt = sprot.load(swissp_src, database, "protein")
    trembl_cnt = sprot.load(trembl_src, database, "protein")
    return database, swissp_cnt, trembl_cnt

def _insert_all(con: cx_Oracle.Connection, db: ProteinDatabase) -> int:
    cur = con.cursor()
    orautils.drop_table(cur, "INTERPRO", "PROTEIN_STG", purge=True)
    cur.execute(
        """
        CREATE TABLE INTERPRO.PROTEIN_STG (
          PROTEIN_AC VARCHAR2(15) NOT NULL,
          NAME VARCHAR2(16) NOT NULL,
          DBCODE CHAR(1) NOT NULL,
          CRC64 VARCHAR2(16) NOT NULL,
          LEN NUMBER(5) NOT NULL,
          FRAGMENT CHAR(1) NOT NULL,
          TAX_ID NUMBER(15) NOT NULL
        ) NOLOGGING
        """
    )
    cur.close()

    query = """
      INSERT /*+ APPEND */ INTO INTERPRO.PROTEIN_STG
      VALUES(:1, :2, :3, :4, :5, :6, :7)
    """
    table = orautils.TablePopulator(con, query, autocommit=True)

    for ac, name, is_rev, crc64, length, is_frag, taxid in db:
        table.insert((ac, name, 'S' if is_rev else 'T', crc64, length,
                     'Y' if is_frag else 'N', taxid))

    table.close()

    cur = con.cursor()
    cur.execute(
        """
        CREATE UNIQUE INDEX UI_PROTEIN_STG
        ON INTERPRO.PROTEIN_STG (PROTEIN_AC) NOLOGGING
        """
    )
    cur.close()

    return table.rowcount


def _update_annotations(con: cx_Oracle.Connection) -> int:
    cur = con.cursor()
    cur.execute(
        """
        SELECT
          PS.NAME, PS.DBCODE, PS.LEN, PS.FRAGMENT, PS.TAX_ID, PS.PROTEIN_AC
        FROM INTERPRO.PROTEIN_STG PS
        INNER JOIN INTERPRO.PROTEIN P
          ON PS.PROTEIN_AC = P.PROTEIN_AC
        WHERE PS.CRC64 = P.CRC64
        AND (
          PS.NAME != P.NAME
          OR PS.DBCODE != P.DBCODE
          OR PS.LEN != P.LEN
          OR PS.FRAGMENT != P.FRAGMENT
          OR PS.TAX_ID != P.TAX_ID
        )
        """
    )
    query = """
        UPDATE INTERPRO.PROTEIN
        SET
          NAME = :1, DBCODE = :2,  LEN = :3, TIMESTAMP = SYSDATE,
          USERSTAMP = USER, FRAGMENT = :4, TAX_ID = :5
        WHERE PROTEIN_AC = :6
    """
    table = orautils.TablePopulator(con, query)
    for row in cur:
        table.update(row)
    table.close()
    cur.close()
    return table.rowcount


def _update_sequences(con: cx_Oracle.Connection) -> int:
    cur = con.cursor()
    orautils.drop_table(cur, "INTERPRO", "PROTEIN_CHANGES", purge=True)
    cur.execute(
        """
        CREATE TABLE INTERPRO.PROTEIN_CHANGES
        (PROTEIN_AC VARCHAR2(15) PRIMARY KEY NOT NULL)
        """
    )
    cur.execute(
        """
        SELECT
          PS.NAME, PS.DBCODE, PS.CRC64, PS.LEN,
          PS.FRAGMENT, PS.TAX_ID, PS.PROTEIN_AC
        FROM INTERPRO.PROTEIN_STG PS
        INNER JOIN INTERPRO.PROTEIN P
          ON PS.PROTEIN_AC = P.PROTEIN_AC AND PS.CRC64 != P.CRC64
        """
    )

    query = """
        INSERT /*+ APPEND */ INTO INTERPRO.PROTEIN_CHANGES
        VALUES (:1)
    """
    table1 = orautils.TablePopulator(con, query, autocommit=True)
    query = """
        UPDATE INTERPRO.PROTEIN
        SET
          NAME = :1, DBCODE = :2, CRC64 = :3,  LEN = :4,
          TIMESTAMP = SYSDATE, USERSTAMP = USER, FRAGMENT = :5,
          TAX_ID = :6
        WHERE PROTEIN_AC = :7
    """
    table2 = orautils.TablePopulator(con, query)

    for row in cur:
        table1.insert((row[6],))
        table2.update(row)

    table1.close()
    table2.close()
    cur.close()
    return table1.rowcount


def _track_deleted(con) -> int:
    cur = con.cursor()
    orautils.drop_table(cur, "INTERPRO", "PROTEIN_TO_DELETE", purge=True)
    cur.execute(
        """
        CREATE TABLE INTERPRO.PROTEIN_TO_DELETE
        (ID NUMBER NOT NULL, PROTEIN_AC VARCHAR2(15) NOT NULL) NOLOGGING
        """
    )

    query = "INSERT INTO INTERPRO.PROTEIN_TO_DELETE VALUES (:1, :2)"
    table = orautils.TablePopulator(con, query)

    cur.execute(
        """
        SELECT PROTEIN_AC FROM INTERPRO.PROTEIN
        MINUS
        SELECT PROTEIN_AC FROM INTERPRO.PROTEIN_STG
        """
    )
    for i, (protein_ac,) in enumerate(cur):
        table.insert((i+1, protein_ac))

    table.close()
    cur.execute(
        """
        CREATE UNIQUE INDEX UI_PROTEIN_TO_DELETE
        ON INTERPRO.PROTEIN_TO_DELETE (ID) NOLOGGING
        """
    )
    cur.close()
    return table.rowcount


def _insert_new(con: cx_Oracle.Connection) -> int:
    cur = con.cursor()
    orautils.drop_table(cur, "INTERPRO", "PROTEIN_NEW", purge=True)
    cur.execute(
        """
        CREATE TABLE INTERPRO.PROTEIN_NEW NOLOGGING
        AS SELECT *
        FROM INTERPRO.PROTEIN_STG
        WHERE PROTEIN_AC NOT IN (SELECT PROTEIN_AC FROM INTERPRO.PROTEIN)
        """
    )

    cur.execute("SELECT COUNT(*) FROM INTERPRO.PROTEIN_NEW")
    rowcount = cur.fetchone()[0]

    cur.execute(
        """
        INSERT /*+ APPEND */ INTO INTERPRO.PROTEIN_CHANGES
        SELECT PROTEIN_AC FROM INTERPRO.PROTEIN_NEW
        """
    )

    cur.execute(
        """
        INSERT INTO INTERPRO.PROTEIN
        SELECT PROTEIN_AC, NAME, DBCODE, CRC64, LEN, SYSDATE, USER,
               FRAGMENT, 'N', TAX_ID
        FROM INTERPRO.PROTEIN_NEW
        """
    )

    orautils.drop_table(cur, "INTERPRO", "PROTEIN_NEW", purge=True)
    cur.close()
    return rowcount


def update(user: str, dsn: str, version: str, date: str):
    _delete_obsolete(user, dsn, truncate=True)
    _update_database_info(user, dsn, version, date)


def _delete_obsolete(user: str, dsn: str, truncate: bool=False):
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
        partitions = orautils.get_partitions(cur, to, tn)
        if partitions:
            for p in partitions:
                jobs.append((tn, tc, p["name"]))
        else:
            jobs.append((tn, tc, None))

    with futures.ThreadPoolExecutor(max_workers=len(jobs)) as executor:
        fs = {}

        for tn, tc, pn in jobs:
            f = executor.submit(_delete_iter, url, tn, tc, count, 10000, pn)
            fs[f] = (tn, tc, pn)

        num_errors = 0
        for f in futures.as_completed(fs):
            tn, tc, pn = fs[f]
            name = tn if pn is None else "{} ({})".format(tn, pn)

            try:
                f.result()  # returns None
            except Exception as exc:
                logger.debug("{}: exited ({}: {})".format(name, exc.__class__.__name__, exc))
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
        logger.debug("enabling: {}".format(tc))

        if not orautils.toggle_constraint(cur, to, tn, tc, True):
            logger.error("could not enable {}".format(tc))
            num_errors += 1

    cur.close()
    con.close()
    if num_errors:
        raise RuntimeError("{} constraints "
                           "could not be disabled".format(num_errors))

    logger.info("complete")


def _update_database_info(user: str, dsn: str, version: str, date: str):
    date = datetime.strptime(date, "%d-%b-%Y")
    url = orautils.make_connect_string(user, dsn)

    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.execute("SELECT COUNT(*) FROM INTERPRO.PROTEIN WHERE DBCODE = 'S'")
    n_swissprot = cur.fetchone()[0]

    cur.execute("SELECT COUNT(*) FROM INTERPRO.PROTEIN WHERE DBCODE = 'T'")
    n_trembl = cur.fetchone()[0]

    cur.executemany(
        """
        UPDATE INTERPRO.DB_VERSION
        SET
          VERSION = :1,
          ENTRY_COUNT = :2,
          FILE_DATE = :3,
          LOAD_DATE = SYSDATE
          WHERE DBCODE = :4
        """,
        [
            (version, n_swissprot, date, 'S'),
            (version, n_trembl, date, 'T'),
            (version, n_swissprot + n_trembl, date, 'u')
        ]
    )
    con.commit()
    cur.close()
    con.close()


def find_protein_to_refresh(user: str, dsn: str):
    logger.info("PROTEIN_TO_SCAN: refreshing")
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    orautils.drop_table(cur, "INTERPRO", "PROTEIN_TO_SCAN", purge=True)

    # # Assume CRC64 have been checked and that no mismatches were found
    # cur.execute(
    #     """
    #     CREATE TABLE INTERPRO.PROTEIN_TO_SCAN NOLOGGING
    #     AS
    #     SELECT P.NEW_PROTEIN_AC AS PROTEIN_AC, X.UPI
    #     FROM INTERPRO.PROTEIN_CHANGES P
    #     INNER JOIN UNIPARC.XREF X
    #       ON P.NEW_PROTEIN_AC = X.AC
    #     WHERE X.DBID IN (2, 3)
    #     AND X.DELETED = 'N'
    #     """
    # )

    # Assume CRC64 have been checked and that no mismatches were found
    cur.execute(
        """
        CREATE TABLE INTERPRO.PROTEIN_TO_SCAN NOLOGGING
        AS
        SELECT P.PROTEIN_AC, X.UPI
        FROM INTERPRO.PROTEIN_CHANGES P
        INNER JOIN UNIPARC.XREF X
          ON P.PROTEIN_AC = X.AC
        WHERE X.DBID IN (2, 3)
        AND X.DELETED = 'N'
        """
    )

    orautils.gather_stats(cur, "INTERPRO", "PROTEIN_TO_SCAN")

    cur.execute(
        """
        CREATE UNIQUE INDEX PK_PROTEIN_TO_SCAN
        ON INTERPRO.PROTEIN_TO_SCAN (PROTEIN_AC) NOLOGGING
        """
    )

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
        logger.error("{}: {} / {}".format(*row))
        num_errors += 1

    cur.close()
    con.close()

    if num_errors:
        raise RuntimeError("{} CRC64 mismatches".format(num_errors))
    else:
        logger.info("no CRC64 mismatches")


def _count_proteins_to_delete(cur: cx_Oracle.Cursor) -> int:
    cur.execute(
        """
        SELECT COUNT(*)
        FROM INTERPRO.PROTEIN_TO_DELETE
        """
    )
    return cur.fetchone()[0]


def _delete_iter(url: str, table: str, column: str, stop: int, step: int,
                 partition: Optional[str]=None):
    con = cx_Oracle.connect(url)
    cur = con.cursor()

    if partition:
        dst = "{} PARTITION ({})".format(table, partition)
    else:
        dst = "{}".format(table)

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
        con.commit()
        orautils.gather_stats(cur, "INTERPRO", table, partition)

    cur.close()
    con.close()
