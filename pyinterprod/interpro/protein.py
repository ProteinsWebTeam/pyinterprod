# -*- coding: utf-8 -*-

import os
import sqlite3
from concurrent import futures
from tempfile import mkstemp
from typing import Optional, Tuple

import cx_Oracle

from pyinterprod import logger
from pyinterprod.utils import Table, oracle as ora
from pyinterprod.uniprot import flatfile


def export_proteins(url: str, database: str) -> Tuple[int, int]:
    con1 = sqlite3.connect(database)
    with Table(con1, "INSERT INTO protein VALUES (?, ?, ?, ?, ?, ?, ?)") as t:
        swissp_cnt = trembl_cnt = 0

        con2 = cx_Oracle.connect(url)
        cur2 = con2.cursor()
        cur2.execute(
            """
            SELECT
              PROTEIN_AC, NAME, DBCODE, CRC64, LEN, FRAGMENT, TAX_ID
            FROM INTERPRO.PROTEIN
            """
        )

        for row in cur2:
            if row[2] == 'S':
                swissp_cnt += 1
                is_reviewed = 1
            else:
                trembl_cnt += 1
                is_reviewed = 0

            t.insert((
                row[0],
                row[1],
                is_reviewed,
                row[3],
                row[4],
                1 if row[5] == 'Y' else 0,
                row[6]
            ))

        cur2.close()
        con2.close()

    con1.commit()
    con1.close()

    return swissp_cnt, trembl_cnt


class Entry(object):
    def __init__(self, database):
        self.gen = self.iter(database)
        self.record = None
        self.ok = True
        self.count = 0
        self.next()

    @staticmethod
    def iter(database):
        con = sqlite3.connect(database)

        try:
            # accession, identifier, is_reviewed, crc64, length, is_fragment, taxon_id
            for row in con.execute("SELECT * FROM protein ORDER BY accession"):
                yield row
        finally:
            con.close()

    def next(self):
        if not self.ok:
            return

        try:
            record = next(self.gen)
        except StopIteration:
            self.ok = False
        else:
            self.record = record
            self.count += 1

    def __eq__(self, other):
        return self.record[0] == other.record[0]

    def __lt__(self, other):
        return self.record[0] < other.record[0]

    def diff_sequence(self, other):
        return self.record[3] != other.record[3]

    def diff_annotation(self, other):
        for i in (1, 2, 4, 5, 6):
            if self.record[i] != other.record[i]:
                return True
        return False

    @property
    def accession(self):
        return self.record[0]

    @property
    def annotation(self):
        return (
            self.record[1],
            'S' if self.record[2] else 'T',
            self.record[4],
            'Y' if self.record[5] else 'N',
            self.record[6],
            self.record[0]
        )

    @property
    def sequence(self):
        return (
            self.record[1],
            'S' if self.record[2] else 'T',
            self.record[3],
            self.record[4],
            'Y' if self.record[5] else 'N',
            self.record[6],
            self.record[0]
        )

    @property
    def new(self):
        return (
            self.record[0],
            self.record[1],
            'S' if self.record[2] else 'T',
            self.record[3],
            self.record[4],
            'Y' if self.record[5] else 'N',
            self.record[6]
        )


def init_database(tmpdir: Optional[str] = None) -> str:
    fd, database = mkstemp(suffix=".sqlite", dir=tmpdir)
    os.close(fd)
    os.remove(database)

    con = sqlite3.connect(database)
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
    con.close()

    return database


def track_changes(url: str, swissp: str, trembl: str,
                  tmpdir: Optional[str] = None):
    logger.info("starting")
    database_old = init_database(dir=tmpdir)
    database_new = init_database(dir=tmpdir)

    failed = False
    with futures.ProcessPoolExecutor(max_workers=2) as executor:
        fs = {}
        f = executor.submit(export_proteins, url, database_old)
        fs[f] = "current"

        f = executor.submit(flatfile.load, swissp, trembl, database_new)
        fs[f] = "next release"

        for f in futures.as_completed(fs):
            try:
                swissp_cnt, trembl_cnt = f.result()
            except Exception as exc:
                logger.error(f"{fs[f]}: failed ({exc})")
                failed = True
            else:
                logger.info(f"{fs[f]}: Swiss-Prot: {swissp_cnt}, "
                            f"TrEMBL: {trembl_cnt}")

    if failed:
        os.remove(database_old)
        os.remove(database_new)
        raise RuntimeError("failed to track changes between UniProt releases")

    size = os.path.getsize(database_old) + os.path.getsize(database_new)
    logger.info(f"disk usage: {size/1024**2:.0f} MB")

    con = cx_Oracle.connect(url)
    cur = con.cursor()
    ora.truncate_table(cur, "INTERPRO.PROTEIN_CHANGES")
    ora.truncate_table(cur, "INTERPRO.PROTEIN_TO_DELETE")
    cur.close()

    # Annotation changes
    sql = """
        UPDATE INTERPRO.PROTEIN
        SET NAME = :1, DBCODE = :2, LEN = :3, TIMESTAMP = SYSDATE, 
            USERSTAMP = USER, FRAGMENT = :4, TAX_ID = :5
        WHERE PROTEIN_AC = :6
    """
    ann_table = Table(con, sql)

    # Sequence changes
    sql = """
        UPDATE INTERPRO.PROTEIN
        SET NAME = :1, DBCODE = :2, CRC64 = :3, LEN = :4, TIMESTAMP = SYSDATE, 
            USERSTAMP = USER, FRAGMENT = :5, TAX_ID = :6
        WHERE PROTEIN_AC = :7
    """
    seq_table = Table(con, sql)

    # Obsolete proteins
    sql = "INSERT INTO INTERPRO.PROTEIN_TO_DELETE VALUES (:1, :2)"
    del_table = Table(con, sql)

    # New proteins
    sql = """
        INSERT INTO INTERPRO.PROTEIN
        VALUES (:1, :2, :3, :4, :5, SYSDATE, USER, :6, 'N', :7)
    """
    new_table = Table(con, sql)

    # All changes
    all_table = Table(con, "INSERT INTO INTERPRO.PROTEIN_CHANGES VALUES (:1)")

    current = Entry(database_old)
    forthcoming = Entry(database_new)

    while True:
        if current.ok:
            if forthcoming.ok:
                # Still entries in both databases: compare
                if current == forthcoming:
                    if current.diff_sequence(forthcoming):
                        seq_table.update(forthcoming.sequence)
                        all_table.insert((forthcoming.accession,))
                    elif current.diff_annotation(forthcoming):
                        ann_table.update(forthcoming.annotation)

                    current.next()
                    forthcoming.next()
                elif current < forthcoming:
                    # Entry in current but not in forthcoming: deleted
                    del_table.insert(
                        (del_table.count+1, current.accession))
                    current.next()
                else:
                    # Entry in forthcoming but not in current: new
                    new_table.insert(forthcoming.new)
                    all_table.insert((forthcoming.accession,))
                    forthcoming.next()
            else:
                # Still entries in current, but not in forthcoming: all deleted
                del_table.insert((del_table.count+1, current.accession))
                current.next()
        elif forthcoming.ok:
            # Still entries in forthcoming, but not in current: all new
            new_table.insert(forthcoming.new)
            all_table.insert((forthcoming.accession,))
            forthcoming.next()
        else:
            break

    ann_table.close()
    seq_table.close()
    del_table.close()
    new_table.close()
    all_table.close()
    con.commit()

    os.remove(database_old)
    os.remove(database_new)

    logger.info(f"annotations updated: {ann_table.count:>10}")
    logger.info(f"sequences updated:   {seq_table.count:>10}")
    logger.info(f"obsolete entries:    {del_table.count:>10}")
    logger.info(f"new entries:         {new_table.count:>10}")

    cur = con.cursor()
    ora.gather_stats(cur, "INTERPRO", "PROTEIN_CHANGES")
    ora.gather_stats(cur, "INTERPRO", "PROTEIN_TO_DELETE")
    cur.close()
    con.close()


def delete_obsoletes(url: str, truncate: bool = False, threads: int = 8,
                     step: int = 10000):
    con = cx_Oracle.connect(url)
    cur = con.cursor()

    # Count the protein to delete
    cur.execute("SELECT COUNT(*) FROM INTERPRO.PROTEIN_TO_DELETE")
    stop, = cur.fetchone()

    logger.info(f"{stop:,} proteins to delete")

    if not stop:
        # Nothing to do: exit
        cur.close()
        con.close()
        return

    # Get tables with a FOREIGN KEY to INTERPRO.PROTEIN.PROTEIN_AC
    tables = []
    for table, constraint, column in ora.get_child_tables(cur, "INTERPRO", "PROTEIN"):
        if truncate and table.startswith("MV_") or table.endswith("_NEW"):
            logger.info(f"truncating {table}")
            ora.truncate_table(cur, table)
        else:
            tables.append((table, constraint, column))

    # Add  INTERPRO.PROTEIN as we want also to delete rows in this table
    tables.append(("PROTEIN", None, "PROTEIN_AC"))

    logger.info("disabling referential constraints")
    num_errors = 0
    for table, constraint, column in tables:
        if not constraint:
            continue

        try:
            ora.toggle_constraint(cur, table, constraint, False)
        except cx_Oracle.DatabaseError as exc:
            logger.error(exc)
            num_errors += 1

    if num_errors:
        cur.close()
        con.close()
        raise RuntimeError(f"{num_errors} constraints could not be disabled")

    # Find partitions to run a DELETE statement for each partition
    tasks = []
    for table, constraint, column in tables:
        partitions = ora.get_partitions(cur, "INTERPRO", table)
        if partitions:
            for p in partitions:
                tasks.append((table, p["name"], column))
        else:
            tasks.append((table, None, column))
    cur.close()
    con.close()

    logger.info(f"{len(tasks)} tables to update")
    with futures.ThreadPoolExecutor(max_workers=threads) as executor:
        fs = {}

        for table, partition, column in tasks:
            f = executor.submit(iterative_delete, url, table, partition,
                                column, step, stop)
            fs[f] = (table, partition)

        num_errors = 0
        for f in futures.as_completed(fs):
            table, partition = fs[f]
            if partition:
                name = f"{table} ({partition})"
            else:
                name = table

            try:
                f.result()
            except Exception as exc:
                logger.info(f"{name}: failed ({exc})")
                num_errors += 1
            else:
                logger.info(f"{name}: done")

        if num_errors:
            raise RuntimeError(f"{num_errors} tables failed")

    logger.info("enabling referential constraints")
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    num_errors = 0
    constraints = set()
    for table, constraint, column in tables:
        if not constraint or constraint in constraints:
            """
            Either no constraint
            or prevent the same constrain to be enabled several times
            """
            continue

        constraints.add(constraint)

        try:
            ora.toggle_constraint(cur, table, constraint, True)
        except cx_Oracle.DatabaseError as exc:
            logger.error(exc)
            num_errors += 1

    if num_errors:
        cur.close()
        con.close()
        raise RuntimeError(f"{num_errors} constraints could not be enabled")

    for table, constraint, column in tables:
        for index in ora.get_indexes(cur, "INTERPRO", table):
            if index["unusable"]:
                logger.info(f"rebuilding index {index['name']}")
                ora.rebuild_index(cur, index["name"])

    cur.close()
    con.close()
    logger.info("complete")


def iterative_delete(url: str, table: str, partition: Optional[str],
                     column: str, step: int, stop: int):
    con = cx_Oracle.connect(url)
    cur = con.cursor()

    if partition:
        _table = f"{table} PARTITION ({partition})"
    else:
        _table = table

    cur.execute(
        f"""
        SELECT COUNT(*)
        FROM {_table}
        WHERE {column} IN (
            SELECT PROTEIN_AC 
            FROM INTERPRO.PROTEIN_TO_DELETE
        ) 
        """
    )
    num_rows, = cur.fetchone()

    if not num_rows:
        cur.close()
        con.close()
        return

    for i in range(1, stop, step):
        cur.execute(
            f"""
            DELETE FROM INTERPRO.{_table}
            WHERE {column} IN (
              SELECT PROTEIN_AC
              FROM INTERPRO.PROTEIN_TO_DELETE
              WHERE ID BETWEEN :1 and :2
            )
            """, (i, i + step - 1)
        )

    con.commit()
    ora.gather_stats(cur, "INTERPRO", table, partition)
    cur.close()
    con.close()


def check_proteins(cur: cx_Oracle.Cursor) -> int:
    num_errors = 0
    cur.execute(
        """
        SELECT IP.PROTEIN_AC
        FROM INTERPRO.PROTEIN IP
        LEFT OUTER JOIN UNIPARC.XREF UX
          ON IP.PROTEIN_AC = UX.AC
          AND UX.DBID IN (2, 3)   -- Swiss-Prot/TrEMBL
          AND UX.DELETED = 'N'    -- Not deleted  
        WHERE UX.AC IS NULL 
        """
    )

    for accession, in cur:
        logger.error(f"missing: {accession}")
        num_errors += 1

    cur.execute(
        """
        SELECT IP.PROTEIN_AC, IP.CRC64, UP.CRC64
        FROM INTERPRO.PROTEIN IP
        INNER JOIN UNIPARC.XREF UX
          ON IP.PROTEIN_AC = UX.AC
          AND UX.DBID IN (2, 3)   -- Swiss-Prot/TrEMBL
          AND UX.DELETED = 'N'    -- Not deleted
        INNER JOIN UNIPARC.PROTEIN UP
          ON UX.UPI = UP.UPI
        WHERE IP.CRC64 != UP.CRC64     
        """
    )

    for accession, crc64_1, crc64_2 in cur:
        logger.error(f"CRC64 mismatch: {accession}: {crc64_1} / {crc64_2}")
        num_errors += 1

    return num_errors


def check_proteins_to_scan(url: str):
    con = cx_Oracle.connect(url)
    cur = con.cursor()

    logger.info("checking for errors")
    num_errors = check_proteins(cur)
    if num_errors:
        cur.close()
        con.close()
        raise RuntimeError(f"{num_errors} errors found")

    logger.info("no error found")

    ora.truncate_table(cur, "INTERPRO.PROTEIN_TO_SCAN")
    cur.execute(
        """
        INSERT INTO INTERPRO.PROTEIN_TO_SCAN (PROTEIN_AC, UPI) 
        SELECT P.PROTEIN_AC, X.UPI
        FROM INTERPRO.PROTEIN_CHANGES P
        INNER JOIN UNIPARC.XREF X 
        ON P.PROTEIN_AC = X.AC
        AND X.DBID IN (2, 3)   -- Swiss-Prot/TrEMBL
        AND X.DELETED = 'N'    -- Not deleted  
        """
    )

    con.commit()

    ora.gather_stats(cur, "INTERPRO", "PROTEIN_TO_SCAN")

    cur.close()
    con.close()

    logger.info("complete")
