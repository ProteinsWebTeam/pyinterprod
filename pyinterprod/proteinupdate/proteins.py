# -*- coding: utf-8 -*-

import os
import sqlite3
from concurrent import futures
from datetime import datetime
from tempfile import mkstemp
from typing import Optional, Tuple

import cx_Oracle

from . import sprot
from .. import logger, orautils


def load(user: str, dsn: str, swissp_src: str, trembl_src: str,
         tmpdir: Optional[str]=None):
    url = orautils.make_connect_string(user, dsn)
    database_old = database_new = None
    with futures.ProcessPoolExecutor(max_workers=2) as executor:
        fs = {
            executor.submit(_export_proteins, url, tmpdir): "old",
            executor.submit(_load_flat_files, swissp_src, trembl_src, tmpdir): "new"
        }

        for f in futures.as_completed(fs):
            try:
                database, swissp_cnt, trembl_cnt = f.result()
            except Exception as exc:
                logger.error(f"{fs[f]}: exited ({exc})")
            else:
                if fs[f] == "old":
                    database_old = database
                else:
                    database_new = database

                logger.info(f"{fs[f]} counts: {swissp_cnt} (Swiss-Prot), {trembl_cnt} (TrEMBL)")

    if not database_old or not database_new:
        if database_old:
            os.remove(database_old)
        if database_new:
            os.remove(database_new)
        raise RuntimeError("failed")

    size = os.path.getsize(database_old) + os.path.getsize(database_new)
    logger.info("disk space used: {:.0f}MB".format(size/1024**2))

    _diff_databases(url, database_old, database_new)
    os.remove(database_old)
    os.remove(database_new)


def _export_proteins(url: str, tmpdir: Optional[str]=None) -> Tuple[str, int, int]:
    fd, database = mkstemp(dir=tmpdir)
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


def _load_flat_files(swissp_src: str, trembl_src: str, tmpdir: Optional[str]=None) -> Tuple[str, int, int]:
    fd, database = mkstemp(dir=tmpdir)
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


def _diff_databases(url: str, database_old: str, database_new: str):
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    orautils.drop_table(cur, "INTERPRO", "PROTEIN_CHANGES", purge=True)
    orautils.drop_table(cur, "INTERPRO", "PROTEIN_TO_DELETE", purge=True)
    cur.execute(
        """
        CREATE TABLE INTERPRO.PROTEIN_CHANGES
        (PROTEIN_AC VARCHAR2(15) PRIMARY KEY NOT NULL)
        """
    )
    cur.execute(
        """
        CREATE TABLE INTERPRO.PROTEIN_TO_DELETE
        (ID NUMBER NOT NULL, PROTEIN_AC VARCHAR2(15) NOT NULL) NOLOGGING
        """
    )
    cur.close()

    # Annotation changes
    query = """
        UPDATE INTERPRO.PROTEIN
        SET
          NAME = :1, DBCODE = :2,  LEN = :3, TIMESTAMP = SYSDATE,
          USERSTAMP = USER, FRAGMENT = :4, TAX_ID = :5
        WHERE PROTEIN_AC = :6
    """
    ann_changes = orautils.TablePopulator(con, query)

    # Sequence changes
    query = """
        UPDATE INTERPRO.PROTEIN
        SET
          NAME = :1, DBCODE = :2, CRC64 = :3,  LEN = :4,
          TIMESTAMP = SYSDATE, USERSTAMP = USER, FRAGMENT = :5,
          TAX_ID = :6
        WHERE PROTEIN_AC = :7
    """
    seq_changes = orautils.TablePopulator(con, query)

    # Obsolete proteins
    query = "INSERT INTO INTERPRO.PROTEIN_TO_DELETE VALUES (:1, :2)"
    del_changes = orautils.TablePopulator(con, query)

    # New proteins
    query = """
        INSERT INTO INTERPRO.PROTEIN
        VALUES (:1, :2, :3, :4, :5, SYSDATE, USER, :6, 'N', :7)
    """
    new_proteins = orautils.TablePopulator(con, query)

    # Proteins to scan
    query = "INSERT INTO INTERPRO.PROTEIN_CHANGES VALUES (:1)"
    all_changes = orautils.TablePopulator(con, query)

    it1 = iter(_iter_proteins(database_old))
    it2 = iter(_iter_proteins(database_new))
    row1 = next(it1)
    row2 = next(it2)
    is_alive1 = is_alive2 = True

    while is_alive1 or is_alive2:
        """
        row:
            - accession (str)
            - name (str)
            - is_reviewed (int: 0/1)
            - crc64 (str)
            - length (int)
            - is_fragment (int: 0/1)
            - taxon_id (int)
        """
        acc1 = row1[0]
        acc2 = row2[0]

        if acc1 == acc2:
            if row1[3] == row2[3]:
                # Same CRC64
                for i in (1, 2, 4, 5, 6):
                    if row1[i] != row2[i]:
                        ann_changes.update((
                            row2[1],
                            'S' if row2[2] else 'T',
                            row2[4],
                            'Y' if row2[5] else 'N',
                            row2[6],
                            acc2
                        ))
                        break
            else:
                # Sequence changed
                seq_changes.update((
                    row2[1],
                    'S' if row2[2] else 'T',
                    row2[3],
                    row2[4],
                    'Y' if row2[5] else 'N',
                    row2[6],
                    acc2
                ))
                all_changes.insert((acc2,))

            try:
                row1 = next(it1)
            except StopIteration:
                is_alive1 = False

            try:
                row2 = next(it2)
            except StopIteration:
                is_alive2 = False
        elif acc1 < acc2 and is_alive1:
            del_changes.insert((del_changes.rowcount, acc1))
            try:
                row1 = next(it1)
            except StopIteration:
                is_alive1 = False
        elif acc1 > acc2 and is_alive2:
            new_proteins.insert((
                acc2,
                row2[1],
                'S' if row2[2] else 'T',
                row2[3],
                row2[4],
                'Y' if row2[5] else 'N',
                row2[6]
            ))
            all_changes.insert((acc2,))
            try:
                row2 = next(it2)
            except StopIteration:
                is_alive2 = False
        elif is_alive1:
            try:
                row1 = next(it1)
            except StopIteration:
                is_alive1 = False
        elif is_alive2:
            try:
                row2 = next(it2)
            except StopIteration:
                is_alive2 = False

    ann_changes.close()
    seq_changes.close()
    del_changes.close()
    new_proteins.close()
    all_changes.close()
    con.commit()
    cur = con.cursor()
    cur.execute(
        """
        CREATE UNIQUE INDEX UI_PROTEIN_TO_DELETE
        ON INTERPRO.PROTEIN_TO_DELETE (ID) NOLOGGING
        """
    )
    cur.close()
    con.close()

    logger.info(f"annotation changes: {ann_changes.rowcount:>10}")
    logger.info(f"sequence changes:   {seq_changes.rowcount:>10}")
    logger.info(f"obsolete proteins:  {del_changes.rowcount:>10}")
    logger.info(f"new proteins:       {new_proteins.rowcount:>10}")


def _iter_proteins(database: str):
    with sqlite3.connect(database) as con:
        for row in con.execute("SELECT * FROM protein ORDER BY accession"):
            yield row


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
                logger.info("{}: exited ({})".format(name, exc))
                num_errors += 1
            else:
                logger.info("{}: done".format(name))

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
