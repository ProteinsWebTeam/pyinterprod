# -*- coding: utf-8 -*-

import os
import pickle
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Optional, Sequence, Tuple

import cx_Oracle

from pyinterprod import logger
from pyinterprod.pronto.signature import get_swissprot_descriptions
from pyinterprod.utils import Table, oracle as ora
from . import contrib
from .database import Database
from .match import MATCH_PARTITIONS, SITE_PARTITIONS


FILE_DB_SIG_DIFF = "signatures.anno.changes.pickle"
FILE_SIG_DESCR = "signatures.descr.pickle"


def export_swissprot_descriptions(pg_url, data_dir: str):
    with open(os.path.join(data_dir, FILE_SIG_DESCR), "wb") as fh:
        pickle.dump(get_swissprot_descriptions(pg_url), fh)


def add_staging(url: str, update: Sequence[Tuple[Database, str]]):
    con = cx_Oracle.connect(url)
    cur = con.cursor()

    ora.truncate_table(cur, "METHOD_STG")

    sql = """
        INSERT INTO INTERPRO.METHOD_STG 
        VALUES (:1, :2, :3, :4, :5, :6) 
    """
    with Table(con, sql) as table:
        errors = 0
        for db, src in update:
            if db.identifier == 'Q':
                # HAMAP
                signatures = contrib.hamap.parse_signatures(src)
            elif db.identifier == 'V':
                # PANTHER
                signatures = contrib.panther.parse_signatures(src)
            else:
                logger.error(f"{db.name}: unsupported member database")
                errors += 1
                continue

            for m in signatures:
                table.insert((
                    m.accession,
                    db.identifier,
                    m.name,
                    m.description,
                    m.sig_type,
                    m.abstract
                ))

    if errors:
        cur.close()
        con.close()
        raise RuntimeError(f"{errors} errors occurred")

    con.commit()

    code2name = {db.identifier: db.name for db, _ in update}
    cur.execute(
        """
        SELECT DBCODE, COUNT(*)
        FROM INTERPRO.METHOD_STG
        GROUP BY DBCODE
        """
    )
    for dbcode, cnt in cur:
        logger.info(f"{code2name[dbcode]:<30} {cnt:>10,} signatures")

    cur.close()
    con.close()


def track_signature_changes(url: str, databases: Sequence[Database], data_dir: str):
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    results = {}
    for db in databases:
        # cur.execute(
        #     """
        #     SELECT COUNT(*)
        #     FROM INTERPRO.MATCH
        #     WHERE DBCODE = :1
        #     """, (dbcode,)
        # )
        # num_matches, = cur.fetchone()
        #
        # cur.execute(
        #     """
        #     SELECT COUNT(*)
        #     FROM INTERPRO.MATCH
        #     WHERE DBCODE = :dbcode
        #     AND METHOD_AC IN (
        #         SELECT METHOD_AC
        #         FROM INTERPRO.METHOD
        #         WHERE DBCODE = :dbcode
        #         MINUS
        #         SELECT METHOD_AC
        #         FROM INTERPRO.METHOD_STG
        #         WHERE DBCODE = :dbcode
        #     )
        #     """, dbcode=dbcode
        # )
        # num_deleted_matches, = cur.fetchone()

        cur.execute(
            """
            SELECT M.METHOD_AC, M.NAME, M.DESCRIPTION, M.SIG_TYPE, EM.ENTRY_AC
            FROM INTERPRO.METHOD M
            LEFT OUTER JOIN INTERPRO.ENTRY2METHOD EM 
            ON M.METHOD_AC = EM.METHOD_AC
            WHERE DBCODE = :1
            """, (db.identifier,)
        )
        old_signatures = {row[0]: row[1:] for row in cur}

        cur.execute(
            """
            SELECT METHOD_AC, NAME, DESCRIPTION, SIG_TYPE
            FROM INTERPRO.METHOD_STG
            WHERE DBCODE = :1
            """, (db.identifier,)
        )
        new_signatures = {row[0]: row[1:] for row in cur}

        deleted = []
        name_changes = []
        descr_changes = []
        type_changes = []
        for acc in sorted(old_signatures):
            old_name, old_descr, old_type, entry_acc = old_signatures[acc]

            try:
                new_name, new_descr, new_type = new_signatures.pop(acc)
            except KeyError:
                deleted.append((acc, old_name, old_descr, entry_acc))
                continue

            if old_name != new_name:
                name_changes.append((acc, old_name, new_name))

            if old_descr != new_descr:
                descr_changes.append((acc, old_descr, new_descr))

            if old_type != new_type:
                type_changes.append((acc, old_type, new_type))

        results[db.identifier] = {
            "new": [
                (acc, *new_signatures[acc])
                for acc in sorted(new_signatures)
            ],
            "deleted": deleted,
            "changes": {
                "names": name_changes,
                "descriptions": descr_changes,
                "types": type_changes
            }
        }

    cur.close()
    con.close()

    with open(os.path.join(data_dir, FILE_DB_SIG_DIFF), "wb") as fh:
        pickle.dump(results, fh)


def delete_from_table(url: str, table: str, partition: Optional[str],
                      column: str, step: int, stop: int) -> int:
    con = cx_Oracle.connect(url)
    cur = con.cursor()

    if partition:
        db_obj = f"{table} PARTITION ({partition})"
    else:
        db_obj = table

    cur.execute(
        f"""
        SELECT COUNT(*)
        FROM INTERPRO.{db_obj}
        WHERE {column} IN (
            SELECT METHOD_AC 
            FROM INTERPRO.METHOD_TO_DELETE
        ) 
        """
    )
    num_rows, = cur.fetchone()

    if not num_rows:
        cur.close()
        con.close()
        return num_rows

    for i in range(1, stop, step):
        cur.execute(
            f"""
            DELETE FROM INTERPRO.{db_obj}
            WHERE {column} IN (
              SELECT METHOD_AC
              FROM INTERPRO.METHOD_TO_DELETE
              WHERE ID BETWEEN :1 and :2
            )
            """, (i, i + step - 1)
        )

    con.commit()
    ora.gather_stats(cur, "INTERPRO", table, partition)
    cur.close()
    con.close()
    return num_rows


def delete_obsoletes(url: str, databases: Sequence[Database], **kwargs):
    step = kwargs.get("step", 10000)
    threads = kwargs.get("threads", 8)

    con = cx_Oracle.connect(url)
    cur = con.cursor()

    # track signatures that need to be deleted
    ora.drop_table(cur, "METHOD_TO_DELETE")
    cur.execute(
        """
        CREATE TABLE INTERPRO.METHOD_TO_DELETE (
            ID NUMBER NOT NULL,
            METHOD_AC VARCHAR2(25) NOT NULL
        )
        """
    )

    for db in databases:
        cur.execute(
            """
            INSERT INTO INTERPRO.METHOD_TO_DELETE (ID, METHOD_AC)
            SELECT ROWNUM, METHOD_AC
            FROM (
                SELECT METHOD_AC
                FROM INTERPRO.METHOD
                WHERE DBCODE = :dbcode
                MINUS
                SELECT METHOD_AC
                FROM INTERPRO.METHOD_STG
                WHERE DBCODE = :dbcode
            )
            """, dbcode=db.identifier
        )

    con.commit()
    cur.execute(
        """
        CREATE UNIQUE INDEX UI_METHOD_TO_DELETE
        ON INTERPRO.METHOD_TO_DELETE (ID)
        """
    )

    cur.execute("SELECT COUNT(*) FROM INTERPRO.METHOD_TO_DELETE")
    stop, = cur.fetchone()

    logger.info(f"{stop:,} signatures to delete")

    if not stop:
        # Nothing to delete
        cur.close()
        con.close()
        return

    # Get tables with a FOREIGN KEY to INTERPRO.METHOD
    tables = []
    child_tables = ora.get_child_tables(cur, "INTERPRO", "METHOD")
    for table, constraint, column in child_tables:
        tables.append((table, constraint, column))

    # Add INTERPRO.METHOD as we want also to delete rows in this table
    tables.append(("METHOD", None, "METHOD_AC"))

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

    tasks = []
    for table, constraint, column in tables:
        if table == "MATCH":
            for db in databases:
                partition = MATCH_PARTITIONS[db.identifier]
                tasks.append((table, partition, column))
        elif table == "SITE_MATCH":
            for db in databases:
                partition = SITE_PARTITIONS[db.identifier]
                tasks.append((table, partition, column))
        else:
            tasks.append((table, None, column))

    cur.close()
    con.close()

    with ThreadPoolExecutor(max_workers=threads) as executor:
        fs = {}

        for table, partition, column in tasks:
            args = (url, table, partition, column, step, stop)
            f = executor.submit(delete_from_table, *args)
            fs[f] = (table, partition)

        num_errors = 0
        for f in as_completed(fs):
            table, partition = fs[f]
            if partition:
                name = f"{table} ({partition})"
            else:
                name = table

            try:
                num_rows = f.result()
            except Exception as exc:
                logger.info(f"{name}: failed ({exc})")
                num_errors += 1
            else:
                logger.info(f"{name}: {num_rows:,} rows deleted")

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

        logger.debug(f"enabling: {table}.{constraint}")
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


def update_signatures(url: str):
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.execute(
        """
        MERGE INTO INTERPRO.METHOD M
        USING INTERPRO.METHOD_STG S
          ON (M.METHOD_AC = S.METHOD_AC)
        WHEN MATCHED THEN 
          UPDATE SET M.NAME = S.NAME,
                     M.DESCRIPTION = S.DESCRIPTION,
                     M.SIG_TYPE = S.SIG_TYPE,
                     M.ABSTRACT = S.ABSTRACT,
                     M.TIMESTAMP = SYSDATE
        WHEN NOT MATCHED THEN
          INSERT (METHOD_AC, NAME, DBCODE, CANDIDATE, DESCRIPTION, SIG_TYPE, ABSTRACT)
          VALUES (S.METHOD_AC, S.NAME, S.DBCODE, 'Y', S.DESCRIPTION, S.SIG_TYPE, S.ABSTRACT)
        """
    )
    con.commit()
    cur.close()
    con.close()
