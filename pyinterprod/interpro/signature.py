# -*- coding: utf-8 -*-

import os
import pickle
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Dict, Optional, Sequence, Set

import cx_Oracle
import psycopg2

from pyinterprod import logger
from pyinterprod.utils import Table, oracle as ora, pg
from . import contrib


SIGNATURES_DESCR_FILE = "swissprot_descr.dat"


def get_swissprot_descriptions(url: str) -> Dict[str, Set[str]]:
    con = psycopg2.connect(**pg.url2dict(url))
    cur = con.cursor("spdecur")
    cur.execute(
        """
        SELECT DISTINCT s2p.signature_acc, pn.text
        FROM interpro.signature2protein s2p
        INNER JOIN interpro.protein_name pn ON s2p.name_id = pn.name_id
        WHERE s2p.is_reviewed            
        """
    )
    signatures = {}
    for acc, text in cur:
        try:
            signatures[acc].add(text)
        except KeyError:
            signatures[acc] = {text}

    cur.close()
    con.close()
    return signatures


def export_swissprot_description(url: str, data_dir: str):
    signatures = get_swissprot_descriptions(url)

    with open(os.path.join(data_dir, SIGNATURES_DESCR_FILE), "wb") as fh:
        pickle.dump(signatures, fh)


def get_key2db(cur: cx_Oracle.Cursor, databases: Sequence[str]) -> dict:
    cur.execute(
        """
        SELECT LOWER(DBSHORT), DBCODE, DBNAME
        FROM INTERPRO.CV_DATABASE
        """
    )
    all_key2db = {row[0]: row[1:] for row in cur}

    key2db = {}
    errors = []
    for dbkey in databases:
        try:
            dbcode, dbname = all_key2db[dbkey.lower()]
        except KeyError:
            errors.append(dbkey)
        else:
            key2db[dbkey.lower()] = (dbcode, dbname)

    if errors:
        cur.close()
        cur.connection.close()
        raise ValueError(f"Unknown database(s): {', '.join(errors)}")

    return key2db


def add_staging(url: str, databases: dict):
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    key2db = get_key2db(cur, databases)
    ora.truncate_table(cur, "METHOD_STG")

    sql = """
        INSERT INTO INTERPRO.METHOD_STG 
        VALUES (:1, :2, :3, :4, :5, :6) 
    """
    with Table(con, sql) as table:
        errors = 0
        for dbkey, src in databases.items():
            dbkey = dbkey.lower()
            dbcode, dbname = key2db[dbkey]
            if dbkey == "panther":
                signatures = contrib.panther.parse_signatures(src)
            else:
                logger.error(f"{dbkey}: unsupported member database")
                errors += 1
                continue

            for m in signatures:
                table.insert((
                    m.accession,
                    dbcode,
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

    code2name = dict(key2db.values())
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


def make_pre_report(url: str, databases: Sequence[str], output: str):
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    key2db = get_key2db(cur, databases)
    results = {}
    for dbcode, dbname in key2db.values():
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
            """, (dbcode,)
        )
        old_signatures = {row[0]: row[1:] for row in cur}

        cur.execute(
            """
            SELECT METHOD_AC, NAME, DESCRIPTION, SIG_TYPE
            FROM INTERPRO.METHOD_STG
            WHERE DBCODE = :1
            """, (dbcode,)
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

        results[dbcode] = {
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

    with open(output, "wb") as fh:
        pickle.dump(results, fh)


def delete_from_table(url: str, table: str, partition: Optional[str],
                      column: str):
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
                SELECT METHOD_AC 
                FROM INTERPRO.METHOD_TO_DELETE
            ) 
            """
    )
    num_rows, = cur.fetchone()
    cur.close()
    con.close()
    return num_rows

    # if not num_rows:
    #     cur.close()
    #     con.close()
    #     return


def delete_obsolete(url: str, databases: Sequence[str], **kwargs):
    threads = kwargs.get("threads", 8)

    con = cx_Oracle.connect(url)
    cur = con.cursor()
    key2db = get_key2db(cur, databases)

    # track signatures that need to be deleted
    ora.drop_table(cur, "METHOD_TO_DELETE")
    cur.execute(
        """
        CREATE TABLE INTERPRO.METHOD_TO_DELETE
        AS
        SELECT METHOD_AC 
        FROM INTERPRO.METHOD
        WHERE 1 = 0
        """
    )

    for dbcode, dbname in key2db.values():
        cur.execute(
            """
            INSERT INTO INTERPRO.METHOD_TO_DELETE (METHOD_AC)
            SELECT METHOD_AC
            FROM INTERPRO.METHOD
            WHERE DBCODE = :dbcode
            MINUS
            SELECT METHOD_AC
            FROM INTERPRO.METHOD_STG
            WHERE DBCODE = :dbcode
            """, dbcode=dbcode
        )

    con.commit()
    cur.execute("SELECT COUNT(*) FROM INTERPRO.METHOD_TO_DELETE")
    cnt, = cur.fetchone()

    if not cnt:
        # Nothing to delete
        cur.close()
        con.close()
        return

    # Get tables with a FOREIGN KEY to INTERPRO.METHOD
    tables = []
    child_tables = ora.get_child_tables(cur, "INTERPRO", "METHOD")
    for table, constraint, column in child_tables:
        tables.append((table, constraint, column))

    # Add  INTERPRO.METHOD as we want also to delete rows in this table
    tables.append(("METHOD", None, "METHOD_AC"))

    logger.info("disabling referential constraints")
    num_errors = 0
    for table, constraint, column in tables:
        if not constraint:
            continue

        # try:
        #     ora.toggle_constraint(cur, table, constraint, False)
        # except cx_Oracle.DatabaseError as exc:
        #     logger.error(exc)
        #     num_errors += 1

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
    with ThreadPoolExecutor(max_workers=threads) as executor:
        fs = {}

        for table, partition, column in tasks:
            f = executor.submit(delete_from_table, url, table, partition,
                                column)
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
