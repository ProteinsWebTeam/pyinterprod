# -*- coding: utf-8 -*-

import os
import pickle
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Dict, Mapping, Sequence, Set, Tuple

import cx_Oracle
import psycopg2

from pyinterprod import logger
from pyinterprod.utils import Table, oracle as ora, pg
from . import contrib
from .database import Database


SIGNATURES_DESCR_FILE = "swissprot_descr.dat"
MATCH_PARTITIONS = {
    'B': 'MATCH_DBCODE_B',  # SFLD
    'F': 'MATCH_DBCODE_F',
    'H': 'MATCH_DBCODE_H',  # Pfam
    'J': 'MATCH_DBCODE_J',  # CDD
    'M': 'MATCH_DBCODE_M',
    'N': 'MATCH_DBCODE_N',
    'P': 'MATCH_DBCODE_P',
    'Q': 'MATCH_DBCODE_Q',
    'R': 'MATCH_DBCODE_R',
    'U': 'MATCH_DBCODE_U',
    'V': 'MATCH_DBCODE_V',  # PANTHER
    'X': 'MATCH_DBCODE_X',
    'Y': 'MATCH_DBCODE_Y',
}
SITE_PARTITIONS = {
    'B': 'SFLD',
    'J': 'CDD',
}


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
            if db.identifier == 'V':
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


def make_pre_report(url: str, databases: Sequence[Database], output: str):
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

    with open(output, "wb") as fh:
        pickle.dump(results, fh)


def delete_from_table(url: str, table: str, column: str, step: int,
                      stop: int) -> int:
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.execute(
        f"""
        SELECT COUNT(*)
        FROM {table}
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
            DELETE FROM INTERPRO.{table}
            WHERE {column} IN (
              SELECT METHOD_AC
              FROM INTERPRO.METHOD_TO_DELETE
              WHERE ID BETWEEN :1 and :2
            )
            """, (i, i + step - 1)
        )

    con.commit()
    ora.gather_stats(cur, "INTERPRO", table)
    cur.close()
    con.close()
    return num_rows


def create_and_exchange(url: str, table: str, partition: str, column: str) -> int:
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.execute(
        f"""
        SELECT COUNT(*)
        FROM {table} PARTITION ({partition})
        WHERE {column} IN (
            SELECT METHOD_AC 
            FROM INTERPRO.METHOD_TO_DELETE
        ) 
        """
    )
    num_rows, = cur.fetchone()

    logger.debug(f"{table} ({partition}): {num_rows:,} rows to delete")
    if not num_rows:
        cur.close()
        con.close()
        return num_rows

    # Create table without obsolete signatures
    logger.debug(f"{table} ({partition}): creating temporary table")
    tmp_table = f"{partition}_TMP"
    ora.drop_table(cur, tmp_table)
    cur.execute(
        f"""
        CREATE TABLE {tmp_table} NOLOGGING
        AS  
        SELECT *
        FROM {table} PARTITION ({partition})
        WHERE {column} NOT IN (
            SELECT METHOD_AC
            FROM INTERPRO.METHOD_TO_DELETE
        )
        """
    )

    # Add constraints/indexes to be able to exchange partition
    logger.debug(f"{table} ({partition}): creating indexes and constraints")
    cur.execute(f"CREATE INDEX {tmp_table}$D ON {tmp_table} (DBCODE) NOLOGGING")
    cur.execute(f"CREATE INDEX {tmp_table}$E ON {tmp_table} (EVIDENCE) NOLOGGING")
    cur.execute(f"CREATE INDEX {tmp_table}$S ON {tmp_table} (STATUS) NOLOGGING")
    cur.execute(f"CREATE INDEX {tmp_table}$M ON {tmp_table} (METHOD_AC) NOLOGGING")
    cur.execute(
        f"""
        ALTER TABLE {tmp_table}
        ADD CONSTRAINT {tmp_table}$CK1 
        CHECK ( POS_FROM >= 1 ) 
        """
    )
    cur.execute(
        f"""
        ALTER TABLE {tmp_table}
        ADD CONSTRAINT {tmp_table}$CK2
        CHECK ( POS_TO - POS_FROM > 0 ) 
        """
    )
    cur.execute(
        f"""
        ALTER TABLE {tmp_table}
        ADD CONSTRAINT {tmp_table}$CK3
        CHECK ( STATUS != 'N' OR (STATUS = 'N' AND DBCODE IN ('P', 'M', 'Q')) )
        """
    )
    cur.execute(
        f"""
        ALTER TABLE {tmp_table}
        ADD CONSTRAINT {tmp_table}$PK
        PRIMARY KEY (PROTEIN_AC, METHOD_AC, POS_FROM, POS_TO)
        """
    )
    cur.execute(
        f"""
        ALTER TABLE {tmp_table}
        ADD CONSTRAINT {tmp_table}$FK1 FOREIGN KEY (DBCODE) 
        REFERENCES INTERPRO.CV_DATABASE (DBCODE)
        """
    )
    cur.execute(
        f"""
        ALTER TABLE {tmp_table}
        ADD CONSTRAINT {tmp_table}$FK2 FOREIGN KEY (EVIDENCE) 
        REFERENCES INTERPRO.CV_EVIDENCE (CODE)
        """
    )
    """
    Do not create a FK to INTERPRO.METHOD as creating one would raise
        ORA-14128 (FOREIGN KEY constraint mismatch 
                   in ALTER TABLE EXCHANGE PARTITION)
    when exchanging partition.
    INTERPRO.MATCH *has* a FK to INTERPRO.METHOD, but we disabled it before
    starting to delete obsolete signatures
    """
    # cur.execute(
    #     f"""
    #     ALTER TABLE {tmp_table}
    #     ADD CONSTRAINT {tmp_table}$FK3 FOREIGN KEY (METHOD_AC)
    #     REFERENCES INTERPRO.METHOD (METHOD_AC)
    #     """
    # )
    cur.execute(
        f"""
        ALTER TABLE {tmp_table}
        ADD CONSTRAINT {tmp_table}$FK4 FOREIGN KEY (PROTEIN_AC) 
        REFERENCES INTERPRO.PROTEIN (PROTEIN_AC)
        """
    )
    cur.execute(
        f"""
        ALTER TABLE {tmp_table}
        ADD CONSTRAINT {tmp_table}$FK5 FOREIGN KEY (STATUS) 
        REFERENCES INTERPRO.CV_STATUS (CODE)
        """
    )
    cur.execute(
        f"""
        ALTER TABLE {tmp_table}
        ADD CONSTRAINT {tmp_table}$CK4 
        CHECK (PROTEIN_AC IS NOT NULL )
        """
    )
    cur.execute(
        f"""
        ALTER TABLE {tmp_table}
        ADD CONSTRAINT {tmp_table}$CK5 
        CHECK (METHOD_AC IS NOT NULL )
        """
    )

    logger.debug(f"{table} ({partition}): exchanging partition")
    cur.execute(
        f"""
        ALTER TABLE {table} 
        EXCHANGE PARTITION ({partition}) 
        WITH TABLE {tmp_table}
        WITHOUT VALIDATION
        """
    )
    ora.drop_table(cur, tmp_table, purge=True)

    cur.close()
    con.close()
    return num_rows


def delete_obsolete(url: str, databases: Sequence[Database], **kwargs):
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

    # MATCH and SITE_MATCH are too big for DELETE: exchange partitions
    delete_tasks = []
    exchange_tasks = []
    for table, constraint, column in tables:
        if table == "MATCH":
            for db in databases:
                partition = MATCH_PARTITIONS[db.identifier]
                exchange_tasks.append((table, partition, column))
        elif table == "SITE_MATCH":
            for db in databases:
                partition = SITE_PARTITIONS[db.identifier]
                exchange_tasks.append((table, partition, column))
        else:
            # Normal DELETE statement
            delete_tasks.append((table, column))

    cur.close()
    con.close()

    with ThreadPoolExecutor(max_workers=threads) as executor:
        fs = {}

        for table, column in delete_tasks:
            args = (url, table, column, step, stop)
            f = executor.submit(delete_from_table, *args)
            fs[f] = (table, None)

        for table, partition, column in exchange_tasks:
            args = (url, table, partition, column)
            f = executor.submit(create_and_exchange, *args)
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
