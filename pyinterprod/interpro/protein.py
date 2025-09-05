import os
import pickle
import sqlite3
from concurrent import futures
from dataclasses import dataclass
from tempfile import mkdtemp, mkstemp

import oracledb

from pyinterprod import logger
from pyinterprod.utils import Table, oracle as ora
from pyinterprod.uniprot import sprot
from .match import export_entries_protein_counts


@dataclass
class Sequence:
    accession: str
    identifier: str
    is_reviewed: bool
    crc64: str
    length: int
    is_fragment: bool
    taxon_id: int

    def __post_init__(self):
        if not isinstance(self.is_reviewed, bool):
            self.is_reviewed = self.is_reviewed in (1, 'S')

        if not isinstance(self.is_fragment, bool):
            self.is_fragment = self.is_fragment in (1, 'Y')

    @property
    def annotation(self):
        return (self.identifier, self.is_reviewed, self.length,
                self.is_fragment, self.taxon_id)

    def asdict(self):
        return dict(acc=self.accession,
                    identifer=self.identifier,
                    dbcode='S' if self.is_reviewed else 'T',
                    crc64=self.crc64,
                    length=self.length,
                    fragment='Y' if self.is_fragment else 'N',
                    taxid=self.taxon_id)


def export_proteins(url: str, outdir: str, buffer_size: int = 1000000) -> list[str]:
    con = oracledb.connect(url)
    cur = con.cursor()
    cur.execute(
        """
        SELECT PROTEIN_AC, NAME, DBCODE, CRC64, LEN, FRAGMENT, TAX_ID
        FROM INTERPRO.PROTEIN
        ORDER BY PROTEIN_AC
        """
    )

    buffer = {}
    files = []
    for row in cur:
        seq = Sequence(*row)
        buffer[seq.accession] = seq

        if len(buffer) == buffer_size:
            fd, file = mkstemp(dir=outdir)
            with open(fd, "wb") as fh:
                pickle.dump(buffer, fh)

            buffer = {}
            files.append(file)

    if buffer:
        fd, file = mkstemp(dir=outdir)
        with open(fd, "wb") as fh:
            pickle.dump(buffer, fh)

        files.append(file)

    cur.close()
    con.close()

    return files


def track_changes(url: str, swissp: str, trembl: str, version: str, date: str,
                  data_dir: str, tmpdir: str | None = None):
    os.makedirs(data_dir, exist_ok=True)
    workdir = mkdtemp(dir=tmpdir)

    con = oracledb.connect(url)
    cur = con.cursor()
    cur.execute("SELECT VERSION FROM INTERPRO.DB_VERSION WHERE DBCODE = 'u'")
    old_version, = cur.fetchone()

    logger.info(f"exporting protein counts per entry (UniProt {old_version})")
    export_entries_protein_counts(cur, data_dir)

    cur.close()
    con.close()

    logger.info(f"dumping UniProt {old_version} proteins")
    files = export_proteins(url, workdir)

    logger.info(f"loading UniProt {version} proteins")
    fd, database = mkstemp(dir=workdir)
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

    sprot.load(swissp, database, "protein")
    sprot.load(trembl, database, "protein")

    size = os.path.getsize(database) + sum(map(os.path.getsize, files))

    logger.info("tracking changes")
    con = oracledb.connect(url)
    cur = con.cursor()
    ora.truncate_table(cur, "INTERPRO.PROTEIN_CHANGES")
    ora.truncate_table(cur, "INTERPRO.PROTEIN_TO_DELETE")
    cur.close()

    # New proteins
    sql = """
        INSERT INTO INTERPRO.PROTEIN
        VALUES (:acc, :identifer, :dbcode, :crc64, :length, SYSDATE, USER,
                :fragment, 'N', :taxid)
    """
    new_proteins = Table(con, sql)

    # Annotation/sequence changes
    sql = """
        UPDATE INTERPRO.PROTEIN
        SET NAME = :identifer, DBCODE = :dbcode, CRC64 = :crc64, LEN = :length,
            TIMESTAMP = SYSDATE, USERSTAMP = USER, FRAGMENT = :fragment,
            TAX_ID = :taxid
        WHERE PROTEIN_AC = :acc
    """
    existing_proteins = Table(con, sql)

    # Obsolete proteins
    sql = "INSERT INTO INTERPRO.PROTEIN_TO_DELETE VALUES (:1, :2)"
    obsolete_proteins = Table(con, sql)

    # Proteins to track
    sql = "INSERT INTO INTERPRO.PROTEIN_CHANGES VALUES (:1)"
    track_proteins = Table(con, sql)

    old_reviewed = old_unreviewed = 0
    new_reviewed = new_unreviewed = 0

    min_acc = max_acc = None

    con2 = sqlite3.connect(database)
    cur2 = con2.cursor()
    for i, file in enumerate(files):
        with open(file, "rb") as fh:
            old_proteins = pickle.load(fh)

        os.remove(file)

        start = min(old_proteins)
        stop = max(old_proteins)
        if i == 0:
            cur2.execute(
                """
                SELECT *
                FROM protein
                WHERE accession >= ? AND accession <= ?
                """, (start, stop)
            )
            min_acc = start
        else:
            cur2.execute(
                """
                SELECT *
                FROM protein
                WHERE accession > ? AND accession <= ?
                """, (max_acc, stop)
            )
                       
        max_acc = stop

        for row in cur2:
            new_seq = Sequence(*row)

            if new_seq.is_reviewed:
                new_reviewed += 1
            else:
                new_unreviewed += 1

            try:
                old_seq = old_proteins.pop(new_seq.accession)
            except KeyError:
                new_proteins.insert(new_seq.asdict())
                track_proteins.insert((new_seq.accession,))
                continue

            if old_seq.is_reviewed:
                old_reviewed += 1
            else:
                old_unreviewed += 1

            if new_seq.crc64 != old_seq.crc64:
                # Sequence update
                existing_proteins.update(new_seq.asdict())

                # Track the protein (sequence change -> match changes)
                track_proteins.insert((new_seq.accession,))
            elif new_seq.annotation != old_seq.annotation:
                # Annotation update
                existing_proteins.update(new_seq.asdict())

        for old_seq in old_proteins.values():
            obsolete_proteins.insert((
                obsolete_proteins.count + 1,
                old_seq.accession
            ))

            if old_seq.is_reviewed:
                old_reviewed += 1
            else:
                old_unreviewed += 1

    """
    If there is a new protein with an accession lower than the lowest accession
    of the old proteins, or with an an accession greater than the greatest
    accession of the old proteins, it has not been considered until now
    """
    cur2.execute(
        """
        SELECT *
        FROM protein
        WHERE accession < ? OR accession > ?
        """, (min_acc, max_acc)
    )
    for row in cur2:
        new_seq = Sequence(*row)

        if new_seq.is_reviewed:
            new_reviewed += 1
        else:
            new_unreviewed += 1

        new_proteins.insert(new_seq.asdict())
        track_proteins.insert((new_seq.accession,))

    cur2.close()
    con2.close()
    os.remove(database)
    os.rmdir(workdir)

    new_proteins.close()
    existing_proteins.close()
    obsolete_proteins.close()
    track_proteins.close()

    cur = con.cursor()
    cur.executemany(
        """
        UPDATE INTERPRO.DB_VERSION
        SET
          VERSION = :1,
          ENTRY_COUNT = :2,
          FILE_DATE = TO_DATE(:3, 'DD-Mon-YYYY'),
          LOAD_DATE = SYSDATE
          WHERE DBCODE = :4
        """, [
            (version, new_reviewed, date, 'S'),
            (version, new_unreviewed, date, 'T'),
            (version, new_reviewed + new_unreviewed, date, 'u')
        ]
    )
    con.commit()

    ora.gather_stats(cur, "INTERPRO", "PROTEIN_CHANGES")
    ora.gather_stats(cur, "INTERPRO", "PROTEIN_TO_DELETE")
    cur.close()
    con.close()

    logger.info(f"Reviewed (before):     {old_reviewed:>12,}")
    logger.info(f"Reviewed (now):        {new_reviewed:>12,}")
    logger.info(f"Unreviewed (before):   {old_unreviewed:>12,}")
    logger.info(f"Unreviewed (now):      {new_unreviewed:>12,}")
    logger.info(f"New proteins:          {new_proteins.count:>12,}")
    logger.info(f"Updated proteins:      {existing_proteins.count:>12,}")
    logger.info(f"Obsolete proteins:     {obsolete_proteins.count:>12,}")
    logger.info(f"Disk usage:            {size / 1024 ** 2:>9,.0f} MB")

    if new_reviewed == 0 or new_unreviewed == 0:
        raise RuntimeError("No reviewed or unreviewed UniProtKB entries")


def delete_obsoletes(url: str, truncate: bool = False, threads: int = 8,
                     step: int = 10000):
    con = oracledb.connect(url)
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
        except oracledb.DatabaseError as exc:
            logger.error(exc)
            num_errors += 1

    if num_errors:
        cur.close()
        con.close()
        raise RuntimeError(f"{num_errors} constraints could not be disabled")

    # Rebuild indexes (DELETE fails if an index is unusable)
    for table, constraint, column in tables:
        for index in ora.get_indexes(cur, "INTERPRO", table):
            if index["is_unusable"]:
                logger.info(f"rebuilding index {index['name']}")
                ora.rebuild_index(cur, index["name"])

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
                num_rows = f.result()
            except Exception as exc:
                logger.info(f"{name}: failed ({exc})")
                num_errors += 1
            else:
                logger.info(f"{name}: {num_rows:,} rows deleted")

        if num_errors:
            raise RuntimeError(f"{num_errors} tables failed")

    logger.info("enabling referential constraints")
    con = oracledb.connect(url)
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
        except oracledb.DatabaseError as exc:
            logger.error(exc)
            num_errors += 1

    if num_errors:
        cur.close()
        con.close()
        raise RuntimeError(f"{num_errors} constraints could not be enabled")

    for table, constraint, column in tables:
        for index in ora.get_indexes(cur, "INTERPRO", table):
            if index["is_unusable"]:
                logger.info(f"rebuilding index {index['name']}")
                ora.rebuild_index(cur, index["name"])

    cur.close()
    con.close()
    logger.info("complete")


def iterative_delete(url: str, table: str, partition: str | None,
                     column: str, step: int, stop: int,
                     gather_stats: bool = True) -> int:
    con = oracledb.connect(url)
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
        return num_rows

    for i in range(1, stop, step):
        cur.execute(
            f"""
            DELETE FROM INTERPRO.{_table}
            WHERE {column} IN (
              SELECT PROTEIN_AC
              FROM INTERPRO.PROTEIN_TO_DELETE
              WHERE ID BETWEEN :1 and :2
            )
            """,
            [i, i + step - 1]
        )
        con.commit()
        logger.debug(f"{_table}: {i + step - 1:,} / {stop:,}")

    if gather_stats:
        ora.gather_stats(cur, "INTERPRO", table, partition)

    cur.close()
    con.close()
    return num_rows


def check_proteins(cur: oracledb.Cursor) -> int:
    num_errors = 0
    cur.execute(
        """
        SELECT PROTEIN_AC
        FROM INTERPRO.PROTEIN
        MINUS
        SELECT AC
        FROM UNIPARC.XREF
        WHERE DBID IN (2, 3)   -- Swiss-Prot/TrEMBL
        AND DELETED = 'N'      -- Not deleted
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
    con = oracledb.connect(url)
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
