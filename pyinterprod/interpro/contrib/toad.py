import heapq
import os
import shutil
import tarfile
from tempfile import mkdtemp

import oracledb

from pyinterprod import logger
from pyinterprod.utils import Table
from pyinterprod.utils.io import dump, iter_until_eof
from pyinterprod.utils.oracle import drop_table, get_partitions


def load_matches(uri: str, databases: dict[str, str], **kwargs):
    con = oracledb.connect(uri)
    cur = con.cursor()

    partitions = {}
    for p in get_partitions(cur, "INTERPRO", "TOAD_MATCH"):
        name = p["name"]
        dbcode = p["value"][1:-1]  # 'X' -> X
        partitions[dbcode] = name

    for i, (dbshort, filepath) in enumerate(databases.items()):
        cur.execute(
            """
            SELECT DBCODE, DBNAME 
            FROM INTERPRO.CV_DATABASE 
            WHERE DBSHORT = :1
            """,
            [dbshort.upper()]
        )
        row = cur.fetchone()
        if row is None:
            cur.close()
            con.close()
            raise ValueError(f"No database found for {dbshort}")

        dbcode, dbname = row

        try:
            partition = partitions[dbcode]
        except KeyError:
            cur.close()
            con.close()
            err = f"No partition in TOAD_MATCH for database {dbshort}"
            raise KeyError(err)

        logger.info(f"loading {dbname} matches")
        last = i + 1 == len(databases)
        load_database_matches(cur, partition, filepath,
                              tmpdir=kwargs.get("tmpdir"),
                              purge=last)

    cur.close()
    con.close()


def load_database_matches(cur: oracledb.Cursor, partition: str, filepath: str,
                          tmpdir: str | None = None, purge: bool = False):
    logger.info(f"\tloading matches from {filepath}")
    # Insert "raw" matches (as provided by DeepMind)
    drop_table(cur, "INTERPRO.TOAD_MATCH_NEW", purge=purge)
    cur.execute(
        """
        CREATE TABLE INTERPRO.TOAD_MATCH_NEW NOLOGGING
        AS SELECT * FROM INTERPRO.TOAD_MATCH WHERE 1 = 0
        """
    )

    query = """
        INSERT /*+ APPEND */ 
        INTO INTERPRO.TOAD_MATCH_NEW 
        VALUES (:1, :2, :3, :4, :5, :6, :7)
    """
    with Table(con=cur.connection, query=query, autocommit=True) as table:
        for item in process_tarfile(filepath, tmpdir=tmpdir):
            uniprot_acc, method_acc, pos_from, pos_to, group_id, score = item
            table.insert((
                uniprot_acc,
                method_acc,
                '-',  # DBCODE (cannot be null but value isn't important)
                pos_from,
                pos_to,
                group_id,
                score
            ))
    logger.info(f"\t{table.count:,} matches inserted")

    # Filter matches (deleted proteins/signatures or out-of-bounds)
    logger.info("\tfiltering matches")
    drop_table(cur, "INTERPRO.TOAD_MATCH_TMP", purge=purge)
    cur.execute(
        """
        CREATE TABLE INTERPRO.TOAD_MATCH_TMP NOLOGGING
        AS SELECT * FROM INTERPRO.TOAD_MATCH WHERE 1 = 0
        """
    )
    cur.execute(
        """
        INSERT /*+ APPEND */ INTO INTERPRO.TOAD_MATCH_TMP
        SELECT *
        FROM (
            SELECT T.PROTEIN_AC,
                   T.METHOD_AC,
                   M.DBCODE,
                   T.POS_FROM,
                   CASE WHEN T.POS_TO <= P.LEN
                        THEN T.POS_TO 
                        ELSE P.LEN 
                        END POS_TO,
                   T.GROUP_ID,
                   T.SCORE 
            FROM INTERPRO.TOAD_MATCH_NEW T
            INNER JOIN INTERPRO.PROTEIN P ON T.PROTEIN_AC = P.PROTEIN_AC
            INNER JOIN INTERPRO.METHOD M ON T.METHOD_AC = M.METHOD_AC
        )
        WHERE POS_FROM <= POS_TO
        """,
    )
    cnt = cur.rowcount
    cur.connection.commit()

    drop_table(cur, "INTERPRO.TOAD_MATCH_NEW", purge=purge)
    logger.info(f"\t{cnt:,} matches kept")

    logger.info("\tcreating indexes and constraints")
    cur.execute(
        """
        ALTER TABLE INTERPRO.TOAD_MATCH_TMP
        ADD CONSTRAINT CK_TOAD_MATCH_TMP$FROM 
        CHECK (POS_FROM >= 1)
        """
    )
    cur.execute(
        """
        ALTER TABLE INTERPRO.TOAD_MATCH_TMP
        ADD CONSTRAINT CK_TOAD_MATCH_TMP$TO 
        CHECK (POS_FROM <= POS_TO)
        """
    )
    cur.execute(
        """
        ALTER TABLE INTERPRO.TOAD_MATCH_TMP
        ADD CONSTRAINT PK_TOAD_MATCH_TMP
        PRIMARY KEY (PROTEIN_AC, METHOD_AC, DBCODE, POS_FROM, POS_TO)
        """
    )
    cur.execute(
        """
        ALTER TABLE INTERPRO.TOAD_MATCH_TMP
        ADD CONSTRAINT FK_TOAD_MATCH_TMP$PROTEIN 
        FOREIGN KEY (PROTEIN_AC) REFERENCES INTERPRO.PROTEIN (PROTEIN_AC)
        """
    )
    cur.execute(
        """
        ALTER TABLE INTERPRO.TOAD_MATCH_TMP
        ADD CONSTRAINT FK_TOAD_MATCH_TMP$METHOD 
        FOREIGN KEY (METHOD_AC) REFERENCES INTERPRO.METHOD (METHOD_AC)
        """
    )
    cur.execute(
        """
        ALTER TABLE INTERPRO.TOAD_MATCH_TMP
        ADD CONSTRAINT FK_TOAD_MATCH_TMP$DBCODE 
        FOREIGN KEY (DBCODE) REFERENCES INTERPRO.CV_DATABASE (DBCODE)
        """
    )

    logger.info("\texchanging partition")
    cur.execute(
        f"""
        ALTER TABLE INTERPRO.TOAD_MATCH
        EXCHANGE PARTITION ({partition})
        WITH TABLE INTERPRO.TOAD_MATCH_TMP
        """
    )

    drop_table(cur, "INTERPRO.TOAD_MATCH_TMP", purge=purge)


def process_tarfile(filepath: str, tmpdir: str | None, **kwargs):
    if tmpdir:
        os.makedirs(tmpdir, exist_ok=True)

    outdir = mkdtemp(dir=tmpdir)
    files = parse_matches(filepath, outdir, **kwargs)

    for protein_acc, matches in process_matches(files):
        for signature_acc, locations in matches.items():
            for i, (_, _, fragments, score) in enumerate(locations):
                for pos_from, pos_to in fragments:
                    yield (
                        protein_acc,
                        signature_acc,
                        pos_from,
                        pos_to,
                        i + 1,
                        score
                    )

    shutil.rmtree(outdir)


def parse_matches(filepath: str, outdir: str | None,
                  buffersize: int = 1000000) -> list[str]:
    files = []
    proteins = {}
    i = 0
    for uniprot_acc, method_acc, fragments, score in iter_matches(filepath):
        try:
            matches = proteins[uniprot_acc]
        except KeyError:
            matches = proteins[uniprot_acc] = []

        matches.append((method_acc, fragments, score))

        i += 1
        if i % buffersize == 0:
            file = dump(proteins, outdir)
            files.append(file)
            proteins.clear()

    if proteins:
        file = dump(proteins, outdir)
        files.append(file)
        proteins.clear()

    return files


def process_matches(files: list[str]):
    iterable = [iter_until_eof(f) for f in files]
    protein_acc = None
    matches = []
    for key, value in heapq.merge(*iterable, key=lambda x: x[0]):
        if key != protein_acc:
            if matches:
                yield protein_acc, filter_matches(matches)

            protein_acc = key
            matches.clear()

        matches += value

    if matches:
        yield protein_acc, filter_matches(matches)


def filter_matches(matches: list[tuple]) -> dict[str, list[tuple]]:
    profiles = {}
    for acc, fragments, score in matches:
        try:
            hits = profiles[acc]
        except KeyError:
            hits = profiles[acc] = []

        hits.append((fragments, score))

    for acc, hits in profiles.items():
        locations = []
        for fragments, score in hits:
            pos_from = fragments[0][0]
            pos_to = max(end for start, end in fragments)
            locations.append((pos_from, pos_to, fragments, score))

        # Sort by score (descending)
        locations.sort(key=lambda x: -x[3])

        filtered_locations = []
        for loc in sorted(locations, key=lambda x: -x[3]):
            for other in filtered_locations:
                overlap = min(loc[1], other[1]) - max(loc[0], other[0])
                if overlap >= 0:
                    # Do not keep (overlaps with one having a better score)
                    break
            else:
                filtered_locations.append(loc)

        profiles[acc] = filtered_locations

    return profiles


def iter_matches(filepath: str):
    """Iterate TOAD inferences
    """
    with tarfile.open(filepath, mode="r") as tar:
        for member in tar:
            if member.name.endswith(".tsv"):
                br = tar.extractfile(member)
                lines = br.read().decode("utf-8").splitlines(keepends=False)

                # First line is a header
                for line in lines[1:]:
                    values = line.split("\t")
                    if len(values) == 5:
                        # No discontinuous domains
                        uniprot_acc, signature_acc, start, end, score = values
                        yield (uniprot_acc, signature_acc,
                               [(int(start), int(end))], float(score))
                    elif len(values) == 23:
                        # Discontinuous domains (up to ten fragments)
                        uniprot_acc, signature_acc = values[:2]
                        fragments = []
                        for i in range(10):
                            start = values[2+i*2]
                            end = values[3+i*2]
                            if start == "NULL":
                                break
                            else:
                                fragments.append((int(start), int(end)))

                        score = float(values[-1])
                        fragments.sort()
                        yield uniprot_acc, signature_acc, fragments, score
                    else:
                        err = f"Unexpected number of columns: {values}"
                        raise ValueError(err)
