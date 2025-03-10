import heapq
import os
import shutil
import tarfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from logging import DEBUG
from tempfile import mkdtemp

import oracledb

from pyinterprod import logger
from pyinterprod.utils import Table
from pyinterprod.utils.io import dump, iter_until_eof
from pyinterprod.utils.oracle import drop_table, get_partitions


def load_matches(
    uri: str, databases: dict[str, str], tmpdir: str | None = None, processes: int = 1
):
    """
    :param uri: Oracle database connection string
    :param databases: dictionary of databases to update
                      key -> dbcode
                      value -> path to tar archive
    :param tmpdir: directory for temporary files
    :param processes: number of parallel processes
    """
    con = oracledb.connect(uri)
    cur = con.cursor()
    partitions = {}
    for p in get_partitions(cur, "INTERPRO", "TOAD_MATCH"):
        name = p["name"]
        dbcode = p["value"][1:-1]  # 'X' -> X
        partitions[dbcode] = name

    cur.close()
    con.close()

    tasks = []
    for dbcode, filepath in databases.items():
        try:
            partition = partitions[dbcode]
        except KeyError:
            raise KeyError(f"No partition for database with dbcode '{dbcode}'")
        else:
            tasks.append((dbcode, partition, filepath))

    if processes > 1:
        logger.setLevel(DEBUG)

        for i, (dbcode, partition, filepath) in enumerate(tasks):
            logger.info(f"updating partition: {partition}")
            last = i + 1 == len(databases)
            load_database_matches(uri, partition, filepath, tmpdir=tmpdir, purge=last)
    else:
        with ProcessPoolExecutor(max_workers=processes) as executor:
            fs = {}
            for dbcode, partition, filepath in tasks:
                f = executor.submit(
                    load_database_matches,
                    uri,
                    partition,
                    filepath,
                    tmpdir=tmpdir,
                    suffix=f"_{dbcode}",
                )
                fs[f] = partition

            errors = 0
            for f in as_completed(fs):
                partition = fs[f]

                try:
                    f.result()
                except Exception as exc:
                    logger.error(f"{partition}: {exc}")
                    errors += 1
                else:
                    logger.info(f"{partition}: done")

            if errors:
                raise RuntimeError(f"{errors} errors occurred")


def load_database_matches(
    uri: str,
    partition: str,
    filepath: str,
    tmpdir: str | None = None,
    suffix: str = "",
    purge: bool = False,
):
    if tmpdir:
        os.makedirs(tmpdir, exist_ok=True)

    outdir = mkdtemp(dir=tmpdir)

    logger.debug(f"\textracting matches from {filepath}")
    files = extract_matches(filepath, outdir)

    logger.debug(f"\tinsert matches")
    con = oracledb.connect(uri)
    cur = con.cursor()
    drop_table(cur, f"INTERPRO.TOAD_MATCH_NEW{suffix}", purge=purge)
    cur.execute(
        f"""
        CREATE TABLE INTERPRO.TOAD_MATCH_NEW{suffix} NOLOGGING
        AS SELECT * FROM INTERPRO.TOAD_MATCH WHERE 1 = 0
        """
    )

    query = f"""
        INSERT /*+ APPEND */ 
        INTO INTERPRO.TOAD_MATCH_NEW{suffix}
        VALUES (:1, :2, :3, :4, :5, :6, :7)
    """
    with Table(
        con=cur.connection, query=query, autocommit=True, buffer_size=1000
    ) as table:
        for uniprot_acc, method_acc, pos_from, pos_to, group_id, score in iter_matches(
            files
        ):
            table.insert(
                (
                    uniprot_acc,
                    method_acc,
                    "-",  # DBCODE (cannot be null but value isn't important)
                    pos_from,
                    pos_to,
                    group_id,
                    score,
                )
            )
    logger.debug(f"\t{table.count:,} matches inserted")

    shutil.rmtree(outdir)

    # Filter matches (deleted proteins/signatures or out-of-bounds)
    logger.debug("\tfiltering matches")
    drop_table(cur, f"INTERPRO.TOAD_MATCH_TMP{suffix}", purge=purge)
    cur.execute(
        f"""
        CREATE TABLE INTERPRO.TOAD_MATCH_TMP{suffix} NOLOGGING
        AS SELECT * FROM INTERPRO.TOAD_MATCH WHERE 1 = 0
        """
    )
    cur.execute(
        f"""
        INSERT /*+ APPEND */ INTO INTERPRO.TOAD_MATCH_TMP{suffix}
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
            FROM INTERPRO.TOAD_MATCH_NEW{suffix} T
            INNER JOIN INTERPRO.PROTEIN P ON T.PROTEIN_AC = P.PROTEIN_AC
            INNER JOIN INTERPRO.METHOD M ON T.METHOD_AC = M.METHOD_AC
        )
        WHERE POS_FROM <= POS_TO
        """,
    )
    cnt = cur.rowcount
    cur.connection.commit()

    drop_table(cur, f"INTERPRO.TOAD_MATCH_NEW{suffix}", purge=purge)
    logger.debug(f"\t{cnt:,} matches kept")

    logger.debug("\tcreating indexes and constraints")
    cur.execute(
        f"""
        ALTER TABLE INTERPRO.TOAD_MATCH_TMP{suffix}
        ADD CONSTRAINT CK_TOAD_MATCH_TMP{suffix}$FROM
        CHECK (POS_FROM >= 1)
        """
    )
    cur.execute(
        f"""
        ALTER TABLE INTERPRO.TOAD_MATCH_TMP{suffix}
        ADD CONSTRAINT CK_TOAD_MATCH_TMP{suffix}$TO 
        CHECK (POS_FROM <= POS_TO)
        """
    )
    cur.execute(
        f"""
        ALTER TABLE INTERPRO.TOAD_MATCH_TMP{suffix}
        ADD CONSTRAINT PK_TOAD_MATCH_TMP{suffix}
        PRIMARY KEY (PROTEIN_AC, METHOD_AC, DBCODE, POS_FROM, POS_TO)
        """
    )
    cur.execute(
        f"""
        ALTER TABLE INTERPRO.TOAD_MATCH_TMP{suffix}
        ADD CONSTRAINT FK_TOAD_MATCH_TMP{suffix}$PROTEIN 
        FOREIGN KEY (PROTEIN_AC) REFERENCES INTERPRO.PROTEIN (PROTEIN_AC)
        """
    )
    cur.execute(
        f"""
        ALTER TABLE INTERPRO.TOAD_MATCH_TMP{suffix}
        ADD CONSTRAINT FK_TOAD_MATCH_TMP{suffix}$METHOD 
        FOREIGN KEY (METHOD_AC) REFERENCES INTERPRO.METHOD (METHOD_AC)
        """
    )
    cur.execute(
        f"""
        ALTER TABLE INTERPRO.TOAD_MATCH_TMP{suffix}
        ADD CONSTRAINT FK_TOAD_MATCH_TMP{suffix}$DBCODE 
        FOREIGN KEY (DBCODE) REFERENCES INTERPRO.CV_DATABASE (DBCODE)
        """
    )

    logger.debug("\texchanging partition")
    cur.execute(
        f"""
        ALTER TABLE INTERPRO.TOAD_MATCH
        EXCHANGE PARTITION ({partition})
        WITH TABLE INTERPRO.TOAD_MATCH_TMP{suffix}
        """
    )

    drop_table(cur, f"INTERPRO.TOAD_MATCH_TMP{suffix}", purge=purge)
    logger.debug("\tdone")


def extract_matches(
    filepath: str, outdir: str | None, buffersize: int = 1000000
) -> list[str]:
    files = []
    proteins = {}
    i = 0
    for uniprot_acc, method_acc, fragments, score in iter_tarfile(filepath):
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


def iter_tarfile(filepath: str):
    """
    Extract TOAD matches from TAR archive
    :param filepath: Path to TAR archive
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
                        yield (
                            uniprot_acc,
                            signature_acc,
                            [(int(start), int(end))],
                            float(score),
                        )
                    elif len(values) == 23:
                        # Discontinuous domains (up to ten fragments)
                        uniprot_acc, signature_acc = values[:2]
                        fragments = []
                        for i in range(10):
                            start = values[2 + i * 2]
                            end = values[3 + i * 2]
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


def iter_matches(files: list[str]):
    """
    Process and iterate over TOAD matches
    :param files: list of file paths to extracted TOAD matches
    """
    for protein_acc, matches in process_matches(files):
        for signature_acc, locations in matches.items():
            for i, (_, _, fragments, score) in enumerate(locations):
                for pos_from, pos_to in fragments:
                    yield (protein_acc, signature_acc, pos_from, pos_to, i + 1, score)


def process_matches(files: list[str]):
    """
    Process extracted TOAD matches to group them by UniProt accession
    :param files: list of file paths to extracted TOAD matches
    """
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
    """
    Evaluate TOAD matches to keep the best non-overlapping predictions
    :param matches: list of predicted matches
    :return: dictionary of matches group by member database accession
    """
    profiles = {}
    for acc, fragments, score in matches:
        try:
            hits = profiles[acc]
        except KeyError:
            hits = profiles[acc] = []

        hits.append((fragments, score))

    for acc, hits in profiles.items():
        # Create locations from fragments
        locations = []
        for fragments, score in hits:
            # fragments are sorted when extracting matches from tar archive,
            # see `iter_tarfile()`
            pos_from = fragments[0][0]
            pos_to = max(end for start, end in fragments)
            locations.append((pos_from, pos_to, fragments, score))

        # Sort by score (descending)
        locations.sort(key=lambda x: -x[3])

        # When two locations overlap, keep the one with the highest score
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
