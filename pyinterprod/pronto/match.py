import hashlib
import heapq
import os
import pickle
import shutil
from multiprocessing import Process, Queue
from tempfile import mkstemp
from typing import Optional

import cx_Oracle
import psycopg2

from pyinterprod import logger
from pyinterprod.utils import pg
from pyinterprod.utils.kvdb import KVdb


"""
At least 50% of the residues of the shortest signature
  must overlap the other signature
(shorted signature = signature with the least residues in the protein)
"""
_MIN_OVERLAP = 0.5

"""
One of the signatures must hit at least 50% of the proteins hit by the other
signature
"""
_MIN_COLLOCATION = 0.5

# Threshold for Jaccard index/coefficients
_MIN_SIMILARITY = 0.75

# Domain org.: introduce a gap when distance between two positions > 20 aa
_MAX_GAP = 20

_INDEX_SUFFIX = ".i"


def export(url: str, output: str, cachesize: int = 10000000,
           tmpdir: Optional[str] = None):
    if tmpdir:
        os.makedirs(tmpdir, exist_ok=True)

    logger.info("exporting matches")
    files = _export_matches(url, cachesize, tmpdir)

    logger.info("sorting proteins")
    index = []
    with open(output, "wb") as fh:
        i = o = 0
        for protein in _merge_files(files):
            pickle.dump(protein, fh)

            i += 1
            if i % cachesize == 0:
                index.append((o, cachesize))
                o = fh.tell()

            if i % 10e6 == 0:
                logger.info(f"{i:>15,}")

        if i % cachesize != 0:
            index.append((o, i % cachesize))

    size = 0
    for file in files:
        size += os.path.getsize(file)
        os.unlink(file)

    logger.info(f"temporary files: {size / 1024 ** 2:.0f} MB")

    # Write index
    with open(f"{output}{_INDEX_SUFFIX}", "wb") as fh:
        pickle.dump(index, fh)

    logger.info("done")


def _export_matches(url: str, cachesize: int,
                    tmpdir: Optional[str] = None) -> list[str]:
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.execute(
        """
        SELECT P.PROTEIN_AC, P.DBCODE, P.FRAGMENT, E.LEFT_NUMBER, 
               LOWER(D.DBSHORT), M.METHOD_AC, M.FRAGMENTS, M.POS_FROM, 
               M.POS_TO, M.IS_METHOD
        FROM (
            SELECT PROTEIN_AC, METHOD_AC, DBCODE,
                   POS_FROM, POS_TO, FRAGMENTS, 'Y' IS_METHOD
            FROM INTERPRO.MATCH
            UNION ALL
            SELECT PROTEIN_AC, METHOD_AC, DBCODE, 
                   POS_FROM, POS_TO, NULL, 'N' IS_METHOD
            FROM INTERPRO.FEATURE_MATCH FM
            WHERE FM.DBCODE = 'a'  -- AntiFam
        ) M
        INNER JOIN INTERPRO.CV_DATABASE D 
            ON M.DBCODE = D.DBCODE
        INNER JOIN INTERPRO.PROTEIN P 
            ON P.PROTEIN_AC = M.PROTEIN_AC
        INNER JOIN INTERPRO.ETAXI E
            ON E.TAX_ID = P.TAX_ID
        """
    )

    i = 0
    cache = {}
    files = []
    for row in cur:
        protein_acc = row[0]

        try:
            obj = cache[protein_acc]
        except KeyError:
            obj = cache[protein_acc] = (
                row[1] == "S",  # is reviewed
                row[2] == "N",  # is complete
                row[3],  # taxon left number
                []  # matches
            )

        if row[6]:
            fragments = row[6]
        else:
            # Single/continuous fragment using match start/end positions
            fragments = f"{row[7]}-{row[8]}-S"

        obj[2].append((
            row[5],  # match accession
            row[4],  # match DB
            row[9] == "Y",  # is signature?
            fragments
        ))

        i += 1
        if i % cachesize == 0:
            fd, file = mkstemp(dir=tmpdir)
            _dump(cache, fd)
            cache.clear()
            files.append(file)

        if i % 1e8 == 0:
            logger.info(f"{i:>15,}")

    cur.close()
    con.close()

    fd, file = mkstemp(dir=tmpdir)
    _dump(cache, fd)
    cache.clear()
    files.append(file)

    logger.info(f"{i:>15,}")

    return files


def _dump(data: dict, fd: int):
    with open(fd, "wb") as fh:
        for key in sorted(data):
            pickle.dump((key, data[key]), fh)


def _merge_files(files: list[str]):
    iterable = [_load(file) for file in files]
    protein_acc = is_reviewed = is_complete = left_number = None
    matches = {}
    for key, value in heapq.merge(*iterable, key=lambda x: x[0]):
        if key != protein_acc:
            if protein_acc:
                yield (protein_acc, is_reviewed, is_complete, left_number,
                       matches)

            protein_acc = key
            is_reviewed, is_complete, left_number, _ = value
            matches = {}

        for match_acc, match_db, is_signature, fragments in value[3]:
            if match_acc in matches:
                matches[match_acc][2].append(fragments)
            else:
                matches[match_acc] = (match_db, is_signature, [fragments])

    yield protein_acc, is_reviewed, is_complete, left_number, matches


def _load(file: str):
    with open(file, "rb") as fh:
        while True:
            try:
                yield pickle.load(fh)
            except EOFError:
                break


def insert_signature2protein(url: str, names_db: str, matches_file: str,
                             processes: int = 1, tmpdir: Optional[str] = None):
    if tmpdir:
        os.makedirs(tmpdir, exist_ok=True)

    logger.info("copying SQLite database")
    fd, tmp_db = mkstemp(dir=tmpdir)
    os.close(fd)
    shutil.copyfile(names_db, tmp_db)

    logger.info("creating signature2protein")
    con = psycopg2.connect(**pg.url2dict(url))
    with con.cursor() as cur:
        cur.execute("DROP TABLE IF EXISTS signature2protein")
        cur.execute(
            """
            CREATE TABLE signature2protein (
                signature_acc VARCHAR(25) NOT NULL,
                protein_acc VARCHAR(15) NOT NULL,
                is_reviewed BOOLEAN NOT NULL,
                taxon_left_num INTEGER NOT NULL,
                name_id INTEGER NOT NULL,
                md5 VARCHAR(32) NOT NULL
            )
            """
        )
        con.commit()

    con.close()

    # Load jobs to send to workers
    with open(f"{matches_file}{_INDEX_SUFFIX}", "rb") as fh:
        index = pickle.load(fh)

    logger.info("populating signature2protein")
    inqueue = Queue()
    outqueue = Queue()
    workers = []
    for _ in range(max(1, processes - 1)):
        p = Process(target=_populate_signature2protein,
                    args=(url, tmp_db, matches_file, inqueue, outqueue))
        p.start()
        workers.append(p)

    total = 0
    for offset, count in index:
        inqueue.put((offset, count))
        total += count

    for _ in workers:
        inqueue.put(None)

    progress = 0
    for _ in index:
        count = outqueue.get()
        progress += count
        logger.info(f"{progress / total * 100:>7.0f}%")

    for p in workers:
        p.join()

    logger.info(f"temporary files: "
                f"{os.path.getsize(tmp_db) / 1024 ** 2:.0f} MB")
    os.unlink(tmp_db)

    con = psycopg2.connect(**pg.url2dict(url))
    with con.cursor() as cur:
        logger.info("indexing")
        cur.execute(
            """
            CREATE INDEX signature2protein_protein_idx
            ON signature2protein (protein_acc)
            """
        )
        cur.execute(
            """
            CREATE INDEX signature2protein_signature_idx
            ON signature2protein (signature_acc)
            """
        )
        con.commit()

    con.close()

    logger.info("done")


def _populate_signature2protein(url: str, names_db: str, matches_file: str,
                                src: Queue, dst: Queue):
    con = psycopg2.connect(**pg.url2dict(url))
    with con.cursor() as cur:
        proteins = _iter_proteins(names_db, matches_file, src, dst)
        cur.copy_from(file=pg.CsvIO(proteins, sep='|'),
                      table="signature2protein",
                      sep='|')
        con.commit()

    con.close()


def _iter_proteins(names_db: str, matches_file: str, src: Queue, dst: Queue):
    with KVdb(names_db) as names, open(matches_file, "rb") as fh:
        for offset, count in iter(src.get, None):
            fh.seek(offset)

            for _ in range(count):
                prot_acc, is_rev, is_comp, left_num, matches = pickle.load(fh)
                if not is_comp:
                    # Ignore fragmented proteins
                    continue

                name_id = names[prot_acc]
                md5 = _hash_matches(matches)

                for sig_acc in matches:
                    yield sig_acc, prot_acc, is_rev, left_num, name_id, md5

            dst.put(count)


def _hash_matches(matches: dict[str, tuple[str, bool, list[str]]]) -> str:
    # Flatten all matches
    locations = []

    for match_acc, (match_db, is_signature, hits) in matches.items():
        if is_signature:
            # Only consider signatures, not sequence features
            for start, end in _merge_overlapping(hits):
                locations.append((start, match_acc))
                locations.append((end, match_acc))

    """
    Evaluate the protein's match structure,
        i.e. how signatures match the proteins
    -----------------------------   Protein
     <    >                         Signature 1
       <    >                       Signature 2
                  < >               Signature 3
    Flattened:
    -----------------------------   Protein
     < <  > >     < >
     1 2  1 2     3 3
    Structure, with '-' representing a "gap"
        (more than N bp between two positions):
    1212-33
    """

    # Sort locations by position
    locations.sort()

    positions = iter(locations)

    """
    Do not set last_pos to 0, but to the first position:
    if two proteins have the same structure, but the distance between 0 
    and the first position is > max_gap in the first protein, 
    and <= max_gap in the second, a gap will be used for the first protein 
    and not for the other, which will result in two different hashes
    """
    last_pos, signature_acc = next(positions)
    # Overall domain organisation
    dom_org = []
    # Local domain organisation (positions within MAX_GAP)
    loc_org = [signature_acc]

    for pos, signature_acc in positions:
        if pos - last_pos > _MAX_GAP:
            dom_org += loc_org
            dom_org.append('')  # Add a gap
            loc_org.clear()

        loc_org.append(signature_acc)

    dom_org += loc_org

    return hashlib.md5('/'.join(dom_org).encode("utf-8")).hexdigest()


def _merge_overlapping(hits: list[str]) -> list[tuple[int, int]]:
    merged = []
    pos_start = pos_end = None
    for start, end in sorted(_flatten_hits(hits)):
        if pos_start is None:
            # Leftmost match
            pos_start = start
            pos_end = end
        elif start > pos_end:
            """
              pos_end
                ----] [----
                      start

            Gap: new location
            """
            merged.append((pos_start, pos_end))
            pos_start = start
            pos_end = end
        elif end > pos_end:
            """
                    pos_end
                ----]
                  ------]
                        end

            Extend current location
            """
            pos_end = end

    merged.append((pos_start, pos_end))
    return merged


def _flatten_hits(hits: list[str]) -> list[tuple[int, int]]:
    locations = []

    for fragments in hits:
        for fragment in fragments.split(","):
            for start, end, _ in fragment.split("-"):
                locations.append((int(start), int(end)))

    return locations


def insert_matches(uri: str, matches_file: str):
    logger.info("populating")
    con = psycopg2.connect(**pg.url2dict(uri))
    cur = con.cursor()
    cur.execute("DROP TABLE IF EXISTS match")
    cur.execute(
        """
        CREATE TABLE match (
            protein_acc VARCHAR(15) NOT NULL,
            signature_acc VARCHAR(25) NOT NULL,
            --database_id INTEGER NOT NULL,
            fragments TEXT NOT NULL
        )
        """
    )

    cur.execute("SELECT name, id FROM database")
    databases = dict(cur.fetchall())

    matches = _iter_matches(matches_file, databases)
    cur.copy_from(file=pg.CsvIO(matches, sep="|"),
                  table="match",
                  sep="|")
    con.commit()

    logger.info("indexing")
    cur.execute(
        """
        CREATE INDEX match_protein_idx
        ON interpro.match (protein_acc)
        """
    )
    con.commit()

    cur.close()
    con.close()

    logger.info("done")


def _iter_matches(matches_file: str, databases: dict[str, int]):
    for prot_acc, is_rev, is_comp, left_num, matches in _load(matches_file):
        for match_acc, (match_db, is_signature, hits) in matches.items():
            # sig_db_id = databases[sig_db]

            for fragments in hits:
                # yield prot_acc, sig_acc, sig_db_id, fragments
                yield prot_acc, match_acc, fragments
