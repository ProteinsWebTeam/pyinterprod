import gzip
import hashlib
import heapq
import math
import os
import pickle
import shutil
from multiprocessing import Process, Queue
from tempfile import mkstemp

import oracledb
import psycopg

from pyinterprod import logger
from pyinterprod.utils import pg
from pyinterprod.utils.kvdb import KVdb


# Domain org.: introduce a gap when distance between two positions > 20 aa
_MAX_GAP = 20

INDEX_SUFFIX = ".i"


def export(url: str, output: str, cachesize: int = 10000000,
           tmpdir: str | None = None):
    if tmpdir:
        os.makedirs(tmpdir, exist_ok=True)

    logger.info("exporting matches")
    files = _export_matches(url, cachesize, tmpdir)

    logger.info("sorting proteins")
    index = []
    with open(output, "wb") as fh:
        i = o = 0
        for protein in _merge_matches(files):
            pickle.dump(protein, fh)

            i += 1
            if i % cachesize == 0:
                index.append((o, cachesize))
                o = fh.tell()

            if i % 1e8 == 0:
                logger.info(f"{i:>15,}")

        if i % cachesize != 0:
            index.append((o, i % cachesize))

        logger.info(f"{i:>15,}")

    size = 0
    for file in files:
        size += os.path.getsize(file)
        os.unlink(file)

    logger.info(f"temporary files: {size / 1024 ** 2:.0f} MB")

    # Write index
    with open(f"{output}{INDEX_SUFFIX}", "wb") as fh:
        pickle.dump(index, fh)

    logger.info("done")


def _export_matches(url: str, cachesize: int,
                    tmpdir: str | None = None) -> list[str]:
    con = oracledb.connect(url)
    cur = con.cursor()

    # Loading databases
    cur.execute(
        """
        SELECT DBCODE, LOWER(DBSHORT)
        FROM INTERPRO.CV_DATABASE
        """
    )
    databases = dict(cur.fetchall())

    # Loading taxonomy
    cur.execute(
        """
        SELECT TAX_ID, LEFT_NUMBER
        FROM INTERPRO.ETAXI
        """
    )
    taxonomy = dict(cur.fetchall())

    cur.execute(
        """
        SELECT P.PROTEIN_AC, P.DBCODE, P.FRAGMENT, P.TAX_ID, M.METHOD_AC, 
            M.DBCODE, M.FRAGMENTS, M.POS_FROM, M.POS_TO, M.MODEL_AC
        FROM INTERPRO.MATCH M
        INNER JOIN INTERPRO.PROTEIN P 
            ON P.PROTEIN_AC = M.PROTEIN_AC
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
                taxonomy[row[3]],  # taxon left number
                []  # matches
            )

        if row[6]:
            fragments = row[6]
        else:
            # Single/continuous fragment using match start/end positions
            fragments = f"{row[7]}-{row[8]}-S"

        obj[3].append((
            row[4],  # match accession
            databases[row[5]],  # match DB
            fragments,
            row[9] if row[5] == 'V' else row[4]  # used for PANTHER subfamilies
        ))

        i += 1
        if i % cachesize == 0:
            files.append(_dump(cache, tmpdir))

        if i % 1e8 == 0:
            logger.info(f"{i:>15,}")

    cur.close()
    con.close()

    files.append(_dump(cache, tmpdir))
    logger.info(f"{i:>15,}")

    return files


def _dump(data: dict, tmpdir: str | None = None) -> str:
    fd, file = mkstemp(dir=tmpdir)
    os.close(fd)
    with gzip.open(file, "wb", compresslevel=6) as fh:
        for key in sorted(data):
            pickle.dump((key, data[key]), fh)

    data.clear()
    return file


def _merge_matches(files: list[str]):
    iterable = [iter_util_eof(file, compressed=True) for file in files]
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

        for signature_acc, signature_db, fragments, model_acc in value[3]:
            try:
                models = matches[signature_acc]
            except KeyError:
                models = matches[signature_acc] = {}

            try:
                models[model_acc][1].append(fragments)
            except KeyError:
                models[model_acc] = (signature_db, [fragments])

    yield protein_acc, is_reviewed, is_complete, left_number, matches


def iter_util_eof(file: str, compressed: bool):
    if compressed:
        fh = gzip.open(file, "rb")
    else:
        fh = open(file, "rb")

    try:
        while True:
            try:
                yield pickle.load(fh)
            except EOFError:
                break
    finally:
        fh.close()


def insert_signature2protein(url: str, names_db: str, matches_file: str,
                             processes: int = 1, tmpdir: str | None = None):
    if tmpdir:
        os.makedirs(tmpdir, exist_ok=True)

    logger.info("copying SQLite database")
    fd, tmp_db = mkstemp(dir=tmpdir)
    os.close(fd)
    shutil.copyfile(names_db, tmp_db)

    logger.info("creating signature2protein")
    con = psycopg.connect(**pg.url2dict(url))
    with con.cursor() as cur:
        cur.execute("DROP TABLE IF EXISTS signature2protein")
        cur.execute(
            """
            CREATE TABLE signature2protein (
                signature_acc VARCHAR(25) NOT NULL,
                model_acc VARCHAR(25),
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
    with open(f"{matches_file}{INDEX_SUFFIX}", "rb") as fh:
        index = pickle.load(fh)

    logger.info("populating")
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

    done = 0
    milestone = step = 10
    for _ in index:
        count = outqueue.get()
        done += count
        progress = math.floor(done / total * 100)
        if progress >= milestone:
            logger.info(f"{progress}%")
            milestone += step

    for p in workers:
        p.join()

    logger.info(f"temporary files: "
                f"{os.path.getsize(tmp_db) / 1024 ** 2:.0f} MB")
    os.unlink(tmp_db)

    logger.info("done")


def _populate_signature2protein(url: str, names_db: str, matches_file: str,
                                src: Queue, dst: Queue):
    con = psycopg.connect(**pg.url2dict(url))
    with con.cursor() as cur:
        proteins = _iter_proteins(names_db, matches_file, src, dst)

        sql = """
            INSERT INTO signature2protein 
                (signature_acc, model_acc, protein_acc, is_reviewed, 
                taxon_left_num, name_id, md5) 
            VALUES (%s, %s, %s, %s, %s, %s, %s)
        """

        records = []
        for row in proteins:
            records.append(row)
            if len(records) == 1000:
                cur.executemany(sql, records)
                con.commit()
                records.clear()
        if records:
            cur.executemany(sql, records)
            con.commit()
            records.clear()

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
                    for model_acc in matches[sig_acc]:
                        if matches[sig_acc][model_acc][0] == 'panther':
                            yield (sig_acc, model_acc, prot_acc, is_rev, 
                                   left_num, name_id, md5)
                        else:
                            yield (sig_acc, None, prot_acc, is_rev, 
                                   left_num, name_id, md5)

            dst.put(count)


def _hash_matches(matches: dict[str, tuple[str, list[str]]]) -> str:
    # Flatten all matches
    locations = []

    for signature_acc in matches:
        for model_acc, (_, hits) in matches[signature_acc].items():
            for start, end in merge_overlapping(hits):
                locations.append((start, signature_acc))
                locations.append((end, signature_acc))

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


def merge_overlapping(hits: list[str]) -> list[tuple[int, int]]:
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
            start, end, _ = fragment.split("-")
            locations.append((int(start), int(end)))

    return locations


def finalize_signature2protein(uri: str):
    con = psycopg.connect(**pg.url2dict(uri))
    with con.cursor() as cur:
        logger.info("indexing")

        logger.info("\tprotein_acc")
        cur.execute(
            """
            CREATE INDEX IF NOT EXISTS signature2protein_protein_idx
            ON signature2protein (protein_acc)
            """
        )
        con.commit()

        logger.info("\tsignature_acc")
        cur.execute(
            """
            CREATE INDEX IF NOT EXISTS signature2protein_signature_idx
            ON signature2protein (signature_acc)
            """
        )
        con.commit()

        logger.info("\tclustering")
        cur.execute("CLUSTER signature2protein "
                    "USING signature2protein_signature_idx")
        con.commit()

    con.close()
    logger.info("done")


def create_match_table(uri: str):
    con = psycopg.connect(**pg.url2dict(uri))
    with con.cursor() as cur:
        cur.execute("DROP TABLE IF EXISTS match")
        cur.execute(
            """
            CREATE TABLE match (
                protein_acc VARCHAR(15) NOT NULL,
                signature_acc VARCHAR(30) NOT NULL,
                database_id INTEGER NOT NULL,
                fragments TEXT NOT NULL
            )
            """
        )
        con.commit()

    con.close()
    logger.info("done")


def insert_fmatches(ora_uri: str, pg_uri: str):
    logger.info("populating")

    con = psycopg.connect(**pg.url2dict(pg_uri))
    with con.cursor() as cur:
        cur.execute("SELECT name, id FROM database")
        name2id = dict(cur.fetchall())

        sql = """
            INSERT INTO match 
                (protein_acc, signature_acc, database_id, fragments) 
            VALUES (%s, %s, %s, %s)
        """

        records = []
        for row in _get_fmatches(ora_uri, name2id):
            records.append(row)
            if len(records) == 1000:
                cur.executemany(sql, records)
                con.commit()
                records.clear()
        if records:
            cur.executemany(sql, records)
            con.commit()
            records.clear()

    con.close()
    logger.info("done")


def _get_fmatches(uri: str, name2id: dict[str, int]):
    con = oracledb.connect(uri)
    cur = con.cursor()
    cur.execute(
        """
        SELECT FM.PROTEIN_AC, FM.METHOD_AC, LOWER(D.DBSHORT), POS_FROM, POS_TO
        FROM INTERPRO.FEATURE_MATCH FM
        INNER JOIN INTERPRO.CV_DATABASE D ON FM.DBCODE = D.DBCODE
        WHERE FM.DBCODE IN ('a', 'f')
        """
    )

    for protein_acc, feature_acc, database, start, end in cur:
        database_id = name2id[database]

        # No fragment column in FEATURE_MATCH, so pretend we have one
        fragments = f"{start}-{end}-S"

        yield protein_acc, feature_acc, database_id, fragments

    cur.close()
    con.close()


def insert_matches(uri: str, matches_file: str, processes: int = 1):
    # Load jobs to send to workers
    with open(f"{matches_file}{INDEX_SUFFIX}", "rb") as fh:
        index = pickle.load(fh)

    logger.info("populating")
    inqueue = Queue()
    outqueue = Queue()
    workers = []
    for _ in range(max(1, processes - 1)):
        p = Process(target=_populate_matches,
                    args=(uri, matches_file, inqueue, outqueue))
        p.start()
        workers.append(p)

    total = 0
    for offset, count in index:
        inqueue.put((offset, count))
        total += count

    for _ in workers:
        inqueue.put(None)

    done = 0
    milestone = step = 10
    for _ in index:
        count = outqueue.get()
        done += count
        progress = math.floor(done / total * 100)
        if progress >= milestone:
            logger.info(f"{progress}%")
            milestone += step

    for p in workers:
        p.join()

    logger.info("done")


def _populate_matches(url: str, matches_file: str, src: Queue, dst: Queue):
    con = psycopg.connect(**pg.url2dict(url))
    with con.cursor() as cur:
        cur.execute("SELECT name, id FROM database")
        name2id = dict(cur.fetchall())
        matches = _iter_matches(matches_file, name2id, src, dst)

        sql = """
              INSERT INTO match (protein_acc, signature_acc, database_id, fragments) 
              VALUES (%s, %s, %s, %s)
              """

        records = []
        for row in matches:
            records.append(row)
            if len(records) == 1000:
                cur.executemany(sql, records)
                con.commit()
                records.clear()
        if records:
            cur.executemany(sql, records)
            con.commit()
            records.clear()

    con.close()


def _iter_matches(matches_file: str, name2id: dict[str, int],
                  src: Queue, dst: Queue):
    with open(matches_file, "rb") as fh:
        for offset, count in iter(src.get, None):
            fh.seek(offset)

            for _ in range(count):
                prot_acc, is_rev, is_comp, left_num, matches = pickle.load(fh)

                for sig_acc in matches:
                    for model_acc, (sig_db, hits) in matches[sig_acc].items():
                        sig_db_id = name2id[sig_db]

                        for fragments in hits:
                            yield prot_acc, sig_acc, sig_db_id, fragments

            dst.put(count)


def finalize_match_table(uri: str):
    con = psycopg.connect(**pg.url2dict(uri))
    with con.cursor() as cur:

        logger.info("indexing")
        cur.execute(
            """
            CREATE INDEX IF NOT EXISTS match_protein_idx
            ON interpro.match (protein_acc)
            """
        )
        con.commit()

        logger.info("clustering")
        cur.execute("CLUSTER match USING match_protein_idx")
        con.commit()

    con.close()
    logger.info("done")
