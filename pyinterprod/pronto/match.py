# -*- coding: utf-8 -*-

import heapq
import os
import pickle
import shutil
from multiprocessing import Process, Queue
from tempfile import mkdtemp, mkstemp
from typing import Dict, List, Optional, Sequence

import cx_Oracle
import psycopg2

from pyinterprod import logger
from pyinterprod.utils.kvdb import KVdb
from pyinterprod.utils.pg import CsvIO, drop_index, url2dict


"""
At least 50% of the residues of the shortest signature 
  must overlap the other signature
(shorted signature = signature with the least residues in the protein)
"""
MIN_OVERLAP = 0.5

"""
One of the signatures must hit at least 50% of the proteins hit by the other
signature 
"""
MIN_COLLOCATION = 0.5

# Threshold for Jaccard index/coefficients
MIN_SIMILARITY = 0.75


def _dump_signatures(signatures: Dict[str, Sequence],
                     tmpdir: Optional[str]) -> str:
    fd, filepath = mkstemp(dir=tmpdir)
    os.close(fd)
    with open(filepath, "wb") as fh:
        for key in sorted(signatures):
            pickle.dump((key, set(signatures[key])), fh)

    return filepath


def _iter_matches(url: str, databases: Dict[str, int], outdir: str):
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.execute(
        """
        SELECT M.PROTEIN_AC, M.METHOD_AC, LOWER(D.DBSHORT),
               M.POS_FROM, M.POS_TO, M.FRAGMENTS
        FROM INTERPRO.MATCH M
        INNER JOIN INTERPRO.CV_DATABASE D ON M.DBCODE = D.DBCODE
        UNION ALL
        SELECT FM.PROTEIN_AC, FM.METHOD_AC, LOWER(D.DBSHORT),
               FM.POS_FROM, FM.POS_TO, NULL
        FROM INTERPRO.FEATURE_MATCH FM
        INNER JOIN INTERPRO.CV_DATABASE D ON FM.DBCODE = D.DBCODE
        WHERE FM.DBCODE = 'g'
        """
    )

    i = 0
    signatures = {}
    for row in cur:
        yield (
            row[0],
            row[1],
            databases[row[2]],
            # Format: start-end-type
            # with type=S (continuous single chain domain)
            row[5] if row[5] else f"{row[3]}-{row[4]}-S"
        )

        try:
            signatures[row[1]].append(row[0])
        except KeyError:
            signatures[row[1]] = [row[0]]

        i += 1
        if not i % 1000000:
            _dump_signatures(signatures, outdir)
            signatures = {}

        if not i % 100000000:
            logger.info(f"{i:>13,}")

    cur.close()
    con.close()

    _dump_signatures(signatures, outdir)
    logger.info(f"{i:>13,}")


def iterfile(path: str):
    with open(path, "rb") as fh:
        while True:
            try:
                obj = pickle.load(fh)
            except EOFError:
                break
            else:
                yield obj


def _agg_signatures(paths: Sequence[str]):
    files = [iterfile(path) for path in paths]
    signature = None
    proteins = set()
    for key, values in heapq.merge(*files, key=lambda x: x[0]):
        if key != signature:
            if signature:
                yield signature, len(proteins)

            signature = key
            proteins = set()

        proteins |= values

    if signature:
        yield signature, len(proteins)


def import_matches(ora_url: str, pg_url: str, output: str,
                   tmpdir: Optional[str] = None):
    tmpdir = mkdtemp(dir=tmpdir)

    logger.info("populating")
    pg_con = psycopg2.connect(**url2dict(pg_url))
    with pg_con.cursor() as pg_cur:
        pg_cur.execute("TRUNCATE TABLE match")
        pg_con.commit()

        drop_index(pg_con, "match_protein_idx")

        pg_cur.execute("SELECT name, id FROM database")
        databases = dict(pg_cur.fetchall())

        gen = _iter_matches(ora_url, databases, tmpdir)
        pg_cur.copy_from(file=CsvIO(gen, sep='|'), table="match", sep='|')

        logger.info("indexing")
        pg_cur.execute("ANALYZE match")
        pg_cur.execute(
            """
            CREATE INDEX match_protein_idx
            ON match (protein_acc)
            """
        )
        pg_con.commit()

    pg_con.close()

    logger.info("counting proteins per signature")
    paths = [os.path.join(tmpdir, name) for name in os.listdir(tmpdir)]
    with open(output, "wb") as fh:
        pickle.dump(dict(_agg_signatures(paths)), fh)

    size = 0
    for path in paths:
        size += os.path.getsize(path)
        os.remove(path)

    os.rmdir(tmpdir)
    logger.info(f"disk usage: {size/1024/1024:,.0f} MB")
    logger.info("complete")


def export_comp_seq_matches(url: str, filepath: str):
    logger.info("exporting matches")
    with open(filepath, "wb") as fh:
        i = 0
        for row in _iter_comp_seq_matches(url):
            pickle.dump(row, fh)
            i += 1
            if not i % 100000000:
                logger.info(f"{i:>13,}")

    logger.info(f"{i:>13,}")


def _iter_comp_seq_matches(url: str, filepath: Optional[str] = None):
    if filepath:
        return iterfile(filepath)
    else:
        con = cx_Oracle.connect(url)
        cur = con.cursor()
        cur.execute(
            """
            SELECT
                P.PROTEIN_AC, P.DBCODE, E.LEFT_NUMBER,
                M.METHOD_AC, M.POS_FROM, M.POS_TO, M.FRAGMENTS
            FROM INTERPRO.PROTEIN P
            INNER JOIN INTERPRO.ETAXI E
              ON P.TAX_ID = E.TAX_ID
              AND P.FRAGMENT = 'N'
            INNER JOIN INTERPRO.MATCH M
              ON P.PROTEIN_AC = M.PROTEIN_AC
            ORDER BY P.PROTEIN_AC
            """
        )

        protein_acc = None
        is_reviewed = None
        taxon_left_num = None
        matches = []
        for row in cur:
            if row[0] != protein_acc:
                if protein_acc:
                    yield protein_acc, is_reviewed, taxon_left_num, matches

                protein_acc = row[0]
                is_reviewed = row[1] == 'S'
                taxon_left_num = row[2]
                matches = []

            # Format: start-end-type
            # with type=S (continuous single chain domain)
            fragments = row[6] if row[6] else f"{row[4]}-{row[5]}-S"
            matches.append((row[3], fragments))

        cur.close()
        con.close()

        yield protein_acc, is_reviewed, taxon_left_num, matches


def _merge_matches(matches: Sequence[tuple]) -> Dict[str, List[tuple]]:
    # Merge matches by signature
    signatures = {}
    for signature_acc, fragments in matches:
        try:
            locations = signatures[signature_acc]
        except KeyError:
            locations = signatures[signature_acc] = []

        """
        Format: START-END-TYPE
        Types:
            * S: Continuous single chain domain
            * N: N terminus discontinuous
            * C: C terminus discontinuous
            * NC: N and C terminus discontinuous
        """
        for frag in fragments.split(','):
            pos_start, pos_end, pos_type = frag.split('-')
            locations.append((int(pos_start), int(pos_end)))

    # Merge overlapping locations
    for signature_acc in signatures:
        locations = []
        pos_start = pos_end = None
        for start, end in sorted(signatures[signature_acc]):
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
                locations.append((pos_start, pos_end))
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

        locations.append((pos_start, pos_end))
        signatures[signature_acc] = locations

    return signatures


def _process_chunk(url: str, names_db: str, inqueue: Queue, outqueue: Queue):
    signatures = {}  # number of proteins/residues per signature
    comparisons = {}  # collocations/overlaps between signatures

    con = psycopg2.connect(**url2dict(url))
    with con.cursor() as cur, KVdb(names_db) as names:
        for chunk in iter(inqueue.get, None):
            values = []
            for obj in chunk:
                protein_acc = obj[0]
                is_reviewed = obj[1]
                taxon_left_num = obj[2]
                matches = _merge_matches(obj[3])

                name_id = names[protein_acc]
                for signature_acc in matches:
                    values.append((
                        signature_acc,
                        protein_acc,
                        is_reviewed,
                        taxon_left_num,
                        name_id
                    ))

                    # Number of residues covered by signature's matches
                    locs_1 = matches[signature_acc]
                    residues_1 = sum(end - start + 1 for start, end in locs_1)

                    try:
                        obj = signatures[signature_acc]
                    except KeyError:
                        signatures[signature_acc] = [1, residues_1]
                        comparisons[signature_acc] = {}
                    else:
                        obj[0] += 1
                        obj[1] += residues_1

                    for other_acc in matches:
                        if other_acc <= signature_acc:
                            continue

                        locs_2 = matches[other_acc]
                        residues_2 = sum(end - start + 1
                                         for start, end in locs_2)

                        # Check overlapping matches
                        residues = 0
                        i = 0
                        start_2, end_2 = locs_2[i]
                        for start_1, end_1 in locs_1:
                            while end_2 < start_1:
                                i += 1
                                try:
                                    start_2, end_2 = locs_2[i]
                                except IndexError:
                                    break

                            # Overlap (start_1, end1) <-> (start_2, end_2)
                            o = min(end_1, end_2) - max(start_1, start_2) + 1
                            if o > 0:
                                residues += o  # matches overlap

                        # Add collocation, protein overlap, residue overlap
                        try:
                            cmp = comparisons[signature_acc][other_acc]
                        except KeyError:
                            comparisons[signature_acc][other_acc] = [0, 0, 0]
                            cmp = comparisons[signature_acc][other_acc]

                        # collocation
                        cmp[0] += 1

                        # Overlapping proteins
                        shortest = min(residues_1, residues_2)
                        if residues >= MIN_OVERLAP * shortest:
                            cmp[1] += 1

                        # Overlapping residues
                        cmp[2] += residues

            cur.copy_from(file=CsvIO(iter(values), sep='|'),
                          table="signature2protein",
                          sep='|')

    con.commit()
    con.close()
    outqueue.put((signatures, comparisons))


def proc_comp_seq_matches(ora_url: str, pg_url: str, database: str,
                          output: str, **kwargs):
    tmpdir = kwargs.get("tmpdir")
    matches_dat = kwargs.get("matches")
    processes = kwargs.get("processes", 4)

    logger.info("copying database")
    fd, tmp_database = mkstemp(dir=tmpdir)
    os.close(fd)
    shutil.copyfile(database, tmp_database)

    logger.info("populating: signature2protein")
    pg_con = psycopg2.connect(**url2dict(pg_url))
    with pg_con.cursor() as pg_cur:
        pg_cur.execute("TRUNCATE TABLE signature2protein")
        pg_con.commit()
        drop_index(pg_con, "signature2protein_protein_idx")
        drop_index(pg_con, "signature2protein_signature_idx")
    pg_con.close()

    inqueue = Queue(maxsize=1)
    outqueue = Queue()
    workers = []
    for _ in range(max(1, processes-1)):
        p = Process(target=_process_chunk,
                    args=(pg_url, tmp_database, inqueue, outqueue))
        p.start()
        workers.append(p)

    chunk = []
    i = 0
    for obj in _iter_comp_seq_matches(ora_url, matches_dat):
        chunk.append(obj)
        if len(chunk) == 1000:
            inqueue.put(chunk)
            chunk = []

        i += 1
        if not i % 10000000:
            logger.info(f"{i:>12,}")
    inqueue.put(chunk)

    for _ in workers:
        inqueue.put(None)

    logger.info(f"{i:>12,}")

    logger.info("aggregating signature comparisons")
    signatures = {}
    comparisons = {}
    for _ in workers:
        _signatures, _comparisons = outqueue.get()
        for signature_acc, [num_proteins, num_residues] in _signatures.items():
            try:
                obj = signatures[signature_acc]
            except KeyError:
                signatures[signature_acc] = [num_proteins, num_residues]
            else:
                obj[0] += num_proteins
                obj[1] += num_residues

        for signature_acc, others in _comparisons.items():
            try:
                cmp = comparisons[signature_acc]
            except KeyError:
                comparisons[signature_acc] = others
            else:
                for other_acc, values in others.items():
                    if other_acc in cmp:
                        for i, v in enumerate(values):
                            cmp[other_acc][i] += v
                    else:
                        cmp[other_acc] = list(values)

    for p in workers:
        p.join()

    with open(output, "wb") as fh:
        pickle.dump(signatures, fh)

    os.remove(tmp_database)

    pg_con = psycopg2.connect(**url2dict(pg_url))
    with pg_con.cursor() as pg_cur:
        logger.info("populating: comparison")
        pg_cur.execute("TRUNCATE TABLE comparison")
        pg_con.commit()
        drop_index(pg_con, "comparison_idx")
        gen = _iter_comparisons(comparisons)
        pg_cur.copy_from(file=CsvIO(gen, sep='|'),
                         table="comparison",
                         sep='|')
        pg_con.commit()

        logger.info("optimizing: comparison")
        pg_cur.execute("ANALYZE comparison")
        pg_cur.execute(
            """
            CREATE INDEX comparison_idx
            ON comparison (signature_acc_1, signature_acc_2)
            """
        )
        pg_cur.execute(
            """
            CLUSTER comparison 
            USING comparison_idx
            """
        )
        pg_con.commit()

        logger.info("populating: prediction")
        pg_cur.execute("TRUNCATE TABLE prediction")
        pg_con.commit()
        drop_index(pg_con, "prediction_idx")
        gen = _iter_predictions(signatures, comparisons)
        pg_cur.copy_from(file=CsvIO(gen, sep='|'),
                         table="prediction",
                         sep='|')
        pg_con.commit()

        logger.info("optimizing: prediction")
        pg_cur.execute("ANALYZE prediction")
        pg_cur.execute(
            """
            CREATE INDEX prediction_idx
            ON prediction (signature_acc_1)
            """
        )
        pg_cur.execute(
            """
            CLUSTER prediction 
            USING prediction_idx
            """
        )
        pg_con.commit()

        logger.info("optimizing: signature2protein")
        logger.debug("\tANALYZE")
        pg_cur.execute("ANALYZE signature2protein")
        logger.debug("\tsignature2protein_protein_idx")
        pg_cur.execute(
            """
            CREATE INDEX signature2protein_protein_idx
            ON signature2protein (protein_acc)
            """
        )
        logger.debug("\tsignature2protein_signature_idx")
        pg_cur.execute(
            """
            CREATE INDEX signature2protein_signature_idx
            ON signature2protein (signature_acc)
            """
        )
        logger.debug("\tCLUSTER")
        pg_cur.execute(
            """
            CLUSTER signature2protein 
            USING signature2protein_signature_idx
            """
        )
        pg_con.commit()

    pg_con.close()
    logger.info("complete")


def _iter_comparisons(comparisons: dict):
    for acc1, others in comparisons.items():
        for acc2, [collocs, prot_overlaps, res_overlaps] in others.items():
            yield acc1, acc2, collocs, prot_overlaps
            yield acc2, acc1, collocs, prot_overlaps


def _iter_predictions(signatures: dict, comparisons: dict):
    for acc1, others in comparisons.items():
        for acc2, [collocs, prot_overlaps, res_overlaps] in others.items():
            num_proteins1, num_residues_1 = signatures[acc1]
            num_proteins2, num_residues_2 = signatures[acc2]

            if collocs / min(num_proteins1, num_proteins2) >= MIN_COLLOCATION:
                yield acc1, acc2, collocs, prot_overlaps, res_overlaps
                yield acc2, acc1, collocs, prot_overlaps, res_overlaps
