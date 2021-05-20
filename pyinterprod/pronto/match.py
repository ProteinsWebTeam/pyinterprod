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
from pyinterprod.utils.pg import CsvIO, url2dict


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


class MatchIterator:
    def __init__(self, url: str, databases: Dict[str, int],
                 tmpdir: Optional[str] = None):
        self.url = url
        self.databases = databases
        self.tmpdir = tmpdir
        self.files = []

    def __iter__(self):
        con = cx_Oracle.connect(self.url)
        cur = con.cursor()

        # Get accession for Swiss-Prot proteins
        cur.execute(
            """
            SELECT PROTEIN_AC 
            FROM INTERPRO.PROTEIN 
            WHERE DBCODE = 'S'
            """
        )
        reviewed = {acc for acc, in cur.fetchall()}

        cur.execute(
            """
            SELECT M.PROTEIN_AC, M.METHOD_AC, LOWER(D.DBSHORT),
                   M.POS_FROM, M.POS_TO, M.FRAGMENTS
            FROM INTERPRO.MATCH M
            INNER JOIN INTERPRO.CV_DATABASE D ON M.DBCODE = D.DBCODE
            --UNION ALL
            --SELECT FM.PROTEIN_AC, FM.METHOD_AC, LOWER(D.DBSHORT),
            --       FM.POS_FROM, FM.POS_TO, NULL
            --FROM INTERPRO.FEATURE_MATCH FM
            --INNER JOIN INTERPRO.CV_DATABASE D ON FM.DBCODE = D.DBCODE
            --WHERE FM.DBCODE = 'g'
            """
        )

        i = 0
        signatures = {}
        for row in cur:
            prot_acc = row[0]
            sign_acc = row[1]

            yield (
                prot_acc,
                sign_acc,
                self.databases[row[2]],
                # Format: start-end-type
                # with type=S (continuous single chain domain)
                row[5] if row[5] else f"{row[3]}-{row[4]}-S"
            )

            try:
                matches = signatures[sign_acc]
            except KeyError:
                matches = signatures[sign_acc] = []

            matches.append((prot_acc, prot_acc in reviewed))

            i += 1
            if not i % 1000000:
                self.dump(signatures)
                signatures = {}

            if not i % 100000000:
                logger.info(f"{i:>13,}")

        cur.close()
        con.close()

        self.dump(signatures)
        logger.info(f"{i:>13,}")

    def dump(self, signatures: Dict[str, Sequence]):
        fd, path = mkstemp(dir=self.tmpdir)
        with open(fd, "wb") as fh:
            for sign_acc in sorted(signatures):
                reviewed_proteins = set()
                n_reviewed_matches = 0
                unreviewed_proteins = set()

                for prot_acc, is_reviewed in signatures[sign_acc]:
                    if is_reviewed:
                        reviewed_proteins.add(prot_acc)
                        n_reviewed_matches += 1
                    else:
                        unreviewed_proteins.add(prot_acc)

                pickle.dump((
                    sign_acc,
                    reviewed_proteins,
                    n_reviewed_matches,
                    unreviewed_proteins
                ), fh)

        self.files.append(path)

    @staticmethod
    def load(path: str):
        with open(path, "rb") as fh:
            while True:
                try:
                    obj = pickle.load(fh)
                except EOFError:
                    break
                else:
                    yield obj

    def merge(self) -> Dict[str, tuple]:
        counts = {}
        signature = None
        reviewed_proteins = set()
        unreviewed_proteins = set()
        n_reviewed_matches = 0

        iterable = [self.load(path) for path in self.files]
        for item in heapq.merge(*iterable, key=lambda x: x[0]):
            sig_acc, rev_prots, n_rev_matches, unrev_prots = item

            if sig_acc != signature:
                if signature:
                    counts[signature] = (
                        len(reviewed_proteins),
                        n_reviewed_matches,
                        len(unreviewed_proteins)
                    )

                signature = sig_acc
                reviewed_proteins = set()
                unreviewed_proteins = set()
                n_reviewed_matches = 0

            reviewed_proteins |= rev_prots
            n_reviewed_matches += n_rev_matches
            unreviewed_proteins |= unrev_prots

        if signature:
            counts[signature] = (
                len(reviewed_proteins),
                n_reviewed_matches,
                len(unreviewed_proteins)
            )

        return counts

    @property
    def size(self) -> int:
        return sum(map(os.path.getsize, self.files))

    def remove(self):
        for path in self.files:
            os.remove(path)


def import_matches(ora_url: str, pg_url: str, output: str,
                   tmpdir: Optional[str] = None):
    logger.info("populating")
    pg_con = psycopg2.connect(**url2dict(pg_url))
    pg_cur = pg_con.cursor()
    pg_cur.execute("DROP TABLE IF EXISTS match")
    pg_cur.execute(
        """
        CREATE TABLE match (
            protein_acc VARCHAR(15) NOT NULL,
            signature_acc VARCHAR(25) NOT NULL,
            database_id INTEGER NOT NULL,
            fragments TEXT NOT NULL
        )
        """
    )

    pg_cur.execute("SELECT name, id FROM database")
    databases = dict(pg_cur.fetchall())

    matches = MatchIterator(ora_url, databases, tmpdir)
    pg_cur.copy_from(file=CsvIO(iter(matches), sep='|'),
                     table="match",
                     sep='|')

    logger.info("indexing")
    pg_cur.execute(
        """
        CREATE INDEX match_protein_idx
        ON match (protein_acc)
        """
    )
    pg_con.commit()
    pg_cur.close()
    pg_con.close()

    logger.info("counting proteins per signature")
    with open(output, "wb") as fh:
        pickle.dump(matches.merge(), fh)

    logger.info(f"disk usage: {matches.size/1024**2:,.0f} MB")
    matches.remove()
    logger.info("complete")


def export_comp_seq_matches(url: str, filepath: str):
    logger.info("exporting matches")
    with open(filepath, "wb") as fh:
        i = 0
        for row in iter_comp_seq_matches(url):
            pickle.dump(row, fh)
            i += 1
            if not i % 100000000:
                logger.info(f"{i:>13,}")

    logger.info(f"{i:>13,}")


def iter_comp_seq_matches(url: str, filepath: Optional[str] = None):
    if filepath:
        return MatchIterator.load(filepath)
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
            INNER JOIN INTERPRO.MATCH M
              ON P.PROTEIN_AC = M.PROTEIN_AC
            WHERE P.FRAGMENT = 'N'
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


def merge_matches(matches: Sequence[tuple]) -> Dict[str, List[tuple]]:
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


def process_chunk(url: str, names_db: str, inqueue: Queue, outqueue: Queue):
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
                matches = merge_matches(obj[3])

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
        pg_cur.execute("DROP TABLE IF EXISTS signature2protein")
        pg_cur.execute(
            """
            CREATE TABLE signature2protein (
                signature_acc VARCHAR(25) NOT NULL,
                protein_acc VARCHAR(15) NOT NULL,
                is_reviewed BOOLEAN NOT NULL,
                taxon_left_num INTEGER NOT NULL,
                name_id INTEGER NOT NULL
            )
            """
        )
        pg_con.commit()

    pg_con.close()

    inqueue = Queue(maxsize=1)
    outqueue = Queue()
    workers = []
    for _ in range(max(1, processes-1)):
        p = Process(target=process_chunk,
                    args=(pg_url, tmp_database, inqueue, outqueue))
        p.start()
        workers.append(p)

    chunk = []
    i = 0
    for obj in iter_comp_seq_matches(ora_url, matches_dat):
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
        pg_cur.execute("DROP TABLE IF EXISTS comparison")
        pg_cur.execute(
            """
            CREATE TABLE comparison (
                signature_acc_1 VARCHAR(25) NOT NULL,
                signature_acc_2 VARCHAR(25) NOT NULL,
                num_collocations INTEGER NOT NULL,
                num_overlaps INTEGER NOT NULL
            )
            """
        )

        gen = iter_comparisons(comparisons)
        pg_cur.copy_from(file=CsvIO(gen, sep='|'),
                         table="comparison",
                         sep='|')

        logger.info("optimizing: comparison")
        pg_cur.execute(
            """
            CREATE UNIQUE INDEX comparison_idx
            ON comparison (signature_acc_1, signature_acc_2)
            """
        )
        pg_con.commit()

        pg_cur.execute("CLUSTER comparison USING comparison_idx")
        pg_con.commit()

        logger.info("populating: prediction")
        pg_cur.execute("DROP TABLE IF EXISTS prediction")
        pg_cur.execute(
            """
            CREATE TABLE prediction (
                signature_acc_1 VARCHAR(25) NOT NULL,
                signature_acc_2 VARCHAR(25) NOT NULL,
                num_collocations INTEGER NOT NULL,
                num_protein_overlaps INTEGER NOT NULL,
                num_residue_overlaps BIGINT NOT NULL
            )
            """
        )

        gen = iter_predictions(signatures, comparisons)
        pg_cur.copy_from(file=CsvIO(gen, sep='|'),
                         table="prediction",
                         sep='|')

        logger.info("optimizing: prediction")
        pg_cur.execute(
            """
            CREATE INDEX prediction_idx
            ON prediction (signature_acc_1)
            """
        )
        pg_con.commit()
        pg_cur.execute("CLUSTER prediction USING prediction_idx")
        pg_con.commit()

        logger.info("optimizing: signature2protein")
        logger.info("\tsignature2protein_protein_idx")
        pg_cur.execute(
            """
            CREATE INDEX signature2protein_protein_idx
            ON signature2protein (protein_acc)
            """
        )
        logger.info("\tsignature2protein_signature_idx")
        pg_cur.execute(
            """
            CREATE INDEX signature2protein_signature_idx
            ON signature2protein (signature_acc)
            """
        )
        pg_con.commit()

        logger.info("\tCLUSTER")
        pg_cur.execute(
            """
            CLUSTER signature2protein
            USING signature2protein_signature_idx
            """
        )
        pg_con.commit()

    pg_con.close()
    logger.info("complete")


def iter_comparisons(comparisons: dict):
    for acc1, others in comparisons.items():
        for acc2, [collocs, prot_overlaps, res_overlaps] in others.items():
            yield acc1, acc2, collocs, prot_overlaps
            yield acc2, acc1, collocs, prot_overlaps


def iter_predictions(signatures: dict, comparisons: dict):
    for acc1, others in comparisons.items():
        for acc2, [collocs, prot_overlaps, res_overlaps] in others.items():
            num_proteins1, num_residues_1 = signatures[acc1]
            num_proteins2, num_residues_2 = signatures[acc2]

            if collocs / min(num_proteins1, num_proteins2) >= MIN_COLLOCATION:
                yield acc1, acc2, collocs, prot_overlaps, res_overlaps
                yield acc2, acc1, collocs, prot_overlaps, res_overlaps
