import math
import pickle
from multiprocessing import Process, Queue

import oracledb
import psycopg

from pyinterprod import logger
from pyinterprod.utils.oracle import clob_as_str
from pyinterprod.utils.pg import url2dict
from .match import INDEX_SUFFIX, merge_overlapping


"""
At least 50% of the shortest match must overlap with the other match
"""
_MIN_OVERLAP = 0.5

"""
One of the signatures must hit at least 50% of the proteins hit by the other
signature
"""
_MIN_COLLOCATION = 0.5

# Threshold for Jaccard index/coefficients
_MIN_SIMILARITY = 0.75


def _compare_signatures(matches_file: str, src: Queue, dst: Queue):
    signatures = {}
    comparisons = {}

    with open(matches_file, "rb") as fh:
        for offset, count in iter(src.get, None):
            fh.seek(offset)

            for _ in range(count):
                prot_acc, is_rev, is_comp, left_num, matches = pickle.load(fh)

                # Merge overlapping hits
                for signature_acc, models in matches.items():
                    for model_acc, (_, hits) in models.items():
                        hits = sorted(merge_overlapping(hits))
                        matches[signature_acc] = hits
                            
                for signature_acc in matches:
                    """
                    Count the number of proteins,
                    regardless of the sequence status
                    """
                    try:
                        sig = signatures[signature_acc]
                    except KeyError:
                        sig = signatures[signature_acc] = [
                            0,  # number of proteins
                            0,  # number of reviewed proteins
                            0,  # number of complete proteins
                            0,  # number of complete reviewed proteins
                            0,  # number of residues in complete proteins
                        ]
                        comparisons[signature_acc] = {}

                    sig[0] += 1
                    if is_rev:
                        sig[1] += 1

                    if not is_comp:
                        # Skip incomplete/fragment proteins
                        continue

                    sig[2] += 1
                    if is_rev:
                        sig[3] += 1

                    locs_1 = matches[signature_acc]

                    # Number of residues covered by the signature's matches
                    residues_1 = sum(end - start + 1 for start, end in locs_1)
                    sig[4] += residues_1

                    for other_acc in matches:
                        if other_acc <= signature_acc:
                            continue

                        locs_2 = matches[other_acc]
                        residues_2 = sum(
                            end - start + 1 for start, end in locs_2)

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

                        try:
                            cmp = comparisons[signature_acc][other_acc]
                        except KeyError:
                            cmp = comparisons[signature_acc][other_acc] = [
                                0,  # number of collocations
                                0,  # number of proteins with overlaps
                                0,  # number of overlapping residues
                            ]

                        cmp[0] += 1

                        # Overlapping proteins
                        shortest = min(residues_1, residues_2)
                        if residues >= _MIN_OVERLAP * shortest:
                            cmp[1] += 1

                        # Overlapping residues
                        cmp[2] += residues

            dst.put(count)

    dst.put((signatures, comparisons))


def insert_signatures(ora_uri: str, pg_uri: str, matches_file: str,
                      processes: int = 1):
    logger.info("iterating proteins")

    # Load jobs to send to workers
    with open(f"{matches_file}{INDEX_SUFFIX}", "rb") as fh:
        index = pickle.load(fh)

    inqueue = Queue()
    outqueue = Queue()
    workers = []
    for _ in range(max(1, processes - 1)):
        p = Process(target=_compare_signatures,
                    args=(matches_file, inqueue, outqueue))
        p.start()
        workers.append(p)

    total_proteins = 0
    for offset, count in index:
        inqueue.put((offset, count))
        total_proteins += count

    for _ in workers:
        inqueue.put(None)

    # Number of proteins (all, complete only, reviewed only) and residues
    signatures = {}
    # collocations/overlaps between signatures
    comparisons = {}

    done_proteins = done_workers = 0
    milestone = step = 10
    while done_workers < len(workers):
        obj = outqueue.get()
        if isinstance(obj, int):
            # Number of proteins processed by the worker
            done_proteins += obj
            progress = math.floor(done_proteins / total_proteins * 100)
            if progress >= milestone:
                logger.info(f"{progress}%")
                milestone += step
        else:
            # Aggregate results from worker
            done_workers += 1
            _signatures, _comparisons = obj

            for accession, counts in _signatures.items():
                if accession in signatures:
                    for i, count in enumerate(counts):
                        signatures[accession][i] += count
                else:
                    signatures[accession] = counts

            for accession, others in _comparisons.items():
                if accession not in comparisons:
                    comparisons[accession] = others
                    continue

                for other_acc, counts in others.items():
                    if other_acc in comparisons[accession]:
                        for i, count in enumerate(counts):
                            comparisons[accession][other_acc][i] += count
                    else:
                        comparisons[accession][other_acc] = counts

    for p in workers:
        p.join()

    # Load signatures from Oracle
    logger.info("loading signatures")
    con = oracledb.connect(ora_uri)
    cur = con.cursor()
    cur.outputtypehandler = clob_as_str
    cur.execute(
        """
        SELECT
            M.METHOD_AC, LOWER(D.DBSHORT), M.NAME, M.DESCRIPTION, T.ABBREV, 
            M.ABSTRACT, M.ABSTRACT_LONG
        FROM INTERPRO.METHOD M
        INNER JOIN INTERPRO.CV_DATABASE D ON M.DBCODE = D.DBCODE
        INNER JOIN INTERPRO.CV_ENTRY_TYPE T ON M.SIG_TYPE = T.CODE
        WHERE NOT REGEXP_LIKE(M.METHOD_AC, '^PTHR\d+:SF\d+$')
        """
    )
    rows = cur.fetchall()
    cur.close()
    con.close()

    con = psycopg.connect(**url2dict(pg_uri))
    with con.cursor() as cur:
        cur.execute("SELECT name, id FROM database")
        name2id = dict(cur.fetchall())

        logger.info("inserting signatures")
        values = []
        while rows:
            row = rows.pop()

            # Will raise a KeyError if the DB key is not in the database table
            signature_db_id = name2id[row[1]]

            signature_acc = row[0]
            values.append((
                signature_acc,
                signature_db_id,
                row[2],             # name
                row[3],             # description
                row[4],             # type
                row[5] or row[6],   # abstract
                *signatures.get(signature_acc, [0] * 5)
            ))

        cur.execute("DROP TABLE IF EXISTS signature")
        cur.execute(
            """
            CREATE TABLE signature (
                accession VARCHAR(25) NOT NULL 
                    CONSTRAINT signature_pkey PRIMARY KEY,
                database_id INTEGER NOT NULL,
                name VARCHAR(100),
                description VARCHAR(400),
                type VARCHAR(25) NOT NULL,
                abstract TEXT,
                num_sequences INTEGER NOT NULL,
                num_reviewed_sequences INTEGER NOT NULL,
                num_complete_sequences INTEGER NOT NULL,
                num_complete_reviewed_sequences INTEGER NOT NULL,
                num_residues BIGINT NOT NULL
            )
            """
        )

        sql = """
            COPY signature 
                (accession, database_id, name, description, type, abstract, 
                num_sequences, num_reviewed_sequences, num_complete_sequences, 
                num_complete_reviewed_sequences, num_residues) 
            FROM SDTIN
        """

        with cur.copy(sql) as copy:
            for rec in values:
                copy.write_row(rec)

        cur.execute(
            """
            CREATE INDEX signature_database_idx
            ON signature (database_id)
            """
        )
        con.commit()

        logger.info("inserting comparisons")
        cur.execute("DROP TABLE IF EXISTS comparison")
        cur.execute(
            """
            CREATE TABLE comparison (
                signature_acc_1 VARCHAR(25) NOT NULL,
                signature_acc_2 VARCHAR(25) NOT NULL,
                num_collocations INTEGER NOT NULL,
                num_overlaps INTEGER NOT NULL
            )
            """
        )

        sql = """
            COPY comparison 
                (signature_acc_1, signature_acc_2, num_collocations, 
                num_overlaps) 
            FROM STDIN
        """

        with cur.copy(sql) as copy:
            for row in _iter_comparisons(comparisons):
                copy.write_row(row)

        cur.execute(
            """
            CREATE UNIQUE INDEX comparison_idx
            ON comparison (signature_acc_1, signature_acc_2)
            """
        )
        con.commit()

        cur.execute("CLUSTER comparison USING comparison_idx")
        con.commit()

        logger.info("inserting predictions")
        cur.execute("DROP TABLE IF EXISTS prediction")
        cur.execute(
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

        sql = """
            COPY prediction 
                (signature_acc_1, signature_acc_2, num_collocations, 
                num_protein_overlaps, num_residue_overlaps) 
            FROM SDTIN
        """

        with cur.copy(sql) as copy:
            for row in _iter_predictions(comparisons, signatures):
                copy.write_row(row)

        cur.execute(
            """
            CREATE INDEX prediction_idx
            ON prediction (signature_acc_1)
            """
        )
        con.commit()

        cur.execute("CLUSTER prediction USING prediction_idx")
        con.commit()

    con.close()
    logger.info("done")


def _iter_comparisons(comps: dict[str, dict[str, list[int, int, int]]]):
    for acc1, others in comps.items():
        for acc2, (collocations, protein_overlaps, _) in others.items():
            yield acc1, acc2, collocations, protein_overlaps
            yield acc2, acc1, collocations, protein_overlaps


def _iter_predictions(comps: dict[str, dict[str, list[int, int, int]]],
                      sigs: dict[str, list[int, int, int, int, int]]):
    for acc1, others in comps.items():
        for acc2, (collocts, prot_overlaps, res_overlaps) in others.items():
            _, _, num_proteins1, num_reviewed1, num_residues1 = sigs[acc1]
            _, _, num_proteins2, num_reviewed2, num_residues2 = sigs[acc2]

            num_proteins = min(num_proteins1, num_proteins2)
            if collocts / num_proteins >= _MIN_COLLOCATION:
                yield acc1, acc2, collocts, prot_overlaps, res_overlaps
                yield acc2, acc1, collocts, prot_overlaps, res_overlaps


def get_swissprot_descriptions(pg_url: str) -> dict:
    con = psycopg.connect(**url2dict(pg_url))
    with con.cursor() as cur:
        cur.execute(
            """
            SELECT DISTINCT s2p.signature_acc, pn.text
            FROM interpro.signature2protein s2p
            INNER JOIN interpro.protein_name pn ON s2p.name_id = pn.name_id
            WHERE s2p.is_reviewed            
            """
        )

        signatures = {}
        for signature_acc, text in cur:
            try:
                signatures[signature_acc].add(text)
            except KeyError:
                signatures[signature_acc] = {text}

    con.close()

    return signatures
