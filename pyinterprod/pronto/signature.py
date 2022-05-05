import cx_Oracle
import psycopg2
from psycopg2.extras import execute_values

from pyinterprod import logger
from pyinterprod.utils.oracle import clob_as_str
from pyinterprod.utils.pg import url2dict, CsvIO
from .match import iter_util_eof, merge_overlapping


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


def insert_signatures(ora_uri: str, pg_uri: str, matches_file: str):
    logger.info("iterating proteins")

    # Number of proteins (all, complete only, reviewed only) and residues
    signatures = {}
    # collocations/overlaps between signatures
    comparisons = {}
    progress = 0

    for prot_acc, is_rev, is_comp, _, matches in iter_util_eof(matches_file):
        # Merge overlapping hits
        for signature_acc, (_, hits) in matches.items():
            matches[signature_acc] = sorted(merge_overlapping(hits))

        for signature_acc in matches:
            # Count the number of proteins, regardless of the sequence status
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
                residues_2 = sum(end - start + 1 for start, end in locs_2)

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

        progress += 1
        if progress % 1e7 == 0:
            logger.info(f"{progress:>15,}")

    logger.info(f"{progress:>15,}")

    # Load signatures from Oracle
    logger.info("loading signatures")
    con = cx_Oracle.connect(ora_uri)
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
        """
    )
    rows = cur.fetchall()
    cur.close()
    con.close()

    con = psycopg2.connect(**url2dict(pg_uri))
    with con.cursor() as cur:
        cur.execute("SELECT name, id FROM database")
        name2id = dict(cur.fetchall())

        logger.info("inserting signatures")
        values = []
        while rows:
            row = rows.pop()
            signature_acc = row[0]

            values.append((
                signature_acc,
                name2id[row[1]],
                row[2],
                row[3],
                row[4],
                row[5] or row[6],
                *signatures.get(signature_acc, [0] * 5)
            ))

        cur.execute("DROP TABLE IF EXISTS signature")
        cur.execute(
            """
            CREATE TABLE signature (
                accession VARCHAR(25) NOT NULL 
                    CONSTRAINT signature_pkey PRIMARY KEY,
                database_id INTEGER NOT NULL,
                name VARCHAR(100) NOT NULL,
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
        execute_values(cur, "INSERT INTO signature VALUES %s", values,
                       page_size=1000)
        con.commit()

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

        cur.copy_from(file=CsvIO(_iter_comparisons(comparisons), sep='|'),
                      table="comparison",
                      sep='|')
        con.commit()

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

        records = _iter_predictions(comparisons, signatures)
        cur.copy_from(file=CsvIO(records, sep='|'),
                      table="prediction",
                      sep='|')
        con.commit()

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
    con = psycopg2.connect(**url2dict(pg_url))
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
