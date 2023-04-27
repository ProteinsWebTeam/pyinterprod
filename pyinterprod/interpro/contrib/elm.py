import hashlib
import logging
from datetime import datetime

from .common import Method
from pyinterprod.utils.oracle import drop_table


def parse_instances(cur, signatures_source: str, sequences_source: str) -> list[Method]:
    """
    Parse the ELM instances.tsv file, available at http://elm.eu.org/
    """
    valid_acc = get_updated_sequences(cur, sequences_source)

    instances = set()
    match_data = set()
    with open(signatures_source, "rt") as fh:
        for line in fh:
            if line[0] == '#':
                continue
            else:
                elm_acc, _, elm_id, _, primary_acc, _, start, end, _, methods, inst_logic, _, _ = line.split("\t")
                if inst_logic.strip('"') == 'true positive':
                    if primary_acc in valid_acc:
                        instances.add(Method(elm_id, None))
                        match_data.add((primary_acc, elm_id, methods, int(start), int(end)))
    insert_matches(cur, list(match_data))
    return list(instances)


def get_updated_sequences(cur, filepath: str) -> list[str]:
    fasta_acc = []
    fasta_md5 = {}
    with open(filepath, "rt") as fh:
        for line in fh:
            if line.startswith('>'):
                _, protein_acc, _ = line.split('|')
                protein_acc = protein_acc.split("-")[0]
                fasta_sequence = next(fh)
                fasta_md5[protein_acc] = hashlib.md5(fasta_sequence.encode('utf-8')).hexdigest()
                fasta_acc.append(protein_acc)

    valid_acc = []
    step = 1000
    for i in range(0, len(fasta_acc), step):
        params = fasta_acc[i:i+step]
        args = ",".join([":" + str(i+1) for i in range(len(params))])
        cur.execute(f"""
            SELECT xr.AC, p.MD5
            FROM UNIPARC.XREF xr
            JOIN UNIPARC.PROTEIN p ON xr.UPI = p.UPI
            WHERE xr.DELETED = 'N' 
                AND xr.DBID IN (2, 3)
                AND xr.AC IN ({args})
        """, params)
        protein_md5 = dict(cur.fetchall())

        for acc in fasta_acc[i:i+step]:
            try:
                if fasta_md5[acc].upper() == protein_md5[acc].upper():
                    valid_acc.append(protein_acc)
            except KeyError:
                logging.warning(f"Can't find a valid protein accession '{acc}'")

    return valid_acc


def insert_matches(cur, data: list):
    drop_table(cur, "INTERPRO.ELM_MATCH", purge=True)
    cur.execute(
        """
        CREATE TABLE INTERPRO.ELM_MATCH (
            PROTEIN_ID VARCHAR(15) NOT NULL,
            METHOD_AC VARCHAR2(25) NOT NULL,
            SEQ_FEATURE VARCHAR2(60) NOT NULL,
            POS_FROM NUMBER(5) NOT NULL,
            POS_TO NUMBER(5) NOT NULL
        ) NOLOGGING
        """
    )

    for i in range(0, len(data)):
        cur.executemany(
            """
            INSERT INTO INTERPRO.ELM_MATCH (PROTEIN_ID, METHOD_AC, SEQ_FEATURE, POS_FROM, POS_TO)
            VALUES (:1, :2, :3, :4, :5)
            """,
            data[i]
        )

    cur.connection.commit()
