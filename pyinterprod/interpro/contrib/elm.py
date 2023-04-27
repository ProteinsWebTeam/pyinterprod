import hashlib
from datetime import datetime

from .common import Method
from pyinterprod.utils.oracle import drop_table


def parse_instances(cur, signatures_source: str, sequences_source: str, last_modified: str) -> list[Method]:
    """
    Parse the ELM instances.tsv file, available at http://elm.eu.org/
    """
    valid_acc = get_updated_sequences(cur, sequences_source)

    instances = []
    match_data = []
    date = datetime.strptime(last_modified, "%d-%b-%Y-%H:%M:%S")
    with open(signatures_source, "rt") as fh:
        for line in fh:
            if line[0] == '#':
                continue
            else:
                elm_acc, _, elm_id, _, primary_acc, _, start, end, _, methods, inst_logic, _, _ = line.split("\t")
                if inst_logic.strip('"') == 'true positive':
                    if primary_acc in valid_acc:
                        instances.append(Method(elm_id, None, None, None, None, date))
                        match_data.append((primary_acc, elm_id, methods, int(start), int(end)))
    insert_matches(cur, match_data)
    return instances


def get_updated_sequences(cur, filepath: str) -> list[str]:
    cur.execute("""
        SELECT xr.AC, p.MD5
        FROM UNIPARC.XREF xr, UNIPARC.PROTEIN p
        WHERE xr.UPI = p.UPI 
            AND xr.DELETED = 'N' 
            AND xr.DBID IN (2, 3)
    """)
    sequences_md5 = dict(cur.fetchall())

    valid_acc = []
    with open(filepath, "rt") as fh:
        for line in fh:
            if line.startswith('>'):
                db, protein_acc, protein_id = line.split('|')
                fasta_sequence = next(fh)
                fasta_md5 = hashlib.md5(fasta_sequence.encode('utf-8')).hexdigest()

            if fasta_md5.upper() == sequences_md5[protein_acc].upper():
                valid_acc.append(protein_acc)
    return valid_acc


def insert_matches(cur, data: list):
    drop_table(cur, "INTERPRO.ELM_MATCH", purge=True)
    cur.execute(
        """
        CREATE TABLE INTERPRO.ELM_MATCH (
            PROTEIN_ID VARCHAR(15) NOT NULL,
            METHOD_AC VARCHAR2(25) NOT NULL,
            POS_FROM NUMBER(5) NOT NULL,
            POS_TO NUMBER(5) NOT NULL
        ) NOLOGGING
        """
    )

    for i in range(0, len(data)):
        cur.executemany(
            """
            INSERT INTO INTERPRO.ELM_MATCH (PROTEIN_ID, METHOD_AC, POS_FROM, POS_TO)
            VALUES (:1, :2, :3, :4)
            """,
            data[i]
        )

    cur.connection.commit()
