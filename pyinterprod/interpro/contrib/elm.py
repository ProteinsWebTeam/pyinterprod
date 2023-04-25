import hashlib
from .common import Method


def parse_instances(cur, signatures_source: str, sequences_source: str) -> list[Method]:
    """
    Parse the ELM instances.tsv file, available at http://elm.eu.org/
    """
    valid_acc = get_updated_sequences(cur, sequences_source)

    instances = []
    with open(signatures_source, "rt") as fh:
        date = None
        for line in fh:
            if line[0] == '#':
                continue
            else:
                cols = line.split("\t")
                if cols[10].strip('"') == 'true positive':

                    instances.append(Method(cols[0].strip('"'), None, cols[9].strip('"'), None, None, date))

    return instances


def get_updated_sequences(cur, filepath: str) -> list[str]:
    cur.execute("""
        SELECT xr.AC, p.MD5
        FROM UNIPARC.XREF xr, UNIPARC.PROTEIN p
        WHERE xr.UPI = p.UPI 
            AND xr.DELETED = 'N' 
            AND xr.DBID IN (2, 3);
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
