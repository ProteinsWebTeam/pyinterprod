import csv
import hashlib

from pyinterprod.utils.oracle import drop_table

from .common import Method


def ignore_comments(file: str):
    with open(file, "rt", newline='') as fh:
        for line in fh:
            if line and line[0] != "#":
                yield line


def load_classes(cur, classes_file: str, instances_file: str,
                 sequences_file: str) -> list[Method]:
    methods = {}

    # Load classes from TSV
    content = ignore_comments(classes_file)
    reader = csv.DictReader(content, quotechar='"', delimiter="\t")
    for row in reader:
        m = Method(accession=row["accession"],
                   sig_type=None,
                   name=row["ELMIdentifier"],
                   description=row["Description"])
        methods[m.name] = m

    # Load the MD5 hash of sequences used in ELM
    checksums = {}
    with open(sequences_file, "rt") as fh:
        for line in fh:
            if line and line[0] == ">":
                source, accession, identifier = line[1:].split("|")
                # Some accessions are suffixed with what I guess is the version
                accession = accession.split("-")[0]
                sequence = next(fh).strip().encode("utf-8")
                checksums[accession] = hashlib.md5(sequence).hexdigest()

    # Check whether the sequences have changed
    sequences = list(checksums.keys())
    sequences_ok = set()
    step = 1000
    for i in range(0, len(sequences), step):
        params = sequences[i:i + step]
        args = ",".join([":" + str(i + 1) for i in range(len(params))])

        cur.execute(
            f"""
            SELECT X.AC, P.MD5
            FROM UNIPARC.XREF X
            INNER JOIN UNIPARC.PROTEIN P ON X.UPI = P.UPI
            WHERE X.AC IN ({args})
              AND X.DBID IN (2, 3)             
              AND X.DELETED = 'N'  
            """,
            params
        )

        for accession, md5 in cur.fetchall():
            if md5 == checksums[accession]:
                sequences_ok.add(accession)

    # Load instances (i.e. protein matches)
    matches = []
    content = ignore_comments(instances_file)
    reader = csv.DictReader(content, quotechar='"', delimiter="\t")
    for row in reader:
        elm_name = row["ELMIdentifier"]
        uniprot_acc = row["ProteinName"]
        pos_from = int(row["Start"])
        pos_to = int(row["End"])
        evidences = row["Methods"]
        logic = row["InstanceLogic"]

        if elm_name not in methods:
            continue
        elif uniprot_acc not in sequences_ok:
            continue
        elif logic != "true positive":
            continue

        matches.append((
            uniprot_acc,
            methods[elm_name].accession,
            # evidences,
            pos_from,
            pos_to
        ))

    # Populate ELM_MATCH
    drop_table(cur, "INTERPRO.ELM_MATCH", purge=True)
    cur.execute(
        """
        CREATE TABLE INTERPRO.ELM_MATCH (
            PROTEIN_ID VARCHAR(15) NOT NULL,
            METHOD_AC VARCHAR2(25) NOT NULL,
            -- EVIDENCES VARCHAR2(500),
            POS_FROM NUMBER(5) NOT NULL,
            POS_TO NUMBER(5) NOT NULL
        ) NOLOGGING
        """
    )
    cur.executemany(
        """
        INSERT INTO INTERPRO.ELM_MATCH
        VALUES (:1, :2, :3, :4)
        """,
        matches
    )
    cur.connection.commit()

    return list(methods.values())
