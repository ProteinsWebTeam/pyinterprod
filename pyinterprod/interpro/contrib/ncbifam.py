from argparse import ArgumentParser
from configparser import ConfigParser
import re

import oracledb

from pyinterprod.utils.oracle import drop_table
from .common import Method


def parse_tsv(tsvfile: str):
    with open(tsvfile, "rt") as fh:
        for line in fh:
            if line[0] != "#":
                yield line.split("\t")


def get_signatures(cur: oracledb.Cursor,
                   tsvfile: str,
                   triagefile: str,
                   keep_integrated: bool = True) -> list[Method]:
    """Parse NCBIfam TSV file
    https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.tsv
    and returns signatures to be imported in InterPro
    :param cur: Oracle cursor
    :param tsvfile: Path to TSV file
    :param triagefile: Path to TXT file of families to import
    :param keep_integrated: If True, keep integrated signatures
           even if they are not in the list of models to import
    :return:
    """
    cur.execute("SELECT METHOD_AC FROM INTERPRO.ENTRY2METHOD")
    integrated = {signature_acc for signature_acc, in cur.fetchall()}

    to_import = set()
    with open(triagefile, "rt") as fh:
        for model_acc in map(str.rstrip, fh):
            to_import.add(model_acc)

    methods = []
    for values in parse_tsv(tsvfile):
        model_acc = values[0]
        signature_acc, model_version = model_acc.split(".")

        if model_acc not in to_import:
            if signature_acc not in integrated or not keep_integrated:
                continue

        name = values[2]
        _type = values[6]
        # for_AMRFinder = values[9] == "Y"

        if values[14]:
            references = list(map(int, set(values[14].split(","))))
        else:
            references = []

        source = values[19]
        hmm_name = values[21]
        comment = values[22].strip()

        if not comment or comment.lower() == "null":
            comment = None

        if _type == "repeat":
            _type = "R"
        elif _type == "domain":
            _type = "D"
        else:
            # Defaults to family
            _type = "F"

        m = Method(
            accession=signature_acc,
            sig_type=_type,
            name=name,
            description=hmm_name,
            abstract=comment,
            references=references,
            model=model_acc
        )
        methods.append(m)

    return methods


def create_custom_hmm(cur: oracledb.Cursor, tsvfile: str, triagefile: str,
                      hmmfile: str):
    models = {m.model: m for m in get_signatures(cur, tsvfile, triagefile)}

    reg_acc = re.compile(r"^ACC\s+(\w+\.\d+)$", flags=re.M)
    reg_desc = re.compile(r"^(DESC\s+)(.+)$", flags=re.M)

    with open(hmmfile, "rt") as fh:
        buffer = ""
        for line in fh:
            buffer += line

            if line[:2] == "//":
                model_acc = reg_acc.search(buffer).group(1)
                if model_acc in models:
                    desc = models[model_acc].description
                    buffer = reg_desc.sub(replace_desc(desc), buffer)
                    print(buffer, end="")

                buffer = ""

        if buffer:
            model_acc = reg_acc.search(buffer).group(1)
            if model_acc in models:
                desc = models[model_acc].description
                buffer = reg_desc.sub(replace_desc(desc), buffer)
                print(buffer, end="")


def replace_desc(new_desc: str):
    def repl(match: re.Match):
        return f"{match.group(1)}{new_desc}"

    return repl


def update_go_terms(uri: str, tsvfile: str):
    con = oracledb.connect(uri)
    cur = con.cursor()
    cur.execute("SELECT METHOD_AC FROM INTERPRO.METHOD WHERE DBCODE = 'N'")
    signatures = {acc for acc, in cur.fetchall()}

    drop_table(cur, "INTERPRO.NCBIFAM2GO", purge=True)
    cur.execute(
        """
        CREATE TABLE INTERPRO.NCBIFAM2GO
        (
            METHOD_AC VARCHAR2(25) NOT NULL
                CONSTRAINT FK_NCBIFAM2GO
                REFERENCES INTERPRO.METHOD (METHOD_AC) ON DELETE CASCADE,
            GO_ID VARCHAR2(10) NOT NULL,
            CONSTRAINT PK_NCBIFAM2GO
            PRIMARY KEY (METHOD_AC, GO_ID)
        ) NOLOGGING
        """
    )

    records = []
    for values in parse_tsv(tsvfile):
        accession = values[0].split(".")[0]

        if accession in signatures and values[13]:
            for go_id in set(values[13].split(",")):
                records.append((accession, go_id))

    for i in range(0, len(records), 1000):
        cur.executemany("INSERT INTO INTERPRO.NCBIFAM2GO VALUES (:1, :2)",
                        records[i:i+1000])

    con.commit()
    cur.close()
    con.close()


def main():
    parser = ArgumentParser(description="create a custom NCBIfam "
                                        "HMM file for InterProScan")
    parser.add_argument("config", metavar="pyinterprod.conf",
                        help="configuration file")
    parser.add_argument("tsv", metavar="hmm_PGAP.tsv",
                        help="signatures information, in TSV format")
    parser.add_argument("triage", metavar="for-interpro.tsv",
                        help="triage file of models to import")
    parser.add_argument("hmm", metavar="hmm_PGAP.LIB",
                        help="file of concatenated HMM profiles")
    args = parser.parse_args()

    config = ConfigParser()
    config.read(args.config)
    uri = config["oracle"]["ipro-interpro"]

    con = oracledb.connect(uri)
    cur = con.cursor()

    try:
        create_custom_hmm(cur, args.tsv, args.triage, args.hmm)
    finally:
        cur.close()
        con.close()


if __name__ == "__main__":
    main()
