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


def get_signatures(tsvfile: str,
                   cur: oracledb.Cursor | None = None) -> list[Method]:
    """Parse NCBIFAM TSV file
    https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.tsv
    and returns signatures to be imported in InterPro
    :param tsvfile: Path to TSV file
    :param cur: Oracle cursor or None
    :return:
    """
    methods = []
    amr_models = []
    for values in parse_tsv(tsvfile):
        model_acc = values[0]
        signature_acc, model_version = model_acc.split(".")

        name = values[2]
        _type = values[6]
        for_AMRFinder = values[9] == "Y"

        if values[15]:
            references = list(map(int, set(values[15].split(","))))
        else:
            references = []

        source = values[20]
        hmm_name = values[22]
        comment = values[23].strip()

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

        if for_AMRFinder:
            amr_models.append((signature_acc,))

    if cur:
        drop_table(cur, "INTERPRO.NCBIFAM_AMR", purge=True)
        cur.execute(
            """
            CREATE TABLE INTERPRO.NCBIFAM_AMR
            (
                METHOD_AC VARCHAR2(25) NOT NULL
                CONSTRAINT PK_NCBIFAM_AMR PRIMARY KEY
            )
            """
        )
        cur.execute("GRANT SELECT ON INTERPRO.NCBIFAM_AMR TO INTERPRO_SELECT")
        if amr_models:
            cur.executemany("INSERT INTO INTERPRO.NCBIFAM_AMR VALUES (:1)",
                            amr_models)
            cur.connection.commit()

    return methods


def create_custom_hmm(tsvfile: str, hmmfile: str):
    models = {m.model: m for m in get_signatures(tsvfile)}

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

        if accession in signatures and values[14]:
            for go_id in set(values[14].split(",")):
                records.append((accession, go_id))

    for i in range(0, len(records), 1000):
        cur.executemany("INSERT INTO INTERPRO.NCBIFAM2GO VALUES (:1, :2)",
                        records[i:i+1000])

    con.commit()
    cur.close()
    con.close()


def main():
    parser = ArgumentParser(description="create a custom NCBIFAM "
                                        "HMM file for InterProScan")
    parser.add_argument("tsv", metavar="hmm_PGAP.tsv",
                        help="signatures information, in TSV format")
    parser.add_argument("hmm", metavar="hmm_PGAP.LIB",
                        help="file of concatenated HMM profiles")
    args = parser.parse_args()
    create_custom_hmm(args.tsv, args.hmm)


if __name__ == "__main__":
    main()
