import gzip
import json
import os
import re
import tarfile
from io import BytesIO

import oracledb

from .common import Clan, Method
from pyinterprod import logger
from pyinterprod.utils import Table
from pyinterprod.utils.oracle import drop_table


_TYPES = {
    "Domain": 'D',
    "Family": 'F',
    "Repeat": 'R',
    "Coiled-coil": 'I',
    "Disordered": 'O',
    "Motif": 'C'
}
_STO_MAIN_FIELDS = {"AC", "ID", "DE", "AU", "BM", "SM", "GA", "TC", "TP", "SQ"}
_STO_OTHER_FIELDS = {"DR", "CC", "WK", "CL", "MB"}
_STO_REFERENCE_FIELDS = {"RC", "RM", "RT", "RA", "RL"}


def iter_interpro_n_matches(file: str):
    """Iterate Pfam inferences from the InterPro-N archive
    """
    with tarfile.open(file, mode="r") as tar:
        for member in tar:
            if member.name.endswith(".tsv"):
                br = tar.extractfile(member)
                lines = br.read().decode("utf-8").splitlines(keepends=False)

                for line in lines[1:]:
                    uniprot_acc, pfam_acc, start, end, score = line.split("\t")
                    yield uniprot_acc, pfam_acc, int(start), int(end)


def get_pfam_n_entries(cur: oracledb.Cursor, file: str) -> list[Method]:
    drop_table(cur, "INTERPRO.PFAMN_MATCH", purge=True)
    cur.execute(
        """
        CREATE TABLE INTERPRO.PFAMN_MATCH (
            PROTEIN_ID VARCHAR(15) NOT NULL,
            METHOD_AC VARCHAR2(25) NOT NULL,
            POS_FROM NUMBER(5) NOT NULL,
            POS_TO NUMBER(5) NOT NULL
        ) NOLOGGING
        """
    )

    query = """
        INSERT /*+ APPEND */ 
        INTO INTERPRO.PFAMN_MATCH 
        VALUES (:1, :2, :3, :4)
    """

    with Table(con=cur.connection, query=query, autocommit=True) as table:
        accessions = set()
        for uniprot_acc, pfam_acc, start, end in iter_interpro_n_matches(file):
            accessions.add(pfam_acc)
            table.insert((uniprot_acc, pfam_acc, start, end))

    return [Method(pfam_acc, sig_type=None) for pfam_acc in accessions]


def get_signatures(pfam_seed_file: str) -> list[Method]:
    methods = []
    for entry in parse_sto(pfam_seed_file):
        accession, version = entry.features["AC"].split(".")
        comment = entry.features.get("CC")
        rn2pmid = {}
        if comment:
            for i, obj in enumerate(entry.features.get("RN", [])):
                rn2pmid[i+1] = int(" ".join(obj["RM"]))

            abstract = _repl_references(accession, comment, rn2pmid)
        else:
            abstract = None

        methods.append(
            Method(
                accession,
                _TYPES[entry.features["TP"]],
                entry.features["ID"],
                entry.features["DE"],
                abstract,
                list(rn2pmid.values())
            )
        )

    return methods


def _repl_references(acc: str, text: str, references: dict[int, int]):
    def _repl(match: re.Match) -> str:
        refs = []
        for ref_num in map(int, map(str.strip, match.group(1).split(','))):
            try:
                pmid = references[ref_num]
            except KeyError:
                logger.warning(f"{acc}: missing PubMed ID "
                               f"for reference #{ref_num}")
            else:
                refs.append(f"PMID:{pmid}")

        return f"[{','.join(refs)}]"

    return re.sub(r"\[([\d\s,]+)]", _repl, text)


class StockholdMSA:
    def __init__(self):
        self.lines = ""
        self.features = {}
        self.sequences = []
        self.complete = False

    def add_line(self, line: str):
        self.lines += line
        if line == "//":
            self.close()
        elif line.startswith("#=GF"):
            # Generic per-File
            _, field, value = line.split(maxsplit=2)
            self.add_feature(field, value)
        elif line.startswith("#=GS"):
            # Generic per-Sequence
            pass
        elif line.startswith("#=GR"):
            # Generic per-Residue
            pass
        elif line.startswith("#=GC"):
            # Generic per-Column
            pass
        elif line[0] != "#":
            # Sequence line
            self.sequences.append(line)

    def add_feature(self, name: str, value: str):
        if name == "RN":
            try:
                self.features[name].append({})
            except KeyError:
                self.features[name] = [{}]
        elif name in _STO_REFERENCE_FIELDS:
            if "RN" not in self.features:
                # logger.warning(f"{name} {value} ignored (preceding RN field)")
                return

            ref_dict = self.features["RN"][-1]
            try:
                ref_dict[name].append(value)
            except KeyError:
                ref_dict[name] = [value]
        elif name in _STO_MAIN_FIELDS or name in _STO_OTHER_FIELDS:
            try:
                self.features[name].append(value)
            except KeyError:
                self.features[name] = [value]

    def close(self):
        if self.complete:
            return

        self.complete = True
        for key, values in self.features.items():
            if key in ("BM", "SM", "CC", "DE"):
                self.features[key] = " ".join(values)
            elif key in ("AU", "WK", "RN", "DR", "MB"):
                pass
            elif len(values) != 1:
                # Other fields should have one value only
                raise ValueError(f"more than one value "
                                 f"for field {key}: {values}")
            elif key == "SQ":
                num_sequences = int(values.pop())
                if num_sequences == len(self.sequences):
                    self.features[key] = num_sequences
                else:
                    raise ValueError(f"inconsistent number of sequences: "
                                     f"{num_sequences} != {len(self.sequences)}"
                                     f" ({self.features})")
            else:
                self.features[key] = values.pop()

    def compress(self, compresslevel: int = 6) -> bytes:
        with BytesIO() as bs:
            with gzip.GzipFile(mode="wb",
                               compresslevel=compresslevel,
                               fileobj=bs) as gz:
                for line in self.lines:
                    gz.write((line + "\n").encode("utf-8"))

            data = bs.getvalue()
        return data

    def is_empty(self) -> bool:
        return len(self.lines) == 0


def parse_sto(file: str):
    """Parse a Gzip-compressed Pfam file in the Stockhold format

    :param file: string representing the path to the file.
    """
    if file.endswith(".gz"):
        _open = gzip.open
    else:
        _open = open

    with _open(file, "rb") as fh:
        entry = StockholdMSA()
        for line in map(_decode, fh):
            entry.add_line(line)
            if entry.complete:
                yield entry
                entry = StockholdMSA()

    assert entry.is_empty() is True


def parse_fasta(file: str) -> dict[str, int]:
    """Count the number of sequences per Pfam family

    :param file: string representation of the path to Pfam-A.fasta.gz
    :return: dictionary of Pfam accession -> number of sequences
    """
    if file.endswith(".gz"):
        _open = gzip.open
    else:
        _open = open

    counts = {}
    with _open(file, "rt") as fh:
        for line in map(str.rstrip, fh):
            if line[0] == ">":
                # >A0A8S0GYS1_9PSED/45-323 A0A8S0GYS1.1 PF00389.35;2-Hacid_dh;
                _, _, entry = line.split(None, maxsplit=2)
                accession, version = entry.split(".")
                try:
                    counts[accession] += 1
                except KeyError:
                    counts[accession] = 1

    return counts


def _decode(b: bytes) -> str:
    """Decode bytes to a string using UTF-8 or Latin-1

    :param b: bytes to decode
    :return: a string, with any trailing whitespace trimmed
    """
    try:
        s = b.decode("utf-8")
    except UnicodeDecodeError:
        pass
    else:
        return s.rstrip()

    return b.decode("latin-1").rstrip()


def persist_pfam_a(uri: str, pfama_seed: str, pfama_full: str):
    """Extract information about Pfam families and persist it in Oracle.

    :param uri: Oracle connection string
    :param pfama_seed: string representation of the path to Pfam-A.seed[.gz]
    :param pfama_full: string representation of the path to Pfam-A.full[.gz]
    """
    logger.info(f"parsing {os.path.basename(pfama_seed)}")
    seeds = {}
    for entry in parse_sto(pfama_seed):
        seeds[entry.features["AC"]] = entry

    con = oracledb.connect(uri)
    cur = con.cursor()
    drop_table(cur, "INTERPRO.PFAM_A", purge=True)
    cur.execute(
        """
        CREATE TABLE INTERPRO.PFAM_A (
            ACCESSION VARCHAR2(7),
            VERSION NUMBER NOT NULL,
            SEQ_ONTOLOGY_ID VARCHAR2(20),
            BUILD_CMD VARCHAR2(255) NOT NULL,
            SEARCH_CMD VARCHAR2(255) NOT NULL,
            SEQ_GA FLOAT NOT NULL,
            DOM_GA FLOAT NOT NULL,
            SEED_NUM NUMBER NOT NULL,
            SEED_ALN BLOB NOT NULL,
            FULL_NUM NUMBER NOT NULL,
            FULL_ALN BLOB NOT NULL,
            AUTHORS VARCHAR2(1000) NOT NULL,
            WIKIPEDIA VARCHAR2(255) NOT NULL,
            CONSTRAINT PK_PFAM_A PRIMARY KEY (ACCESSION)
        ) NOLOGGING
        """
    )

    logger.info(f"parsing {os.path.basename(pfama_full)}")
    for full_entry in parse_sto(pfama_full):
        seed_entry = seeds.pop(full_entry.features["AC"])

        accession, version = full_entry.features["AC"].split(".")
        seq_ontology = None
        for ref in full_entry.features.get("DR", []):
            m = re.match(r"SO;\s*(\d+);", ref)
            if m:
                seq_ontology = m.group(1)
                break

        seq_ga_str, dom_ga_str = full_entry.features["GA"].rstrip(";").split()
        seq_ga = float(seq_ga_str)
        dom_ga = float(dom_ga_str)
        authors = []
        for author in full_entry.features["AU"]:
            name, orcid = author.split(";")
            authors.append({
                "author": author,
                "orcid": orcid or None
            })

        cur.execute(
            """
            INSERT /*+ APPEND */  INTO INTERPRO.PFAM_A
            VALUES (:1, :2, :3, :4, :5, :6, :7, :8, :9, :10, :11, :12, :13)
            """,
            [
                accession,
                version,
                seq_ontology,
                full_entry.features["BM"],
                full_entry.features["SM"],
                seq_ga,
                dom_ga,
                seed_entry.features["SQ"],
                seed_entry.compress(),
                full_entry.features["SQ"],
                full_entry.compress(),
                json.dumps(authors),
                json.dumps(full_entry.features.get("WK", []))
            ]
        )

    con.commit()
    cur.close()
    con.close()

    logger.info("done")


def get_clans(pfam_c: str, pfama_fa: str, pfama_full: str) -> list[Clan]:
    """Extract information about Pfam clans and their members

    :param pfam_c: string representation of the path to Pfam-C[.gz]
    :param pfama_fa: string representation of the path to Pfam-A.fasta[.gz]
    :param pfama_full: string representation of the path to Pfam-A.full[.gz]
    :return: list of Clan objects
    """
    logger.info(f"parsing {os.path.basename(pfama_fa)}")
    num_seqs = parse_fasta(pfama_fa)

    logger.info(f"parsing {os.path.basename(pfama_full)}")
    num_full = {}
    for entry in parse_sto(pfama_full):
        accession, version = entry.features["AC"].split(".")
        num_full[accession] = entry.features["SQ"]

    logger.info(f"parsing {os.path.basename(pfam_c)}")
    clans = []
    for entry in parse_sto(pfam_c):
        accession, version = entry.features["AC"].split(".")
        name = entry.features["ID"]
        description = entry.features["DE"]

        total = 0
        members = []
        for member in entry.features.get("MB", []):
            member = member.rstrip(";")
            total += num_seqs.get(member, 0)
            members.append(member)

        clans.append(Clan(accession, name, description))
        for member in members:
            clans[-1].members.append({
                "accession": member,
                "score": num_full.get(member, 0) / total if total > 0 else 0
            })

    return clans


def persist_pfam_c(uri: str, pfam_c: str):
    """Extract Clan information from Pfam-C and persist it in Oracle.

    :param uri: Oracle connection string
    :param pfam_c: string representation of the path to Pfam-C[.gz]
    """
    con = oracledb.connect(uri)
    cur = con.cursor()
    drop_table(cur, "INTERPRO.PFAM_C", purge=True)
    cur.execute(
        """
        CREATE TABLE INTERPRO.PFAM_C (
            ACCESSION VARCHAR2(6),
            NAME VARCHAR2(40) NOT NULL,
            DESCRIPTION VARCHAR2(100) NOT NULL,
            ABSTRACT VARCHAR2(4000),
            AUTHORS VARCHAR2(255) NOT NULL,
            REFERENCES VARCHAR2(4000) NOT NULL,
            CONSTRAINT PK_PFAM_C PRIMARY KEY (ACCESSION)
        ) NOLOGGING
        """
    )

    logger.info(f"parsing {os.path.basename(pfam_c)}")
    for entry in parse_sto(pfam_c):
        accession, version = entry.features["AC"].split(".")
        name = entry.features["ID"]
        description = entry.features["DE"]

        authors = []
        for author in entry.features["AU"]:
            authors += [e.strip() for e in author.split(",")]

        references = []
        for ref_dict in entry.features.get("RN", []):
            references.append({
                "PMID": int(" ".join(ref_dict["RM"])),
                "title": " ".join(ref_dict["RT"]),
                "authors": list(
                    map(
                        str.strip,
                        " ".join(ref_dict["RA"]).rstrip(";").split(",")
                    )
                ),
                "journal": " ".join(ref_dict["RL"])
            })

        cur.execute(
            """
            INSERT /*+ APPEND */  INTO INTERPRO.PFAM_C
            VALUES (:1, :2, :3, :4, :5, :6)
            """,
            [
                accession,
                name,
                description,
                entry.features.get("CC"),
                json.dumps(authors),
                json.dumps(references)
            ]
        )

    con.commit()
    cur.close()
    con.close()
