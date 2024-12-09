import gzip
import json
import os
import re
import tarfile

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
    for entry, _ in parse_sto(pfam_seed_file):
        accession, version = entry.features["AC"].split(".")
        comment = entry.features.get("CC")
        rn2pmid = {}
        if comment:
            for i, obj in enumerate(entry.features.get("RN", [])):
                rn2pmid[i+1] = int(" ".join(obj["RM"]))

            abstract = _repl_references(accession, comment, rn2pmid)

            abstract = _repl_xreferences(abstract)
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
        for ref_num in _expand_range(match.group(1)):
            try:
                pmid = references[ref_num]
            except KeyError:
                logger.warning(f"{acc}: missing PubMed ID "
                               f"for reference #{ref_num}")
            else:
                refs.append(f"PMID:{pmid}")

        return f"[{','.join(refs)}]"

    return re.sub(r"\[([\d\s,-]+)]", _repl, text)


def _expand_range(s: str) -> list[int]:
    r = []
    for i in s.split(','):
        values = i.split("-")
        if len(values) == 1:
            # single reference number
            r.append(int(i))
        else:
            # range of reference numbers, e.g. 1-4
            low, high = map(int, values)
            r += list(range(low, high + 1))

    return r


def _repl_xreferences(text: str) -> str:
    # EC number
    text = re.sub(r"EC:(\d+\.\d+\.\d+\.(\d+|-))", r"[ec:\1]", text)

    # Pfam
    text = re.sub(r"Pfam:(PF\d{5})", r"[pfam:\1]", text)

    # Swiss-Prot
    text = re.sub(r"Swiss:([A-Z0-9]+)", r"[swissprot:\1]", text)

    return text


class StockholdMSA:
    def __init__(self):
        self.features = {}
        self.sequences = []
        self.complete = False

    def process_line(self, line: str):
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
        raw = ""
        for line in map(_decode, fh):
            raw += line + "\n"
            entry.process_line(line)
            if entry.complete:
                yield entry, raw
                entry = StockholdMSA()
                raw = ""

    assert len(raw) == 0


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
    for entry, raw in parse_sto(pfama_seed):
        seeds[entry.features["AC"]] = (entry, raw)

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
    progress = 0
    for full_entry, full_raw in parse_sto(pfama_full):
        seed_entry, seed_raw = seeds.pop(full_entry.features["AC"])

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
                gzip.compress(seed_raw.encode("utf-8"), compresslevel=6),
                full_entry.features["SQ"],
                gzip.compress(full_raw.encode("utf-8"), compresslevel=6),
                json.dumps(authors),
                json.dumps(full_entry.features.get("WK", []))
            ]
        )
        progress += 1
        if progress % 100 == 0:
            logger.info(f"{progress:>15,}")

    logger.info(f"{progress:>15,}")
    con.commit()
    cur.close()
    con.close()

    logger.info("done")


def get_clans(pfam_c: str, pfama_full: str) -> list[Clan]:
    """Extract information about Pfam clans and their members

    :param pfam_c: string representation of the path to Pfam-C[.gz]
    :param pfama_full: string representation of the path to Pfam-A.full[.gz]
    :return: list of Clan objects
    """
    logger.info(f"parsing {os.path.basename(pfama_full)}")
    num_full = {}
    for entry, _ in parse_sto(pfama_full):
        accession, version = entry.features["AC"].split(".")
        num_full[accession] = entry.features["SQ"]

    logger.info(f"parsing {os.path.basename(pfam_c)}")
    clans = []
    for entry, _ in parse_sto(pfam_c):
        accession, version = entry.features["AC"].split(".")
        name = entry.features["ID"]
        description = entry.features["DE"]

        total = 0
        members = []
        for member in entry.features.get("MB", []):
            member = member.rstrip(";")
            num_seqs = num_full.get(member, 0)
            if num_seqs > total:
                total = num_seqs
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
    for entry, _ in parse_sto(pfam_c):
        accession, version = entry.features["AC"].split(".")
        name = entry.features["ID"]
        description = entry.features["DE"]

        authors = []
        for author in entry.features["AU"]:
            authors += [e.strip() for e in author.split(",")]

        rn2pmid = {}
        references = []
        for i, ref_dict in enumerate(entry.features.get("RN", [])):
            pmid = int(" ".join(ref_dict["RM"]))
            rn2pmid[i+1] = pmid
            references.append({
                "PMID": pmid,
                "title": " ".join(ref_dict["RT"]),
                "authors": list(
                    map(
                        str.strip,
                        " ".join(ref_dict["RA"]).rstrip(";").split(",")
                    )
                ),
                "journal": " ".join(ref_dict["RL"])
            })

        comment = entry.features.get("CC")
        if comment:
            comment = _repl_references(accession, comment, rn2pmid)
            comment = _repl_xreferences(comment)

        cur.execute(
            """
            INSERT /*+ APPEND */  INTO INTERPRO.PFAM_C
            VALUES (:1, :2, :3, :4, :5, :6)
            """,
            [
                accession,
                name,
                description,
                comment,
                json.dumps(authors),
                json.dumps(references)
            ]
        )

    con.commit()
    cur.close()
    con.close()
