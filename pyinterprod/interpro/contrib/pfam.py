import gzip
import os
import re

import MySQLdb
import MySQLdb.cursors
import oracledb

from .common import Clan, Method
from pyinterprod import logger
from pyinterprod.utils import Table
from pyinterprod.utils.oracle import drop_table
from pyinterprod.utils.pg import url2dict


_RECORD_BREAK = "//"
_TYPES = {
    "Domain": 'D',
    "Family": 'F',
    "Repeat": 'R',
    "Coiled-coil": 'I',
    "Disordered": 'O',
    "Motif": 'C'
}


class AbstractFormatter:
    def __init__(self,
        references: dict[int, int],
    ):
        self.refs = references
        self.acc = None

    def update(self, acc: str, abstract: str) -> str:
        self.acc = acc
        # abstract = re.sub(r"\[(\d+-\d+)]", self.split_refs, abstract)
        abstract = re.sub(r"\[([\d\s,]+)]", self.repl_refs, abstract)
        return abstract

    @staticmethod
    def split_refs(match: re.Match) -> str:
        start, end = map(int, match.group(1).split("-"))
        refs = []
        for ref_pos in range(start, end + 1):
            refs.append(str(ref_pos))

        return f"[{', '.join(refs)}]"

    def repl_refs(self, match: re.Match) -> str:
        refs = []
        for ref_pos in map(int, map(str.strip, match.group(1).split(','))):
            try:
                pmid = self.refs[ref_pos]
            except KeyError:
                logger.error(
                    "%s: no PMID for %s -- %s",
                    self.acc,
                    ref_pos,
                    self.refs,
                )
                continue

            refs.append(f"PMID:{pmid}")

        return f"[{','.join(refs)}]" if refs else ""


def get_default_values(
    signatures=False,
    clans=False
) -> tuple:
    """Return default values for Pfam variables

    :param signatures: bool - get default values for a signature
    :param clans: bool - get default values for a clan
    
    In order, for signatures return:
    Accession [str]
    name [int]
    Description [str]
    Long type (e.g. Family, Domain, etc.) [str]
    Abstract [empty str that is added]
    References {order added/pos ref [int]: pmid [int]}
    sequence ontology: None
    Author info [(author, orcid)]
    Build method: None
    Search method: None
    Sequence gathering threshold: None
    Domain gather threshold: None
    Wikipedia article: None
    Version: None

    For clans return:
    Clan accessio number, excluding version num [int]
    Clan id [int]
    Description [str]
    List of member accessions [List[str]]
    """
    if signatures:
        return None, None, None, None, "", {}, None, [], None, None, None, None, None, None
    if clans:
        return None, None, None, []


def get_signatures(pfam_path: str, persist_pfam=False) -> list[Method] | dict:
    """Parse Pfam-A.seed.gz file and extract signatures.

    :param pfam_path: str, path to Pfam-A.seed.gz file
    :param update_sigs: bool, parser called during member database
        update to persist the Pfam data in the oracle db
    
    Return 
    * list of Method objs
    * dict of entries
    """
    signatures = []
    entries = {}

    try:
        with gzip.open(pfam_path, 'rt') as fh:
            (
                accession,
                name,
                description,
                long_type,
                abstract,
                references,  # {order added/pos ref [int]: pmid [int]}
                sequence_ontology, author_info,
                build_method, search_method,
                sequence_ga, domain_ga,
                wiki, version,
            ) = get_default_values(signatures=True)

            for line in fh:
                if line.strip() == _RECORD_BREAK:
                    # create signature record from previous Pfam record
                    formatter = AbstractFormatter(references)

                    if persist_pfam:
                        entries[accession] = {
                            "curation": {
                                "sequence_ontology": sequence_ontology,
                                "authors": author_info,
                            },
                            "hmm": {
                                "commands": {
                                    "build": build_method,
                                    "search": search_method,
                                },
                                "cutoffs": {
                                    "gathering": {
                                        "sequence": sequence_ga,
                                        "domain": domain_ga
                                    }
                                },
                                "version": version,
                            },
                            "wiki": wiki,
                        }

                    else:
                        signatures.append(
                            Method(
                                accession,
                                _TYPES[long_type],  # convert to single letter
                                name,
                                description,
                                abstract=formatter.update(accession, abstract) if abstract else None,
                                references=list(references.values())
                            )
                        )

                    # start afresh for the next record/signature
                    (
                        accession,
                        name,
                        description,
                        long_type,
                        abstract,
                        references,  # {order added/pos ref [int]: pmid [int]}
                        sequence_ontology, author_info,
                        build_method, search_method,
                        sequence_ga, domain_ga,
                        wiki, version,
                    ) = get_default_values(signatures=True)
                
                elif line.startswith("#=GF AC"):
                    accession = line.split(" ")[-1].strip()
                    version = accession.split(".")[-1]
                
                elif line.startswith("#=GF ID"):
                    name = line[7:].strip()   # pfam id

                elif line.startswith("#=GF DE"):
                    description = line[7:].strip()

                elif line.startswith("#=GF TP"):
                    long_type = line[7:].strip()

                elif line.startswith("#=GF CC"):
                    abstract += line[7:].strip()

                elif line.startswith("#=GF RN"):
                    current_ref_pos = int(line.split("[")[-1][:-2])
                
                elif line.startswith("#=GF RM"):
                    references[current_ref_pos] = int(line[7:].strip())  # = pmid

                elif line.startswith("#=GF DR"):
                    sequence_ontology = line.split(";")[1].strip()

                elif line.startswith("#=GF AU"):
                    author_info.append(
                        (
                            line[7:].split(";")[0].strip(),  # author
                            line.split(";")[1].strip(),  # orcidID
                        )
                    )

                elif line.startswith("#=GF BM"):
                    build_method = line[7:].strip()

                elif line.startswith("#=GF SM"):
                    search_method = line[7:].strip()

                elif line.startswith("#=GF GA"):  # GA = gathering threshold
                    sequence_ga = line[7:].strip().split(" ")[0].strip()
                    domain_ga = line[7:].strip().split(" ")[1].replace(";","").strip()

                elif line.startswith("#=GF WK"):
                    wiki = line[7:].strip().replace(" ", "_")   # wikipedia article

    except UnicodeDecodeError:
        pass

    except FileNotFoundError:
        logger.error(
            (
                "Could not find Pfam-A-seed (seed alignment) file at %s\n"
                "Not retrieving signature data"
            ), pfam_path
        )
        if persist_pfam:
            return entries
        return signatures

    if persist_pfam:
        return entries
    return signatures


def get_fam_seq_counts(
    pfam_fasta_path: str
) -> dict[str: int] | None:
    """Parse PfamA-FASTA and get the number of seqs per PFAM family.

    :param pfam_fasta_path: path to Pfam-A.fasta.gz

    Return {fam accession [str] : number of seqs in fam [int]}
    """
    fams = {}

    try:
        with gzip.open(pfam_fasta_path, 'rt') as fh:
            for line in fh:
                if line.startswith(">"):
                    pfamA_acc = line.strip().split(" ")[-1].split(";")[0].split(".")[0]
                    try:
                        fams[pfamA_acc] += 1
                    except KeyError:
                        fams[pfamA_acc] = 1
    except UnicodeDecodeError:
        pass

    except FileNotFoundError:
        return

    return fams


def get_num_full(
    pfam_full_path: str,
) -> dict[str: int] | None:
    """Calculate num_full values by parsing Pfam-A.full.gz alignment file
    
    :param pfam_full_path: path to pfam-a-full alignemnt file
    """
    try:
        with gzip.open(pfam_full_path, 'rt') as fh:
            num_fulls = {}  # {acc [str]: num_full [int]}
            num_full_count = 0

            for line in fh:

                if line.strip() == _RECORD_BREAK:
                    # store the previous record
                    num_fulls[accession] = num_full_count
                    
                    # start afresh for the next record
                    num_full_count = 0

                elif line.startswith("#=GF AC"):
                    accession = line.split(" ")[-1].split(".")[0].strip()
                    version = accession.split(".")[-1]
                
                elif re.search(r"^#=GS\s+\S+\s+AC\s+", line.strip()):
                    num_full_count += 1
                    
    except UnicodeDecodeError:
        pass
    
    except FileNotFoundError:
        return

    return num_fulls


def get_clans(
    pfam_clan_path: str,
    pfam_fasta_path: str,
    pfam_full_path: str,
) -> list[Clan]:
    """Retrieve clans and clan members from the Pfam file Pfam-C
    
    :param pfam_clan_path: path to Pfam-C.tsv.gz file
    :param pfam_fasta_path: path to Pfam-A.fasta.gz file
    :param pfam_full_path: path to Pfam-A-full.gz alignment file
    """
    clans = []

    clan_file_found, fasta_file_found, full_file_found = True, True, True
    if os.path.isfile(pfam_clan_path) is False:
        logger.error("Could not find Pfam-C (clan file) at %s", pfam_clan_path)
        clan_file_found = False
    if os.path.isfile(pfam_fasta_path) is False:
        logger.error("Could not find Pfam-A-fasta file at %s", pfam_fasta_path)
        fasta_file_found = False
    if os.path.isfile(pfam_full_path) is False:
        logger.error("Could not find Pfam-A-full file at %s", pfam_full_path)
    if any(_ is False for _ in (clan_file_found, fasta_file_found, full_file_found)):
        logger.error("Check Pfam file paths are correct.\nNot retrieving clan data")
        return clans

    logger.info("Getting num_full values")
    num_fulls = get_num_full(pfam_full_path)
    if num_fulls is None:
        logger.error(
            (
                "Could not retrieve num_full values from Pfam-A-full (full alignment file) at %s\n"
                "Not retrieving clan data"
            ), pfam_full_path
        )
        return clans

    logger.info("Getting Pfam family seq counts")
    fam_seq_counts = get_fam_seq_counts(pfam_fasta_path)
    if fam_seq_counts is None:
        logger.error(
            (
                "Could not parse Pfam FASTA file at %s\n"
                "Not retrieving clan data"
            ), pfam_fasta_path
        )
        return clans

    logger.info("Getting Clans")
    try:
        with gzip.open(pfam_clan_path, 'rt') as fh:
            (
                accession,
                clan_id,
                description,
                members,
            ) = get_default_values(clans=True)

            for line in fh:

                if line.strip() == _RECORD_BREAK:
                    # store the previous record
                    clan = Clan(
                            accession,
                            clan_id,  # referred to as clan name in interpro
                            description,
                        )
                    
                    # totalSeq = totalSeq + pfam.pfama_acc.num_full; --> Pfam/PfamWeb/root/rest/clan/entry_xml.tt
                    clan_seq = 0
                    for pfamA_acc in fam_seq_counts:
                        try:
                            clan_seq += fam_seq_counts[pfamA_acc]
                        except KeyError:
                            logger.error(
                                (
                                    "Could not retrieve the number of sequences for Pfam family %s in clan %s -"
                                    "Setting number as 0"
                                ), pfamA_acc, accession
                            )
                            pass
                    
                    for pfamA_acc in members:
                        try:
                            num_full = num_fulls[pfamA_acc]  # num of seq accessions listed in full alignment (pfam-a.full)
                        except KeyError:
                            if pfamA_acc in dead_fams:
                                continue
                                
                            num_full = 0
                            logger.warning("Could not find num_full for %s in clan %s -> Setting score as '0'", pfamA_acc, accession)

                        clan.members.append(
                            {
                                "accession": pfamA_acc,
                                "score": num_full / clan_seq if clan_seq > 0 else 0
                            }
                        )

                    clans.append(clan)
                    
                    # start afresh for the next record
                    (
                        accession,
                        clan_id,
                        description,
                        members
                    ) = get_default_values(clans=True)
                
                elif line.startswith("#=GF AC"):
                    accession = line.split(" ")[-1].split(".")[0].strip()
                    version = accession.split(".")[-1]
                
                elif line.startswith("#=GF ID"):
                    clan_id = line[7:].strip()

                elif line.startswith("#=GF DE"):
                    description = line[7:].strip()

                elif line.startswith("#=GF MB"):
                    members.append(line[7:].strip().replace(";",""))

    except UnicodeDecodeError:
        pass
    
    except FileNotFoundError:
        logger.error(
            (
                "Could not parse Pfam-C (clan) file at %s\n"
                "Not retrieving clan data"
            ), pfam_clan_path
        )

    return list(clans.values())


def iter_protenn_matches(file: str):
    """Iterate ProtENN matches from the Pfam-N domain calls TSV file
    """
    with open(file, "rt") as fh:
        for line in fh:
            sequence_id, pfam_acc, start, end = line.rstrip().split("\t")
            if re.fullmatch(r"PF\d+", pfam_acc):
                yield sequence_id, pfam_acc, int(start), int(end)


def get_protenn_entries(cur, file: str) -> list[Method]:
    drop_table(cur, "INTERPRO.PFAMN_MATCH_TMP", purge=True)
    cur.execute(
        """
        CREATE TABLE INTERPRO.PFAMN_MATCH_TMP (
            PROTEIN_ID VARCHAR(15) NOT NULL,
            METHOD_AC VARCHAR2(25) NOT NULL,
            POS_FROM NUMBER(5) NOT NULL,
            POS_TO NUMBER(5) NOT NULL
        ) NOLOGGING
        """
    )

    query = """
        INSERT /*+ APPEND */ 
        INTO INTERPRO.PFAMN_MATCH_TMP 
        VALUES (:1, :2, :3, :4)
    """

    with Table(con=cur.connection, query=query, autocommit=True) as table:
        it = iter_protenn_matches(file)

        sequence_id, pfam_acc, start, end = next(it)
        accessions = {pfam_acc}
        uniparc = sequence_id.startswith("UPI")
        table.insert((sequence_id, pfam_acc, start, end))

        for sequence_id, pfam_acc, start, end in it:
            accessions.add(pfam_acc)
            table.insert((sequence_id, pfam_acc, start, end))

    drop_table(cur, "INTERPRO.PFAMN_MATCH", purge=True)
    if uniparc:
        cur.execute(
            """
            CREATE TABLE INTERPRO.PFAMN_MATCH NOLOGGING
            AS SELECT * FROM INTERPRO.PFAMN_MATCH_TMP WHERE 1 = 0
            """
        )

        cur.execute(
            """
            INSERT /*+ APPEND */ INTO INTERPRO.PFAMN_MATCH
            SELECT X.AC, M.METHOD_AC, M.POS_FROM, M.POS_TO
            FROM INTERPRO.PFAMN_MATCH_TMP M
            INNER JOIN UNIPARC.XREF X ON M.PROTEIN_ID = X.UPI
            WHERE X.DBID IN (2, 3) AND X.DELETED = 'N'
            """
        )

        cur.connection.commit()
        drop_table(cur, "INTERPRO.PFAMN_MATCH_TMP", purge=True)
    else:
        cur.execute("RENAME PFAMN_MATCH_TMP TO PFAMN_MATCH")

    return [Method(pfam_acc, sig_type=None) for pfam_acc in accessions]


def persist_extra_pfam_data(
    db_props: dict[str, str],
    ora_url: str,
):
    """
    Persist additional Pfam data from Pfam release files that is not captured in the signatures
        or clan tables. Primarily persists data required by interpro7.

    :param db_props: db properties from members.conf file
    :param con: oracle db connection
    """
    signatures = get_signatures(
        pfam_path=db_props["seed"],
        persist_pfam=True,
    )

    pfam_query = """
        INSERT /*+ APPEND */ 
        INTO INTERPRO.PFAM_DATA 
        VALUES (:1, :2, :3, :4, :5, :6, :7)
    """
    author_query = """
        INSERT /*+ APPEND */ 
        INTO INTERPRO.PFAM_AUTHOR
        VALUES (:1, :2, :3)
    """
    wiki_query = """
        INSERT /*+ APPEND */ 
        INTO INTERPRO.PFAM_WIKIPEDIA 
        VALUES (:1, :2)
    """

    con = oracledb.connect(ora_url)
    with con.cursor() as cur:
        drop_table(cur, "INTERPRO.PFAM_DATA", purge=True)
        cur.execute(
            """
            CREATE TABLE INTERPRO.PFAM_DATA (
                accession VARCHAR2(25) PRIMARY KEY,
                seq_ontology NUMBER,
                hmm_build VARCHAR2(255),
                hmm_search VARCHAR2(255),
                seq_gathering FLOAT,
                domain_gathering FLOAT,
                version NUMBER
            )
            """
        )

        drop_table(cur, "INTERPRO.PFAM_AUTHOR", purge=True)
        cur.execute(
            """
            CREATE TABLE INTERPRO.PFAM_DATA (
                accession VARCHAR2(25) PRIMARY KEY,
                author VARCHAR2(225),
                orcid VARCHAR2(225)
            )
            """
        )

        drop_table(cur, "INTERPRO.PFAM_WIKIPEDIA", purge=True)
        cur.execute(
            """
            CREATE TABLE INTERPRO.PFAM_WIKIPEDIA (
                accession VARCHAR2(25) PRIMARY KEY,
                title VARCHAR2(225)
            )
            """
        )
        
        for pfam_acc in signatures:
            cur.execute(
                pfam_query,
                [
                    pfam_acc,
                    signatures[pfam_acc]["curation"]["sequence_ontology"],
                    signatures[pfam_acc]["hmm"]["commands"]["build"],
                    signatures[pfam_acc]["hmm"]["commands"]["search"],
                    signatures[pfam_acc]["hmm"]["cutoffs"]["gathering"]["sequence"],
                    signatures[pfam_acc]["hmm"]["cutoffs"]["gathering"]["domain"],
                    signatures[pfam_acc]["hmm"]["version"],
                ]
            )

            for author_info in signatures[pfam_acc]["curation"]["authors"]:
                cur.execute(
                    author_query,
                    author_info
                )
            
            cur.execute(
                wiki_query,
                signatures[pfam_acc]["wiki"]
            )
