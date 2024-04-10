import gzip
import os
import re

from pathlib import Path

import oracledb

from .common import Clan, Method
from pyinterprod import logger
from pyinterprod.utils import Table
from pyinterprod.utils.oracle import drop_table


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
    clans=False,
    clans_lit=False,
) -> tuple:
    """Return default values for Pfam variables

    :param signatures: bool - get default values for a signature
    :param clans: bool - get default values for a clan
    :param clans_lit: bool - get clan literature data
    
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
    Wikipedia article: [str] -- can be multiple wiki articles
    Version: None

    For clans return:
    Clan accessio number, excluding version num [int]
    Clan id [int]
    Description [str]
    List of member accessions [List[str]]

    For clans lit:
    Clan accession (excluding version) number [int]
    Authors [str]
    References [dict]
    """
    if signatures:
        return None, None, None, None, "", {}, None, [], None, None, None, None, [], None
    if clans:
        return None, None, None, []
    if clans_lit:
        return None, None, {}


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
    line_count = 0

    try:
        with gzip.open(pfam_path, 'rb') as fh:
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

            for _line in fh:
                line_count += 1
                try:
                    line = _line.decode('utf-8')
                except UnicodeDecodeError:
                    try:
                        line = _line.decode('latin-1')
                    except UnicodeDecodeError:
                        logger.error(
                            "UnicodeDecodeError encountered on PFAM-A.seed line %s. Skipping line.",
                            line_count
                        )
                        continue

                if line.strip() == _RECORD_BREAK:
                    # create signature record from previous Pfam record
                    formatter = AbstractFormatter(references)

                    if persist_pfam:
                        entries[accession] = {
                            "name": name,
                            "description": description,
                            "type": _TYPES[long_type],
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
                            "wiki": wiki,  # list, can be multiple articles
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
                    wiki.append(line[7:].strip().replace(" ", "_"))   # wikipedia article

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


def get_fam_seq_counts(
    pfam_fasta_path: str
) -> dict[str: int] | None:
    """Parse PfamA-FASTA and get the number of seqs per PFAM family.

    :param pfam_fasta_path: path to Pfam-A.fasta.gz

    Return {fam accession [str] : number of seqs in fam [int]}
    """
    fams = {}
    line_count = 0

    try:
        with gzip.open(pfam_fasta_path, 'rb') as fh:
            for _line in fh:
                line_count += 1
                try:
                    line = _line.decode('utf-8')
                except UnicodeDecodeError:
                    try:
                        line = _line.decode('latin-1')
                    except UnicodeDecodeError:
                        logger.error(
                            "UnicodeDecodeError encountered on PFAM-A.fasta line %s. Skipping line.",
                            line_count
                        )
                        continue
                if line.startswith(">"):
                    pfamA_acc = line.strip().split(" ")[-1].split(";")[0].split(".")[0]
                    try:
                        fams[pfamA_acc] += 1
                    except KeyError:
                        fams[pfamA_acc] = 1

    except FileNotFoundError:
        return

    return fams


def get_num_full(
    pfam_full_path: str,
    seed=False,
) -> dict[str: int] | None:
    """Calculate num_full values by parsing Pfam-A.full.gz alignment file
    
    :param pfam_full_path: path to pfam-a-full alignemnt file
    """
    line_count = 0
    try:
        with gzip.open(pfam_full_path, 'rb') as fh:
            num_fulls = {}  # {acc [str]: num_full [int]}
            num_full_count = 0

            for _line in fh:
                line_count += 1
                try:
                    line = _line.decode('utf-8')
                except UnicodeDecodeError:
                    try:
                        line = _line.decode('latin-1')
                    except UnicodeDecodeError:
                        logger.error(
                            "UnicodeDecodeError encountered on PFAM-A.%s line %s. Skipping line.",
                            "seed" if seed else "full", line_count
                        )
                        continue

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
    line_count = 0
    try:
        with gzip.open(pfam_clan_path, 'rb') as fh:
            (
                accession,
                clan_id,
                description,
                members,
            ) = get_default_values(clans=True)

            for _line in fh:
                try:
                    line = _line.decode('utf-8')
                except UnicodeDecodeError:
                    try:
                        line = _line.decode('latin-1')
                    except UnicodeDecodeError:
                        logger.error(
                            "UnicodeDecodeError encountered on PFAM-C line %s. Skipping line.",
                            line_count
                        )
                        continue
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
    
    except FileNotFoundError:
        logger.error(
            (
                "Could not parse Pfam-C (clan) file at %s\n"
                "Not retrieving clan data"
            ), pfam_clan_path
        )

    return clans


def get_clan_literature(pfam_clan_path: str) -> dict:
    """Get clan authors and literature data from Pfam-C file

    :param pfam_clan_path: path to Pfam-C release file.
    """
    clans = {}
    line_count = 0
    try:
        with gzip.open(pfam_clan_path, 'rb') as fh:
            accession, authors, references = get_default_values(clans_lit=True)

            for _line in fh:
                line_count += 1
                try:
                    line = _line.decode('utf-8')
                except UnicodeDecodeError:
                    try:
                        line = _line.decode('latin-1')
                    except UnicodeDecodeError:
                        logger.error(
                            "UnicodeDecodeError encountered on PFAM-C line %s. Skipping line.",
                            line_count
                        )
                        continue
            
                if line.strip() == _RECORD_BREAK:
                    clans[accession] = {
                        'authors': authors,
                        'references': references,
                    }
                    # set to default values for next accession
                    accession, authors, references = get_default_values(clans_lit=True)

                elif line.startswith("#=GF AC"):
                    accession = line.split(" ")[-1].split(".")[0].strip()

                elif line.startswith('#=GF AU'):
                    authors = line[7:].split(",")
                    authors = [_.strip() for _ in authors]

                elif line.startswith("#=GF RN"):
                    current_ref_pos = int(line.split("[")[-1][:-2])
                    references[current_ref_pos] = {
                        'pmid': None,
                        'authors': '',
                        'title': '',
                        'journal': '',
                    }

                elif line.startswith("#=GF RM"):
                    references[current_ref_pos]['pmid'] = int(line[7:].strip())

                elif line.startswith("#=GF RA"):
                    references[current_ref_pos]['authors'] += line[7:-1].strip()

                elif line.startswith("#=GF RT"):
                    references[current_ref_pos]['title'] += line[7:].strip()

                elif line.startswith("#=GF RL"):
                    references[current_ref_pos]['journal'] += line[7:].split(";")[0][:-4].strip()

    except FileNotFoundError:
        logger.error(
            (
                "Could not parse Pfam-C (clan) file at %s\n"
                "Not retrieving clan data"
            ), pfam_clan_path
        )
    return clans


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
    :param ora_url: connection str for oracle db
    """
    def add_num_val(acc, num_dict, _record):
        """Add alignment count if available to record, else add None"""
        try:
            _record.append(num_dict[acc])
        except KeyError:
            _record.append(None)
        return record

    signatures = get_signatures(
        pfam_path=db_props["seed"],
        persist_pfam=True,
    )
    clans = get_clan_literature(db_props['clan'])

    logger.info("Getting pfam-A.seed alignment counts")
    seed_num = get_num_full(db_props["seed"], seed=True)
    if seed_num is None:
        logger.error("Could not find file pfam-A.seed at %s\nNot persisting seed_num values", db_props["seed"])
        seed_num = {}
    logger.info("Getting pfam-A.full alignment counts")
    full_num = get_num_full(db_props["full"])
    if full_num is None:
        logger.error("Could not find file pfam-A.seed at %s\nNot persisting full_num values", db_props["full"])
        full_num = {}
    logger.info("Getting RP15 alignment counts")
    rp15_num = get_alignment_counts(db_props["rp15"], "rp15")
    logger.info("Getting RP35 alignment counts")
    rp35_num = get_alignment_counts(db_props["rp35"], "rp35")
    logger.info("Getting RP55 alignment counts")
    rp55_num = get_alignment_counts(db_props["rp55"], "rp55")
    logger.info("Getting RP75 alignment counts")
    rp75_num = get_alignment_counts(db_props["rp75"], "rp75")
    logger.info("Getting UniProt alignment counts")
    uniprot_num = get_alignment_counts(db_props["uniprot"], "uniprot")

    logger.info("Persisting data for %s signatures", len(signatures))
    logger.info("Persisting data for %s clans", len(clans))

    pfam_query = """
        INSERT /*+ APPEND */ 
        INTO INTERPRO.PFAM_DATA 
        VALUES (:1, :2, :3, :4, :5, :6, :7, :8, :9, :10, :11, :12, :13, :14, :15, :16, :17)
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
    clan_au_query = """
        INSERT /*+ APPEND */ 
        INTO INTERPRO.PFAM_CLAN_AUTHOR
        VALUES (:1, :2)
    """
    clan_lit_query = """
        INSERT /*+ APPEND */ 
        INTO INTERPRO.PFAM_CLAN_LITERATURE
        VALUES (:1, :2, :3, :4, :5, :6)
    """

    con = oracledb.connect(ora_url)
    with con.cursor() as cur:
        drop_table(cur, "INTERPRO.PFAM_DATA", purge=True)
        cur.execute(
            """
            CREATE TABLE INTERPRO.PFAM_DATA (
                accession VARCHAR2(25) PRIMARY KEY,
                name VARCHAR2(255),
                description VARCHAR2(500),
                type VARCHAR2(25),
                seq_ontology NUMBER,
                hmm_build VARCHAR2(255),
                hmm_search VARCHAR2(255),
                seq_gathering FLOAT,
                domain_gathering FLOAT,
                version NUMBER,
                seed_num NUMBER, full_num NUMBER,
                rp15_num NUMBER, rp35_num NUMBER,
                rp55_num NUMBER, rp75_num NUMBER,
                uniprot_num NUMBER
            )
            """
        )

        drop_table(cur, "INTERPRO.PFAM_AUTHOR", purge=True)
        cur.execute(
            """
            CREATE TABLE INTERPRO.PFAM_AUTHOR (
                accession VARCHAR2(25),
                author VARCHAR2(225),
                orcid VARCHAR2(225)
            )
            """
        )

        drop_table(cur, "INTERPRO.PFAM_WIKIPEDIA", purge=True)
        cur.execute(
            """
            CREATE TABLE INTERPRO.PFAM_WIKIPEDIA (
                accession VARCHAR2(25),
                title VARCHAR2(225)
            )
            """
        )

        drop_table(cur, "INTERPRO.PFAM_CLAN_AUTHOR", purge=True)
        cur.execute(
            """
            CREATE TABLE INTERPRO.PFAM_CLAN_AUTHOR (
                clan_id VARCHAR2(25),
                author VARCHAR2(225)
            )
            """
        )

        drop_table(cur, "INTERPRO.PFAM_CLAN_LITERATURE", purge=True)
        cur.execute(
            """
            CREATE TABLE INTERPRO.PFAM_CLAN_LITERATURE (
                clan_id VARCHAR2(25),
                pubmed_id NUMBER,
                title VARCHAR2(500),
                author VARCHAR2(500),
                journal VARCHAR2(225),
                order_added NUMBER
            )
            """
        )
        
        for pfam_acc in signatures:
            record = [
                pfam_acc.split(".")[0],
                signatures[pfam_acc]['name'],
                signatures[pfam_acc]['description'],
                signatures[pfam_acc]['type'],
                signatures[pfam_acc]["curation"]["sequence_ontology"],
                signatures[pfam_acc]["hmm"]["commands"]["build"],
                signatures[pfam_acc]["hmm"]["commands"]["search"],
                signatures[pfam_acc]["hmm"]["cutoffs"]["gathering"]["sequence"],
                signatures[pfam_acc]["hmm"]["cutoffs"]["gathering"]["domain"],
                signatures[pfam_acc]["hmm"]["version"],
            ]
            record = add_num_val(pfam_acc.split(".")[0], seed_num, record)
            record = add_num_val(pfam_acc.split(".")[0], full_num, record)
            record = add_num_val(pfam_acc.split(".")[0], rp15_num, record)
            record = add_num_val(pfam_acc.split(".")[0], rp35_num, record)
            record = add_num_val(pfam_acc.split(".")[0], rp55_num, record)
            record = add_num_val(pfam_acc.split(".")[0], rp75_num, record)
            record = add_num_val(pfam_acc.split(".")[0], uniprot_num, record)

            cur.execute(
                pfam_query,
                record,
            )

            for author_info in signatures[pfam_acc]["curation"]["authors"]:
                cur.execute(
                    author_query,
                    (pfam_acc.split(".")[0],) + author_info
                )
            
            for wiki_title in signatures[pfam_acc]["wiki"]:
                cur.execute(
                    wiki_query,
                    [pfam_acc.split(".")[0], wiki_title]
                )
    
        for clan_id in clans:
            for author in clans[clan_id]['authors']:
                cur.execute(clan_au_query, [clan_id, author])
            
            for order_added in clans[clan_id]['references']:
                try:
                    cur.execute(
                        clan_lit_query,
                        [
                            clan_id,
                            clans[clan_id]['references'][order_added]['pmid'],
                            _check_str_len(clan_id, clans[clan_id]['references'][order_added]['title'], "TITLE", 500),
                            _check_str_len(clan_id, clans[clan_id]['references'][order_added]['authors'][:-1], "AUTHORS", 500),
                            clans[clan_id]['references'][order_added]['journal'],
                            order_added
                        ]
                    )
                except oracledb.DatabaseError as err:
                    logger.error(
                        "Could not insert citation (order_added:%s) for clan %s\nERROR:%s",
                        order_added, clan_id, err
                    )

    con.commit()

    try:
        db_props["alignment"]
        persist_alignments(db_props, con)
    except KeyError:
        pass

    con.close()


def _check_str_len(clan_id: str, field: str, field_name: str, max_len: int) -> str:
    """Check if the len of a field is too long for the db. If too long, flag, and truncate.

    :param clan_id: clan accession
    :param field: str to be inserted into db
    :param field_name: name of column to be inserted into
    :param max_len: int, max number of accepted chars
    """
    if len(field) > max_len:
        logger.error(
            (
                "The %s for clan %s is too long (actual: %s, maximum: %s)\n"
                "Truncating %s to %s chars"
            ), field_name, clan_id, len(field), max_len, field_name, max_len
        )
        return field[:max_len]
    return field


def get_alignment_counts(pfam_file: str, alignment_name: str) -> dict[str, int]:
    """Count the number of hits in the alignemnt file for each accession

    :param pfam_file: path to rp15, rp35, rp55, rp75, uniprot
    :param alignment_name: name of alignment: rp15, rp35, rp55, rp75, uniprot
    
    Return {accession: count [int]}
    """
    line_count, count_dict = 0, {}
    try:
        with gzip.open(pfam_file, 'rb') as fh:
            count, acc = 0, ""

            for _line in fh:
                line_count += 1
                try:
                    line = _line.decode("utf-8")
                except UnicodeDecodeError:
                    try:
                        line = _line.decode('latin-1')
                    except UnicodeDecodeError:
                        logger.error(
                            (
                                "Unicode Decode Error encountered in PFAM %s file line %s\n"
                                "Not parsing this line and continuing to the next"
                            ),
                            alignment_name, line_count
                        )
                        continue

                if line.strip() == _RECORD_BREAK:
                    count_dict[acc] = count
                    count, acc = 0, ""

                elif line.startswith("#=GF AC"):
                    acc = line.split(" ")[-1].split(".")[0].strip()

                elif line.startswith("#") is False:
                    count += 1

    except FileNotFoundError:
        logger.error("Could not find Pfam file (%s) at %s", alignment_name, pfam_file)
    
    return count_dict


def persist_alignments(db_props: dict[str, str], con) -> None:
    """Persist Pfam alignments in the oracle db.

    Alignments types to be stored are defined by assigning a value to the 
    corresponding file/alignment type in members.conf:
        seed, full, rpXX (where XX is a number), uniprot

    :param db_props: db properties from members.conf file
    :param con: open oracle db connection
    """
    def get_ali_paths(directory: Path):
        return (entry for entry in Path(directory).iterdir() if entry.is_file())

    types_to_store = [
        _ for _ in db_props 
        if db_props[_] is not None 
        and _ not in ['hmm', 'members', 'signatures', 'clan', 'alignment']
    ]
    if len(types_to_store) == 0:
        logger.warning("No alignment types specified. Not storing alignments")
        return
    
    logger.info("Storing alignments: %s", types_to_store)
    sql_insert = "INSERT INTO INTERPRO.PFAM_ALIGNMENTS (accession, type, alignment) VALUES (:1, :2, :3)"

    with con.cursor() as cur:
        drop_table(cur, "INTERPRO.PFAM_ALIGNMENTS", purge=True)
        cur.execute(
            """
            CREATE TABLE INTERPRO.PFAM_ALIGNMENTS (
                accession VARCHAR2(25),
                type VARCHAR2(35),
                alignment BLOB
            )
            """
        )

        dir_paths = (entry for entry in Path(db_props['alignment']).iterdir() if entry.is_dir())
        for _path in dir_paths:
            pfam_acc = _path.name
            for alignment_path in get_ali_paths(_path):
                ali_type = alignment_path.name.split(".")[1]
                if ali_type in types_to_store:
                    with gzip.open(alignment_path, 'rb') as fh:
                        alignment = fh.read()

                        cur.execute(sql_insert, [pfam_acc, ali_type, alignment])

    con.commit()
