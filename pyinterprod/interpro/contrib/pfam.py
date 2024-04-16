import gzip
import os
import re
from io import BytesIO
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
    def __init__(self, references: dict[int, int]):
        self.refs = references
        self.acc = None

    def update(self, acc: str, abstract: str) -> str:
        self.acc = acc
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


def _decode_line(_line):
    try:
        line = _line.decode('utf-8')
    except UnicodeDecodeError:
        try:
            line = _line.decode('latin-1')
        except UnicodeDecodeError:
            line = None
    return line


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
                line = _decode_line(_line)
                if not line:
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
                line = _decode_line(_line)
                if not line:
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
                line_count += 1
                line = _decode_line(_line)
                if not line:
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
                line = _decode_line(_line)
                if not line:
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

    try:
        logger.info("Getting signatures from %s", Path(db_props["seed"]).name)
    except KeyError:
        return

    signatures = get_signatures(
        pfam_path=db_props["seed"],
        persist_pfam=True,
    )
    clans = get_clan_literature(db_props['clan'])

    logger.info("Getting Pfam-A.seed alignment counts")
    seed_num = get_num_full(db_props["seed"], seed=True)
    if seed_num is None:
        logger.error("Could not find file pfam-A.seed at %s\nNot persisting seed_num values", db_props["seed"])
        seed_num = {}

    try:
        logger.info("Getting %s alignment counts", Path(db_props["full"]).name)
        full_num = get_num_full(db_props["full"])
        if full_num is None:
            logger.error("Could not find file pfam-A.seed at %s\nNot persisting full_num values", db_props["full"])
            full_num = {}
    except KeyError:
        logger.error(
            (
                f"Pfam file path for full alignments (Pfam-A.full) is not defined in "
                "the members.config file.\n"
                f"Skipping retrieving full (Pfam-A.full) alignemnts."
            ),
        )
        full_num = {}

    logger.info("Getting RP15 alignment counts")
    rp15_num = get_alignment_counts(db_props, "rp15")
    logger.info("Getting RP35 alignment counts")
    rp35_num = get_alignment_counts(db_props, "rp35")
    logger.info("Getting RP55 alignment counts")
    rp55_num = get_alignment_counts(db_props, "rp55")
    logger.info("Getting RP75 alignment counts")
    rp75_num = get_alignment_counts(db_props, "rp75")
    logger.info("Getting UniProt alignment counts")
    uniprot_num = get_alignment_counts(db_props, "uniprot")

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
                author VARCHAR2(750),
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
                            _check_str_len(clan_id, clans[clan_id]['references'][order_added]['authors'][:-1], "AUTHORS", 750),
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

    n = persist_alignments(db_props["alignment"], con)
    logger.info(f"{n:,} alignments persisted")
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


def get_alignment_counts(db_props: dict[str, str], alignment_name: str) -> dict[str, int]:
    """Count the number of hits in the alignemnt file for each accession

    :param db_props: db properties from members.conf file
    :param alignment_name: name of alignment: rp15, rp35, rp55, rp75, uniprot
    
    Return {accession: count [int]}
    """
    line_count, count_dict = 0, {}
    try:
        pfam_file = db_props[alignment_name]
    except KeyError:
        logger.error(
            (
                f"Pfam file path for {alignment_name} alignments is not defined in "
                "the members.config file.\n"
                f"Skipping retrieving {alignment_name} alignemnts."
            ),
        )
        return count_dict

    try:
        with gzip.open(pfam_file, 'rb') as fh:
            count, acc = 0, ""

            for _line in fh:
                line_count += 1
                line = _decode_line(_line)
                if not line:
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


def persist_alignments(alndir: str, con: oracledb.Connection) -> int:
    """Persist Pfam alignments in Oracle

    :param alndir: string representing the path to the directory of alignments
    :param con: open oracle db connection
    """
    count = 0
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

        for accession, aln_type, aln_bytes in find_alignments(alndir):
            cur.execute(
                """
                INSERT INTO INTERPRO.PFAM_ALIGNMENTS
                VALUES (:1, :2, :3)
                """,
                [accession, aln_type, aln_bytes]
            )
            count += 1

    con.commit()
    return count


def find_alignments(root: str):
    for path in Path(root).rglob("*.gz"):
        accession, aln_type, _ = path.name.split(".")
        if aln_type in _ALN_TYPES:
            with gzip.open(path, "rb") as fh:
                yield accession, aln_type, fh.read()


def get_signatures_alt(pfam_seed_file: str) -> list[Method]:
    methods = []
    for entry in parse_sto(pfam_seed_file):
        comment = entry.features.get("CC")
        rn2pmid = {}
        if comment:
            for i, obj in enumerate(entry.references):
                rn2pmid[i+1] = int(obj["RM"])

            abstract = _repl_references(comment, rn2pmid)
        else:
            abstract = None

        accession, version = entry.features["AC"].split(".")
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


def _repl_references(text: str, references: dict[int, int]):
    def _repl(match: re.Match) -> str:
        refs = []
        for ref_num in map(int, map(str.strip, match.group(1).split(','))):
            pmid = references[ref_num]
            refs.append(f"PMID:{pmid}")

        return f"[{','.join(refs)}]"

    return re.sub(r"\[([\d\s,]+)]", _repl, text)


def _init_entry() -> dict:
    return {
        "accession": None,
        "name": None,
        "description": None,
        "type": None,
        "curation": {
            "sequence_ontology": None,
            "authors": [],
            "comment": None,
            "references": {}
        },
        "hmm": {
            "commands": {
                "build": None,
                "search": None,
            },
            "cutoffs": {
                "gathering": {
                    "sequence": None,
                    "domain": None
                }
            },
            "version": None,
        },
        "wikipedia": [],
        "alignments": [],
    }


_KEPT_FIELDS = {"AC", "ID", "DE", "AU", "BM", "SM", "GA", "TC", "TP", "SQ",
                "DR", "CC", "WK", "CL", "MB"}
_REFERENCE_FIELDS = {"RC", "RM", "RT", "RA", "RL"}


class StockholdMSA:
    def __init__(self):
        self.features = {}
        self.references = []
        self.sequences = []

    def add_feature(self, name: str, value: str):
        if name == "RN":
            self.references.append({})
        elif name in _REFERENCE_FIELDS:
            obj = self.references[-1]
            obj[name] = self.concat(obj.get(name, ""), value)
        elif name in _KEPT_FIELDS:
            self.features[name] = self.concat(self.features.get(name, ""),
                                              value)

    @staticmethod
    def concat(s: str, s2: str) -> str:
        if s and s[-1].isspace():
            return s + s2
        else:
            return s + " " + s2

    def get_alignments(self, compresslevel: int = 9) -> bytes:
        with BytesIO() as bs:
            with gzip.GzipFile(mode="wb",
                               compresslevel=compresslevel,
                               fileobj=bs) as gz:
                for line in self.sequences:
                    gz.write((line + "\n").encode("utf-8"))

            data = bs.getvalue()
        return data

    def is_empty(self) -> bool:
        return (not self.features and
                not self.references and
                not self.sequences)


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
            if line == "//":
                yield entry
                entry = StockholdMSA()
            elif line.startswith("#=GF"):
                # Generic per-File
                _, field, value = line.split(maxsplit=2)
                entry.add_feature(field, value)
                # # Compulsory fields
                # if field == "ID":
                #     entry["name"] = value
                # elif field == "AC":
                #     accession, version = value.split(".")
                #     entry["accession"] = accession
                #     entry["hmm"]["version"] = version
                # elif field == "DE":
                #     entry["description"] = value
                # elif field == "AU":
                #     author, orcid = value.split(";")
                #     entry["curation"]["authors"].append({
                #         "author": author,
                #         "orcid": orcid or None
                #     })
                # elif field == "BM":
                #     entry["hmm"]["commands"] = value
                # elif field == "SM":
                #     entry["hmm"]["search"] = value
                # elif field == "GA":
                #     seq_ga, dom_ga = value.rstrip(";").split()
                #     entry["hmm"]["cutoffs"]["sequence"] = float(seq_ga)
                #     entry["hmm"]["cutoffs"]["domain"] = float(dom_ga)
                # elif field == "TP":
                #     entry["type"] = value
                # # Optional fields
                # elif field == "DR":
                #     m = re.match(r"SO;\s*(\d+);", value)
                #     if m:
                #         entry["curation"]["sequence_ontology"] = m.group(1)
                # elif field == "RN":
                #     # expects [1], [2], etc.
                #     current_ref = int(value[1:-1])
                #     entry["curation"]["references"][current_ref] = None
                # elif field == "RM":
                #     references = entry["curation"]["references"]
                #     # Ensure there is only one Medline number per reference
                #     assert references[current_ref] is None
                #     references[current_ref] = int(value)
                # elif field == "CC":
                #     if entry["curation"]["comment"] is None:
                #         entry["curation"]["references"] = ""
                #     else:
                #         entry["curation"]["references"] += " "
                #
                #     entry["curation"]["references"] += value
                # elif field == "WK":
                #     entry["wikipedia"].append(value.replace(" ", "_"))
            elif line.startswith("#=GS"):
                # Generic per-Sequence
                pass
            elif line.startswith("#=GR"):
                # Generic per-Residue
                pass
            elif line.startswith("#=GC"):
                # Generic per-Column
                pass
            else:
                # Sequence line
                entry.sequences.append(line)

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
                pfam_acc = entry.split(".")[0]
                try:
                    counts[pfam_acc] += 1
                except KeyError:
                    counts[pfam_acc] = 1

    return counts


def _decode(b: bytes) -> str:
    try:
        s = b.decode("utf-8")
    except UnicodeDecodeError:
        pass
    else:
        return s.rstrip()

    return b.decode("latin-1").rstrip()


def _is_empty(obj) -> bool:
    if isinstance(obj, dict):
        for v in obj.values():
            if not _is_empty(v):
                return False

        return True
    elif isinstance(obj, list):
        return len(obj) == 0
    else:
        return obj is None


def compress_alignments(alignments: list[str]) -> bytes:
    with BytesIO() as bs:
        with gzip.GzipFile(mode="wb", compresslevel=9, fileobj=bs) as gz:
            for line in alignments:
                gz.write((line + "\n").encode("utf-8"))

        data = bs.getvalue()
    return data
