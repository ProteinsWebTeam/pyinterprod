import re

import MySQLdb
import MySQLdb.cursors

from .common import Clan, Method
from pyinterprod import logger
from pyinterprod.utils import Table
from pyinterprod.utils.oracle import drop_table
from pyinterprod.utils.pg import url2dict


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

    For clans return:
    Clan accessio number, excluding version num [int]
    Clan id [int]
    Description [str]
    List of member accessions [List[str]]
    """
    if signatures:
        return None, None, None, None, "", {}
    if clans:
        return None, None, None, []


def get_signatures(pfam_path: str) -> list[Method]:
    """Parse Pfam-A.seed.gz file and extract signatures.

    :param pfam_path: str, path to Pfam-A.seed.gz file
    
    Return list of Method objs
    """
    signatures =[]

    try:
        with gzip.open(pfam_path, 'rt') as fh:
            (
                accession,
                name,
                description,
                long_type,
                abstract,
                references,  # {order added/pos ref [int]: pmid [int]}
            ) = get_default_values(signatures=True)

            for line in fh:
                if line.strip() == RECORD_BREAK:
                    # create signature record from previous Pfam record
                    formatter = AbstractFormatter(references)

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

                i += 1
                if limit is not None:
                    if i > limit:
                        break

    except UnicodeDecodeError:
        pass

    except FileNotFoundError:
        logger.error(
            (
                "Could not find Pfam-A-seed (seed alignment) file at %s\n"
                "Not retrieving signature data"
            ), pfam_path
        )
        return signatures

    return signatures




def get_clans(uri: str) -> list[Clan]:
    con = connect_mysql(uri)
    cur = MySQLdb.cursors.SSCursor(con)
    cur.execute(
        """
        SELECT
          c.clan_acc, c.clan_id, c.clan_description,
          c.number_sequences, m.pfamA_acc, f.num_full
        FROM clan c
        INNER JOIN clan_membership m ON c.clan_acc = m.clan_acc
        INNER JOIN pfamA f ON m.pfamA_acc = f.pfamA_acc
        """
    )

    clans = {}
    for row in cur:
        clan_acc = row[0]

        try:
            clan = clans[clan_acc]
        except KeyError:
            clan = clans[clan_acc] = Clan(clan_acc, row[1], row[2])

        clan.members.append({
            "accession": row[4],
            "score": row[5] / row[3] if row[3] > 0 else 0
        })

    cur.close()
    con.close()

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
