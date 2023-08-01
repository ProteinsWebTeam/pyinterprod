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


def connect_mysql(url: str) -> MySQLdb.Connection:
    obj = url2dict(url)
    # MySQLdb and psycopg use different keyword params for password/database
    obj.update({
        "passwd": obj.pop("password"),
        "db": obj.pop("dbname")
    })
    return MySQLdb.connect(**obj)


def get_clans(url: str) -> list[Clan]:
    con = connect_mysql(url)
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


class AbstractFormater:
    def __init__(self, references, citations):
        self.refs = references
        self.cite = citations
        self.acc = None

    def update(self, acc, abstract):
        self.acc = acc
        return '<p>' + re.sub(r"\[([\d\s,]+)\]", self.sub, abstract) + '</p>'

    def sub(self, match):
        refs = []
        for ref_pos in map(int, map(str.lower, match.group(1).split(','))):
            try:
                ref_id = self.refs[self.acc][ref_pos]
            except KeyError:
                logger.error(f"{self.acc}: no reference for {ref_pos}")
                refs.append("PMID:")
                continue

            try:
                pmid = self.cite[ref_id]
            except KeyError:
                logger.error(f"{ref_id}: no PMID")
                refs.append("PMID:")
                continue
            else:
                refs.append(f"PMID:{pmid}")

        return f"[{', '.join(refs)}]"


def get_signatures(url: str) -> list[Method]:
    con = connect_mysql(url)
    cur = con.cursor()
    cur.execute(
        """
        SELECT auto_lit, pmid
        FROM literature_reference
        """
    )
    citations = dict(cur.fetchall())

    cur.execute(
        """
        SELECT pfamA_acc, auto_lit, order_added
        FROM pfamA_literature_reference
        """
    )
    references = {}
    for acc, ref_id, ref_pos in cur:
        if acc in references:
            references[acc][ref_pos] = ref_id
        else:
            references[acc] = {ref_pos: ref_id}

    cur.execute(
        """
        SELECT pfamA_acc, pfamA_id, description, type, comment
        FROM pfamA
        """
    )
    entries = {row[0]: row[1:] for row in cur}
    cur.close()
    con.close()

    formater = AbstractFormater(references, citations)
    signatures = []
    for acc in sorted(entries):
        name, description, long_type, abstract = entries[acc]
        short_type = _TYPES[long_type]

        if abstract:
            abstract = formater.update(acc, abstract)
        else:
            abstract = ''

        m = Method(acc, short_type, name, description, abstract)
        signatures.append(m)

    return signatures


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
        cur.execute("RENAME INTERPRO.PFAMN_MATCH_TMP TO INTERPRO.PFAMN_MATCH")

    return [Method(pfam_acc, sig_type=None) for pfam_acc in accessions]
