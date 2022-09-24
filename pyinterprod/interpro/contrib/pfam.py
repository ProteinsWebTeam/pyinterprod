import re

import MySQLdb
import MySQLdb.cursors

from .common import Clan, Method
from pyinterprod import logger
from pyinterprod.utils.pg import url2dict


_TYPES = {
    "Domain": 'D',
    "Family": 'F',
    "Repeat": 'R',
    "Coiled-coil": 'U',
    "Disordered": 'U',
    "Motif": 'C'
}


def connect_mysql(url: str) -> MySQLdb.Connection:
    obj = url2dict(url)
    # MySQLdb and psycopg2 use different keyword params for password/database
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
            "score": row[5] / row[3]
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
            uniprot_acc, pfam_acc, start, end = line.rstrip().split("\t")
            if re.fullmatch(r"PF\d+", pfam_acc):
                yield uniprot_acc, pfam_acc, int(start), int(end)


def get_protenn_entries(file: str) -> list[Method]:
    accessions = set()

    for _, pfam_acc, _, _ in iter_protenn_matches(file):
        pfam_acc.add(pfam_acc)

    return [Method(pfam_acc, sig_type=None) for pfam_acc in accessions]
