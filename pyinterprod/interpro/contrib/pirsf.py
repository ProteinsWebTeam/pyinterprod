import re

from oracledb import Cursor

from .common import Clan, Method


def get_clans(file: str) -> list[Clan]:
    """
    Get the PIRSF clans and their members.
    The file is available on PIRSF's FTP:
        ftp://ftp.pir.georgetown.edu/databases/pirsf/pirsfinfo.dat

    :param file: pirsfinfo.dat file
    :return: list of clans
    """
    clans = {}
    with open(file, "rt") as fh:
        for line in fh:
            if line[0] != '>':
                continue

            m = re.search(r">(PIRSF\d+)\s*\(.+?\)\s*(.+)?", line.rstrip())
            member_acc, string = m.groups()
            if not string:
                continue

            m = re.search(r"\[Parent=(PIRSF\d+)\]", string)
            if not m:
                continue

            name = string[:m.start(0)].rstrip()
            clan_acc = m.group(1)

            try:
                clan = clans[clan_acc]
            except KeyError:
                clan = clans[clan_acc] = Clan(clan_acc, name)

            clan.members.append({
                "accession": member_acc,
                "score": 1
            })

    return list(clans.values())


def get_signatures(cur: Cursor) -> list[Method]:
    """
    We don't expect to ever do an actual PIRSF update,
    as PIR hasn't created new families in ages.
    If we need to do a data (matches) updates, we simply get all existing
    signatures.
    """
    cur.execute(
        """
        SELECT METHOD_AC, NAME, DESCRIPTION, SIG_TYPE, ABSTRACT, ABSTRACT_LONG
        FROM INTERPRO.METHOD
        WHERE DBCODE = 'U'
        """
    )

    methods = []
    for row in cur:
        m = Method(
            accession=row[0],
            sig_type=row[3],
            name=row[1],
            description=row[2],
            abstract=row[5].read() if row[5] else row[4]
        )
        methods.append(m)

    return methods
