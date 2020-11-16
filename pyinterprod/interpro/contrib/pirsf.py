# -*- coding: utf-8 -*-

import re
from typing import List

from .model import Clan


def get_clans(file: str) -> List[Clan]:
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

    return [c for c in clans.values() if len(c.members) > 0]
