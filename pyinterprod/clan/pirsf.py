# -*- coding: utf-8 -*-

import os
import re
from typing import List

from .utils import download


def get_clans() -> List[dict]:
    filepath = download("ftp://ftp.pir.georgetown.edu/databases/pirsf/pirsfinfo.dat")

    clans = {}
    with open(filepath, "rt") as fh:
        for line in fh:
            if line[0] != '>':
                continue

            m = re.search(">(PIRSF\d+)\s*\(.+?\)\s*(.+)?", line.rstrip())
            member_acc, string = m.groups()
            if not string:
                continue

            m = re.search("\[Parent=(PIRSF\d+)\]", string)
            if not m:
                continue

            name = string[:m.start(0)].rstrip()
            clan_acc = m.group(1)

            try:
                clan = clans[clan_acc]
            except KeyError:
                clan = clans[clan_acc] = {
                    "accession": clan_acc,
                    "name": name,
                    "description": None,
                    "members": []
                }

            clan["members"].append({
                "accession": member_acc,
                "score": 1
            })

    os.remove(filepath)
    return list(clans.values())
