# -*- coding: utf-8 -*-

import re
from typing import List

from cx_Oracle import Connection


def get_clans(con: Connection) -> List[dict]:
    cur = con.cursor()
    cur.execute(
        """
        SELECT METHOD_AC 
        FROM INTERPRO.METHOD 
        WHERE DBCODE = 'V'
        """
    )

    clans = {}
    for row in cur:
        member_acc, = row

        match = re.match("(PTHR\d+):SF\d+", member_acc)
        if not match:
            continue

        clan_acc = match.group(1)

        try:
            clan = clans[clan_acc]
        except KeyError:
            clan = clans[clan_acc] = {
                "accession": clan_acc,
                "name": None,
                "description": None,
                "members": []
            }

        clan["members"].append({
            "accession": member_acc,
            "score": 1
        })

    cur.close()
    return list(clans.values())
