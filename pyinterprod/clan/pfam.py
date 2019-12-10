# -*- coding: utf-8 -*-

from typing import List

import MySQLdb
import MySQLdb.cursors

from pyinterprod.orautils import parse_url

DBCODE = 'H'


def get_clans(url) -> List[dict]:
    con = MySQLdb.connect(**parse_url(url))
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
            clan = clans[clan_acc] = {
                "accession": clan_acc,
                "name": row[1],
                "description": row[2],
                "members": []
            }

        clan["members"].append({
            "accession": row[4],
            "score": row[5] / row[3]
        })

    cur.close()
    con.close()

    return list(clans.values())
