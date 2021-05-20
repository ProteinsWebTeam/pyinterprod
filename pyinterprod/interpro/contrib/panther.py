# -*- coding: utf-8 -*-

import re
from typing import List

import cx_Oracle

from .common import Clan, Method

_TYPE = 'F'


def get_clans(url: str) -> List[Clan]:
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.execute(
        """
        SELECT METHOD_AC
        FROM INTERPRO.METHOD 
        WHERE DBCODE = 'V'
        """
    )

    reg_subfam = re.compile(r"(PTHR\d+):SF\d+")
    clans = {}
    for member_acc, in cur:
        match = reg_subfam.match(member_acc)
        if not match:
            continue

        clan_acc = match.group(1)
        try:
            clan = clans[clan_acc]
        except KeyError:
            clan = clans[clan_acc] = Clan(clan_acc)

        clan.members.append({
            "accession": member_acc,
            "score": 1
        })

    return list(clans.values())


def parse_signatures(filepath: str) -> List[Method]:
    """
    Parse the names.tab file distributed with PANTHER releases

    :param filepath:
    :return:
    """
    signatures = []
    prog = re.compile(r"(PTHR\d+)\.(SF\d+)?")
    with open(filepath, "rt") as fh:
        for line in fh:
            cols = line.rstrip().split(sep='\t', maxsplit=1)
            try:
                accession, name = cols
            except ValueError:
                accession, = cols
                name = None

            fam_id, sub_id = prog.search(accession).groups()
            if sub_id:
                accession = f"{fam_id}:{sub_id}"
            else:
                accession = fam_id

            if name is None or name.upper().endswith("NOT NAMED"):
                descr = None
            else:
                descr = name

            m = Method(accession, _TYPE, description=descr)
            signatures.append(m)

    return signatures
