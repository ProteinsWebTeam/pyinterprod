import json
import os
import re

import cx_Oracle

from pyinterprod.utils.oracle import drop_table
from .common import Clan, Method

_TYPE = 'F'


def get_clans(uri: str) -> list[Clan]:
    con = cx_Oracle.connect(uri)
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


def parse_signatures(db_sources: dict) -> list[Method]:
    """
    Parse the names.tab file distributed with PANTHER releases

    :param db_sources:
    :return:
    """
    signatures = []
    prog = re.compile(r"(PTHR\d+)\.(SF\d+)?")
    with open(db_sources["sig_source"], "rt") as fh:
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


def update_go_terms(uri: str, root: str):
    con = cx_Oracle.connect(uri)
    cur = con.cursor()
    drop_table(cur, "INTERPRO.PANTHER2GO", purge=True)
    cur.execute(
        """
        CREATE TABLE INTERPRO.PANTHER2GO
        (
            METHOD_AC VARCHAR2(25) NOT NULL
                CONSTRAINT FK_PANTHER2GO
                REFERENCES INTERPRO.METHOD (METHOD_AC) ON DELETE CASCADE,
            AN_ID VARCHAR(10) NOT NULL,
            GO_ID VARCHAR2(10) NOT NULL,
            CONSTRAINT PK_PANTHER2GO
            PRIMARY KEY (METHOD_AC, AN_ID, GO_ID)
        ) NOLOGGING
        """
    )

    sql = """
        INSERT /*+ APPEND */ 
        INTO INTERPRO.PANTHER2GO
        VALUES (:1, :2, :3)
    """

    records = []
    for name in os.listdir(root):
        if name.endswith(".json"):
            with open(os.path.join(root, name), "rt") as fh:
                data = json.load(fh)

            for ancestral_node_id, (subfam_id, go_terms, _, _) in data.items():
                if subfam_id and go_terms:
                    for go_id in set(go_terms.split(",")):
                        records.append((subfam_id, ancestral_node_id, go_id))

                        if len(records) == 1000:
                            cur.executemany(sql, records)
                            con.commit()
                            records.clear()

    if records:
        cur.executemany(sql, records)
        con.commit()
        records.clear()

    cur.close()
    con.close()
