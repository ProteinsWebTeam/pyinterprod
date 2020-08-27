# -*- coding: utf-8 -*-

from dataclasses import dataclass
from datetime import datetime
from typing import Dict, Sequence

import cx_Oracle


@dataclass
class Database:
    identifier: str
    name: str
    version: str
    date: datetime
    analysis_id: int


def get_databases(url: str, names: Sequence[str]) -> Dict[str, Database]:
    names = set(names)
    con = cx_Oracle.connect(url)
    cur = con.cursor()

    args = [':' + str(i+1) for i in range(len(names))]
    cur.execute(
        f"""
        SELECT D.DBCODE, LOWER(D.DBSHORT), D.DBNAME, V.VERSION, V.FILE_DATE, 
               I2C.IPRSCAN_SIG_LIB_REL_ID
        FROM INTERPRO.CV_DATABASE D
        LEFT OUTER JOIN INTERPRO.DB_VERSION V ON D.DBCODE = V.DBCODE
        LEFT OUTER JOIN INTERPRO.IPRSCAN2DBCODE I2C ON D.DBCODE = I2C.DBCODE 
        WHERE LOWER(D.DBSHORT) IN ({','.join(args)})
        """, tuple(map(str.lower, names))
    )

    databases = {}
    for code, key, name, version, date, analysis_id in cur:
        databases[key] = Database(code, name, version, date, analysis_id)

    cur.close()
    con.close()

    for db in databases.values():
        if db.analysis_id is None:
            raise RuntimeError(f"{db.name} not in INTERPRO.IPRSCAN2DBCODE")

    return databases
