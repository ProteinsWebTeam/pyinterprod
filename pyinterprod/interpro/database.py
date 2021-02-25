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
    has_site_matches: bool


def get_databases(url: str, names: Sequence[str], expects_new: bool = False) -> Dict[str, Database]:
    con = cx_Oracle.connect(url)
    cur = con.cursor()

    args = [':' + str(i+1) for i in range(len(names))]
    cur.execute(
        f"""
        SELECT D.DBCODE, LOWER(D.DBSHORT), D.DBNAME, V.VERSION, V.FILE_DATE,
                I2C.IPRSCAN_SIG_LIB_REL_ID, I2C.PREV_ID, T.SITE_TABLE
        FROM INTERPRO.CV_DATABASE D
        INNER JOIN INTERPRO.DB_VERSION V 
            ON D.DBCODE = V.DBCODE
        INNER JOIN INTERPRO.IPRSCAN2DBCODE I2C 
            ON D.DBCODE = I2C.DBCODE
        INNER JOIN IPRSCAN.IPM_ANALYSIS@ISPRO A 
            ON I2C.IPRSCAN_SIG_LIB_REL_ID = A.ANALYSIS_ID
        INNER JOIN IPRSCAN.IPM_ANALYSIS_MATCH_TABLE@ISPRO T 
            ON A.ANALYSIS_MATCH_TABLE_ID = T.ID
        WHERE LOWER(D.DBSHORT) IN ({','.join(args)})
        """, tuple(map(str.lower, names))
    )

    databases = {}
    not_ready = []
    for row in cur:
        db = Database(identifier=row[0],
                      name=row[2],
                      version=row[3],
                      date=row[4],
                      analysis_id=row[5],
                      has_site_matches=row[7] is not None)

        databases[row[1]] = db

        if expects_new and db.analysis_id == row[6]:
            not_ready.append(db.name)

    cur.close()
    con.close()

    unknown = set(names) - set(databases.keys())
    if unknown:
        raise RuntimeError(f"Unknown databases: {', '.join(unknown)}")
    elif not_ready:
        raise RuntimeError(f"Database(s) not ready: "
                           f"{', '.join(not_ready)}")

    return databases


def update_database(url: str, name: str, version: str, date: str,
                    confirm: bool = True):
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.execute(
        """
        SELECT D.DBCODE, D.DBNAME, V.VERSION, V.FILE_DATE, 
               I.IPRSCAN_SIG_LIB_REL_ID, I.PREV_ID
        FROM INTERPRO.CV_DATABASE D
        INNER JOIN INTERPRO.DB_VERSION V ON D.DBCODE = V.DBCODE
        INNER JOIN INTERPRO.IPRSCAN2DBCODE I ON D.DBCODE = I.DBCODE
        WHERE LOWER(D.DBSHORT) = LOWER(:1)
        """, (name,)
    )
    row = cur.fetchone()

    if not row:
        cur.close()
        con.close()
        raise RuntimeError(f"No database matching '{name}'")

    dbcode, name, current_version, current_date, current_id, prev_id = row

    # Find the 'active' analysis in ISPRO
    cur.execute(
        """
        SELECT ANALYSIS_ID
        FROM IPM_ANALYSIS@ISPRO
        WHERE ANALYSIS_MATCH_TABLE_ID = (
            SELECT ANALYSIS_MATCH_TABLE_ID
            FROM IPM_ANALYSIS@ISPRO
            WHERE ANALYSIS_ID = :1
        )
        AND ACTIVE = 1
        """, (current_id,)
    )
    row = cur.fetchone()
    active_id, = row

    print(f"Updating {name}")
    print(f"  Currently: {current_version} ({current_date:%Y-%m-%d})")
    print(f"  Update to: {version} ({date})")

    if confirm and input("Do you want to continue? [y/N] ").lower() != 'y':
        cur.close()
        con.close()
        print("Abort.")
        return

    try:
        cur.execute(
            """
                UPDATE INTERPRO.DB_VERSION 
                SET VERSION = :1,
                    FILE_DATE = TO_DATE(:2, 'YYYY-MM-DD'),
                    LOAD_DATE = SYSDATE
                WHERE DBCODE = :3
            """, (version, date, dbcode)
        )

        if active_id != current_id:
            # Update IPRSCAN2DBCODE
            cur.execute(
                """
                UPDATE INTERPRO.IPRSCAN2DBCODE
                SET IPRSCAN_SIG_LIB_REL_ID = :1
                WHERE DBCODE = :2 AND IPRSCAN_SIG_LIB_REL_ID = :3
                """, (active_id, dbcode, current_id)
            )
    except cx_Oracle.Error as exc:
        con.rollback()
        print(f"{type(exc)}: {exc}")
    else:
        con.commit()
        print("Success.")
    finally:
        cur.close()
        con.close()
