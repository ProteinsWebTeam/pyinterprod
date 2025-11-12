from dataclasses import dataclass
from datetime import datetime

import oracledb
import psycopg

from pyinterprod.utils.oracle import get_partitions
from pyinterprod.utils.pg import url2dict


_NOT_IN_ISPRO = ["ELM", "Pfam-N"]


@dataclass
class Database:
    identifier: str
    name: str
    version: str
    date: datetime
    analysis_id: int
    has_site_matches: bool
    is_member_db: bool
    is_feature_db: bool


def get_databases(ora_uri: str, pg_uri: str, names: list[str],
                  expects_new: bool = False) -> dict[str, Database]:
    names = {name.lower(): name for name in names}

    con = psycopg.connect(**url2dict(pg_uri))
    with con.cursor() as cur:
        cur.execute(
            """
            SELECT a.id
            FROM iprscan.analysis a
            INNER JOIN iprscan.analysis_tables t ON a.name = t.name
            WHERE t.site_table IS NOT NULL
            """
        )
        with_sites = {row[0] for row in cur.fetchall()}
    con.close()

    con = oracledb.connect(ora_uri)
    cur = con.cursor()

    cur.execute(
        """
        SELECT DBCODE, COUNT(*)
        FROM INTERPRO.METHOD
        GROUP BY DBCODE
        """
    )
    cnt_signatures = dict(cur.fetchall())

    cur.execute(
        """
        SELECT DBCODE, COUNT(*)
        FROM INTERPRO.FEATURE_METHOD
        GROUP BY DBCODE
        """
    )
    cnt_features = dict(cur.fetchall())

    match_databases = set()
    for p in get_partitions(cur, "INTERPRO", "MATCH"):
        dbcode = p["value"][1:-1]  # 'X' -> X
        match_databases.add(dbcode)

    fmatch_databases = set()
    for p in get_partitions(cur, "INTERPRO", "FEATURE_MATCH"):
        dbcode = p["value"][1:-1]  # 'X' -> X
        fmatch_databases.add(dbcode)

    args = [':' + str(i+1) for i in range(len(names))]
    cur.execute(
        f"""
        SELECT D.DBCODE, LOWER(D.DBSHORT), D.DBNAME, V.VERSION, V.FILE_DATE,
                I2C.IPRSCAN_SIG_LIB_REL_ID, I2C.PREV_ID
        FROM INTERPRO.CV_DATABASE D
        LEFT OUTER JOIN INTERPRO.DB_VERSION V 
            ON D.DBCODE = V.DBCODE
        LEFT OUTER JOIN INTERPRO.IPRSCAN2DBCODE I2C 
            ON D.DBCODE = I2C.DBCODE
        WHERE LOWER(D.DBSHORT) IN ({','.join(args)})
        """,
        list(names.keys())
    )

    databases = {}
    missing = []
    not_ready = []
    for row in cur:
        db = Database(identifier=row[0],
                      name=row[2],
                      version=row[3],
                      date=row[4],
                      analysis_id=row[5],
                      has_site_matches=row[5] in with_sites,
                      is_member_db=(cnt_signatures.get(row[0], 0) > 0
                                    or row[0] in match_databases),
                      is_feature_db=(cnt_features.get(row[0], 0) > 0
                                     or row[0] in fmatch_databases))

        databases[row[1]] = db

        if row[3] is None or row[4] is None:
            missing.append(db.name)
        elif db.name in _NOT_IN_ISPRO:
            pass
        elif any(e is None for e in row[5:7]):
            missing.append(db.name)
        elif expects_new and db.analysis_id == row[6]:
            not_ready.append(db.name)

        try:
            del names[row[1]]
        except KeyError:
            pass

    cur.close()
    con.close()

    if names:
        names = sorted(names.values())
        raise RuntimeError(f"Unknown databases: {', '.join(names)}. "
                           f"Check CV_DATABASE.")
    elif missing:
        raise RuntimeError(f"Database(s) not initialised: "
                           f"{', '.join(missing)}. Run ipr-pre-memdb.")
    elif not_ready:
        raise RuntimeError(f"Database(s) outdated in IPRSCAN2DBCODE: "
                           f"{', '.join(not_ready)}. Run ipr-pre-memdb.")

    return databases


def update_database(ora_uri: str, pg_uri: str,
                    name: str, version: str, date: str,
                    by_name: bool = False, confirm: bool = True):
    con = oracledb.connect(ora_uri)
    cur = con.cursor()
    cur.execute(
        """
        SELECT D.DBCODE, D.DBNAME, V.VERSION, V.FILE_DATE, 
               I.IPRSCAN_SIG_LIB_REL_ID, I.PREV_ID
        FROM INTERPRO.CV_DATABASE D
        LEFT OUTER JOIN INTERPRO.DB_VERSION V ON D.DBCODE = V.DBCODE
        LEFT OUTER JOIN INTERPRO.IPRSCAN2DBCODE I ON D.DBCODE = I.DBCODE
        WHERE LOWER(D.DBSHORT) = LOWER(:1)
        """,
        [name]
    )
    row = cur.fetchone()

    if not row:
        cur.close()
        con.close()
        raise RuntimeError(f"No database matching '{name}'")

    dbcode, name, current_version, current_date, current_id, prev_id = row
    cur.close()
    con.close()

    # Find the 'active' analysis
    con = psycopg.connect(**url2dict(pg_uri))
    with con.cursor() as cur:
        if by_name or current_id is None:
            # Use name instead of analysis ID
            cur.execute(
                """
                SELECT id, name, version
                FROM iprscan.analysis
                WHERE LOWER(name) = LOWER(%s)
                  AND active = true
                ORDER BY id
                """,
                [name]
            )
        else:
            cur.execute(
                """
                SELECT id, name, version
                FROM iprscan.analysis
                WHERE name = (
                    SELECT name
                    FROM iprscan.analysis
                    WHERE id = %s
                )
                  AND active = true
                ORDER BY id
                """,
                [current_id]
            )

        rows = cur.fetchall()

    con.close()

    if name in _NOT_IN_ISPRO:
        id_to_use = current_id
    elif not rows:
        raise ValueError("Missing analysis")
    else:
        id_to_use = None
        for _id, _name, _version in rows:
            if _name.lower() == name.lower() and _version == version:
                id_to_use = _id
                break

        if id_to_use is None:
            raise ValueError(f"No active analysis for {name}-{version}")

    print(f"Updating {name}")
    if current_version and current_date:
        print(f"  Currently: {current_version} ({current_date:%Y-%m-%d})")
    else:
        print(f"  Currently: N/A")

    print(f"  Update to: {version} ({date})")

    if confirm and input("Do you want to continue? [y/N] ").lower() != 'y':
        print("Abort.")
        return

    con = oracledb.connect(ora_uri)
    cur = con.cursor()

    if (current_version, current_date) != (version, date):
        # Update DB_VERSION
        cur.execute(
            """
            SELECT COUNT(*) 
            FROM INTERPRO.DB_VERSION 
            WHERE DBCODE = :1
            """,
            [dbcode]
        )
        cnt, = cur.fetchone()
        if cnt:
            cur.execute(
                """
                UPDATE INTERPRO.DB_VERSION 
                SET VERSION = :1,
                    FILE_DATE = TO_DATE(:2, 'YYYY-MM-DD'),
                    LOAD_DATE = SYSDATE
                WHERE DBCODE = :3
                """,
                [version, date, dbcode]
            )
        else:
            cur.execute(
                """
                INSERT INTO INTERPRO.DB_VERSION (
                    DBCODE, VERSION, ENTRY_COUNT, FILE_DATE
                ) VALUES (:1, :2, 0, TO_DATE(:3, 'YYYY-MM-DD'))
                """,
                [dbcode, version, date]
            )

    if id_to_use != current_id:
        # Update IPRSCAN2DBCODE
        cur.execute(
            """
            SELECT COUNT(*) 
            FROM INTERPRO.IPRSCAN2DBCODE
            WHERE DBCODE = :1
            """,
            [dbcode]
        )
        cnt, = cur.fetchone()
        if cnt:
            cur.execute(
                """
                UPDATE INTERPRO.IPRSCAN2DBCODE
                SET IPRSCAN_SIG_LIB_REL_ID = :1
                WHERE DBCODE = :2
                """,
                [id_to_use, dbcode]
            )
        else:
            cur.execute(
                """
                SELECT CODE, DESCRIPTION
                FROM INTERPRO.CV_EVIDENCE
                """
            )

            print("Available evidences:")
            evidences = set()
            for code, description in cur:
                print(f"  {code:<10}{description}")
                evidences.add(code)

            evidence_to_use = None
            while evidence_to_use not in evidences:
                evidence_to_use = input("Enter evidence to use: ")

            cur.execute(
                """
                INSERT INTO INTERPRO.IPRSCAN2DBCODE (
                    IPRSCAN_SIG_LIB_REL_ID, DBCODE, EVIDENCE, PREV_ID
                ) VALUES (:1, :2, :3, 0)
                """,
                [id_to_use, dbcode, evidence_to_use]
            )

    con.commit()
    cur.close()
    con.close()
    print("Success.")
