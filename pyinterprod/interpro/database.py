from dataclasses import dataclass
from datetime import datetime

import cx_Oracle

from .match import MATCH_PARTITIONS, FEATURE_MATCH_PARTITIONS


_NOT_IN_ISPRO = ["Pfam-N"]


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


def get_databases(uri: str, names: list[str],
                  expects_new: bool = False) -> dict[str, Database]:
    names = {name.lower(): name for name in names}

    con = cx_Oracle.connect(uri)
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

    args = [':' + str(i+1) for i in range(len(names))]
    cur.execute(
        f"""
        SELECT D.DBCODE, LOWER(D.DBSHORT), D.DBNAME, V.VERSION, V.FILE_DATE,
                I2C.IPRSCAN_SIG_LIB_REL_ID, I2C.PREV_ID, T.SITE_TABLE
        FROM INTERPRO.CV_DATABASE D
        LEFT OUTER JOIN INTERPRO.DB_VERSION V 
            ON D.DBCODE = V.DBCODE
        LEFT OUTER JOIN INTERPRO.IPRSCAN2DBCODE I2C 
            ON D.DBCODE = I2C.DBCODE
        LEFT OUTER JOIN IPRSCAN.ANALYSIS@ISPRO A 
            ON I2C.IPRSCAN_SIG_LIB_REL_ID = A.ID
        LEFT OUTER JOIN IPRSCAN.ANALYSIS_TABLES@ISPRO T 
            ON A.NAME = T.NAME
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
                      has_site_matches=row[7] is not None,
                      is_member_db=(cnt_signatures.get(row[0], 0) > 0
                                    or row[0] in MATCH_PARTITIONS),
                      is_feature_db=(cnt_features.get(row[0], 0) > 0
                                     or row[0] in FEATURE_MATCH_PARTITIONS))

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


def update_database(uri: str, name: str, version: str, date: str,
                    confirm: bool = True):
    con = cx_Oracle.connect(uri)
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

    # Find the 'active' analysis in ISPRO
    if current_id is not None:
        cur.execute(
            """
            SELECT ID, NAME, VERSION
            FROM IPRSCAN.ANALYSIS@ISPRO
            WHERE NAME = (
                SELECT NAME
                FROM IPRSCAN.ANALYSIS@ISPRO
                WHERE ID = :1
            )
            AND ACTIVE = 'Y'
            ORDER BY ID
            """,
            [current_id]
        )
    else:
        # Not in IPRSCAN2DBCODE: fallback to name
        cur.execute(
            """
            SELECT ID, NAME, VERSION
            FROM IPRSCAN.ANALYSIS@ISPRO
            WHERE LOWER(NAME) = LOWER(:1)
            AND ACTIVE = 'Y'
            ORDER BY ID
            """,
            [name]
        )

    rows = cur.fetchall()
    cur.close()
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

    con = cx_Oracle.connect(uri)
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
