# -*- coding: utf-8 -*-

import cx_Oracle
import psycopg2

from pyinterprod import logger
from pyinterprod.utils.pg import url2dict


DATABASES = [
    'B',    # SFLD
    'F',    # PRINTS
    'H',    # Pfam
    'I',    # InterPro
    'J',    # CDD
    'M',    # PROSITE profiles
    'N',    # TIGRFAMs
    'P',    # PROSITE patterns
    'Q',    # HAMAP
    'R',    # SMART
    'U',    # PIRSF
    'V',    # PANTHER
    'X',    # CATH-Gene3D
    'Y',    # SUPERFAMILY
    'g',    # MobiDB Lite
    'u',    # UniProtKB
]


def import_databases(ora_url: str, pg_url: str):
    logger.info("populating")
    pg_con = psycopg2.connect(**url2dict(pg_url))
    with pg_con.cursor() as pg_cur:
        pg_cur.execute("TRUNCATE TABLE database RESTART IDENTITY")

        ora_con = cx_Oracle.connect(ora_url)
        ora_cur = ora_con.cursor()
        ora_cur.execute(
            f"""
            SELECT D.DBCODE, LOWER(D.DBSHORT), D.DBNAME, V.VERSION, V.FILE_DATE
            FROM INTERPRO.CV_DATABASE D
            LEFT OUTER JOIN INTERPRO.DB_VERSION V
              ON D.DBCODE = V.DBCODE
            """
        )

        for row in ora_cur:
            if row[0] in DATABASES:
                pg_cur.execute(
                    """
                    INSERT INTO database (name, name_long, version, updated)
                    VALUES (%s, %s, %s, %s)
                    """, row[1:]
                )

        ora_cur.close()
        ora_con.close()

        pg_con.commit()
        pg_cur.execute("ANALYZE database")

    pg_con.close()
    logger.info("complete")
