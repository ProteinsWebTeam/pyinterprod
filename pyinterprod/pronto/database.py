import oracledb
import psycopg

from pyinterprod import logger
from pyinterprod.utils.pg import url2dict


DATABASES = [
    'B',    # SFLD
    'F',    # PRINTS
    'H',    # Pfam
    'I',    # InterPro
    'J',    # CDD
    'M',    # PROSITE profiles
    'N',    # NCBIfam
    'P',    # PROSITE patterns
    'Q',    # HAMAP
    'R',    # SMART
    'U',    # PIRSF
    'V',    # PANTHER
    'X',    # CATH-Gene3D
    'Y',    # SUPERFAMILY
    'a',    # AntiFam
    'd',    # Pfam-N
    'f',    # FunFam
    'g',    # MobiDB
    'l',    # ELM
    'u',    # UniProtKB
]


def import_databases(ora_url: str, pg_url: str):
    logger.info("populating")
    pg_con = psycopg.connect(**url2dict(pg_url))
    with pg_con.cursor() as pg_cur:
        pg_cur.execute("DROP TABLE IF EXISTS database")
        pg_cur.execute(
            """
            CREATE TABLE database (
                id SERIAL NOT NULL 
                    CONSTRAINT database_pkey PRIMARY KEY,
                name VARCHAR(50) NOT NULL,
                name_long VARCHAR(50) NOT NULL,
                version VARCHAR(50) NOT NULL,
                updated DATE NOT NULL
            )
            """
        )

        ora_con = oracledb.connect(ora_url)
        ora_cur = ora_con.cursor()
        ora_cur.execute(
            """
            UPDATE INTERPRO.PRONTO_STATES
            SET ACTIVE = 'Y'
            WHERE NAME = 'UPDATING'
            """
        )
        ora_con.commit()

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

        pg_cur.execute(
            """
            CREATE UNIQUE INDEX database_name_idx
            ON database (name)
            """
        )
        pg_con.commit()

    pg_con.close()
    logger.info("complete")


def set_ready(ora_url: str):
    logger.info("updating status")
    ora_con = oracledb.connect(ora_url)
    ora_cur = ora_con.cursor()
    ora_cur.execute(
        """
        UPDATE INTERPRO.PRONTO_STATES
        SET ACTIVE = 'N'
        WHERE NAME = 'UPDATING'
        """
    )
    ora_con.commit()
    ora_cur.close()
    ora_con.close()
    logger.info("complete")
