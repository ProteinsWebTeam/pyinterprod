# -*- coding: utf-8 -*-

import cx_Oracle
import psycopg2
from psycopg2.extras import execute_values

from pyinterprod import logger
from pyinterprod.utils.kvdb import KVdb
from pyinterprod.utils.pg import url2dict


def import_similarity_comments(ora_url: str, pg_url: str):
    logger.info("populating")
    pg_con = psycopg2.connect(**url2dict(pg_url))
    with pg_con.cursor() as pg_cur:
        pg_cur.execute("TRUNCATE TABLE protein_similarity")

        ora_con = cx_Oracle.connect(ora_url)
        ora_cur = ora_con.cursor()
        ora_cur.execute(
            """
            SELECT E.ACCESSION, NVL(B.TEXT, SS.TEXT)
            FROM SPTR.DBENTRY@SWPREAD E
            INNER JOIN SPTR.COMMENT_BLOCK@SWPREAD B
              ON E.DBENTRY_ID = B.DBENTRY_ID
              AND B.COMMENT_TOPICS_ID = 34        -- SIMILARITY comments
            LEFT OUTER JOIN SPTR.COMMENT_STRUCTURE@SWPREAD S
              ON B.COMMENT_BLOCK_ID = S.COMMENT_BLOCK_ID
              AND S.CC_STRUCTURE_TYPE_ID = 1      -- TEXT structure
            LEFT OUTER JOIN SPTR.COMMENT_SUBSTRUCTURE@SWPREAD SS
              ON S.COMMENT_STRUCTURE_ID = SS.COMMENT_STRUCTURE_ID
            WHERE E.ENTRY_TYPE = 0                -- Swiss-Prot
              AND E.MERGE_STATUS != 'R'           -- not 'Redundant'
              AND E.DELETED = 'N'                 -- not deleted
              AND E.FIRST_PUBLIC IS NOT NULL      -- published
            """
        )

        sql = "INSERT INTO protein_similarity VALUES %s"
        execute_values(pg_cur, sql, ora_cur, page_size=1000)

        ora_cur.close()
        ora_con.close()

        pg_con.commit()
        pg_cur.execute("ANALYZE protein_similarity")

    pg_con.close()
    logger.info("complete")


def import_protein_names(ora_url: str, pg_url: str):
    logger.info("populating")
    pg_con = psycopg2.connect(**url2dict(pg_url))
    with pg_con.cursor() as pg_cur:
        pg_cur.execute("TRUNCATE TABLE protein_name")
        pg_cur.execute("TRUNCATE TABLE protein2name")

        ora_con = cx_Oracle.connect(ora_url)
        ora_cur = ora_con.cursor()
        ora_cur.execute(
            """
            SELECT ACCESSION, DESCR
            FROM (
              SELECT
                E.ACCESSION,
                D.DESCR,
                ROW_NUMBER() OVER (
                  PARTITION BY E.ACCESSION
                  ORDER BY CV.DESC_ID,    -- 1=RecName, 2=AltName, 3=SubName
                  CV.ORDER_IN,            -- Swiss-Prot manual order
                  D.DESCR                 -- TrEMBL alphabetic order
              ) R
              FROM SPTR.DBENTRY@SWPREAD E
              INNER JOIN SPTR.DBENTRY_2_DESC@SWPREAD D
                ON E.DBENTRY_ID = D.DBENTRY_ID
                AND D.DESC_ID IN (1,4,11,13,16,23,25,28,35)  --Full description section
              INNER JOIN SPTR.CV_DESC@SWPREAD CV
                ON D.DESC_ID = CV.DESC_ID
              WHERE E.ENTRY_TYPE IN (0, 1)          -- Swiss-Prot/TrEMBL
                AND E.MERGE_STATUS != 'R'           -- not 'Redundant'
                AND E.DELETED = 'N'                 -- not deleted
                AND E.FIRST_PUBLIC IS NOT NULL      -- published
            )
            WHERE R = 1                             -- one name per protein                            --
            """
        )

        names = {}
        values = []
        for protein_acc, text in ora_cur:
            try:
                name_id = names[text]
            except KeyError:
                name_id = names[text] = len(names) + 1

            values.append((protein_acc, name_id))
            if len(values) == 1000:
                execute_values(pg_cur, "INSERT INTO protein2name VALUES %s",
                               values, page_size=1000)
                values = []

        ora_cur.close()
        ora_con.close()

        if values:
            execute_values(pg_cur, "INSERT INTO protein2name VALUES %s",
                           values, page_size=1000)

        execute_values(pg_cur, "INSERT INTO protein_name VALUES %s",
                       ((name_id, text) for text, name_id in names.items()),
                       page_size=1000)

        pg_con.commit()
        pg_cur.execute("ANALYZE protein2name")
        pg_cur.execute("ANALYZE protein_name")

    pg_con.close()
    logger.info("complete")


def import_proteins(ora_url: str, pg_url: str):
    logger.info("populating")
    pg_con = psycopg2.connect(**url2dict(pg_url))
    with pg_con.cursor() as pg_cur:
        pg_cur.execute("TRUNCATE TABLE protein")

        ora_con = cx_Oracle.connect(ora_url)
        ora_cur = ora_con.cursor()
        ora_cur.execute(
            """
            SELECT PROTEIN_AC, NAME, LEN, TAX_ID, FRAGMENT, DBCODE
            FROM INTERPRO.PROTEIN                 --
            """
        )

        sql = "INSERT INTO protein VALUES %s"
        execute_values(pg_cur, sql, ((
            row[0],
            row[1],
            row[2],
            row[3],
            row[4] == 'Y',
            row[5] == 'S'
        ) for row in ora_cur), page_size=1000)

        ora_cur.close()
        ora_con.close()

        pg_con.commit()
        pg_cur.execute("ANALYZE protein")

    pg_con.close()
    logger.info("complete")


def export_names(url: str, database: str):
    logger.info("exporting protein names")
    con = psycopg2.connect(**url2dict(url))
    with con.cursor(name="names") as cur, KVdb(database, True) as names:
        cur.itersize = 1000000
        cur.execute("SELECT protein_acc, name_id FROM protein2name")
        for i, (key, value) in enumerate(cur):
            if not i % 1000000:
                names.sync()

            names[key] = value

    con.close()
    logger.info("complete")
