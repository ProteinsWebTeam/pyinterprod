# -*- coding: utf-8 -*-

import os
import shutil
from tempfile import mkstemp
from typing import Optional

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
        pg_cur.execute("DROP TABLE IF EXISTS protein_similarity")
        pg_cur.execute(
            """
            CREATE TABLE protein_similarity (
                comment_id INTEGER NOT NULL,
                comment_text TEXT NOT NULL,
                protein_acc VARCHAR(15) NOT NULL
            )
            """
        )

        ora_con = cx_Oracle.connect(ora_url)
        ora_cur = ora_con.cursor()
        ora_cur.execute(
            """
            SELECT 
                DENSE_RANK() OVER (ORDER BY TEXT),
                TEXT,
                ACCESSION
            FROM (
                SELECT E.ACCESSION, NVL(B.TEXT, SS.TEXT) AS TEXT
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
            )
            """
        )

        sql = "INSERT INTO protein_similarity VALUES %s"
        execute_values(pg_cur, sql, ora_cur, page_size=1000)

        ora_cur.close()
        ora_con.close()

        pg_cur.execute(
            """
            CREATE INDEX protein_similarity_comment_idx
            ON protein_similarity (comment_id)
            """
        )
        pg_cur.execute(
            """
            CREATE INDEX protein_similarity_protein_idx
            ON protein_similarity (protein_acc)
            """
        )
        pg_con.commit()

    pg_con.close()
    logger.info("complete")


def import_protein_names(ora_url: str, pg_url: str, database: str,
                         tmpdir: Optional[str] = None):
    logger.info("populating protein2name")
    fd, tmp_database = mkstemp(dir=tmpdir)
    os.close(fd)
    os.remove(tmp_database)

    pg_con = psycopg2.connect(**url2dict(pg_url))
    with pg_con.cursor() as pg_cur:
        pg_cur.execute("DROP TABLE IF EXISTS protein_name")
        pg_cur.execute("DROP TABLE IF EXISTS protein2name")
        pg_cur.execute(
            """
            CREATE TABLE protein_name (
                name_id INTEGER NOT NULL 
                    CONSTRAINT protein_name_pkey PRIMARY KEY,
                text TEXT NOT NULL
            )
            """
        )
        pg_cur.execute(
            """
            CREATE TABLE protein2name (
                protein_acc VARCHAR(15) NOT NULL 
                    CONSTRAINT protein2name_pkey PRIMARY KEY,
                name_id INTEGER NOT NULL
            )
            """
        )

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
            WHERE R = 1                             -- one name per protein
            """
        )

        names = {}
        values = []
        i = 0
        with KVdb(tmp_database, True) as namesdb:
            for protein_acc, text in ora_cur:
                try:
                    name_id = names[text]
                except KeyError:
                    name_id = names[text] = len(names) + 1

                values.append((protein_acc, name_id))
                namesdb[protein_acc] = name_id
                i += 1

                if not i % 100000:
                    namesdb.sync()
                    execute_values(
                        cur=pg_cur,
                        sql="INSERT INTO protein2name VALUES %s",
                        argslist=values,
                        page_size=1000
                    )
                    values = []

                    if not i % 10000000:
                        logger.info(f"{i:>12,}")

        ora_cur.close()
        ora_con.close()
        logger.info(f"{i:>12,}")

        if values:
            execute_values(
                cur=pg_cur,
                sql="INSERT INTO protein2name VALUES %s",
                argslist=values,
                page_size=1000
            )

        logger.info("populating protein_name")
        execute_values(
            cur=pg_cur,
            sql="INSERT INTO protein_name VALUES %s",
            argslist=((name_id, text) for text, name_id in names.items()),
            page_size=1000
        )

        logger.info("analyzing tables")
        pg_cur.execute("ANALYZE protein2name")
        pg_cur.execute("ANALYZE protein_name")
        pg_con.commit()

    pg_con.close()

    logger.info("copying database")
    shutil.copyfile(tmp_database, database)
    logger.info(f"disk usage: {os.path.getsize(tmp_database)/1024**2:.0f} MB")
    os.remove(tmp_database)
    logger.info("complete")


def import_proteins(ora_url: str, pg_url: str):
    logger.info("populating")
    pg_con = psycopg2.connect(**url2dict(pg_url))
    with pg_con.cursor() as pg_cur:
        pg_cur.execute("DROP TABLE IF EXISTS protein")
        pg_cur.execute(
            """
            CREATE TABLE protein (
                accession VARCHAR(15) NOT NULL
                    CONSTRAINT protein_pkey PRIMARY KEY,
                identifier VARCHAR(16) NOT NULL,
                length INTEGER NOT NULL,
                taxon_id INTEGER NOT NULL,
                is_fragment BOOLEAN NOT NULL,
                is_reviewed BOOLEAN NOT NULL
            )
            """
        )

        ora_con = cx_Oracle.connect(ora_url)
        ora_cur = ora_con.cursor()
        ora_cur.execute(
            """
            SELECT PROTEIN_AC, NAME, LEN, TAX_ID, FRAGMENT, DBCODE
            FROM INTERPRO.PROTEIN
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

        pg_cur.execute(
            """
            CREATE UNIQUE INDEX protein_identifier_uidx
            ON protein (identifier)
            """
        )
        pg_cur.execute(
            """
            CREATE INDEX protein_reviewed_idx
            ON protein (is_reviewed)
            """
        )
        pg_cur.execute(
            """
            CREATE INDEX protein_taxon_idx 
            ON protein (taxon_id)
            """
        )
        pg_con.commit()

    pg_con.close()
    logger.info("complete")
