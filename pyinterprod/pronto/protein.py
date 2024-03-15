import os
import re
import shutil
from tempfile import mkstemp

import oracledb
import psycopg

from pyinterprod import logger
from pyinterprod.utils.kvdb import KVdb
from pyinterprod.utils.pg import url2dict


def import_similarity_comments(swp_url: str, ipr_url: str):
    logger.info("populating")
    pg_con = psycopg.connect(**url2dict(ipr_url))
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

        ora_con = oracledb.connect(swp_url)
        ora_cur = ora_con.cursor()
        ora_cur.execute(
            """
            SELECT 
                DENSE_RANK() OVER (ORDER BY TEXT),
                TEXT,
                ACCESSION
            FROM (
                SELECT E.ACCESSION, NVL(B.TEXT, SS.TEXT) AS TEXT
                FROM SPTR.DBENTRY E
                INNER JOIN SPTR.COMMENT_BLOCK B
                  ON E.DBENTRY_ID = B.DBENTRY_ID
                  AND B.COMMENT_TOPICS_ID = 34        -- SIMILARITY comments
                LEFT OUTER JOIN SPTR.COMMENT_STRUCTURE S
                  ON B.COMMENT_BLOCK_ID = S.COMMENT_BLOCK_ID
                  AND S.CC_STRUCTURE_TYPE_ID = 1      -- TEXT structure
                LEFT OUTER JOIN SPTR.COMMENT_SUBSTRUCTURE SS
                  ON S.COMMENT_STRUCTURE_ID = SS.COMMENT_STRUCTURE_ID
                WHERE E.ENTRY_TYPE = 0                -- Swiss-Prot
                  AND E.MERGE_STATUS != 'R'           -- not 'Redundant'
                  AND E.DELETED = 'N'                 -- not deleted
                  AND E.FIRST_PUBLIC IS NOT NULL      -- published            
            )
            """
        )

        records = []
        sql = """
            INSERT INTO protein_similarity 
                (comment_id, comment_text, protein_acc) 
            VALUES (%s, %s, %s)
        """

        for rec in ora_cur:
            records.append(rec)
            if len(records) == 1000:
                pg_cur.executemany(sql, records)
                pg_con.commit()
                records.clear()

        if records:
            pg_cur.executemany(sql, records)
            pg_con.commit()
            records.clear()

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


def import_protein_names(swp_url: str, ipr_url: str, database: str,
                         tmpdir: str | None = None):
    os.makedirs(os.path.dirname(database), exist_ok=True)

    logger.info("populating protein2name")
    fd, tmp_database = mkstemp(dir=tmpdir)
    os.close(fd)
    os.remove(tmp_database)

    pg_con = psycopg.connect(**url2dict(ipr_url))
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

        ora_con = oracledb.connect(swp_url)
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
              FROM SPTR.DBENTRY E
              INNER JOIN SPTR.DBENTRY_2_DESC D
                ON E.DBENTRY_ID = D.DBENTRY_ID
                AND D.DESC_ID IN (1,4,11,13,16,23,25,28,35)  --Full description section
              INNER JOIN SPTR.CV_DESC CV
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

        sql = """
            COPY protein2name (protein_acc, name_id) 
            FROM STDIN
        """

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
                    with pg_cur.copy(sql) as copy:
                        for rec in values:
                            copy.write_row(rec)
                    values = []

                    if not i % 10000000:
                        logger.info(f"{i:>12,}")

        ora_cur.close()
        ora_con.close()
        logger.info(f"{i:>12,}")

        if values:
            with pg_cur.copy(sql) as copy:
                for rec in values:
                    copy.write_row(rec)

        logger.info("populating protein_name")

        sql = """
            COPY protein_name (name_id, text) 
            FROM STDIN
        """

        with pg_cur.copy(sql) as copy:
            for text, name_id in names.items():
                copy.write_row((name_id, text))

        logger.info("analyzing tables")
        pg_cur.execute("ANALYZE protein2name")
        pg_cur.execute("ANALYZE protein_name")
        pg_con.commit()

    pg_con.close()

    logger.info("copying database")
    shutil.copyfile(tmp_database, database)
    logger.info(f"disk usage: "
                f"{os.path.getsize(tmp_database) / 1024 ** 2:.0f} MB")
    os.remove(tmp_database)
    logger.info("complete")


def import_proteins(ora_url: str, pg_url: str):
    logger.info("populating")
    pg_con = psycopg.connect(**url2dict(pg_url))
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

        ora_con = oracledb.connect(ora_url)
        ora_cur = ora_con.cursor()
        ora_cur.execute(
            """
            SELECT PROTEIN_AC, NAME, LEN, TAX_ID, FRAGMENT, DBCODE
            FROM INTERPRO.PROTEIN
            """
        )

        sql = """
            COPY protein 
                (accession, identifier, length, taxon_id, is_fragment, 
                is_reviewed) 
            FROM STDIN
        """

        with pg_cur.copy(sql) as copy:
            for row in ora_cur:
                copy.write_row((row[0],
                                row[1],
                                row[2],
                                row[3],
                                row[4] == 'Y',
                                row[5] == 'S'
                                ))

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


def import_protein_pubmed(swp_url: str, pg_url: str):
    logger.info("populating protein2publication")
    pg_con = psycopg.connect(**url2dict(pg_url))
    with pg_con.cursor() as pg_cur:
        pg_cur.execute("DROP TABLE IF EXISTS protein2publication")
        pg_cur.execute(
            """
            CREATE TABLE protein2publication (
                PROTEIN_ACC VARCHAR(15) NOT NULL,
                PUBMED_ID INTEGER NOT NULL
            )
            """
        )

        ora_con = oracledb.connect(swp_url)
        ora_cur = ora_con.cursor()
        swp2pmid = set()

        ora_cur.execute(
            """
                SELECT E.ACCESSION, NVL(B.TEXT, SS.TEXT) AS TEXT
                FROM SPTR.DBENTRY E
                INNER JOIN SPTR.COMMENT_BLOCK B ON E.DBENTRY_ID = B.DBENTRY_ID
                    AND B.COMMENT_TOPICS_ID = 2
                INNER JOIN SPTR.CV_COMMENT_TOPICS CT 
                    ON B.COMMENT_TOPICS_ID = CT.COMMENT_TOPICS_ID
                LEFT OUTER JOIN SPTR.COMMENT_STRUCTURE S
                    ON B.COMMENT_BLOCK_ID = S.COMMENT_BLOCK_ID
                    AND S.CC_STRUCTURE_TYPE_ID = 1
                LEFT OUTER JOIN SPTR.COMMENT_SUBSTRUCTURE SS
                    ON S.COMMENT_STRUCTURE_ID = SS.COMMENT_STRUCTURE_ID
                WHERE E.ENTRY_TYPE = 0
                  AND E.MERGE_STATUS != 'R'           
                  AND E.DELETED = 'N'
            """
        )
        protein_text = {acc: text for acc, text in ora_cur.fetchall()}

        sql = """
            INSERT INTO protein2publication
                (protein_acc, pubmed_id)
            VALUES (%s, %s)
        """
        reg_pubmed = re.compile(r"(PubMed:(\d+))")
        for acc, text in protein_text.items():
            pmids = reg_pubmed.findall(text)
            if pmids:
                for pmid in pmids:
                    swp2pmid.add((acc, pmid[1]))
                    if len(swp2pmid) == 1000:
                        pg_cur.executemany(sql, swp2pmid)
                        pg_con.commit()
                        swp2pmid.clear()
        if swp2pmid:
            pg_cur.executemany(sql, swp2pmid)
            pg_con.commit()
            swp2pmid.clear()

        pg_cur.execute(
            """
            CREATE INDEX protein2publication_idx
            ON protein2publication (protein_acc)
            """
        )

    pg_con.commit()
    pg_con.close()
    ora_con.close()
    logger.info("done")


if __name__ == '__main__':
    import configparser

    config = configparser.ConfigParser()
    config.read("/Users/lcf/PycharmProjects/pyinterprod/test_data/pyinterprod.config")

    import_protein_pubmed(
        config["oracle"]["unpr-swpread"],
        config["postgresql"]["pronto"]
    )
