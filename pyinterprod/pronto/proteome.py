import cx_Oracle
import psycopg2
from psycopg2.extras import execute_values

from pyinterprod import logger
from pyinterprod.utils.pg import url2dict


class ProteomeInterator:
    def __init__(self, url: str):
        self.url = url
        self.proteomes = {}

    def __iter__(self):
        con = cx_Oracle.connect(self.url)
        cur = con.cursor()
        cur.execute(
            """
            SELECT DISTINCT P.UPID, P.PROTEOME_NAME, P.PROTEOME_TAXID,
                            P2U.ACCESSION
            FROM SPTR.PROTEOME2UNIPROT P2U
            INNER JOIN SPTR.PROTEOME P ON P2U.PROTEOME_ID = P.PROTEOME_ID
            WHERE P2U.ENTRY_TYPE IN (0, 1)  -- Swiss-Prot, TrEMBL
              AND P2U.MERGE_STATUS != 'R'   -- Not 'Redundant'
              AND P.IS_REFERENCE = 1
            """
        )

        for upid, name, taxid, accession in cur:
            if upid not in self.proteomes:
                self.proteomes[upid] = (upid, name, taxid)

            yield upid, accession

        cur.close()
        con.close()


def import_proteomes(ora_url: str, pg_url: str):
    logger.info("inserting reference proteomes info")

    pg_con = psycopg2.connect(**url2dict(pg_url))
    with pg_con.cursor() as pg_cur:
        pg_cur.execute("DROP TABLE IF EXISTS proteome")
        pg_cur.execute(
            """
            CREATE TABLE proteome (
                id VARCHAR(20) NOT NULL
                    CONSTRAINT proteome_pkey PRIMARY KEY,
                name VARCHAR(200) NOT NULL,
                taxon_id INTEGER NOT NULL
            )
            """
        )

        pg_cur.execute("DROP TABLE IF EXISTS proteome2protein")
        pg_cur.execute(
            """
            CREATE TABLE proteome2protein (
                id VARCHAR(20) NOT NULL,            
                protein_acc VARCHAR(15) NOT NULL            
            )
            """
        )

        sql = "INSERT INTO proteome2protein VALUES %s"
        iterator = ProteomeInterator(ora_url)
        execute_values(pg_cur, sql, iterator, page_size=1000)

        sql = "INSERT INTO proteome VALUES %s"
        execute_values(pg_cur, sql, list(iterator.proteomes.values()),
                       page_size=1000)

        logger.info("indexing")
        pg_cur.execute(
            """
            CREATE INDEX proteome_name_idx
            ON proteome (name)
            """
        )
        pg_cur.execute(
            """
            CREATE INDEX proteome2protein_id
            ON proteome2protein (id)
            """
        )
        pg_con.commit()

    pg_con.close()
    logger.info("complete")
