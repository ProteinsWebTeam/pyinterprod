import cx_Oracle
import psycopg2
from psycopg2.extras import execute_values

from pyinterprod import logger
from pyinterprod.utils.pg import url2dict


def import_proteomes(ora_url: str, pg_url: str):
    logger.info("inserting reference proteomes info")

    pg_con = psycopg2.connect(**url2dict(pg_url))
    with pg_con.cursor() as pg_cur:
        pg_cur.execute("DROP TABLE IF EXISTS proteome2protein")
        pg_cur.execute(
            """
            CREATE TABLE proteome2protein (
                id VARCHAR(20) NOT NULL
                    CONSTRAINT proteome2protein_pkey PRIMARY KEY,
                name VARCHAR(200) NOT NULL,
                taxon_id INTEGER NOT NULL,            
                protein_acc VARCHAR(15) NOT NULL            
            )
            """
        )

        ora_con = cx_Oracle.connect(ora_url)
        ora_cur = ora_con.cursor()
        ora_cur.execute(
            """
            SELECT DISTINCT P.UPID, P.PROTEOME_NAME, P.PROTEOME_TAXID,
                            P2U.ACCESSION
            FROM SPTR.PROTEOME2UNIPROT P2U
            INNER JOIN SPTR.PROTEOME P ON P2U.PROTEOME_ID = P.PROTEOME_ID
            WHERE P2U.ENTRY_TYPE IN (1, 2)
              AND P2U.MERGE_STATUS != 'R'
              AND P.IS_REFERENCE = 1
            """
        )

        sql = "INSERT INTO proteome2protein VALUES %s"
        execute_values(pg_cur, sql, ora_cur, page_size=1000)

        ora_cur.close()
        ora_con.close()

        pg_cur.execute(
            """
            CREATE INDEX proteome2protein_taxon_idx
            ON proteome2protein (taxon_id)
            """
        )
        pg_con.commit()

    pg_con.close()
    logger.info("complete")
