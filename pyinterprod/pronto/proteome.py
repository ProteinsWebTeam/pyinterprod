import oracledb
import psycopg

from pyinterprod import logger
from pyinterprod.utils.pg import url2dict


class ProteomeIterator:
    def __init__(self, url: str):
        self.url = url
        self.proteomes = {}

    def __iter__(self):
        con = oracledb.connect(self.url)
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
            try:
                self.proteomes[upid][3] += 1
            except KeyError:
                self.proteomes[upid] = [upid, name, taxid, 1]

            yield upid, accession

        cur.close()
        con.close()


def import_proteomes(ora_url: str, pg_url: str):
    logger.info("inserting reference proteomes info")

    pg_con = psycopg.connect(**url2dict(pg_url))
    with pg_con.cursor() as pg_cur:
        pg_cur.execute("DROP TABLE IF EXISTS proteome")
        pg_cur.execute(
            """
            CREATE TABLE proteome (
                id VARCHAR(20) NOT NULL
                    CONSTRAINT proteome_pkey PRIMARY KEY,
                name VARCHAR(200) NOT NULL,
                taxon_id INTEGER NOT NULL,
                num_proteins INTEGER NOT NULL
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

        records = []
        iterator = ProteomeIterator(ora_url)

        sql = """
              INSERT INTO proteome2protein (id, protein_acc) 
              VALUES (%s, %s)
              """

        for rec in iterator:
            records.append(rec)
            if len(records) == 1000:
                pg_cur.executemany(sql, records)
                pg_con.commit()
                records.clear()
        if records:
            pg_cur.executemany(sql, records)
            pg_con.commit()
            records.clear()

        sql = """
              INSERT INTO proteome (id, name, taxon_id, num_proteins) 
              VALUES (%s, %s, %s, %s)
              """

        for rec in list(iterator.proteomes.values()):
            records.append(rec)
            if len(records) == 1000:
                pg_cur.executemany(sql, records)
                pg_con.commit()
                records.clear()
        if records:
            pg_cur.executemany(sql, records)
            pg_con.commit()
            records.clear()

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
