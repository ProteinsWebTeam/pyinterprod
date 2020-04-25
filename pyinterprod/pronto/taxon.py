# -*- coding: utf-8 -*-

import cx_Oracle
import psycopg2
from psycopg2.extras import execute_values

from pyinterprod import logger
from pyinterprod.utils.pg import url2dict


RANKS = {"superkingdom", "kingdom", "phylum", "class", "order", "family",
         "genus", "species"}


def iter_lineage(taxa: dict):
    for tax_id, (name, rank, left_num, right_num, parent_id) in taxa.items():
        if rank in RANKS:
            yield tax_id, tax_id, rank

        while parent_id:
            node = taxa[parent_id]
            rank = node[2]

            if rank in RANKS:
                yield tax_id, parent_id, rank

            parent_id = node[0]


def import_taxonomy(ora_url: str, pg_url: str):
    logger.info("loading taxonomy info")
    ora_con = cx_Oracle.connect(ora_url)
    ora_cur = ora_con.cursor()
    ora_cur.execute(
        """
        SELECT TAX_ID, SCIENTIFIC_NAME, RANK, LEFT_NUMBER, RIGHT_NUMBER, 
               PARENT_ID
        FROM INTERPRO.ETAXI
        """
    )
    taxa = {}
    for tax_id, name, rank, left_num, right_num, parent_id in ora_cur:
        if tax_id in (1, 131567):
            """
            Skip root and meta-superkingdom (131567) which contains:
                * Bacteria (2)
                * Archaea (2157)
                * Eukaryota (2759)
            """
            continue
        elif parent_id in (1, 131567):
            rank = "superkingdom"
            parent_id = None

        taxa[tax_id] = (name, rank, left_num, right_num, parent_id)

    ora_cur.close()
    ora_con.close()

    pg_con = psycopg2.connect(**url2dict(pg_url))
    with pg_con.cursor() as pg_cur:
        pg_cur.execute("TRUNCATE TABLE taxon")
        pg_cur.execute("TRUNCATE TABLE lineage")

        logger.info("populating: taxon")
        execute_values(pg_cur, "INSERT INTO taxon VALUES %s", (
            (tax_id, name, rank, left_num, right_num, parent_id)
            for tax_id, (name, rank, left_num, right_num, parent_id)
            in taxa.items()
        ), page_size=1000)

        logger.info("populating: lineage")
        execute_values(pg_cur, "INSERT INTO lineage VALUES %s",
                       iter_lineage(taxa),
                       page_size=1000)

        pg_con.commit()

        pg_cur.execute("ANALYZE taxon")
        pg_cur.execute("ANALYZE lineage")

    pg_con.close()
    logger.info("complete")
