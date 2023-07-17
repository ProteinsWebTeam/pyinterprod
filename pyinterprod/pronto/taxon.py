# -*- coding: utf-8 -*-

import json

import oracledb
import psycopg

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
            rank = node[1]

            if rank in RANKS:
                yield tax_id, parent_id, rank

            parent_id = node[4]


def get_lineage(taxa: dict, tax_id: int):
    path = [tax_id]
    node_id = taxa[tax_id][4]

    while node_id:
        node = taxa[node_id]
        rank = node[1]
        if rank in RANKS:
            path.append(node_id)

        node_id = node[4]

    return path[::-1]


def import_taxonomy(ora_url: str, pg_url: str):
    logger.info("loading taxonomy info")
    ora_con = oracledb.connect(ora_url)
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

    pg_con = psycopg.connect(**url2dict(pg_url))
    with pg_con.cursor() as pg_cur:
        pg_cur.execute("DROP TABLE IF EXISTS taxon")
        pg_cur.execute("DROP TABLE IF EXISTS lineage")
        pg_cur.execute(
            """
            CREATE TABLE taxon (
                id INTEGER NOT NULL 
                    CONSTRAINT taxon_id_pkey PRIMARY KEY,
                name VARCHAR(255) NOT NULL,
                rank VARCHAR(50) NOT NULL,
                left_number INTEGER NOT NULL,
                right_number INTEGER NOT NULL,
                parent_id INTEGER,
                lineage TEXT NOT NULL
            )
            """
        )
        pg_cur.execute(
            """
            CREATE TABLE lineage (
                child_id INTEGER NOT NULL,
                parent_id INTEGER NOT NULL,
                parent_rank VARCHAR(255) NOT NULL
            )
            """
        )

        logger.info("populating: taxon")
        records = []
        sql = """
              INSERT INTO taxon (id, name, rank, left_number, right_number, parent_id, lineage)
              VALUES (%s, %s, %s, %s, %s, %s, %s)
              """

        for tax_id, value in taxa.items():
            (name, rank, left_num, right_num, parent_id) = value
            json_lineage = json.dumps(get_lineage(taxa, tax_id))
            records.append((tax_id, name, rank, left_num,
                            right_num, parent_id, json_lineage))

            if len(records) == 1000:
                pg_cur.executemany(sql, records)
                pg_con.commit()
                records.clear()

        if records:
            pg_cur.executemany(sql, records)
            pg_con.commit()
            records.clear()

        pg_cur.execute(
            """
            CREATE INDEX taxon_left_number_idx
            ON taxon (left_number)
            """
        )

        logger.info("populating: lineage")
        sql = """
              INSERT INTO lineage (child_id, parent_id, parent_rank) 
              VALUES (%s, %s, %s)
              """

        for rec in iter_lineage(taxa):
            records.append(rec)

            if len(records) == 1000:
                pg_cur.executemany(sql, records)
                pg_con.commit()
                records.clear()

        if records:
            pg_cur.executemany(sql, records)
            pg_con.commit()
            records.clear()

        pg_cur.execute(
            """
            CREATE UNIQUE INDEX lineage_child_parent_uidx
            ON lineage (child_id, parent_id)
            """
        )
        pg_cur.execute(
            """
            CREATE INDEX lineage_child_rank_idx
            ON lineage (child_id, parent_rank)
            """
        )

        pg_con.commit()

    pg_con.close()
    logger.info("complete")
