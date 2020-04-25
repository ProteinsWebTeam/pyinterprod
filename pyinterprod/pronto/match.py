# -*- coding: utf-8 -*-

import os
import pickle
from tempfile import mkstemp
from typing import Dict, Optional

import cx_Oracle
import psycopg2

from pyinterprod import logger
from pyinterprod.utils.kvdb import KVdb
from pyinterprod.utils.pg import CsvIO, drop_index, url2dict
from .protein import export_names


def _iter_matches(url: str, databases: Dict[str, int]):
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.execute(
        """
        SELECT M.PROTEIN_AC, M.METHOD_AC, D.DBNAME,
               M.POS_FROM, M.POS_TO, M.FRAGMENTS
        FROM INTERPRO.MATCH M
        INNER JOIN INTERPRO.CV_DATABASE D ON M.DBCODE = D.DBCODE
        UNION ALL
        SELECT FM.PROTEIN_AC, FM.METHOD_AC, D.DBNAME,
               FM.POS_FROM, FM.POS_TO, NULL
        FROM INTERPRO.FEATURE_MATCH FM
        INNER JOIN INTERPRO.CV_DATABASE D ON FM.DBCODE = D.DBCODE
        WHERE FM.DBCODE = 'g'
        """
    )

    i = 0
    for row in cur:
        yield (
            row[0],
            row[1],
            databases[row[2]],
            row[5] if row[5] else f"{row[3]}-{row[4]}-S"
        )

        i += 1
        if not i % 100000000:
            logger.info(f"{i:>13,}")

    cur.close()
    con.close()
    logger.info(f"{i:>13,}")


def import_matches(ora_url: str, pg_url: str):
    logger.info("populating: match")
    pg_con = psycopg2.connect(**url2dict(pg_url))
    with pg_con.cursor() as pg_cur:
        pg_cur.execute("TRUNCATE TABLE match")

        drop_index(pg_con, "match_protein_idx")
        # drop_index(pg_con, "match_signature_idx")

        pg_cur.execute("SELECT name, id FROM database")
        databases = dict(pg_cur.fetchall())

        logger.info("starting")
        gen = _iter_matches(ora_url, databases)
        pg_cur.copy_from(file=CsvIO(gen, sep='|'), table="match", sep='|')

        pg_con.commit()

        logger.info("analyze")
        pg_cur.execute("ANALYZE match")

        logger.info("index")
        pg_cur.execute(
            """
            CREATE INDEX match_protein_idx
            ON match (protein_acc)
            """
        )
        # pg_cur.execute(
        #     """
        #     CREATE INDEX match_signature_idx
        #     ON match (signature_acc)
        #     """
        # )

    pg_con.close()
    logger.info("complete")


def export_complete_sequence_matches(url: str, filepath: str):
    logger.info("exporting matches")
    with open(filepath, "wb") as fh:
        i = 0
        for row in _iter_complete_sequence_matches(url):
            pickle.dump(row, fh)
            i += 1
            if not i % 100000000:
                logger.info(f"{i:>13,}")

    logger.info(f"{i:>13,}")


def _iter_complete_sequence_matches(url: str, filepath: Optional[str]=None):
    if filepath:
        with open(filepath, "rb") as fh:
            while True:
                try:
                    row = pickle.load(fh)
                except EOFError:
                    break
                else:
                    yield row
    else:
        con = cx_Oracle.connect(url)
        cur = con.cursor()
        cur.execute(
            """
            SELECT
                P.PROTEIN_AC, P.LEN, P.DBCODE, P.TAX_ID, E.LEFT_NUMBER,
                M.METHOD_AC, M.POS_FROM, M.POS_TO, M.FRAGMENTS
            FROM INTERPRO.PROTEIN P
            INNER JOIN INTERPRO.ETAXI E
              ON P.TAX_ID = E.TAX_ID
              AND P.FRAGMENT = 'N'
            INNER JOIN INTERPRO.MATCH M
              ON P.PROTEIN_AC = M.PROTEIN_AC
            ORDER BY P.PROTEIN_AC
            """
        )
        yield from cur
        cur.close()
        con.close()


def process_complete_sequence_matches(ora_url: str, pg_url: str, **kwargs):
    names_db = kwargs.get("names")
    if names_db:
        keep_db = True
    else:
        keep_db = False
        fd, names_db = mkstemp(dir=kwargs.get("dir"))
        os.close(fd)
        os.remove(names_db)
        export_names(pg_url, names_db)

    logger.info("loading reviewed proteins")
    pg_con = psycopg2.connect(**url2dict(pg_url))
    with pg_con.cursor() as pg_cur:
        pg_cur.execute(
            """
            SELECT accession
            FROM protein
            WHERE is_reviewed IS TRUE
            """
        )
        reviewed = {row[0] for row in pg_cur}
    pg_con.close()

    logger.info("starting")
    with KVdb(names_db) as names:
        it = _iter_complete_sequence_matches(ora_url, kwargs.get("matches"))
        for obj in it:
            name = names[obj[0]]
            print(obj)
            print(name)
            print(obj[0] in reviewed)
            break

    if not keep_db:
        os.remove(names_db)
