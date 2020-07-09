# -*- coding: utf-8 -*-

import cx_Oracle
import psycopg2
from psycopg2.extras import execute_values

from pyinterprod import logger
from pyinterprod.utils.pg import url2dict


def import_annotations(ora_url: str, pg_url: str):
    pg_con = psycopg2.connect(**url2dict(pg_url))
    with pg_con.cursor() as pg_cur:
        pg_cur.execute("TRUNCATE TABLE protein2go")
        pg_cur.execute("TRUNCATE TABLE publication")
        pg_cur.execute("TRUNCATE TABLE term")

        ora_con = cx_Oracle.connect(ora_url)
        ora_cur = ora_con.cursor()

        logger.info("populating: protein2go")
        """
        Filtering on length:
        Some annotations are not on proteins, but on post-translation modifications or processing events.
          e.g. P27958:PRO_0000037566 (protein: P27958; chain: PRO_0000037573)

        Protein accessions are 15 characters long (max), so anything longer than 15 characters cannot be an accession.
        A better (but heavier) approach would be to join with our PROTEIN table.
        """
        ora_cur.execute(
            """
            SELECT A.ENTITY_ID, A.GO_ID, A.REF_DB_CODE, A.REF_DB_ID
            FROM GO.ANNOTATIONS@GOAPRO A
            INNER JOIN GO.ECO2EVIDENCE@GOAPRO E ON A.ECO_ID = E.ECO_ID
            INNER JOIN GO.CV_SOURCES@GOAPRO S ON S.CODE = A.SOURCE
            WHERE A.ENTITY_TYPE = 'protein'
            AND LENGTH(A.ENTITY_ID) <= 15
            AND E.GO_EVIDENCE != 'IEA'
            AND S.IS_PUBLIC = 'Y'
            """
        )

        sql = "INSERT INTO protein2go VALUES %s"
        execute_values(pg_cur, sql, ora_cur, page_size=1000)

        logger.info("populating: publication")
        ora_cur.execute(
            """
            SELECT ID, TITLE, FIRST_PUBLISH_DATE
            FROM GO.PUBLICATIONS@GOAPRO
            WHERE ID IN (
              SELECT DISTINCT A.REF_DB_ID
              FROM GO.ANNOTATIONS@GOAPRO A
              INNER JOIN GO.ECO2EVIDENCE@GOAPRO E ON A.ECO_ID = E.ECO_ID
              INNER JOIN GO.CV_SOURCES@GOAPRO S ON S.CODE = A.SOURCE
              WHERE A.ENTITY_TYPE = 'protein'
              AND LENGTH(A.ENTITY_ID) <= 15
              AND E.GO_EVIDENCE != 'IEA'
              AND S.IS_PUBLIC = 'Y'
              AND A.REF_DB_CODE = 'PMID'
            )
            """
        )
        sql = "INSERT INTO publication VALUES %s"
        execute_values(pg_cur, sql, ora_cur, page_size=1000)

        logger.info("populating: term")
        ora_cur.execute(
            """
            SELECT CHILD_ID, PARENT_ID
            FROM GO.ANCESTORS@GOAPRO
            WHERE CHILD_ID != PARENT_ID
            """
        )
        ancestors = {}
        for term_id, parent_id in ora_cur:
            try:
                ancestors[term_id].add(parent_id)
            except KeyError:
                ancestors[term_id] = {parent_id}

        ora_cur.execute(
            """
            SELECT DISTINCT GO_ID, CONSTRAINT_ID
            FROM GO.TERM_TAXON_CONSTRAINTS@GOAPRO
            """
        )
        constraints = {}
        for term_id, constraint_id in ora_cur:
            try:
                constraints[term_id].add(constraint_id)
            except KeyError:
                constraints[term_id] = {constraint_id}

        ora_cur.execute(
            """
            SELECT T.GO_ID, T.NAME, C.TERM_NAME, T.IS_OBSOLETE,
                   D.DEFINITION, NULL
            FROM GO.TERMS@GOAPRO T
            INNER JOIN GO.DEFINITIONS@GOAPRO D ON T.GO_ID = D.GO_ID
            INNER JOIN GO.CV_CATEGORIES@GOAPRO C ON T.CATEGORY = C.CODE
            UNION ALL
            SELECT S.SECONDARY_ID, T.NAME, C.TERM_NAME, T.IS_OBSOLETE,
                   D.DEFINITION, T.GO_ID
            FROM GO.SECONDARIES@GOAPRO S
            INNER JOIN GO.TERMS@GOAPRO T ON S.GO_ID = T.GO_ID
            INNER JOIN GO.CV_CATEGORIES@GOAPRO C ON T.CATEGORY = C.CODE
            INNER JOIN GO.DEFINITIONS@GOAPRO D ON T.GO_ID = D.GO_ID
            """
        )

        sql = "INSERT INTO term VALUES %s"
        execute_values(pg_cur, sql, ((
            row[0],
            row[1],
            row[2],
            len(_get_constraints(row[0], ancestors, constraints)),
            row[3] == 'Y',
            row[4],
            row[5]
        ) for row in ora_cur), page_size=1000)

        ora_cur.close()
        ora_con.close()

        pg_con.commit()

        pg_cur.execute("ANALYZE protein2go")
        pg_cur.execute("ANALYZE publication")
        pg_cur.execute("ANALYZE term")

    pg_con.close()
    logger.info("complete")


def _get_constraints(term_id: str, ancestors: dict, constraints: dict) -> set:
    result = constraints.get(term_id, set())

    for parent_id in ancestors.get(term_id, []):
        result |= _get_constraints(parent_id, ancestors, constraints)

    return result
