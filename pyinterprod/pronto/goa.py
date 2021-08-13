# -*- coding: utf-8 -*-

import cx_Oracle
import psycopg2
from psycopg2.extras import execute_values

from pyinterprod import logger
from pyinterprod.utils.pg import url2dict


def import_annotations(ora_url: str, pg_url: str):
    pg_con = psycopg2.connect(**url2dict(pg_url))
    with pg_con.cursor() as pg_cur:
        for name in ("protein2go", "publication", "term"):
            pg_cur.execute(f"DROP TABLE IF EXISTS {name}")

        pg_cur.execute(
            """
            CREATE TABLE protein2go (
                protein_acc VARCHAR(15) NOT NULL,
                term_id VARCHAR(10) NOT NULL,
                ref_db_code VARCHAR(10) NOT NULL,
                ref_db_id VARCHAR(60) NOT NULL
            )
            """
        )
        pg_cur.execute(
            """
            CREATE TABLE publication (
                id VARCHAR(25) NOT NULL 
                    CONSTRAINT publication_pkey PRIMARY KEY,
                title VARCHAR(1500) NOT NULL,
                published DATE NOT NULL
            )
            """
        )
        pg_cur.execute(
            """
            CREATE TABLE term (
                id VARCHAR(10) NOT NULL 
                    CONSTRAINT term_pkey PRIMARY KEY,
                name VARCHAR(200) NOT NULL,
                category VARCHAR(25) NOT NULL,
                num_constraints INTEGER NOT NULL,
                is_obsolete BOOLEAN NOT NULL,
                definition VARCHAR NOT NULL,
                replaced_id VARCHAR(10)
            )
            """
        )

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
            FROM GO.ANNOTATIONS A
            INNER JOIN GO.ECO2EVIDENCE E ON A.ECO_ID = E.ECO_ID
            INNER JOIN GO.CV_SOURCES S ON S.CODE = A.SOURCE
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
            FROM GO.PUBLICATIONS
            WHERE ID IN (
              SELECT DISTINCT A.REF_DB_ID
              FROM GO.ANNOTATIONS A
              INNER JOIN GO.ECO2EVIDENCE E ON A.ECO_ID = E.ECO_ID
              INNER JOIN GO.CV_SOURCES S ON S.CODE = A.SOURCE
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
            FROM GO.ANCESTORS
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
            FROM GO.TERM_TAXON_CONSTRAINTS
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
            FROM GO.TERMS T
            INNER JOIN GO.DEFINITIONS D ON T.GO_ID = D.GO_ID
            INNER JOIN GO.CV_CATEGORIES C ON T.CATEGORY = C.CODE
            UNION ALL
            SELECT S.SECONDARY_ID, T.NAME, C.TERM_NAME, T.IS_OBSOLETE,
                   D.DEFINITION, T.GO_ID
            FROM GO.SECONDARIES S
            INNER JOIN GO.TERMS T ON S.GO_ID = T.GO_ID
            INNER JOIN GO.CV_CATEGORIES C ON T.CATEGORY = C.CODE
            INNER JOIN GO.DEFINITIONS D ON T.GO_ID = D.GO_ID
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

        pg_cur.execute(
            """
            CREATE INDEX protein2go_protein_idx
            ON protein2go (protein_acc)
            """
        )
        pg_cur.execute(
            """
            CREATE INDEX protein2go_term_idx
            ON protein2go (term_id)
            """
        )
        pg_cur.execute("ANALYZE publication")
        pg_cur.execute("ANALYZE term")
        pg_con.commit()

    pg_con.close()
    logger.info("complete")


def _get_constraints(term_id: str, ancestors: dict, constraints: dict) -> set:
    result = constraints.get(term_id, set())

    for parent_id in ancestors.get(term_id, []):
        result |= _get_constraints(parent_id, ancestors, constraints)

    return result
