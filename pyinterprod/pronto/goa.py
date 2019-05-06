# -*- coding: utf-8 -*-

import cx_Oracle

from .. import orautils


def load_annotations(user: str, dsn: str):
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    orautils.drop_table(cur, owner, "PROTEIN2GO")
    cur.execute(
        """
        CREATE TABLE {}.PROTEIN2GO
        (
            PROTEIN_AC VARCHAR2(15) NOT NULL,
            GO_ID VARCHAR2(10) NOT NULL,
            EVIDENCE VARCHAR2(100) NOT NULL,
            REF_DB_CODE VARCHAR2(10),
            REF_DB_ID VARCHAR2(60)
        ) NOLOGGING
        """.format(owner)
    )

    """
    Filtering on length:
    Some annotations are not on proteins, but on post-translation modifications or processing events.
        e.g. P27958:PRO_0000037566 (protein: P27958; chain: PRO_0000037573)

    Protein accessions are 15 characters long (max), so anything above 15 characters cannot be an accession.
    A better (but heavier) approach would be to join with our PROTEIN table
    """
    cur.execute(
        """
        INSERT /*+ APPEND */ INTO {}.PROTEIN2GO (PROTEIN_AC, GO_ID, EVIDENCE, REF_DB_CODE, REF_DB_ID)
        SELECT A.ENTITY_ID, A.GO_ID, E.GO_EVIDENCE, A.REF_DB_CODE, A.REF_DB_ID
        FROM GO.ANNOTATIONS@GOAPRO A
        INNER JOIN GO.ECO2EVIDENCE@GOAPRO E ON A.ECO_ID = E.ECO_ID
        INNER JOIN GO.CV_SOURCES@GOAPRO S ON S.CODE = A.SOURCE
        WHERE A.ENTITY_TYPE = 'protein'
        AND LENGTH(A.ENTITY_ID) <= 15
        AND E.GO_EVIDENCE != 'IEA'
        AND S.IS_PUBLIC = 'Y'
        ORDER BY A.ENTITY_ID
        """.format(owner)
    )
    con.commit()

    cur.execute(
        """
        CREATE INDEX I_PROTEIN2GO$P$G
        ON {}.PROTEIN2GO (PROTEIN_AC, GO_ID) NOLOGGING
        """.format(owner)
    )
    cur.execute(
        """
        CREATE INDEX I_PROTEIN2GO$E
        ON {}.PROTEIN2GO (EVIDENCE) NOLOGGING
        """.format(owner)
    )
    cur.execute(
        """
        CREATE INDEX I_PROTEIN2GO$RC
        ON {}.PROTEIN2GO (REF_DB_CODE) NOLOGGING
        """.format(owner)
    )

    orautils.gather_stats(cur, owner, "PROTEIN2GO")
    orautils.grant(cur, owner, "PROTEIN2GO", "SELECT", "INTERPRO_SELECT")
    cur.close()
    con.close()


def load_publications(user: str, dsn: str):
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    orautils.drop_table(cur, owner, "PUBLICATION")
    cur.execute(
        """
        CREATE TABLE {}.PUBLICATION
        (
            ID VARCHAR2(25) NOT NULL,
            TITLE VARCHAR2(1500),
            FIRST_PUBLISHED_DATE DATE
        ) NOLOGGING
        """.format(owner)
    )

    cur.execute(
        """
        INSERT /*+APPEND*/ INTO {}.PUBLICATION (ID, TITLE, FIRST_PUBLISHED_DATE)
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
        """.format(owner)
    )
    con.commit()

    cur.execute(
        """
        ALTER TABLE {}.PUBLICATION
        ADD CONSTRAINT PK_PUBLICATION PRIMARY KEY (ID)
        """.format(owner)
    )

    orautils.gather_stats(cur, owner, "PUBLICATION")
    orautils.grant(cur, owner, "PUBLICATION", "SELECT", "INTERPRO_SELECT")
    cur.close()
    con.close()


def _traverse_ancestors(go_id: str, ancestors: dict, constraints: dict,
                        dst: set):
    dst |= constraints.get(go_id, set())

    for parent_id in ancestors.get(go_id, []):
        _traverse_ancestors(parent_id, ancestors, constraints, dst)


def load_terms(user: str, dsn: str):
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    cur.execute(
        """
        SELECT CHILD_ID, PARENT_ID
        FROM GO.ANCESTORS@GOAPRO
        WHERE CHILD_ID != PARENT_ID
        """
    )
    ancestors = {}
    for child_id, parent_id in cur:
        if child_id in ancestors:
            ancestors[child_id].add(parent_id)
        else:
            ancestors[child_id] = {parent_id}

    cur.execute(
        """
        SELECT DISTINCT GO_ID, CONSTRAINT_ID
        FROM GO.TERM_TAXON_CONSTRAINTS@GOAPRO
        """
    )
    constraints = {}
    for go_id, constraint_id in cur:
        if go_id in constraints:
            constraints[go_id].add(constraint_id)
        else:
            constraints[go_id] = {constraint_id}

    orautils.drop_table(cur, owner, "TERM")
    cur.execute(
        """
        CREATE TABLE {}.TERM
        (
            GO_ID VARCHAR2(10) NOT NULL,
            NAME VARCHAR2(200) NOT NULL,
            CATEGORY CHAR(1) NOT NULL ,
            NUM_CONSTRAINTS NUMBER NOT NULL,
            IS_OBSOLETE CHAR(1) NOT NULL ,
            DEFINITION VARCHAR2(4000),
            REPLACED_BY VARCHAR2(10)
        ) NOLOGGING
        """.format(owner)
    )

    cur.execute(
        """
        SELECT
          T.GO_ID, T.NAME, T.CATEGORY,
          T.IS_OBSOLETE, D.DEFINITION, NULL
        FROM GO.TERMS@GOAPRO T
        INNER JOIN GO.DEFINITIONS@GOAPRO D
          ON T.GO_ID = D.GO_ID
        UNION ALL
        SELECT
          S.SECONDARY_ID, T.NAME, T.CATEGORY,
          T.IS_OBSOLETE, D.DEFINITION, T.GO_ID
        FROM GO.SECONDARIES@GOAPRO S
        INNER JOIN GO.TERMS@GOAPRO T
          ON S.GO_ID = T.GO_ID
        INNER JOIN GO.DEFINITIONS@GOAPRO D
          ON T.GO_ID = D.GO_ID
        """
    )

    query = """
        INSERT /*+ APPEND */ INTO {}.TERM
        VALUES (:1, :2, :3, :4, :5, :6, :7)
    """.format(owner)
    table = orautils.TablePopulator(con, query, autocommit=True)
    for row in cur:
        term_constraints = set()
        _traverse_ancestors(row[0], ancestors, constraints, term_constraints)
        table.insert((row[0], row[1], row[2], len(term_constraints),
                     row[3], row[4], row[5]))

    table.close()

    cur.execute(
        """
        ALTER TABLE {}.TERM
        ADD CONSTRAINT PK_TERM PRIMARY KEY (GO_ID)
        """.format(owner)
    )

    orautils.gather_stats(cur, owner, "TERM")
    orautils.grant(cur, owner, "TERM", "SELECT", "INTERPRO_SELECT")
    cur.close()
    con.close()
