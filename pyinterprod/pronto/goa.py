import oracledb
import psycopg

from pyinterprod import logger
from pyinterprod.utils.pg import url2dict


def import_annotations(ora_url: str, pg_url: str):
    pg_con = psycopg.connect(**url2dict(pg_url))
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

        ora_con = oracledb.connect(ora_url)
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

        records = []
        sql = """
            INSERT INTO protein2go 
                (protein_acc, term_id, ref_db_code, ref_db_id) 
            VALUES (%s, %s, %s, %s)
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

        records = []
        sql = """
            INSERT INTO publication (id, title, published) 
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

        records = []
        sql = """
            INSERT INTO term 
                (id, name, category, num_constraints, is_obsolete, definition, 
                replaced_id) 
            VALUES (%s, %s, %s, %s, %s, %s, %s)
        """

        for row in ora_cur:
            records.append(
                (
                    row[0],
                    row[1],
                    row[2],
                    len(_get_constraints(row[0], ancestors, constraints)),
                    row[3] == "Y",
                    row[4],
                    row[5],
                )
            )

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


def import_go_constraints(go_url: str, pg_url: str):
    ora_con = oracledb.connect(go_url)
    ora_cur = ora_con.cursor()

    ora_cur.execute(
        """
        SELECT UNION_ID, LISTAGG(TAXON_ID, ',')
        FROM GO.NCBI_TAXON_UNIONS
        GROUP BY UNION_ID
        """
    )
    union2taxa = {}
    for union_id, taxon_ids in ora_cur:
        union2taxa[union_id] = set(map(int, taxon_ids.split(",")))

    ora_cur.execute(
        """
        SELECT T.GO_ID, TC.RELATIONSHIP, TC.TAX_ID_TYPE, TC.TAX_ID
        FROM GO.TERMS T
        INNER JOIN GO.ANCESTORS A 
                ON T.GO_ID = A.CHILD_ID
        INNER JOIN GO.TERMS AT 
                ON A.PARENT_ID = AT.GO_ID
        INNER JOIN GO.TERM_TAXON_CONSTRAINTS TTC 
                ON AT.GO_ID = TTC.GO_ID
        INNER JOIN GO.TAXON_CONSTRAINTS TC 
                ON TTC.CONSTRAINT_ID = TC.CONSTRAINT_ID
        """
    )

    go2constraints = {}
    for go_id, relationship, tax_id_type, tax_id in ora_cur:
        if go_id not in go2constraints:
            go2constraints[go_id] = {}

        if relationship not in go2constraints[go_id]:
            go2constraints[go_id][relationship] = set()

        if tax_id_type == "NCBITaxon_Union":
            go2constraints[go_id][relationship] |= set(union2taxa[tax_id])
        elif tax_id_type == "NCBITaxon":
            go2constraints[go_id][relationship].add(tax_id)
        else:
            raise ValueError(f"Unknown TAX_ID_TYPE: {tax_id_type}")

    ora_cur.close()
    ora_con.close()

    pg_con = psycopg.connect(**url2dict(pg_url))
    with pg_con.cursor() as pg_cur:
        pg_cur.execute(f"DROP TABLE IF EXISTS go2constraints")

        pg_cur.execute(
            """
            CREATE TABLE go2constraints (
                go_id VARCHAR(15) NOT NULL,
                relationship VARCHAR(50) NOT NULL,
                taxon INTEGER NOT NULL
            )
            """
        )

        sql = """
            INSERT INTO GO2CONSTRAINTS (go_id, relationship, taxon)
            VALUES (%s, %s, %s)
        """

        records = []
        for go_id, relat2const in go2constraints.items():
            for relation, constraints in relat2const.items():
                for taxon_id in constraints:
                    if taxon_id == 131567:
                        """
                        meta-superkingdom (131567) includes three superkingdoms
                            * Bacteria (2)
                            * Archaea (2157)
                            * Eukaryota (2759)
                        """
                        records += [
                            (go_id, relation, 2),
                            (go_id, relation, 2157),
                            (go_id, relation, 2759)
                        ]
                    else:
                        records.append((go_id, relation, taxon_id))

                    if len(records) >= 1000:
                        pg_cur.executemany(sql, records)
                        pg_con.commit()
                        records.clear()
        if records:
            pg_cur.executemany(sql, records)
            pg_con.commit()
            records.clear()

        pg_cur.execute(
            """
            CREATE INDEX go2constraints_go_idx
            ON go2constraints (go_id)
            """
        )
        pg_con.commit()

    pg_con.close()
