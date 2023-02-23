import os
import pickle
import re
import warnings
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Optional

import cx_Oracle

from pyinterprod import logger
from pyinterprod.pronto.signature import get_swissprot_descriptions
from pyinterprod.utils import Table, oracle as ora
from . import contrib
from .database import Database
from .match import FEATURE_MATCH_PARTITIONS, MATCH_PARTITIONS, SITE_PARTITIONS
from .match import get_sig_protein_counts


FILE_DB_SIG = "signatures.update.pickle"
FILE_SIG_DESCR = "signatures.descr.pickle"


def export_swissprot_descriptions(pg_uri, data_dir: str):
    with open(os.path.join(data_dir, FILE_SIG_DESCR), "wb") as fh:
        pickle.dump(get_swissprot_descriptions(pg_uri), fh)


def add_staging(uri: str, update: list[tuple[Database, list[str]]]):
    con = cx_Oracle.connect(uri)
    cur = con.cursor()

    mapping_pmid_pubid = get_mapping_pmid_pubid(cur)

    ora.drop_table(cur, "METHOD_STG", purge=True)
    cur.execute(
        """
        CREATE TABLE INTERPRO.METHOD_STG (
            METHOD_AC VARCHAR2(25) NOT NULL,
            NAME VARCHAR2(100),
            DBCODE CHAR(1) NOT NULL,
            DESCRIPTION VARCHAR2(400),
            SIG_TYPE CHAR(1) NOT NULL,
            ABSTRACT VARCHAR2(4000),
            ABSTRACT_LONG CLOB
        )
        """
    )

    sql = """
        INSERT INTO INTERPRO.METHOD_STG
        VALUES (:1, :2, :3, :4, :5, :6, :7)
    """
    with Table(con, sql) as table:
        errors = 0
        for db, db_sources in update:
            if db.identifier == 'H':
                # Pfam
                signatures = contrib.pfam.get_signatures(*db_sources)
            elif db.identifier == 'J':
                # CDD
                signatures = contrib.cdd.parse_signatures(*db_sources)
            elif db.identifier == 'M':
                # PROSITE profiles
                signatures = contrib.prosite.parse_profiles(*db_sources)
            elif db.identifier == 'N':
                # NCBIFam
                signatures = contrib.ncbifam.get_signatures(*db_sources)
            elif db.identifier == 'P':
                # PROSITE patterns
                signatures = contrib.prosite.parse_patterns(*db_sources)
            elif db.identifier == 'Q':
                # HAMAP
                signatures = contrib.hamap.parse_signatures(*db_sources)
            elif db.identifier == 'R':
                # SMART
                signatures = contrib.smart.parse_signatures(*db_sources)
            elif db.identifier == 'V':
                # PANTHER
                signatures = contrib.panther.parse_signatures(*db_sources)
            elif db.identifier == 'X':
                # CATH-Gene3D
                signatures = contrib.cath.parse_superfamilies(*db_sources)
            else:
                logger.error(f"{db.name}: unsupported member database")
                errors += 1
                continue

            for m in signatures:
                m.abstract = update_abstract_pmid(cur, m.abstract, m.accession, mapping_pmid_pubid)
                if m.abstract is None:
                    abstract = abstract_long = None
                elif len(m.abstract) <= 4000:
                    abstract = m.abstract
                    abstract_long = None
                else:
                    abstract = None
                    abstract_long = m.abstract

                if m.name:
                    # Sanitize name (strip + remove multi-spaces)
                    name = re.sub(r"\s{2,}", " ", m.name.strip())
                else:
                    name = None

                if m.description:
                    descr = re.sub(r"\s{2,}", " ", m.description.strip())
                else:
                    descr = None

                table.insert((
                    m.accession,
                    name,
                    db.identifier,
                    descr,
                    m.sig_type,
                    abstract,
                    abstract_long
                ))

    if errors:
        cur.close()
        con.close()
        raise RuntimeError(f"{errors} errors occurred")

    con.commit()

    code2name = {db.identifier: db.name for db, _ in update}
    cur.execute(
        """
        SELECT DBCODE, COUNT(*)
        FROM INTERPRO.METHOD_STG
        GROUP BY DBCODE
        """
    )
    for dbcode, cnt in cur:
        logger.info(f"{code2name[dbcode]:<30} {cnt:>10,} signatures")

    cur.close()
    con.close()


def track_signature_changes(ora_uri: str, pg_uri: str,
                            databases: list[Database], data_dir: str):
    os.makedirs(data_dir, exist_ok=True)

    # First, get the SwissProt descriptions (before the update)
    all_sig2descs = get_swissprot_descriptions(pg_uri)

    con = cx_Oracle.connect(ora_uri)
    cur = con.cursor()
    results = {}
    for db in databases:
        cur.execute(
            """
            SELECT M.METHOD_AC, M.NAME, M.DESCRIPTION, M.SIG_TYPE, EM.ENTRY_AC
            FROM INTERPRO.METHOD M
            LEFT OUTER JOIN INTERPRO.ENTRY2METHOD EM
            ON M.METHOD_AC = EM.METHOD_AC
            WHERE DBCODE = :1
            """, (db.identifier,)
        )
        old_signatures = {row[0]: row[1:] for row in cur}

        cur.execute(
            """
            SELECT METHOD_AC, NAME, DESCRIPTION, SIG_TYPE
            FROM INTERPRO.METHOD_STG
            WHERE DBCODE = :1
            """, (db.identifier,)
        )
        new_signatures = {row[0]: row[1:] for row in cur}

        deleted = []
        name_changes = []
        descr_changes = []
        type_changes = []
        for acc in sorted(old_signatures):
            old_name, old_descr, old_type, entry_acc = old_signatures[acc]

            try:
                new_name, new_descr, new_type = new_signatures.pop(acc)
            except KeyError:
                deleted.append((acc, old_name, old_descr, entry_acc))
                continue

            if old_name != new_name:
                name_changes.append((acc, old_name, new_name))

            if old_descr != new_descr:
                descr_changes.append((acc, old_descr, new_descr))

            if old_type != new_type:
                type_changes.append((acc, old_type, new_type))

        results[db.identifier] = {
            "new": [
                (acc, *new_signatures[acc])
                for acc in sorted(new_signatures)
            ],
            "deleted": deleted,
            "changes": {
                "names": name_changes,
                "descriptions": descr_changes,
                "types": type_changes
            },
            "proteins": get_sig_protein_counts(cur, db.identifier),
            "descriptions": {
                acc: all_sig2descs[acc]
                for acc in all_sig2descs if acc in old_signatures
            }
        }

        logger.info(db.name)
        logger.info(f"  new:     {len(results[db.identifier]['new']):>6}")
        logger.info(f"  deleted: {len(results[db.identifier]['deleted']):>6}")

    cur.close()
    con.close()

    with open(os.path.join(data_dir, FILE_DB_SIG), "wb") as fh:
        pickle.dump(results, fh)


def delete_from_table(uri: str, table: str, partition: Optional[str],
                      column: str, step: int, stop: int) -> int:
    con = cx_Oracle.connect(uri)
    cur = con.cursor()

    if partition:
        db_obj = f"{table} PARTITION ({partition})"
    else:
        db_obj = table

    cur.execute(
        f"""
        SELECT COUNT(*)
        FROM INTERPRO.{db_obj}
        WHERE {column} IN (
            SELECT METHOD_AC
            FROM INTERPRO.METHOD_TO_DELETE
        )
        """
    )
    num_rows, = cur.fetchone()

    if not num_rows:
        cur.close()
        con.close()
        return num_rows

    for i in range(1, stop, step):
        cur.execute(
            f"""
            DELETE FROM INTERPRO.{db_obj}
            WHERE {column} IN (
              SELECT METHOD_AC
              FROM INTERPRO.METHOD_TO_DELETE
              WHERE ID BETWEEN :1 and :2
            )
            """, (i, i + step - 1)
        )

    con.commit()
    ora.gather_stats(cur, "INTERPRO", table, partition)
    cur.close()
    con.close()
    return num_rows


def delete_obsoletes(uri: str, databases: list[Database], **kwargs):
    step = kwargs.get("step", 10000)
    threads = kwargs.get("threads", 8)

    con = cx_Oracle.connect(uri)
    cur = con.cursor()

    # track signatures that need to be deleted
    ora.drop_table(cur, "METHOD_TO_DELETE")
    cur.execute(
        """
        CREATE TABLE INTERPRO.METHOD_TO_DELETE (
            ID NUMBER NOT NULL,
            METHOD_AC VARCHAR2(25) NOT NULL
        )
        """
    )

    for db in databases:
        cur.execute(
            """
            INSERT INTO INTERPRO.METHOD_TO_DELETE (ID, METHOD_AC)
            SELECT ROWNUM, METHOD_AC
            FROM (
                SELECT METHOD_AC
                FROM INTERPRO.METHOD
                WHERE DBCODE = :dbcode
                MINUS
                SELECT METHOD_AC
                FROM INTERPRO.METHOD_STG
                WHERE DBCODE = :dbcode
            )
            """, dbcode=db.identifier
        )

    con.commit()
    cur.execute(
        """
        CREATE UNIQUE INDEX UI_METHOD_TO_DELETE
        ON INTERPRO.METHOD_TO_DELETE (ID)
        """
    )

    cur.execute("SELECT COUNT(*) FROM INTERPRO.METHOD_TO_DELETE")
    stop, = cur.fetchone()

    logger.info(f"{stop:,} signatures to delete")

    if not stop:
        # Nothing to delete
        cur.close()
        con.close()
        return

    # Get tables with a FOREIGN KEY to INTERPRO.METHOD
    tables = []
    child_tables = ora.get_child_tables(cur, "INTERPRO", "METHOD")
    for table, constraint, column in child_tables:
        tables.append((table, constraint, column))

    # Add INTERPRO.METHOD as we want also to delete rows in this table
    tables.append(("METHOD", None, "METHOD_AC"))

    logger.info("disabling referential constraints")
    num_errors = 0
    for table, constraint, column in tables:
        if not constraint:
            continue

        try:
            ora.toggle_constraint(cur, table, constraint, False)
        except cx_Oracle.DatabaseError as exc:
            logger.error(exc)
            num_errors += 1

    if num_errors:
        cur.close()
        con.close()
        raise RuntimeError(f"{num_errors} constraints could not be disabled")

    tasks = []
    for table, constraint, column in tables:
        if table == "MATCH":
            for db in databases:
                partition = MATCH_PARTITIONS[db.identifier]
                logger.info(f"truncating {table} ({partition})")
                ora.truncate_partition(cur, table, partition)
        elif table == "SITE_MATCH":
            for db in databases:
                partition = SITE_PARTITIONS[db.identifier]
                logger.info(f"truncating {table} ({partition})")
                ora.truncate_partition(cur, table, partition)
        else:
            tasks.append((table, None, column))

    cur.close()
    con.close()

    with ThreadPoolExecutor(max_workers=threads) as executor:
        fs = {}

        for table, partition, column in tasks:
            args = (uri, table, partition, column, step, stop)
            f = executor.submit(delete_from_table, *args)
            fs[f] = (table, partition)

        num_errors = 0
        for f in as_completed(fs):
            table, partition = fs[f]
            if partition:
                name = f"{table} ({partition})"
            else:
                name = table

            try:
                num_rows = f.result()
            except Exception as exc:
                logger.info(f"{name}: failed ({exc})")
                num_errors += 1
            else:
                logger.info(f"{name}: {num_rows:,} rows deleted")

        if num_errors:
            raise RuntimeError(f"{num_errors} tables failed")

    logger.info("enabling referential constraints")
    con = cx_Oracle.connect(uri)
    cur = con.cursor()
    num_errors = 0
    constraints = set()
    for table, constraint, column in tables:
        if not constraint or constraint in constraints:
            """
            Either no constraint
            or prevent the same constrain to be enabled several times
            """
            continue

        logger.debug(f"enabling: {table}.{constraint}")
        constraints.add(constraint)
        try:
            ora.toggle_constraint(cur, table, constraint, True)
        except cx_Oracle.DatabaseError as exc:
            logger.error(exc)
            num_errors += 1

    if num_errors:
        cur.close()
        con.close()
        raise RuntimeError(f"{num_errors} constraints could not be enabled")

    for table, constraint, column in tables:
        for index in ora.get_indexes(cur, "INTERPRO", table):
            if index["is_unusable"]:
                logger.info(f"rebuilding index {index['name']}")
                ora.rebuild_index(cur, index["name"])

    cur.close()
    con.close()
    logger.info("complete")


def update_signatures(uri: str, go_sources: list[tuple[str, str]]):
    con = cx_Oracle.connect(uri)
    cur = con.cursor()
    cur.execute(
        """
        SELECT COUNT(*), DBCODE
        FROM INTERPRO.METHOD_STG
        GROUP BY DBCODE
        """
    )
    counts = cur.fetchall()

    cur.execute(
        """
        MERGE INTO INTERPRO.METHOD M
        USING INTERPRO.METHOD_STG S
          ON (M.METHOD_AC = S.METHOD_AC)
        WHEN MATCHED THEN
          UPDATE SET M.NAME = S.NAME,
                     M.DESCRIPTION = S.DESCRIPTION,
                     M.SIG_TYPE = S.SIG_TYPE,
                     M.ABSTRACT = S.ABSTRACT,
                     M.ABSTRACT_LONG = S.ABSTRACT_LONG,
                     M.TIMESTAMP = SYSDATE
        WHEN NOT MATCHED THEN
          INSERT (METHOD_AC, NAME, DBCODE, CANDIDATE, DESCRIPTION,
                  SIG_TYPE, ABSTRACT, ABSTRACT_LONG)
          VALUES (S.METHOD_AC, S.NAME, S.DBCODE, 'Y', S.DESCRIPTION,
                  S.SIG_TYPE, S.ABSTRACT, S.ABSTRACT_LONG)
        """
    )

    cur.executemany(
        """
        UPDATE INTERPRO.DB_VERSION
        SET ENTRY_COUNT = :1
        WHERE DBCODE = :2
        """, counts
    )

    con.commit()
    cur.close()
    con.close()

    for dbname, source in go_sources:
        if dbname == "panther":
            logger.info("updating PANTHER GO terms")
            contrib.panther.update_go_terms(uri, source)


def update_features(uri: str, update: list[tuple[Database, list[str]]]):
    con = cx_Oracle.connect(uri)
    cur = con.cursor()

    for db, db_sources in update:
        try:
            partition = FEATURE_MATCH_PARTITIONS[db.identifier]
        except KeyError:
            cur.close()
            con.close()
            raise ValueError(f"{db.name}: not a sequence feature database")

        # Remove matches
        ora.truncate_partition(cur, "INTERPRO.FEATURE_MATCH", partition)

        # Remove features
        cur.execute(
            """
            DELETE FROM INTERPRO.FEATURE_METHOD
            WHERE DBCODE = :1
            """,
            (db.identifier,)
        )

        if db.identifier == 'a':
            # AntiFam
            features = contrib.antifam.parse_models(*db_sources)
        elif db.identifier == 'd':
            # Pfam-N
            features = contrib.pfam.get_protenn_entries(cur, *db_sources)
        elif db.identifier == 'f':
            # FunFams
            features = contrib.cath.parse_functional_families(*db_sources)
        else:
            cur.close()
            con.close()
            raise ValueError(f"{db.name}: unsupported member database")

        # Add features
        params = []
        for f in features:
            params.append((
                f.accession,
                f.name if f.name else None,
                db.identifier,
                f.description,
                f.abstract[:4000] if f.abstract else None
            ))

        cur.executemany(
            """
            INSERT INTO INTERPRO.FEATURE_METHOD 
                (METHOD_AC, NAME, DBCODE, METHOD_DATE, TIMESTAMP, USERSTAMP, 
                DESCRIPTION, ABSTRACT) 
            VALUES (:1, :2, :3, SYSDATE, SYSDATE, USER, :4, :5)
            """,
            params
        )

        cur.execute(
            """
            UPDATE INTERPRO.DB_VERSION
            SET ENTRY_COUNT = :1
            WHERE DBCODE = :2
            """,
            (len(params), db.identifier)
        )

    con.commit()
    cur.close()
    con.close()


def get_mapping_pmid_pubid(cur: cx_Oracle.Cursor) -> dict:
    pmid_pubid = {}
    cur.execute(
        """
        SELECT DISTINCT PUB_ID, PUBMED_ID 
        FROM interpro.CITATION
        """
    )
    for pub_id, pubmed_id in cur:
        pmid_pubid.update({pubmed_id: pub_id})

    return pmid_pubid


def update_abstract_pmid(cur: cx_Oracle.Cursor, abstract: str, method_ac: str, pubid_pubmed: dict) -> str:
    method2pub = []
    pmids = re.findall('PMID:\s*[0-9]+', str(abstract))
    for pmid in pmids:
        pmid = int(pmid.split(':')[1].strip())
        try:
            pub_id = pubid_pubmed[pmid]
        except KeyError:
            pub_id = update_citation(cur, pmid)
        if pub_id:
            abstract = abstract.replace(f'PMID:{pmid}', f'[cite:{pub_id}]')
            method2pub.append(f"{pub_id}:{method_ac}")
    update_method2pub(cur, method2pub)
    return abstract


def update_method2pub(cur: cx_Oracle.Cursor, new_entries: list):
    current_method2pub = []
    cur.execute(
        """
        SELECT DISTINCT PUB_ID, METHOD_AC FROM interpro.METHOD2PUB
        """
    )
    for pub_id, method_ac in cur:
        current_method2pub.append(f"{pub_id}:{method_ac}")

    no_entries_method2pub = set(new_entries).difference(set(current_method2pub))
    for entry in no_entries_method2pub:
        pub_method = entry.split(':')
        cur.execute(
            """
            INSERT INTO interpro.METHOD2PUB (PUB_ID, METHOD_AC)
            VALUES (:1, :2)
            """,
            [pub_method[0], pub_method[1]]
        )


def update_citation(cur: cx_Oracle.Cursor, pmid: int):
    cur.execute(
        f"""
            SELECT
              C.EXTERNAL_ID, I.VOLUME, I.ISSUE, I.PUBYEAR, C.TITLE,
              C.PAGE_INFO, J.MEDLINE_ABBREVIATION, J.ISO_ABBREVIATION,
              A.AUTHORS, U.URL
            FROM CDB.CITATIONS@LITPUB C
              LEFT OUTER JOIN CDB.JOURNAL_ISSUES@LITPUB I
                ON C.JOURNAL_ISSUE_ID = I.ID
              LEFT JOIN CDB.CV_JOURNALS@LITPUB J
                ON I.JOURNAL_ID = J.ID
              LEFT OUTER JOIN CDB.FULLTEXT_URL@LITPUB U
                ON (
                    C.EXTERNAL_ID = U.EXTERNAL_ID AND
                    U.DOCUMENT_STYLE  ='DOI' AND
                    U.SOURCE = 'MED'
                )
              LEFT OUTER JOIN CDB.AUTHORS@LITPUB A
                ON (
                  C.ID = A.CITATION_ID AND
                  A.HAS_SPECIAL_CHARS = 'N'
                )
            WHERE C.EXTERNAL_ID = '{pmid}'
            """
    )

    citation = cur.fetchone()
    if citation:
        if len(citation[4]) > 740:
            citation[4] = citation[4][:737] + "..."
        pub_id = cur.var(cx_Oracle.STRING)
        cur.execute(
            """
            INSERT INTO INTERPRO.CITATION (
              PUB_ID, PUB_TYPE, PUBMED_ID, VOLUME, ISSUE,
              YEAR, TITLE, RAWPAGES, MEDLINE_JOURNAL,
              ISO_JOURNAL, AUTHORS, DOI_URL
            ) VALUES (
              INTERPRO.NEW_PUB_ID(), 'J', :1, :2, :3, :4, :5,
              :6, :7, :8, :9, :10
            )
            RETURNING PUB_ID INTO :11
            """,
            (*citation, pub_id)
        )
    else:
        warnings.warn(f"Citation related to PMID {pmid} not found.")
        return None
    return pub_id
