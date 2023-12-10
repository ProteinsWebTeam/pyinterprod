import os
import pickle
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime

import oracledb

from pyinterprod import logger
from pyinterprod.pronto.signature import get_swissprot_descriptions
from pyinterprod.utils import Table
from pyinterprod.utils import oracle as ora

from . import contrib
from .contrib.common import Method
from .database import Database
from .match import (
    FEATURE_MATCH_PARTITIONS,
    MATCH_PARTITIONS,
    SITE_PARTITIONS,
    get_sig_protein_counts,
)

FILE_DB_SIG = "signatures.update.pickle"
FILE_SIG_DESCR = "signatures.descr.pickle"


def export_swissprot_descriptions(pg_uri, data_dir: str):
    with open(os.path.join(data_dir, FILE_SIG_DESCR), "wb") as fh:
        pickle.dump(get_swissprot_descriptions(pg_uri), fh)


def add_staging(uri: str, update: list[tuple[Database, dict[str, str]]]):
    con = oracledb.connect(uri)
    cur = con.cursor()

    pmid2pubid = get_pmid2pubid(cur)
    method2pub = get_method2pub(cur)

    ora.drop_table(cur, "METHOD2PUB_STG", purge=True)
    cur.execute(
        """
        CREATE TABLE INTERPRO.METHOD2PUB_STG
        AS SELECT * FROM INTERPRO.METHOD2PUB WHERE 1 = 0        
        """
    )

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
        for db, db_props in update:
            if db.identifier == 'H':
                # Pfam
                signatures = contrib.pfam.get_signatures(
                    db_props["signatures"]
                )
            elif db.identifier == 'J':
                # CDD
                signatures = contrib.cdd.parse_signatures(
                    db_props["signatures"]
                )
            elif db.identifier == 'M':
                # PROSITE profiles
                signatures = contrib.prosite.parse_profiles(
                    db_props["signatures"]
                )
            elif db.identifier == 'N':
                # NCBIFam
                signatures = contrib.ncbifam.get_signatures(
                    cur, db_props["signatures"], db_props["triage"]
                )
            elif db.identifier == 'P':
                # PROSITE patterns
                signatures = contrib.prosite.parse_patterns(
                    db_props["signatures"]
                )
            elif db.identifier == 'Q':
                # HAMAP
                signatures = contrib.hamap.parse_signatures(
                    db_props["signatures"]
                )
            elif db.identifier == 'R':
                # SMART
                signatures = contrib.smart.parse_signatures(
                    db_props["signatures"]
                )
            elif db.identifier == 'V':
                # PANTHER
                signatures = contrib.panther.parse_signatures(
                    db_props["signatures"]
                )
            elif db.identifier == 'X':
                # CATH-Gene3D
                signatures = contrib.cath.parse_superfamilies(
                    db_props["signatures"]
                )
            else:
                logger.error(f"{db.name}: unsupported member database")
                errors += 1
                continue

            for m in signatures:
                method2pub[m.accession] = update_references(cur, m, pmid2pubid)
                abstract = abstract_long = None
                if m.abstract is not None:
                    abstract = enclose_paragraph(m.abstract)
                    if len(abstract) > 4000:
                        abstract_long = abstract
                        abstract = None

                if m.name:
                    # Sanitize name (strip + remove multi-spaces)
                    name = re.sub(r"\s{2,}", " ", m.name.strip())
                    if name.lower() == m.accession.lower():
                        name = None
                else:
                    name = None

                if m.description:
                    descr = re.sub(r"\s{2,}", " ", m.description.strip())
                    if descr.lower() == m.accession.lower():
                        descr = None
                else:
                    descr = None

                table.insert(
                    (
                        m.accession,
                        name,
                        db.identifier,
                        descr,
                        m.sig_type,
                        abstract,
                        abstract_long,
                    )
                )

    populate_method2pub_stg(cur, method2pub)

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


def enclose_paragraph(abstract: str) -> str:
    if not abstract.lower().startswith("<p>"):
        abstract = f"<p>{abstract}"

    if not abstract.lower().endswith("</p>"):
        abstract += "</p>"

    return abstract


def track_signature_changes(
    ora_uri: str, pg_uri: str, databases: list[Database], data_dir: str
):
    os.makedirs(data_dir, exist_ok=True)

    # First, get the SwissProt descriptions (before the update)
    all_sig2descs = get_swissprot_descriptions(pg_uri)

    con = oracledb.connect(ora_uri)
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
            """,
            (db.identifier,),
        )
        old_signatures = {row[0]: row[1:] for row in cur}

        cur.execute(
            """
            SELECT METHOD_AC, NAME, DESCRIPTION, SIG_TYPE
            FROM INTERPRO.METHOD_STG
            WHERE DBCODE = :1
            """,
            (db.identifier,),
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
            "new": [(acc, *new_signatures[acc]) for acc in sorted(new_signatures)],
            "deleted": deleted,
            "changes": {
                "names": name_changes,
                "descriptions": descr_changes,
                "types": type_changes,
            },
            "proteins": get_sig_protein_counts(cur, db.identifier),
            "descriptions": {
                acc: all_sig2descs[acc]
                for acc in all_sig2descs
                if acc in old_signatures
            },
        }

        logger.info(db.name)
        logger.info(f"  new:     {len(results[db.identifier]['new']):>6}")
        logger.info(f"  deleted: {len(results[db.identifier]['deleted']):>6}")

    cur.close()
    con.close()

    with open(os.path.join(data_dir, FILE_DB_SIG), "wb") as fh:
        pickle.dump(results, fh)


def delete_from_table(uri: str, table: str, partition: str | None,
                      column: str, step: int, stop: int) -> int:
    con = oracledb.connect(uri)
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
    (num_rows,) = cur.fetchone()

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
            """,
            (i, i + step - 1),
        )

    con.commit()
    ora.gather_stats(cur, "INTERPRO", table, partition)
    cur.close()
    con.close()
    return num_rows


def delete_obsoletes(uri: str, databases: list[Database], **kwargs):
    step = kwargs.get("step", 1000)
    threads = kwargs.get("threads", 8)
    truncate_match_tables = kwargs.get("truncate_match_tables", True)

    con = oracledb.connect(uri)
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

    i = 0

    for db in databases:
        cur.execute(
            """
            SELECT METHOD_AC
            FROM INTERPRO.METHOD
            WHERE DBCODE = :dbcode
            MINUS
            SELECT METHOD_AC
            FROM INTERPRO.METHOD_STG
            WHERE DBCODE = :dbcode
            """,
            dbcode=db.identifier,
        )

        params = []
        for accession, in cur.fetchall():
            i += 1
            params.append((i, accession))

        for j in range(0, len(params), 1000):
            cur.executemany(
                """
                INSERT INTO INTERPRO.METHOD_TO_DELETE (ID, METHOD_AC)
                VALUES (:1, :2)
                """,
                params[j:j+1000]
            )

    con.commit()
    cur.execute(
        """
        CREATE UNIQUE INDEX UI_METHOD_TO_DELETE
        ON INTERPRO.METHOD_TO_DELETE (ID)
        """
    )

    cur.execute("SELECT COUNT(*) FROM INTERPRO.METHOD_TO_DELETE")
    (stop,) = cur.fetchone()

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

    """
    Any row deleted in METHOD2PUB should be deleted from METHOD2PUB_STG
    (which does not have a FK to METHOD)
    """
    tables.append(("METHOD2PUB_STG", None, "METHOD_AC"))

    logger.info("disabling referential constraints")
    num_errors = 0
    for table, constraint, column in tables:
        if not constraint:
            continue

        try:
            ora.toggle_constraint(cur, table, constraint, False)
        except oracledb.DatabaseError as exc:
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

                if truncate_match_tables:
                    logger.info(f"truncating {table} ({partition})")
                    ora.truncate_partition(cur, table, partition)
                else:
                    tasks.append((table, partition, column))
        elif table == "SITE_MATCH":
            for db in databases:
                partition = SITE_PARTITIONS[db.identifier]

                if truncate_match_tables:
                    logger.info(f"truncating {table} ({partition})")
                    ora.truncate_partition(cur, table, partition)
                else:
                    tasks.append((table, partition, column))
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
    con = oracledb.connect(uri)
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
        except oracledb.DatabaseError as exc:
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
    con = oracledb.connect(uri)
    cur = con.cursor()
    ora.truncate_table(cur, "INTERPRO.METHOD2PUB", reuse_storage=True)

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
        """,
        counts,
    )

    cur.execute(
        """
        INSERT INTO INTERPRO.METHOD2PUB 
        SELECT * 
        FROM INTERPRO.METHOD2PUB_STG
        """
    )

    con.commit()
    cur.close()
    con.close()

    for dbname, source in go_sources:
        if dbname == "ncbifam":
            logger.info("updating NCBIFAM GO terms")
            contrib.ncbifam.update_go_terms(uri, source)
        elif dbname == "panther":
            logger.info("updating PANTHER GO terms")
            contrib.panther.update_go_terms(uri, source)


def update_features(uri: str, update: list[tuple[Database, dict[str, str]]]):
    con = oracledb.connect(uri)
    cur = con.cursor()

    for db, db_props in update:
        try:
            partition = FEATURE_MATCH_PARTITIONS[db.identifier]
        except KeyError:
            cur.close()
            con.close()
            raise ValueError(f"{db.name}: not a sequence feature database")

        # Remove matches
        ora.truncate_partition(cur, "INTERPRO.FEATURE_MATCH", partition)

        if db.identifier == "a":
            # AntiFam
            features = contrib.antifam.parse_models(db_props["signatures"])
        elif db.identifier == 'd':
            # Pfam-N
            file = db_props["signatures"]
            features = contrib.pfam.get_protenn_entries(cur, file)
        elif db.identifier == 'f':
            # FunFams
            file = db_props["signatures"]
            features = contrib.cath.parse_functional_families(file)
        elif db.identifier == 'l':
            # ELM
            features = contrib.elm.load_classes(cur,
                                                db_props["classes"],
                                                db_props["instances"],
                                                db_props["sequences"])
        elif db.identifier == 'x':
            # COILS
            features = []
            cur.execute(
                """
                SELECT METHOD_AC, NAME, METHOD_DATE, DESCRIPTION, ABSTRACT
                FROM INTERPRO.FEATURE_METHOD
                WHERE DBCODE = :1
                """,
                [db.identifier]
            )

            for row in cur.fetchall():
                features.append(Method(accession=row[0],
                                       sig_type=None,
                                       name=row[1],
                                       description=row[3],
                                       abstract=row[4],
                                       date=row[2]))
        else:
            cur.close()
            con.close()
            raise ValueError(f"{db.name}: unsupported member database")

        # Remove existing features
        cur.execute(
            """
            DELETE FROM INTERPRO.FEATURE_METHOD
            WHERE DBCODE = :1
            """,
            [db.identifier]
        )

        # Add new features
        params = []
        for f in features:
            params.append(
                (
                    f.accession,
                    f.name or None,
                    db.identifier,
                    f.date or datetime.now(),
                    f.description,
                    f.abstract[:4000] if f.abstract else None,
                )
            )

        cur.executemany(
            """
            INSERT INTO INTERPRO.FEATURE_METHOD (
                METHOD_AC, NAME, DBCODE, METHOD_DATE, TIMESTAMP, USERSTAMP, 
                DESCRIPTION, ABSTRACT
            ) VALUES (:1, :2, :3, :4, SYSDATE, USER, :5, :6)
            """,
            params
        )

        cur.execute(
            """
            UPDATE INTERPRO.DB_VERSION
            SET ENTRY_COUNT = :1
            WHERE DBCODE = :2
            """,
            [len(params), db.identifier]
        )

    con.commit()
    cur.close()
    con.close()


def get_pmid2pubid(cur: oracledb.Cursor) -> dict[int, str]:
    cur.execute(
        """
        SELECT DISTINCT PUBMED_ID, PUB_ID 
        FROM INTERPRO.CITATION
        """
    )
    return dict(cur.fetchall())


def get_method2pub(cur: oracledb.Cursor) -> dict[str, set]:
    cur.execute("SELECT METHOD_AC, PUB_ID FROM INTERPRO.METHOD2PUB")

    method2pub = {}
    for signature_acc, pub_id in cur:
        try:
            method2pub[signature_acc].add(pub_id)
        except KeyError:
            method2pub[signature_acc] = {pub_id}

    return method2pub


def update_references(cur: oracledb.Cursor, method: Method,
                      pmid2pubid: dict[int, str]) -> set[str]:
    pub_ids = set()
    if method.abstract is not None:
        text = method.abstract

        # Case found in PR01452
        text = re.sub(r"\s*\[PMID:\s*NOT FOUND]", r"", text, flags=re.I)

        pmids = re.finditer(r"PMID:\s*\.?([0-9]+)", text, flags=re.I)
        for match in pmids:
            pmid = int(match.group(1))
            try:
                pub_id = pmid2pubid[pmid]
            except KeyError:
                pub_id = update_citation(cur, pmid)
                pmid2pubid[pmid] = pub_id

            if pub_id:
                text = text.replace(match.group(0), f"[cite:{pub_id}]")
                pub_ids.add(pub_id)
            else:
                logger.warning(f"{method.accession}: "
                               f"no citation found with PubMed ID {pmid}")
                # text = text.replace(match.group(0), "")

        method.abstract = text

    # Required for PMIDs not in the abstract
    for pmid in method.references:
        try:
            pub_id = pmid2pubid[pmid]
        except KeyError:
            pub_id = update_citation(cur, pmid)
            pmid2pubid[pmid] = pub_id

        if pub_id:
            pub_ids.add(pub_id)

    return pub_ids


def populate_method2pub_stg(cur: oracledb.Cursor, method2pub: dict[str, set]):
    data = [(pmid, acc) for acc, pmids in method2pub.items() for pmid in pmids]
    step = 1000
    for i in range(0, len(data), step):
        cur.executemany(
            """
            INSERT INTO INTERPRO.METHOD2PUB_STG (PUB_ID, METHOD_AC)
            VALUES (:1, :2)
            """,
            data[i:i+step],
        )


def update_citation(cur: oracledb.Cursor, pmid: int) -> str | None:
    cur.execute(
        """
        SELECT EXTERNAL_ID, VOLUME, ISSUE, YEAR, TITLE, RAWPAGES, 
               MEDLINE_JOURNAL, ISO_JOURNAL, AUTHORS, DOI_URL
        FROM (
            SELECT
                C.EXTERNAL_ID AS EXTERNAL_ID, I.VOLUME AS VOLUME, 
                I.ISSUE AS ISSUE, I.PUBYEAR AS YEAR, C.TITLE AS TITLE, 
                C.PAGE_INFO AS RAWPAGES, 
                J.MEDLINE_ABBREVIATION AS MEDLINE_JOURNAL, 
                J.ISO_ABBREVIATION AS ISO_JOURNAL,
                A.AUTHORS AS AUTHORS, U.URL AS DOI_URL,
                ROW_NUMBER() OVER (
                      PARTITION BY C.EXTERNAL_ID
                      ORDER BY U.DATE_UPDATED DESC
                  ) R
            FROM CDB.CITATIONS@LITPUB C
            LEFT OUTER JOIN CDB.JOURNAL_ISSUES@LITPUB I
                ON C.JOURNAL_ISSUE_ID = I.ID
            LEFT JOIN CDB.CV_JOURNALS@LITPUB J
                ON I.JOURNAL_ID = J.ID
            LEFT OUTER JOIN CDB.FULLTEXT_URL_MEDLINE@LITPUB U
                ON (
                    C.EXTERNAL_ID = U.EXTERNAL_ID AND
                    UPPER(U.SITE) = 'DOI'
                )
            LEFT OUTER JOIN CDB.AUTHORS@LITPUB A
                ON (
                    C.ID = A.CITATION_ID AND
                    A.HAS_SPECIAL_CHARS = 'N'
                ) 
            WHERE C.EXTERNAL_ID = :1
        ) 
        WHERE R = 1
        """,
        (str(pmid),),
    )
    citation = cur.fetchone()

    if citation:
        citation = list(citation)
        citation[0] = int(citation[0])
        if len(citation[4]) > 740:
            citation[4] = citation[4][:737] + "..."

        pub_id = cur.var(oracledb.STRING)
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
            (*citation, pub_id),
        )
        return pub_id.getvalue()[0]
    else:
        return None


def update_citations(cur: oracledb.Cursor):
    cur.execute(
        """
        DELETE FROM INTERPRO.CITATION
        WHERE PUB_ID NOT IN 
        (
            SELECT PUB_ID FROM INTERPRO.ENTRY2PUB
            UNION
            SELECT PUB_ID FROM INTERPRO.METHOD2PUB
            UNION
            SELECT PUB_ID FROM INTERPRO.SUPPLEMENTARY_REF
        )
        """
    )

    pmids = _get_used_citations_pmids(cur)
    step = 1000
    for i in range(0, len(pmids), step):
        params = pmids[i:i+step]
        args = ",".join([":" + str(i + 1) for i in range(len(params))])

        cur.execute(
            f"""
            SELECT VOLUME, ISSUE, YEAR, TITLE, RAWPAGES, MEDLINE_JOURNAL, 
                   ISO_JOURNAL, AUTHORS, DOI_URL, EXTERNAL_ID
            FROM (
                SELECT
                    C.EXTERNAL_ID AS EXTERNAL_ID, I.VOLUME AS VOLUME, 
                    I.ISSUE AS ISSUE, I.PUBYEAR AS YEAR, C.TITLE AS TITLE, 
                    C.PAGE_INFO AS RAWPAGES, 
                    J.MEDLINE_ABBREVIATION AS MEDLINE_JOURNAL, 
                    J.ISO_ABBREVIATION AS ISO_JOURNAL,
                    A.AUTHORS AS AUTHORS, U.URL AS DOI_URL,
                    ROW_NUMBER() OVER (
                          PARTITION BY C.EXTERNAL_ID
                          ORDER BY U.DATE_UPDATED DESC
                      ) R
                FROM CDB.CITATIONS@LITPUB C
                LEFT OUTER JOIN CDB.JOURNAL_ISSUES@LITPUB I
                    ON C.JOURNAL_ISSUE_ID = I.ID
                LEFT JOIN CDB.CV_JOURNALS@LITPUB J
                    ON I.JOURNAL_ID = J.ID
                LEFT OUTER JOIN CDB.FULLTEXT_URL_MEDLINE@LITPUB U
                    ON (
                        C.EXTERNAL_ID = U.EXTERNAL_ID AND
                        UPPER(U.SITE) = 'DOI'
                    )
                LEFT OUTER JOIN CDB.AUTHORS@LITPUB A
                    ON (
                        C.ID = A.CITATION_ID AND
                        A.HAS_SPECIAL_CHARS = 'N'
                    )
                WHERE C.EXTERNAL_ID IN ({args})
            ) WHERE R = 1
            """,
            params,
        )
        citations = cur.fetchall()

        for citation in citations:
            citation = list(citation)
            if len(citation[3]) > 740:
                citation[3] = citation[3][:737] + "..."
            if citation[8] is not None:
                citation[8] = citation[8].lower().replace(" ", "")
            citation[9] = int(citation[9])

        cur.executemany(
            """
            UPDATE INTERPRO.CITATION
                SET VOLUME= :1, 
                    ISSUE= :2, 
                    YEAR= :3, 
                    TITLE= :4, 
                    RAWPAGES= :5, 
                    MEDLINE_JOURNAL= :6, 
                    ISO_JOURNAL= :7, 
                    AUTHORS= :8, 
                    DOI_URL= :9
            WHERE PUBMED_ID = :10
            """, 
            citations
        )


def _get_used_citations_pmids(cur: oracledb.Cursor) -> list[str]:
    cur.execute(
        """
        SELECT PUBMED_ID
        FROM INTERPRO.CITATION
        WHERE PUB_ID IN (
            SELECT PUB_ID FROM INTERPRO.ENTRY2PUB
            UNION
            SELECT PUB_ID FROM INTERPRO.SUPPLEMENTARY_REF
            UNION
            SELECT PUB_ID FROM INTERPRO.METHOD2PUB
        )
        """
    )
    return [str(pmid) for pmid, in cur.fetchall()]
