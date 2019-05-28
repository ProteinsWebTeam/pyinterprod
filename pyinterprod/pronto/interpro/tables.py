# -*- coding: utf-8 -*-

import os
import pickle
from multiprocessing import Process, Queue
from typing import Optional

import cx_Oracle

from . import proteins, RANKS
from .utils import merge_organizers
from ... import logger, orautils


def load_databases(user: str, dsn: str):
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    orautils.drop_table(cur, owner, "CV_DATABASE")
    cur.execute(
        """
        CREATE TABLE {}.CV_DATABASE
        (
            DBCODE VARCHAR2(10) NOT NULL,
            DBNAME VARCHAR2(50) NOT NULL,
            DBSHORT VARCHAR2(10) NOT NULL,
            VERSION VARCHAR2(20),
            FILE_DATE DATE,
            IS_READY CHAR(1) DEFAULT 'N',
            CONSTRAINT PK_DATABASE PRIMARY KEY (DBCODE)
        ) NOLOGGING
        """.format(owner)
    )

    cur.execute(
        """
        INSERT /*+ APPEND */ INTO {}.CV_DATABASE (
            DBCODE, DBNAME, DBSHORT, VERSION, FILE_DATE
        )
        SELECT DB.DBCODE, DB.DBNAME, DB.DBSHORT, V.VERSION, V.FILE_DATE
        FROM INTERPRO.CV_DATABASE DB
        LEFT OUTER JOIN INTERPRO.DB_VERSION V
          ON DB.DBCODE = V.DBCODE
        """.format(owner)
    )
    con.commit()

    orautils.gather_stats(cur, owner, "CV_DATABASE")
    orautils.grant(cur, owner, "CV_DATABASE", "SELECT", "INTERPRO_SELECT")
    cur.close()
    con.close()


def load_matches(user: str, dsn: str):
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    orautils.drop_table(cur, owner, "MATCH")
    cur.execute(
        """
        CREATE TABLE {}.MATCH
        (
            PROTEIN_AC VARCHAR2(15) NOT NULL,
            METHOD_AC VARCHAR2(25) NOT NULL,
            DBCODE CHAR(1) NOT NULL,
            MODEL_AC VARCHAR2(25),
            POS_FROM NUMBER(5) NOT NULL,
            POS_TO NUMBER(5) NOT NULL,
            FRAGMENTS VARCHAR2(200) DEFAULT NULL
        ) NOLOGGING
        """.format(owner)
    )

    """
    Insert directly signature matches and MobiDB-lite matches
    For some signature matches, the HMM model accession
        is the signature accession itself: in such cases, use NULL
    """
    cur.execute(
        """
        INSERT /*+ APPEND */ INTO {}.MATCH
        SELECT
          PROTEIN_AC, METHOD_AC, DBCODE,
          CASE
            WHEN MODEL_AC IS NOT NULL AND MODEL_AC != METHOD_AC
            THEN MODEL_AC
            ELSE NULL
          END AS MODEL_AC,
          POS_FROM, POS_TO, FRAGMENTS
        FROM INTERPRO.MATCH
        UNION ALL
        SELECT
          PROTEIN_AC, METHOD_AC, DBCODE,
          NULL,
          POS_FROM, POS_TO, NULL
        FROM INTERPRO.FEATURE_MATCH
        WHERE DBCODE = 'g'
        """.format(owner)
    )
    con.commit()

    orautils.gather_stats(cur, owner, "MATCH")
    orautils.grant(cur, owner, "MATCH", "SELECT", "INTERPRO_SELECT")

    cur.execute(
        """
        CREATE INDEX I_MATCH$PROTEIN
        ON {}.MATCH (PROTEIN_AC) NOLOGGING
        """.format(owner)
    )
    cur.execute(
        """
        CREATE INDEX I_MATCH$METHOD
        ON {}.MATCH (METHOD_AC) NOLOGGING
        """.format(owner)
    )
    cur.execute(
        """
        CREATE INDEX I_MATCH$DBCODE
        ON {}.MATCH (DBCODE) NOLOGGING
        """.format(owner)
    )

    cur.execute(
        """
        SELECT M.METHOD_AC, COUNT(DISTINCT M.PROTEIN_AC)
        FROM INTERPRO.PROTEIN P
        INNER JOIN {}.MATCH M
          ON P.PROTEIN_AC = M.PROTEIN_AC
        WHERE P.FRAGMENT = 'N'
        GROUP BY M.METHOD_AC
        """.format(owner)
    )
    signatures = cur.fetchall()

    for acc, num_proteins in signatures:
        cur.execute(
            """
            UPDATE {0}.METHOD
            SET PROTEIN_COUNT = :1
            WHERE METHOD_AC = :2
            """.format(owner), (num_proteins, acc)
        )
    con.commit()
    cur.close()
    con.close()


def load_signatures(user: str, dsn: str):
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    orautils.drop_table(cur, owner, "METHOD")
    cur.execute(
        """
        CREATE TABLE {}.METHOD
        (
            METHOD_AC VARCHAR2(25) NOT NULL,
            NAME VARCHAR2(100),
            DBCODE CHAR(1) NOT NULL,
            CANDIDATE CHAR(1) NOT NULL,
            DESCRIPTION VARCHAR2(220),
            SIG_TYPE CHAR(1),
            ABSTRACT VARCHAR2(4000),
            ABSTRACT_LONG CLOB,
            PROTEIN_COUNT NUMBER(8) DEFAULT 0 NOT NULL
        ) NOLOGGING
        """.format(owner)
    )

    cur.execute(
        """
        INSERT /*+ APPEND */ INTO {}.METHOD (
            METHOD_AC, NAME, DBCODE, CANDIDATE,
            DESCRIPTION, SIG_TYPE, ABSTRACT, ABSTRACT_LONG
        )
        SELECT
            METHOD_AC, NAME, DBCODE, CANDIDATE,
            DESCRIPTION, SIG_TYPE, ABSTRACT, ABSTRACT_LONG
        FROM INTERPRO.METHOD
        """.format(owner)
    )
    con.commit()

    cur.execute(
        """
        ALTER TABLE {}.METHOD
        ADD CONSTRAINT PK_METHOD PRIMARY KEY (METHOD_AC)
        """.format(owner)
    )
    cur.execute(
        """
        CREATE INDEX I_METHOD$DBCODE
        ON {}.METHOD (DBCODE) NOLOGGING
        """.format(owner)
    )

    orautils.gather_stats(cur, owner, "METHOD")
    orautils.grant(cur, owner, "METHOD", "SELECT", "INTERPRO_SELECT")
    cur.close()
    con.close()


def load_taxa(user: str, dsn: str):
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    orautils.drop_table(cur, owner, "ETAXI")
    cur.execute(
        """
        CREATE TABLE {}.ETAXI
        NOLOGGING
        AS
        SELECT
          TAX_ID, PARENT_ID, SCIENTIFIC_NAME, RANK, FULL_NAME
        FROM INTERPRO.ETAXI
        """.format(owner)
    )
    orautils.gather_stats(cur, owner, "ETAXI")
    orautils.grant(cur, owner, "ETAXI", "SELECT", "INTERPRO_SELECT")
    cur.execute(
        """
        ALTER TABLE {}.ETAXI
        ADD CONSTRAINT PK_ETAXI PRIMARY KEY (TAX_ID)
        """.format(owner)
    )

    orautils.drop_table(cur, owner, "LINEAGE")
    cur.execute(
        """
        CREATE TABLE {}.LINEAGE
        (
            TAX_ID NUMBER(10) NOT NULL,
            RANK VARCHAR2(50) NOT NULL,
            RANK_TAX_ID NUMBER(10)
        ) NOLOGGING
        """.format(owner)
    )

    """
    taxID 131567 (cellular organisms) contains three superkingdoms:
        * Bacteria (2)
        * Archaea (2157)
        * Eukaryota (2759)

    therefore it is not needed (we don't want a meta-superkingdom)
    """
    cur.execute(
        """
        SELECT TAX_ID, PARENT_ID, RANK
        FROM {}.ETAXI
        WHERE TAX_ID != 131567
        """.format(owner)
    )
    taxa = {}
    for tax_id, parent_id, rank in cur:
        if parent_id == 1:
            taxa[tax_id] = ("superkingdom", parent_id)
        else:
            taxa[tax_id] = (rank, parent_id)

    table = orautils.TablePopulator(con,
                                    query="INSERT /*+ APPEND */ "
                                          "INTO {}.LINEAGE "
                                          "VALUES (:1, :2, :3)".format(owner),
                                    autocommit=True)
    for tax_id in taxa:
        rank, parent_id = taxa[tax_id]
        if rank in RANKS:
            table.insert((tax_id, rank, tax_id))

        while parent_id in taxa:
            rank_tax_id = parent_id
            rank, parent_id = taxa[rank_tax_id]
            if rank in RANKS:
                table.insert((tax_id, rank, rank_tax_id))
    table.close()

    orautils.gather_stats(cur, owner, "LINEAGE")
    orautils.grant(cur, owner, "LINEAGE", "SELECT", "INTERPRO_SELECT")
    cur.execute(
        """
        CREATE INDEX I_LINEAGE
        ON {}.LINEAGE (TAX_ID, RANK)
        NOLOGGING
        """.format(owner)
    )

    cur.close()
    con.close()


def load_proteins(user: str, dsn: str):
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    orautils.drop_table(cur, owner, "PROTEIN")
    cur.execute(
        """
        CREATE TABLE {}.PROTEIN
        NOLOGGING
        AS
        SELECT PROTEIN_AC, NAME, DBCODE, LEN, FRAGMENT, TAX_ID
        FROM INTERPRO.PROTEIN
        """.format(owner)
    )
    orautils.gather_stats(cur, owner, "PROTEIN")
    orautils.grant(cur, owner, "PROTEIN", "SELECT", "INTERPRO_SELECT")
    cur.execute(
        """
        ALTER TABLE {}.PROTEIN
        ADD CONSTRAINT PK_PROTEIN PRIMARY KEY (PROTEIN_AC)
        """.format(owner)
    )
    cur.execute(
        """
        CREATE INDEX I_PROTEIN$DBCODE
        ON {}.PROTEIN (DBCODE) NOLOGGING
        """.format(owner)
    )
    cur.execute(
        """
        CREATE INDEX I_PROTEIN$NAME
        ON {}.PROTEIN (NAME) NOLOGGING
        """.format(owner)
    )

    cur.close()
    con.close()


def export_matches(user: str, dsn: str, dst: str, buffer_size: int=1000000):
    with open(dst, "wb") as fh:
        buffer = []
        for row in _get_matches(user, dsn):
            buffer.append(row)
            if len(buffer) == buffer_size:
                pickle.dump(buffer, fh)
                buffer = []

        pickle.dump(buffer, fh)


def _get_matches(user: str, dsn: str, filepath: Optional[str]=None):
    if filepath is not None:
        with open(filepath, "rb") as fh:
            while True:
                try:
                    buffer = pickle.load(fh)
                except EOFError:
                    break
                else:
                    for row in buffer:
                        yield row
    else:
        owner = user.split('/')[0]
        con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
        cur = con.cursor()
        cur.execute(
            """
            SELECT
                P.PROTEIN_AC, P.LEN, P.DBCODE, P.TAX_ID, D.DESC_ID,
                MA.METHOD_AC, MA.POS_FROM, MA.POS_TO, MA.FRAGMENTS
            FROM INTERPRO.PROTEIN P
            INNER JOIN {}.PROTEIN_DESC D
              ON P.PROTEIN_AC = D.PROTEIN_AC
            INNER JOIN INTERPRO.MATCH MA
              ON P.PROTEIN_AC = MA.PROTEIN_AC
            WHERE P.FRAGMENT = 'N'
            ORDER BY P.PROTEIN_AC
            """.format(owner)
        )

        for row in cur:
            yield row

        cur.close()
        con.close()


def load_signature2protein(user: str, dsn: str, processes: int=1,
                           tmpdir: Optional[str]=None, chunk_size: int=1000,
                           filepath: Optional[str]=None):
    if tmpdir is not None:
        os.makedirs(tmpdir, exist_ok=True)

    task_queue = Queue(maxsize=1)
    done_queue = Queue()
    consumers = []
    for _ in range(max(processes-1, 1)):
        p = Process(target=proteins.consume_proteins,
                    args=(user, dsn, task_queue, done_queue, tmpdir))
        consumers.append(p)
        p.start()

    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    orautils.drop_table(cur, owner, "METHOD2PROTEIN")
    cur.execute(
        """
        CREATE TABLE {}.METHOD2PROTEIN
        (
            METHOD_AC VARCHAR2(25) NOT NULL,
            PROTEIN_AC VARCHAR2(15) NOT NULL,
            DBCODE CHAR(1) NOT NULL,
            MD5 VARCHAR(32) NOT NULL,
            LEN NUMBER(5) NOT NULL,
            TAX_ID NUMBER(10) NOT NULL,
            DESC_ID NUMBER(10) NOT NULL
        )
        PARTITION BY LIST (DBCODE) (
          PARTITION M2P_SWISSP VALUES ('S'),
          PARTITION M2P_TREMBL VALUES ('T')
        ) NOLOGGING
        """.format(owner)
    )
    cur.close()
    con.close()

    chunk = []
    matches = []
    _protein_acc = dbcode = length = descid = taxid = None
    num_proteins = 0
    for row in _get_matches(user, dsn, filepath):
        protein_acc = row[0]

        if protein_acc != _protein_acc:
            if _protein_acc:
                chunk.append((_protein_acc, dbcode, length, taxid, descid,
                              matches))

                if len(chunk) == chunk_size:
                    task_queue.put(chunk)
                    chunk = []

                num_proteins += 1
                if not num_proteins % 1e7:
                    logger.debug("proteins: {:,}".format(num_proteins))

            elif _protein_acc is None:
                logger.debug("processing proteins")

            _protein_acc = protein_acc
            length = row[1]
            dbcode = row[2]
            taxid = row[3]
            descid = row[4]
            matches = []

        signature_acc = row[5]
        pos_start = row[6]
        pos_end = row[7]
        fragments = row[8]
        if fragments is not None:
            for frag in fragments.split(','):
                """
                Format: START-END-STATUS
                Types:
                    * S: Continuous single chain domain
                    * N: N terminus discontinuous
                    * C: C terminus discontinuous
                    * NC: N and C terminus discontinuous
                """
                pos_start, pos_end, _ = frag.split('-')
                matches.append((signature_acc, int(pos_start), int(pos_end)))
        else:
            matches.append((signature_acc, pos_start, pos_end))

    if _protein_acc:
        # Last protein
        chunk.append((_protein_acc, dbcode, length, taxid, descid,
                      matches))
        task_queue.put(chunk)
        chunk = []
        matches = []
        num_proteins += 1
        logger.debug("proteins: {:,}".format(num_proteins))

    for _ in consumers:
        task_queue.put(None)

    names = []
    taxa = []
    terms = []
    comparators = []
    sizes = [0, 0, 0, 0]
    for _ in consumers:
        _names, _taxa, _terms, comparator, _sizes = done_queue.get()
        names.append(_names)
        taxa.append(_taxa)
        terms.append(_terms)
        comparators.append(comparator)

        for i, size in enumerate(_sizes):
            sizes[i] += size

    for p in consumers:
        p.join()

    logger.debug("disk usage:")
    logger.debug("\tdescr.: {:.0f} MB".format(sizes[0]/1024**2))
    logger.debug("\ttaxa: {:.0f} MB".format(sizes[1]/1024**2))
    logger.debug("\tGO terms: {:.0f} MB".format(sizes[2]/1024**2))
    logger.debug("\toverlaps: {:.0f} MB".format(sizes[3]/1024**2))
    logger.debug("\ttotal: {:.0f} MB".format(sum(sizes)/1024**2))

    consumers = [
        Process(target=_load_description_counts, args=(user, dsn, names)),
        Process(target=_load_taxonomy_counts, args=(user, dsn, taxa)),
        Process(target=_load_term_counts, args=(user, dsn, terms)),
        Process(target=_load_comparisons, args=(user, dsn, comparators))
    ]

    for p in consumers:
        p.start()

    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    orautils.gather_stats(cur, owner, "METHOD2PROTEIN")
    orautils.grant(cur, owner, "METHOD2PROTEIN", "SELECT", "INTERPRO_SELECT")
    cur.execute(
        """
        CREATE UNIQUE INDEX UI_METHOD2PROTEIN
        ON {}.METHOD2PROTEIN (METHOD_AC, PROTEIN_AC)
        NOLOGGING
        """.format(owner)
    )
    cur.execute(
        """
        CREATE INDEX I_METHOD2PROTEIN$M
        ON {}.METHOD2PROTEIN (METHOD_AC)
        NOLOGGING
        """.format(owner)
    )
    cur.execute(
        """
        CREATE INDEX I_METHOD2PROTEIN$P
        ON {}.METHOD2PROTEIN (PROTEIN_AC)
        NOLOGGING
        """.format(owner)
    )
    cur.execute(
        """
        CREATE INDEX I_METHOD2PROTEIN$T
        ON {}.METHOD2PROTEIN (TAX_ID)
        NOLOGGING
        """.format(owner)
    )
    logger.debug("METHOD2PROTEIN ready")
    cur.close()
    con.close()

    for p in consumers:
        p.join()


def _load_comparisons(user: str, dsn: str, comparators: list):
    signatures = {}
    comparisons = {}
    for c in comparators:
        for acc_1, cnts, _comparisons in c:
            if acc_1 in signatures:
                for i, cnt in enumerate(cnts):
                    signatures[acc_1][i] += cnt

                for acc_2, cnts in _comparisons.items():
                    if acc_2 in comparisons[acc_1]:
                        for i, cnt in enumerate(cnts):
                            comparisons[acc_1][acc_2][i] += cnt
                    else:
                        comparisons[acc_1][acc_2] = cnts
            else:
                signatures[acc_1] = cnts
                comparisons[acc_1] = _comparisons

        c.remove()

    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    orautils.drop_table(cur, owner, "METHOD_COUNT")
    orautils.drop_table(cur, owner, "METHOD_COMPARISON")
    cur.execute(
        """
        CREATE TABLE {}.METHOD_COUNT
        (
            METHOD_AC VARCHAR2(25) NOT NULL,
            PROTEIN_COUNT NUMBER NOT NULL,
            RESIDUE_COUNT NUMBER NOT NULL
        ) NOLOGGING
        """.format(owner)
    )
    cur.execute(
        """
        CREATE TABLE {}.METHOD_COMPARISON
        (
            METHOD_AC1 VARCHAR2(25) NOT NULL,
            METHOD_AC2 VARCHAR2(25) NOT NULL,
            COLLOCATION NUMBER NOT NULL,
            PROTEIN_MATCH_OVERLAP_ NUMBER NOT NULL,
            PROTEIN_RESIDUE_OVERLAP NUMBER NOT NULL,
            RESIDUE_OVERLAP NUMBER NOT NULL
        ) NOLOGGING
        """.format(owner)
    )
    table = orautils.TablePopulator(con,
                                    query="INSERT /*+ APPEND */ "
                                          "INTO {}.METHOD_COUNT "
                                          "VALUES (:1, :2, :3)".format(owner),
                                    autocommit=True)
    for acc_1, cnts in signatures.items():
        table.insert((acc_1, *cnts))
    table.close()
    signatures = None

    table = orautils.TablePopulator(con,
                                    query="INSERT /*+ APPEND */ "
                                          "INTO {}.METHOD_COMPARISON "
                                          "VALUES (:1, :2, :3, :4, :5, :6)".format(owner),
                                    autocommit=True)
    for acc_1 in comparisons:
        for acc_2, cnts in comparisons[acc_1].items():
            table.insert((acc_1, acc_2, *cnts))
    table2.close()
    comparisons = None

    for table in ("METHOD_COUNT", "METHOD_COMPARISON"):
        orautils.gather_stats(cur, owner, table)
        orautils.grant(cur, owner, table, "SELECT", "INTERPRO_SELECT")

    cur.execute(
        """
        CREATE UNIQUE INDEX UI_METHOD_COUNT
        ON {}.METHOD_COUNT (METHOD_AC) NOLOGGING
        """.format(owner)
    )
    cur.execute(
        """
        CREATE INDEX I_METHOD_COMPARISON$AC1
        ON {}.METHOD_COMPARISON (METHOD_AC1) NOLOGGING
        """.format(owner)
    )
    cur.execute(
        """
        CREATE INDEX I_METHOD_COMPARISON$AC2
        ON {}.METHOD_COMPARISON (METHOD_AC2) NOLOGGING
        """.format(owner)
    )
    logger.debug("METHOD_COUNT and METHOD_COMPARISON ready")
    cur.close()
    con.close()


def _load_term_counts(user: str, dsn: str, organizers: list):
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    orautils.drop_table(cur, owner, "METHOD_TERM")
    cur.execute(
        """
        CREATE TABLE {}.METHOD_TERM
        (
            METHOD_AC VARCHAR2(25) NOT NULL,
            GO_ID VARCHAR2(10) NOT NULL,
            PROTEIN_COUNT NUMBER(10) NOT NULL
        ) NOLOGGING
        """.format(owner)
    )
    cur.close()

    table = orautils.TablePopulator(con,
                                    query="INSERT /*+ APPEND */ "
                                          "INTO {}.METHOD_TERM "
                                          "VALUES (:1, :2, :3)".format(owner),
                                    autocommit=True)

    for acc, terms in merge_organizers(organizers):
        counts = {}
        for go_id in terms:
            if go_id in _terms:
                counts[go_id] += 1
            else:
                counts[go_id] = 1

        for go_id, count in counts.items():
            table.insert((acc, go_id, count))

    table.close()

    for o in organizers:
        o.remove()

    cur = con.cursor()
    cur.execute(
        """
        CREATE INDEX I_METHOD_TERM
        ON {}.METHOD_TERM (METHOD_AC) NOLOGGING
        """.format(owner)
    )
    orautils.gather_stats(cur, owner, "METHOD_TERM")
    orautils.grant(cur, owner, "METHOD_TERM", "SELECT", "INTERPRO_SELECT")
    logger.debug("METHOD_TERM ready")
    cur.close()
    con.close()

def _load_description_counts(user: str, dsn: str, organizers: list):
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    orautils.drop_table(cur, owner, "METHOD_DESC")
    cur.execute(
        """
        CREATE TABLE {}.METHOD_DESC
        (
            METHOD_AC VARCHAR2(25) NOT NULL,
            DESC_ID NUMBER(10) NOT NULL,
            REVIEWED_COUNT NUMBER(10) NOT NULL,
            UNREVIEWED_COUNT NUMBER(10) NOT NULL
        ) NOLOGGING
        """.format(owner)
    )
    cur.close()

    table = orautils.TablePopulator(con,
                                    query="INSERT /*+ APPEND */ "
                                          "INTO {}.METHOD_DESC "
                                          "VALUES (:1, :2, :3, :4)".format(owner),
                                    autocommit=True)

    for acc, descriptions in merge_organizers(organizers):
        counts = {}
        for descid, dbcode in descriptions:
            if descid in counts:
                d = counts[descid]
            else:
                d = counts[descid] = {'S': 0, 'T': 0}

            counts[dbcode] += 1

        for descid, dbcodes in counts.items():
            table.insert((acc, descid, dbcodes['S'], dbcodes['T']))

    table.close()

    for o in organizers:
        o.remove()

    cur = con.cursor()
    cur.execute(
        """
        CREATE INDEX I_METHOD_DESC
        ON {}.METHOD_DESC (METHOD_AC) NOLOGGING
        """.format(owner)
    )
    orautils.gather_stats(cur, owner, "METHOD_DESC")
    orautils.grant(cur, owner, "METHOD_DESC", "SELECT", "INTERPRO_SELECT")
    logger.debug("METHOD_DESC ready")
    cur.close()
    con.close()


def _load_taxonomy_counts(user: str, dsn: str, organizers: list):
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    orautils.drop_table(cur, owner, "METHOD_TAXA")
    cur.execute(
        """
        CREATE TABLE {}.METHOD_TAXA
        (
            METHOD_AC VARCHAR2(25) NOT NULL,
            RANK VARCHAR2(50) NOT NULL,
            TAX_ID NUMBER(10) NOT NULL,
            PROTEIN_COUNT NUMBER(10) NOT NULL
        ) NOLOGGING
        """.format(owner)
    )

    orautils.drop_table(cur, owner, "METHOD_TAXA_TMP")
    cur.execute(
        """
        CREATE TABLE {}.METHOD_TAXA_TMP
        (
            METHOD_AC VARCHAR2(25) NOT NULL,
            TAX_ID NUMBER(10) NOT NULL
        ) NOLOGGING
        """.format(owner)
    )

    # Get lineages for the METHOD_TAXA table
    cur.execute(
        """
        SELECT TAX_ID, RANK, RANK_TAX_ID
        FROM {}.LINEAGE
        """.format(owner)
    )
    lineages = {}
    for tax_id, rank, rank_tax_id in cur:
        if tax_id in lineages:
            lineages[tax_id][rank] = rank_tax_id
        else:
            lineages[tax_id] = {rank: rank_tax_id}

    table1 = orautils.TablePopulator(con,
                                     query="INSERT /*+ APPEND */ "
                                           "INTO {}.METHOD_TAXA "
                                           "VALUES (:1, :2, :3, :4)".format(owner),
                                     autocommit=True)

    table2 = orautils.TablePopulator(con,
                                     query="INSERT /*+ APPEND */ "
                                           "INTO {}.METHOD_TAXA_TMP "
                                           "VALUES (:1, :2)".format(owner),
                                     autocommit=True)

    invalid_taxa = set()
    for acc, tax_ids in merge_organizers(organizers):
        counts = {}
        taxa = set()
        for tax_id in tax_ids:
            if tax_id in lineages:
                taxa.add(tax_id)  # we want unique Ids

                for rank, rank_tax_id in lineages[tax_id].items():
                    if rank in counts:
                        if rank_tax_id in counts[rank]:
                            counts[rank][rank_tax_id] += 1
                        else:
                            counts[rank][rank_tax_id] = 1
                    else:
                        counts[rank] = {rank_tax_id: 1}
            else:
                invalid_taxa.add(tax_id)

        for rank in counts:
            for rank_tax_id, count in counts[rank].items():
                table1.insert((acc, rank, rank_tax_id, count))

        for tax_id in taxa:
            table2.insert((acc, tax_id))

    table1.close()
    table2.close()

    if invalid_taxa:
        logger.warning("{} invalid tax IDs".format(len(invalid_taxa)))

    for o in organizers:
        o.remove()

    cur = con.cursor()
    cur.execute(
        """
        CREATE INDEX I_METHOD_TAXA
        ON {}.METHOD_TAXA (METHOD_AC, RANK) NOLOGGING
        """.format(owner)
    )
    orautils.gather_stats(cur, owner, "METHOD_TAXA")
    orautils.grant(cur, owner, "METHOD_TAXA", "SELECT", "INTERPRO_SELECT")
    logger.debug("METHOD_TAXA ready")
    cur.close()
    con.close()


def copy_schema(user_src: str, user_dst: str, dsn: str):
    _enable_schema(user_src, dsn)
    orautils.clear_schema(user_dst, dsn)


def _enable_schema(user: str, dsn: str):
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    cur.execute(
        """
        UPDATE {}.CV_DATABASE
        SET IS_READY = 'Y'
        """.format(owner)
    )
    con.commit()
    cur.close()
    con.close()
