#!/usr/bin/env python

import os
import pickle
from concurrent.futures import as_completed, ThreadPoolExecutor
from multiprocessing import Process, Queue
from typing import Optional

import cx_Oracle

from .. import logger, orautils
from . import prediction, protein
from .utils import merge_organizers


def load_matches(user: str, dsn: str):
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    orautils.drop_table(cur, owner, "MATCH", purge=True)
    cur.execute(
        """
        CREATE TABLE {}.MATCH NOLOGGING
        AS
        SELECT
          PROTEIN_AC, METHOD_AC, DBCODE, POS_FROM, POS_TO, FRAGMENTS
        FROM INTERPRO.MATCH
        UNION ALL
        SELECT
          PROTEIN_AC, METHOD_AC, DBCODE, POS_FROM, POS_TO, NULL
        FROM INTERPRO.FEATURE_MATCH
        WHERE DBCODE = 'g'
        """.format(owner)
    )
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
    cur.close()
    con.close()


def load_signatures(user: str, dsn: str):
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    orautils.drop_table(cur, owner, "METHOD", purge=True)
    cur.execute(
        """
        CREATE TABLE {}.METHOD
        (
            METHOD_AC VARCHAR2(25) NOT NULL,
            NAME VARCHAR2(100),
            DBCODE CHAR(1) NOT NULL,
            CANDIDATE CHAR(1) NOT NULL,
            DESCRIPTION VARCHAR2(400),
            SIG_TYPE CHAR(1),
            ABSTRACT VARCHAR2(4000),
            ABSTRACT_LONG CLOB,
            PROTEIN_COUNT NUMBER DEFAULT 0 NOT NULL,
            FULL_SEQ_COUNT NUMBER DEFAULT 0 NOT NULL,
            RESIDUE_COUNT NUMBER DEFAULT 0 NOT NULL,
            DESC_COUNT NUMBER DEFAULT 0 NOT NULL,
            RANK_COUNT VARCHAR2(250) DEFAULT '{{}}' NOT NULL,
            TERM_COUNT NUMBER DEFAULT 0 NOT NULL
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


def update_signatures(user: str, dsn: str):
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    cur.execute(
        """
        MERGE INTO {0}.METHOD ME
        USING (
            SELECT METHOD_AC, COUNT(DISTINCT PROTEIN_AC) PROTEIN_COUNT
            FROM {0}.MATCH
            GROUP BY METHOD_AC
        ) MA
        ON (ME.METHOD_AC = MA.METHOD_AC)
        WHEN MATCHED THEN UPDATE SET ME.PROTEIN_COUNT = MA.PROTEIN_COUNT
        """.format(owner)
    )
    con.commit()
    cur.close()
    con.close()


def _export_matches(user: str, dsn: str, dst: str, buffer_size: int=1000000):
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
    pool = []
    for _ in range(max(1, processes-1)):
        p = Process(target=protein.consume_proteins,
                    args=(user, dsn, task_queue, done_queue, tmpdir))
        pool.append(p)
        p.start()

    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    orautils.drop_table(cur, owner, "METHOD2PROTEIN", purge=True)
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
    _protein_acc = dbcode = length = descid = tax_id = None
    num_proteins = 0
    for row in _get_matches(user, dsn, filepath):
        protein_acc = row[0]

        if protein_acc != _protein_acc:
            if _protein_acc:
                chunk.append((_protein_acc, dbcode, length, tax_id,
                              descid, matches))

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
            tax_id = row[3]
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
        chunk.append((_protein_acc, dbcode, length, tax_id,
                      descid, matches))
        task_queue.put(chunk)
        chunk = []
        matches = []
        num_proteins += 1

    for _ in pool:
        task_queue.put(None)

    names = []
    taxa = []
    terms = []
    comparators = []
    size = 0
    for _ in pool:
        # 1 comparator, 3 organizers, size
        obj = done_queue.get()
        comparators.append(obj[0])
        names.append(obj[1])
        taxa.append(obj[2])
        terms.append(obj[3])
        size += obj[4]

    for p in pool:
        p.join()

    logger.debug("proteins: {:,}".format(num_proteins))

    with ThreadPoolExecutor() as executor:
        fs = [executor.submit(_finalize_method2protein, user, dsn)]

        with open(os.path.join(tmpdir, "comparators.p"), "wb") as fh:
            pickle.dump(comparators, fh)
        # prediction.load_comparators(user, dsn, comparators)

        _create_method_term(user, dsn, terms)
        _create_method_desc(user, dsn, names)
        _create_method_taxa(user, dsn, taxa)
        terms = names = taxa = None

        # prediction.cmp_terms(user, dsn, processes, tmpdir)
        # prediction.cmp_descriptions(user, dsn, processes, tmpdir)
        # prediction.cmp_taxa(user, dsn, processes, tmpdir)

        for f in as_completed(fs):
            exc = f.exception()
            if exc is not None:
                raise exc

    logger.info("disk usage: {:.0f}MB".format(size/1024**2))


    if num_errors:
        raise RuntimeError("one or more tables could not be created")


def _create_method_desc(user: str, dsn: str, organizers: list):
    logger.debug("creating METHOD_DESC")
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    orautils.drop_table(cur, owner, "METHOD_DESC", purge=True)
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
    cur.execute(
        """
        SELECT DESC_ID
        FROM {0}.DESC_VALUE
        WHERE TEXT LIKE 'Predicted protein%'
          OR TEXT LIKE 'Uncharacterized protein%'
        """.format(owner)
    )
    excluded_descr = {row[0] for row in cur}
    cur.close()

    table = orautils.TablePopulator(con,
                                    query="INSERT /*+ APPEND */ "
                                          "INTO {}.METHOD_DESC "
                                          "VALUES (:1, :2, :3, :4)".format(owner),
                                    autocommit=True)

    for acc, descriptions in merge_organizers(organizers, remove=True):
        counts = {}
        for descid, dbcode in descriptions:
            if descid in counts:
                dbcodes = counts[descid]
            else:
                dbcodes = counts[descid] = {'S': 0, 'T': 0}

            dbcodes[dbcode] += 1

        descriptions = set()
        for descid, dbcodes in counts.items():
            table.insert((acc, descid, dbcodes['S'], dbcodes['T']))

            if descid not in excluded_descr:
                descriptions.add(descid)

    table.close()
    con.close()
    _finalize_method_desc(user, dsn)


def _finalize_method_desc(user: str, dsn: str):
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    orautils.grant(cur, owner, "METHOD_DESC", "SELECT", "INTERPRO_SELECT")
    orautils.gather_stats(cur, owner, "METHOD_DESC")
    cur.execute(
        """
        CREATE INDEX I_METHOD_DESC
        ON {}.METHOD_DESC (METHOD_AC) NOLOGGING
        """.format(owner)
    )
    cur.close()
    con.close()
    logger.debug("METHOD_DESC ready")


def _create_method_taxa(user: str, dsn: str, organizers: list):
    logger.debug("creating METHOD_TAXA")
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    orautils.drop_table(cur, owner, "METHOD_TAXA", purge=True)
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

    taxa = {}
    cur.execute(
        """
        SELECT TAX_ID, RANK, RANK_TAX_ID
        FROM {}.LINEAGE
        """.format(owner)
    )
    for tax_id, rank, rank_tax_id in cur:
        if tax_id in taxa:
            taxa[tax_id][rank] = rank_tax_id
        else:
            taxa[tax_id] = {rank: rank_tax_id}
    cur.close()

    table = orautils.TablePopulator(con,
                                    query="INSERT /*+ APPEND */ "
                                          "INTO {}.METHOD_TAXA "
                                          "VALUES (:1, :2, :3, :4)".format(owner),
                                    autocommit=True)

    for acc, tax_ids in merge_organizers(organizers, remove=True):
        counts = {}
        for tax_id in tax_ids:
            try:
                ranks = taxa[tax_id]
            except KeyError:
                continue  # should never happend (or ETAXI is incomplete)

            for rank, rank_tax_id in ranks.items():
                if rank in counts:
                    if rank_tax_id in counts[rank]:
                        counts[rank][rank_tax_id] += 1
                    else:
                        counts[rank][rank_tax_id] = 1
                else:
                    counts[rank] = {rank_tax_id: 1}

        for rank in counts:
            for tax_id, count in counts[rank].items():
                table.insert((acc, rank, tax_id, count))

    table.close()
    con.close()
    _finalize_method_taxa(user, dsn)


def _finalize_method_taxa(user: str, dsn: str):
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    orautils.grant(cur, owner, "METHOD_TAXA", "SELECT", "INTERPRO_SELECT")
    orautils.gather_stats(cur, owner, "METHOD_TAXA")
    cur.execute(
        """
        CREATE INDEX I_METHOD_TAXA
        ON {}.METHOD_TAXA (METHOD_AC, RANK) NOLOGGING
        """.format(owner)
    )
    cur.close()
    con.close()
    logger.debug("METHOD_TAXA ready")


def _create_method_term(user: str, dsn: str, organizers: list):
    logger.debug("creating METHOD_TERM")
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    orautils.drop_table(cur, owner, "METHOD_TERM", purge=True)
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
    for acc, terms in merge_organizers(organizers, remove=True):
        counts = {}
        for go_id in terms:
            if go_id in counts:
                counts[go_id] += 1
            else:
                counts[go_id] = 1

        for go_id, count in counts.items():
            table.insert((acc, go_id, count))

    table.close()
    con.close()
    _finalize_method_term(user, dsn)


def _finalize_method_term(user: str, dsn: str):
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    orautils.grant(cur, owner, "METHOD_TERM", "SELECT", "INTERPRO_SELECT")
    orautils.gather_stats(cur, owner, "METHOD_TERM")
    cur.execute(
        """
        CREATE INDEX I_METHOD_TERM
        ON {}.METHOD_TERM (METHOD_AC) NOLOGGING
        """.format(owner)
    )
    cur.close()
    con.close()
    logger.debug("METHOD_TERM ready")


def _finalize_method2protein(user: str, dsn: str):
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    orautils.grant(cur, owner, "METHOD2PROTEIN", "SELECT", "INTERPRO_SELECT")
    orautils.gather_stats(cur, owner, "METHOD2PROTEIN")
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
    cur.close()
    con.close()
    logger.debug("METHOD2PROTEIN ready")
