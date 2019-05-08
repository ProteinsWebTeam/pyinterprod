# -*- coding: utf-8 -*-

import heapq
import os
from multiprocessing import Process, Queue
from typing import Optional

import cx_Oracle

from ... import orautils
from . import proteins


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
        SELECT METHOD_AC, COUNT(DISTINCT PROTEIN_AC)
        FROM {}.MATCH
        GROUP BY METHOD_AC
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
          TAX_ID, PARENT_ID, SCIENTIFIC_NAME, RANK,
          LEFT_NUMBER, RIGHT_NUMBER, FULL_NAME
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
            LEFT_NUMBER NUMBER NOT NULL,
            TAX_ID NUMBER(10) NOT NULL,
            RANK VARCHAR2(50)
        ) NOLOGGING
        """.format(owner)
    )

    taxons = {}
    cur.execute(
        """
        SELECT TAX_ID, PARENT_ID, LEFT_NUMBER, RANK
        FROM {}.ETAXI
        """.format(owner)
    )
    for tax_id, parent_id, left_num, rank in cur:
        if tax_id != 131567:
            """
            taxID 131567 (cellular organisms) contains three superkingdoms:
                * Bacteria (2)
                * Archaea (2157)
                * Eukaryota (2759)

            therefore it is not needed (we don't want a meta-superkingdom)
            """
            taxons[tax_id] = {
                "parent": parent_id,
                "left_num": left_num,
                "rank": "superkingdom" if parent_id == 1 else rank
            }

    lineage = []
    for tax_id in taxons:
        t = taxons[tax_id]
        left_num = t["left_num"]
        if not left_num:
            continue
        elif t["rank"] != "no rank":
            lineage.append((left_num, tax_id, t["rank"]))

        parent_id = t["parent"]
        while parent_id:
            if parent_id not in taxons:
                # taxID 131567 missing from dictionary
                break

            t = taxons[parent_id]
            if t["rank"] != "no rank":
                lineage.append((left_num, parent_id, t["rank"]))

            parent_id = t["parent"]

    table = orautils.TablePopulator(con,
                                    query="INSERT /*+ APPEND */ "
                                          "INTO {}.LINEAGE "
                                          "VALUES (:1, :2, :3)".format(owner),
                                    autocommit=True)
    for record in lineage:
        table.insert(record)
    table.close()

    orautils.gather_stats(cur, owner, "LINEAGE")
    orautils.grant(cur, owner, "LINEAGE", "SELECT", "INTERPRO_SELECT")
    cur.execute(
        """
        CREATE INDEX I_LINEAGE$L$R
        ON {}.LINEAGE (LEFT_NUMBER, RANK)
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


def load_signature2protein(user: str, dsn: str, processes: int=1,
                           tmpdir: Optional[str]=None, chunk_size: int=1000):
    if tmpdir is not None:
        os.makedirs(tmpdir, exist_ok=True)

    task_queue = Queue(maxsize=processes)
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
            LEFT_NUMBER NUMBER NOT NULL,
            DESC_ID NUMBER(10) NOT NULL
        )
        PARTITION BY LIST (DBCODE) (
          PARTITION M2P_SWISSP VALUES ('S'),
          PARTITION M2P_TREMBL VALUES ('T')
        ) NOLOGGING
        """.format(owner)
    )

    cur.execute(
        """
        SELECT
            P.PROTEIN_AC, P.LEN, P.DBCODE, D.DESC_ID,
            NVL(E.LEFT_NUMBER, 0),
            MA.METHOD_AC, MA.POS_FROM, MA.POS_TO, MA.FRAGMENTS
        FROM INTERPRO.PROTEIN P
        INNER JOIN INTERPRO.ETAXI E
          ON P.TAX_ID = E.TAX_ID
        INNER JOIN {}.PROTEIN_DESC D
          ON P.PROTEIN_AC = D.PROTEIN_AC
        INNER JOIN INTERPRO.MATCH MA
          ON P.PROTEIN_AC = MA.PROTEIN_AC
        WHERE P.FRAGMENT = 'N'
        ORDER BY P.PROTEIN_AC
        """.format(owner)
    )
    chunk = []
    matches = []
    _protein_acc = dbcode= length = descid = left_num = None
    num_proteins = 0
    for row in cur:
        protein_acc = row[0]

        if protein_acc != _protein_acc:
            if _protein_acc:
                # acc, dbcode, length, descid, leftnum, matches
                chunk.append((_protein_acc, dbcode, length, descid, left_num,
                              matches))

                if len(chunk) == chunk_size:
                    task_queue.put(chunk)
                    chunk = []

                num_proteins += 1
                if not num_proteins % 1e7:
                    logger.debug("proteins: {:,}".format(num_proteins))

            else:
                logger.debug("processing proteins")

            length = row[1]
            dbcode = row[2]
            descid = row[3]
            left_num = row[4]

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
        chunk.append((_protein_acc, dbcode, length, descid, left_num,
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
    size = 0
    signatures = {}
    collocations = {}
    match_overlaps = {}
    residue_overlaps = {}
    for _ in consumers:
        _names, _taxa, comparator, _size = done_queue.get()
        names.append(_names)
        taxa.append(_taxa)
        size += _size

        for acc, s in comparator.signatures.items():
            if acc in signatures:
                for k, v in s.items():
                    signatures[acc][k] += v
            else:
                signatures[acc] = s

        for acc_1, s in comparator.collocations.items():
            if acc_1 in collocations:
                for acc_2, v in s.items():
                    if acc_2 in collocations[acc_1]:
                        collocations[acc_1][acc_2] += v
                    else:
                        collocations[acc_1][acc_2] = v
            else:
                collocations[acc1] = s

        for acc_1, s in comparator.match_overlaps.items():
            if acc_1 in match_overlaps:
                for acc_2, v in s.items():
                    if acc_2 in match_overlaps[acc_1]:
                        match_overlaps[acc_1][acc_2] += v
                    else:
                        match_overlaps[acc_1][acc_2] = v
            else:
                match_overlaps[acc1] = s

        for acc_1, s in comparator.residue_overlaps.items():
            if acc_1 in residue_overlaps:
                for acc_2, v in s.items():
                    if acc_2 in residue_overlaps[acc_1]:
                        residue_overlaps[acc_1][acc_2] += v
                    else:
                        residue_overlaps[acc_1][acc_2] = v
            else:
                residue_overlaps[acc1] = s

    for p in consumers:
        p.join()

    consumers = [
        Process(target=_load_description_counts, args=(user, dsn, names)),
        Process(target=_load_taxonomy_counts, args=(user, dsn, taxa)),
        Process(target=_load_comparisons, args=(user, dsn, signatures,
                                                collocations, match_overlaps,
                                                residue_overlaps))
    ]

    for p in consumers:
        p.start()

    orautils.gather_stats(cur, owner, "METHOD2PROTEIN")
    cur.close()
    con.close()
    logger.debug("METHOD2PROTEIN ready")

    for p in consumers:
        p.join()


def _load_comparisons(user: str, dsn: str, signatures: dict, collocs: dict,
                      match_ovrlps: dict, residue_ovrlps: dict):
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    orautils.drop_table(cur, owner, "METHOD_COUNT")
    cur.execute(
        """
        CREATE TABLE {}.METHOD_COUNT
        (
            METHOD_AC VARCHAR2(25) NOT NULL,
            PROTEIN_COUNT NUMBER NOT NULL,
            MATCH_COUNT NUMBER NOT NULL,
            RESIDUE_COUNT NUMBER NOT NULL
        ) NOLOGGING
        """.format(owner)
    )
    table = orautils.TablePopulator(con,
                                    query="INSERT /*+ APPEND */ "
                                          "INTO {}.METHOD_COUNT "
                                          "VALUES (:1, :2, :3, :4)".format(owner),
                                    autocommit=True)

    for acc, s in signatures.items():
        table.insert((acc, s["proteins"], s["matches"], s["residues"]))

    table.close()
    orautils.gather_stats(cur, owner, "METHOD_COUNT")
    cur.execute(
        """
        CREATE UNIQUE INDEX UI_METHOD_COUNT
        ON {}.METHOD_COUNT (METHOD_AC) NOLOGGING
        """.format(owner)
    )
    logger.debug("METHOD_COUNT ready")

    orautils.drop_table(cur, owner, "METHOD_COMPARISON")
    cur.execute(
        """
        CREATE TABLE {}.METHOD_COMPARISON
        (
            METHOD_AC1 VARCHAR2(25) NOT NULL,
            METHOD_AC2 VARCHAR2(25) NOT NULL,
            PROTEIN_COUNT NUMBER NOT NULL,
            MATCH_COUNT NUMBER NOT NULL,
            RESIDUE_COUNT NUMBER NOT NULL
        ) NOLOGGING
        """.format(owner)
    )
    table = orautils.TablePopulator(con,
                                    query="INSERT /*+ APPEND */ "
                                          "INTO {}.METHOD_COMPARISON "
                                          "VALUES (:1, :2, :3, :4, :5)".format(owner),
                                    autocommit=True)

    for acc_1 in collocs:
        for acc_2, num_collocs in collocs[acc_1].items():
            num_matches = match_ovrlps[acc_1][acc_2]
            num_residues = residue_ovrlps[acc_1][acc_2]
            table.insert((acc_1, acc_2, num_collocs, num_matches, num_residues))

    table.close()
    orautils.gather_stats(cur, owner, "METHOD_COMPARISON")
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
    logger.debug("METHOD_COMPARISON ready")
    cur.close()
    con.close()


def _load_description_counts(user: str, dsn: str, organisers: list):
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

    _acc = None
    _descriptions = {}
    for acc, descriptions in heapq.merge(*organisers):
        if acc != _acc:
            for descid, counts in _descriptions.items():
                table.insert((_acc, descid, counts['S'], counts['T']))

            _acc = acc
            _descriptions = {}

        for descid, dbcode in descriptions:
            if descid in _descriptions:
                d = _descriptions[descid]
            else:
                d = _descriptions[descid] = {'S': 0, 'T': 0}

            d[dbcode] += 1

    for descid, counts in _descriptions.items():
        table.insert((_acc, descid, counts['S'], counts['T']))

    table.close()

    for o in organisers:
        o.remove()

    cur = con.cursor()
    orautils.gather_stats(cur, owner, "METHOD_DESC")
    cur.close()
    con.close()
    logger.debug("METHOD_DESC ready")


def _load_taxonomy_counts(user: str, dsn: str, organisers: list):
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

    ranks = ["superkingdom", "kingdom", "phylum", "class", "order",
             "family", "genus", "species"]

    # Get lineages for the METHOD_TAXA table
    cur.execute(
        """
        SELECT LEFT_NUMBER, TAX_ID, RANK
        FROM {}.LINEAGE
        WHERE RANK IN ({})
        """.format(owner, ','.join(':'+str(i+1) for i in range(len(ranks)))),
        ranks
    )
    lineages = {}
    for left_num, tax_id, rank in cur:
        if left_num in lineages:
            lineages[left_num][rank] = tax_id
        else:
            lineages[left_num] = {rank: tax_id}

    table = orautils.TablePopulator(con,
                                    query="INSERT /*+ APPEND */ "
                                          "INTO {}.METHOD_TAXA "
                                          "VALUES (:1, :2, :3, :4)".format(owner),
                                    autocommit=True)

    _acc = None
    _ranks = {}
    for acc, left_numbers in heapq.merge(*organisers):
        if acc != _acc:
            for rank in ranks:
                for tax_id in _ranks.get(rank, {}):
                    table.insert((_acc, rank, tax_id, _ranks[tax_id]))

            _acc = acc
            _ranks = {}

        for left_num in left_numbers:
            for rank, tax_id in lineages.get(left_num, {}).items():
                if rank in _ranks:
                    r = _ranks[rank]
                else:
                    r = _ranks[rank] = {}

                if tax_id in r:
                    r[tax_id] += 1
                else:
                    r[tax_id] = 1

    table.close()

    for o in organisers:
        o.remove()

    cur = con.cursor()
    orautils.gather_stats(cur, owner, "METHOD_TAXA")
    cur.close()
    con.close()
    logger.debug("METHOD_TAXA ready")


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
