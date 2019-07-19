# -*- coding: utf-8 -*-

import time
from multiprocessing import Process, Queue
from typing import Dict, List, Optional, Tuple

import cx_Oracle

from .. import logger, orautils
from .utils import merge_comparators, Kvdb, PersistentBuffer

JACCARD_THRESHOLD = 0.5
PROGRESS_SECONDS = 3600


def load_comparators(user: str, dsn: str, comparators: Optional[list]=None,
                     remove: bool=True):
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    orautils.drop_table(cur, owner, "METHOD_SIMILARITY", purge=True)
    cur.execute(
        """
        CREATE TABLE {}.METHOD_SIMILARITY
        (
            METHOD_AC1 VARCHAR2(25) NOT NULL,
            METHOD_AC2 VARCHAR2(25) NOT NULL,
            COLL_INDEX BINARY_DOUBLE DEFAULT NULL,
            COLL_CONT1 BINARY_DOUBLE DEFAULT NULL,
            COLL_CONT2 BINARY_DOUBLE DEFAULT NULL,
            POVR_INDEX BINARY_DOUBLE DEFAULT NULL,
            POVR_CONT1 BINARY_DOUBLE DEFAULT NULL,
            POVR_CONT2 BINARY_DOUBLE DEFAULT NULL,
            ROVR_INDEX BINARY_DOUBLE DEFAULT NULL,
            ROVR_CONT1 BINARY_DOUBLE DEFAULT NULL,
            ROVR_CONT2 BINARY_DOUBLE DEFAULT NULL,
            DESC_INDEX BINARY_DOUBLE DEFAULT NULL,
            DESC_CONT1 BINARY_DOUBLE DEFAULT NULL,
            DESC_CONT2 BINARY_DOUBLE DEFAULT NULL,
            TAXA_INDEX BINARY_DOUBLE DEFAULT NULL,
            TAXA_CONT1 BINARY_DOUBLE DEFAULT NULL,
            TAXA_CONT2 BINARY_DOUBLE DEFAULT NULL,
            TERM_INDEX BINARY_DOUBLE DEFAULT NULL,
            TERM_CONT1 BINARY_DOUBLE DEFAULT NULL,
            TERM_CONT2 BINARY_DOUBLE DEFAULT NULL,
            CONSTRAINT PK_METHOD_SIMILARITY PRIMARY KEY (METHOD_AC1, METHOD_AC2)
        ) NOLOGGING
        """.format(owner)
    )
    cur.close()

    if comparators:
        counts, comparisons = merge_comparators(comparators, remove=remove)

        table = orautils.TablePopulator(
            con=con,
            query="""
                    INSERT /*+ APPEND */ INTO {}.METHOD_SIMILARITY (
                      METHOD_AC1, METHOD_AC2, COLL_INDEX, COLL_CONT1, COLL_CONT2,
                      POVR_INDEX, POVR_CONT1, POVR_CONT2,
                      ROVR_INDEX, ROVR_CONT1, ROVR_CONT2
                    )
                    VALUES (:1, :2, :3, :4, :5, :6, :7, :8, :9, :10, :11)
                """.format(owner),
            autocommit=True
        )
        for acc1, (n_prot1, n_res1) in counts.items():
            for acc2 in comparisons[acc1]:
                n_prot2, n_res2 = counts[acc2]
                n_col, n_prot_over, n_res_over = comparisons[acc1][acc2]
                table.insert((
                    acc1, acc2,
                    # Collocation
                    n_col / (n_prot1 + n_prot2 - n_col),
                    n_col / n_prot1,
                    n_col / n_prot2,
                    # Protein overlap
                    n_prot_over / (n_prot1 + n_prot2 - n_prot_over),
                    n_prot_over / n_prot1,
                    n_prot_over / n_prot2,
                    # Residue overlap
                    n_res_over / (n_res1 + n_res2 - n_res_over),
                    n_res_over / n_res1,
                    n_res_over / n_res2,
                ))

        table.close()
        con.commit()

    cur = con.cursor()
    orautils.grant(cur, owner, "METHOD_SIMILARITY", "SELECT", "INTERPRO_SELECT")
    cur.execute(
        """
        CREATE INDEX I_METHOD_SIMILARITY$AC1
        ON {}.METHOD_SIMILARITY (METHOD_AC1) NOLOGGING
        """.format(owner)
    )
    cur.execute(
        """
        CREATE INDEX I_METHOD_SIMILARITY$AC2
        ON {}.METHOD_SIMILARITY (METHOD_AC2) NOLOGGING
        """.format(owner)
    )
    orautils.gather_stats(cur, owner, "METHOD_SIMILARITY")
    cur.close()
    con.close()


def _process(kvdb: Kvdb, task_queue: Queue, done_queue: Queue,
                dir: Optional[str]):
    signatures = {}
    with PersistentBuffer(dir=dir) as buffer:
        for acc_1, values_1 in iter(task_queue.get, None):
            counts = {}
            gen = kvdb.range(acc_1)
            next(gen)
            for acc_2, values_2 in gen:
                counts[acc_2] = len(values_1 & values_2)

            buffer.add((acc_1, counts))
            signatures[acc_1] = len(values_1)

    done_queue.put((signatures, buffer))


def _calc_similarity(counts: Dict[str, int], src: PersistentBuffer,
                     queue: Queue, dir: Optional[str]=None):
    with PersistentBuffer(dir=dir) as dst:
        for acc1, cmps in src:
            cnt1 = counts[acc1]
            val = {}
            for acc2, intersect in cmps.items():
                cnt2 = counts[acc2]
                # J(A, B) = intersection(A, B) / union(A, B)
                idx = intersect / (cnt1 + cnt2 - intersect)
                # C(A, B) = intersection(A, B) / size(A)
                #       --> larger values indicate more of A lying in B
                ct1 = intersect / cnt1
                ct2 = intersect / cnt2
                if any([sim >= JACCARD_THRESHOLD for sim in (idx, ct1, ct2)]):
                    val[acc2] = (idx, ct1, ct2)

            if val:
                dst.add((acc1, val))

    queue.put(dst)


def _compare(kvdb: Kvdb, processes: int, dir: Optional[str]) -> Tuple[List[PersistentBuffer], int]:
    logger.debug("comparing")
    pool = []
    task_queue = Queue(maxsize=1)
    done_queue = Queue()
    for _ in range(max(1, processes-1)):
        p = Process(target=_process, args=(kvdb, task_queue, done_queue, dir))
        p.start()
        pool.append(p)

    s = len(kvdb)               # matrix of shape (s, s)
    n = (pow(s, 2) - s) // 2    # number of items (half-matrix - diagonal)
    c = 0                       # number of items enqueued
    ts1 = time.time()
    for i, (acc_1, taxids_1) in enumerate(kvdb):
        task_queue.put((acc_1, taxids_1))
        c += s - (i + 1)
        ts2 = time.time()
        if ts2 > ts1 + PROGRESS_SECONDS:
            ts1 = ts2
            logger.debug(f"{c/n*100:>3.0f}%")

    for _ in pool:
        task_queue.put(None)

    buffers = []
    signatures = {}
    for _ in pool:
        _signatures, buffer = done_queue.get()
        signatures.update(_signatures)
        buffers.append(buffer)

    for p in pool:
        p.join()

    logger.debug("calculating similarities")
    pool = []
    for buffer in buffers:
        p = Process(target=_calc_similarity,
                    args=(signatures, buffer, done_queue, dir))
        p.start()
        pool.append(p)

    buffers2 = []
    for _ in pool:
        buffers2.append(done_queue.get())

    for p in pool:
        p.join()

    size = 0
    for buffer in buffers:
        size += buffer.size
        buffer.remove()

    return buffers2, size


def _export_signatures(cur: cx_Oracle.Cursor, dir: Optional[str]) -> Kvdb:
    logger.debug("exporting")
    with Kvdb(dir=dir, insertonly=True) as kvdb:
        values = set()
        _acc = None
        for acc, val in cur:
            if acc != _acc:
                if _acc:
                    kvdb[_acc] = values

                _acc = acc
                values = set()

            values.add(val)

        if _acc:
            kvdb[_acc] = values

    return kvdb


def _load_comparisons(user: str, dsn: str, column: str, buffers: List[PersistentBuffer], remove: bool=True):
    logger.debug(f"updating METHOD_SIMILARITY ({column})")
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    table = orautils.TablePopulator(
        con=con,
        query="""
                MERGE INTO {0}.METHOD_SIMILARITY
                USING DUAL ON (METHOD_AC1 = :ac1 AND METHOD_AC2 = :ac2)
                WHEN MATCHED THEN
                  UPDATE SET {1}_INDEX=:idx,
                             {1}_CONT1=:ct1,
                             {1}_CONT2=:ct2
                WHEN NOT MATCHED THEN
                  INSERT (
                    METHOD_AC1, METHOD_AC2, {1}_INDEX, {1}_CONT1, {1}_CONT2
                  ) VALUES (:ac1, :ac2, :idx, :ct1, :ct2)
            """.format(owner, column)
    )

    ts1 = time.time()
    i = 0
    for buffer in buffers:
        for acc1, cmps in buffer:
            for acc2, (idx, ct1, ct2) in cmps.items():
                table.insert({
                    "ac1": acc1,
                    "ac2": acc2,
                    "idx": idx,
                    "ct1": ct1,
                    "ct2": ct2
                })
                    
            i += 1
            ts2 = time.time()
            if ts2 > ts1 + PROGRESS_SECONDS:
                ts1 = ts2
                logger.debug(f"{i:>10}")

        if remove:
            buffer.remove()

    table.close()
    con.commit()
    con.close()


def cmp_descriptions(user: str, dsn: str, processes: int=1,
                     dir: Optional[str]=None, remove: bool=True):
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    cur.execute(
        """
            SELECT METHOD_AC, DESC_ID
            FROM {0}.METHOD_DESC
            WHERE DESC_ID NOT IN (
              SELECT DESC_ID
              FROM {0}.DESC_VALUE
              WHERE TEXT IN ('Uncharacterized protein', 'Predicted protein')
            )
            ORDER BY METHOD_AC
        """.format(owner)
    )
    kvdb = _export_signatures(cur, dir)
    cur.close()
    con.close()
    buffers, size = _compare(kvdb, processes, dir)
    size += kvdb.size
    kvdb.remove()
    return buffers, size


def cmp_taxa(user: str, dsn: str, processes: int=1,
             dir: Optional[str]=None, remove: bool=True):
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    cur.execute(
        """
        SELECT METHOD_AC, TAX_ID
        FROM {}.METHOD_TAXA
        WHERE RANK IN ('superkingdom','kingdom','class','family','species')
        ORDER BY METHOD_AC
        """.format(owner)
    )
    kvdb = _export_signatures(cur, dir)
    cur.close()
    con.close()
    buffers, size = _compare(kvdb, processes, dir)
    size += kvdb.size
    kvdb.remove()
    return buffers, size


def cmp_terms(user: str, dsn: str, processes: int=1,
              dir: Optional[str]=None, remove: bool=True):
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    cur.execute(
        """
        SELECT METHOD_AC, GO_ID
        FROM {0}.METHOD_TERM
        WHERE GO_ID NOT IN (
          SELECT GO_ID
          FROM {0}.TERM
          WHERE NAME IN (
            'protein binding', 'molecular_function', 'biological_process',
            'cellular_component'
          )
        )
        ORDER BY METHOD_AC
        """.format(owner)
    )
    kvdb = _export_signatures(cur, dir)
    cur.close()
    con.close()
    buffers, size = _compare(kvdb, processes, dir)
    size += kvdb.size
    kvdb.remove()
    return buffers, size
