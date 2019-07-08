# -*- coding: utf-8 -*-

import math
from multiprocessing import Process, Queue
from typing import Dict, List, Optional, Tuple

import cx_Oracle

from .. import logger, orautils
from .utils import merge_comparators, Kvdb, PersistentBuffer


def load_comparators(user: str, dsn: str, comparators: list):
    counts, comparisons = merge_comparators(comparators, remove=True)
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    orautils.drop_table(cur, owner, "METHOD_SIMILARITY")
    cur.execute(
        """
        CREATE TABLE {}.METHOD_SIMILARITY
        (
            METHOD_AC1 VARCHAR2(25) NOT NULL,
            METHOD_AC2 VARCHAR2(25) NOT NULL,
            DESC_INDEX BINARY_DOUBLE DEFAULT NULL,
            DESC_CONT1 BINARY_DOUBLE DEFAULT NULL,
            DESC_CONT2 BINARY_DOUBLE DEFAULT NULL,
            TAXA_INDEX BINARY_DOUBLE DEFAULT NULL,
            TAXA_CONT1 BINARY_DOUBLE DEFAULT NULL,
            TAXA_CONT2 BINARY_DOUBLE DEFAULT NULL,
            TERM_INDEX BINARY_DOUBLE DEFAULT NULL,
            TERM_CONT1 BINARY_DOUBLE DEFAULT NULL,
            TERM_CONT2 BINARY_DOUBLE DEFAULT NULL,
            COLL_INDEX BINARY_DOUBLE DEFAULT NULL,
            COLL_CONT1 BINARY_DOUBLE DEFAULT NULL,
            COLL_CONT2 BINARY_DOUBLE DEFAULT NULL,
            POVR_INDEX BINARY_DOUBLE DEFAULT NULL,
            POVR_CONT1 BINARY_DOUBLE DEFAULT NULL,
            POVR_CONT2 BINARY_DOUBLE DEFAULT NULL,
            ROVR_INDEX BINARY_DOUBLE DEFAULT NULL,
            ROVR_CONT1 BINARY_DOUBLE DEFAULT NULL,
            ROVR_CONT2 BINARY_DOUBLE DEFAULT NULL,
            CONSTRAINT PK_METHOD_SIMILARITY PRIMARY KEY (METHOD_AC1, METHOD_AC2)
        ) NOLOGGING
        """.format(owner)
    )
    cur.close()

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
    con.close()

    orautils.grant(cur, owner, "METHOD_SIMILARITY", "SELECT", "INTERPRO_SELECT")
    orautils.gather_stats(cur, owner, "METHOD2PROTEIN")


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


def _compare(kvdb: Kvdb, processes: int, dir: Optional[str]) -> Tuple[Dict[str, int], List[PersistentBuffer]]:
    logger.debug("compare")
    pool = []
    task_queue = Queue(maxsize=1)
    done_queue = Queue()
    for _ in range(max(1, processes-1)):
        p = Process(target=_process, args=(kvdb, task_queue, done_queue, dir))
        p.start()
        pool.append(p)

    s = len(kvdb)               # matrix of shape (s, s)
    n = (pow(s, 2) - s) // 2    # number of items (half-matrix - diagonal)
    c = pc = _pc = 0            # c: number of items done, pc/_pc: percent
    for i, (acc_1, taxids_1) in enumerate(kvdb):
        task_queue.put((acc_1, taxids_1))
        c += s - (i + 1)
        pc = math.floor(c*100/n*10) / 10
        if pc > _pc:
            logger.debug(f"{i+1:>9}/{s:}{pc:>10}%")
            _pc = pc

    logger.debug(f"{i+1:>9}/{s:}{pc:>10}%")

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

    return signatures, buffers


def export_signatures(cur: cx_Oracle.Cursor, dir: Optional[str]) -> Kvdb:
    logger.debug("export")
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


def _load_comparisons(user: str, dsn: str, column: str, counts: Dict[str: int],
                      buffers: List[PersistentBuffer], remove: bool=True):
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    table = orautils.TablePopulator(
        con=con,
        query="""
                MERGE INTO {0}.METHOD_SIMILARITY
                USING DUAL ON (METHOD_AC1 = :ac1 AND METHOD_AC2 = :ac2)
                WHEN MATCHED THEN 
                  UPDATE SET TERM_INDEX=:idx, 
                             TERM_CONT1=:ct1, 
                             TERM_CONT2=:ct2
                WHEN NOT MATCHED THEN 
                  INSERT (
                    METHOD_AC1, METHOD_AC2, {1}_INDEX, {1}_CONT1, {1}_CONT2
                  ) VALUES (:ac1, :ac2, :idx, :ct1, :ct2)
            """.format(owner, column)
    )

    for buffer in buffers:
        for acc1, cmps in buffer:
            cnt1 = counts[acc1]
            for acc2, i in cmps.items():
                cnt2 = counts[acc2]
                table.insert({
                    "ac1": acc1,
                    "ac2": acc2,
                    "idx": i / (cnt1 + cnt2 - i),
                    # J(A, B) = intersection(A, B) / union(A, B)
                    "ct1": i / cnt1,
                    # C(A, B) = intersection(A, B) / size(A) --> larger values indicate more of A lying in B
                    "ct2": i / cnt2
                })

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
    kvdb = export_signatures(cur, dir)
    cur.close()
    con.close()
    signatures, buffers = _compare(kvdb, processes, dir)
    size = kvdb.size + sum([b.size for b in buffers])
    kvdb.remove()
    if remove:
        _load_comparisons(user, dsn, "DESC", signatures, buffers)
        return size
    else:
        return signatures, buffers, size


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
    kvdb = export_signatures(cur, dir)
    cur.close()
    con.close()
    signatures, buffers = _compare(kvdb, processes, dir)
    size = kvdb.size + sum([b.size for b in buffers])
    kvdb.remove()
    if remove:
        _load_comparisons(user, dsn, "TAXA", signatures, buffers)
        return size
    else:
        return signatures, buffers, size


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
    kvdb = export_signatures(cur, dir)
    cur.close()
    con.close()
    signatures, buffers = _compare(kvdb, processes, dir)
    size = kvdb.size + sum([b.size for b in buffers])
    kvdb.remove()
    if remove:
        _load_comparisons(user, dsn, "TERM", signatures, buffers)
        return size
    else:
        return signatures, buffers, size
