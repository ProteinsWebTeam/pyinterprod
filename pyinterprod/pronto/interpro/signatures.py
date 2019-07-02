# -*- coding: utf-8 -*-

from multiprocessing import Process, Queue
from typing import List, Optional, Tuple

import cx_Oracle

from ... import logger
from .utils import Kvdb, Organizer, PersistentBuffer


def process(dir: Optional[str], task_queue: Queue, max_items: int,
            done_queue: Queue):
    with PersistentBuffer(dir) as buffer:
        counts = {}
        comparisons = {}
        num_items = 0
        for accessions in iter(task_queue.get, None):
            for i, acc_1 in enumerate(accessions):
                try:
                    counts[acc_1] += 1
                except KeyError:
                    counts[acc_1] = 1

                try:
                    cmp = comparisons[acc_1]
                except KeyError:
                    cmp = comparisons[acc_1] = {}

                for acc_2 in accessions[i+1:]:
                    try:
                        cmp[acc_2] += 1
                    except KeyError:
                        cmp[acc_2] = 1
                        num_items += 1

                if max_items and num_items >= max_items:
                    buffer.add(comparisons)
                    comparisons = {}
                    num_items = 0

        buffer.add(comparisons)

    done_queue.put((counts, buffer))


def compare(cur: cx_Oracle.Cursor, processes: int,
            dir: Optional[str], max_items: int, chunk_size: int) -> Tuple[Kvdb, int]:
    pool = []
    task_queue = Queue(maxsize=max(1, processes-1))
    done_queue = Queue()
    for _ in range(max(1, processes-1)):
        p = Process(target=process, args=(dir, task_queue, max_items, done_queue))
        p.start()
        pool.append(p)

    _key = None
    accessions = []
    i = 0
    for key, acc in cur:
        if key != _key:
            task_queue.put(accessions)
            accessions = []
            _key = key

        accessions.append(acc)
        i += 1
        if not i % 10000000:
            logger.debug("{:>12,}".format(i))

    task_queue.put(accessions)
    accessions = []
    logger.debug("{:>12,}".format(i))

    for p in pool:
        task_queue.put(None)

    logger.debug("collecting")
    counts = {}
    buffers = []
    for p in pool:
        _counts, buffer = done_queue.get()
        for acc, cnt in _counts.items():
            try:
                counts[acc] += cnt
            except KeyError:
                counts[acc] = cnt

        buffers.append(buffer)

    for p in pool:
        p.join()

    logger.debug("organizing")
    keys = sorted(counts.keys())
    keys = [keys[i] for i in range(0, len(keys), chunk_size)]
    organizer = Organizer(keys, dir=dir)
    size = 0
    for buffer in buffers:
        size += buffer.size
        for comparisons in buffer:
            for acc, cmps in comparisons.items():
                organizer.add(acc, cmps)

            organizer.flush()
        buffer.remove()

    logger.debug("merging")
    size += organizer.merge(processes, fn=sum_counts)

    with Kvdb(dir=dir) as kvdb:
        for acc, (cmps,) in organizer:
            kvdb[acc] = cmps

    organizer.remove()
    return kvdb, size


def sum_counts(values):
    counts = {}
    for cmps in values:
        for acc, cnt in cmps.items():
            if acc in counts:
                counts[acc] += cnt
            else:
                counts[acc] = cnt

    return [counts]


def cmp_descriptions(user: str, dsn: str, processes: int=1,
                     dir: Optional[str]=None, max_items: int=10000000,
                     chunk_size: int=10) -> Tuple[Kvdb, int]:
    con = cx_Oracle.connect(user + '@' + dsn)
    cur = con.cursor()
    cur.execute(
        """
            SELECT DESC_ID, METHOD_AC
            FROM {0}.METHOD_DESC
            WHERE DESC_ID NOT IN (
              SELECT DESC_ID
              FROM {0}.DESC_VALUE
              WHERE TEXT IN ('Uncharacterized protein', 'Predicted protein')
            )
            ORDER BY DESC_ID, METHOD_AC
        """.format(user.split('/')[0])
    )
    obj = compare(cur, processes, dir, max_items, chunk_size)
    cur.close()
    con.close()
    return obj


def cmp_taxa(user: str, dsn: str, processes: int=1,
             dir: Optional[str]=None, max_items: int=10000000,
             chunk_size: int=10, rank: Optional[str]=None) -> Tuple[Kvdb, int]:
    con = cx_Oracle.connect(user + '@' + dsn)
    cur = con.cursor()
    if rank:
        cur.execute(
            """
                SELECT TAX_ID, METHOD_AC
                FROM {}.METHOD_TAXA
                WHERE RANK = :1
                ORDER BY TAX_ID, METHOD_AC
            """.format(user.split('/')[0]),
            (rank,)
        )
    else:
        cur.execute(
            """
                SELECT TAX_ID, METHOD_AC
                FROM {}.METHOD_TAXA
                ORDER BY TAX_ID, METHOD_AC
            """.format(user.split('/')[0])
        )
    obj = compare(cur, processes, dir, max_items, chunk_size)
    cur.close()
    con.close()
    return obj


def cmp_terms(user: str, dsn: str, processes: int=1,
              dir: Optional[str]=None, max_items: int=10000000,
              chunk_size: int=10) -> Tuple[Kvdb, int]:
    con = cx_Oracle.connect(user + '@' + dsn)
    cur = con.cursor()
    cur.execute(
        """
        SELECT GO_ID, METHOD_AC
        FROM {0}.METHOD_TERM
        WHERE GO_ID NOT IN (
          SELECT GO_ID
          FROM {0}.TERM
          WHERE NAME IN (
            'protein binding', 'molecular_function', 'biological_process',
            'cellular_component'
          )
        )
        ORDER BY GO_ID, METHOD_AC
        """.format(user.split('/')[0])
    )
    obj = compare(cur, processes, dir, max_items, chunk_size)
    cur.close()
    con.close()
    return obj


def process_new(kvdb: Kvdb, task_queue: Queue, done_queue: Queue,
                dir: Optional[str]):
    with PersistentBuffer(dir=dir) as buffer:
        for acc_1, values_1 in iter(task_queue.get, None):
            counts = {}
            gen = kvdb.range(acc_1)
            next(gen)
            for acc_2, values_2 in gen:
                counts[acc_2] = len(values_1 & values_2)

            buffer.add((acc_1, len(values_1), counts))

    done_queue.put(buffer)


def compare_new(kvdb: Kvdb, processes: int, dir: Optional[str]):
    logger.debug("compare")
    pool = []
    task_queue = Queue(maxsize=1)
    done_queue = Queue()
    for _ in range(max(1, processes-1)):
        p = Process(target=process_new, args=(kvdb, task_queue, done_queue, dir))
        p.start()
        pool.append(p)

    i = j = 0
    s = len(kvdb)
    n = (pow(s, 2) - s) // 2
    for acc_1, taxids_1 in kvdb:
        task_queue.put((acc_1, taxids_1))
        i += 1
        j += (s - i)
        logger.debug(f"{acc_1:<30}{i:>9,}{s:,}   {j:>15,}{n:,}")

    logger.debug(f"{acc_1:<30}{i:>9,}{s:,}   {j:>15,}{n:,}")

    for _ in pool:
        task_queue.put(None)

    buffers = [done_queue.get() for _ in pool]

    for p in pool:
        p.join()

    return buffers


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


def cmp_descriptions_new(user: str, dsn: str, processes: int=1, dir: Optional[str]=None):
    con = cx_Oracle.connect(user + '@' + dsn)
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
        """.format(user.split('/')[0])
    )
    kvdb_tmp = export_signatures(cur, dir)
    cur.close()
    con.close()
    buffers = compare_new(kvdb_tmp, processes, dir)
    size = kvdb_tmp.size
    kvdb_tmp.remove()
    return buffers, size


def cmp_taxa_new(user: str, dsn: str, processes: int=1, dir: Optional[str]=None, rank: Optional[str]=None):
    con = cx_Oracle.connect(user + '@' + dsn)
    cur = con.cursor()
    if rank:
        cur.execute(
            """
                SELECT METHOD_AC, TAX_ID
                FROM {}.METHOD_TAXA
                WHERE RANK = :1
                ORDER BY METHOD_AC
            """.format(user.split('/')[0]),
            (rank,)
        )
    else:
        cur.execute(
            """
                SELECT METHOD_AC, TAX_ID
                FROM {}.METHOD_TAXA
                ORDER BY METHOD_AC
            """.format(user.split('/')[0])
        )
    kvdb_tmp = export_signatures(cur, dir)
    cur.close()
    con.close()
    buffers = compare_new(kvdb_tmp, processes, dir)
    size = kvdb_tmp.size
    kvdb_tmp.remove()
    return buffers, size


def cmp_terms_new(user: str, dsn: str, processes: int=1, dir: Optional[str]=None):
    con = cx_Oracle.connect(user + '@' + dsn)
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
        """.format(user.split('/')[0])
    )
    kvdb_tmp = export_signatures(cur, dir)
    cur.close()
    con.close()
    buffers = compare_new(kvdb_tmp, processes, dir)
    size = kvdb_tmp.size
    kvdb_tmp.remove()
    return buffers, size


# def merge_comparisons(user: str, dsn: str, kvdbs: tuple, comparators: list,
#                       ranks: List[str], **kwargs):
#     logger.debug("collecting descriptions")
#     de_kvdb = compare_descriptions(user, dsn, kvdbs[0], **kwargs)
#
#     logger.debug("collecting terms")
#     te_kvdb = compare_terms(user, dsn, kvdbs[1], **kwargs)
#
#     logger.debug("collecting ranks")
#     ta_kvdb = compare_taxa(user, dsn, kvdbs[2], ranks, **kwargs)
#
#     logger.debug("collecting overlaps")
#     signatures, comparisons = merge_comparators(comparators, remove=False)
#
#     for acc_1, (n_seq, n_res) in signatures.items():
#         try:
#             n_desc, d_comp = d_kvdb[acc_1]
#         except KeyError:
#             n_desc = 0
#             d_comp = {}
#
#         try:
#             n_term, te_comp = te_kvdb[acc_1]
#         except KeyError:
#             n_term = 0
#             te_comp = {}
#
#         try:
#             n_ranks, ta_comp = ta_kvdb[acc_1]
#         except KeyError:
#             n_ranks = [0] * len(ranks)
#             ta_comp = {}
#         finally:
#             n_ranks = dict((zip(ranks, n_ranks)))
#
#         row1 = (acc_1, n_seq, n_res, n_desc, json.dumps(n_ranks), n_term)
#         rows = []
#         comp = comparisons[acc_1]
#         for acc_2, (n_colloc, n_prot_over, n_res_over) in comp.items():
#             n_desc = d_comp.get(acc_2, 0)
#             n_rank = dict(zip(ranks, ta_comp.get(acc_2, [0] * len(ranks))))
#             n_term = te_comp.get(acc_2, 0)
#             rows.append((acc_1, acc_2, n_colloc, n_prot_over, n_res, n_desc,
#                          json.dumps(n_rank), n_term))
#
#         yield row1, rows
#
#     for kvdb in (de_kvdb, te_kvdb, ta_kvdb):
#         kvdb.remove()
