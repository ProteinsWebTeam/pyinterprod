# -*- coding: utf-8 -*-

import json
from multiprocessing import Process, Queue
from typing import Callable, List, Optional, Tuple

from ... import logger
from .utils import Kvdb, PersistentBuffer


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


def compare(user: str, dsn: str, query: str, processes: int,
            dir: Optional[str], max_items: int, chunk_size: int):
    pool = []
    task_queue = Queue(maxsize=max(1, processes-1))
    done_queue = Queue()
    for _ in range(max(1, processes-1)):
        p = Process(target=process, args=(dir, task_queue, max_items, done_queue))
        p.start()
        pool.append(p)

    logger.debug("starting")
    con = cx_Oracle.connect(user + '@' + dsn)
    cur = con.cursor()
    cur.execute(query)
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
    cur.close()
    con.close()

    for _ in pool:
        task_queue.put(None)

    logger.debug("collecting")
    counts = {}
    buffers = []
    for i, p in enumerate(pool):
        _counts, buffer = done_queue.get()
        for acc, cnt in _counts.items():
            try:
                counts[acc] += cnt
            except KeyError:
                counts[acc] = cnt

        buffers.append(buffer)
        logger.debug("collected: {}/{} ({} chunks)".format(i+1, len(pool), len(buffer)))

    for p in pool:
        p.join()

    logger.debug("aggregating")
    keys = sorted(counts.keys())
    keys = [keys[i] for i in range(0, len(keys), chunk_size)]
    organizer = utils.Organizer(keys, dir=dir)
    size = 0
    num_items = 0
    for i, buffer in enumerate(buffers):
        for comparisons in buffer:
            for acc, cmps in comparisons.items():
                organizer.add(acc, cmps)
                num_items += len(cmps)

                if num_items >= max_items * processes:
                    organizer.flush()
                    num_items = 0

        size += buffer.size
        buffer.remove()
        logger.debug("merged: {}/{}".format(i+1, len(buffers)))

    logger.debug("merging")
    size += organizer.merge(processes, fn=sum_counts)

    logger.debug("total: {:.0f} MB".format(size/1024**2))


def cmp_descriptions(user: str, dsn: str, processes: int=1,
                     dir: Optional[str]=None, max_items: int=10000000,
                     chunk_size: int=100):
    return _compare(_process, *args, **kwargs)


def cmp_taxa(user: str, dsn: str, processes: int=1,
             dir: Optional[str]=None, max_items: int=10000000,
             chunk_size: int=100):
    return _compare(_process2, *args, **kwargs)


def cmp_terms(user: str, dsn: str, processes: int=1,
              dir: Optional[str]=None, max_items: int=10000000,
              chunk_size: int=100):
    return _compare(_process, *args, **kwargs)


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
