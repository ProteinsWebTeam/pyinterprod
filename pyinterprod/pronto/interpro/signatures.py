# -*- coding: utf-8 -*-

import json
from multiprocessing import Process, Queue
from typing import Callable, List, Optional

from .utils import merge_comparators, merge_kvdbs, Kvdb
from ... import logger


def _process2(database: str, dir: Optional[str], task_queue: Queue, done_queue: Queue):
    with Kvdb(database) as src, Kvdb(dir=dir) as dst:
        for acc_1, values_1 in iter(task_queue.get, None):
            comparisons = {}
            counts = {k: len(values) for k, values in values_1.items()}
            it = src.range(acc_1)
            next(it)
            for acc_2, values_2 in it:
                comparisons[acc_2] = {}
                for k, values in values_1.items():
                    if k in values_2:
                        comparisons[acc_2][k] = len(values & values_2[k])
                    else:
                        comparisons[acc_2][k] = 0

            dst[acc_1] = (counts, comparisons)

    done_queue.put(dst.filepath)


def _process(database: str, dir: Optional[str], task_queue: Queue, done_queue: Queue):
    with Kvdb(database) as src, Kvdb(dir=dir) as dst:
        for acc_1, values_1 in iter(task_queue.get, None):
            comparisons = {}
            it = src.range(acc_1)
            next(it)
            for acc_2, values_2 in it:
                comparisons[acc_2] = len(values_1 & values_2)

            dst[acc_1] = (len(values_1), comparisons)

    done_queue.put(dst)


def _compare(database: str, processes: int, dir: Optional[str], fn: Callable) -> Kvdb:
    task_queue = Queue(maxsize=1)
    done_queue = Queue()
    pool = []
    for _ in range(max(1, processes-1)):
        p = Process(target=fn, args=(database, dir, task_queue, done_queue))
        p.start()
        pool.append(p)

    with Kvdb(database) as kvdb:
        for acc_1, values_1 in kvdb:
            task_queue.put((acc_1, values_1))

        size = kvdb.size
        #kvdb.remove()

    for _ in pool:
        task_queue.put(None)

    kvdbs = []
    for _ in pool:
        kvdb = done_queue.get()
        kvdbs.append(kvdb)
        size += kvdb.size

    for p in pool:
        p.join()

    with Kvdb(dir=dir) as kvdb:
        for key, values in merge_kvdbs(kvdbs, remove=False):
            kvdb[key], = values  # one single item

    return kvdb


def compare_descriptions(database: str, processes: int, dir: Optional[str]) -> Kvdb:
    return _compare(database, processes, dir, _process)


def compare_taxa(database: str, processes: int, dir: Optional[str]) -> Kvdb:
    return _compare(database, processes, dir, _process2)


def compare_terms(database: str, processes: int, dir: Optional[str]) -> Kvdb:
    return _compare(database, processes, dir, _process)


def merge_comparisons(user: str, dsn: str, kvdbs: tuple, comparators: list,
                      ranks: List[str], **kwargs):
    logger.debug("collecting descriptions")
    de_kvdb = compare_descriptions(user, dsn, kvdbs[0], **kwargs)

    logger.debug("collecting terms")
    te_kvdb = compare_terms(user, dsn, kvdbs[1], **kwargs)

    logger.debug("collecting ranks")
    ta_kvdb = compare_taxa(user, dsn, kvdbs[2], ranks, **kwargs)

    logger.debug("collecting overlaps")
    signatures, comparisons = merge_comparators(comparators, remove=False)

    for acc_1, (n_seq, n_res) in signatures.items():
        try:
            n_desc, d_comp = d_kvdb[acc_1]
        except KeyError:
            n_desc = 0
            d_comp = {}

        try:
            n_term, te_comp = te_kvdb[acc_1]
        except KeyError:
            n_term = 0
            te_comp = {}

        try:
            n_ranks, ta_comp = ta_kvdb[acc_1]
        except KeyError:
            n_ranks = [0] * len(ranks)
            ta_comp = {}
        finally:
            n_ranks = dict((zip(ranks, n_ranks)))

        row1 = (acc_1, n_seq, n_res, n_desc, json.dumps(n_ranks), n_term)
        rows = []
        comp = comparisons[acc_1]
        for acc_2, (n_colloc, n_prot_over, n_res_over) in comp.items():
            n_desc = d_comp.get(acc_2, 0)
            n_rank = dict(zip(ranks, ta_comp.get(acc_2, [0] * len(ranks))))
            n_term = te_comp.get(acc_2, 0)
            rows.append((acc_1, acc_2, n_colloc, n_prot_over, n_res, n_desc,
                         json.dumps(n_rank), n_term))

        yield row1, rows

    for kvdb in (de_kvdb, te_kvdb, ta_kvdb):
        kvdb.remove()
