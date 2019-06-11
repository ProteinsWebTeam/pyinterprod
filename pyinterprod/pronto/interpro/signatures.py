# -*- coding: utf-8 -*-

import json
from multiprocessing import Process, Queue
from typing import Callable, List, Optional, Tuple

import cx_Oracle

from .utils import merge_comparators, merge_kvdbs, Kvdb
from ... import logger, orautils


def compare(task_queue: Queue, done_queue: Queue, dir: Optional[str]=None):
    signatures = {}
    comparisons = {}
    for values in iter(task_queue.get, None):
        accessions = sorted({acc for val in values for acc in val})
        for i, acc_1 in enumerate(accessions):
            if acc_1 in signatures:
                signatures[acc_1] += 1
            else:
                signatures[acc_1] = 1
                comparisons[acc_1] = {}

            for acc_2 in accessions[i:]:
                if acc_2 in comparisons[acc_1]:
                    comparisons[acc_1][acc_2] += 1
                else:
                    comparisons[acc_1][acc_2] = 1

    with Kvdb(dir=dir, cache=False) as kvdb:
        for acc, count in signatures.items():
            kvdb[acc] = (count, comparisons[acc])

    done_queue.put(kvdb)


def compare2(ranks: List[str], task_queue: Queue, done_queue: Queue, dir: Optional[str]=None):
    indices = {rank: i for i, rank in enumerate(ranks)}
    signatures = {}
    comparisons = {}
    with Kvdb(dir=dir, cache=False) as kvdb:
        for ranks, values in iter(task_queue.get, None):
            accessions = sorted({acc for val in values for acc in val})
            for j, acc_1 in enumerate(accessions):
                try:
                    counts, comparisons = kvdb[acc_1]
                except KeyError:
                    counts = [0] * len(indices)
                    comparisons = {}

                for rank in ranks:
                    i = indices[rank]
                    counts[i] += 1

                for acc_2 in accessions[j:]:
                    if acc_2 not in comparisons:
                        comparisons[acc_2] = [0] * len(indices)

                    for rank in ranks:
                        i = indices[rank]
                        comparisons[acc_2][i] += 1

                kvdb[acc_1] = (counts, comparisons)

    done_queue.put(kvdb)


def collect_counts(pool: List[Process], queue: Queue, dir: Optional[str]=None) -> Kvdb:
    kvdbs = [queue.get() for _ in pool]
    for p in pool:
        p.join()

    logger.debug("\tstoring counts")
    with Kvdb(dir=dir, cache=False) as kvdb:
        for acc_1, items in merge_kvdbs(kvdbs):
            count = None
            comp = {}
            for cnt, cmp in items:
                if count is None:
                    count = cnt
                elif isinstance(count, list):
                    for i, c in enumerate(cnt):
                        count[i] += c
                else:
                    count += cnt

                for acc_2, cnt in cmp.items():
                    if acc_2 in comp:
                        if isinstance(cnt, list):
                            for i, c in enumerate(cnt):
                                comp[acc_2][i] += c
                        else:
                            comp[acc_2] += cnt
                    else:
                        comp[acc_2] = cnt

            kvdb[acc_1] = (count, comp)

    # for db in kvdbs:
    #     db.remove()

    return kvdb


def init_pool(size: int, func: Callable, args: Tuple) -> List[Process]:
    pool = []
    for _ in range(max(1, size-1)):
        p = Process(target=func, args=args)
        p.start()
        pool.append(p)

    return pool


def merge_comparisons(user: str, dsn: str, comparators: list, kvdbs: tuple,
                      ranks: List[str], processes: int=1,
                      dir: Optional[str]=None):
    task_queue = Queue(maxsize=1)
    done_queue = Queue()

    logger.debug("collecting descriptions")
    pool = init_pool(processes, compare, (task_queue, done_queue, dir))
    for desc_id, values in merge_kvdbs(kvdbs[0]):
        task_queue.put(values)
    for _ in pool:
        task_queue.put(None)
    d_kvdb = collect_counts(pool, done_queue, dir=dir)

    logger.debug("collecting terms")
    pool = init_pool(processes, compare, (task_queue, done_queue, dir))
    for go_id, values in merge_kvdbs(kvdbs[1]):
        task_queue.put(values)
    for _ in pool:
        task_queue.put(None)
    te_kvdb = collect_counts(pool, done_queue, dir=dir)

    logger.debug("collecting ranks")
    pool = init_pool(processes, compare2, (ranks, task_queue, done_queue, dir))
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    taxa = {}
    cur.execute("SELECT TAX_ID, RANK FROM {}.LINEAGE".format(owner))
    for tax_id, rank in cur:
        key = str(tax_id)
        if key in taxa:
            taxa[key].append(rank)
        else:
            taxa[key] = [rank]
    cur.close()
    con.close()
    for tax_id, values in merge_kvdbs(kvdbs[2]):
        try:
            ranks = taxa[tax_id]
        except KeyError:
            continue  # if it happens, ETAXI is incomplete
        task_queue.put((ranks, values))
    for _ in pool:
        task_queue.put(None)
    ta_kvdb = collect_counts(pool, done_queue, dir=dir)

    signatures, comparisons = merge_comparators(comparators)

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

    # d_kvdb.remove()
    # te_kvdb.remove()
    # ta_kvdb.remove()
