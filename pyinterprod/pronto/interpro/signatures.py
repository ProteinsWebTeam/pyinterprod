# -*- coding: utf-8 -*-

import json
import os
import time
from multiprocessing import Process, Queue
from typing import Callable, Dict, List, Optional, Tuple

import cx_Oracle

from .utils import merge_comparators, merge_kvdbs, merge_organizers, Kvdb, Organizer
from ... import logger, orautils


def chunk_accessions(user: str, dsn: str, bucket_size: int) -> List[str]:
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    cur.execute("SELECT METHOD_AC FROM {}.METHOD".format(owner))
    keys = sorted([row[0] for row in cur])
    cur.close()
    con.close()
    return [keys[i] for i in range(0, len(keys), bucket_size)]


def collect(keys: List[str], dir: Optional[str], task_queue: Queue, fn: Callable, done_queue: Queue):
    organizer = Organizer(keys, dir=dir)
    for signatures, comparisons in iter(task_queue.get, None):
        for acc, cnt in signatures.items():
            organizer[acc] = (cnt, comparisons.pop(acc))
        organizer.flush()
    size = organizer.merge(fn=fn)
    done_queue.put((organizer, size))


def init_pool(size: int, fn: Callable, args: Tuple) -> List[Process]:
    pool = []
    for _ in range(max(1, size)):
        p = Process(target=fn, args=args)
        p.start()
        pool.append(p)

    return pool


def close_pool(pool: List[Process], inqueue: Queue, outqueue: Optional[Queue]=None):
    for _ in pool:
        inqueue.put(None)

    if outqueue:
        result = [outqueue.get() for _ in pool]
    else:
        result = None

    for p in pool:
        p.join()
        # p.close()  # Python >= 3.7

    return result


def sum_counts(values: List[Tuple[int, Dict[str, int]]]) -> Tuple[int, Dict[str, int]]:
    count = 0
    comparisons = {}
    for cnt, cmp in values:
        count += cnt
        for acc, cnt2 in cmp.items():
            if acc in comparisons:
                comparisons[acc] += cnt2
            else:
                comparisons[acc] = cnt2

    return count, comparisons


def sum_ranks(values: List[Tuple[List[int], Dict[str, List[int]]]]) -> Tuple[List[int], Dict[str, List[int]]]:
    counts = None
    comparisons = {}
    for cnts, cmp in values:
        if counts is None:
            counts = cnts
        else:
            for i, c in enumerate(cnts):
                counts[i] += c

        for acc, cnts2 in cmp.items():
            if acc in comparisons:
                for i, c in enumerate(cnts2):
                    comparisons[acc][i] += c
            else:
                comparisons[acc] = cnts2

    return counts, comparisons


def _compare_simple(user: str, dsn: str, kvdbs: List[Kvdb], **kwargs) -> Tuple[Kvdb, int]:
    bucket_size = kwargs.get("bucket_size", 1000)
    dir = kwargs.get("dir")
    max_items = kwargs.get("max_items", 1000000)
    processes = kwargs.get("processes", 1)

    task_queue = Queue(maxsize=1)
    done_queue = Queue(maxsize=1)
    keys = chunk_accessions(user, dsn, bucket_size)
    pool = init_pool(processes-1, collect, (keys, dir, task_queue, sum_counts, done_queue))

    cnt = 0
    num_items = 0
    signatures = {}
    comparisons = {}
    for key, values in merge_kvdbs(kvdbs, remove=False):
        accessions = sorted({acc for val in values for acc in val})
        for i, acc_1 in enumerate(accessions):
            if acc_1 in signatures:
                signatures[acc_1] += 1
                cmp = comparisons[acc_1]
            else:
                signatures[acc_1] = 1
                cmp = comparisons[acc_1] = {}

            for acc_2 in accessions[i:]:
                if acc_2 in cmp:
                    cmp[acc_2] += 1
                else:
                    cmp[acc_2] = 1
                    num_items += 1

            if num_items >= max_items:
                task_queue.put((signatures, comparisons))
                signatures = {}
                comparisons = {}
                logger.debug("{:<10}{:>15}: enqueued {} items ({} / {})".format(key, cnt, num_items, i, len(accessions)))
                num_items = 0

        cnt += 1
    task_queue.put((signatures, comparisons))
    organizers = []
    size = 0
    for o, s in close_pool(pool, task_queue, done_queue):
        organizers.append(o)
        size += s

    logger.debug("\tdisk space: {} MB".format(size/1024**2))
    with Kvdb(dir=dir) as kvdb:
        i = 0
        for acc, items in merge_organizers(organizers, remove=True):
            kvdb[acc_1] = sum_counts(items)

            i += 1
            if not i % 1000:
                logger.debug("{:>15}".format(i))

    return kvdb, size + kvdb.size


def compare_terms(*args, **kwargs) -> Tuple[Kvdb, int]:
    return _compare_simple(*args, **kwargs)


def compare_descriptions(*args, **kwargs) -> Tuple[Kvdb, int]:
    return _compare_simple(*args, **kwargs)


def get_lineages(user: str, dsn: str) -> Dict[str, List[str]]:
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
    return taxa


def compare_taxa(user: str, dsn: str, kvdbs: List[Kvdb], ranks: List[str], **kwargs) -> Tuple[Kvdb, int]:
    bucket_size = kwargs.get("bucket_size", 1000)
    dir = kwargs.get("dir")
    max_items = kwargs.get("max_items", 1000000)
    processes = kwargs.get("processes", 1)

    task_queue = Queue(maxsize=1)
    done_queue = Queue(maxsize=1)
    keys = chunk_accessions(user, dsn, bucket_size)
    pool = init_pool(processes-1, collect, (keys, dir, task_queue, sum_ranks, done_queue))

    taxa = get_lineages(user, dsn)

    cnt = 0
    num_items = 0
    signatures = {}
    comparisons = {}
    for key, values in merge_kvdbs(kvdbs, remove=False):
        try:
            _ranks = taxa[key]
        except KeyError:
            continue  # ETAXI is incomplete

        accessions = sorted({acc for val in values for acc in val})
        values = [1 if r in _ranks else 0 for r in ranks]
        for i, acc_1 in enumerate(accessions):
            if acc_1 in signatures:
                for x, v in enumerate(values):
                    signatures[acc_1][x] += v
                cmp = comparisons[acc_1]
            else:
                signatures[acc_1] = values
                cmp = comparisons[acc_1] = {}

            for acc_2 in accessions[i:]:
                if acc_2 in cmp:
                    for x, v in enumerate(values):
                        cmp[acc_2][x] += v
                else:
                    cmp[acc_2] = values
                    num_items += 1

            if num_items >= max_items:
                task_queue.put((signatures, comparisons))
                signatures = {}
                comparisons = {}
                logger.debug("{:<10}{:>15}: enqueued {} items ({} / {})".format(key, cnt, num_items, i, len(accessions)))
                num_items = 0

        cnt += 1
    task_queue.put((signatures, comparisons))
    organizers = []
    size = 0
    for o, s in close_pool(pool, task_queue, done_queue):
        organizers.append(o)
        size += s

    logger.debug("\tdisk space: {} MB".format(size/1024**2))
    with Kvdb(dir=dir) as kvdb:
        i = 0
        for acc, items in merge_organizers(organizers, remove=True):
            kvdb[acc_1] = sum_ranks(items)

            i += 1
            if not i % 1000:
                logger.debug("{:>15}".format(i))

    return kvdb, size + kvdb.size


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
