# -*- coding: utf-8 -*-

import json
import time
from multiprocessing import Process, Queue
from typing import Callable, Dict, List, Optional, Tuple

import cx_Oracle

from .utils import merge_comparators, merge_kvdbs, merge_organizers, Kvdb, Organizer
from ... import logger, orautils


# def compare(keys: List[str], task_queue: Queue, done_queue: Queue,
#             dir: Optional[str]=None, buffer_size: int=0):
#     organizer = Organizer(keys, dir=dir, buffer_size=buffer_size)
#     signatures = {}
#     for values in iter(task_queue.get, None):
#         accessions = sorted({acc for val in values for acc in val})
#         for i, acc_1 in enumerate(accessions):
#             if acc_1 in signatures:
#                 signatures[acc_1] += 1
#             else:
#                 signatures[acc_1] = 1
#
#             for acc_2 in accessions[i:]:
#                 organizer.add(acc_1, acc_2)
#
#     size = organizer.merge(fn=count_accessions)
#     done_queue.put((signatures, organizer, size))


def compare(task_queue: Queue, done_queue: Queue, dir: Optional[str]=None,
            sync_freq: int=1):
    with Kvdb(dir=dir, buffer_size=-1) as kvdb:
        n_iter = 0
        for values in iter(task_queue.get, None):
            accessions = sorted({acc for val in values for acc in val})
            for i, acc_1 in enumerate(accessions):
                try:
                    count, comparisons = kvdb[acc_1]
                except KeyError:
                    count = 0
                    comparisons = {}

                count += 1
                for acc_2 in accessions[i:]:
                    try:
                        comparisons[acc_2] += 1
                    except KeyError:
                        comparisons[acc_2] = 1

                kvdb[acc_1] = (count, comparisons)
            n_iter += 1
            if not n_iter %  sync_freq:
                kvdb.sync()

    done_queue.put(kvdb)


def count_accessions(values: List[str]) -> List[Tuple[str, int]]:
    counts = {}
    for acc in values:
        if acc in counts:
            counts[acc] += 1
        else:
            counts[acc] = 1
    return list(counts.items())


# def collect_counts(pool: List[Process], queue: Queue,
#                    dir: Optional[str]=None) -> Tuple[Kvdb, int]:
#     signatures = {}
#     organizers = []
#     size = 0
#     for _ in pool:
#         counts, organizer, _size = queue.get()
#         organizers.append(organizer)
#         size += _size
#         for acc, cnt in counts.items():
#             if acc in signatures:
#                 signatures[acc] += cnt
#             else:
#                 signatures[acc] = cnt
#     for p in pool:
#         p.join()
#
#     with Kvdb(dir=dir) as kvdb:
#         for acc_1, cmp in merge_organizers(organizers, remove=True):
#             comparisons = {}
#             for acc_2, cnt in cmp:
#                 if acc_2 in comparisons:
#                     comparisons[acc_2] += cnt
#                 else:
#                     comparisons[acc_2] = cnt
#
#             kvdb[acc_1] = (signatures[acc_1], comparisons)
#
#     size += kvdb.size
#     return kvdb, size


def collect_counts(pool: List[Process], queue: Queue,
                   dir: Optional[str]=None) -> Tuple[Kvdb, int]:
    kvdbs = []
    for _ in pool:
        kvdbs.append(queue.get())
    for p in pool:
        p.join()

    size = sum([kvdb.size for kvdb in kvdbs])
    with Kvdb(dir=dir) as kvdb:
        for acc_1, items in merge_kvdbs(kvdbs, remove=True):
            count = 0
            comparisons = {}
            for cnt, cmps in items:
                count += cnt
                for acc_2, cnt in cmps.items():
                    try:
                        comparisons[acc_2] += cnt
                    except KeyError:
                        comparisons[acc_2] = cnt

            kvdb[acc_1] = (count, comparisons)

    size += kvdb.size
    return kvdb, size


def compare2(keys: List[str], ranks: List[str], task_queue: Queue,
             done_queue: Queue, dir: Optional[str]=None, buffer_size: int=0):
    organizer = Organizer(keys, dir=dir, buffer_size=buffer_size)
    signatures = {}
    for ranks_to_incr, values in iter(task_queue.get, None):
        incr = [1 if r in ranks_to_incr else 0 for r in ranks]
        accessions = sorted({acc for val in values for acc in val})
        for i, acc_1 in enumerate(accessions):
            if acc_1 in signatures:
                for x, c in enumerate(incr):
                    signatures[acc_1][x] += c
            else:
                signatures[acc_1] = incr

            for acc_2 in accessions[i:]:
                organizer.add(acc_1, (acc_2, incr))

    size = organizer.merge(fn=count_accessions_per_rank)
    done_queue.put((signatures, organizer, size))


def count_accessions_per_rank(values: List[Tuple[str, list]]) -> List[Tuple[str, list]]:
    counts = {}
    for acc, incr in values:
        if acc in counts:
            for i, c in enumerate(incr):
                counts[acc][i] += c
        else:
            counts[acc] = incr
    return list(counts.items())


def collect_counts2(pool: List[Process], queue: Queue,
                    dir: Optional[str]=None) -> Tuple[Kvdb, int]:
    logger.debug("\tgathering")
    signatures = {}
    organizers = []
    size = 0
    for _ in pool:
        counts, organizer, _size = queue.get()
        organizers.append(organizer)
        size += _size
        for acc, incr in counts.items():
            if acc in signatures:
                for i, c in enumerate(incr):
                    signatures[acc][i] += c
            else:
                signatures[acc] = incr
    for p in pool:
        p.join()

    logger.debug("\tdisk space: {} MB".format(size/1024**2))
    with Kvdb(dir=dir) as kvdb:
        for acc_1, cmp in merge_organizers(organizers, remove=True):
            comparisons = {}
            for acc_2, incr in cmp:
                if acc_2 in comparisons:
                    for i, c in enumerate(incr):
                        comparisons[acc_2][i] += c
                else:
                    comparisons[acc_2] = incr

            kvdb[acc_1] = (signatures[acc_1], comparisons)

    size += kvdb.size
    return kvdb, size


def init_pool(size: int, func: Callable, args: Tuple) -> List[Process]:
    pool = []
    for _ in range(max(1, size-1)):
        p = Process(target=func, args=args)
        p.start()
        pool.append(p)

    return pool


def chunk_accessions(user: str, dsn: str, bucket_size: int) -> List[str]:
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    cur.execute("SELECT METHOD_AC FROM {}.METHOD".format(owner))
    keys = sorted([row[0] for row in cur])
    keys = [keys[i] for i in range(0, len(keys), bucket_size)]
    cur.close()
    con.close()
    return keys


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


def compare_descriptions(user: str, dsn: str, kvdbs: List[Kvdb], **kwargs) -> Kvdb:
    bucket_size = kwargs.get("bucket_size", 100)
    buffer_size = kwargs.get("buffer_size", 0)
    dir = kwargs.get("dir")
    processes = kwargs.get("processes", 1)

    keys = chunk_accessions(user, dsn, bucket_size)
    task_queue = Queue(maxsize=1)
    done_queue = Queue()
    pool = init_pool(processes, compare, (keys, task_queue, done_queue, dir, buffer_size))
    i = 0
    for desc_id, values in merge_kvdbs(kvdbs, remove=False):
        task_queue.put(values)
        i += 1
        if not i % 100000:
            logger.debug("{:>15}".format(i))
    for _ in pool:
        task_queue.put(None)
    return collect_counts(pool, done_queue, dir=dir)


# def compare_terms(user: str, dsn: str, kvdbs: List[Kvdb], **kwargs) -> Kvdb:
#     bucket_size = kwargs.get("bucket_size", 100)
#     buffer_size = kwargs.get("buffer_size", 0)
#     dir = kwargs.get("dir")
#     processes = kwargs.get("processes", 1)
#
#     keys = chunk_accessions(user, dsn, bucket_size)
#     task_queue = Queue(maxsize=1)
#     done_queue = Queue()
#     pool = init_pool(processes, compare, (keys, task_queue, done_queue, dir, buffer_size))
#     i = 0
#     for desc_id, values in merge_kvdbs(kvdbs, remove=False):
#         task_queue.put(values)
#         i += 1
#         if not i % 1000:
#             logger.debug("{:>15}".format(i))
#     for _ in pool:
#         task_queue.put(None)
#     return collect_counts(pool, done_queue, dir=dir)


def compare_terms(user: str, dsn: str, kvdbs: List[Kvdb], **kwargs) -> Tuple[Kvdb, int]:
    dir = kwargs.get("dir")
    max_items = kwargs.get("max_items", 0)
    processes = kwargs.get("processes", 1)
    task_queue = Queue(maxsize=1)
    done_queue = Queue()
    pool = init_pool(processes, compare, (task_queue, done_queue, dir, max_items))
    i = 0
    for go_id, values in merge_kvdbs(kvdbs, remove=False):
        task_queue.put(values)
        i += 1
        if not i % 1000:
            logger.debug("{:>15}".format(i))
    for _ in pool:
        task_queue.put(None)
    return collect_counts(pool, done_queue, dir=dir)


def compare_taxa(user: str, dsn: str, kvdbs: List[Kvdb], ranks: List[str], **kwargs) -> Kvdb:
    bucket_size = kwargs.get("bucket_size", 100)
    buffer_size = kwargs.get("buffer_size", 0)
    dir = kwargs.get("dir")
    processes = kwargs.get("processes", 1)

    keys = chunk_accessions(user, dsn, bucket_size)
    task_queue = Queue(maxsize=1)
    done_queue = Queue()
    pool = init_pool(processes, compare2, (keys, ranks, task_queue, done_queue, dir, buffer_size))
    taxa = get_lineages(user, dsn)
    i = 0
    for tax_id, values in merge_kvdbs(kvdbs, remove=False):
        try:
            ranks = taxa[tax_id]
        except KeyError:
            pass  # ETAXI is incomplete
        else:
            task_queue.put((ranks, values))
            i += 1
            if not i % 1000000:
                logger.debug("{:>15}".format(i))
    for _ in pool:
        task_queue.put(None)
    return collect_counts2(ranks, pool, done_queue, dir=dir)


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
