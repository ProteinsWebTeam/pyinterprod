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


def collect(keys: List[str], task_queue: Queue, done_queue: Queue, dir: Optional[str]):
    organizer = Organizer(keys, dir=dir)
    for signatures, comparisons in iter(task_queue.get, None):
        for acc, cnt in signatures.items():
            organizer[acc] = (cnt, comparisons[acc])
        organizer.flush()
    size = organizer.merge(fn=sum_counts)
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
        p.close()

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


def _compare_simple(user: str, dsn: str, kvdbs: List[Kvdb], **kwargs) -> Tuple[Kvdb, int]:
    bucket_size = kwargs.get("bucket_size", 1000)
    dir = kwargs.get("dir")
    max_items = kwargs.get("max_items", 1000000)
    processes = kwargs.get("processes", 1)

    task_queue = Queue(maxsize=1)
    done_queue = Queue(maxsize=1)
    keys = chunk_accessions(user, dsn, bucket_size)
    pool = init_pool(processes-1, collect, (keys, task_queue, done_queue, dir))

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


def feed_processes(processes: int, fn: Callable, kvdbs: List[Kvdb], done_queue: Queue):
    task_queue = Queue(maxsize=1)

    # Init consumers; -2: -1 (parent proc) + -1 (this proc)
    pool = []
    for _ in range(max(1, processes-2)):
        p = Process(target=fn, args=(task_queue, done_queue))
        p.start()
        pool.append(p)

    # Submit tasks to pool of workers
    cnt = 0
    for key, values in merge_kvdbs(kvdbs, remove=False):
        task_queue.put(values)
        accessions = sorted({acc for val in values for acc in val})
        cnt += 1
        if not cnt % 1000:
            logger.debug("{:>15}".format(cnt))

    # Indicate to workers that that no more data will be sent
    for _ in pool:
        queue.put(None)

    # Join workers
    for p in pool:
        p.join()

    # Release workers
    for p in pool:
        p.close()

    # Indicate to parent process that no more data will be sent
    done_queue.put(None)


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
