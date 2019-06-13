# -*- coding: utf-8 -*-

import json
import time
from multiprocessing import Process, Queue
from typing import Callable, List, Optional, Tuple

import cx_Oracle

from .utils import merge_comparators, merge_kvdbs, merge_organizers, Kvdb, Organizer
from ... import logger, orautils


def compare(keys: List[str], task_queue: Queue, done_queue: Queue, dir: Optional[str]=None):
    organizer = Organizer(keys, dir=dir)
    signatures = {}
    for values in iter(task_queue.get, None):
        accessions = sorted({acc for val in values for acc in val})
        for i, acc_1 in enumerate(accessions):
            if acc_1 in signatures:
                signatures[acc_1] += 1
            else:
                signatures[acc_1] = 1

            for acc_2 in accessions[i:]:
                organizer.add(acc_1, acc_2)

        organizer.flush()
    organizer.merge()
    done_queue.put((signatures, organizer))


def collect_counts(pool: List[Process], queue: Queue, dir: Optional[str]=None) -> Kvdb:
    signatures = {}
    organizers = []
    logger.debug("\tgathering results")
    for _ in pool:
        counts, organizer = queue.get()
        organizers.append(organizer)
        for acc, cnt in counts.items():
            if acc in signatures:
                signatures[acc] += cnt
            else:
                signatures[acc]  = cnt
    for p in pool:
        p.join()

    logger.debug("\tstoring counts")
    with Kvdb(dir=dir) as kvdb:
        for acc_1, accessions in merge_organizers(organizers):
            comparisons = {}
            for acc_2 in accessions:
                if acc_2 in comparisons:
                    comparisons[acc_2] += 1
                else:
                    comparisons[acc_2] = 1

            kvdb[acc_1] = (signatures[acc_1], comparisons)

    # for o in organizers:
    #     o.remove()

    return kvdb


def compare2(keys: List[str], ranks: List[str], task_queue: Queue, done_queue: Queue, dir: Optional[str]=None):
    organizer = Organizer(keys, dir=dir)
    signatures = {}
    indices = {rank: i for i, rank in enumerate(ranks)}
    for ranks, values in iter(task_queue.get, None):
        ranks = [indices[rank] for rank in ranks]
        accessions = sorted({acc for val in values for acc in val})
        for i, acc_1 in enumerate(accessions):
            if acc_1 not in signatures:
                signatures[acc_1] = [0] * len(indices)

            for x in ranks:
                signatures[acc_1][x] += 1

            for acc_2 in accessions[i:]:
                organizer.add(acc_1, (acc_2, ranks))

        organizer.flush()
    organizer.merge()
    done_queue.put((signatures, organizer))


def collect_counts2(pool: List[Process], queue: Queue, dir: Optional[str]=None) -> Kvdb:
    signatures = {}
    organizers = []
    num_items = None
    logger.debug("\tgathering results")
    for _ in pool:
        counts, organizer = queue.get()
        organizers.append(organizer)
        for acc, val in counts.items():
            if acc in signatures:
                for i, v in enumerate(val):
                    signatures[acc][i] += v
            else:
                signatures[acc] = val
                num_items = len(val)  # we need to know the number of ranks
    for p in pool:
        p.join()

    logger.debug("\tstoring counts")
    with Kvdb(dir=dir) as kvdb:
        for acc_1, items in merge_organizers(organizers):
            comparisons = {}
            for acc_2, indices in items:
                if acc_2 not in comparisons:
                    comparisons[acc_2] = [0] * num_items

                for x in indices:
                    comparisons[acc_2][x] += 1

            kvdb[acc_1] = (signatures[acc_1], comparisons)

    # for o in organizers:
    #     o.remove()

    return kvdb


def init_pool(size: int, func: Callable, args: Tuple) -> List[Process]:
    pool = []
    for _ in range(max(1, size-1)):
        p = Process(target=func, args=args)
        p.start()
        pool.append(p)

    return pool


def merge_comparisons(user: str, dsn: str, comparators: list, kvdbs: tuple,
                      ranks: List[str], **kwargs):
    buffer_size = kwargs.get("buffer_size", 0)
    bucket_size = kwargs.get("bucket_size", 100)
    dir = kwargs.get("dir")
    processes = kwargs.get("processes", 1)

    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    cur.execute("SELECT METHOD_AC FROM {}.METHOD".format(owner))
    keys = sorted([row[0] for row in cur])
    keys = [keys[i] for i in range(0, len(keys), bucket_size)]
    cur.close()
    con.close()

    task_queue = Queue(maxsize=1)
    done_queue = Queue()

    logger.debug("collecting descriptions")
    pool = init_pool(processes, compare, (keys, task_queue, done_queue, dir))
    i = 0
    for desc_id, values in merge_kvdbs(kvdbs[0]):
        task_queue.put(values)
        i += 1
        if not i % 1000:
            logger.debug("{:>15}".format(i))
    for _ in pool:
        task_queue.put(None)
    d_kvdb = collect_counts(pool, done_queue, dir=dir)

    logger.debug("collecting terms")
    pool = init_pool(processes, compare, (keys, task_queue, done_queue, dir))
    for go_id, values in merge_kvdbs(kvdbs[1]):
        task_queue.put(values)
    for _ in pool:
        task_queue.put(None)
    te_kvdb = collect_counts(pool, done_queue, dir=dir)

    logger.debug("collecting ranks")
    pool = init_pool(processes, compare2, (keys, ranks, task_queue, done_queue, dir))
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
    i = 0
    for tax_id, values in merge_kvdbs(kvdbs[2]):
        try:
            ranks = taxa[tax_id]
        except KeyError:
            continue  # if it happens, ETAXI is incomplete
        task_queue.put((ranks, values))
        i += 1
        if not i % 1000:
            logger.debug("{:>15}".format(i))
    for _ in pool:
        task_queue.put(None)
    ta_kvdb = collect_counts2(ranks, pool, done_queue, dir=dir)

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
