# -*- coding: utf-8 -*-

from multiprocessing import Process, Queue
from typing import Dict, List, Optional, Tuple

import cx_Oracle

from .utils import Organizer, merge_organizers
from ... import logger, orautils


def _cmp_descriptions(keys: List[str], task_queue: Queue,
                      done_queue: Queue, dir: Optional[str]=None):
    o = Organizer(keys, dir)
    for accessions in iter(task_queue.get, None):
        accessions.sort()
        for i, acc_1 in enumerate(accessions):
            for acc_2 in accessions[i+1:]:
                o.add(acc_1, acc_2)

        o.flush()

    size = o.merge()
    done_queue.put((o, size))


def get_descriptions(user: str, dsn: str, processes: int=4,
                     bucket_size: int=100,
                     dir: Optional[str]=None) -> Tuple[Dict, Organizer]:
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    cur.execute(
        """
        SELECT METHOD_AC
        FROM {}.METHOD
        ORDER BY METHOD_AC
        """.format(owner)
    )
    keys = []
    for i, row in enumerate(cur):
        if not i % bucket_size:
            keys.append(row[0])

    pool = []
    task_queue = Queue()
    done_queue = Queue()
    for _ in range(processes):
        p = Process(target=_cmp_descriptions,
                    args=(keys, task_queue, done_queue, dir))
        p.start()
        pool.append(p)

    descriptions = {}
    signatures = {}
    cur.execute(
        """
        SELECT MD.METHOD_AC, MD.DESC_ID
        FROM {0}.METHOD_DESC MD
        INNER JOIN {0}.DESC_VALUE DV
          ON MD.DESC_ID = DV.DESC_ID
        WHERE DV.TEXT NOT LIKE 'Predicted protein%'
          AND DV.TEXT NOT LIKE 'Uncharacterized protein%'
        """.format(owner)
    )
    for acc, descid in cur:
        if descid in descriptions:
            descriptions[descid].append(acc)
        else:
            descriptions[descid] = [acc]

        if acc in signatures:
            signatures[acc] += 1
        else:
            signatures[acc] = 1

    for accessions in descriptions.values():
        task_queue.put(accessions)

    descriptions = None

    for _ in pool:
        task_queue.put(None)

    organizers = []
    size = 0
    for _ in pool:
        o, s = done_queue.get()
        organizers.append(o)
        size += s

    for p in pool:
        p.join()

    organizer = Organizer(keys, dir)
    for acc_1, accessions in merge_organizers(organizers):
        counts = {}
        for acc_2 in accessions:
            if acc_2 in counts:
                counts[acc_2] += 1
            else:
                counts[acc_2] = 1
        organizer.write(acc_1, counts)

    size += organizer.size
    logger.debug("disk usage (description): {:.0f} MB".format(size/1024**2))

    for o in organizers:
        o.remove()

    return signatures, organizer


def _cmp_taxa(keys: List[str], task_queue: Queue, done_queue: Queue,
              dir: Optional[str]=None):
    o = Organizer(keys, dir)
    for rank, taxa in iter(task_queue.get, None):
        for taxid, accessions in taxa.items():
            accessions.sort()
            for i, acc_1 in enumerate(accessions):
                for acc_2 in accessions[i+1:]:
                    o.add(acc_1, (rank, acc_2))

        o.flush()

    size = o.merge()
    done_queue.put((o, size))


def get_taxa(user: str, dsn: str, processes: int=4, bucket_size: int=100,
             dir: Optional[str]=None) -> Tuple[Dict, Organizer]:
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    cur.execute(
        """
        SELECT METHOD_AC
        FROM {}.METHOD
        ORDER BY METHOD_AC
        """.format(owner)
    )
    keys = []
    for i, row in enumerate(cur):
        if not i % bucket_size:
            keys.append(row[0])

    pool = []
    task_queue = Queue()
    done_queue = Queue()
    for _ in range(processes):
        p = Process(target=_cmp_taxa,
                    args=(keys, task_queue, done_queue, dir))
        p.start()
        pool.append(p)

    signatures = {}
    taxa = {}
    cur.execute(
        """
        SELECT METHOD_AC, RANK, TAX_ID
        FROM INTERPRO_ANALYSIS.METHOD_TAXA
        """
    )
    for acc, rank, taxid in cur:
        if rank in taxa:
            if taxid in taxa[rank]:
                taxa[rank][taxid].append(acc)
            else:
                taxa[rank][taxid] = [acc]
        else:
            taxa[rank] = {taxid: [acc]}

        if acc in signatures:
            if rank in signatures[acc]:
                signatures[acc][rank] += 1
            else:
                signatures[acc][rank] = 1
        else:
            signatures[acc] = {rank: 1}

    for rank in taxa:
        task_queue.put((rank, taxa[rank]))

    taxa = None

    for _ in pool:
        task_queue.put(None)

    organizers = []
    size = 0
    for _ in pool:
        o, s = done_queue.get()
        organizers.append(o)
        size += s

    for p in pool:
        p.join()

    organizer = Organizer(keys, dir)
    for acc_1, items in merge_organizers(organizers):
        counts = {}
        for rank, acc_2 in items:
            if rank in counts:
                if acc_2 in counts[rank]:
                    counts[rank][acc_2] += 1
                else:
                    counts[rank][acc_2] = 1
            else:
                counts[rank] = {acc_2: 1}

        organizer.write(acc_1, counts)

    size += organizer.size
    logger.debug("disk usage (taxonomy): {:.0f} MB".format(size/1024**2))

    for o in organizers:
        o.remove()

    return signatures, organizer


def get_matches(user: str, dsn: str) -> Tuple[Dict, Dict]:
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    cur.execute(
        """
        SELECT METHOD_AC, PROTEIN_COUNT, RESIDUE_COUNT
        FROM {}.METHOD_COUNT
        """.format(owner)
    )
    signatures = {}
    for row in cur:
        signatures[row[0]] = (row[1], row[2])

    cur.execute(
        """
        SELECT METHOD_AC1, METHOD_AC2, COLLOCATION, PROTEIN_OVERLAP, RESIDUE_OVERLAP
        FROM {}.METHOD_COMPARISON
        """.format(owner)
    )
    comparisons = {}
    for row in cur:
        acc_1 = row[0]
        acc_2 = row[1]
        if acc_1 in comparisons:
            comparisons[acc_1][acc_2] = row[2:]
        else:
            comparisons[acc_1] = {acc_2: row[2:]}

    cur.close()
    con.close()

    return signatures, comparisons


def get_terms(user: str, dsn: str) -> Tuple[Dict, Dict]:
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    terms = {}
    signatures = {}
    cur.execute(
        """
        SELECT METHOD_AC, GO_ID
        FROM {}.METHOD_TERM
        """.format(owner)
    )
    for acc, goid in cur:
        if goid in terms:
            terms[goid].append(acc)
        else:
            terms[goid] = [acc]

        if acc in signatures:
            signatures[acc] += 1
        else:
            signatures[acc] = 1

    comparisons = {}
    for goid, accessions in terms.items():
        accessions.sort()
        for i, acc_1 in enumerate(accessions):
            for acc_2 in accessions[i+1:]:
                if acc_1 in comparisons:
                    if acc_2 in comparisons[acc_1]:
                        comparisons[acc_1][acc_2] += 1
                    else:
                        comparisons[acc_1][acc_2] = 1
                else:
                    comparisons[acc_1] = {acc_2: 1}

    return signatures, comparisons


def make_predictions(user: str, dsn: str, processes: int=4,
                     dir: Optional[str]=None):
    sd, od = get_descriptions(user, dsn, processes=processes, dir=dir)
    st, ot = get_taxa(user, dsn, processes=processes, dir=dir)
    sm, cm = get_matches(user, dsn)
    sg, cg = get_terms(user, dsn)

    accessions = set()
    for d in (sd, st, sm, sg):
        accessions |= set(d.keys())

    accessions = sorted(accessions)
