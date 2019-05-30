# -*- coding: utf-8 -*-

from multiprocessing import Process, Queue
from typing import Dict, List, Optional, Tuple

import cx_Oracle

from .utils import Organizer, merge_organizers
from ... import logger, orautils


def _cmp_descriptions(keys: List[str], task_queue: Queue,
                      done_queue: Queue, dir: Optional[str]=None):
    o = Organizer(keys, dir=dir)
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
    task_queue = Queue(1)
    done_queue = Queue()
    for _ in range(max(1, processes-1)):
        p = Process(target=_cmp_descriptions,
                    args=(keys, task_queue, done_queue, dir))
        p.start()
        pool.append(p)

    cur.execute(
        """
        SELECT MD.METHOD_AC, MD.DESC_ID
        FROM {0}.METHOD_DESC MD
        INNER JOIN {0}.DESC_VALUE DV
          ON MD.DESC_ID = DV.DESC_ID
        WHERE DV.TEXT NOT LIKE 'Predicted protein%'
          AND DV.TEXT NOT LIKE 'Uncharacterized protein%'
        ORDER BY MD.DESC_ID
        """.format(owner)
    )

    signatures = {}
    accessions = []
    _descid = None
    for acc, descid in cur:
        if acc in signatures:
            signatures[acc] += 1
        else:
            signatures[acc] = 1

        if descid != _descid:
            task_queue.put(accessions)
            _descid = descid
            accessions = []

        accessions.append(acc)

    task_queue.put(accessions)
    accessions = []

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

    organizer = Organizer(keys, dir=dir)
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
    o = Organizer(keys, dir=dir)
    for tax_id, accessions in iter(task_queue.get, None):
        accessions.sort()
        for i, acc_1 in enumerate(accessions):
            for acc_2 in accessions[i+1:]:
                o.add(acc_1, (acc_2, tax_id))

        o.flush()

    size = o.merge()
    done_queue.put((o, size))


def get_taxa(user: str, dsn: str, processes: int=4, bucket_size: int=20,
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
    task_queue = Queue(1)
    done_queue = Queue()
    for _ in range(max(1, processes-1)):
        p = Process(target=_cmp_taxa,
                    args=(keys, task_queue, done_queue, dir))
        p.start()
        pool.append(p)

    cur.execute(
        """
        SELECT METHOD_AC, TAX_ID
        FROM {}.METHOD_TAXA_TMP
        ORDER BY TAX_ID
        """.format(owner)
    )
    signatures = {}
    accessions = []
    n = -1
    _tax_id = None
    for acc, tax_id in cur:
        if acc in signatures:
            if tax_id in signatures[acc]:
                signatures[acc][tax_id] += 1
            else:
                signatures[acc][tax_id] = 1
        else:
            signatures[acc] = {tax_id: 1}

        if tax_id != _tax_id:
            task_queue.put((_tax_id, accessions))
            n += 1
            logger.debug("{:>10}\t{:<20}".format(n, _tax_id))
            _tax_id = tax_id
            accessions = []

        accessions.append(acc)

    task_queue.put((_tax_id, accessions))
    n += 1
    logger.debug("{:>10}\t{:<20}".format(n, _tax_id))
    accessions = []
    cur.close()
    con.close()

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

    organizer = Organizer(keys, dir=dir)
    for acc_1, items in merge_organizers(organizers):
        counts = {}
        for acc_2, tax_id in items:
            if acc_2 in counts:
                if tax_id in counts[acc_2]:
                    counts[acc_2][tax_id] += 1
                else:
                    counts[acc_2][tax_id] = 1
            else:
                counts[acc_2] = {tax_id: 1}

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
