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
        for acc_1 in accessions:
            for acc_2 in accessions:
                if acc_1 < acc_2:
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
            for acc_1 in accessions:
                for acc_2 in accessions:
                    if acc_1 < acc_2:
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
