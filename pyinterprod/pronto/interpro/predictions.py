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


def _agg_taxa(user: str, dsn: str, keys: List[str], task_queue: Queue,
              dir: Optional[str]=None):
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    cur.execute("SELECT TAX_ID, RANK FROM {}.LINEAGE".format(owner))
    ranks = {}
    for tax_id, rank in cur:
        if tax_id in ranks:
            ranks[tax_id].append(rank)
        else:
            ranks[tax_id] = [rank]
    cur.close()
    con.close()

    o = Organizer(keys, dir=dir)
    for acc_1, acc_2, taxa in iter(task_queue.get, None):
        counts = {}
        for tax_id in taxa:
            for rank in ranks[tax_id]:
                if rank in counts:
                    counts[rank] += 1
                else:
                    counts[rank] = 1

        o.write(acc_1, (acc_2, counts))
        logger.debug("{:<30}\t{:<30}".format(acc_1, acc_2))


def _cmp_taxa(keys: List[str], src: str, task_queue: Queue, done_queue: Queue):
    o = Organizer(keys, exists=True, dir=src)

    for acc_1, taxa_1 in iter(task_queue.get, None):
        for acc_2, items in o:
            taxa_2 = items[0]
            done_queue.put((acc_1, acc_2, taxa_1 & taxa_2))


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

    logger.debug("exporting")
    cur.execute(
        """
        SELECT METHOD_AC, TAX_ID
        FROM INTERPRO_ANALYSIS_LOAD.METHOD_TAXA_TMP
        ORDER BY METHOD_AC
        """
    )
    _acc = None
    taxa = set()
    organizer = Organizer(keys, dir=dir)
    for acc, tax_id in cur:
        if acc != _acc:
            if _acc:
                organizer.write(_acc, taxa)

            _acc = acc
            taxa = set()

        taxa.add(tax_id)

    if _acc:
        organizer.write(_acc, taxa)
        taxa = set()

    cur.close()
    con.close()

    logger.debug("disk usage (taxonomy): {:.0f} MB".format(organizer.size/1024**2))

    pool = []
    task_queue = Queue(1)
    done_queue = Queue(1)
    for _ in range(max(1, processes - 2)):
        p = Process(target=_cmp_taxa,
                    args=(keys, organizer.dir, task_queue, done_queue))
        p.start()
        pool.append(p)

    agg = Process(target=_agg_taxa, args=(user, dsn, keys, done_queue, dir))
    agg.start()

    logger.debug("comparing")
    for acc, items in organizer:
        task_queue.put((acc, items[0]))

    for _ in pool:
        task_queue.put(None)

    for p in pool:
        p.join()

    done_queue.put(None)
    agg.join()

    signatures = {}
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
