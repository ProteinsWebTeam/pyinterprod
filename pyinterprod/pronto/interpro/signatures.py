# -*- coding: utf-8 -*-

from multiprocessing import Process, Queue

from .utils import merge_comparators, merge_kvdbs


def agg_kvdb(task_queue: Queue, done_queue: Queue):
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

    done_queue.put((signatures, comparisons))


def agg_kvdb2(list_size: int, task_queue: Queue, done_queue: Queue):
    signatures = {}
    comparisons = {}
    for rank, values in iter(task_queue.get, None):
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


def load_comparisons(user: str, dsn: str, comparators: list, kvdbs: tuple, processes: int=1):
    # m_signatures, m_comparisons = merge_comparators(comparators)

    # kvdbs: descriptions, taxa, terms
    d_signatures = {}
    d_comparisons = {}
    for desc_id, values in merge_kvdbs(kvdbs[0]):



    return