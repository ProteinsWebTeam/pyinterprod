# -*- coding: utf-8 -*-

import os
import shutil
import time
from multiprocessing import Process, Queue
from tempfile import mkstemp
from typing import Dict, Optional

import cx_Oracle
from mundone import Task

from .. import logger, orautils
from .utils import merge_comparators, Kvdb, PersistentBuffer

JACCARD_THRESHOLD = 0.5
PROGRESS_SECONDS = 3600


def load_comparators(user: str, dsn: str, comparators: list):
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    orautils.drop_table(cur, owner, "METHOD_SIMILARITY", purge=True)
    cur.execute(
        """
        CREATE TABLE {}.METHOD_SIMILARITY
        (
            METHOD_AC1 VARCHAR2(25) NOT NULL,
            METHOD_AC2 VARCHAR2(25) NOT NULL,
            COLL_INDEX BINARY_DOUBLE DEFAULT NULL,
            COLL_CONT1 BINARY_DOUBLE DEFAULT NULL,
            COLL_CONT2 BINARY_DOUBLE DEFAULT NULL,
            POVR_INDEX BINARY_DOUBLE DEFAULT NULL,
            POVR_CONT1 BINARY_DOUBLE DEFAULT NULL,
            POVR_CONT2 BINARY_DOUBLE DEFAULT NULL,
            ROVR_INDEX BINARY_DOUBLE DEFAULT NULL,
            ROVR_CONT1 BINARY_DOUBLE DEFAULT NULL,
            ROVR_CONT2 BINARY_DOUBLE DEFAULT NULL,
            DESC_INDEX BINARY_DOUBLE DEFAULT NULL,
            DESC_CONT1 BINARY_DOUBLE DEFAULT NULL,
            DESC_CONT2 BINARY_DOUBLE DEFAULT NULL,
            TAXA_INDEX BINARY_DOUBLE DEFAULT NULL,
            TAXA_CONT1 BINARY_DOUBLE DEFAULT NULL,
            TAXA_CONT2 BINARY_DOUBLE DEFAULT NULL,
            TERM_INDEX BINARY_DOUBLE DEFAULT NULL,
            TERM_CONT1 BINARY_DOUBLE DEFAULT NULL,
            TERM_CONT2 BINARY_DOUBLE DEFAULT NULL,
            CONSTRAINT PK_METHOD_SIMILARITY PRIMARY KEY (METHOD_AC1, METHOD_AC2)
        ) NOLOGGING
        """.format(owner)
    )
    cur.close()

    similarities, num_full_sequences = merge_comparators(comparators)
    table = orautils.TablePopulator(
        con=con,
        query="""
                INSERT /*+ APPEND */ INTO {}.METHOD_SIMILARITY (
                  METHOD_AC1, METHOD_AC2, COLL_INDEX, COLL_CONT1, COLL_CONT2,
                  POVR_INDEX, POVR_CONT1, POVR_CONT2,
                  ROVR_INDEX, ROVR_CONT1, ROVR_CONT2
                )
                VALUES (:1, :2, :3, :4, :5, :6, :7, :8, :9, :10, :11)
            """.format(owner),
        autocommit=True
    )
    for acc1 in similarities:
        for acc2, values in similarities[acc1].items():
            table.insert((acc1, acc2) + values)
    table.close()

    table = orautils.TablePopulator(
        con=con,
        query="""
                UPDATE {}.METHOD
                SET FULL_SEQ_COUNT = :1
                WHERE METHOD_AC = :2
            """.format(owner),
        autocommit=True
    )
    for acc, n in num_full_sequences.items():
        table.update((n, acc))
    table.close()
    con.commit()
    con.close()


def _process(kvdb: Kvdb, start: str, stop: str, task_queue: Queue,
             done_queue: Queue, dir: Optional[str]=None):
    with PersistentBuffer(dir=dir) as buffer:
        for acc_1, values_1 in iter(task_queue.get, None):
            intersections = {}
            for acc_2, values_2 in kvdb.range(start, stop):
                if acc_1 < acc_2:
                    intersections[acc_2] = len(values_1 & values_2)

            buffer.add((acc_1, intersections))

    done_queue.put(buffer)


def _load_similarities(user: str, dsn: str, column: str, queue: Queue):
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    table = orautils.TablePopulator(
        con=con,
        query="""
              UPDATE {0}.METHOD_SIMILARITY
              SET {1}_INDEX=:idx, {1}_CONT1=:ct1, {1}_CONT2=:ct2
              WHERE METHOD_AC1 = :ac1 AND METHOD_AC2 = :ac2
              """.format(owner, column),
        autocommit=True
    )

    with PersistentBuffer() as buffer:
        buffer.close()
        buffer.remove()

        for filename in iter(queue.get, None):
            buffer.filename = filename

            for similarities in buffer:
                for acc_1, obj in similarities.items():
                    for acc_2, (idx, ct1, ct2) in obj.items():
                        table.update({
                            "ac1": acc_1,
                            "ac2": acc_2,
                            "idx": idx,
                            "ct1": ct1,
                            "ct2": ct2
                        })

            buffer.remove()

    table.close()
    con.close()


def _calc_similarities(buffer: PersistentBuffer, counts: Dict[str, int],
                       queue: Queue):
    results = {}
    for acc1, cmps in buffer:
        cnt1 = counts[acc1]
        val = {}
        for acc2, intersect in cmps.items():
            cnt2 = counts[acc2]
            # J(A, B) = intersection(A, B) / union(A, B)
            idx = intersect / (cnt1 + cnt2 - intersect)
            # C(A, B) = intersection(A, B) / size(A)
            #       --> larger values indicate more of A lying in B
            ct1 = intersect / cnt1
            ct2 = intersect / cnt2
            if any([sim >= JACCARD_THRESHOLD for sim in (idx, ct1, ct2)]):
                val[acc2] = (idx, ct1, ct2)

        if val:
            results[acc1] = val

    buffer.remove()
    queue.put(results)


def _compare(kvdb_src: str, i_start: str, i_stop: str, j_start: str,
             j_stop: str, outdir: str, processes: int=8,
             tmpdir: Optional[str]=None):
    # Copy Kvdb file locally
    if tmpdir:
        os.makedirs(tmpdir, exist_ok=True)
    fd, kvdb_dst = mkstemp(dir=tmpdir)
    os.close(fd)
    os.remove(kvdb_dst)
    shutil.copy(kvdb_src, kvdb_dst)

    counts = {}  # required counts to compute similarities
    pool = []
    task_queue = Queue(maxsize=1)
    done_queue = Queue()
    with Kvdb(kvdb_dst) as kvdb:
        for _ in range(max(1, processes-1)):
            p = Process(target=_process,
                        args=(kvdb, j_start, j_stop, task_queue, done_queue,
                              tmpdir))
            p.start()
            pool.append(p)

        for acc_1, taxids_1 in kvdb.range(i_start, i_stop):
            task_queue.put((acc_1, taxids_1))
            counts[acc_1] = len(taxids_1)

        if i_start != j_start:
            # Get remaining required counts
            for acc_1, taxids_1 in kvdb.range(j_start, j_stop):
                counts[acc_1] = len(taxids_1)

    for _ in pool:
        task_queue.put(None)

    # We don't need the local Kvdb any more
    size = kvdb.size
    kvdb.remove()

    # Temporary buffers containing intersections
    tmp_buffers = []
    for _ in pool:
        buffer = done_queue.get()
        tmp_buffers.append(buffer)
        size += buffer.size

    for p in pool:
        p.join()

    pool = []
    for buffer in tmp_buffers:
        # Temporary buffers are deleted in _calc_similarities()
        p = Process(target=_calc_similarities,
                    args=(buffer, counts, done_queue))
        p.start()
        pool.append(p)

    with PersistentBuffer(dir=outdir) as buffer:
        for _ in pool:
            buffer.add(done_queue.get())

    for p in pool:
        p.join()

    logger.info(f"disk usage: {size/1024**2:.0f}MB")
    return buffer.filename


def _export_signatures(cur: cx_Oracle.Cursor, dst: str):
    with Kvdb(insertonly=True) as kvdb:
        values = set()
        _acc = None
        for acc, val in cur:
            if acc != _acc:
                if _acc:
                    kvdb[_acc] = values

                _acc = acc
                values = set()

            values.add(val)

        if _acc:
            kvdb[_acc] = values

    shutil.copy(kvdb.filepath, dst)
    kvdb.remove()


def _chunk_jobs(cur: cx_Oracle.Cursor, schema: str, chunk_size: int):
    cur.execute(
        f"""
        SELECT METHOD_AC
        FROM {schema}.METHOD
        ORDER BY METHOD_AC
        """
    )
    accessions = [row[0] for row in cur]
    chunks = []
    for i in range(0, len(accessions), chunk_size):
        start = accessions[i]
        try:
            stop = accessions[i+chunk_size-1]
        except IndexError:
            stop = accessions[-1]

        chunks.append((start, stop))

    jobs = []
    for i, (row_start, row_stop) in enumerate(chunks):
        for col_start, col_stop in chunks[i:]:
            jobs.append((row_start, row_stop, col_start, col_stop))

    return jobs


def _run_comparisons(user: str, dsn: str, query: str, column: str,
                     outdir: str, processes: int, chunk_size: int,
                     max_jobs: int, job_processes: int,
                     job_tmpdir: Optional[str], job_queue: Optional[str]):
    os.makedirs(outdir, exist_ok=True)
    fd, kvdb_path = mkstemp(suffix=".db", dir=outdir)
    os.close(fd)
    os.remove(kvdb_path)

    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    logger.info(f"{column}: exporting")
    cur.execute(query)
    _export_signatures(cur, kvdb_path)
    logger.info(f"{column}:     {os.path.getsize(kvdb_path)/1024**2:.0f}MB")

    queue = Queue()
    loaders = []
    for _ in range(max(1, processes-1)):
        p = Process(target=_load_similarities,
                    args=(user, dsn, column, queue))
        p.start()
        loaders.append(p)

    logger.info(f"{column}: comparing")
    pending = _chunk_jobs(cur, owner, chunk_size)
    cur.close()
    con.close()

    running = []
    submitted = 0
    done = 0
    failed = 0
    while pending or running:
        if max_jobs and pending:
            for _ in range(max_jobs - len(running)):
                # Can submit a new job
                row_start, row_stop, col_start, col_stop = pending.pop()
                t = Task(
                    name=f"pronto-cmp-{submitted}",
                    fn=_compare,
                    args=(kvdb_path, row_start, row_stop, col_start, col_stop,
                          outdir, job_processes, job_tmpdir),
                    scheduler=dict(queue=job_queue, cpu=job_processes,
                                   mem=4000, scratch=5000)
                )
                t.run(workdir=outdir)
                running.append(t)
                submitted += 1
        elif pending:
            # Submit all jobs at once
            for row_start, row_stop, col_start, col_stop in pending:
                t = Task(
                    name=f"pronto-cmp-{submitted}",
                    fn=_compare,
                    args=(kvdb_path, row_start, row_stop, col_start, col_stop,
                          outdir, job_processes, job_tmpdir),
                    scheduler=dict(queue=job_queue, cpu=job_processes,
                                   mem=4000, scratch=5000)
                )
                t.run(workdir=outdir)
                running.append(t)
                submitted += 1
            pending = []

        # Check completed tasks
        _running = []
        for task in running:
            if task.done():
                logger.debug(f"{task.stdout}\n{task.stderr}")
                if task.successful():
                    # todo: change this when updating `mundone`
                    queue.put(task.output.read())
                    done += 1
                else:
                    pending = []  # cancel pending tasks
                    failed += 1

                logger.info(f"{column}:     {done} done; {failed} failed")
            else:
                _running.append(task)

        running = _running
        time.sleep(10)

    for _ in loaders:
        queue.put(None)

    os.remove(kvdb_path)

    for p in loaders:
        p.join()

    if failed:
        raise RuntimeError(f"{failed} task(s) failed")

    logger.info(f"{column}: complete")


def cmp_descriptions(user: str, dsn: str, outdir: str, processes: int=4,
                     chunk_size: int=10000, max_jobs: int=0,
                     job_processes: int=8, job_tmpdir: Optional[str]=None,
                     job_queue: Optional[str]=None):
    query = """
        SELECT METHOD_AC, DESC_ID
        FROM {0}.METHOD_DESC
        WHERE DESC_ID NOT IN (
          SELECT DESC_ID
          FROM {0}.DESC_VALUE
          WHERE TEXT IN ('Uncharacterized protein', 'Predicted protein')
        )
        ORDER BY METHOD_AC
    """.format(user.split('/')[0])
    _run_comparisons(user, dsn, query, "DESC", outdir, processes, chunk_size,
                     max_jobs, job_processes, job_tmpdir, job_queue)


def cmp_taxa(user: str, dsn: str, outdir: str, processes: int=4,
             chunk_size: int=10000, max_jobs: int=0,
             job_processes: int=8, job_tmpdir: Optional[str]=None,
             job_queue: Optional[str]=None):
    query = """
        SELECT METHOD_AC, TAX_ID
        FROM {}.METHOD_TAXA
        ORDER BY METHOD_AC
    """.format(user.split('/')[0])
    _run_comparisons(user, dsn, query, "TAXA", outdir, processes, chunk_size,
                     max_jobs, job_processes, job_tmpdir, job_queue)


def cmp_terms(user: str, dsn: str, outdir: str, processes: int=4,
              chunk_size: int=10000, max_jobs: int=0,
              job_processes: int=8, job_tmpdir: Optional[str]=None,
              job_queue: Optional[str]=None):
    query = """
        SELECT METHOD_AC, GO_ID
        FROM {0}.METHOD_TERM
        WHERE GO_ID NOT IN (
          SELECT GO_ID
          FROM {0}.TERM
          WHERE NAME IN (
            'protein binding', 'molecular_function', 'biological_process',
            'cellular_component'
          )
        )
        ORDER BY METHOD_AC
    """.format(user.split('/')[0])
    _run_comparisons(user, dsn, query, "TERM", outdir, processes, chunk_size,
                     max_jobs, job_processes, job_tmpdir, job_queue)


def compare(user: str, dsn: str, outdir: str, processes: int=4,
            chunk_size: int=10000, max_jobs: int=0, job_processes: int=8,
            job_tmpdir: Optional[str]=None, job_queue: Optional[str]=None):
    cmp_terms(user, dsn, outdir, processes, chunk_size, max_jobs,
              job_processes, job_tmpdir, job_queue)

    cmp_descriptions(user, dsn, outdir, processes, chunk_size, max_jobs,
                     job_processes, job_tmpdir, job_queue)

    cmp_taxa(user, dsn, outdir, processes, chunk_size, max_jobs,
             job_processes, job_tmpdir, job_queue)

    logger.info("indexing/anayzing METHOD_SIMILARITY")
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    orautils.grant(cur, owner, "METHOD_SIMILARITY", "SELECT", "INTERPRO_SELECT")
    orautils.drop_index(cur, owner, "I_METHOD_SIMILARITY$AC1")
    orautils.drop_index(cur, owner, "I_METHOD_SIMILARITY$AC2")
    cur.execute(
        """
        CREATE INDEX I_METHOD_SIMILARITY$AC1
        ON {}.METHOD_SIMILARITY (METHOD_AC1) NOLOGGING
        """.format(owner)
    )
    cur.execute(
        """
        CREATE INDEX I_METHOD_SIMILARITY$AC2
        ON {}.METHOD_SIMILARITY (METHOD_AC2) NOLOGGING
        """.format(owner)
    )
    orautils.gather_stats(cur, owner, "METHOD_SIMILARITY")
    cur.close()
    con.close()
    logger.info("METHOD_SIMILARITY is ready")
