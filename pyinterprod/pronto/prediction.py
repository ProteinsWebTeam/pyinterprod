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
from .utils import Kvdb, PersistentBuffer as Buffer


PROGRESS_SECONDS = 3600
SIMILARITY_THRESHOLD = 0.75


def _compare_chunk(database: str, start1: str, stop1: str, start2: str,
                   stop2: str, dst: str, tmpdir: Optional[str]):
    fd, _database = mkstemp(dir=tmpdir)
    os.close(fd)
    os.remove(_database)
    shutil.copy(database, _database)

    counts = {}
    with Buffer(dir=tmpdir) as tmp, Buffer(filepath=dst) as buffer:
        with Kvdb(_database) as kvdb:
            for acc_1, values_1 in kvdb.range(start1, stop1):
                counts[acc_1] = len(values_1)
                intersections = {}

                for acc_2, values_2 in kvdb.range(start2, stop2):
                    if acc_1 == acc_2:
                        continue  # skip diagonal

                    intersection = len(values_1 & values_2)
                    if intersection:
                        intersections[acc_2] = intersection

                if intersections:
                    tmp.add((acc_1, intersections))

            if start1 != start2:
                for acc_1, values_1 in kvdb.range(start2, stop2):
                    counts[acc_1] = len(values_1)

        size = kvdb.size + tmp.size
        kvdb.remove()

        for acc_1, intersections in tmp:
            cnt_1 = counts[acc_1]
            similarities = {}

            for acc_2, intersection in intersections.items():
                cnt_2 = counts[acc_2]

                # J(A, B) = intersection(A, B) / union(A, B)
                idx = intersection / (cnt_1 + cnt_2 - intersection)

                # C(A, B) = intersection(A, B) / size(A)
                #       --> larger values indicate more of A lying in B
                ct1 = intersection / cnt_1
                ct2 = intersection / cnt_2

                similarities[acc_2] = (idx, ct1, ct2)

            buffer.add((acc_1, similarities))

        tmp.remove()

    logger.info(f"disk usage: {size/1024**2:.0f}MB")


def _process(kvdb: Kvdb, start: str, stop: str, task_queue: Queue,
             done_queue: Queue, dir: Optional[str]=None):
    with PersistentBuffer(dir=dir) as buffer:
        for acc_1, values_1 in iter(task_queue.get, None):
            intersections = {}
            for acc_2, values_2 in kvdb.range(start, stop):
                if acc_1 < acc_2:
                    intersection = len(values_1 & values_2)
                    if intersection:
                        intersections[acc_2] = intersection

            if intersections:
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


def _export_signatures(cur: cx_Oracle.Cursor, database: str):
    with Kvdb(database, insertonly=True) as kvdb:
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

    for i, (row_start, row_stop) in enumerate(chunks):
        for col_start, col_stop in chunks[i:]:
            yield row_start, row_stop, col_start, col_stop


def _run_comparisons(user: str, dsn: str, query: str, column: str,
                     outdir: str, chunk_size: int, job_processes: int,
                     job_tmpdir: Optional[str], job_queue: Optional[str]):
    os.makedirs(job_tmpdir, exist_ok=True)
    fd, kvdb_tmp = mkstemp(suffix=".db", dir=job_tmpdir)
    os.close(fd)
    os.remove(kvdb_tmp)

    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    logger.info(f"{column}: exporting")
    cur.execute(query)
    _export_signatures(cur, kvdb_tmp)
    logger.info(f"{column}:     {os.path.getsize(kvdb_tmp)/1024**2:.0f}MB")

    os.makedirs(outdir, exist_ok=True)
    kvdb_path = os.path.join(outdir, os.path.basename(kvdb_tmp))
    try:
        os.remove(kvdb_path)
    except FileNotFoundError:
        pass
    finally:
        shutil.move(kvdb_tmp, kvdb_path)

    queue = Queue()
    loader = Process(target=_load_similarities,
                     args=(user, dsn, column, queue))
    loader.start()

    # Submit all jobs at once
    logger.info(f"{column}: comparing")
    running = []
    submitted = 0
    for start1, stop1, start2, stop2 in _chunk_jobs(cur, owner, chunk_size):
        t = Task(
            name=f"pronto-cmp-{submitted}",
            fn=_compare,
            args=(kvdb_path, start1, stop1, start2, stop2,
                  outdir, job_processes, job_tmpdir),
            scheduler=dict(queue=job_queue, cpu=job_processes,
                           mem=4000, scratch=5000)
        )
        t.run(workdir=outdir)
        running.append(t)
        submitted += 1

    cur.close()
    con.close()

    done = 0
    failed = 0
    while running:
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
                    failed += 1
            else:
                _running.append(task)

        running = _running
        time.sleep(10)

    queue.put(None)
    os.remove(kvdb_path)
    loader.join()

    if failed:
        raise RuntimeError(f"{failed} task(s) failed")

    logger.info(f"{column}: complete")


def _predict(user: str, dsn: str, queue: str, num_producers: int):



def cmp_descriptions(user: str, dsn: str, outdir: str, chunk_size: int=10000,
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
    _run_comparisons(user, dsn, query, "DESC", outdir, chunk_size,
                     job_processes, job_tmpdir, job_queue)


def cmp_taxa(user: str, dsn: str, outdir: str, chunk_size: int=10000,
             job_processes: int=8, job_tmpdir: Optional[str]=None,
             job_queue: Optional[str]=None):
    query = """
        SELECT METHOD_AC, TAX_ID
        FROM {}.METHOD_TAXA
        ORDER BY METHOD_AC
    """.format(user.split('/')[0])
    _run_comparisons(user, dsn, query, "TAXA", outdir, chunk_size,
                     job_processes, job_tmpdir, job_queue)


def cmp_terms(user: str, dsn: str, outdir: str, chunk_size: int=10000,
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
    _run_comparisons(user, dsn, query, "TERM", outdir, chunk_size,
                     job_processes, job_tmpdir, job_queue)


def _compare_signatures(user: str, dsn: str, query: str, source: str,
                        outdir: str, queue: Queue, chunk_size: int,
                        tmpdir: Optional[str], job_queue: Optional[str]):
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()

    # Export signatures to local database
    fd, database = mkstemp(suffix=".db", dir=tmpdir)
    os.close(fd)
    os.remove(database)
    with Kvdb(database, insertonly=True) as kvdb:
        cur.execute(query)
        values = set()
        _acc = None
        for acc, val in cur:
            if acc != _acc:
                if _acc:
                    kvdb[_acc] = values

                _acc = acc
                values.clear()

            values.add(val)

        if _acc:
            kvdb[_acc] = values
            values.clear()

    # Move database
    dst = os.path.join(outdir, os.path.basename(database))
    try:
        os.remove(dst)
    except FileNotFoundError:
        pass
    finally:
        shutil.move(database, dst)
        database = dst

    # Submit all jobs
    owner = user.split('/')[0]
    tasks = {}
    for start1, stop1, start2, stop2 in _chunk_jobs(cur, owner, chunk_size):
        fd, filepath = mkstemp(dir=outdir)
        os.close(fd)
        t = Task(
            fn=_compare_chunk,
            args=(database, start1, stop1, start2, stop2, filepath, tmpdir),
            scheduler=dict(queue=job_queue, mem=4000)
        )
        t.run(workdir=outdir)
        tasks[t] = (start1, start2, filepath)

    cur.close()
    con.close()
    running = len(tasks)
    while running:
        for t in tasks:
            if t.done():
                logger.debug(f"{task.stdout}\n{task.stderr}")
                running -= 1
                if t.successful():
                    start1, start2, filepath = tasks[t]
                    queue.put(((start1, start2), source, filepath))

        time.sleep(10)

    queue.put(None)


def compare(user: str, dsn: str, outdir: str, chunk_size: int=10000,
            job_tmpdir: Optional[str]=None, job_queue: Optional[str]=None):
    os.makedirs(outdir, exist_ok=True)
    if job_tmpdir:
        os.makedirs(job_tmpdir, exist_ok=True)

    queue = Queue()

    owner = user.split('/')[0]
    producers = []
    query = f"""
        SELECT METHOD_AC, DESC_ID
        FROM {owner}.METHOD_DESC
        WHERE DESC_ID NOT IN (
          SELECT DESC_ID
          FROM {owner}.DESC_VALUE
          WHERE TEXT IN ('Uncharacterized protein', 'Predicted protein')
        )
        ORDER BY METHOD_AC
    """
    producers.append(Process(target=_compare_signatures,
                             args=(user, dsn, query, "desc", outdir, queue,
                                   chunk_size, job_tmpdir, job_queue)))

    query = f"""
        SELECT METHOD_AC, TAX_ID
        FROM {owner}.METHOD_TAXA
        ORDER BY METHOD_AC
    """
    producers.append(Process(target=_compare_signatures,
                             args=(user, dsn, query, "taxa", outdir, queue,
                                   chunk_size, job_tmpdir, job_queue)))

    query = f"""
        SELECT METHOD_AC, GO_ID
        FROM {owner}.METHOD_TERM
        WHERE GO_ID NOT IN (
          SELECT GO_ID
          FROM {owner}.TERM
          WHERE NAME IN (
            'protein binding', 'molecular_function', 'biological_process',
            'cellular_component'
          )
        )
        ORDER BY METHOD_AC
    """
    producers.append(Process(target=_compare_signatures,
                             args=(user, dsn, query, "term", outdir, queue,
                                   chunk_size, job_tmpdir, job_queue)))

    for p in producers:
        p.start()

    running = len(producers)
    chunks = {}
    while running:
        for key, source, filepath in iter(queue.get, None):
            if key in chunks:
                pass
            else:
                chunks[key] = {
                    "desc": None,
                    "taxa": None,
                    "term": None
                }

        running -= 1

    for p in producers:
        p.join()

    return

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
