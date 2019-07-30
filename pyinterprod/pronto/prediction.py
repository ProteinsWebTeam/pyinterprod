# -*- coding: utf-8 -*-

import os
import shutil
import time
from multiprocessing import Process, Queue
from tempfile import mkstemp
from typing import Dict, List, Optional

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
                if acc_1 != acc_2:
                    intersections[acc_2] = len(values_1 & values_2)

            buffer.add((acc_1, intersections))

    done_queue.put(buffer)


def _calc_similarity(counts: Dict[str, int], src: PersistentBuffer,
                     queue: Queue, dir: Optional[str]=None):
    with PersistentBuffer(dir=dir) as dst:
        for acc1, cmps in src:
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
                dst.add((acc1, val))

    queue.put(dst)


def _compare(src: str, i_start: str, i_stop: str, j_start: str, j_stop: str,
             outdir: str, processes: int=8, tmpdir: Optional[str]=None):
    # Copy Kvdb file locally
    logger.info("copying")
    if tmpdir:
        os.makedirs(tmpdir, exist_ok=True)
    fd, dst = mkstemp(dir=tmpdir)
    os.close(fd)
    os.remove(dst)
    shutil.copy(src, dst)

    logger.info("comparing")
    counts = {}  # required counts to compute similarities
    pool = []
    task_queue = Queue(maxsize=1)
    done_queue = Queue()
    with Kvdb(dst) as kvdb:
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
            # Get remaning required counts
            for acc_1, taxids_1 in kvdb.range(j_start, j_stop):
                counts[acc_1] = len(taxids_1)

    for _ in pool:
        task_queue.put(None)

    # We don't need the local Kvdb any more
    size = kvdb.size
    kvdb.remove()

    # Temporary buffers (with intersections, NOT similarities)
    tmp_buffers = []
    for _ in pool:
        buffer = done_queue.get()
        tmp_buffers.append(buffer)
        size += buffer.size

    for p in pool:
        p.join()

    logger.info("calculating similarities")
    pool = []
    for buffer in tmp_buffers:
        p = Process(target=_calc_similarity,
                    args=(counts, buffer, done_queue, outdir))
        p.start()
        pool.append(p)

    # Final similarity buffers
    sim_buffers = []
    for _ in pool:
        buffer = done_queue.get()
        sim_buffers.append(buffer)

    for p in pool:
        p.join()

    for buffer in tmp_buffers:
        buffer.remove()

    logger.info(f"disk usage: {size/1024**2:.0f}MB")


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
    n = len(chunks)
    for i in range(n):
        row_start, row_stop = chunks[i]
        for j in range(i, n):
            col_start, col_stop = chunks[j]
            jobs.append((row_start, row_stop, col_start, col_stop))

    return jobs


def _load_comparisons(user: str, dsn: str, column: str, files: List[str]):
    logger.info(f"updating METHOD_SIMILARITY ({column})")
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    table = orautils.TablePopulator(
        con=con,
        query="""
              UPDATE {0}.METHOD_SIMILARITY
              SET {1}_INDEX=:idx, {1}_CONT1=:ct1, {1}_CONT2=:ct2
              WHERE METHOD_AC1 = :ac1 AND METHOD_AC2 = :ac2
              """.format(owner, column)
    )

    with PersistentBuffer() as buffer:
        # remove temp file created by new instance
        buffer.close()
        buffer.remove()

        for f in files:
            buffer.filename = f

            for acc1, cmps in buffer:
                for acc2, (idx, ct1, ct2) in cmps.items():
                    table.update({
                        "ac1": acc1,
                        "ac2": acc2,
                        "idx": idx,
                        "ct1": ct1,
                        "ct2": ct2
                    })

            buffer.remove()

    table.close()
    con.commit()
    con.close()
    logger.info(f"METHOD_SIMILARITY ({column}) updated")


def _run_comparisons(user: str, dsn: str, query: str, outdir: str,
                     processes: int=8, max_jobs: int=0,
                     tmpdir: Optional[str]=None, chunk_size: int=10000,
                     job_queue: Optional[str]=None) -> List[str]:
    os.makedirs(outdir, exist_ok=True)
    kvdb_path = os.path.join(outdir, "signatures.db")

    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    cur.execute(query)
    logger.info("exporting")
    _export_signatures(cur, kvdb_path)
    pending = _chunk_jobs(cur, owner, chunk_size)
    cur.close()
    con.close()

    if tmpdir:
        os.makedirs(tmpdir, exist_ok=True)

    logger.info("comparing")
    running = []
    submitted = 0
    failed = False
    while pending or running:
        if max_jobs and pending:
            for _ in range(max_jobs - len(running)):
                # Can submit a new job
                row_start, row_stop, col_start, col_stop = pending.pop()
                t = Task(
                    name=f"pronto-cmp-{submitted}",
                    fn=_compare,
                    args=(kvdb_path, row_start, row_stop, col_start, col_stop,
                          outdir, processes, tmpdir),
                    scheduler=dict(queue=job_queue, cpu=processes, mem=1000,
                                   scratch=5000)
                )
                t.run()
                running.append(t)
                submitted += 1
        elif pending:
            # Submit all jobs at once
            for row_start, row_stop, col_start, col_stop in pending:
                t = Task(
                    name=f"pronto-cmp-{submitted}",
                    fn=_compare,
                    args=(kvdb_path, row_start, row_stop, col_start, col_stop,
                          outdir, processes, tmpdir),
                    scheduler=dict(queue=job_queue, cpu=processes, mem=1000,
                                   scratch=5000)
                )
                t.run()
                running.append(t)
                submitted += 1
            pending = []

        # Check completed tasks
        _running = []
        for task in running:
            if task.done():
                logger.debug(f"{task.stdout}\n{task.stderr}")
                if task.successful():
                    logger.info(f"{task.name} done")
                else:
                    logger.info(f"{task.name} failed")
                    pending = []  # cancel pending tasks
                    failed = True
            else:
                _running.append(task)

        running = _running
        time.sleep(10)

    os.remove(kvdb_path)
    files = [os.path.join(outdir, f) for f in os.listdir(outdir)]
    if failed:
        for path in files:
            os.remove(path)
        raise RuntimeError("one or more tasks failed")

    return files


def cmp_descriptions(user: str, dsn: str, outdir: str, processes: int=8,
                     max_jobs: int=0, tmpdir: Optional[str]=None,
                     chunk_size: int=10000, job_queue: Optional[str]=None):
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
    files = _run_comparisons(user, dsn, query, outdir, processes, max_jobs,
                             tmpdir, chunk_size, job_queue)
    _load_comparisons(user, dsn, "DESC", files)


def cmp_taxa(user: str, dsn: str, outdir: str, processes: int=8,
             max_jobs: int=0, tmpdir: Optional[str]=None,
             chunk_size: int=10000, job_queue: Optional[str]=None):
    query = """
        SELECT METHOD_AC, TAX_ID
        FROM {}.METHOD_TAXA
        ORDER BY METHOD_AC
    """.format(user.split('/')[0])
    files = _run_comparisons(user, dsn, query, outdir, processes, max_jobs,
                             tmpdir, chunk_size, job_queue)
    _load_comparisons(user, dsn, "TAXA", files)


def cmp_terms(user: str, dsn: str, outdir: str, processes: int=8,
              max_jobs: int=0, tmpdir: Optional[str]=None,
              chunk_size: int=10000, job_queue: Optional[str]=None):
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
    files = _run_comparisons(user, dsn, query, outdir, processes, max_jobs,
                             tmpdir, chunk_size, job_queue)
    _load_comparisons(user, dsn, "TERM", files)


def compare(user: str, dsn: str, outdir: str, processes: int=8,
            max_jobs: int = 0, tmpdir: Optional[str] = None,
            chunk_size: int = 10000, job_queue: Optional[str]=None):
    cmp_terms(user, dsn, outdir, processes, max_jobs, tmpdir, chunk_size,
              job_queue)
    cmp_descriptions(user, dsn, outdir, processes, max_jobs, tmpdir,
                     chunk_size, job_queue)
    cmp_taxa(user, dsn, outdir, processes, max_jobs, tmpdir, chunk_size,
             job_queue)

    logger.info("indexing/anayzing METHOD_SIMILARITY")
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    orautils.grant(cur, owner, "METHOD_SIMILARITY", "SELECT", "INTERPRO_SELECT")
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
