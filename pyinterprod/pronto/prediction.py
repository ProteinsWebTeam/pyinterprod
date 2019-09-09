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


SIMILARITY_THRESHOLD = 0.75


def chunk_jobs(cur: cx_Oracle.Cursor, schema: str, chunk_size: int):
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


def _compare_chunk(database: str, start: str, stop: str, task_queue: Queue,
                   done_queue: Queue):
    with Kvdb(database) as kvdb:
        for acc_1, values_1 in iter(task_queue.get, None):
            intersections = {}
            for acc_2, values_2 in kvdb.range(start, stop):
                if acc_1 < acc_2:
                    intersection = len(values_1 & values_2)
                    if intersection:
                        intersections[acc_2] = intersection

            if intersections:
                done_queue.put((acc_1, intersections))


def _persist_comparisons(database: str, start1: str, stop1: str, start2: str,
                         stop2: str, dst: str, queue: Queue):
    with Kvdb(database) as kvdb:
        counts = {}
        for acc, values in kvdb.range(start1, stop1):
            counts[acc] = len(values)

        if start2 != start1:
            for acc, values in kvdb.range(start2, stop2):
                counts[acc] = len(values)

    with Buffer(filepath=dst) as buffer:
        for acc_1, intersections in iter(queue.get, None):
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


def compare_chunk(database: str, start1: str, stop1: str, start2: str,
                  stop2: str, dst: str, processes: int=8,
                  tmpdir: Optional[str]=None):
    fd, _database = mkstemp(dir=tmpdir)
    os.close(fd)
    os.remove(_database)
    shutil.copy(database, _database)
    database = _database

    task_queue = Queue(maxsize=1)
    done_queue = Queue()

    workers = []
    for _ in range(max(1, processes-2)):
        p = Process(target=_compare_chunk,
                    args=(database, start2, stop2, task_queue, done_queue))
        p.start()
        workers.append(p)

    persister = Process(target=_persist_comparisons,
                        args=(database, start1, stop1, start2, stop2, dst,
                              done_queue))
    persister.start()

    with Kvdb(database) as kvdb:
        for acc, values in kvdb.range(start1, stop1):
            task_queue.put((acc, values))

    for _ in workers:
        task_queue.put(None)

    for p in workers:
        p.join()

    done_queue.put(None)
    persister.join()
    os.remove(database)


def compare_signatures(user: str, dsn: str, query: str, outdir: str,
                       queue: Queue, chunk_size: int, tmpdir: Optional[str],
                       job_queue: Optional[str]):
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()

    # Export signatures to local database
    fd, database = mkstemp(dir=tmpdir)
    os.close(fd)
    os.remove(database)
    with Kvdb(database) as kvdb:
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
    for start1, stop1, start2, stop2 in chunk_jobs(cur, owner, chunk_size):
        fd, filepath = mkstemp(dir=outdir)
        os.close(fd)
        t = Task(
            fn=compare_chunk,
            args=(database, start1, stop1, start2, stop2, filepath),
            kwargs=dict(processes=8, tmpdir=tmpdir),
            scheduler=dict(queue=job_queue, cpu=8, mem=4000)
        )
        t.run(workdir=outdir)
        tasks[t] = (start1, stop1, start2, stop2, filepath)

    cur.close()
    con.close()
    while tasks:
        _tasks = {}
        for t in tasks:
            if t.done():
                start1, stop1, start2, stop2, filepath = tasks[t]

                logger.debug(f"{t.stdout}\n{t.stderr}")
                if t.successful():
                    queue.put((query, start1, start2, filepath))
                else:
                    queue.put((None, start1, start2, filepath))
            else:
                _tasks[t] = tasks[t]

        tasks = _tasks
        time.sleep(10)

    queue.put(None)
    os.remove(database)


def compare(user: str, dsn: str, outdir: str, chunk_size: int=10000,
            job_tmpdir: Optional[str]=None, job_queue: Optional[str]=None):
    os.makedirs(outdir, exist_ok=True)
    if job_tmpdir:
        os.makedirs(job_tmpdir, exist_ok=True)

    cmp_queue = Queue()

    owner = user.split('/')[0]
    producers = []
    sources = {}
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
    sources[query] = "desc"
    producers.append(Process(target=compare_signatures,
                             args=(user, dsn, query, outdir, cmp_queue,
                                   chunk_size, job_tmpdir, job_queue)))

    query = f"""
        SELECT METHOD_AC, TAX_ID
        FROM {owner}.METHOD_TAXA
        ORDER BY METHOD_AC
    """
    sources[query] = "taxa"
    producers.append(Process(target=compare_signatures,
                             args=(user, dsn, query, outdir, cmp_queue,
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
    sources[query] = "terms"
    producers.append(Process(target=compare_signatures,
                             args=(user, dsn, query, outdir, cmp_queue,
                                   chunk_size, job_tmpdir, job_queue)))

    for p in producers:
        p.start()

    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    orautils.drop_table(cur, owner, "METHOD_SIMILARITY", purge=True)
    cur.execute(
        f"""
        CREATE TABLE {owner}.METHOD_SIMILARITY
        (
            METHOD_AC1 VARCHAR2(25) NOT NULL,
            METHOD_AC2 VARCHAR2(25) NOT NULL,
            DBCODE1 CHAR(1) NOT NULL,
            DBCODE2 CHAR(1) NOT NULL,
            PROT_COUNT1 NUMBER(*) NOT NULL,
            PROT_COUNT2 NUMBER(*) NOT NULL,
            COLL_COUNT NUMBER(*) NOT NULL,
            PROT_OVER_COUNT NUMBER(*) NOT NULL,
            PROT_SIM NUMBER(*) NOT NULL,
            PROT_PRED CHAR(1),
            RESI_PRED CHAR(1),
            DESC_PRED CHAR(1),
            TAXA_PRED CHAR(1),
            TERM_PRED CHAR(1)
        ) NOLOGGING
        """
    )
    cur.close()
    con.close()

    chunk_queue = Queue()
    loaders = []
    for _ in range(4):
        p = Process(target=load_similarities,
                    args=(user, dsn, chunk_queue))
        p.start()
        loaders.append(p)

    running = len(producers)
    chunks = {}
    failed = 0
    while running:
        for query, start1, start2, filepath in iter(cmp_queue.get, None):
            if query is None:
                failed += 1
                continue

            key = (start1, start2)
            if key not in chunks:
                chunks[key] = {k: None for k in sources.values()}

            source = sources[query]
            chunks[key][source] = filepath

            if all(chunks[key].values()):
                chunk_queue.put(chunks.pop(key))

        running -= 1

    for p in producers:
        p.join()

    for _ in loaders:
        chunk_queue.put(None)

    for p in loaders:
        p.join()

    if failed:
        raise RuntimeError(f"{failed} jobs failed")

    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    orautils.grant(cur, owner, "METHOD_SIMILARITY", "SELECT",
                   "INTERPRO_SELECT")
    cur.execute(
        f"""
        CREATE INDEX I_METHOD_SIMILARITY$AC1
        ON {owner}.METHOD_SIMILARITY (METHOD_AC1) NOLOGGING
        """
    )
    cur.execute(
        f"""
        CREATE INDEX I_METHOD_SIMILARITY$AC2
        ON {owner}.METHOD_SIMILARITY (METHOD_AC2) NOLOGGING
        """
    )
    cur.close()
    con.close()


def load_similarities(user: str, dsn: str, queue: Queue):
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    get_query = f"""
        SELECT
          M1.DBCODE, M2.DBCODE, MS.PROT_COUNT1, MS.PROT_COUNT2, 
          MS.COLL_COUNT, MS.PROT_OVER_COUNT,
          MS.RES_COUNT1, MS.RES_COUNT2, MS.RES_OVER_COUNT
        FROM {owner}.METHOD_OVERLAP MS
        INNER JOIN {owner}.METHOD M1 
          ON MS.METHOD_AC1 = M1.METHOD_AC
        INNER JOIN {owner}.METHOD M2 
          ON MS.METHOD_AC2 = M2.METHOD_AC
        WHERE MS.METHOD_AC1 = :1 AND MS.METHOD_AC2 = :2
        """

    table = orautils.TablePopulator(
        con=con,
        query=f"""
            INSERT /*+ APPEND */ INTO {owner}.METHOD_SIMILARITY
            VALUES (
                :1, :2, :3, :4, :5, :6, :7, :8, :9, :10, :11, :12, :13, :14
            )
        """,
        autocommit=True
    )

    indices = {
        "desc": 9,
        "taxa": 10,
        "terms": 11
    }

    for chunk in iter(queue.get, None):
        rows = {}

        for source, filepath in chunk.items():
            with Buffer(filepath) as buffer:
                for acc_1, similarities in buffer:
                    for acc_2, (sim, ct1, ct2) in similarities.items():
                        if acc_1 in rows and acc_2 in rows[acc_1]:
                            row = rows[acc_1][acc_2]
                        else:
                            if acc_1 not in rows:
                                rows[acc_1] = {}

                            cur.execute(get_query, (acc_1, acc_2))
                            row = cur.fetchone()

                            if row is None:
                                continue

                            dbcode_1 = row[0]
                            dbcode_2 = row[1]
                            prot_cnt_1 = row[2]
                            prot_cnt_2 = row[3]
                            coll_cnt = row[4]
                            prot_over_cnt = row[5]
                            res_cnt_1 = row[6]
                            res_cnt_2 = row[7]
                            res_over_cnt = row[8]

                            # Protein overlap similarity
                            psim = prot_over_cnt / (prot_cnt_1 + prot_cnt_2
                                                    - prot_over_cnt)
                            pct1 = prot_over_cnt / prot_cnt_1
                            pct2 = prot_over_cnt / prot_cnt_2

                            # Residue overlap similarity
                            rsim = res_over_cnt / (res_cnt_1 + res_cnt_2
                                                   - res_over_cnt)
                            rct1 = res_over_cnt / res_cnt_1
                            rct2 = res_over_cnt / res_cnt_2

                            # Protein/residue predictions
                            ppred = _predict(psim, pct1, pct2)
                            rpred = _predict(rsim, rct1, rct2)

                            row = rows[acc_1][acc_2] = [
                                dbcode_1, dbcode_2, prot_cnt_1, prot_cnt_2,
                                coll_cnt, prot_over_cnt, psim,
                                ppred, rpred, None, None, None
                            ]

                        i = indices[source]
                        row[i] = _predict(sim, ct1, ct2)

                buffer.remove()

        for acc_1 in rows:
            for acc_2, row in rows[acc_1].items():
                table.insert((acc_1, acc_2) + tuple(row))

    cur.close()
    table.close()
    con.close()


def _predict(similarity: float, containment1: float, containment2: float) -> Optional[str]:
    if similarity >= SIMILARITY_THRESHOLD:
        return 'S'  # Similar
    elif containment1 >= SIMILARITY_THRESHOLD:
        if containment2 >= SIMILARITY_THRESHOLD:
            return 'R'  # Related
        else:
            return 'C'  # Child
    elif containment2 >= SIMILARITY_THRESHOLD:
        return 'P'  # Parent
    else:
        return None  # Dissimilar
