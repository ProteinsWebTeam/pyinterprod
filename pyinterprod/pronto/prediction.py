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


def _compare_signatures(user: str, dsn: str, query: str, outdir: str,
                        queue: Queue, chunk_size: int, tmpdir: Optional[str],
                        job_queue: Optional[str]):
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
        tasks[t] = (start1, stop1, start2, stop2, filepath)

    cur.close()
    con.close()
    running = len(tasks)
    while running:
        for t in tasks:
            if t.done():
                start1, stop1, start2, stop2, filepath = tasks[t]
                running -= 1

                logger.debug(f"{task.stdout}\n{task.stderr}")
                if t.successful():
                    queue.put((query, start1, stop1, start2, stop2, filepath))
                else:
                    queue.put((None, start1, stop1, start2, stop2, filepath))

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
    producers.append(Process(target=_compare_signatures,
                             args=(user, dsn, query, outdir, queue,
                                   chunk_size, job_tmpdir, job_queue)))

    query = f"""
        SELECT METHOD_AC, TAX_ID
        FROM {owner}.METHOD_TAXA
        ORDER BY METHOD_AC
    """
    sources[query] = "taxa"
    producers.append(Process(target=_compare_signatures,
                             args=(user, dsn, query, outdir, queue,
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
    producers.append(Process(target=_compare_signatures,
                             args=(user, dsn, query, outdir, queue,
                                   chunk_size, job_tmpdir, job_queue)))

    for p in producers:
        p.start()

    buffer = Buffer(os.path.join(outdir, "results"))

    running = len(producers)
    chunks = {}
    failed = 0
    while running:
        for query, start1, stop1, start2, stop2, filepath in iter(queue.get, None):
            if query is None:
                failed += 1
                continue

            source = sources[query]
            key = (start1, stop1, start2, stop2)
            try:
                c = chunks[key]
            except KeyError:
                c = chunks[key] = {
                    "desc": None,
                    "taxa": None,
                    "term": None
                }

            c[source] = filepath
            if all(c.values()):
                buffer.add((start1, stop1, start2, stop2, c))

        running -= 1

    for p in producers:
        p.join()

    logger.warning(f"{failed} jobs failed")

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


class Predictor(object):
    def __init__(self, user: str, dsn: str):
        owner = user.split('/')[0]
        self.con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
        cur = self.con.cursor()
        orautils.drop_table(cur, owner, "METHOD_SIMILARITY", purge=True)
        cur.execute(
            f"""
            CREATE TABLE {owner}.METHOD_SIMILARITY
            (
                METHOD_AC1 VARCHAR2(25) NOT NULL,
                PROT_COUNT1 NUMBER(*) NOT NULL,
                METHOD_AC2 VARCHAR2(25) NOT NULL,
                PROT_COUNT2 NUMBER(*) NOT NULL,
                COLL_COUNT NUMBER(*) NOT NULL,
                PROT_OVER_COUNT NUMBER(*) NOT NULL,
                PROT_PRED CHAR(1) NOT NULL,
                DESC_PRED CHAR(1) NOT NULL,
                TAXA_PRED CHAR(1) NOT NULL,
                TERM_PRED CHAR(1) NOT NULL
            ) NOLOGGING
            """
        )
        cur.close()

        query = f"""
            INSERT /*+ APPEND */ INTO {owner}.METHOD_SIMILARITY
            VALUES (;1, :2, :3, :4, :5, :6, :7, :8, :9, :10)
        """
        self.table = orautils.TablePopulator(self.con, query, autocommit=True)

        self.query = f"""
            SELECT 
              METHOD_AC1, METHOD_AC2, PROT_COUNT1, RES_COUNT1, PROT_COUNT2, 
              RES_COUNT2, COLL_COUNT, PROT_OVER_COUNT, RES_OVER_COUNT
            FROM {owner}.METHOD_OVERLAP
            WHERE METHOD_AC1 BETWEEN :1 AND :2
              AND METHOD_AC2 BETWEEN :3 AND :4
        """

    def digest(self, start1: str, stop1: str, start2: str, stop2: str, chunk: Dict[str, str]):
        cur = self.con.cursor()
        cur.execute(self.query, (start1, stop1, start2, stop2))

        overlaps = {}
        for row in cur:
            acc_1 = row[0]
            acc_2 = row[1]

            if acc_1 in overlaps:
                overlaps[acc_1][acc_2] = row[2:]
            else:
                overlaps[acc_1] = {acc_2: row[2:]}

        cur.close()

        data = {}
        for source, filepath in chunk.items():
            with Buffer(filepath) as buffer:
                for acc_1, similarities in buffer:
                    pass
