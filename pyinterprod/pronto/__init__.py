#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import json
import os
import sys
from typing import List, Optional

import cx_Oracle
from mundone import Task, Workflow

from .. import orautils
from . import go, prediction, protein, signature


def load_databases(user: str, dsn: str):
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    orautils.drop_table(cur, owner, "CV_DATABASE", purge=True)
    cur.execute(
        """
        CREATE TABLE {}.CV_DATABASE
        (
            DBCODE VARCHAR2(10) NOT NULL,
            DBNAME VARCHAR2(50) NOT NULL,
            DBSHORT VARCHAR2(10) NOT NULL,
            VERSION VARCHAR2(20),
            FILE_DATE DATE,
            IS_READY CHAR(1) DEFAULT 'N',
            CONSTRAINT PK_DATABASE PRIMARY KEY (DBCODE)
        ) NOLOGGING
        """.format(owner)
    )

    cur.execute(
        """
        INSERT /*+ APPEND */ INTO {}.CV_DATABASE (
            DBCODE, DBNAME, DBSHORT, VERSION, FILE_DATE
        )
        SELECT DB.DBCODE, DB.DBNAME, DB.DBSHORT, V.VERSION, V.FILE_DATE
        FROM INTERPRO.CV_DATABASE DB
        LEFT OUTER JOIN INTERPRO.DB_VERSION V
          ON DB.DBCODE = V.DBCODE
        """.format(owner)
    )
    con.commit()

    orautils.gather_stats(cur, owner, "CV_DATABASE")
    orautils.grant(cur, owner, "CV_DATABASE", "SELECT", "INTERPRO_SELECT")
    cur.close()
    con.close()


def load_taxa(user: str, dsn: str, **kwargs):
    ranks = kwargs.get("ranks", ("superkingdom", "kingdom", "phylum", "class",
                                 "order", "family", "genus", "species"))
    ranks = set(ranks)

    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    orautils.drop_table(cur, owner, "ETAXI", purge=True)
    cur.execute(
        """
        CREATE TABLE {}.ETAXI
        NOLOGGING
        AS
        SELECT
          TAX_ID, PARENT_ID, SCIENTIFIC_NAME, RANK, FULL_NAME
        FROM INTERPRO.ETAXI
        """.format(owner)
    )
    orautils.gather_stats(cur, owner, "ETAXI")
    orautils.grant(cur, owner, "ETAXI", "SELECT", "INTERPRO_SELECT")
    cur.execute(
        """
        ALTER TABLE {}.ETAXI
        ADD CONSTRAINT PK_ETAXI PRIMARY KEY (TAX_ID)
        """.format(owner)
    )

    orautils.drop_table(cur, owner, "LINEAGE", purge=True)
    cur.execute(
        """
        CREATE TABLE {}.LINEAGE
        (
            TAX_ID NUMBER(10) NOT NULL,
            RANK VARCHAR2(50) NOT NULL,
            RANK_TAX_ID NUMBER(10)
        ) NOLOGGING
        """.format(owner)
    )

    """
    taxID 131567 (cellular organisms) contains three superkingdoms:
        * Bacteria (2)
        * Archaea (2157)
        * Eukaryota (2759)

    therefore it is not needed (we don't want a meta-superkingdom)
    """
    cur.execute(
        """
        SELECT TAX_ID, PARENT_ID, RANK
        FROM {}.ETAXI
        WHERE TAX_ID != 131567
        """.format(owner)
    )
    taxa = {}
    for tax_id, parent_id, rank in cur:
        if parent_id == 1:
            taxa[tax_id] = ("superkingdom", parent_id)
        else:
            taxa[tax_id] = (rank, parent_id)

    table = orautils.TablePopulator(con,
                                    query="INSERT /*+ APPEND */ "
                                          "INTO {}.LINEAGE "
                                          "VALUES (:1, :2, :3)".format(owner),
                                    autocommit=True)
    for tax_id in taxa:
        rank, parent_id = taxa[tax_id]
        if rank in ranks:
            table.insert((tax_id, rank, tax_id))

        while parent_id in taxa:
            rank_tax_id = parent_id
            rank, parent_id = taxa[rank_tax_id]
            if rank in ranks:
                table.insert((tax_id, rank, rank_tax_id))
    table.close()

    orautils.gather_stats(cur, owner, "LINEAGE")
    orautils.grant(cur, owner, "LINEAGE", "SELECT", "INTERPRO_SELECT")
    cur.execute(
        """
        CREATE INDEX I_LINEAGE
        ON {}.LINEAGE (TAX_ID, RANK)
        NOLOGGING
        """.format(owner)
    )

    cur.close()
    con.close()


def report_description_changes(user: str, dsn: str, dst: str):
    try:
        os.remove(dst)
    except FileNotFoundError:
        pass

    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()

    cur.execute("SELECT ENTRY_AC, CHECKED FROM INTERPRO.ENTRY")
    entries = dict(cur.fetchall())

    cur.execute(
        """
        SELECT DISTINCT EM.ENTRY_AC, M.DESCRIPTION
        FROM INTERPRO.ENTRY2METHOD EM
        INNER JOIN INTERPRO.METHOD2SWISS_DE M
          ON EM.METHOD_AC = M.METHOD_AC
        WHERE EM.ENTRY_AC IN (
          SELECT ENTRY_AC FROM INTERPRO.ENTRY WHERE ENTRY_TYPE='F'
        )
        """
    )
    entries_then = {}
    for acc, description in cur:
        if acc in entries_then:
            entries_then[acc].add(description)
        else:
            entries_then[acc] = {description}

    cur.execute(
        """
        SELECT DISTINCT EM.ENTRY_AC, D.TEXT
        FROM {0}.METHOD2PROTEIN PARTITION(M2P_SWISSP) M
        INNER JOIN {0}.DESC_VALUE D
          ON M.DESC_ID = D.DESC_ID
        INNER JOIN INTERPRO.ENTRY2METHOD EM
          ON M.METHOD_AC = EM.METHOD_AC
       WHERE EM.ENTRY_AC IN (
         SELECT ENTRY_AC FROM INTERPRO.ENTRY WHERE ENTRY_TYPE='F'
       )
       """.format(owner)
    )
    entries_now = {}
    for acc, description in cur:
        if acc in entries_now:
            entries_now[acc].add(description)
        else:
            entries_now[acc] = {description}

    cur.close()
    con.close()

    changes = {}  # key: accession, value: (gained, lost)
    for acc, descs_now in entries_now.items():
        try:
            descs_then = entries_then.pop(acc)
        except KeyError:
            # This entry did not have descriptions then
            changes[acc] = (descs_now, [])
        else:
            changes[acc] = (descs_now - descs_then, descs_then - descs_now)

    # Entries that no longer have descriptions
    for acc, descs_then in entries_then.items():
        changes[acc] = ([], descs_then)

    with open(dst, "wt") as fh:
        fh.write("Accession\tChecked\t# Lost\t# Gained\tLost\tGained\n")
        for acc in sorted(changes):
            lost, gained = changes[acc]
            fh.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(
                acc, entries[acc], len(lost), len(gained), " | ".join(lost),
                " | ".join(gained)
            ))


def copy_schema(user_src: str, user_dst: str, dsn: str):
    owner = user_src.split('/')[0]

    con = cx_Oracle.connect(orautils.make_connect_string(user_src, dsn))
    cur = con.cursor()
    cur.execute("UPDATE {}.CV_DATABASE SET IS_READY = 'Y'".format(owner))

    tables = []
    for t in orautils.get_tables(cur, owner):
        tables.append({
            "name": t,
            "grants": orautils.get_grants(cur, owner, t),
            "constraints": orautils.get_constraints(cur, owner, t),
            "indexes": orautils.get_indices(cur, owner, t),
            "partitions": orautils.get_partitions(cur, owner, t)
        })
    #orautils.clear_schema(user_dst, dsn)


def _get_tasks(**kwargs):
    user1 = kwargs.get("user1")
    user2 = kwargs.get("user2")
    dsn = kwargs.get("dsn")
    outdir = kwargs.get("outdir", os.getcwd())
    queue = kwargs.get("queue")
    report_dst = kwargs.get("report", "swiss_de_families.tsv")

    return [
        Task(
            name="annotations",
            fn=go.load_annotations,
            args=(user1, dsn),
            scheduler=dict(queue=queue, mem=500)
        ),
        Task(
            name="publications",
            fn=go.load_publications,
            args=(user1, dsn),
            scheduler=dict(queue=queue, mem=500)
        ),
        Task(
            name="terms",
            fn=go.load_terms,
            args=(user1, dsn),
            scheduler=dict(queue=queue, mem=500)
        ),
        Task(
            name="databases",
            fn=load_databases,
            args=(user1, dsn),
            scheduler=dict(queue=queue, mem=500)
        ),
        Task(
            name="taxa",
            fn=load_taxa,
            args=(user1, dsn),
            scheduler=dict(queue=queue, mem=1000)
        ),
        Task(
            name="comments",
            fn=protein.load_comments,
            args=(user1, dsn),
            scheduler=dict(queue=queue, mem=2000)
        ),
        Task(
            name="descriptions",
            fn=protein.load_descriptions,
            args=(user1, dsn),
            scheduler=dict(queue=queue, mem=1000)
        ),
        Task(
            name="enzymes",
            fn=protein.load_enzymes,
            args=(user1, dsn),
            scheduler=dict(queue=queue, mem=500)
        ),
        Task(
            name="proteins",
            fn=protein.load_proteins,
            args=(user1, dsn),
            scheduler=dict(queue=queue, mem=500)
        ),

        Task(
            name="signatures",
            fn=signature.load_signatures,
            args=(user1, dsn),
            scheduler=dict(queue=queue, mem=500)
        ),
        Task(
            name="matches",
            fn=signature.load_matches,
            args=(user1, dsn),
            scheduler=dict(queue=queue, mem=500),
            requires=["signatures"]
        ),
        Task(
            name="signatures-proteins",
            fn=signature.load_signature2protein,
            args=(user1, dsn),
            kwargs=dict(processes=16, tmpdir="/scratch/"),
            scheduler=dict(queue=queue, cpu=16, mem=24000, scratch=32000),
            requires=["descriptions", "signatures", "taxa", "terms"]
        ),
        Task(
            name="index",
            fn=signature.finalize_method2protein,
            args=(user1, dsn),
            scheduler=dict(queue=queue, mem=500),
            requires=["signatures-proteins"]
        ),
        Task(
            name="compare",
            fn=prediction.compare,
            args=(user1, dsn, outdir),
            kwargs=dict(chunk_size=5000, max_jobs=0, job_processes=8,
                        job_tmpdir="/scratch/", job_queue=queue),
            scheduler=dict(queue=queue, cpu=2, mem=2000, scratch=5000),
            requires=["signatures-proteins"]
        )

        # todo: add report and copy tasks
    ]


def run(config_path: str, task_names: Optional[List[str]]=None,
        report: Optional[str]=None, raise_on_error: bool=True):

    with open(config_path, "rt") as fh:
        config = json.load(fh)

    outdir = config["paths"]["results"]

    if not report:
        report = os.path.join(outdir, "swiss_de_families.tsv")

    tasks = _get_tasks(
        user1=config["database"]["users"]["pronto-pri"],
        user2=config["database"]["users"]["pronto-sec"],
        dsn=config["database"]["dsn"],
        outdir=outdir,
        queue=config["workflow"]["lsf-queue"],
        report=report
    )

    wdir = config["workflow"]["dir"]
    wdb = os.path.join(wdir, "pronto.db")
    with Workflow(tasks, db=wdb, dir=wdir) as w:
        status = w.run(task_names, dependencies=False)

    if not status and raise_on_error:
        raise RuntimeError("")

    return status


def main():
    parser = argparse.ArgumentParser(description="Pronto schema update")
    parser.add_argument("config", metavar="CONFIG.JSON",
                        help="config JSON file")
    parser.add_argument("-t", "--tasks", nargs="+",
                        choices=[t.name for t in _get_tasks()], default=None,
                        help="tasks to run")
    parser.add_argument("-o", "--output", default="swiss_de_families.tsv",
                        help=("output report for curators (default: swiss_de_families.tsv)"))
    args = parser.parse_args()

    success = run(args.config, args.tasks, args.output, raise_on_error=False)
    sys.exit(0 if success else 1)
