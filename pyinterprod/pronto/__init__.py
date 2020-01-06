#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import json
import os
import subprocess
import sys
from typing import List, Optional

import cx_Oracle
from mundone import Task, Workflow

from .. import logger, orautils
from . import go, prediction, protein, signature


_DATA_PUMP_DIR = "PANDA_DATA_PUMP_DIR"
_REPORT = "swiss_de_families.tsv"


def load_databases(user: str, dsn: str):
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    orautils.drop_table(cur, owner, "CV_DATABASE", purge=True)
    cur.execute(
        f"""
        CREATE TABLE {owner}.CV_DATABASE
        (
            DBCODE VARCHAR2(10) NOT NULL,
            DBNAME VARCHAR2(50) NOT NULL,
            DBSHORT VARCHAR2(10) NOT NULL,
            VERSION VARCHAR2(20),
            FILE_DATE DATE,
            IS_READY CHAR(1) DEFAULT 'N',
            CONSTRAINT PK_DATABASE PRIMARY KEY (DBCODE)
        ) NOLOGGING
        """
    )

    cur.execute(
        f"""
        INSERT /*+ APPEND */ INTO {owner}.CV_DATABASE (
            DBCODE, DBNAME, DBSHORT, VERSION, FILE_DATE
        )
        SELECT DB.DBCODE, DB.DBNAME, DB.DBSHORT, V.VERSION, V.FILE_DATE
        FROM INTERPRO.CV_DATABASE DB
        LEFT OUTER JOIN INTERPRO.DB_VERSION V
          ON DB.DBCODE = V.DBCODE
        """
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
        f"""
        CREATE TABLE {owner}.ETAXI
        NOLOGGING
        AS
        SELECT TAX_ID, PARENT_ID, SCIENTIFIC_NAME, RANK, LEFT_NUMBER,
               RIGHT_NUMBER, FULL_NAME
        FROM INTERPRO.ETAXI
        """
    )
    orautils.gather_stats(cur, owner, "ETAXI")
    orautils.grant(cur, owner, "ETAXI", "SELECT", "INTERPRO_SELECT")
    cur.execute(
        f"""
        ALTER TABLE {owner}.ETAXI
        ADD CONSTRAINT PK_ETAXI PRIMARY KEY (TAX_ID)
        """
    )

    orautils.drop_table(cur, owner, "LINEAGE", purge=True)
    cur.execute(
        f"""
        CREATE TABLE {owner}.LINEAGE
        (
            TAX_ID NUMBER(10) NOT NULL,
            RANK VARCHAR2(50) NOT NULL,
            RANK_TAX_ID NUMBER(10)
        ) NOLOGGING
        """
    )

    """
    taxID 131567 (cellular organisms) contains three superkingdoms:
        * Bacteria (2)
        * Archaea (2157)
        * Eukaryota (2759)

    therefore it is not needed (we don't want a meta-superkingdom)
    """
    cur.execute(
        f"""
        SELECT TAX_ID, PARENT_ID, RANK
        FROM {owner}.ETAXI
        WHERE TAX_ID != 131567
        """
    )
    taxa = {}
    for tax_id, parent_id, rank in cur:
        if parent_id == 1:
            taxa[tax_id] = ("superkingdom", parent_id)
        else:
            taxa[tax_id] = (rank, parent_id)

    table = orautils.TablePopulator(con,
                                    query=f"INSERT /*+ APPEND */ "
                                          f"INTO {owner}.LINEAGE "
                                          f"VALUES (:1, :2, :3)",
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
        f"""
        CREATE INDEX I_LINEAGE
        ON {owner}.LINEAGE (TAX_ID, RANK)
        NOLOGGING
        """
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

    cur.execute("SELECT ENTRY_AC, ENTRY_TYPE, CHECKED FROM INTERPRO.ENTRY")
    entries = {row[0]: (row[1], row[2]) for row in cur}

    cur.execute(
        """
        SELECT DISTINCT EM.ENTRY_AC, M.DESCRIPTION
        FROM INTERPRO.ENTRY2METHOD EM
        INNER JOIN INTERPRO.METHOD2SWISS_DE M
          ON EM.METHOD_AC = M.METHOD_AC
        """
    )
    entries_then = {}
    for acc, description in cur:
        if acc in entries_then:
            entries_then[acc].add(description)
        else:
            entries_then[acc] = {description}

    cur.execute(
        f"""
        SELECT DISTINCT EM.ENTRY_AC, D.TEXT
        FROM {owner}.METHOD2PROTEIN PARTITION(M2P_SWISSP) M
        INNER JOIN {owner}.DESC_VALUE D
          ON M.DESC_ID = D.DESC_ID
        INNER JOIN INTERPRO.ENTRY2METHOD EM
          ON M.METHOD_AC = EM.METHOD_AC
       """
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
        fh.write("Accession\tType\tChecked\t# Lost\t# Gained\tLost\tGained\n")
        for acc in sorted(changes):
            gained, lost = changes[acc]
            if lost or gained:
                entry_type, is_checked = entries[acc]
                fh.write(f"{acc}\t{entry_type}\t{is_checked}\t{len(lost)}\t"
                         f"{len(gained)}\t{' | '.join(lost)}\t"
                         f"{' | '.join(gained)}\n")


def copy_tables(user_src: str, user_dst: str, dsn: str, set_status: bool=False):
    if user_src == user_dst:
        logger.warning("identical source and target schemas")
        return

    schema_src = user_src.split('/')[0]
    schema_dst = user_dst.split('/')[0]
    dumpfile = schema_src.upper()

    con = cx_Oracle.connect(orautils.make_connect_string(user_src, dsn))
    cur = con.cursor()
    tables = [f"{schema_src}.{t}"
              for t in orautils.get_tables(cur, schema_src)]

    if set_status:
        cur.execute(f"UPDATE {schema_src}.CV_DATABASE SET IS_READY = 'N'")
        con.commit()

    cur.close()
    con.close()

    returncode = subprocess.call(["expdp",
                                  orautils.make_connect_string(user_src, dsn),
                                  f"TABLES={','.join(tables)}",
                                  f"DIRECTORY={_DATA_PUMP_DIR}",
                                  f"DUMPFILE={dumpfile}.dmp",
                                  f"LOGFILE={dumpfile}-exp.log",
                                  "REUSE_DUMPFILES=YES"])

    if returncode:
        raise RuntimeError(f"expdp exited with code {returncode}")

    con = cx_Oracle.connect(orautils.make_connect_string(user_src, dsn))
    cur = con.cursor()
    cur.execute(f"UPDATE {schema_src}.CV_DATABASE SET IS_READY = 'Y'")
    con.commit()
    cur.close()
    con.close()

    orautils.drop_all(user_dst, dsn)
    returncode = subprocess.call(["impdp",
                                  orautils.make_connect_string(user_dst, dsn),
                                  f"TABLES={','.join(tables)}",
                                  f"DIRECTORY={_DATA_PUMP_DIR}",
                                  f"DUMPFILE={dumpfile}.dmp",
                                  f"LOGFILE={dumpfile}-imp.log",
                                  f"REMAP_SCHEMA={schema_src}:{schema_dst}"])

    if returncode:
        raise RuntimeError(f"impdp exited with code {returncode}")

    con = cx_Oracle.connect(orautils.make_connect_string(user_dst, dsn))
    cur = con.cursor()
    cur.execute(f"UPDATE {schema_dst}.CV_DATABASE SET IS_READY = 'Y'")
    con.commit()
    cur.close()
    con.close()


def run(config_path: str, **kwargs):
    raise_on_error = kwargs.get("raise_on_error", False)
    report = kwargs.get("report", _REPORT)
    task_names = kwargs.get("tasks")
    exclude = kwargs.get("exclude", [])

    with open(config_path, "rt") as fh:
        config = json.load(fh)

    user1 = config["database"]["users"]["pronto-pri"]
    user2 = config["database"]["users"]["pronto-sec"]
    dsn = config["database"]["dsn"]
    outdir = config["paths"]["results"]
    queue = config["workflow"]["lsf-queue"]

    tasks = [
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
            scheduler=dict(queue=queue, cpu=16, mem=32000, scratch=32000),
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
            kwargs=dict(job_tmpdir="/scratch/", job_queue=queue),
            scheduler=dict(queue=queue, cpu=8, mem=16000, scratch=10000),
            requires=["signatures-proteins"]
        ),
        Task(
            name="report",
            fn=report_description_changes,
            args=(user1, dsn, report),
            scheduler=dict(queue=queue, mem=2000),
            requires=["index"]
        ),
        Task(
            name="copy",
            fn=copy_tables,
            args=(user1, user2, dsn),
            scheduler=dict(queue=queue, mem=500),
            requires=["annotations", "publications", "databases", "comments",
                      "enzymes", "proteins", "matches", "index", "compare"]
        )
    ]

    if not task_names and exclude:
        task_names = [t.name for t in tasks if t.name not in exclude]

    wdir = os.path.join(config["workflow"]["dir"], config["release"]["version"])
    wdb = os.path.join(wdir, "pronto.db")
    with Workflow(tasks, db=wdb, dir=wdir) as w:
        status = w.run(task_names, dependencies=False)

    if not status and raise_on_error:
        raise RuntimeError("one or several tasks failed")

    return status


def main():
    parser = argparse.ArgumentParser(description="Pronto schema update")
    parser.add_argument("config", metavar="CONFIG.JSON",
                        help="config JSON file")
    parser.add_argument("-t", "--tasks", nargs="*",
                        help="tasks to run (default: all)")
    parser.add_argument("-o", "--output", default=_REPORT,
                        help=f"output report for curators (default: {_REPORT})")
    args = parser.parse_args()

    success = run(args.config, report=args.output, tasks=args.tasks)
    sys.exit(0 if success else 1)
