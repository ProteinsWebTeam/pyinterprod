#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import json
import logging
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from tempfile import gettempdir

import cx_Oracle

from .. import logger, orautils
from . import go, protein, signature


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
    then = {}
    for acc, description in cur:
        if acc in then:
            then[acc].add(description)
        else:
            then[acc] = {description}

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
    now = {}
    for acc, description in cur:
        if acc in now:
            now[acc].add(description)
        else:
            now[acc] = {description}

    cur.close()
    con.close()

    changes = {}
    for acc, descs_then in then.items():
        try:
            descs_now = now.pop(acc)
        except KeyError:
            # Lost all descriptions
            changes[acc] = ([], descs_then)
        else:
            changes[acc] = (descs_then-descs_now, descs_now-descs_then)

    for acc, descs_now in now.items():
        changes[acc] = (descs_now, [])

    dst_tmp = dst + ".tmp"
    with open(dst_tmp, "wt") as fh:
        fh.write("Accession\tChecked\t# Lost\t# Gained\tLost\tGained\n")
        for acc in sorted(changes):
            lost, gained = changes[acc]
            fh.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(
                acc, entries[acc], len(lost), len(gained), " | ".join(lost),
                " | ".join(gained)
            ))

    os.rename(dst_tmp, dst)


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


def _default_steps() -> list:
    return sorted(k for k, v in _get_steps().items() if not v.get("skip"))


def _get_steps() -> dict:
    return {
        "databases": {
            "func": load_databases
        },
        "taxa": {
            "func": load_taxa
        },
        "annotations": {
            "func": go.load_annotations
        },
        "terms": {
            "func": go.load_terms
        },
        "publications": {
            "func": go.load_publications
        },
        "comments": {
            "func": protein.load_comments
        },
        "descriptions": {
            "func": protein.load_descriptions
        },
        "enzymes": {
            "func": protein.load_enzymes
        },
        "proteins": {
            "func": protein.load_proteins
        },

        "signatures": {
            "func": signature.load_signatures
        },
        "matches": {
            "func": signature.load_matches
        },
        "signatures2": {
            "func": signature.update_signatures,
            "requires": ("signatures", "matches")
        },
        "signature2protein": {
            "func": signature.load_signature2protein,
            "requires": ("descriptions", "signatures", "taxa", "terms")
        },
        "copy": {
            "func": copy_schema,
            "requires": ("annotations", "comments", "databases", "enzymes",
                         "proteins", "publications", "signatures2",
                         "signature2protein")
        },
        "report": {
            "func": report_description_changes,
            "requires": ("copy",)
        }
    }


def run(user1: str, user2: str, dsn: str, **kwargs):
    report_dst = kwargs.get("report", "swiss_de_families.tsv")
    level = kwargs.get("level", logging.INFO)
    processes = kwargs.get("processes", 1)
    steps = kwargs.get("steps", _get_steps())
    tmpdir = kwargs.get("tmpdir", gettempdir())

    logger.setLevel(level)

    for name, step in steps.items():
        if name == "signature2protein":
            step["args"] = (user1, dsn, processes, tmpdir)
        elif name == "copy":
            step["args"] = (user1, user2, dsn)
        elif name == "report":
            step["args"] = (user1, dsn, report_dst)
        else:
            step["args"] = (user1, dsn)

    for step in steps.values():
        step["requires"] = list(step.get("requires", []))

    with ThreadPoolExecutor(max_workers=processes) as executor:
        pending = {}
        running = {}
        for name, step in steps.items():
            for req_name in step["requires"]:
                if req_name in steps:
                    pending[name] = step
                    break
            else:
                # no requirement scheduled to run
                running[name] = step

        fs = {}
        for name, step in running.items():
            logger.info("{:<20}running".format(name))
            f = executor.submit(step["func"], *step["args"])
            fs[f] = name

        done = set()
        failed = set()
        while fs or pending:
            for f in as_completed(fs):
                name = fs[f]
                try:
                    f.result()
                except Exception as exc:
                    logger.error("{:<19}failed ({})".format(name, exc))
                    failed.add(name)
                else:
                    logger.info("{:<20}done".format(name))
                    done.add(name)

                del running[name]

                # Look if any pending step can be submitted/cancelled
                num_submitted = 0
                for pend_name in list(pending.keys()):
                    pend_step = pending[pend_name]
                    tmp = []
                    for req_name in pend_step["requires"]:
                        if req_name in failed:
                            # cancel step (one dependency failed)
                            del pending[pend_name]
                            break
                        elif req_name not in done:
                            tmp.append(req_name)
                    else:
                        # No dependency failed
                        # Update list of dependencies still to run
                        pend_step["requires"] = tmp
                        if not tmp:
                            logger.info("{:<20}running".format(pend_name))
                            del pending[pend_name]
                            running[pend_name] = pend_step
                            f = executor.submit(pend_step["func"], *pend_step["args"])
                            fs[f] = pend_name
                            num_submitted += 1

                if num_submitted:
                    break

            fs = {f: n for f, n in fs.items() if n in running}

    if failed:
        raise RuntimeError("one or more step failed")


def main():
    parser = argparse.ArgumentParser(description="Pronto schema update")
    parser.add_argument("config", metavar="CONFIG.JSON",
                        help="config JSON file")
    parser.add_argument("-s", "--steps", nargs="+", choices=_get_steps(),
                        default=_default_steps(),
                        help="steps to run (default: all)")
    parser.add_argument("-t", "--tmp", metavar="DIRECTORY",
                        help="temporary directory", default=gettempdir())
    parser.add_argument("-p", "--processes", type=int, default=1,
                        help="number of processes (default: 1)")
    parser.add_argument("-o", "--output", default="swiss_de_families.tsv",
                        help="output report for curators "
                             "(default: swiss_de_families.tsv)")
    parser.add_argument("--verbose", action="store_const",
                        const=logging.DEBUG, default=logging.INFO,
                        help="display additional logging messages")
    args = parser.parse_args()

    try:
        with open(args.config, "rt") as fh:
            config = json.load(fh)
    except FileNotFoundError:
        parser.error("{}: no such file or directory".format(args.config))
    except json.JSONDecodeError:
        parser.error("{}: not a valid JSON file".format(args.config))

    try:
        os.makedirs(args.tmp, exist_ok=True)
    except Exception as e:
        parser.error(e)

    run(config["database"]["users"]["pronto_main"],
        config["database"]["users"]["pronto_alt"],
        config["database"]["dsn"],
        steps={k: v for k, v in _get_steps().items() if k in args.steps},
        tmpdir=args.tmp,
        processes=args.processes,
        level=args.verbose,
        report=args.output)
