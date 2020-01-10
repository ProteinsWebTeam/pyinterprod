# -*- coding: utf-8 -*-


import argparse
import json
import os
import sys

from mundone import Task, Workflow

from . import (interproscan, matches, misc, proteins, signatures,
               uniparc, uniprot)
from .. import __version__, pronto


def run(config_path: str, **kwargs):
    raise_on_error = kwargs.get("raise_on_error", False)
    args_tasks = kwargs.get("tasks")
    resume = kwargs.get("resume")
    dry_run = kwargs.get("dry_run")
    submit = kwargs.get("submit")

    exclude = kwargs.get("exclude", [])

    with open(config_path, "rt") as fh:
        config = json.load(fh)

    db_info = config["database"]
    db_dsn = db_info["dsn"]
    db_users = db_info["users"]
    paths = config["paths"]
    queue = config["workflow"]["lsf-queue"]
    notify = config["email-notifications"]
    os.makedirs(paths["results"], exist_ok=True)

    tasks = [
        Task(
            name="load-proteins",
            fn=proteins.load,
            args=(
                db_users["interpro"], db_dsn,
                paths["flat-files"]["swissprot"],
                paths["flat-files"]["trembl"]
            ),
            kwargs=dict(tmpdir="/scratch/"),
            scheduler=dict(queue=queue, cpu=2, mem=500, scratch=40000)
        ),
        Task(
            name="update-proteins",
            fn=proteins.update,
            args=(
                db_users["interpro"], db_dsn,
                config["release"]["version"], config["release"]["date"]
            ),
            scheduler=dict(queue=queue, mem=500),
            requires=["load-proteins"]
        ),
        Task(
            name="update-uniparc",
            fn=uniparc.update,
            args=(db_users["uniparc"], db_dsn),
            scheduler=dict(queue=queue, mem=500)
        ),
        Task(
            name="check-crc64",
            fn=proteins.check_crc64,
            args=(db_users["interpro"], db_dsn),
            scheduler=dict(queue=queue, mem=500),
            requires=["update-proteins", "update-uniparc"]
        ),
        Task(
            name="proteins2scan",
            fn=proteins.find_protein_to_refresh,
            args=(db_users["interpro"], db_dsn),
            scheduler=dict(queue=queue, mem=500),
            requires=["check-crc64"]
        ),
        Task(
            name="import-matches",
            fn=interproscan.import_matches,
            args=(db_users["iprscan"], db_dsn),
            kwargs=dict(max_workers=8),
            scheduler=dict(queue=queue, mem=500)
        ),
        Task(
            name="import-sites",
            fn=interproscan.import_sites,
            args=(db_users["iprscan"], db_dsn),
            scheduler=dict(queue=queue, mem=500)
        ),
        Task(
            name="update-variants",
            fn=matches.update_variant_matches,
            args=(db_users["interpro"], db_dsn),
            scheduler=dict(queue=queue, mem=500),
            requires=["import-matches", "update-uniparc"]
        ),
        Task(
            name="update-matches",
            fn=matches.update_matches,
            args=(db_users["interpro"], db_dsn, paths["results"]),
            scheduler=dict(queue=queue, mem=500),
            requires=["import-matches", "proteins2scan"]
        ),
        Task(
            name="update-feature-matches",
            fn=matches.update_feature_matches,
            args=(db_users["interpro"], db_dsn),
            kwargs=dict(drop_indices=True),
            scheduler=dict(queue=queue, mem=500),
            requires=["import-matches", "proteins2scan"]
        ),
        Task(
            name="update-sites",
            fn=matches.update_site_matches,
            args=(db_users["interpro"], db_dsn),
            kwargs=dict(drop_indices=True),
            scheduler=dict(queue=queue, mem=500),
            requires=["import-sites", "update-matches"]
        ),
        Task(
            name="aa-iprscan",
            fn=uniprot.build_aa_iprscan,
            args=(db_users["iprscan"], db_dsn),
            scheduler=dict(queue=queue, mem=500),
            requires=["update-matches"]
        ),
        Task(
            name="report-unintegrated",
            fn=uniprot.report_integration_changes,
            kwargs=dict(notify=notify),
            args=(db_users["interpro"], db_dsn),
            scheduler=dict(queue=queue, mem=2000),
            requires=["update-matches"]
        ),
        Task(
            name="xref-summary",
            fn=uniprot.build_xref_summary,
            args=(db_users["interpro"], db_dsn),
            scheduler=dict(queue=queue, mem=500),
            requires=["report-unintegrated"]
        ),
        Task(
            name="xref-condensed",
            fn=uniprot.build_xref_condensed,
            args=(db_users["interpro"], db_dsn),
            scheduler=dict(queue=queue, mem=500),
            requires=["update-matches"]
        ),
        Task(
            name="alert-interpro",
            fn=uniprot.ask_to_snapshot,
            args=(db_users["interpro"], db_dsn),
            kwargs=dict(notify=notify),
            scheduler=dict(queue=queue, mem=500),
            requires=["aa-iprscan", "xref-summary", "xref-condensed",
                      "update-feature-matches"]
        ),
        Task(
            name="dump-xrefs",
            fn=uniprot.export_databases,
            args=(db_users["interpro"], db_dsn, paths["xrefs"]),
            kwargs=dict(notify=notify),
            scheduler=dict(queue=queue, mem=500),
            requires=["xref-summary"]
        ),
        Task(
            name="match-counts",
            fn=misc.refresh_mviews,
            args=(db_users["interpro"], db_dsn),
            kwargs=dict(notify=notify),
            scheduler=dict(queue=queue, mem=500),
            requires=["update-matches"]
        ),
        Task(
            name="export-sib",
            fn=misc.export_for_sib,
            args=(db_users["interpro"], db_dsn),
            kwargs=dict(notify=notify),
            scheduler=dict(queue=queue, mem=500),
            requires=["update-matches"]
        ),
        Task(
            name="taxonomy",
            fn=misc.refresh_taxonomy,
            args=(db_users["interpro"], db_dsn),
            scheduler=dict(queue=queue, mem=500)
        ),
        Task(
            name="signatures-descriptions",
            fn=signatures.update_method2descriptions,
            args=(db_users["interpro"], db_dsn),
            scheduler=dict(queue=queue, mem=500)
        ),
        Task(
            name="pronto",
            fn=pronto.run,
            args=(config_path,),
            kwargs=dict(
                raise_on_error=True,
                report=os.path.join(paths["results"], "swiss_de_changes.tsv"),
                exclude=["copy"]
            ),
            scheduler=dict(queue=queue, mem=500),
            requires=["update-matches", "update-feature-matches", "taxonomy",
                      "signatures-descriptions"]
        ),
        Task(
            name="pronto-copy",
            fn=pronto.run,
            args=(config_path,),
            kwargs=dict(
                raise_on_error=True,
                tasks=["copy"]
            ),
            scheduler=dict(queue=queue, mem=500),
            requires=["pronto"]
        ),
        Task(
            name="report-curators",
            fn=misc.report_to_curators,
            args=(db_users["interpro"], db_dsn, paths["results"]),
            kwargs=dict(notify=notify),
            scheduler=dict(queue=queue, mem=500),
            requires=["pronto"]
        ),
        Task(
            name="unirule",
            fn=uniprot.import_unirules,
            args=(db_users["interpro"], db_dsn, paths["unirule"]),
            scheduler=dict(queue=queue, mem=500),
            requires=["pronto"]
        ),
    ]

    task_names = [t.name for t in tasks]

    if args_tasks:
        for arg in args_tasks:
            if arg not in task_names:
                print(
                    "argument -t/--tasks: "
                    "invalid choice: '{}' (choose from {})\n".format(
                        arg, ", ".join(args_tasks)
                    )
                )
                sys.exit(1)

    if not args_tasks and exclude:
        task_names = [t.name for t in tasks if t.name not in exclude]

    wdir = os.path.join(config["workflow"]["dir"],
                        config["release"]["version"])
    wdb = os.path.join(wdir, "proteinupdate.db")
    wname = "Protein Update"
    with Workflow(tasks, db=wdb, dir=wdir, name=wname) as w:
        status = w.run(task_names, resume=resume, dry=dry_run,
                       secs=submit)

    if not status and raise_on_error:
        raise RuntimeError("one or several tasks failed")

    return status


def main():

    parser = argparse.ArgumentParser(description="InterPro protein update")
    parser.add_argument("config", metavar="CONFIG.JSON",
                        help="config JSON file")
    parser.add_argument("-t", "--tasks", nargs="*",
                        metavar="TASK", help="tasks to run (default: all)")
    parser.add_argument("--dry-run", action="store_true", default=False,
                        help="list tasks to run and exit (default: off)")
    parser.add_argument("--resume", action="store_true", default=False,
                        help="skip completed tasks (default: off)")
    parser.add_argument("--submit", action="store_const",
                        const=0, default=60,
                        help="submit tasks to run and exit (default: off)")
    parser.add_argument("-v", "--version", action="version",
                        version="%(prog)s {}".format(__version__),
                        help="show the version and exit")
    args = parser.parse_args()

    if not os.path.isfile(args.config):
        parser.error(f"{args.config}: no such file or directory")
    success = run(args.config, tasks=args.tasks,
                  dry_run=args.dry_run, resume=args.resume, submit=args.submit)

    sys.exit(0 if success else 1)
