# -*- coding: utf-8 -*-


import argparse
import json
import os
import sys


def main():
    from mundone import Task, Workflow

    from . import (interproscan, matches, misc, proteins, signatures,
                   uniparc, uniprot)
    from .. import __version__, pronto

    parser = argparse.ArgumentParser(description="InterPro protein update")
    parser.add_argument("config", metavar="CONFIG.JSON",
                        help="config JSON file")
    parser.add_argument("-t", "--tasks", nargs="*", default=None,
                        metavar="TASK", help="tasks to run")
    parser.add_argument("--dry-run", action="store_true", default=False,
                        help="list tasks to run and exit")
    parser.add_argument("--resume", action="store_true", default=False,
                        help="skip completed tasks")
    parser.add_argument("--submit", action="store_const",
                        const=0, default=60,
                        help="submit tasks to run and exit")
    parser.add_argument("-v", "--version", action="version",
                        version="%(prog)s {}".format(__version__),
                        help="show the version and exit")
    args = parser.parse_args()

    try:
        with open(args.config, "rt") as fh:
            config = json.load(fh)
    except FileNotFoundError:
        parser.error("{}: no such file or directory".format(args.config))
    except json.JSONDecodeError:
        parser.error("{}: not a valid JSON file".format(args.config))

    db_info = config["database"]
    db_dsn = db_info["dsn"]
    db_users = db_info["users"]
    export_dirs = config["export"]
    queue = config["workflow"]["lsf-queue"]

    tasks = [
        Task(
            name="load-proteins",
            # fn=proteins.load,
            fn=proteins.load2,
            args=(
                db_users["interpro"], db_dsn,
                config["flat_files"]["swissprot"],
                config["flat_files"]["trembl"]
            ),
            kwargs=dict(dir="/scratch/"),
            #scheduler=dict(queue=queue, mem=500, scratch=16000)
            scheduler=dict(queue=queue, cpu=2, mem=8000, scratch=40000)
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
            name="update-matches",
            fn=matches.update_matches,
            args=(db_users["interpro"], db_dsn, export_dirs["matches"]),
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
            scheduler=dict(queue=queue, mem=500),
            requires=["aa-iprscan", "xref-summary", "xref-condensed",
                      "update-feature-matches"]
        ),
        Task(
            name="dump-xrefs",
            fn=uniprot.export_databases,
            args=(db_users["interpro"], db_dsn, export_dirs["xrefs"]),
            scheduler=dict(queue=queue, mem=500),
            requires=["xref-summary"]
        ),
        Task(
            name="match-counts",
            fn=misc.refresh_mviews,
            args=(db_users["interpro"], db_dsn),
            # TODO: set to `True` after all DB links are OK
            kwargs=dict(plsql=False),
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
            args=(db_users["pronto_main"], db_users["pronto_alt"], db_dsn),
            kwargs=dict(
                processes=16,
                report=os.path.join(export_dirs["matches"], "swiss_de_families.tsv"),
                tmpdir="/scratch/"
            ),
            scheduler=dict(queue=queue, cpu=16, mem=32000, scratch=32000),
            requires=["update-matches", "update-feature-matches", "taxonomy",
                      "signatures-descriptions"]
        ),
        Task(
            name="report-curators",
            fn=misc.report_to_curators,
            args=(db_users["interpro"], db_dsn, [
                os.path.join(export_dirs["matches"], "entries_changes.tsv"),
                os.path.join(export_dirs["matches"], "swiss_de_families.tsv")
            ]),
            scheduler=dict(queue=queue, mem=500),
            requires=["pronto"]
        )
    ]

    task_names = [t.name for t in tasks]

    if args.tasks:
        for arg in args.tasks:
            if arg not in task_names:
                parser.error(
                    "argument -t/--tasks: "
                    "invalid choice: '{}' (choose from {})\n".format(
                        arg,
                        ", ".join(map("'{}'".format, task_names))
                    )
                )

    wdir = config["workflow"]["dir"]
    wname = "Protein Update"
    with Workflow(tasks, name=wname, dir=wdir) as w:
        success = w.run(args.tasks, resume=args.resume, dry=args.dry_run,
                        secs=args.submit)

    sys.exit(0 if success else 1)
