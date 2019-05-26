#!/usr/bin/env python
# -*- coding: utf-8 -*-


def main():
    import argparse
    import json
    import os
    from tempfile import gettempdir

    from mundone import Task, Workflow

    from . import (export, interproscan, matches, proteins, signatures,
                   taxonomy, uniparc)
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
    queue = config["workflow"]["lsf-queue"]

    tasks = [
        Task(
            name="update-proteins",
            fn=proteins.update,
            args=(
                db_users["interpro"], db_dsn,
                config["flat_files"]["swissprot"],
                config["flat_files"]["trembl"],
                config["release"]["version"], config["release"]["date"]
            ),
            kwargs=dict(dir="/scratch"),
            scheduler=dict(queue=queue, mem=500, scratch=32000)
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
            kwargs=dict(max_attempts=3*24, secs=3600, max_workers=4),
            scheduler=dict(queue=queue, mem=500)
        ),
        Task(
            name="import-sites",
            fn=interproscan.import_sites,
            args=(db_users["iprscan"], db_dsn),
            kwargs=dict(max_attempts=3*24, secs=3600),
            scheduler=dict(queue=queue, mem=500)
        ),
        Task(
            name="update-matches",
            fn=matches.update_matches,
            args=(db_users["interpro"], db_dsn, config["export"]["matches"]),
            kwargs=dict(drop_indices=True),
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
            requires=["import-sites", "proteins2scan"]
        ),
        Task(
            name="aa-iprscan",
            fn=export.build_aa_iprscan,
            args=(db_users["iprscan"], db_dsn),
            scheduler=dict(queue=queue, mem=500),
            requires=["import-matches"]
        ),
        Task(
            name="xref-summary",
            fn=export.build_xref_summary,
            args=(db_users["interpro"], db_dsn),
            scheduler=dict(queue=queue, mem=500),
            requires=["update-matches"]
        ),
        Task(
            name="xref-condensed",
            fn=export.build_xref_condensed,
            args=(db_users["interpro"], db_dsn),
            scheduler=dict(queue=queue, mem=500),
            requires=["update-matches"]
        ),
        Task(
            name="dump-xrefs",
            fn=export.export_databases,
            args=(db_users["interpro"], db_dsn, config["export"]["xrefs"]),
            scheduler=dict(queue=queue, mem=500),
            requires=["xref-summary"]
        ),
        Task(
            name="taxonomy",
            fn=taxonomy.update_taxonomy,
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
            args=(db_dsn, db_users["pronto_main"], db_users["pronto_alt"]),
            kwargs=dict(tmpdir="/scratch", processes=16),
            scheduler=dict(queue=queue, cpu=16, mem=32000, scratch=32000),
            requires=["update-matches", "update-feature-matches", "taxonomy",
                      "signatures-descriptions"]
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
    with Workflow(tasks, name=wname, dir=wdir, mail=None) as w:
        w.run(args.tasks, resume=args.resume, dry=args.dry_run)
