import argparse
import json
import os
import sys


def main():
    from mundone import Task, Workflow
    from . import (
        methods,
        generate_old_new_stats,
        report,
        delete_dead_signatures,
        match_tmp,
        data_exchange,
        submit_countsupdate,
    )

    from .. import __version__, proteinupdate
    from ..proteinupdate import proteins

    parser = argparse.ArgumentParser(
        description="InterPro member database update")
    parser.add_argument("config", metavar="CONFIG_MEMBER.JSON", help="config JSON file")
    parser.add_argument("-t", "--tasks", nargs="*", metavar="TASK", help="tasks to run (default: all)")
    parser.add_argument(
        "--dry-run",
        action="store_true",
        default=False,
        help="list tasks to run and exit (default: off)",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        default=False,
        help="skip completed tasks (default: off)",
    )
    parser.add_argument(
        "--submit",
        action="store_const",
        const=0,
        default=60,
        help="submit tasks to run and exit (default: off)",
    )

    args = parser.parse_args()

    if not os.path.isfile(args.config):
        parser.error(f"{args.config}: no such file or directory")

    with open(args.config, "rt") as fh:
        config = json.load(fh)

    db_info = config["database"]
    db_dsn = db_info["dsn"]
    db_users = db_info["users"]
    paths = config["paths"]
    queue = config["workflow"]["lsf-queue"]
    notify = config["email-notifications"]
    #email_sender= config["email-sender"]
    email_receiver=config["email-receiver"]
    iprrelease_date = config["release"]["date"]

    os.makedirs(paths["results"], exist_ok=True)

    wdir = os.path.join(config["workflow"]["dir"],
                        config["release"]["version"])

    memberdb = config["member_databases"]
    memberdb = methods.get_dbcodes_memberdb(
        db_users["interpro"], db_dsn, memberdb)
    print(memberdb)

    tasks = [
        # Task(
        #     name="prepare-method-dat",
        #     fn=methods.create_dat,
        #     args=(
        #         db_users["interpro"], db_dsn,
        #         paths["flat-files"]["method_dat"]
        #     ),
        #     kwargs=dict(tmpdir="/scratch/"),
        #     scheduler=dict(queue=queue, cpu=2, mem=500, scratch=40000)
        # ),
        Task(
            name="populate-method-stg", #populate interpro.method_stg table => ok
            fn=methods.populate_method_stg,
            args=(
                db_users["interpro"],
                db_dsn,
                paths["flat-files"]["method_dat"],
                memberdb,
            ),
            scheduler=dict(queue=queue, mem=500),
        ),
        Task(
            name="generate-old-report", #generate old and new stats reports
            fn=generate_old_new_stats.generate_report,
            args=(db_users["interpro"], db_dsn, wdir, memberdb),
            scheduler=dict(queue=queue, mem=500),
            requires=["populate-method-stg"],
        ),
        Task(
            name="check-crc64",#check crc64 before populate_protein_to_scan => not sure it is needed
            fn=proteins.check_crc64,
            args=(db_users["interpro"], db_dsn),
            scheduler=dict(queue=queue, mem=500)
        ),
        Task(
            name="proteins2scan", #populate_protein_to_scan
            fn=methods.update_proteins2scan,
            args=(db_users["interpro"], db_dsn),
            scheduler=dict(queue=queue, mem=500),
            requires=["check-crc64"],
        ),
        Task(
            name="delete-dead-signatures", #delete obsoleted signatures,execute only if something returned from generate_new_report function
            fn=delete_dead_signatures.delete_dead_signatures,
            args=(db_users["interpro"], db_dsn, memberdb),
            scheduler=dict(queue=queue, mem=500),
        ),
        Task(
            name="update-method", # => ok
            fn=methods.update_method,
            args=(db_users["interpro"], db_dsn),
            scheduler=dict(queue=queue, mem=500),
            requires=["populate-method-stg"],
        ),
        Task(
            name="update-iprscan2dbcode", #update IPRSCAN2DBCODE table
            fn=methods.update_iprscan2dbcode,
            args=(db_users["interpro"], db_dsn, memberdb),
            scheduler=dict(queue=queue, mem=500),
        ),
        Task(
            name="create-match-tmp", #what if new member db? Need to add count of match and match_tmp? => ok
            fn=match_tmp.create_match_tmp,
            args=(db_users["interpro"], db_dsn, memberdb, wdir),
            scheduler=dict(queue=queue, mem=500),
        ),
        Task(
            name="update-match",
            fn=data_exchange.exchange_data,
            args=(db_users["interpro"], db_dsn, memberdb, "MATCH"),
            scheduler=dict(queue=queue, mem=500),
        ),
        Task(
            name="update-db-version", # => ok
            fn=methods.update_db_version,
            args=(db_users["interpro"], db_dsn, memberdb),
            scheduler=dict(queue=queue, mem=500),
            requires=["update-method"],
        ),
        Task(
            name="update-site-match",
            fn=match_tmp.refresh_site_matches,
            args=(db_users["interpro"], db_dsn, memberdb),
            scheduler=dict(queue=queue, mem=500),
        ),
        Task(
            name="drop-match-tmp",
            fn=match_tmp.drop_match_tmp,
            args=(db_users["interpro"], db_dsn),
            scheduler=dict(queue=queue, mem=500),
        ),
        Task(
            name="generate-counts",
            fn=submit_countsupdate.generate_counts,
            args=(db_users["interpro"], db_dsn, memberdb, email_receiver, wdir),
            scheduler=dict(queue=queue, mem=500),
        ),
        Task(
            name="report-changes",
            fn=report.report_swissprot_changes,
            args=(db_users["interpro"], db_dsn, memberdb),
            scheduler=dict(queue=queue, mem=500),
        ),
    ]

    task_names = [t.name for t in tasks]

    if args.tasks:
        for arg in args.tasks:
            if arg not in task_names:
                parser.error(
                    "argument -t/--tasks: "
                    "invalid choice: '{}' (choose from {})\n".format(
                        arg, ", ".join(task_names)
                    )
                )

    # methods.update_db_version("interpro", "IPTST", memberdb)
    wdir = os.path.join(config["workflow"]["dir"],
                        config["release"]["version"])
    wdb = os.path.join(wdir, "memberdbupdate.log")
    wname = "Member database Update"
    with Workflow(tasks, db=wdb, dir=wdir, name=wname) as w:
        success = w.run(
            args.tasks, resume=args.resume, dry=args.dry_run, secs=args.submit
        )

    sys.exit(0 if success else 1)
