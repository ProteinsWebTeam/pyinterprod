# -*- coding: utf-8 -*-

import argparse
import configparser
import os

from mundone import Task, Workflow

from pyinterprod import __version__, iprscan
from pyinterprod.interpro import match, protein
from pyinterprod import uniparc


def run_protein_update():
    parser = argparse.ArgumentParser(description="InterPro protein update")
    parser.add_argument("config",
                        metavar="config.ini",
                        help="configuration file")
    parser.add_argument("-t", "--tasks",
                        nargs="*",
                        default=None,
                        metavar="TASK",
                        help="tasks to run")
    parser.add_argument("--dry-run",
                        action="store_true",
                        default=False,
                        help="list tasks to run and exit")
    parser.add_argument("--detach",
                        action="store_true",
                        help="enqueue tasks to run and exit")
    parser.add_argument("-v", "--version", action="version",
                        version=f"%(prog)s {__version__}",
                        help="show the version and exit")
    args = parser.parse_args()

    if not os.path.isfile(args.config):
        parser.error(f"cannot open '{args.config}': no such file or directory")

    config = configparser.ConfigParser()
    config.read(args.config)

    dsn = config["oracle"]["dsn"]
    interpro_url = f"interpro/{config['oracle']['interpro']}@{dsn}"
    iprscan_url = f"iprscan/{config['oracle']['iprscan']}@{dsn}"
    uniparc_url = f"uniparc/{config['oracle']['uniparc']}@{dsn}"

    uniprot_version = config["uniprot"]["version"]

    data_dir = config["misc"]["data_dir"]
    lsf_queue = config["misc"]["lsf_queue"]
    workflow_dir = config["misc"]["workflow_dir"]

    os.makedirs(data_dir, exist_ok=True)
    tasks = [
        Task(
            fn=uniparc.update,
            args=(uniparc_url,),
            name="update-uniparc",
            scheduler=dict(queue=lsf_queue)
        ),
        Task(
            fn=protein.track_changes,
            args=(interpro_url,
                  config["uniprot"]["swiss-prot"],
                  config["uniprot"]["trembl"]),
            kwargs=dict(dir="/scratch/"),
            name="update-proteins",
            scheduler=dict(cpu=2, queue=lsf_queue),
        ),
        Task(
            fn=protein.delete_obsoletes,
            args=(interpro_url,),
            kwargs=dict(truncate=True),
            name="delete-proteins",
            scheduler=dict(queue=lsf_queue),
            requires=["update-proteins"]
        ),
        Task(
            fn=protein.check_proteins_to_scan,
            args=(interpro_url,),
            name="check-proteins",
            scheduler=dict(queue=lsf_queue),
            requires=["update-proteins", "update-uniparc"]
        ),
        Task(
            fn=iprscan.import_matches,
            args=(iprscan_url,),
            kwargs=dict(threads=8),
            name="import-matches",
            scheduler=dict(queue=lsf_queue),
            requires=["update-uniparc"]
        ),
        Task(
            fn=match.update_matches,
            args=(interpro_url, data_dir),
            name="update-matches",
            scheduler=dict(queue=lsf_queue),
            requires=["check-proteins", "import-matches"]
        ),
    ]

    database = os.path.join(workflow_dir, f"{uniprot_version}.sqlite")
    with Workflow(tasks, dir=workflow_dir, database=database) as wf:
        wf.run(args.tasks, dry_run=args.dry_run, monitor=not args.detach)
