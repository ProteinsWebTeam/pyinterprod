# -*- coding: utf-8 -*-

import argparse
import configparser
import os

from mundone import Task, Workflow

from pyinterprod import __version__, iprscan
from pyinterprod.interpro import protein
from pyinterprod.uniprot import flatfile
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
    interpro_url = f"{config['oracle']['interpro']}@{dsn}"
    iprscan_url = f"{config['oracle']['iprscan']}@{dsn}"
    uniparc_url = f"{config['oracle']['uniparc']}@{dsn}"

    uniprot_version = config["uniprot"]["version"]

    data_dir = config["misc"]["data_dir"]
    lsf_queue = config["misc"]["lsf_queue"]
    workflow_dir = config["misc"]["lsf_queue"]

    os.makedirs(data_dir, exist_ok=True)
    tasks = [
        Task(
            fn=uniparc.update,
            args=(uniparc_url,),
            name="update-uniparc",
            scheduler=dict(queue=lsf_queue)
        ),
        Task(
            fn=protein.export_proteins,
            args=(interpro_url,
                  os.path.join(data_dir, "proteins.old.sqlite")),
            name="export-old-proteins",
            scheduler=dict(queue=lsf_queue)
        ),
        Task(
            fn=flatfile.load,
            args=(config["uniprot"]["swiss-prot"],
                  config["uniprot"]["trembl"],
                  os.path.join(data_dir, "proteins.new.sqlite")),
            name="import-new-proteins",
            scheduler=dict(queue=lsf_queue)
        ),
        Task(
            fn=protein.track_changes,
            args=(interpro_url,
                  os.path.join(data_dir, "proteins.old.sqlite"),
                  os.path.join(data_dir, "proteins.new.sqlite")),
            name="update-proteins",
            scheduler=dict(queue=lsf_queue),
            requires=["export-old-proteins", "import-new-proteins"]
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
    ]

    database = os.path.join(workflow_dir, f"{uniprot_version}.sqlite")
    with Workflow(tasks, dir=workflow_dir, database=database) as wf:
        wf.run(args.tasks, dry_run=args.dry_run, monitor=not args.detach)
