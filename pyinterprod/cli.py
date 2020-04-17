# -*- coding: utf-8 -*-

import argparse
import configparser
import os

from mundone import Task, Workflow

from pyinterprod import __version__, iprscan
from pyinterprod import interpro, uniparc, uniprot


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
    xrefs_dir = config["uniprot"]["xrefs"]

    data_dir = config["misc"]["data_dir"]
    lsf_queue = config["misc"]["lsf_queue"]
    send_emails = config.getboolean("misc", "send_emails")
    workflow_dir = config["misc"]["workflow_dir"]

    os.makedirs(data_dir, exist_ok=True)
    tasks = [
        # Data from UAREAD
        Task(
            fn=uniparc.update,
            args=(uniparc_url,),
            name="update-uniparc",
            scheduler=dict(queue=lsf_queue)
        ),

        # Data from SWPREAD
        Task(
            fn=interpro.taxonomy.refresh_taxonomy,
            args=(interpro_url,),
            name="taxonomy",
            scheduler=dict(queue=lsf_queue)
        ),

        # Data from ISPRO
        Task(
            fn=iprscan.import_matches,
            args=(iprscan_url,),
            kwargs=dict(threads=8),
            name="import-matches",
            scheduler=dict(queue=lsf_queue),
            requires=["update-uniparc"]
        ),
        Task(
            fn=iprscan.import_sites,
            args=(iprscan_url,),
            name="import-sites",
            scheduler=dict(queue=lsf_queue),
            requires=["update-uniparc"]
        ),

        # Data from flat files
        Task(
            fn=interpro.protein.track_changes,
            args=(interpro_url,
                  config["uniprot"]["swiss-prot"],
                  config["uniprot"]["trembl"]),
            kwargs=dict(dir="/scratch/"),
            name="update-proteins",
            scheduler=dict(cpu=2, queue=lsf_queue, scratch=40000),
        ),

        # Update IPPRO
        Task(
            fn=interpro.protein.delete_obsoletes,
            args=(interpro_url,),
            kwargs=dict(truncate=True),
            name="delete-proteins",
            scheduler=dict(queue=lsf_queue),
            requires=["update-proteins"]
        ),
        Task(
            fn=interpro.protein.check_proteins_to_scan,
            args=(interpro_url,),
            name="check-proteins",
            scheduler=dict(queue=lsf_queue),
            requires=["delete-proteins", "update-uniparc"]
        ),
        Task(
            fn=interpro.match.update_matches,
            args=(interpro_url, data_dir),
            name="update-matches",
            scheduler=dict(queue=lsf_queue),
            requires=["check-proteins", "import-matches"]
        ),
        Task(
            fn=interpro.match.update_feature_matches,
            args=(interpro_url,),
            name="update-fmatches",
            scheduler=dict(queue=lsf_queue),
            requires=["check-proteins", "import-matches"]
        ),

        # Data for UniProt
        Task(
            fn=uniprot.exchange.export_sib,
            args=(interpro_url, send_emails),
            name="export-sib",
            scheduler=dict(queue=lsf_queue),
            requires=["update-matches"]
        ),
        Task(
            fn=uniprot.uniprot.report_integration_changes,
            args=(interpro_url, send_emails),
            name="report-changes",
            scheduler=dict(mem=2000, queue=lsf_queue),
            requires=["update-matches"]
        ),
        Task(
            fn=uniprot.uniprot.build_aa_iprscan,
            args=(iprscan_url,),
            name="aa-iprscan",
            scheduler=dict(queue=lsf_queue),
            # Actually depends on impart-matches
            requires=["update-matches"]
        ),
        Task(
            fn=uniprot.uniprot.build_xref_condensed,
            args=(interpro_url,),
            name="xref-condensed",
            scheduler=dict(queue=lsf_queue),
            requires=["update-matches"]
        ),
        Task(
            fn=uniprot.uniprot.build_xref_summary,
            args=(interpro_url,),
            name="xref-summary",
            scheduler=dict(queue=lsf_queue),
            # `report-changes` uses XREF_SUMMARY so we need to wait
            # until it completes before re-creating the table
            requires=["report-changes"]
        ),
        Task(
            fn=uniprot.exchange.export_xrefs,
            args=(interpro_url, xrefs_dir, send_emails),
            name="export-xrefs",
            scheduler=dict(queue=lsf_queue),
            requires=["xref-summary"]
        ),

        Task(
            fn=uniprot.aa.ask_to_snapshot,
            args=(interpro_url, send_emails),
            name="notify-interpro",
            scheduler=dict(queue=lsf_queue),
            requires=["aa-iprscan", "xref-condensed", "xref-summary",
                      "update-fmatches"]
        ),

        # Not urgent tasks (can be run after everything else)
        Task(
            fn=interpro.match.update_variant_matches,
            args=(interpro_url,),
            name="update-varsplic",
            scheduler=dict(queue=lsf_queue),
            requires=["import-matches"]
        ),
        Task(
            fn=interpro.match.update_site_matches,
            args=(interpro_url,),
            name="update-sites",
            scheduler=dict(queue=lsf_queue),
            requires=["import-sites", "update-matches"]
        ),
        Task(
            fn=interpro.legacy.refresh_mviews,
            args=(interpro_url, send_emails),
            name="resfresh-mv",
            scheduler=dict(queue=lsf_queue),
            requires=["update-matches"]
        ),
    ]

    database = os.path.join(workflow_dir, f"{uniprot_version}.sqlite")
    with Workflow(tasks, dir=workflow_dir, database=database) as wf:
        wf.run(args.tasks, dry_run=args.dry_run, monitor=not args.detach)
