# -*- coding: utf-8 -*-

import os
from argparse import ArgumentParser
from configparser import ConfigParser
from typing import List, Sequence

from mundone import Task, Workflow

from pyinterprod import __version__, iprscan
from pyinterprod import interpro, pronto, uniparc, uniprot


def get_pronto_tasks(interpro_url: str, pronto_url: str, data_dir: str,
                     lsf_queue: str) -> List[Task]:
    return [
        Task(
            fn=pronto.database.import_databases,
            args=(interpro_url, pronto_url),
            name="databases",
            scheduler=dict(mem=100, queue=lsf_queue)
        ),

        # Data from GOAPRO
        Task(
            fn=pronto.goa.import_annotations,
            args=(interpro_url, pronto_url),
            name="annotations",
            scheduler=dict(mem=500, queue=lsf_queue)
        ),

        # Data from SWPREAD
        Task(
            fn=pronto.protein.import_similarity_comments,
            args=(interpro_url, pronto_url),
            name="proteins-similarities",
            scheduler=dict(mem=100, queue=lsf_queue),
        ),
        Task(
            fn=pronto.protein.import_protein_names,
            args=(interpro_url, pronto_url,
                  os.path.join(data_dir, "names.sqlite")),
            kwargs=dict(tmpdir="/scratch/"),
            name="proteins-names",
            scheduler=dict(mem=8000, scratch=30000, queue=lsf_queue),
        ),

        # Data from IPPRO
        Task(
            fn=pronto.protein.import_proteins,
            args=(interpro_url, pronto_url),
            name="proteins",
            scheduler=dict(mem=100, queue=lsf_queue),
        ),
        Task(
            fn=pronto.match.import_matches,
            args=(interpro_url, pronto_url,
                  os.path.join(data_dir, "allseqs.dat")),
            kwargs=dict(tmpdir="/scratch/"),
            name="matches",
            scheduler=dict(mem=10000, scratch=20000, queue=lsf_queue),
            requires=["databases"]
        ),
        Task(
            fn=pronto.match.proc_comp_seq_matches,
            args=(interpro_url, pronto_url,
                  os.path.join(data_dir, "names.sqlite"),
                  os.path.join(data_dir, "compseqs.dat")),
            kwargs=dict(tmpdir="/scratch/", processes=8),
            name="signature2proteins",
            scheduler=dict(cpu=8, mem=16000, scratch=30000, queue=lsf_queue),
            requires=["proteins-names"]
        ),
        Task(
            fn=pronto.signature.import_signatures,
            args=(interpro_url, pronto_url,
                  os.path.join(data_dir, "allseqs.dat"),
                  os.path.join(data_dir, "compseqs.dat")),
            name="signatures",
            scheduler=dict(mem=4000, queue=lsf_queue),
            requires=["matches", "signature2proteins"]
        ),
        Task(
            fn=pronto.taxon.import_taxonomy,
            args=(interpro_url, pronto_url),
            name="taxonomy",
            scheduler=dict(mem=2000, queue=lsf_queue),
        ),
    ]


def prep_email(emails: dict, to: Sequence[str], **kwargs) -> dict:
    emails.update({
        "Server": emails["server"],
        "Sender": emails["sender"],
        "To": set(),
        "Cc": set(),
        "Bcc": set(),
    })

    _to = []
    for key in to:
        try:
            addr = emails[key]
        except KeyError:
            continue
        else:
            if addr:
                emails["To"].add(addr)

    _cc = []
    for key in kwargs.get("cc", []):
        try:
            addr = emails[key]
        except KeyError:
            continue
        else:
            if addr and addr not in emails["To"]:
                emails["Cc"].add(addr)

    _bcc = []
    for key in kwargs.get("bcc", []):
        try:
            addr = emails[key]
        except KeyError:
            continue
        else:
            if addr and addr not in emails["To"] and addr not in emails["Cc"]:
                emails["Bcc"].add(addr)

    return emails


def run_protein_update():
    parser = ArgumentParser(description="InterPro protein update")
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

    config = ConfigParser()
    config.read(args.config)

    dsn = config["oracle"]["dsn"]
    interpro_url = f"interpro/{config['oracle']['interpro']}@{dsn}"
    iprscan_url = f"iprscan/{config['oracle']['iprscan']}@{dsn}"
    uniparc_url = f"uniparc/{config['oracle']['uniparc']}@{dsn}"
    pronto_url = config["postgresql"]["pronto"]

    uniprot_version = config["uniprot"]["version"]
    xrefs_dir = config["uniprot"]["xrefs"]

    emails = dict(config["emails"])

    pronto_app = config["misc"]["pronto_url"]
    data_dir = config["misc"]["data_dir"]
    lsf_queue = config["misc"]["lsf_queue"]
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
            kwargs=dict(threads=8),
            name="import-sites",
            scheduler=dict(queue=lsf_queue),
            requires=["update-uniparc"]
        ),

        # Data from flat files
        Task(
            fn=interpro.protein.track_changes,
            args=(interpro_url,
                  config["uniprot"]["swiss-prot"],
                  config["uniprot"]["trembl"],
                  uniprot_version,
                  config["uniprot"]["date"]),
            kwargs=dict(tmpdir="/scratch/"),
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
            args=(interpro_url, prep_email(emails, to=["interpro"])),
            name="export-sib",
            scheduler=dict(queue=lsf_queue),
            requires=["update-matches"]
        ),
        Task(
            fn=uniprot.unirule.report_integration_changes,
            args=(interpro_url, prep_email(emails,
                                           to=["aa_dev"],
                                           cc=["unirule", "interpro"])),
            name="report-changes",
            scheduler=dict(mem=2000, queue=lsf_queue),
            requires=["update-matches"]
        ),
        Task(
            fn=uniprot.unirule.build_aa_iprscan,
            args=(iprscan_url,),
            name="aa-iprscan",
            scheduler=dict(queue=lsf_queue),
            # Actually depends on impart-matches
            requires=["update-matches"]
        ),
        Task(
            fn=uniprot.unirule.build_xref_condensed,
            args=(interpro_url,),
            name="xref-condensed",
            scheduler=dict(queue=lsf_queue),
            requires=["update-matches"]
        ),
        Task(
            fn=uniprot.unirule.build_xref_summary,
            args=(interpro_url,),
            name="xref-summary",
            scheduler=dict(queue=lsf_queue),
            # `report-changes` uses XREF_SUMMARY so we need to wait
            # until it completes before re-creating the table
            requires=["report-changes"]
        ),
        Task(
            fn=uniprot.exchange.export_xrefs,
            args=(interpro_url, xrefs_dir,
                  prep_email(emails, to=["uniprot_db"],
                             cc=["uniprot_prod", "interpro"])
                  ),
            name="export-xrefs",
            scheduler=dict(queue=lsf_queue),
            requires=["xref-summary"]
        ),
        Task(
            fn=uniprot.unirule.ask_to_snapshot,
            args=(interpro_url, prep_email(emails, to=["interpro"])),
            name="notify-interpro",
            scheduler=dict(queue=lsf_queue),
            requires=["aa-iprscan", "xref-condensed", "xref-summary",
                      "update-fmatches"]
        ),

        # Copy SwissProt descriptions
        Task(
            fn=interpro.signature.export_swissprot_descriptions,
            args=(pronto_url, data_dir),
            name="swissprot-de",
            scheduler=dict(queue=lsf_queue),
        )
    ]

    # Adding Pronto tasks
    for t in get_pronto_tasks(interpro_url, pronto_url, data_dir, lsf_queue):
        # Adding 'pronto-' prefix
        t.name = f"pronto-{t.name}"
        if t.requires:
            t.requires = {f"pronto-{r}" for r in t.requires}
        else:
            # Task without dependency:
            # add some so it's submitted at the end of the protein update
            t.requires = {"swissprot-de", "taxonomy", "update-matches",
                          "update-fmatches"}

        tasks.append(t)

    # Generate and send report to curators
    tasks.append(Task(
        fn=interpro.report.send_report,
        args=(interpro_url, pronto_url, data_dir, pronto_app,
              prep_email(emails, to=["interpro"])),
        name="send-report",
        scheduler=dict(mem=4000, queue=lsf_queue),
        requires=["pronto-annotations", "pronto-proteins-similarities",
                  "pronto-proteins", "pronto-signatures",
                  "pronto-taxonomy"]
    ))

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

    database = os.path.join(workflow_dir, f"{uniprot_version}.sqlite")
    with Workflow(tasks, dir=workflow_dir, database=database) as wf:
        wf.run(args.tasks, dry_run=args.dry_run, monitor=not args.detach)


def run_pronto_update():
    parser = ArgumentParser(description="InterPro Pronto update")
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

    config = ConfigParser()
    config.read(args.config)

    dsn = config["oracle"]["dsn"]
    interpro_url = f"interpro/{config['oracle']['interpro']}@{dsn}"
    pronto_url = config["postgresql"]["pronto"]
    uniprot_version = config["uniprot"]["version"]
    data_dir = config["misc"]["data_dir"]
    lsf_queue = config["misc"]["lsf_queue"]
    workflow_dir = config["misc"]["workflow_dir"]

    os.makedirs(data_dir, exist_ok=True)
    tasks = get_pronto_tasks(interpro_url, pronto_url, data_dir, lsf_queue)

    database = os.path.join(workflow_dir, f"{uniprot_version}.pronto.sqlite")
    with Workflow(tasks, dir=workflow_dir, database=database) as wf:
        wf.run(args.tasks, dry_run=args.dry_run, monitor=not args.detach)
