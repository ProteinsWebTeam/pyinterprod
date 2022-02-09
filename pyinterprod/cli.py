import os
import sys
from argparse import ArgumentParser
from configparser import ConfigParser
from typing import List

from mundone import Task, Workflow

from pyinterprod import __version__, iprscan
from pyinterprod import interpro, pronto, uniparc, uniprot
from pyinterprod.interpro.clan import update_cdd_clans, update_hmm_clans


def get_pronto_tasks(ora_ipr_url: str, ora_swp_url: str, ora_goa_url: str,
                     pg_ipr_url: str, data_dir: str, lsf_queue: str) -> List[Task]:
    return [
        Task(
            fn=pronto.database.import_databases,
            args=(ora_ipr_url, pg_ipr_url),
            name="databases",
            scheduler=dict(mem=100, queue=lsf_queue)
        ),

        # Data from GOAPRO
        Task(
            fn=pronto.goa.import_annotations,
            args=(ora_goa_url, pg_ipr_url),
            name="annotations",
            scheduler=dict(mem=500, queue=lsf_queue)
        ),

        # Data from SWPREAD
        Task(
            fn=pronto.protein.import_similarity_comments,
            args=(ora_swp_url, pg_ipr_url),
            name="proteins-similarities",
            scheduler=dict(mem=100, queue=lsf_queue),
        ),
        Task(
            fn=pronto.protein.import_protein_names,
            args=(ora_swp_url, pg_ipr_url,
                  os.path.join(data_dir, "names.sqlite")),
            kwargs=dict(tmpdir="/scratch/"),
            name="proteins-names",
            scheduler=dict(mem=8000, scratch=30000, queue=lsf_queue),
        ),

        # Data from IPPRO
        Task(
            fn=pronto.protein.import_proteins,
            args=(ora_ipr_url, pg_ipr_url),
            name="proteins",
            scheduler=dict(mem=100, queue=lsf_queue),
        ),
        Task(
            fn=pronto.match.import_matches,
            args=(ora_ipr_url, pg_ipr_url,
                  os.path.join(data_dir, "allseqs.dat")),
            kwargs=dict(tmpdir="/scratch/"),
            name="matches",
            scheduler=dict(mem=10000, scratch=20000, queue=lsf_queue),
            requires=["databases"]
        ),
        Task(
            fn=pronto.match.create_signature2protein,
            args=(ora_ipr_url, pg_ipr_url,
                  os.path.join(data_dir, "names.sqlite"),
                  os.path.join(data_dir, "compseqs.dat")),
            kwargs=dict(tmpdir="/scratch/", processes=8),
            name="signature2proteins",
            scheduler=dict(cpu=8, mem=16000, scratch=30000, queue=lsf_queue),
            requires=["proteins-names"]
        ),
        Task(
            fn=pronto.signature.import_signatures,
            args=(ora_ipr_url, pg_ipr_url,
                  os.path.join(data_dir, "allseqs.dat"),
                  os.path.join(data_dir, "compseqs.dat")),
            name="signatures",
            scheduler=dict(mem=4000, queue=lsf_queue),
            requires=["matches", "signature2proteins"]
        ),
        Task(
            fn=pronto.taxon.import_taxonomy,
            args=(ora_ipr_url, pg_ipr_url),
            name="taxonomy",
            scheduler=dict(mem=2000, queue=lsf_queue),
        ),
    ]


def check_ispro():
    parser = ArgumentParser(description="Check matches/sites status in ISPRO")
    parser.add_argument("config",
                        metavar="config.ini",
                        help="configuration file")
    parser.add_argument("-t", "--type",
                        default="matches",
                        choices=("matches", "sites"),
                        help="type of data to check (default: matches)")
    parser.add_argument("-s", "--status",
                        default="production",
                        choices=("active", "all", "production"),
                        help="status of analyses to check "
                             "(default: production)")
    args = parser.parse_args()

    if not os.path.isfile(args.config):
        parser.error(f"cannot open '{args.config}': "
                     f"no such file or directory")

    config = ConfigParser()
    config.read(args.config)

    ora_iprscan_url = config["oracle"]["iprscan"]
    iprscan.check_ispro(ora_iprscan_url, args.type, args.status)


def run_clan_update():
    parser = ArgumentParser(description="Update clans and run "
                                        "profile-profile alignments")
    parser.add_argument("config",
                        metavar="config.ini",
                        help="global configuration file")
    parser.add_argument("databases",
                        metavar="DATABASE",
                        nargs="+",
                        help="database(s) to update")
    parser.add_argument("-t", "--threads", type=int,
                        help="number of alignment workers")
    parser.add_argument("-T", "--tempdir",
                        help="directory to use for temporary files")
    parser.add_argument("-v", "--version", action="version",
                        version=f"%(prog)s {__version__}",
                        help="show the version and exit")
    args = parser.parse_args()

    if not os.path.isfile(args.config):
        parser.error(f"cannot open '{args.config}': "
                     f"no such file or directory")

    config = ConfigParser()
    config.read(args.config)

    if not os.path.isfile(config["misc"]["members"]):
        parser.error(f"cannot open '{config['misc']['members']}': "
                     f"no such file or directory")

    options = ConfigParser()
    options.read(config["misc"]["members"])

    ora_interpro_url = config["oracle"]["interpro"]
    db_names = list(set(args.databases))
    databases = interpro.database.get_databases(ora_interpro_url, db_names)

    kwargs = {
        "threads": args.threads,
        "tmpdir": args.tempdir
    }

    for dbname, database in databases.items():
        params = options[dbname]
        if dbname == "cdd":
            update_cdd_clans(ora_interpro_url, database,
                             cddmasters=params["sequences"],
                             cddid=params["summary"],
                             fam2supfam=params["members"],
                             **kwargs)
        elif dbname == "panther":
            update_hmm_clans(ora_interpro_url, database,
                             hmmdb=params["hmm"],
                             **kwargs)
        else:
            update_hmm_clans(ora_interpro_url, database,
                             hmmdb=params["hmm"],
                             source=params["members"],
                             **kwargs)


def run_hmm_update():
    parser = ArgumentParser(description="Update HMMs")
    parser.add_argument("config",
                        metavar="config.ini",
                        help="global configuration file")
    parser.add_argument("databases",
                        metavar="DATABASE",
                        nargs="+",
                        help="database(s) to update")
    parser.add_argument("-v", "--version", action="version",
                        version=f"%(prog)s {__version__}",
                        help="show the version and exit")
    args = parser.parse_args()

    if not os.path.isfile(args.config):
        parser.error(f"cannot open '{args.config}': "
                     f"no such file or directory")

    config = ConfigParser()
    config.read(args.config)

    if not os.path.isfile(config["misc"]["members"]):
        parser.error(f"cannot open '{config['misc']['members']}': "
                     f"no such file or directory")

    options = ConfigParser()
    options.read(config["misc"]["members"])

    db_names = list(set(args.databases))
    ora_interpro_url = config["oracle"]["interpro"]
    databases = interpro.database.get_databases(ora_interpro_url, db_names)

    for dbname, database in databases.items():
        hmmfile = options[dbname]["hmm"]
        mapfile = options[dbname].get("mapping")
        interpro.hmm.update(ora_interpro_url, database, hmmfile, mapfile)


def run_member_db_update():
    parser = ArgumentParser(description="InterPro member database update")
    parser.add_argument("config",
                        metavar="config.ini",
                        help="global configuration file")
    parser.add_argument("databases",
                        metavar="DATABASE",
                        nargs="+",
                        help="database(s) to update")
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
        parser.error(f"cannot open '{args.config}': "
                     f"no such file or directory")

    config = ConfigParser()
    config.read(args.config)

    if not os.path.isfile(config["misc"]["members"]):
        parser.error(f"cannot open '{config['misc']['members']}': "
                     f"no such file or directory")

    options = ConfigParser()
    options.read(config["misc"]["members"])

    db_names = list(set(args.databases))

    ora_interpro_url = config["oracle"]["interpro"]
    ora_iprscan_url = config["oracle"]["iprscan"]
    ora_goa_url = config["oracle"]["goapro"]
    ora_swpread_url = config["oracle"]["swpread"]
    pg_url = config["postgresql"]["pronto"]

    uniprot_version = config["uniprot"]["version"]
    emails = dict(config["emails"])
    pronto_url = config["misc"]["pronto_url"]
    data_dir = config["misc"]["data_dir"]
    lsf_queue = config["misc"]["lsf_queue"]
    workflow_dir = config["misc"]["workflow_dir"]

    databases = interpro.database.get_databases(url=ora_interpro_url,
                                                names=db_names,
                                                expects_new=True)
    mem_updates = []
    non_mem_updates = []
    site_updates = []
    sources = {}
    for dbname, db in databases.items():
        if db.is_member_db or db.is_feature_db:
            # We need a source for signatures
            try:
                props = options[dbname]
            except KeyError:
                parser.error(f"{config['misc']['members']}: "
                             f"missing database '{dbname}'")

            try:
                sig_source = props["signatures"]
            except KeyError:
                sig_source = None

            if not sig_source:
                parser.error(f"{config['misc']['members']}: "
                             f"'signatures' property missing "
                             f"or empty for database '{dbname}'")
            elif db.is_member_db:
                mem_updates.append(db)
                sources[db.identifier] = sig_source
            elif db.is_feature_db:
                non_mem_updates.append(db)
                sources[db.identifier] = sig_source

        if db.has_site_matches:
            site_updates.append(db)

    if not mem_updates and not non_mem_updates and not site_updates:
        parser.error("No database to update")

    tasks = [
        Task(
            fn=iprscan.import_matches,
            args=(ora_iprscan_url,),
            kwargs=dict(databases=mem_updates + non_mem_updates + site_updates,
                        force_import=True, threads=8),
            name="import-matches",
            scheduler=dict(queue=lsf_queue)
        )
    ]

    if mem_updates:
        tasks += [
            Task(
                fn=interpro.signature.add_staging,
                args=(ora_interpro_url, [(db, sources[db.identifier])
                                         for db in mem_updates]),
                name="load-signatures",
                scheduler=dict(queue=lsf_queue)
            ),
            Task(
                fn=interpro.signature.track_signature_changes,
                args=(ora_interpro_url, pg_url, mem_updates, data_dir),
                name="track-changes",
                scheduler=dict(mem=4000, queue=lsf_queue),
                requires=["load-signatures"]
            ),
            Task(
                fn=interpro.signature.delete_obsoletes,
                args=(ora_interpro_url, mem_updates),
                name="delete-obsoletes",
                scheduler=dict(queue=lsf_queue),
                requires=["track-changes"]
            ),
            Task(
                fn=interpro.signature.update_signatures,
                args=(ora_interpro_url,),
                name="update-signatures",
                scheduler=dict(queue=lsf_queue),
                requires=["delete-obsoletes"]
            ),
            Task(
                fn=interpro.match.update_database_matches,
                args=(ora_interpro_url, databases),
                name="update-matches",
                scheduler=dict(queue=lsf_queue),
                requires=["import-matches", "update-signatures"]
            ),
            Task(
                fn=interpro.match.update_variant_matches,
                args=(ora_interpro_url,),
                name="update-varsplic",
                scheduler=dict(queue=lsf_queue),
                requires=["import-matches", "update-signatures"]
            )
        ]

    if non_mem_updates:
        tasks += [
            Task(
                fn=interpro.signature.update_features,
                args=(ora_interpro_url, [(db, sources[db.identifier])
                                         for db in non_mem_updates]),
                name="update-features",
                scheduler=dict(queue=lsf_queue),
                requires=["import-matches"]
            ),
            Task(
                fn=interpro.match.update_database_feature_matches,
                args=(ora_interpro_url, non_mem_updates),
                name="update-fmatches",
                scheduler=dict(queue=lsf_queue),
                requires=["update-features"]
            )
        ]

    if site_updates:
        if mem_updates:
            req = ["import-sites", "update-matches"]
        else:
            req = ["import-sites"]

        tasks += [
            Task(
                fn=iprscan.import_sites,
                args=(ora_iprscan_url,),
                kwargs=dict(databases=site_updates, force_import=True,
                            threads=2),
                name="import-sites",
                scheduler=dict(queue=lsf_queue)
            ),
            Task(
                fn=interpro.match.update_database_site_matches,
                args=(ora_interpro_url, site_updates),
                name="update-sites",
                scheduler=dict(queue=lsf_queue),
                requires=req
            )
        ]

    if mem_updates:
        # Adding Pronto tasks
        after_pronto = set()
        for t in get_pronto_tasks(ora_interpro_url, ora_swpread_url,
                                  ora_goa_url, pg_url, data_dir, lsf_queue):
            # Adding 'pronto-' prefix
            t.name = f"pronto-{t.name}"
            if t.requires:
                t.requires = {f"pronto-{r}" for r in t.requires}
            else:
                # Task without dependency:
                # add one so it's submitted at the end of the protein update
                t.requires = {"update-matches"}

            tasks.append(t)
            after_pronto.add(t.name)

        tasks += [
            Task(
                fn=interpro.report.send_db_update_report,
                args=(ora_interpro_url, pg_url, mem_updates, data_dir,
                      pronto_url, emails),
                name="send-report",
                scheduler=dict(mem=4000, queue=lsf_queue),
                requires=after_pronto
            ),
        ]

    # Base Mundone database on UniProt version and on the name/version
    # of updated member databases
    versions = [uniprot_version]
    for db in sorted(databases.values(), key=lambda k: k.name.lower()):
        versions.append(f"{db.name.lower().replace(' ', '')}{db.version}")

    database = os.path.join(workflow_dir, f"{'_'.join(versions)}.sqlite")
    with Workflow(tasks, dir=workflow_dir, database=database) as wf:
        if wf.run(args.tasks, dry_run=args.dry_run, monitor=not args.detach):
            sys.exit(0)
        else:
            sys.exit(1)


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

    ora_interpro_url = config["oracle"]["interpro"]
    ora_goa_url = config["oracle"]["goapro"]
    ora_swpread_url = config["oracle"]["swpread"]
    pg_url = config["postgresql"]["pronto"]
    uniprot_version = config["uniprot"]["version"]
    data_dir = config["misc"]["data_dir"]
    lsf_queue = config["misc"]["lsf_queue"]
    workflow_dir = config["misc"]["workflow_dir"]

    os.makedirs(data_dir, exist_ok=True)
    tasks = get_pronto_tasks(ora_interpro_url, ora_swpread_url, ora_goa_url,
                             pg_url, data_dir, lsf_queue)

    database = os.path.join(workflow_dir, f"{uniprot_version}.pronto.sqlite")
    with Workflow(tasks, dir=workflow_dir, database=database) as wf:
        if wf.run(args.tasks, dry_run=args.dry_run, monitor=not args.detach):
            sys.exit(0)
        else:
            sys.exit(1)


def run_uniprot_update():
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

    ora_interpro_url = config["oracle"]["interpro"]
    ora_iprscan_url = config["oracle"]["iprscan"]
    ora_uniparc_url = config["oracle"]["uniparc"]
    ora_goa_url = config["oracle"]["goapro"]
    ora_swpread_url = config["oracle"]["swpread"]
    ora_uaread_url = config["oracle"]["uaread"]
    pg_url = config["postgresql"]["pronto"]

    uniprot_version = config["uniprot"]["version"]
    xrefs_dir = config["uniprot"]["xrefs"]

    emails = dict(config["emails"])

    pronto_url = config["misc"]["pronto_url"]
    data_dir = config["misc"]["data_dir"]
    lsf_queue = config["misc"]["lsf_queue"]
    workflow_dir = config["misc"]["workflow_dir"]

    os.makedirs(data_dir, exist_ok=True)
    tasks = [
        # Data from UAREAD
        Task(
            fn=uniparc.update,
            args=(ora_uniparc_url, ora_uaread_url),
            name="update-uniparc",
            scheduler=dict(queue=lsf_queue)
        ),

        # Data from SWPREAD
        Task(
            fn=interpro.taxonomy.refresh_taxonomy,
            args=(ora_interpro_url, ora_swpread_url),
            name="taxonomy",
            scheduler=dict(queue=lsf_queue)
        ),

        # Data from ISPRO
        Task(
            fn=iprscan.import_matches,
            args=(ora_iprscan_url,),
            kwargs=dict(threads=8),
            name="import-matches",
            scheduler=dict(queue=lsf_queue),
            requires=["update-uniparc"]
        ),
        Task(
            fn=iprscan.import_sites,
            args=(ora_iprscan_url,),
            kwargs=dict(threads=2),
            name="import-sites",
            scheduler=dict(queue=lsf_queue),
            requires=["update-uniparc"]
        ),

        # Data from flat files
        Task(
            fn=interpro.protein.track_changes,
            args=(ora_interpro_url,
                  config["uniprot"]["swiss-prot"],
                  config["uniprot"]["trembl"],
                  uniprot_version,
                  config["uniprot"]["date"],
                  data_dir),
            kwargs=dict(tmpdir="/scratch/"),
            name="update-proteins",
            scheduler=dict(mem=4000, queue=lsf_queue, scratch=40000),
        ),

        # Update IPPRO
        Task(
            fn=interpro.protein.delete_obsoletes,
            args=(ora_interpro_url,),
            kwargs=dict(truncate=True),
            name="delete-proteins",
            scheduler=dict(queue=lsf_queue),
            requires=["update-proteins"]
        ),
        Task(
            fn=interpro.protein.check_proteins_to_scan,
            args=(ora_interpro_url,),
            name="check-proteins",
            scheduler=dict(queue=lsf_queue),
            requires=["delete-proteins", "update-uniparc"]
        ),
        Task(
            fn=interpro.match.update_matches,
            args=(ora_interpro_url,),
            name="update-matches",
            scheduler=dict(mem=1000, queue=lsf_queue),
            requires=["check-proteins", "import-matches"]
        ),
        Task(
            fn=interpro.match.update_feature_matches,
            args=(ora_interpro_url,),
            name="update-fmatches",
            scheduler=dict(queue=lsf_queue),
            requires=["update-matches"]
        ),

        # Data for UniProt/SIB
        Task(
            fn=uniprot.exchange.export_sib,
            args=(ora_interpro_url, emails),
            name="export-sib",
            scheduler=dict(queue=lsf_queue),
            requires=["update-matches"]
        ),
        Task(
            fn=uniprot.unirule.report_integration_changes,
            args=(ora_interpro_url, emails),
            name="report-changes",
            scheduler=dict(mem=2000, queue=lsf_queue),
            requires=["update-matches"]
        ),
        Task(
            fn=uniprot.unirule.build_aa_iprscan,
            args=(ora_iprscan_url,),
            name="aa-iprscan",
            scheduler=dict(queue=lsf_queue),
            # Actually depends on import-matches
            requires=["update-matches"]
        ),
        Task(
            fn=uniprot.unirule.build_xref_condensed,
            args=(ora_interpro_url,),
            name="xref-condensed",
            scheduler=dict(queue=lsf_queue),
            requires=["update-matches"]
        ),
        Task(
            fn=uniprot.unirule.build_xref_summary,
            args=(ora_interpro_url,),
            name="xref-summary",
            scheduler=dict(queue=lsf_queue),
            # `report-changes` uses XREF_SUMMARY so we need to wait
            # until it completes before re-creating the table
            requires=["report-changes"]
        ),
        Task(
            fn=uniprot.exchange.export_xrefs,
            args=(ora_interpro_url, xrefs_dir, emails),
            name="export-xrefs",
            scheduler=dict(queue=lsf_queue),
            requires=["xref-summary"]
        ),
        Task(
            fn=uniprot.unirule.ask_to_snapshot,
            args=(ora_interpro_url, emails),
            name="notify-interpro",
            scheduler=dict(queue=lsf_queue),
            requires=["aa-iprscan", "xref-condensed", "xref-summary",
                      "update-fmatches"]
        ),

        # Copy SwissProt descriptions
        Task(
            fn=interpro.signature.export_swissprot_descriptions,
            args=(pg_url, data_dir),
            name="swissprot-de",
            scheduler=dict(queue=lsf_queue),
        ),

        # Update signatures used by UniRule
        Task(
            fn=uniprot.unirule.update_signatures,
            args=(config["uniprot"]["unirule"], ora_interpro_url),
            name="unirule",
            scheduler=dict(queue=lsf_queue),
        )
    ]

    # Adding Pronto tasks
    for t in get_pronto_tasks(ora_interpro_url, ora_swpread_url, ora_goa_url,
                              pg_url, data_dir, lsf_queue):
        # Adding 'pronto-' prefix
        t.name = f"pronto-{t.name}"
        if t.requires:
            t.requires = {f"pronto-{r}" for r in t.requires}
        else:
            # Task without dependency:
            # add some so it's submitted at the end of the protein update
            t.requires = {"swissprot-de", "taxonomy", "unirule",
                          "update-fmatches"}

        tasks.append(t)

    tasks += [
        # Generate and send report to curators
        Task(
            fn=interpro.report.send_prot_update_report,
            args=(ora_interpro_url, pg_url, data_dir, pronto_url, emails),
            name="send-report",
            scheduler=dict(mem=4000, queue=lsf_queue),
            requires=["pronto-annotations", "pronto-proteins-similarities",
                      "pronto-proteins", "pronto-signatures",
                      "pronto-taxonomy"]
        ),

        # Not urgent tasks (can be run after everything else)
        Task(
            fn=interpro.match.update_variant_matches,
            args=(ora_interpro_url,),
            name="update-varsplic",
            scheduler=dict(queue=lsf_queue),
            requires=["import-matches"]
        ),
        Task(
            fn=interpro.match.update_site_matches,
            args=(ora_interpro_url,),
            name="update-sites",
            scheduler=dict(queue=lsf_queue),
            requires=["import-sites", "update-matches"]
        ),
    ]

    database = os.path.join(workflow_dir, f"{uniprot_version}.sqlite")
    with Workflow(tasks, dir=workflow_dir, database=database) as wf:
        if wf.run(args.tasks, dry_run=args.dry_run, monitor=not args.detach):
            sys.exit(0)
        else:
            sys.exit(1)


def update_database():
    parser = ArgumentParser(description="InterPro pre-member database update")
    parser.add_argument("config", metavar="config.ini",
                        help="Configuration file.")
    parser.add_argument("-n", "--name", required=True,
                        help="Name of member database.")
    parser.add_argument("-d", "--date", required=True,
                        help="Release date of member database "
                             "(format: YYYY-MM-DD).")
    parser.add_argument("-v", "--version", required=True,
                        help="Release version of member database.")
    parser.add_argument("-y", "--yes", dest="confirm", action="store_false",
                        help="Do not ask for confirmation.")
    args = parser.parse_args()

    if not os.path.isfile(args.config):
        parser.error(f"cannot open '{args.config}': no such file or directory")

    config = ConfigParser()
    config.read(args.config)

    ora_interpro_url = config["oracle"]["interpro"]
    interpro.database.update_database(url=ora_interpro_url,
                                      name=args.name,
                                      version=args.version,
                                      date=args.date,
                                      confirm=args.confirm)
