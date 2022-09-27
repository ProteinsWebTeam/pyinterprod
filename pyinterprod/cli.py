import os
import sys
from argparse import ArgumentParser
from configparser import ConfigParser

from mundone import Task, Workflow

from pyinterprod import __version__
from pyinterprod import interpro, interproscan, pronto, uniprot
from pyinterprod.interpro.clan import update_cdd_clans, update_hmm_clans


def get_pronto_tasks(ora_ipr_uri: str, ora_swp_uri: str, ora_goa_uri: str,
                     pg_ipr_uri: str, data_dir: str, lsf_queue: str
                     ) -> list[Task]:
    names_db = os.path.join(data_dir, "names.sqlite")
    matches_file = os.path.join(data_dir, "matches")
    return [
        Task(
            fn=pronto.match.create_match_table,
            args=(pg_ipr_uri,),
            name="init-matches",
            scheduler=dict(mem=100, queue=lsf_queue)
        ),

        # Data from GOAPRO
        Task(
            fn=pronto.goa.import_annotations,
            args=(ora_goa_uri, pg_ipr_uri),
            name="annotations",
            scheduler=dict(mem=500, queue=lsf_queue)
        ),

        # Data from SWPREAD
        Task(
            fn=pronto.protein.import_similarity_comments,
            args=(ora_swp_uri, pg_ipr_uri),
            name="proteins-similarities",
            scheduler=dict(mem=100, queue=lsf_queue),
        ),
        Task(
            fn=pronto.protein.import_protein_names,
            args=(ora_swp_uri, pg_ipr_uri, names_db),
            kwargs=dict(tmpdir="/tmp"),
            name="proteins-names",
            scheduler=dict(mem=2000, tmp=15000, queue=lsf_queue),
        ),
        Task(
            fn=pronto.proteome.import_proteomes,
            args=(ora_swp_uri, pg_ipr_uri),
            name="proteomes",
            scheduler=dict(mem=100, queue=lsf_queue),
        ),

        # Data from IPPRO
        Task(
            fn=pronto.database.import_databases,
            args=(ora_ipr_uri, pg_ipr_uri),
            name="databases",
            scheduler=dict(mem=100, queue=lsf_queue)
        ),
        Task(
            fn=pronto.protein.import_proteins,
            args=(ora_ipr_uri, pg_ipr_uri),
            name="proteins",
            scheduler=dict(mem=100, queue=lsf_queue),
        ),
        Task(
            fn=pronto.match.export,
            args=(ora_ipr_uri, matches_file),
            kwargs=dict(tmpdir="/tmp"),
            name="export-matches",
            scheduler=dict(mem=4000, tmp=100000, queue=lsf_queue)
        ),
        Task(
            fn=pronto.match.insert_fmatches,
            args=(ora_ipr_uri, pg_ipr_uri),
            name="insert-fmatches",
            scheduler=dict(mem=1000, queue=lsf_queue),
            requires=["databases", "init-matches"]
        ),
        Task(
            fn=pronto.match.insert_matches,
            args=(pg_ipr_uri, matches_file),
            kwargs=dict(processes=8),
            name="insert-matches",
            scheduler=dict(cpu=8, mem=8000, queue=lsf_queue),
            requires=["databases", "export-matches", "init-matches"]
        ),
        Task(
            fn=pronto.match.finalize_match_table,
            args=(pg_ipr_uri, ),
            name="index-matches",
            scheduler=dict(mem=100, queue=lsf_queue),
            requires=["insert-fmatches", "insert-matches"]
        ),
        Task(
            fn=pronto.match.insert_signature2protein,
            args=(pg_ipr_uri, names_db, matches_file),
            kwargs=dict(processes=8, tmpdir="/tmp"),
            name="insert-signature2proteins",
            scheduler=dict(cpu=8, mem=4000, tmp=15000, queue=lsf_queue),
            requires=["export-matches", "proteins-names"]
        ),
        Task(
            fn=pronto.match.finalize_signature2protein,
            args=(pg_ipr_uri,),
            name="index-signature2proteins",
            scheduler=dict(mem=100, queue=lsf_queue),
            requires=["insert-signature2proteins"]
        ),
        Task(
            fn=pronto.signature.insert_signatures,
            args=(ora_ipr_uri, pg_ipr_uri, matches_file),
            kwargs=dict(processes=8),
            name="signatures",
            scheduler=dict(cpu=8, mem=16000, queue=lsf_queue),
            requires=["databases", "export-matches"]
        ),
        Task(
            fn=pronto.taxon.import_taxonomy,
            args=(ora_ipr_uri, pg_ipr_uri),
            name="taxonomy",
            scheduler=dict(mem=2000, queue=lsf_queue),
        ),
        Task(
            fn=pronto.database.set_ready,
            args=(pg_ipr_uri,),
            name="ready",
            scheduler=dict(queue=lsf_queue),
            requires=["taxonomy","index-signature2proteins","index-matches", "proteins", 
            "proteins-similarities", "proteomes", "annotations", "signatures"]
        ),
    ]


def check_ispro():
    parser = ArgumentParser(description="Check matches/sites status in ISPRO")
    parser.add_argument("config",
                        metavar="main.conf",
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

    ora_iprscan_uri = config["oracle"]["ipro-iprscan"]
    interpro.iprscan.check_ispro(ora_iprscan_uri, args.type, args.status)


def run_clan_update():
    parser = ArgumentParser(description="Update clans and run "
                                        "profile-profile alignments")
    parser.add_argument("config",
                        metavar="main.conf",
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

    ora_interpro_uri = config["oracle"]["ipro-interpro"]
    db_names = list(set(args.databases))
    databases = interpro.database.get_databases(ora_interpro_uri, db_names)

    kwargs = {
        "threads": args.threads,
        "tmpdir": args.tempdir
    }

    for dbname, database in databases.items():
        params = options[dbname]
        if dbname == "cdd":
            update_cdd_clans(ora_interpro_uri, database,
                             cddmasters=params["sequences"],
                             cddid=params["summary"],
                             fam2supfam=params["members"],
                             **kwargs)
        else:
            update_hmm_clans(ora_interpro_uri, database,
                             hmmdb=params["hmm"],
                             source=params["members"],
                             **kwargs)


def run_hmm_update():
    parser = ArgumentParser(description="Update HMMs")
    parser.add_argument("config",
                        metavar="main.conf",
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
    ora_interpro_uri = config["oracle"]["ipro-interpro"]
    databases = interpro.database.get_databases(ora_interpro_uri, db_names)

    for dbname, database in databases.items():
        hmmfile = options[dbname]["hmm"]
        mapfile = options[dbname].get("mapping")
        interpro.hmm.update(ora_interpro_uri, database, hmmfile, mapfile)


def run_member_db_update():
    parser = ArgumentParser(description="InterPro member database update")
    parser.add_argument("config",
                        metavar="main.conf",
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

    ora_interpro_uri = config["oracle"]["ipro-interpro"]
    ora_iprscan_uri = config["oracle"]["ipro-iprscan"]
    ora_goa_uri = config["oracle"]["unpr-goapro"]
    ora_swpread_uri = config["oracle"]["unpr-swpread"]
    pg_uri = config["postgresql"]["pronto"]

    uniprot_version = config["uniprot"]["version"]
    emails = dict(config["emails"])
    pronto_url = config["misc"]["pronto_url"]
    data_dir = config["misc"]["data_dir"]
    lsf_queue = config["misc"]["lsf_queue"]
    temp_dir = config["misc"]["temporary_dir"]
    wflow_dir = config["misc"]["workflows_dir"]

    databases = interpro.database.get_databases(url=ora_interpro_uri,
                                                names=db_names,
                                                expects_new=True)
    mem_updates = []
    non_mem_updates = []
    site_updates = []
    sig_sources = {}
    go_sources = []
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
                sig_sources[db.identifier] = sig_source
            elif db.is_feature_db:
                non_mem_updates.append(db)
                sig_sources[db.identifier] = sig_source

            try:
                go_source = props["go-terms"]
            except KeyError:
                pass
            else:
                go_sources.append((db, go_source))

        if db.has_site_matches:
            site_updates.append(db)

    if not mem_updates and not non_mem_updates and not site_updates:
        parser.error("No database to update")

    """
    Sources of GO terms for all databases, 
    with a flag if they should be updated
    """
    go_sources = {
        db.identifier: (
            src,
            db in mem_updates or db in non_mem_updates
        )
        for db, src in go_sources
    }

    tasks = [
        Task(
            fn=interpro.iprscan.import_tables,
            args=(ora_iprscan_uri, "matches"),
            kwargs=dict(databases=mem_updates + non_mem_updates + site_updates,
                        force=True,
                        threads=8),
            name="import-ipm-matches",
            scheduler=dict(queue=lsf_queue)
        ),
        Task(
            fn=interpro.iprscan.update_partitions,
            args=(ora_iprscan_uri, "matches"),
            kwargs=dict(databases=mem_updates + non_mem_updates + site_updates,
                        force=True,
                        threads=8),
            name="update-ipm-matches",
            scheduler=dict(queue=lsf_queue),
            requires=["import-ipm-matches"]
        ),
    ]

    if mem_updates:
        tasks += [
            Task(
                fn=interpro.signature.add_staging,
                args=(ora_interpro_uri, [(db, sig_sources[db.identifier])
                                         for db in mem_updates]),
                name="load-signatures",
                scheduler=dict(queue=lsf_queue)
            ),
            Task(
                fn=interpro.signature.track_signature_changes,
                args=(ora_interpro_uri, pg_uri, mem_updates, data_dir),
                name="track-changes",
                scheduler=dict(mem=4000, queue=lsf_queue),
                requires=["load-signatures"]
            ),
            Task(
                fn=interpro.signature.delete_obsoletes,
                args=(ora_interpro_uri, mem_updates),
                name="delete-obsoletes",
                scheduler=dict(queue=lsf_queue),
                requires=["track-changes"]
            ),
            Task(
                fn=interpro.signature.update_signatures,
                args=(ora_interpro_uri, go_sources),
                name="update-signatures",
                scheduler=dict(queue=lsf_queue),
                requires=["delete-obsoletes"]
            ),
            Task(
                fn=interpro.match.update_database_matches,
                args=(ora_interpro_uri, mem_updates),
                name="update-matches",
                scheduler=dict(queue=lsf_queue),
                requires=["update-ipm-matches", "update-signatures"]
            ),
            Task(
                fn=interpro.match.update_variant_matches,
                args=(ora_interpro_uri,),
                name="update-varsplic",
                scheduler=dict(queue=lsf_queue),
                requires=["update-ipm-matches", "update-signatures"]
            )
        ]

    if non_mem_updates:
        tasks += [
            Task(
                fn=interpro.signature.update_features,
                args=(ora_interpro_uri, [(db, sig_sources[db.identifier])
                                         for db in non_mem_updates]),
                name="update-features",
                scheduler=dict(queue=lsf_queue),
                requires=["update-ipm-matches"]
            ),
            Task(
                fn=interpro.match.update_database_feature_matches,
                args=(ora_interpro_uri, non_mem_updates),
                name="update-fmatches",
                scheduler=dict(queue=lsf_queue),
                requires=["update-features"]
            )
        ]

    if site_updates:
        if mem_updates:
            req = ["update-ipm-sites", "update-matches"]
        else:
            req = ["update-ipm-sites"]

        tasks += [
            Task(
                fn=interpro.iprscan.import_tables,
                args=(ora_iprscan_uri, "sites"),
                kwargs=dict(databases=site_updates, force=True, threads=2),
                name="import-ipm-sites",
                scheduler=dict(queue=lsf_queue)
            ),
            Task(
                fn=interpro.iprscan.update_partitions,
                args=(ora_iprscan_uri, "sites"),
                kwargs=dict(databases=site_updates, force=True, threads=2),
                name="update-ipm-sites",
                scheduler=dict(queue=lsf_queue),
                requires=["import-ipm-sites"]
            ),
            Task(
                fn=interpro.match.update_database_site_matches,
                args=(ora_interpro_uri, site_updates),
                name="update-sites",
                scheduler=dict(queue=lsf_queue),
                requires=req
            )
        ]

    if mem_updates:
        # Adding Pronto tasks
        after_pronto = []
        for t in get_pronto_tasks(ora_interpro_uri, ora_swpread_uri,
                                  ora_goa_uri, pg_uri, data_dir, lsf_queue):
            # Adding 'pronto-' prefix
            t.name = f"pronto-{t.name}"
            if t.requires:
                t.requires = {f"pronto-{r}" for r in t.requires}
            else:
                # Task without dependency:
                # add one so it's submitted at the end of the protein update
                t.requires = {"update-matches"}

            tasks.append(t)
            after_pronto.append(t.name)

        tasks += [
            Task(
                fn=interpro.report.send_db_update_report,
                args=(ora_interpro_uri, pg_uri, mem_updates, data_dir,
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

    database = os.path.join(wflow_dir, f"{'_'.join(versions)}.sqlite")
    with Workflow(tasks, dir=temp_dir, database=database) as wf:
        if wf.run(args.tasks, dry_run=args.dry_run, monitor=not args.detach):
            sys.exit(0)
        else:
            sys.exit(1)


def run_pronto_update():
    parser = ArgumentParser(description="InterPro Pronto update")
    parser.add_argument("config",
                        metavar="main.conf",
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

    ora_interpro_uri = config["oracle"]["ipro-interpro"]
    ora_goa_uri = config["oracle"]["unpr-goapro"]
    ora_swpread_uri = config["oracle"]["unpr-swpread"]
    pg_uri = config["postgresql"]["pronto"]
    uniprot_version = config["uniprot"]["version"]
    data_dir = config["misc"]["data_dir"]
    lsf_queue = config["misc"]["lsf_queue"]
    temp_dir = config["misc"]["temporary_dir"]
    wflow_dir = config["misc"]["workflows_dir"]

    tasks = get_pronto_tasks(ora_interpro_uri, ora_swpread_uri, ora_goa_uri,
                             pg_uri, data_dir, lsf_queue)

    database = os.path.join(wflow_dir, f"{uniprot_version}_pronto.sqlite")
    with Workflow(tasks, dir=temp_dir, database=database) as wf:
        if wf.run(args.tasks, dry_run=args.dry_run, monitor=not args.detach):
            sys.exit(0)
        else:
            sys.exit(1)


def run_uniprot_update():
    parser = ArgumentParser(description="InterPro protein update")
    parser.add_argument("config",
                        metavar="main.conf",
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

    ora_interpro_uri = config["oracle"]["ipro-interpro"]
    ora_iprscan_uri = config["oracle"]["ipro-iprscan"]
    ora_uniparc_uri = config["oracle"]["ipro-uniparc"]
    ora_goa_uri = config["oracle"]["unpr-goapro"]
    ora_swpread_uri = config["oracle"]["unpr-swpread"]
    ora_uaread_uri = config["oracle"]["unpr-uaread"]
    pg_uri = config["postgresql"]["pronto"]

    uniprot_version = config["uniprot"]["version"]
    xrefs_dir = config["uniprot"]["xrefs"]

    emails = dict(config["emails"])

    pronto_url = config["misc"]["pronto_url"]
    data_dir = config["misc"]["data_dir"]
    lsf_queue = config["misc"]["lsf_queue"]
    temp_dir = config["misc"]["temporary_dir"]
    wflow_dir = config["misc"]["workflows_dir"]

    tasks = [
        # Data from UAREAD
        Task(
            fn=uniprot.uniparc.update_proteins,
            args=(ora_uniparc_uri, ora_uaread_uri),
            kwargs=dict(top_up=True),
            name="update-uniparc-proteins",
            scheduler=dict(queue=lsf_queue)
        ),
        Task(
            fn=uniprot.uniparc.update_xrefs,
            args=(ora_uniparc_uri, ora_uaread_uri),
            name="update-uniparc-xrefs",
            scheduler=dict(queue=lsf_queue)
        ),

        # Data from SWPREAD
        Task(
            fn=interpro.taxonomy.refresh_taxonomy,
            args=(ora_interpro_uri, ora_swpread_uri),
            name="taxonomy",
            scheduler=dict(queue=lsf_queue)
        ),

        # Data from ISPRO
        Task(
            fn=interpro.iprscan.import_tables,
            args=(ora_iprscan_uri, "matches"),
            kwargs=dict(force=True, threads=8),
            name="import-ipm-matches",
            scheduler=dict(queue=lsf_queue),
            requires=["update-uniparc-proteins"]
        ),
        Task(
            fn=interpro.iprscan.update_partitions,
            args=(ora_iprscan_uri, "matches"),
            kwargs=dict(force=True, threads=8),
            name="update-ipm-matches",
            scheduler=dict(queue=lsf_queue),
            requires=["import-ipm-matches"]
        ),
        Task(
            fn=interpro.iprscan.import_tables,
            args=(ora_iprscan_uri, "sites"),
            kwargs=dict(force=True, threads=2),
            name="import-ipm-sites",
            scheduler=dict(queue=lsf_queue),
            requires=["update-uniparc-proteins"]
        ),
        Task(
            fn=interpro.iprscan.update_partitions,
            args=(ora_iprscan_uri, "sites"),
            kwargs=dict(force=True, threads=2),
            name="update-ipm-sites",
            scheduler=dict(queue=lsf_queue),
            requires=["import-ipm-sites"]
        ),

        # Data from flat files
        Task(
            fn=interpro.protein.track_changes,
            args=(ora_interpro_uri,
                  config["uniprot"]["swiss-prot"],
                  config["uniprot"]["trembl"],
                  uniprot_version,
                  config["uniprot"]["date"],
                  data_dir),
            kwargs=dict(tmpdir="/tmp"),
            name="update-proteins",
            scheduler=dict(mem=4000, queue=lsf_queue, tmp=40000),
        ),

        # Update IPPRO
        Task(
            fn=interpro.protein.delete_obsoletes,
            args=(ora_interpro_uri,),
            kwargs=dict(truncate=True),
            name="delete-proteins",
            scheduler=dict(queue=lsf_queue),
            requires=["update-proteins"]
        ),
        Task(
            fn=interpro.protein.check_proteins_to_scan,
            args=(ora_interpro_uri,),
            name="check-proteins",
            scheduler=dict(queue=lsf_queue),
            requires=["delete-proteins", "update-uniparc-proteins",
                      "update-uniparc-xrefs"]
        ),
        Task(
            fn=interpro.match.update_matches,
            args=(ora_interpro_uri,),
            name="update-matches",
            scheduler=dict(mem=1000, queue=lsf_queue),
            requires=["check-proteins", "update-ipm-matches"]
        ),
        Task(
            fn=interpro.match.update_feature_matches,
            args=(ora_interpro_uri,),
            name="update-fmatches",
            scheduler=dict(queue=lsf_queue),
            requires=["update-matches"]
        ),

        # Data for UniProt/SIB
        Task(
            fn=uniprot.exchange.export_sib,
            args=(ora_interpro_uri, emails),
            name="export-sib",
            scheduler=dict(queue=lsf_queue),
            requires=["xref-condensed"]
        ),
        Task(
            fn=uniprot.unirule.report_integration_changes,
            args=(ora_interpro_uri, emails),
            name="report-changes",
            scheduler=dict(mem=2000, queue=lsf_queue),
            requires=["update-matches"]
        ),
        Task(
            fn=uniprot.unirule.build_aa_alignment,
            args=(ora_iprscan_uri,),
            name="aa-alignment",
            scheduler=dict(queue=lsf_queue),
            # Actually depends on update-ipm-matches
            requires=["update-matches"]
        ),
        Task(
            fn=uniprot.unirule.build_aa_iprscan,
            args=(ora_iprscan_uri,),
            name="aa-iprscan",
            scheduler=dict(queue=lsf_queue),
            # Actually depends on update-ipm-matches
            requires=["update-matches"]
        ),
        Task(
            fn=uniprot.unirule.build_xref_condensed,
            args=(ora_interpro_uri,),
            name="xref-condensed",
            scheduler=dict(queue=lsf_queue),
            requires=["update-matches"]
        ),
        Task(
            fn=uniprot.unirule.build_xref_summary,
            args=(ora_interpro_uri,),
            name="xref-summary",
            scheduler=dict(queue=lsf_queue),
            # `report-changes` uses XREF_SUMMARY so we need to wait
            # until it completes before re-creating the table
            requires=["report-changes"]
        ),
        Task(
            fn=uniprot.exchange.export_xrefs,
            args=(ora_interpro_uri, xrefs_dir, emails),
            name="export-xrefs",
            scheduler=dict(queue=lsf_queue),
            requires=["xref-summary"]
        ),
        Task(
            fn=uniprot.unirule.ask_to_snapshot,
            args=(ora_interpro_uri, emails),
            name="notify-interpro",
            scheduler=dict(queue=lsf_queue),
            requires=["aa-alignment", "aa-iprscan", "xref-condensed",
                      "xref-summary", "update-fmatches"]
        ),

        # Copy SwissProt descriptions
        Task(
            fn=interpro.signature.export_swissprot_descriptions,
            args=(pg_uri, data_dir),
            name="swissprot-de",
            scheduler=dict(queue=lsf_queue),
        ),

        # Update signatures used by UniRule
        Task(
            fn=uniprot.unirule.update_signatures,
            args=(config["uniprot"]["unirule"], ora_interpro_uri),
            name="unirule",
            scheduler=dict(queue=lsf_queue),
        )
    ]

    # Adding Pronto tasks
    after_pronto = []
    for t in get_pronto_tasks(ora_interpro_uri, ora_swpread_uri, ora_goa_uri,
                              pg_uri, data_dir, lsf_queue):
        # Adding 'pronto-' prefix
        t.name = f"pronto-{t.name}"
        if t.requires:
            t.requires = {f"pronto-{r}" for r in t.requires}
        else:
            # Task without dependency:
            # add some so it's submitted at the end of the protein update
            t.requires = {"swissprot-de", "taxonomy", "unirule",
                          "update-fmatches"}

        after_pronto.append(t.name)
        tasks.append(t)

    tasks += [
        # Generate and send report to curators
        Task(
            fn=interpro.report.send_prot_update_report,
            args=(ora_interpro_uri, pg_uri, data_dir, pronto_url, emails),
            name="send-report",
            scheduler=dict(mem=4000, queue=lsf_queue),
            requires=after_pronto
        ),

        # Not urgent tasks (can be run after everything else)
        Task(
            fn=interpro.match.update_variant_matches,
            args=(ora_interpro_uri,),
            name="update-varsplic",
            scheduler=dict(queue=lsf_queue),
            requires=["update-ipm-matches"]
        ),
        Task(
            fn=interpro.match.update_site_matches,
            args=(ora_interpro_uri,),
            name="update-sites",
            scheduler=dict(queue=lsf_queue),
            requires=["import-ipm-sites", "update-matches"]
        ),
    ]

    database = os.path.join(wflow_dir, f"{uniprot_version}.sqlite")
    with Workflow(tasks, dir=temp_dir, database=database) as wf:
        if wf.run(args.tasks, dry_run=args.dry_run, monitor=not args.detach):
            sys.exit(0)
        else:
            sys.exit(1)


def update_database():
    parser = ArgumentParser(description="InterPro pre-member database update")
    parser.add_argument("config", metavar="main.conf",
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

    ora_interpro_uri = config["oracle"]["ipro-interpro"]
    interpro.database.update_database(url=ora_interpro_uri,
                                      name=args.name,
                                      version=args.version,
                                      date=args.date,
                                      confirm=args.confirm)


def run_interproscan_manager():
    parser = ArgumentParser(description="InterProScan matches calculation")
    parser.add_argument("config", metavar="main.conf",
                        help="Configuration file.")

    subparsers = parser.add_subparsers(dest="mode", help="mode", required=True)

    parser_import = subparsers.add_parser("import",
                                          help="import sequences from UniParc")
    parser_import.add_argument("--top-up", action="store_true", default=False,
                               help="if used with --import-sequences: "
                                    "only import sequences not already in "
                                    "the InterProScan database (default: off)")

    parser_clean = subparsers.add_parser("clean", help="delete obsolete data")
    parser_clean.add_argument("-a", "--analyses", nargs="*", default=[],
                              type=int, help="ID of analyses to clean "
                                             "(default: all)")

    parser_search = subparsers.add_parser("search", help="search sequences")
    parser_search.add_argument("--dry-run", action="store_true", default=False,
                               help="show the number of jobs to run and exit "
                                    "(default: off)")
    parser_search.add_argument("-l", "--list", action="store_true",
                               default=False, help="list active analyses "
                                                   "and exit (default: off)")
    parser_search.add_argument("-a", "--analyses", nargs="*", default=[],
                               type=int, help="ID of analyses to run "
                                              "(default: all)")
    parser_search.add_argument("-e", "--exclude", nargs="*", default=[],
                               type=int, help="ID of analyses to exclude")
    parser_search.add_argument("-t", "--threads", type=int, default=8,
                               help="number of job monitoring threads "
                                    "(default: 8)")
    parser_search.add_argument("--concurrent-jobs", type=int, default=1000,
                               help="maximum number of concurrent "
                                    "running jobs (default: 1000)")
    parser_search.add_argument("--max-jobs", type=int, default=-1,
                               help="maximum number of job to run "
                                    "per analysis (default: off)")
    parser_search.add_argument("--max-retries", type=int, default=-0,
                               help="maximum number of attempts to re-run "
                                    "a job after it fails (default: 0)")
    parser_search.add_argument("--keep", choices=["none", "failed", "all"],
                               default="none",
                               help="keep jobs' input/output files "
                                    "(default: none)")
    args = parser.parse_args()

    if not os.path.isfile(args.config):
        parser.error(f"cannot open '{args.config}': no such file or directory")

    config = ConfigParser()
    config.read(args.config)

    if not os.path.isfile(config["misc"]["analyses"]):
        parser.error(f"cannot open '{config['misc']['members']}': "
                     f"no such file or directory")

    iscn_iprscan_uri = config["oracle"]["iscn-iprscan"]
    iscn_uniparc_uri = config["oracle"]["iscn-uniparc"]
    unpr_uniparc_uri = config["oracle"]["unpr-uapro"]

    if args.mode == "import":
        interproscan.database.import_uniparc(ispro_uri=iscn_uniparc_uri,
                                             uniparc_uri=unpr_uniparc_uri,
                                             top_up=args.top_up)
    elif args.mode == "clean":
        interproscan.database.clean_tables(iscn_iprscan_uri, args.analyses)

    elif args.mode == "search":
        if args.list:
            analyses = interproscan.database.get_analyses(iscn_iprscan_uri)
            for analysis_id in sorted(analyses,
                                      key=lambda k: (analyses[k]["name"], k)):
                name = analyses[analysis_id]["name"]
                version = analyses[analysis_id]["version"]
                print(f"{analysis_id:>4}\t{name:<30}\t{version}")

            return

        if not args.dry_run:
            interproscan.database.rebuild_indexes(uri=iscn_iprscan_uri,
                                                  analysis_ids=args.analyses)

        analyses_config = ConfigParser()
        analyses_config.read(config["misc"]["analyses"])

        analyses_configs = {}
        for analysis in analyses_config.sections():
            analyses_configs[analysis] = {}

            for option, value in analyses_config.items(analysis):
                analyses_configs[analysis][option] = int(value)

        # Options for analyses without a custom config
        job_cpu = int(analyses_config["DEFAULT"]["job_cpu"])
        job_mem = int(analyses_config["DEFAULT"]["job_mem"])
        job_size = int(analyses_config["DEFAULT"]["job_size"])
        job_timeout = int(analyses_config["DEFAULT"]["job_timeout"])

        interproscan.manager.run(uri=iscn_iprscan_uri,
                                 work_dir=config["misc"]["match_calc_dir"],
                                 temp_dir=config["misc"]["temporary_dir"],
                                 # Default config
                                 job_cpu=job_cpu,
                                 job_mem=job_mem,
                                 job_size=job_size,
                                 job_timeout=job_timeout,
                                 # Custom configs
                                 config=analyses_configs,
                                 # Performs a dry run
                                 dry_run=args.dry_run,
                                 # LSF queue
                                 lsf_queue=config["misc"]["lsf_queue"],
                                 # Resubmit a job if it fails due to memory
                                 infinite_mem=True,
                                 # Attempts to re-run a failed job (non-memory)
                                 max_retries=args.max_retries,
                                 # Concurrent jobs
                                 max_running_jobs=args.concurrent_jobs,
                                 # Max jobs submitted per analysis
                                 max_jobs_per_analysis=args.max_jobs,
                                 # Number of monitoring threads
                                 pool_threads=args.threads,
                                 # Analyses to perform
                                 analyses=args.analyses,
                                 # Analyses to exclude
                                 exclude=args.exclude,
                                 # Debug options
                                 keep_files=args.keep)
