import os
import sys
from argparse import ArgumentParser
from configparser import ConfigParser, NoOptionError

from mundone import Task, Workflow

from pyinterprod import interpro, interproscan, pronto, uniprot
from pyinterprod.interpro.clan import update_cdd_clans, update_hmm_clans


def get_pronto_tasks(ora_ipr_uri: str,
                     ora_swp_uri: str | None,
                     ora_goa_uri: str,
                     ora_pdbe_uri: str,
                     pg_ipr_uri: str,
                     data_dir: str,
                     temp_dir: str,
                     scheduler: str,
                     queue: str) -> list[Task]:
    """
    Create the list of tasks to update the Pronto database
    :param ora_ipr_uri: connection string of InterPro Oracle database
    :param ora_swp_uri: connection string of Swiss-Prot Oracle database
    :param ora_goa_uri: connection string of GOA Oracle database
    :param ora_pdbe_uri: connection string of PDBe Oracle database
    :param pg_ipr_uri: connection string of Pronto PostgreSQL database
    :param data_dir: path to directory to store/load data files
    :param temp_dir: path to temporary directory for transient files
    :param scheduler: job scheduler
    :param queue: job queue/partition
    """
    names_db = os.path.join(data_dir, "names.sqlite")
    matches_file = os.path.join(data_dir, "matches")
    tasks = [
        Task(
            fn=pronto.match.create_match_table,
            args=(pg_ipr_uri,),
            name="init-matches",
            scheduler=dict(type=scheduler, queue=queue, mem=100, hours=1)
        ),

        # Data from GOAPRO
        Task(
            fn=pronto.goa.import_annotations,
            args=(ora_goa_uri, pg_ipr_uri),
            name="go-terms",
            scheduler=dict(type=scheduler, queue=queue, mem=500, hours=1)
        ),
        Task(
            fn=pronto.goa.import_go_constraints,
            args=(ora_goa_uri, pg_ipr_uri),
            name="go-constraints",
            scheduler=dict(type=scheduler, queue=queue, mem=500, hours=1)
        ),

        # Data from PDBE
        Task(
            fn=pronto.match.import_pdb_matches,
            args=(ora_pdbe_uri, ora_ipr_uri, pg_ipr_uri),
            name="structures",
            scheduler=dict(type=scheduler, queue=queue, mem=1000, hours=6)
        ),
    ]

    if ora_swp_uri:
        # Import data from Swiss-Prot database
        swp_tasks = [
            Task(
                fn=pronto.protein.import_similarity_comments,
                args=(ora_swp_uri, pg_ipr_uri),
                name="proteins-similarities",
                scheduler=dict(type=scheduler, queue=queue, mem=100, hours=3)
            ),
            Task(
                fn=pronto.protein.import_protein_names,
                args=(ora_swp_uri, pg_ipr_uri, names_db),
                kwargs=dict(tmpdir=temp_dir),
                name="proteins-names",
                scheduler=dict(type=scheduler, queue=queue, mem=2000, hours=24)
            ),
            Task(
                fn=pronto.protein.import_protein_pubmed,
                args=(ora_swp_uri, pg_ipr_uri),
                name="proteins-pubmed",
                scheduler=dict(type=scheduler, queue=queue, mem=1000, hours=3)
            ),
            Task(
                fn=pronto.proteome.import_proteomes,
                args=(ora_swp_uri, pg_ipr_uri),
                name="proteomes",
                scheduler=dict(type=scheduler, queue=queue, mem=100, hours=6)
            ),
        ]
    else:
        swp_tasks = []

    return tasks + swp_tasks + [
        # Data from IPPRO
        Task(
            fn=pronto.database.import_databases,
            args=(ora_ipr_uri, pg_ipr_uri),
            name="databases",
            scheduler=dict(type=scheduler, queue=queue, mem=100, hours=1)
        ),
        Task(
            fn=pronto.protein.import_proteins,
            args=(ora_ipr_uri, pg_ipr_uri),
            name="proteins",
            scheduler=dict(type=scheduler, queue=queue, mem=1000, hours=6)
        ),
        Task(
            fn=pronto.match.export,
            args=(ora_ipr_uri, matches_file),
            kwargs=dict(tmpdir=temp_dir),
            name="export-matches",
            scheduler=dict(type=scheduler, queue=queue, mem=4000, hours=24)
        ),
        Task(
            fn=pronto.match.insert_fmatches,
            args=(ora_ipr_uri, pg_ipr_uri),
            name="insert-fmatches",
            scheduler=dict(type=scheduler, queue=queue, mem=1000, hours=3),
            requires=["databases", "init-matches"]
        ),
        Task(
            fn=pronto.match.insert_matches,
            args=(pg_ipr_uri, matches_file),
            kwargs=dict(processes=8),
            name="insert-matches",
            scheduler=dict(type=scheduler, queue=queue, cpu=8, mem=8000,
                           hours=6),
            requires=["databases", "export-matches", "init-matches"]
        ),
        Task(
            fn=pronto.match.finalize_match_table,
            args=(pg_ipr_uri,),
            name="index-matches",
            scheduler=dict(type=scheduler, queue=queue, mem=100, hours=18),
            requires=["insert-fmatches", "insert-matches"]
        ),
        Task(
            fn=pronto.match.insert_signature2protein,
            args=(pg_ipr_uri, names_db, matches_file),
            kwargs=dict(processes=8, tmpdir=temp_dir),
            name="insert-signature2proteins",
            scheduler=dict(type=scheduler, queue=queue, cpu=8, mem=4000,
                           hours=12),
            # We only need proteins-names if importing data from Swiss-Prot
            # However, we always need the latest matches file
            requires=(["export-matches"] +
                      (["proteins-names"] if ora_swp_uri else []))
        ),
        Task(
            fn=pronto.match.finalize_signature2protein,
            args=(pg_ipr_uri,),
            name="index-signature2proteins",
            scheduler=dict(type=scheduler, queue=queue, mem=100, hours=12),
            requires=["insert-signature2proteins"]
        ),
        Task(
            fn=pronto.signature.insert_signatures,
            args=(ora_ipr_uri, pg_ipr_uri, matches_file),
            kwargs=dict(processes=8),
            name="signatures",
            scheduler=dict(type=scheduler, queue=queue, cpu=8, mem=16000,
                           hours=6),
            requires=["databases", "export-matches"]
        ),
        Task(
            fn=pronto.taxon.import_taxonomy,
            args=(ora_ipr_uri, pg_ipr_uri),
            name="taxonomy",
            scheduler=dict(type=scheduler, queue=queue, mem=2000, hours=12),
        ),
        Task(
            fn=pronto.database.set_ready,
            args=(ora_ipr_uri,),
            name="ready",
            requires=["taxonomy", "index-signature2proteins", "index-matches",
                      "proteins",  "go-terms", "go-constraints", "signatures",
                      "structures"] + [t.name for t in swp_tasks]
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
                             cddmasters=params["fasta"],
                             cddid=params["summary"],
                             fam2supfam=params["members"],
                             **kwargs)
        else:
            update_hmm_clans(ora_interpro_uri, database, params["hmm"],
                             **dict(params),
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
    ora_pdbe_uri = config["oracle"]["pdbe"]
    pg_uri = config["postgresql"]["pronto"]

    uniprot_version = config["uniprot"]["version"]
    emails = dict(config["emails"])
    pronto_url = config["misc"]["pronto_url"]
    data_dir = config["misc"]["data_dir"]
    scheduler, queue = parse_scheduler(config["misc"]["scheduler"])
    temp_dir = config["misc"]["temporary_dir"]
    wflow_dir = config["misc"]["workflows_dir"]

    databases = interpro.database.get_databases(uri=ora_interpro_uri,
                                                names=db_names,
                                                expects_new=True)
    member_dbs = []
    feature_dbs = []
    non_ipm_dbs = []
    site_dbs = []
    model_sources = {}
    toad_sources = {}
    go_sources = []
    for dbname, db in databases.items():
        if db.is_member_db or db.is_feature_db:
            props = {}

            # We usually need a source for signatures
            if dbname in options:
                props = options[dbname]
            elif dbname in ("coils", "mobidblt"):
                # Exception for feature databases without data files
                pass
            else:
                parser.error(f"{config['misc']['members']}: "
                             f"missing database '{dbname}'")

            model_sources[db.identifier] = props
            try:
                toad_sources[db.identifier] = options.get("toad", dbname)
            except NoOptionError:
                pass

            if db.analysis_id is None:
                # No analysis ID in ISPRO
                non_ipm_dbs.append(db)
            elif db.is_member_db:
                member_dbs.append(db)
            elif db.is_feature_db:
                feature_dbs.append(db)

            try:
                go_source = props["go-terms"]
            except KeyError:
                if dbname in ("ncbifam", "panther"):
                    parser.error(f"Missing go-terms property for {db.name}")
            else:
                go_sources.append((dbname, go_source))

        if db.has_site_matches:
            site_dbs.append(db)

    if len(member_dbs + feature_dbs + non_ipm_dbs + site_dbs) == 0:
        parser.error("No database to update")

    tasks = []

    if member_dbs or feature_dbs or site_dbs:
        tasks += [
            Task(
                fn=interpro.iprscan.import_matches_or_sites,
                args=(ora_iprscan_uri, "matches"),
                kwargs=dict(databases=member_dbs + feature_dbs + site_dbs,
                            force=True,
                            threads=8),
                name="update-ipm-matches",
                scheduler=dict(type=scheduler, queue=queue, mem=100, hours=48)
            )
        ]
        ipm_dependencies = ["update-ipm-matches"]
    else:
        ipm_dependencies = []

    if member_dbs:
        tasks += [
            Task(
                fn=interpro.signature.add_staging,
                args=(ora_interpro_uri, [(db, model_sources[db.identifier])
                                         for db in member_dbs]),
                name="load-signatures",
                scheduler=dict(type=scheduler, queue=queue, mem=500, hours=1),
            ),
            Task(
                fn=interpro.signature.track_signature_changes,
                args=(ora_interpro_uri, pg_uri, member_dbs, data_dir),
                name="track-changes",
                scheduler=dict(type=scheduler, queue=queue, mem=4000, hours=1),
                requires=["load-signatures"]
            ),
            Task(
                fn=interpro.signature.delete_obsoletes,
                args=(ora_interpro_uri, member_dbs),
                name="delete-obsoletes",
                scheduler=dict(type=scheduler, queue=queue, mem=100, hours=24),
                requires=["track-changes"]
            ),
            Task(
                fn=interpro.signature.update_signatures,
                args=(ora_interpro_uri, go_sources),
                name="update-signatures",
                scheduler=dict(type=scheduler, queue=queue, mem=100, hours=2),
                requires=["delete-obsoletes"]
            ),
            Task(
                fn=interpro.match.update_database_matches,
                args=(ora_interpro_uri, member_dbs),
                name="update-matches",
                scheduler=dict(type=scheduler, queue=queue, mem=100, hours=24),
                requires=ipm_dependencies + ["update-signatures"]
            ),
            Task(
                fn=interpro.match.rebuild_indexes,
                args=(ora_interpro_uri, "MATCH"),
                name="index-matches",
                scheduler=dict(type=scheduler, queue=queue, mem=100, hours=12),
                requires=["update-matches"]
            ),
            Task(
                fn=interpro.match.insert_toad_matches,
                args=(ora_interpro_uri, member_dbs, toad_sources),
                name="update-tmatches",
                scheduler=dict(type=scheduler, queue=queue, mem=24000, hours=24),
                requires=["update-signatures"]
            ),
            Task(
                fn=interpro.match.update_variant_matches,
                args=(ora_interpro_uri,),
                name="update-varsplic",
                scheduler=dict(type=scheduler, queue=queue, mem=100, hours=2),
                requires=ipm_dependencies + ["update-signatures"]
            )
        ]

        for db in member_dbs:
            if db.name.lower() == "pfam":
                props = model_sources[db.identifier]
                tasks += [
                    Task(
                        fn=interpro.signature.contrib.pfam.persist_pfam_a,
                        args=(ora_interpro_uri, props["seed"], props["full"]),
                        name="persist-pfam-a",
                        scheduler=dict(type=scheduler, queue=queue, mem=24000,
                                       hours=6),
                        requires=ipm_dependencies + ["update-signatures"]
                    ),
                    Task(
                        fn=interpro.signature.contrib.pfam.persist_pfam_c,
                        args=(ora_interpro_uri, props["clans"]),
                        name="persist-pfam-c",
                        scheduler=dict(type=scheduler, queue=queue, mem=500,
                                       hours=1),
                        requires=ipm_dependencies + ["update-signatures"]
                    ),
                ]
                break

    feature_dbs += non_ipm_dbs
    if feature_dbs:
        tasks += [
            Task(
                fn=interpro.signature.update_features,
                args=(ora_interpro_uri, [(db, model_sources[db.identifier])
                                         for db in feature_dbs]),
                name="update-features",
                scheduler=dict(type=scheduler, queue=queue, mem=1000, hours=2),
                requires=ipm_dependencies
            ),
            Task(
                fn=interpro.match.update_database_feature_matches,
                args=(ora_interpro_uri, feature_dbs),
                name="update-fmatches",
                scheduler=dict(type=scheduler, queue=queue, mem=500, hours=4),
                requires=["update-features"]
            ),
            Task(
                fn=interpro.match.rebuild_indexes,
                args=(ora_interpro_uri, "FEATURE_MATCH"),
                name="index-fmatches",
                scheduler=dict(type=scheduler, queue=queue, mem=100, hours=12),
                requires=["update-fmatches"]
            )
        ]

    if site_dbs:
        if member_dbs:
            req = ["update-ipm-sites", "update-matches"]
        else:
            req = ["update-ipm-sites"]

        tasks += [
            Task(
                fn=interpro.iprscan.import_matches_or_sites,
                args=(ora_iprscan_uri, "sites"),
                kwargs=dict(databases=site_dbs, force=True, threads=2),
                name="update-ipm-sites",
                scheduler=dict(type=scheduler, queue=queue, mem=100, hours=72),
            ),
            Task(
                fn=interpro.match.update_database_site_matches,
                args=(ora_interpro_uri, site_dbs),
                name="update-sites",
                scheduler=dict(type=scheduler, queue=queue, mem=100, hours=12),
                requires=req
            ),
            Task(
                fn=interpro.match.rebuild_indexes,
                args=(ora_interpro_uri, "SITE_MATCH"),
                name="index-sites",
                scheduler=dict(type=scheduler, queue=queue, mem=100, hours=12),
                requires=["update-sites"]
            )
        ]

    if member_dbs:
        # Adding Pronto tasks
        after_pronto = []
        for t in get_pronto_tasks(ora_interpro_uri, None,
                                  ora_goa_uri, ora_pdbe_uri, pg_uri,
                                  data_dir, temp_dir, scheduler, queue):
            # Adding 'pronto-' prefix
            t.name = f"pronto-{t.name}"
            if t.requires:
                t.requires = {f"pronto-{r}" for r in t.requires}
            else:
                # Task without dependency:
                # add one, so it's submitted at the end of the protein update
                t.requires = {"update-matches"}

            tasks.append(t)
            after_pronto.append(t.name)

        tasks += [
            Task(
                fn=interpro.report.send_db_update_report,
                args=(ora_interpro_uri, pg_uri, member_dbs, data_dir,
                      pronto_url, emails),
                name="send-report",
                scheduler=dict(type=scheduler, queue=queue, mem=4000, hours=4),
                requires=after_pronto
            ),
        ]

    # Base Mundone database on UniProt version and on the name/version
    # of updated member databases
    versions = [uniprot_version]
    for db in sorted(databases.values(), key=lambda k: k.name.lower()):
        versions.append(f"{db.name.lower().replace(' ', '')}{db.version}")

    database = os.path.join(wflow_dir, f"{'_'.join(versions)}.sqlite")
    with Workflow(tasks, dir=wflow_dir, database=database) as wf:
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
    args = parser.parse_args()

    if not os.path.isfile(args.config):
        parser.error(f"cannot open '{args.config}': no such file or directory")

    config = ConfigParser()
    config.read(args.config)

    ora_interpro_uri = config["oracle"]["ipro-interpro"]
    ora_goa_uri = config["oracle"]["unpr-goapro"]
    ora_swpread_uri = config["oracle"]["unpr-swpread"]
    ora_pdbe_uri = config["oracle"]["pdbe"]
    pg_uri = config["postgresql"]["pronto"]
    uniprot_version = config["uniprot"]["version"]
    data_dir = config["misc"]["data_dir"]
    scheduler, queue = parse_scheduler(config["misc"]["scheduler"])
    temp_dir = config["misc"]["temporary_dir"]
    wflow_dir = config["misc"]["workflows_dir"]

    tasks = get_pronto_tasks(ora_interpro_uri, ora_swpread_uri, ora_goa_uri,
                             ora_pdbe_uri, pg_uri, data_dir, temp_dir,
                             scheduler, queue)

    database = os.path.join(wflow_dir, f"{uniprot_version}_pronto.sqlite")
    with Workflow(tasks, dir=wflow_dir, database=database) as wf:
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
    ora_pdbe_uri = config["oracle"]["pdbe"]
    pg_uri = config["postgresql"]["pronto"]

    uniprot_version = config["uniprot"]["version"]
    xrefs_dir = config["uniprot"]["xrefs"]

    emails = dict(config["emails"])

    pronto_url = config["misc"]["pronto_url"]
    data_dir = config["misc"]["data_dir"]
    scheduler, queue = parse_scheduler(config["misc"]["scheduler"])
    temp_dir = config["misc"]["temporary_dir"]
    wflow_dir = config["misc"]["workflows_dir"]

    tasks = [
        # Data from UAREAD
        Task(
            fn=uniprot.uniparc.update_proteins,
            args=(ora_uniparc_uri, ora_uaread_uri),
            kwargs=dict(top_up=True),
            name="update-uniparc-proteins",
            scheduler=dict(type=scheduler, queue=queue, mem=100, hours=6),
        ),
        Task(
            fn=uniprot.uniparc.update_xrefs,
            args=(ora_uniparc_uri, ora_uaread_uri),
            name="update-uniparc-xrefs",
            scheduler=dict(type=scheduler, queue=queue, mem=1000, hours=36),
        ),

        # Data from UniProt
        Task(
            fn=interpro.taxonomy.refresh_taxonomy,
            args=(ora_interpro_uri, ora_swpread_uri),
            name="taxonomy",
            scheduler=dict(type=scheduler, queue=queue, mem=500, hours=1),
        ),
        Task(
            fn=interpro.signature.export_swissprot_descriptions,
            args=(pg_uri, data_dir),
            name="swissprot-de",
            # TODO: review and update
            scheduler=dict(type=scheduler, queue=queue, mem=10000, hours=3),
        ),
        Task(
            fn=uniprot.unirule.update_signatures,
            args=(config["uniprot"]["unirule"], ora_interpro_uri),
            name="unirule",
            scheduler=dict(type=scheduler, queue=queue, mem=500, hours=1),
        ),

        # Data from ISPRO
        Task(
            fn=interpro.iprscan.import_matches_or_sites,
            args=(ora_iprscan_uri, "matches"),
            kwargs=dict(force=True, threads=8),
            name="update-ipm-matches",
            scheduler=dict(type=scheduler, queue=queue, mem=100, hours=72),
            requires=["update-uniparc-proteins"]
        ),
        Task(
            fn=interpro.iprscan.import_matches_or_sites,
            args=(ora_iprscan_uri, "sites"),
            kwargs=dict(force=True, threads=2),
            name="update-ipm-sites",
            scheduler=dict(type=scheduler, queue=queue, mem=100, hours=72),
            requires=["update-uniparc-proteins"]
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
            kwargs=dict(tmpdir=temp_dir),
            name="update-proteins",
            scheduler=dict(type=scheduler, queue=queue, mem=4000, hours=15),
        ),

        # Update IPPRO
        Task(
            fn=interpro.protein.delete_obsoletes,
            args=(ora_interpro_uri,),
            kwargs=dict(truncate=True),
            name="delete-proteins",
            scheduler=dict(type=scheduler, queue=queue, mem=100, hours=96),
            requires=["update-proteins"]
        ),
        Task(
            fn=interpro.protein.check_proteins_to_scan,
            args=(ora_interpro_uri,),
            name="check-proteins",
            scheduler=dict(type=scheduler, queue=queue, mem=100, hours=5),
            requires=["delete-proteins", "update-uniparc-proteins",
                      "update-uniparc-xrefs"]
        ),
        Task(
            fn=interpro.match.update_matches,
            args=(ora_interpro_uri,),
            name="update-matches",
            scheduler=dict(type=scheduler, queue=queue, mem=1000, hours=24),
            requires=["check-proteins", "update-ipm-matches"]
        ),
        Task(
            fn=interpro.match.update_feature_matches,
            args=(ora_interpro_uri,),
            name="update-fmatches",
            scheduler=dict(type=scheduler, queue=queue, mem=100, hours=3),
            requires=["update-matches"]
        ),
        Task(
            fn=interpro.match.update_toad_matches,
            args=(ora_interpro_uri,),
            name="update-tmatches",
            scheduler=dict(type=scheduler, queue=queue, mem=100, hours=3),
            requires=["update-matches"]
        ),
    ]

    # Adding Pronto tasks
    after_pronto = []
    for t in get_pronto_tasks(ora_interpro_uri, ora_swpread_uri, ora_goa_uri,
                              ora_pdbe_uri, pg_uri, data_dir, temp_dir,
                              scheduler, queue):
        # Adding 'pronto-' prefix
        t.name = f"pronto-{t.name}"
        if t.requires:
            t.requires = {f"pronto-{r}" for r in t.requires}
        else:
            # Task without dependency:
            # add some, so it's submitted at the end of the protein update
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
            scheduler=dict(type=scheduler, queue=queue, mem=4000, hours=4),
            requires=after_pronto
        ),

        # Data for UniProt
        Task(
            fn=uniprot.unirule.report_integration_changes,
            args=(ora_interpro_uri, emails),
            name="report-changes",
            scheduler=dict(type=scheduler, queue=queue, mem=2000, hours=1),
            requires=["update-matches"]
        ),
        Task(
            fn=uniprot.aa.create_aa_alignment,
            args=(ora_iprscan_uri,),
            name="aa-alignment",
            scheduler=dict(type=scheduler, queue=queue, mem=100, hours=8),
            # Actually depends on update-ipm-matches
            requires=["update-matches"]
        ),
        Task(
            fn=uniprot.aa.create_aa_iprscan,
            args=(ora_iprscan_uri,),
            name="aa-iprscan",
            scheduler=dict(type=scheduler, queue=queue, mem=100, hours=36),
            # Actually depends on update-ipm-matches, but better to wait
            # until update-matches is over
            requires=["update-matches"]
        ),
        Task(
            fn=uniprot.aa.create_xref_summary,
            args=(ora_interpro_uri,),
            name="xref-summary",
            scheduler=dict(type=scheduler, queue=queue, mem=100, hours=8),
            # `report-changes` uses XREF_SUMMARY, so we need to wait
            # until it completes before re-creating the table
            requires=["report-changes"]
        ),
        Task(
            fn=uniprot.xrefs.export,
            args=(ora_interpro_uri, xrefs_dir, emails),
            name="export-xrefs",
            scheduler=dict(type=scheduler, queue=queue, mem=1000, hours=6),
            requires=["xref-summary"]
        ),
        Task(
            fn=uniprot.aa.create_xref_condensed,
            args=(ora_interpro_uri,),
            name="xref-condensed",
            scheduler=dict(type=scheduler, queue=queue, mem=100, hours=12),
            requires=["update-matches"]
        ),
        Task(
            fn=uniprot.sib.export,
            args=(ora_interpro_uri, emails),
            name="export-sib",
            scheduler=dict(type=scheduler, queue=queue, mem=100, hours=3),
            requires=["xref-condensed"]
        ),
        Task(
            fn=uniprot.aa.export_repr_domains,
            args=(ora_interpro_uri,
                  # File created during the Pronto update
                  os.path.join(data_dir, "matches"),
                  os.path.join(xrefs_dir, "representative-domains.tsv"),
                  emails),
            kwargs=dict(processes=8),
            name="repr-domains",
            scheduler=dict(type=scheduler, queue=queue, cpu=8, mem=4000, hours=48),
            requires=["pronto-export-matches"]
        ),
        Task(
            fn=uniprot.unirule.ask_to_snapshot,
            args=(ora_interpro_uri, emails),
            name="notify-interpro",
            scheduler=dict(type=scheduler, queue=queue, mem=100, hours=1),
            requires=["aa-alignment", "aa-iprscan", "xref-condensed",
                      "xref-summary", "update-fmatches"]
        ),

        # Not urgent tasks (can be run after everything else)
        Task(
            fn=interpro.match.update_variant_matches,
            args=(ora_interpro_uri,),
            name="update-varsplic",
            scheduler=dict(type=scheduler, queue=queue, mem=100, hours=1),
            requires=["update-ipm-matches"]
        ),
        Task(
            fn=interpro.match.update_site_matches,
            args=(ora_interpro_uri,),
            name="update-sites",
            scheduler=dict(type=scheduler, queue=queue, mem=100, hours=6),
            requires=["update-ipm-sites", "update-matches"]
        ),
    ]

    database = os.path.join(wflow_dir, f"{uniprot_version}.sqlite")
    with Workflow(tasks, dir=wflow_dir, database=database) as wf:
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
    parser.add_argument("--by-name", action="store_true",
                        help="Query ISPRO using the database name instead of "
                             "the analysis ID. Useful if a database has been "
                             "renamed in IPPRO (e.g. TIGRFAMs -> NCBIFAM).")
    parser.add_argument("-y", "--yes", dest="confirm", action="store_false",
                        help="Do not ask for confirmation.")
    args = parser.parse_args()

    if not os.path.isfile(args.config):
        parser.error(f"cannot open '{args.config}': no such file or directory")

    config = ConfigParser()
    config.read(args.config)

    ora_interpro_uri = config["oracle"]["ipro-interpro"]
    interpro.database.update_database(uri=ora_interpro_uri,
                                      name=args.name,
                                      version=args.version,
                                      date=args.date,
                                      by_name=args.by_name,
                                      confirm=args.confirm)


def run_interproscan_manager():
    parser = ArgumentParser(description="InterProScan matches calculation")
    parser.add_argument("config", metavar="main.conf",
                        help="Configuration file.")

    subparsers = parser.add_subparsers(dest="mode", help="mode", required=True)

    parser_import = subparsers.add_parser("import",
                                          help="import sequences from UniParc")
    parser_import.add_argument("--max-upi", metavar="UPI",
                               help="import sequences with a UniParc ID lesser "
                                    "or equal to UPI (default: off)")
    parser_import.add_argument("--top-up", action="store_true", default=False,
                               help="import new sequences instead of importing"
                                    " all sequences (default: off)")

    parser_clean = subparsers.add_parser("clean", help="delete obsolete data")
    parser_clean.add_argument("-a", "--analyses", nargs="*", default=[],
                              type=int, help="ID of analyses to clean "
                                             "(default: all)")

    parser_search = subparsers.add_parser("search", help="search sequences")
    parser_search.add_argument("--debug", action="store_true", default=False,
                               help="show debug messages (default: off)")
    parser_search.add_argument("-l", "--list", action="store_true",
                               default=False, help="list active analyses "
                                                   "and exit (default: off)")
    parser_search.add_argument("-a", "--analyses", nargs="*", default=[],
                               type=int, help="ID of analyses to run "
                                              "(default: all)")
    parser_search.add_argument("-e", "--exclude", nargs="*", default=[],
                               type=int, help="ID of analyses to exclude")
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
        interproscan.uniparc.import_sequences(ispro_uri=iscn_uniparc_uri,
                                              uniparc_uri=unpr_uniparc_uri,
                                              top_up=args.top_up,
                                              max_upi=args.max_upi)
    elif args.mode == "clean":
        interproscan.utils.clean_tables(iscn_iprscan_uri, args.analyses)

    elif args.mode == "search":
        if args.list:
            analyses = interproscan.analyses.get_analyses(iscn_iprscan_uri)
            for analysis_id in sorted(analyses,
                                      key=lambda k: (analyses[k]["name"], k)):
                name = analyses[analysis_id]["name"]
                version = analyses[analysis_id]["version"]
                print(f"{analysis_id:>4}\t{name:<30}\t{version}")

            return

        interproscan.utils.rebuild_indexes(uri=iscn_iprscan_uri,
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
        scheduler, queue = parse_scheduler(config["misc"]["scheduler"])

        interproscan.manager.run(uri=iscn_iprscan_uri,
                                 work_dir=config["misc"]["match_calc_dir"],
                                 temp_dir=config["misc"]["match_calc_dir"],
                                 # Default config
                                 job_cpu=job_cpu,
                                 job_mem=job_mem,
                                 job_size=job_size,
                                 job_timeout=job_timeout,
                                 # Custom configs
                                 config=analyses_configs,
                                 # Job scheduler/queue
                                 scheduler=scheduler,
                                 queue=queue,
                                 # Re-run jobs that failed due to memory/time
                                 auto_retry=True,
                                 # Attempts to re-run failed jobs
                                 max_retries=args.max_retries,
                                 # Concurrent jobs
                                 max_running_jobs=args.concurrent_jobs,
                                 # Max jobs submitted per analysis
                                 max_jobs_per_analysis=args.max_jobs,
                                 # Analyses to perform
                                 analyses=args.analyses,
                                 # Analyses to exclude
                                 exclude=args.exclude,
                                 # Debug options
                                 keep_files=args.keep,
                                 debug=args.debug)


def parse_scheduler(value: str) -> tuple[str, str | None]:
    values = value.split(":")
    if len(values) == 2:
        return values[0], values[1]
    elif len(values) == 1:
        return values[0], None
    raise ValueError(value)
