#!/usr/bin/env python
# -*- coding: utf-8 -*-


def main():
    import argparse
    import json
    import os
    from tempfile import gettempdir

    from ..orautils import create_db_links
    from . import proteins, matches, uniparc

    parser = argparse.ArgumentParser(description="InterPro protein update")
    parser.add_argument("config", metavar="CONFIG.JSON",
                        help="config JSON file")
    parser.add_argument("-t", "--tmp", metavar="DIRECTORY",
                        help="temporary directory", default=gettempdir())
    args = parser.parse_args()

    try:
        with open(args.config, "rt") as fh:
            config = json.load(fh)
    except FileNotFoundError:
        parser.error("{}: no such file or directory".format(args.config))
    except json.JSONDecodeError:
        parser.error("{}: not a valid JSON file".format(args.config))

    try:
        os.makedirs(args.tmp, exist_ok=True)
    except Exception as e:
        parser.error(e)

    interpro_db = config["databases"]["interpro"]
    dsn = interpro_db["dsn"]
    users = interpro_db["users"]
    create_db_links(users["interpro"], dsn, config["databases"]["others"])

    swissprot_ff = config["flat_files"]["swissprot"]
    trembl_ff = config["flat_files"]["trembl"]
    proteins.insert_new(users["interpro"], dsn, swissprot_ff, trembl_ff, dir=args.tmp)
    proteins.delete_obsolete(users["interpro"], dsn, truncate=True)
    proteins.update_database_info(users["interpro"], dsn,
                                  version=config["release"]["version"],
                                  date=config["release"]["date"])

    uniparc.update(config["databases"]["interpro"]["uniparc"])
    proteins.find_protein_to_refresh(interpro_url)

    # doesnt work on IPTST: use import_mv_iprscan
    # matches.import_ispro(interpro_url)

    matches_dir = config["export"]["matches"]
    matches.prepare_matches(interpro_url)
    matches.check_matches(interpro_url, matches_dir)
    matches.update_matches(interpro_url)
    matches.track_count_changes(interpro_url, matches_dir)
