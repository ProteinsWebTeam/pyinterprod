#!/usr/bin/env python
# -*- coding: utf-8 -*-


def main():
    import argparse
    import json
    import os
    from datetime import datetime
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

    interpro_url = config["databases"]["interpro"]["interpro"]
    uniparc_url = config["databases"]["interpro"]["uniparc"]
    swissprot_ff = config["flat_files"]["swissprot"]
    trembl_ff = config["flat_files"]["trembl"]

    create_db_links(interpro_url, [
        config["databases"]["iprscan"],
        config["databases"]["uniparc"]
    ])
    proteins.track_changes(interpro_url, swissprot_ff, trembl_ff,
                           dir=args.tmp)
    proteins.update(interpro_url,
                    version=config["release"]["version"],
                    date=datetime.strptime(config["release"]["date"],
                                           "%d-%b-%Y"))
    uniparc.update(uniparc_url, interpro_url)
    proteins.refresh_sequences_to_scan(interpro_url)
    matches.import_ispro(interpro_url)
