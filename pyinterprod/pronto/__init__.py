#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import json
import logging
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from tempfile import gettempdir


def _default_steps() -> list:
    return sorted(k for k, v in _get_steps().items() if not v.get("skip"))


def _get_steps() -> dict:
    from . import goa, uniprot
    from .interpro import tables as interpro
    from .. import orautils

    return {
        "clear": {
            "func": orautils.clear_schema,
            "skip": True
        },
        "annotations": {
            "func": goa.load_annotations
        },
        "comments": {
            "func": uniprot.load_comments
        },
        "databases": {
            "func": interpro.load_databases
        },
        "descriptions": {
            "func": uniprot.load_descriptions
        },
        "enzymes": {
            "func": uniprot.load_enzymes
        },
        "matches": {
            "func": interpro.load_matches
        },
        "proteins": {
            "func": interpro.load_proteins
        },
        "publications": {
            "func": goa.load_publications
        },
        "signatures": {
            "func": interpro.load_signatures
        },
        "taxa": {
            "func": interpro.load_taxa
        },
        "terms": {
            "func": goa.load_terms
        },
        "signature2protein": {
            "func": interpro.load_signature2protein
        },
        "copy": {
            "func": interpro.copy_schema
        }
    }


def run(dsn: str, main_user: str, alt_user: str=None, **kwargs):
    level = kwargs.get("level", logging.INFO)
    processes = kwargs.get("processes", 1)
    steps = kwargs.get("steps", _get_steps())
    tmpdir = kwargs.get("tmpdir", gettempdir())

    from .. import logger
    logger.setLevel(level)

    step = steps.pop("clear", None)
    if step is not None:
        logger.info("{:<20}running".format("clear"))
        step["func"](main_user, dsn)
        logger.info("{:<20}done".format("clear"))

    step_s2p = steps.pop("signature2protein", None)
    step_copy = steps.pop("copy", None)

    if steps:
        num_errors = 0
        with ThreadPoolExecutor(max_workers=processes) as executor:
            fs = {}
            running = set()
            for name, step in steps.items():
                step = steps[name]
                f = executor.submit(step["func"], main_user, dsn)
                fs[f] = name
                logger.info("{:<20}running".format(name))
                running.add(name)

            for f in as_completed(fs):
                name = fs[f]
                running.remove(name)
                try:
                    f.result()
                except Exception as exc:
                    logger.error("{:<19}failed ({})".format(name, exc))
                    num_errors += 1
                else:
                    logger.info("{:<20}done".format(name))

                    # TODO: check that 'annotations' is complete as well
                    if name == "descriptions" and step_s2p:
                        # signature2protein can now run
                        fs = {f: name for f, name in fs.items() if name in running}
                        f = executor.submit(step_s2p["func"], main_user, dsn,
                                            processes, tmpdir)
                        fs[f] = "signature2protein"
                        logger.info("{:<20}running".format(fs[f]))
                        break

            for f in as_completed(fs):
                name = fs[f]
                try:
                    f.result()
                except Exception as exc:
                    logger.error("{:<19}failed ({})".format(name, exc))
                    num_errors += 1
                else:
                    logger.info("{:<20}done".format(name))

        if num_errors:
            raise RuntimeError("one or more step failed")
    elif step_s2p:
        logger.info("{:<20}running".format("signature2protein"))
        step_s2p["func"](main_user, dsn, processes=processes, tmpdir=tmpdir)
        logger.info("{:<20}done".format("signature2protein"))

    if step_copy and alt_user:
        logger.info("{:<20}running".format("copy"))
        step_copy["func"](main_user, alt_user, dsn)
        logger.info("{:<20}done".format("copy"))


def main():
    parser = argparse.ArgumentParser(description="Pronto schema update")
    parser.add_argument("config", metavar="CONFIG.JSON",
                        help="config JSON file")
    parser.add_argument("-s", "--steps", nargs="+", choices=_get_steps(),
                        default=_default_steps(),
                        help="steps to run (default: all)")
    parser.add_argument("-t", "--tmp", metavar="DIRECTORY",
                        help="temporary directory", default=gettempdir())
    parser.add_argument("-p", "--processes", type=int, default=1,
                        help="number of processes (default: 1)")
    parser.add_argument("--verbose", action="store_const",
                        const=logging.DEBUG, default=logging.INFO,
                        help="display additional logging messages")
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

    run(dsn=config["databases"]["interpro"]["dsn"],
        main_user=config["databases"]["interpro"]["users"]["pronto_main"],
        alt_user=config["databases"]["interpro"]["users"]["pronto_alt"],
        steps={k: v for k, v in _get_steps().items() if k in args.steps},
        tmpdir=args.tmp,
        processes=args.processes,
        level=args.verbose)
