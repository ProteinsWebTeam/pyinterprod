#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import json
import logging
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from tempfile import gettempdir

from .. import logger


def _default_steps() -> list:
    return sorted(k for k, v in _get_steps().items() if not v.get("skip"))


def _get_steps() -> dict:
    from . import goa, interpro, uniprot

    return {
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
            "func": interpro.load_signature2protein,
            "requires": ("descriptions", "signatures", "taxa", "terms")
        }
    }


def run(dsn: str, main_user: str, **kwargs):
    level = kwargs.get("level", logging.INFO)
    processes = kwargs.get("processes", 1)
    steps = kwargs.get("steps", _get_steps())
    tmpdir = kwargs.get("tmpdir", gettempdir())

    logger.setLevel(level)

    for name, step in steps.items():
        if name == "signature2protein":
            step["args"] = (main_user, dsn, processes, tmpdir)
        else:
            step["args"] = (main_user, dsn)

    for step in steps.values():
        step["requires"] = list(step.get("requires", []))

    with ThreadPoolExecutor(max_workers=processes) as executor:
        pending = {}
        running = {}
        for name, step in steps.items():
            for req_name in step["requires"]:
                if req_name in steps:
                    pending[name] = step
                    break
            else:
                # no requirement scheduled to run
                running[name] = step

        fs = {}
        for name, step in running.items():
            logger.info("{:<20}running".format(name))
            f = executor.submit(step["func"], *step["args"])
            fs[f] = name

        done = set()
        failed = set()
        while fs or pending:
            for f in as_completed(fs):
                name = fs[f]
                try:
                    f.result()
                except Exception as exc:
                    logger.error("{:<19}failed ({})".format(name, exc))
                    failed.add(name)
                else:
                    logger.info("{:<20}done".format(name))
                    done.add(name)

                del running[name]

                # Look if any pending step can be submitted/cancelled
                num_submitted = 0
                for pend_name in pending.keys():
                    pend_step = pending[pend_name]
                    tmp = []
                    for req_name in pend_step["requires"]:
                        if req_name in failed:
                            # cancel step (one dependency failed)
                            del pending[pend_name]
                            break
                        elif req_name not in done:
                            tmp.append(req_name)
                    else:
                        # No dependency failed
                        # Update list of dependencies still to run
                        pend_name["requires"] = tmp
                        if not tmp:
                            logger.info("{:<20}running".format(pend_name))
                            del pending[pend_name]
                            running[pend_name] = pend_step
                            f = executor.submit(pend_step["func"], *pend_step["args"])
                            fs[f] = pend_name
                            num_submitted += 1

                if num_submitted:
                    break

            fs = {f: n for f, n in fs.items() if n in running}

    if failed:
        raise RuntimeError("one or more step failed")


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
