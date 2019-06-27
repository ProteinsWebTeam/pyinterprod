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
    from .. import orautils

    return {
        "clear": {
            "func": orautils.clear_schema,
            "skip": True,
            "pre": True
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
            "func": interpro.load_signature2protein,
            "requires": ("descriptions", "signatures", "taxa", "terms")
        },
        "predictions": {
            "func": interpro.load_predictions,
            "requires": ("signature2protein",)
        },
        "copy": {
            "func": interpro.copy_schema,
            "post": True
        }
    }


def _submit(pool, name, step):
    logger.info("{:<20}running".format(name))
    return pool.submit(step["func"], *step["args"])


def _run(pool: ThreadPoolExecutor, steps: dict, done: set, failed: set):
    running = {}
    pending = {}
    for n, s in steps.items():
        for rn in s["requires"]:
            if rn in steps:
                pending[n] = s
                break
        else:
            # no requirement scheduled to run
            running[n] = s

    fs = {_submit(pool, n, s): n for n, s in running.items()}

    while fs or pending:
        for f in as_completed(fs):
            n = fs[f]
            try:
                f.result()
            except Exception as exc:
                logger.error("{:<19}failed ({})".format(n, exc))
                failed.add(n)
            else:
                logger.info("{:<20}done".format(n))
                done.add(n)

            del running[n]

            # Look if any pending step can be submitted/cancelled
            num_submitted = 0
            for o in pending.keys():
                s = pending[o]
                tmp = []
                for r in s["requires"]:
                    if r in failed:
                        del pending[o]  # cancel step (one dependency failed)
                        break
                    elif r not in done:
                        tmp.append(r)
                else:
                    # No dependency failed
                    s["requires"] = tmp
                    if not tmp:
                        # All dependencies completed
                        del pending[o]
                        running[o] = s
                        f = _submit(pool, o, s)
                        fs[f] = o
                        num_submitted += 1

            if num_submitted:
                break

        fs = {f: n for f, n in fs.items() if n in running}


def run(dsn: str, main_user: str, alt_user: str=None, **kwargs):
    level = kwargs.get("level", logging.INFO)
    processes = kwargs.get("processes", 1)
    steps = kwargs.get("steps", _get_steps())
    tmpdir = kwargs.get("tmpdir", gettempdir())

    logger.setLevel(level)

    for n, s in steps.items():
        if n == "copy":
            s["args"] = (main_user, alt_user, dsn)
        elif n == "signature2protein":
            s["args"] = (main_user, dsn, processes, tmpdir)
        else:
            s["args"] = (main_user, dsn)

    pre = {}
    post = {}
    others = {}
    for n, s in steps.items():
        s["requires"] = list(s.get("requires", []))
        if s.get("pre"):
            pre[n] = s
        elif s.get("post"):
            post[n] = s
        else:
            others[n] = s

    with ThreadPoolExecutor(max_workers=processes) as executor:
        done = set()
        failed = set()
        _run(executor, pre, done, failed)
        if failed:
            raise RuntimeError("one or more step failed")

        _run(executor, others, done, failed)
        if failed:
            raise RuntimeError("one or more step failed")

        _run(executor, post, done, failed)
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
