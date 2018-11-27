#!/usr/bin/env python
# -*- coding: utf-8 -*-

import asyncio
import logging
import os
from threading import Thread
from typing import Union

from . import interpro, io, uniprot

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s: %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)


def load_proteins_from_flat_files(swissprot_path: str, trembl_path: str,
                                  database: io.ProteinDatabase):
    count = database.insert(uniprot.read_flat_file(swissprot_path))
    logging.info("Swiss-Prot: {} proteins".format(count))

    count = database.insert(uniprot.read_flat_file(trembl_path))
    logging.info("TrEMBL: {} proteins".format(count))


def load_proteins_from_database(url: str, database: io.ProteinDatabase):
    count = database.insert(interpro.get_proteins(url))
    logging.info("database: {} proteins".format(count))


def update(url: str, swissprot_path: str, trembl_path: str,
           dir: Union[str, None]=None):
    if dir:
        os.makedirs(dir, exist_ok=True)

    old_db = io.ProteinDatabase(dir=dir)
    new_db = io.ProteinDatabase(dir=dir)

    t1 = Thread(target=load_proteins_from_flat_files,
                 args=(swissprot_path, trembl_path, new_db))
    t1 = Thread(target=load_proteins_from_database,
                 args=(url, old_db))

    t1.start()
    t2.start()

    t1.join()
    t2.join()

    new_db.insert(old_db.iter(), suffix="_old")
    logging.info("databases merged (size needed: {} bytes)".format(
        new_db.size + old_db.size
    ))
    old_db.drop()

    interpro.insert_proteins(url, new_db)
    new_db.drop()

    logging.info("complete")


def count(url: str):
    tables = interpro.get_child_tables(url, "INTERPRO", "PROTEIN")
    aws = [interpro.count_table(url, t["name"]) for t in tables]

    loop = asyncio.get_event_loop()
    future = asyncio.gather(*aws)
    res = loop.run_until_complete(future)
    loop.close()
    print(res)


def delete(url: str):
    tables = interpro.get_child_tables(url, "INTERPRO", "PROTEIN")

    for t in tables:
        interpro.toggle_constraint(url, "INTERPRO", t["name"], t["constraint"], False)

    interpro.delete_proteins(url, table="MATCH", column="PROTEIN_AC")

    for t in tables:
        interpro.toggle_constraint(url, "INTERPRO", t["name"], t["constraint"], True)
