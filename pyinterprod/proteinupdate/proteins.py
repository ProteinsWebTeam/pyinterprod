#!/usr/bin/env python
# -*- coding: utf-8 -*-

import concurrent.futures
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
    t2 = Thread(target=load_proteins_from_database,
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

    logging.info("track changes")
    interpro.insert_proteins(url, new_db)
    new_db.drop()

    logging.info("evaluate tables")
    count(url)

    logging.info("complete")


def count(url: str) -> list:
    tables = interpro.get_child_tables(url, "INTERPRO", "PROTEIN")
    n = len(tables)
    with concurrent.futures.ThreadPoolExecutor(max_workers=n) as executor:
        futures_to_index = {
            executor.submit(
                interpro.count_rows_to_delete,
                url,
                t["name"],
                t["column"]
            ): i
            for i, t in enumerate(tables)
        }

        for future in concurrent.futures.as_completed(futures_to_index):
            i = futures_to_index[future]
            tables[i]["count"] = future.result()

    return tables


def delete(url: str):
    tables = interpro.get_child_tables(url, "INTERPRO", "PROTEIN")

    for t in tables:
        interpro.toggle_constraint(url, "INTERPRO", t["name"], t["constraint"], False)

    interpro.delete_proteins(url, table="MATCH", column="PROTEIN_AC")

    for t in tables:
        interpro.toggle_constraint(url, "INTERPRO", t["name"], t["constraint"], True)
