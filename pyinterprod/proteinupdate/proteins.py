#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
import os
from concurrent import futures
from threading import Thread
from typing import Union

from . import interprodb, io, uniprotdb

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s: %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)


def load_proteins_from_flat_files(swissprot_path: str, trembl_path: str,
                                  database: io.ProteinDatabase):
    count = database.insert(uniprotdb.read_flat_file(swissprot_path))
    logging.info("Swiss-Prot: {} proteins".format(count))

    count = database.insert(uniprotdb.read_flat_file(trembl_path))
    logging.info("TrEMBL: {} proteins".format(count))


def load_proteins_from_database(url: str, database: io.ProteinDatabase):
    count = database.insert(interprodb.get_proteins(url))
    logging.info("InterPro: {} proteins".format(count))


def update(url: str, swissprot_path: str, trembl_path: str,
           dir: Union[str, None]=None):
    if dir:
        os.makedirs(dir, exist_ok=True)

    logging.info("loading proteins")
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

    logging.info("merging databases")
    new_db.insert(old_db.iter(), suffix="_old")
    logging.info("disk space used: {} bytes".format(new_db.size +
                                                    old_db.size))
    old_db.drop()

    logging.info("update changes")
    count = interprodb.update_proteins(url, new_db)
    logging.info("{} proteins updated".format(count))

    logging.info("add new proteins")
    count = interprodb.insert_proteins(url, new_db)
    logging.info("{} proteins added".format(count))

    logging.info("add proteins to delete")
    count = interprodb.prepare_deletion(url, new_db)
    logging.info("{} proteins to delete".format(count))

    new_db.drop()


def _count(url: str) -> list:
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
    count = interprodb.count_proteins_to_delete(url)
    logging.info("{} proteins to delete".format(count))

    tables = interprodb.get_child_tables(url, "INTERPRO", "PROTEIN")

    logging.info("disabling referential constraints")
    for t in tables:
        interprodb.toggle_constraint(url, owner="INTERPRO", table=t["name"],
                                     constraint=t["constraint"], enable=False)

    logging.info("deleting proteins")
    n = len(tables) + 1
    all_done = True
    with futures.ThreadPoolExecutor(max_workers=n) as executor:
        future_to_idx = {}
        for i, t in enumerate(tables):
            future = executor.submit(interprodb.delete_proteins, url,
                                     t["name"], t["column"], count)
            future_to_idx[future] = i

        future_to_idx[-1] = executor.submit(interprodb.delete_proteins, url,
                                            "PROTEIN", "PROTEIN_AC", count)

        for future in futures.as_completed(future_to_idx):
            i = future_to_idx[future]
            name = "PROTEIN" if i == -1 else tables[i]["name"]

            if future.done():
                logging.info("{} table done".format(name))
            else:
                logging.info("{} table exited".format(name))
                all_done = False

    if all_done:
        logging.info("enabling referential constraints")
        for t in tables:
            interprodb.toggle_constraint(url, owner="INTERPRO",
                                         table=t["name"],
                                         constraint=t["constraint"],
                                         enable=True)

