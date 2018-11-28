#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
import os
from concurrent import futures
from multiprocessing import Process
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

    p1 = Process(target=load_proteins_from_flat_files,
                 args=(swissprot_path, trembl_path, new_db))
    p2 = Process(target=load_proteins_from_database,
                 args=(url, old_db))

    p1.start()
    p2.start()

    p1.join()
    p2.join()

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


def delete(url: str):
    tables = interprodb.get_tables_with_proteins_to_delete(url)

    if tables:

        logging.info("disabling referential constraints")
        for t in tables:
            try:
                contraint = t["constraint"]
            except KeyError:
                # The table does not have a constraint to disable
                continue
            else:
                interprodb.toggle_constraint(url,
                                             owner=t["owner"],
                                             table=t["name"],
                                             constraint=contraint,
                                             enable=False)

        logging.info("deleting proteins")
        n = len(tables)
        all_done = True
        with futures.ThreadPoolExecutor(max_workers=n) as executor:
            future_to_idx = {}
            for i, t in enumerate(tables):
                future = executor.submit(interprodb.delete_proteins,
                                         url, t["name"], t["column"])
                future_to_idx[future] = i

                for future in futures.as_completed(future_to_idx):
                    i = future_to_idx[future]
                    name = tables[i]["name"]

                    if future.done():
                        logging.info("{} table done".format(name))
                    else:
                        logging.info("{} table exited".format(name))
                        all_done = False

        if all_done:
            logging.info("enabling referential constraints")
            for t in tables:
                try:
                    contraint = t["constraint"]
                except KeyError:
                    # The table does not have a constraint to enable
                    continue
                else:
                    interprodb.toggle_constraint(url,
                                                 owner=t["owner"],
                                                 table=t["name"],
                                                 constraint=contraint,
                                                 enable=True)
        else:
            raise RuntimeError("some tables were not processed")

