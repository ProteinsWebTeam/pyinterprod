#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
import os
from concurrent import futures
from typing import Union

from . import interprodb, io, sprot, uniprotdb

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

    db = io.ProteinDatabase(dir=dir)
    count = db.insert(interprodb.get_proteins(url), suffix="_old")
    logging.info("InterPro: {} proteins".format(count))

    db.create()
    count = sprot.load(swissprot_path, db.path, "protein")
    logging.info("Swiss-Prot: {} proteins".format(count))

    count = sprot.load(trembl_path, db.path, "protein")
    logging.info("TrEMBL: {} proteins".format(count))

    logging.info("disk space used: {} bytes".format(db.size))

    logging.info("update changes")
    count = interprodb.update_proteins(url, db)
    logging.info("{} proteins updated".format(count))

    logging.info("add new proteins")
    count = interprodb.insert_proteins(url, db)
    logging.info("{} proteins added".format(count))

    logging.info("add proteins to delete")
    count = interprodb.prepare_deletion(url, db)
    logging.info("{} proteins to delete".format(count))

    db.drop()


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
