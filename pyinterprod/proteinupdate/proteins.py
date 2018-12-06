#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
import os
from concurrent import futures
from typing import Optional

from . import interprodb, io, sprot

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s: %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)


def track_changes(url: str, swissprot_path: str, trembl_path: str,
                  dir: Optional[str]=None):
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


def delete(url: str, truncate_mv: bool=False):
    tables = interprodb.get_tables_with_proteins_to_delete(url)

    if tables:
        to_truncate = []
        to_delete = []
        if truncate_mv:
            for t in tables:
                if t["name"].startswith("MV_"):
                    to_truncate.append(t)
                else:
                    to_delete.append(t)
        else:
            to_delete = tables

        logging.info("disabling referential constraints")
        for t in to_delete:
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
        count = interprodb.count_proteins_to_delete(url)
        n = len(tables)
        all_done = True
        with futures.ThreadPoolExecutor(max_workers=n) as executor:
            future_to_idx = {}
            for i, t in enumerate(to_truncate):
                future = executor.submit(interprodb.truncate_table,
                                         url, t["name"])
                future_to_idx[future] = (True, i)

            for i, t in enumerate(to_delete):
                future = executor.submit(interprodb.delete_proteins,
                                         url, t["name"], t["column"], count)
                future_to_idx[future] = (False, i)

            for future in futures.as_completed(future_to_idx):
                b, i = future_to_idx[future]
                if b:
                    name = to_truncate[i]["name"]
                else:
                    name = to_delete[i]["name"]

                if future.done():
                    logging.info("{} table done".format(name))
                else:
                    logging.info("{} table exited".format(name))
                    all_done = False

        if all_done:
            for t in to_delete:
                try:
                    contraint = t["constraint"]
                except KeyError:
                    # The table does not have a constraint to enable
                    continue
                else:
                    logging.info("enabling: {}.{}.{}".format(t["owner"],
                                                             t["name"],
                                                             contraint))
                    interprodb.toggle_constraint(url,
                                                 owner=t["owner"],
                                                 table=t["name"],
                                                 constraint=contraint,
                                                 enable=True)
        else:
            raise RuntimeError("some tables were not processed")
