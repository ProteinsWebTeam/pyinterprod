#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
import os
from multiprocessing import Process
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

    p1 = Process(target=load_proteins_from_flat_files,
                 args=(swissprot_path, trembl_path, new_db))
    p2 = Process(target=load_proteins_from_database,
                 args=(url, old_db))

    p1.start()
    p2.start()

    p1.join()
    p2.join()

    new_db.insert(old_db.iter(), suffix="_old")
    logging.info("databases merged (size needed: {} bytes)".format(
        new_db.size + old_db.size
    ))
    old_db.drop()

    interpro.insert_proteins(url, new_db)
    new_db.drop()

    logging.info("complete")


def delete(url: str):
    interpro.delete_proteins(url, table="MATCH", column="PROTEIN_AC")
