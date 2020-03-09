# -*- coding: utf-8 -*-

import os
import sqlite3

from pyinterprod import logger, sprot


def load(swissp: str, trembl: str, database: str):
    logger.info("loading Swiss-Prot/TrEMBL data")
    try:
        os.remove(database)
    except FileNotFoundError:
        pass
    con = sqlite3.connect(database)
    con.execute(
        """
        CREATE TABLE protein (
          accession TEXT NOT NULL PRIMARY KEY,
          identifier TEXT NOT NULL,
          is_reviewed INTEGER NOT NULL,
          crc64 TEXT NOT NULL,
          length INTEGER NOT NULL,
          is_fragment INTEGER NOT NULL,
          taxon_id INTEGER NOT NULL
        )
        """
    )
    con.close()

    swissp_cnt = sprot.load(swissp, database, "protein")
    logger.info(f"Swiss-Prot: {swissp_cnt} entries")

    trembl_cnt = sprot.load(trembl, database, "protein")
    logger.info(f"Swiss-Prot: {trembl_cnt} entries")
