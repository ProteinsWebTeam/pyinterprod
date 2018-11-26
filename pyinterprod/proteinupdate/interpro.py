#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
import random
import string
from typing import Generator

import cx_Oracle

from .io import ProteinDatabase


def get_proteins(url: str) -> Generator:
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.execute(
        """
        SELECT 
          NAME, PROTEIN_AC, DBCODE, LEN, FRAGMENT, TAX_ID, CRC64
        FROM INTERPRO.PROTEIN
        """
    )

    for row in cur:
        yield (
            row[0],
            row[1],
            1 if row[2] == 'S' else 0,
            row[3],
            1 if row[4] == 'Y' else 0,
            row[5],
            row[6]
        )

    cur.close()
    con.close()


def insert_proteins(url: str, db: ProteinDatabase):
    max_items = 100000

    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.execute("TRUNCATE TABLE INTERPRO.PROTEIN_CHANGES")
    cur.execute("TRUNCATE TABLE INTERPRO.PROTEIN_NEW")

    items = []
    logging.info("track sequence changes")
    for accession in db.get_sequence_changes():
        items.append(('S', accession, accession))

        if len(items) == max_items:
            cur.executemany(
                """
                INSERT INTO INTERPRO.PROTEIN_CHANGES (
                  FLAG, OLD_PROTEIN_AC, NEW_PROTEIN_AC
                ) VALUES (:1, :2, :3)
                """,
                items
            )
            items = []

    logging.info("track annotation changes")
    for accession in db.get_annotation_changes():
        items.append(('A', accession, accession))

        if len(items) == max_items:
            cur.executemany(
                """
                INSERT INTO INTERPRO.PROTEIN_CHANGES (
                  FLAG, OLD_PROTEIN_AC, NEW_PROTEIN_AC
                ) VALUES (:1, :2, :3)
                """,
                items
            )
            items = []

    logging.info("track deleted proteins")
    for accession in db.get_deleted():
        items.append(('D', accession, None))

        if len(items) == max_items:
            cur.executemany(
                """
                INSERT INTO INTERPRO.PROTEIN_CHANGES (
                  FLAG, OLD_PROTEIN_AC, NEW_PROTEIN_AC
                ) VALUES (:1, :2, :3)
                """,
                items
            )
            items = []

    if items:
        cur.executemany(
            """
            INSERT INTO INTERPRO.PROTEIN_CHANGES (
              FLAG, OLD_PROTEIN_AC, NEW_PROTEIN_AC
            ) VALUES (:1, :2, :3)
            """,
            items
        )

    logging.info("track new proteins")
    # TODO: remove datetime when TIMESTAMP is not in the table
    from datetime import datetime
    timestamp = datetime.today()
    items = []
    for row in db.get_new():
        items.append((
            row[0],                     # accession
            row[1],                     # identifier
            'S' if row[2] else 'T',     # dbcode
            'Y' if row[5] else 'N',     # sequence status (fragment)
            row[3],                     # crc64
            row[4],                     # length
            timestamp,                   # timestamp
            row[6]                      # taxon ID
        ))

        if len(items) == max_items:
            cur.executemany(
                """
                INSERT INTO INTERPRO.PROTEIN_NEW
                VALUES (:1, :2, :3, :4, :5, :6, :7, :8)
                """,
                items
            )
            items = []

    if items:
        cur.executemany(
            """
            INSERT INTO INTERPRO.PROTEIN_NEW
            VALUES (:1, :2, :3, :4, :5, :6, :7, :8)
            """,
            items
        )

    con.commit()
    cur.close()
    con.close()


def delete_proteins(url: str, table: str, column: str="PROTEIN_AC"):
    suffix = ''.join(random.choices(string.ascii_uppercase, k=6))

    con = cx_Oracle.connect(url)
    cur = con.cursor()

    logging.info("create INTERPRO.{}_{}".format(table, suffix))
    cur.execute(
        """
        CREATE TABLE INTERPRO.{0}_{1} AS
        SELECT *
        FROM INTERPRO.{0}
        WHERE {2} NOT IN (
          SELECT OLD_PROTEIN_AC
          FROM INTERPRO.PROTEIN_CHANGES
          WHERE FLAG = 'D'
        )
        """.format(table, suffix, column)
    )

    logging.info("delete from INTERPRO.{}".format(table))
    cur.execute(
        """
        DELETE FROM INTERPRO.{}
        WHERE {} IN (
          SELECT OLD_PROTEIN_AC
          FROM INTERPRO.PROTEIN_CHANGES
          WHERE FLAG = 'D'
        )
        """.format(table, column)
    )
    con.commit()

    logging.info("complete")
    cur.close()
    con.close()
