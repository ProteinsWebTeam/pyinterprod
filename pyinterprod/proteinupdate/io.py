#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sqlite3
from tempfile import mkstemp
from typing import Generator, Tuple, Union


class ProteinDatabase(object):
    def __init__(self, path: Union[str, None]=None,
                 dir: Union[str, None]=None):
        if path:
            self.path = path
            self.temporary = False
        else:
            fd, self.path = mkstemp(dir=dir)
            os.close(fd)
            os.remove(self.path)
            self.temporary = True

    def __del__(self):
        if self.temporary:
            self.drop()

    @property
    def size(self) -> int:
        try:
            return os.path.getsize(self.path)
        except FileNotFoundError:
            return 0

    def create(self, suffix: str=""):
        with sqlite3.connect(self.path) as con:
            con.execute(
                """
                CREATE TABLE IF NOT EXISTS protein{} (
                  identifier TEXT NOT NULL,
                  accession TEXT NOT NULL PRIMARY KEY,
                  is_reviewed INTEGER NOT NULL,
                  length INTEGER NOT NULL,
                  is_fragment INTEGER NOT NULL,
                  taxon_id INTEGER NOT NULL,
                  crc64 TEXT NOT NULL
                )
                """.format(suffix)
            )

    def insert(self, src: Generator, suffix: str="",
               max_items: int=100000) -> int:
        self.create(suffix)

        with sqlite3.connect(self.path) as con:
            count = 0
            items = []

            for protein in src:
                items.append(protein)
                count += 1

                if len(items) == max_items:
                    con.executemany(
                        """
                        INSERT INTO protein{}
                        VALUES (?, ?, ?, ?, ?, ?, ?)
                        """.format(suffix),
                        items
                    )
                    items = []

            if items:
                con.executemany(
                    """
                    INSERT INTO protein{}
                    VALUES (?, ?, ?, ?, ?, ?, ?)
                    """.format(suffix),
                    items
                )
            con.commit()

        return count

    def iter(self, suffix: str="") -> Generator[Tuple, None, None]:
        self.create(suffix)

        with sqlite3.connect(self.path) as con:
            for row in con.execute("SELECT * FROM protein{}".format(suffix)):
                yield row

    def get_deleted(self) -> Generator[str, None, None]:
        with sqlite3.connect(self.path) as con:
            cur = con.execute(
                """
                SELECT accession
                FROM protein_old
                EXCEPT
                SELECT accession
                FROM protein
                """
            )

            for row in cur:
                yield row[0]

    def get_new(self) -> Generator[str, None, None]:
        with sqlite3.connect(self.path) as con:
            cur = con.execute(
                """
                SELECT
                  accession,
                  identifier,
                  is_reviewed,
                  crc64,
                  length,
                  is_fragment,
                  taxon_id
                FROM protein
                WHERE accession NOT IN (
                  SELECT accession
                  FROM protein_old
                )
                """
            )

            for row in cur:
                yield row

    def get_changes(self) -> Generator[Tuple, None, None]:
        with sqlite3.connect(self.path) as con:
            cur = con.execute(
                """
                SELECT
                  accession, p1.identifier, p1.is_reviewed, p1.crc64,
                  p1.length, p1.is_fragment, p1.taxon_id
                FROM protein AS p1
                INNER JOIN protein_old AS p2
                  USING (accession)
                WHERE p1.identifier != p2.identifier
                  OR p1.is_reviewed != p2.is_reviewed
                  OR p1.crc64 != p2.crc64
                  OR p1.length != p2.length
                  OR p1.is_fragment != p2.is_fragment
                  OR p1.taxon_id != p2.taxon_id

                """
            )

            for row in cur:
                yield row

    def drop(self):
        try:
            os.remove(self.path)
        except FileNotFoundError:
            pass
