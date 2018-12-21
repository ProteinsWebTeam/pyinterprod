#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
from typing import Generator


class Entry(object):
    def __init__(self):
        self.identifier = None
        self.is_reviewed = None
        self.length = None
        self.accession = None
        self.secondary_accessions = []
        self.date = None
        self.is_fragment = None
        self.taxon_id = None
        self.crc64 = None

        self.prog_id = re.compile("^ID\s+(\w+)\s+([a-zA-Z]+);\s+(\d+)",
                                  flags=re.M)
        self.prog_ac = re.compile("^AC\s+(.+;$)", flags=re.M)
        self.prog_ac2 = re.compile("\w+", flags=re.M)
        self.prog_dt = re.compile(
            "DT\s+(\d+-[A-Z]{3}-\d{4}), sequence version \d+\.",
            flags=re.M
        )
        self.prog_de = re.compile("^DE\s+Flags:\s*Fragment;", flags=re.M)
        self.prog_ox = re.compile("^OX\s+NCBI_TaxID=(\d+)", flags=re.M)
        self.prog_sq = re.compile("^SQ.+?(\w+)\s*CRC64;", flags=re.M)

    def update(self, buffer: str) -> tuple:
        self.identifier, status, length = self.prog_id.search(buffer).groups()
        self.is_reviewed = status == "Reviewed"
        self.length = int(length)

        accessions = []
        for m in self.prog_ac.finditer(buffer):
            for acc in self.prog_ac2.findall(m.group(1)):
                accessions.append(acc)
        self.accession, *self.secondary_accessions = accessions

        # date_string = self.prog_dt.search(buffer).group(1)
        # self.date = datetime.strptime(date_string, "%d-%b-%Y")

        self.is_fragment = self.prog_de.search(buffer) is not None
        self.taxon_id = int(self.prog_ox.search(buffer).group(1))

        self.crc64 = self.prog_sq.search(buffer).group(1)

        return (
            self.accession,
            self.identifier,
            1 if self.is_reviewed else 0,
            self.crc64,
            self.length,
            1 if self.is_fragment else 0,
            self.taxon_id
        )


def read_flat_file(path: str) -> Generator:
    with open(path, "rt") as fh:
        buffer = ""
        e = Entry()
        for line in fh:
            if line[:2] == "//":
                yield e.update(buffer)
                buffer = ""
            else:
                buffer += line

        if buffer:  # in case the last line of the file is not //
            yield e.update(buffer)
