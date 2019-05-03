# -*- coding: utf-8 -*-

import bisect
import hashlib
import os
import pickle
from multiprocessing import Queue
from tempfile import mkdtemp
from typing import Any, List, Optional, Union

import cx_Oracle

from ... import orautils

MAX_GAP = 20


class Organizer(object):
    def __init__(self, keys: List, tmpdir: Optional[str]=None):
        self.keys = keys
        self.dir = mkdtemp(dir=tmpdir)
        self.buckets = [
            {
                "path": os.path.join(self.dir, str(i+1)),
                "data": {}
            }
            for i in range(len(keys))
        ]

    def __iter__(self):
        for b in self.buckets:
            with open(b["path"], "rb") as fh:
                while True:
                    try:
                        key, value = pickle.load(fh)
                    except EOFError:
                        break
                    else:
                        yield key, value

    @property
    def size(self) -> int:
        size = 0
        for b in self.buckets:
            try:
                size += os.path.getsize(b["path"])
            except FileNotFoundError:
                continue

        return size

    def add(self, key: Union[int, str], value: Any):
        i = bisect.bisect_right(self.keys, key)
        if i:
            bucket = self.buckets[i-1]
            if key in bucket["data"]:
                bucket["data"][key].append(value)
            else:
                bucket["data"][key] = [value]
        else:
            raise KeyError(key)

    def flush(self):
        for b in self.buckets:
            if b["data"]:
                with open(b["path"], "ab") as fh:
                    pickle.dump(b["data"], fh)
                b["data"] = {}

    def merge(self) -> int:
        size_before = self.size
        for b in self.buckets:
            self._merge(b["path"])

        return max(size_before, self.size)

    @staticmethod
    def _merge(path: str):
        data = {}
        if os.path.isfile(path):
            with open(path, "rb") as fh:
                while True:
                    try:
                        chunk = pickle.load(fh)
                    except EOFError:
                        break
                    else:
                        for key, value in chunk:
                            if key in data:
                                data[key] += value
                            else:
                                data[key] = value

        with open(path, "wb") as fh:
            for key in sorted(data):
                pickle.dump((key, data[key]), fh)

    def remove(self):
        for b in self.buckets:
            try:
                os.remove(b["path"])
            except FileNotFoundError:
                pass

        os.rmdir(self.dir)


def consume_proteins(user: str, dsn: str, task_queue: Queue, done_queue: Queue,
                     tmpdir: Optional[str]=None, bucket_size: int=100):
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    cur.execute(
        """
        SELECT METHOD_AC, DBCODE, SIG_TYPE
        FROM {}.METHOD
        ORDER BY METHOD_AC
        """.format(owner)
    )
    keys = []
    signatures = {}
    for i, (acc, dbcode, _type) in enumerate(cur):
        signatures[acc] = (dbcode, _type)
        if not i % 1000:
            keys.append(acc)
        i += 1

    names = Organizer(keys, tmpdir)
    taxa = Organizer(keys, tmpdir)
    table = orautils.TablePopulator(
        con=con,
        query="INSERT /*+ APPEND */ "
              "INTO {}.METHOD2PROTEIN "
              "VALUES (:1, :2, :3, :4, :5, :6, :7)".format(owner),
        autocommit=True
    )
    residues = {}
    overlaps = {}
    for chunk in iter(task_queue, None):
        for acc, dbcode, length, descid, leftnum, matches in chunk:
            md5 = hash_protein(matches)
            compare_matches(matches, residues, overlaps)


def hash_protein(matches: List[tuple]) -> str:
    # flatten matches
    locations = []
    for acc, pos_start, pos_end in matches:
        locations.append((pos_start, acc))
        locations.append((pos_end, acc))

    """
    Evaluate the protein's match structure,
        i.e. how signatures match the proteins

    -----------------------------   Protein
     <    >                         Signature 1
       <    >                       Signature 2
                  < >               Signature 3

    Flattened:
    -----------------------------   Protein
     < <  > >     < >
     1 2  1 2     3 3

    Structure, with '-' representing a "gap"
        (more than N bp between two positions):
    1212-33
    """

    # Sort locations by position
    locations.sort()

    """
    Do not set the offset to 0, but to the first position:
    if two proteins have the same structure,
    but the first position of one protein is > max_gap
    while the first position of the other protein is <= max_gap,
    a gap will be used for the first protein and not for the other,
    which will results in two different structures
    """
    offset = locations[0][0]

    # overall match structure
    structure = []
    # close signatures (less than max_gap between two positions)
    signatures = []

    for pos, acc in locations:
        if pos > offset + MAX_GAP:
            for _pos, _ac in signatures:
                structure.append(_ac)

            signatures = []
            structure.append('')  # add a gap

        offset = pos
        signatures.append((pos, acc))

    for _pos, _ac in signatures:
        structure.append(_ac)

    return hashlib.md5(
        '/'.join(structure).encode("utf-8")
    ).hexdigest()


def compare_matches(matches: List[tuple], residues: dict, overlaps: dict):
    # Group signatures
    signatures = {}
    for acc, pos_start, pos_end in matches:
        if acc in signatures:
            signatures[acc].append((pos_start, pos_end))
        else:
            signatures[acc] = [(pos_start, pos_end)]

    # Sort matches by position
    for matches in signatures.values():
        matches.sort()

    # Evaluate overlap between all signatures
    for acc_1 in signatures:
        # number of residues covered by the signature matches
        _residues = sum([e - s + 1 for s, e in signatures[acc_1]])

        if acc_1 in residues:
            residues[acc_1] += _residues
        else:
            residues[acc_1] = _residues
            overlaps[acc_1] = {}

        for acc_2 in signatures:
            if acc_1 >= acc_2:
                continue
            elif acc_2 not in overlaps[acc_1]:
                overlaps[acc_1][acc_2] = 0

            i = 0
            start_2, end_2 = signatures[acc_2][i]
            for start_1, end_1 in signatures[acc_1]:
                while end_2 < start_1:
                    i += 1
                    try:
                        start_2, end_2 = signatures[acc_2][i]
                    except IndexError:
                        break

                o = min(end_1, end_2) - max(start_1, start_2) + 1
                if o > 0:
                    overlaps[acc_1][acc_2] += o
