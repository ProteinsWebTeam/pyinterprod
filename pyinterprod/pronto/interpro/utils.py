# -*- coding: utf-8 -*-

import bisect
import gzip
import heapq
import os
import pickle
import sqlite3
from abc import ABC, abstractmethod
from multiprocessing import Pool
from tempfile import mkdtemp, mkstemp
from typing import Any, Dict, Generator, Iterable, List, Optional, Union

from . import RANKS


MIN_OVERLAP = 0.5   # at least 50% of overlap
COMPRESS_LVL = 6


class Organizer(object):
    def __init__(self, keys: List, exists: bool=False,
                 dir: Optional[str]=None, prefix: Optional[str]=None):
        if exists:
            self.keys = os.listdir(dir)
            self.dir = dir
        else:
            self.keys = keys
            self.dir = mkdtemp(dir=dir, prefix=prefix)

        self.buckets = [
            os.path.join(self.dir, str(i+1))
            for i in range(len(self.keys))
        ]
        self.buffer = {}

        # for the get() method only
        self.it = None
        self.key = None
        self.val = None

    def __iter__(self) -> Generator[tuple, None, None]:
        for bucket in self.buckets:
            if os.path.isfile(bucket):
                with gzip.open(bucket, "rb") as fh:
                    while True:
                        try:
                            key, values = pickle.load(fh)
                        except EOFError:
                            break
                        else:
                            yield key, values

    def get(self, key: Union[int, str], default: Optional[Any]=None):
        if self.it is None:
            self.it = iter(self)
        elif key == self.key:
            return self.val
        elif self.key is None:
            self.key, self.val = next(self.it)

        while self.key < key:
            try:
                self.key, self.val = next(self.it)
            except StopIteration:
                break

        return self.val if key == self.key else default

    @property
    def size(self) -> int:
        size = 0
        for bucket in self.buckets:
            try:
                size += os.path.getsize(bucket)
            except FileNotFoundError:
                continue

        return size

    def add(self, key: Union[int, str], value: Any):
        if key in self.buffer:
            self.buffer[key].append(value)
        else:
            self.buffer[key] = [value]

    def write(self, key: Union[int, str], value: Any):
        i = bisect.bisect_right(self.keys, key)
        if i:
            with gzip.open(self.buckets[i-1], "ab", COMPRESS_LVL) as fh:
                pickle.dump((key, value), fh)
        else:
            raise KeyError(key)

    def flush(self):
        _bucket = None
        data = {}
        for key in sorted(self.buffer):
            i = bisect.bisect_right(self.keys, key)
            if i:
                bucket = self.buckets[i-1]
                if bucket != _bucket:
                    if _bucket:
                        with gzip.open(_bucket, "ab", COMPRESS_LVL) as fh:
                            pickle.dump(data, fh)

                    _bucket = bucket
                    data = {}

                data[key] = self.buffer[key]
            else:
                raise KeyError(key)

        if _bucket:
            with gzip.open(_bucket, "ab", COMPRESS_LVL) as fh:
                pickle.dump(data, fh)

        self.buffer = {}

    def merge(self, processes: int=1) -> int:
        size_before = self.size
        if processes > 1:
            with Pool(processes-1) as pool:
                pool.map(self._merge, self.buckets)
        else:
            for bucket in self.buckets:
                self._merge(bucket)

        return max(size_before, self.size)

    @staticmethod
    def _merge(path: str):
        data = {}
        if os.path.isfile(path):
            with gzip.open(path, "rb") as fh:
                while True:
                    try:
                        chunk = pickle.load(fh)
                    except EOFError:
                        break
                    else:
                        for key, values in chunk.items():
                            if key in data:
                                data[key] += values
                            else:
                                data[key] = values

            with gzip.open(path, "wb", COMPRESS_LVL) as fh:
                for key in sorted(data):
                    pickle.dump((key, data[key]), fh)

    def remove(self):
        for bucket in self.buckets:
            try:
                os.remove(bucket)
            except FileNotFoundError:
                pass

        os.rmdir(self.dir)


def merge_organizers(organizers: Iterable[Organizer]) -> Generator[tuple, None, None]:
    _key = None
    items = []
    for key, values in heapq.merge(*organizers):
        if key != _key:
            if _key is not None:
                yield _key, items
            _key = key
            items = []

        items += values

    if _key is not None:
        yield _key, items


class Comparator(ABC):
    def __init__(self, dir: Optional[str]=None):
        self.signatures = {}
        self.comparisons = {}
        fd, self.path = mkstemp(dir=dir)
        os.close(fd)
        os.remove(self.path)

    def __iter__(self):
        with gzip.open(self.path, "rb") as fh:
            while True:
                try:
                    signatures, comparisons = pickle.load(fh)
                except EOFError:
                    break
                else:
                    for acc_1, cnts in signatures.items():
                        yield acc_1, cnts, comparisons[acc_1]

    @property
    def size(self) -> int:
        return os.path.getsize(self.path)

    def sync(self):
        with gzip.open(self.path, "ab", COMPRESS_LVL) as fh:
            pickle.dump((self.signatures, self.comparisons), fh)
        self.signatures = {}
        self.comparisons = {}

    def remove(self):
        os.remove(self.path)

    @abstractmethod
    def update(self, *args):
        pass


class TaxonomyComparator(Comparator):
    def __init__(self, dir: Optional[str]=None):
        super().__init__(dir)
        self.ranks = {rank: i for i, rank in enumerate(RANKS)}

    def update(self, accessions: List[str], ranks: Iterable[str]):
        for acc_1 in accessions:
            if acc_1 not in self.signatures:
                self.signatures[acc_1] = [0] * len(self.ranks)
                self.comparisons[acc_1] = {}

            for rank in ranks:
                i = self.ranks[rank]
                self.signatures[acc_1][i] += 1

            for acc_2 in accessions:
                if acc_1 >= acc_2:
                    continue
                elif acc_2 not in self.comparisons[acc_1]:
                    self.comparisons[acc_1][acc_2] = [0] * len(self.ranks)

                _accessions.append(acc_2)

            for rank in ranks:
                i = self.ranks[rank]
                for acc_2 in _accessions:
                    self.comparisons[acc_1][acc_2][i] += 1


class MatchComparator(Comparator):
    def update(self, matches: List[tuple]) -> List[str]:
        signatures = self.prepare(matches)

        # Evaluate overlap between all signatures
        for acc_1 in signatures:
            # number of residues covered by the signature matches
            residues_1 = sum([e - s + 1 for s, e in signatures[acc_1]])

            if acc_1 in self.signatures:
                self.signatures[acc_1][0] += 1
                self.signatures[acc_1][1] += residues_1
            else:
                self.signatures[acc_1] = [1, residues_1]
                self.comparisons[acc_1] = {}

            for acc_2 in signatures:
                if acc_1 >= acc_2:
                    continue
                elif acc_2 not in self.comparisons[acc_1]:
                    """
                    * number of collocations
                    * number of overlapping proteins
                        (at least 1 overlap >= shortest * MIN_OVERLAP)
                    * number of overlapping proteins
                        (residue overlap / signature with least residues >= MIN_OVERLAP)
                    * number of overlapping residues
                    """
                    self.comparisons[acc_1][acc_2] = [0, 0, 0, 0]

                residues_2 = sum([e - s + 1 for s, e in signatures[acc_2]])
                has_match_overlap = False
                n_residues = 0
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
                        # o is the number of overlapping residues
                        n_residues += o

                        # Shortest match
                        shortest = min(end_1 - start_1, end_2 - start_2) + 1
                        if o >= shortest * MIN_OVERLAP:
                            has_match_overlap = True

                # collocation
                self.comparisons[acc_1][acc_2][0] += 1

                # overlapping proteins (based on matches)
                if has_match_overlap:
                    self.comparisons[acc_1][acc_2][1] += 1

                # overlapping proteins (based on residues)
                if n_residues / min(residues_1, residues_2) >= MIN_OVERLAP:
                    self.comparisons[acc_1][acc_2][2] += 1

                # overlapping residues
                self.comparisons[acc_1][acc_2][3] += n_residues

        return list(signatures.keys())

    @staticmethod
    def prepare(matches: List[tuple]) -> Dict[str, list]:
        # Group signatures
        signatures = {}
        for acc, pos_start, pos_end in matches:
            if acc in signatures:
                signatures[acc].append((pos_start, pos_end))
            else:
                signatures[acc] = [(pos_start, pos_end)]

        # Sort/merge matches
        for acc in signatures:
            matches = []
            pos_start = pos_end = None
            for start, end in sorted(signatures[acc]):
                if pos_start is None:
                    # Left most match
                    pos_start = start
                    pos_end = end
                elif start > pos_end:
                    """
                      pos_end
                        ----] [----
                              start
                    """
                    matches.append((pos_start, pos_end))
                    pos_start = start
                    pos_end = end
                elif end > pos_end:
                    """
                            pos_end
                        ----]
                          ------]
                                end
                    """
                    pos_end = end

            matches.append((pos_start, pos_end))
            signatures[acc] = matches

        return signatures


class Kvdb(object):
    def __init__(self, dir: Optional[str]=None):
        fd, self.filepath = mkstemp(dir=dir)
        os.close(fd)
        os.remove(self.filepath)
        self.con = sqlite3.connect(self.filepath)
        self.con.execute(
            """
            CREATE TABLE data (
                id TEXT PRIMARY KEY NOT NULL,
                val TEXT NOT NULL
            )
            """
        )
        self.pending = 0

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __del__(self):
        self.close()

    def __setitem__(self, key: str, value: Any):
        self.con.execute(
            "INSERT INTO data (id, val) VALUES (?, ?)",
            (key, pickle.dumps(value))
        )
        self.pending += 1

    def __getitem__(self, key: str) -> dict:
        cur = self.con.execute("SELECT val FROM data WHERE id=?", (key,))
        row = cur.fetchone()
        if row:
            return pickle.loads(row[0])
        else:
            raise KeyError(key)

    def keys(self) -> List[str]:
        return [row[0] for row in self.con.execute("SELECT id FROM data")]

    def commit(self):
        if self.pending:
            self.con.commit()
            self.pending = 0

    def open(self):
        self.close()
        self.con = sqlite3.connect(self.filepath)

    def close(self):
        if self.con is not None:
            self.commit()
            self.con.close()
            self.con = None

    def remove(self):
        self.close()
        os.remove(self.filepath)
