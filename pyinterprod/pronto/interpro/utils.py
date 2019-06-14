# -*- coding: utf-8 -*-

import bisect
import gzip
import heapq
import os
import pickle
import sqlite3
from multiprocessing import Pool
from tempfile import mkdtemp, mkstemp
from typing import *


MIN_OVERLAP = 0.5   # at least 50% of overlap
COMPRESS_LVL = 6


class Organizer(object):
    def __init__(self, keys: List, exists: bool=False, dir: Optional[str]=None,
                 prefix: Optional[str]=None, buffer_size: int=0):
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
        self.buffer_size = 0
        self.buffer_max_size = buffer_size

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

        self.buffer_size += 1
        if self.buffer_size == self.buffer_max_size:
            self.flush()

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
        self.buffer_size = 0

    def merge(self, processes: int=1, fn: Optional[Callable]=None) -> int:
        size_before = self.size
        if processes > 1:
            with Pool(processes-1) as pool:
                pool.map(self._merge, self.buckets, fn)
        else:
            for bucket in self.buckets:
                self._merge(bucket, fn)

        return max(size_before, self.size)

    @staticmethod
    def _merge(path: str, fn: Optional[Callable]=None):
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
                    if fn:
                        pickle.dump((key, fn(data[key])), fh)
                    else:
                        pickle.dump((key, data[key]), fh)

    def remove(self):
        for bucket in self.buckets:
            try:
                os.remove(bucket)
            except FileNotFoundError:
                pass

        os.rmdir(self.dir)


class MatchComparator(object):
    def __init__(self, buffer_size: int=0, dir: Optional[str]=None):
        self.buffer_size = buffer_size
        self.num_items = 0
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
                    for acc_1, val in signatures.items():
                        yield acc_1, val, comparisons[acc_1]

    @property
    def size(self) -> int:
        return os.path.getsize(self.path)

    def sync(self):
        with gzip.open(self.path, "ab", COMPRESS_LVL) as fh:
            pickle.dump((self.signatures, self.comparisons), fh)
        self.signatures = {}
        self.comparisons = {}
        self.num_items = 0

    def remove(self):
        os.remove(self.path)

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
                self.num_items += 1

            for acc_2 in signatures:
                if acc_1 >= acc_2:
                    continue
                elif acc_2 not in self.comparisons[acc_1]:
                    """
                    * number of collocations
                    * number of overlapping proteins
                        (residue overlap / signature with least residues >= MIN_OVERLAP)
                    * number of overlapping residues
                    """
                    self.comparisons[acc_1][acc_2] = [0, 0, 0]
                    self.num_items += 1

                residues_2 = sum([e - s + 1 for s, e in signatures[acc_2]])
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

                # collocation
                self.comparisons[acc_1][acc_2][0] += 1

                # overlapping proteins
                if n_residues / min(residues_1, residues_2) >= MIN_OVERLAP:
                    self.comparisons[acc_1][acc_2][1] += 1

                # overlapping residues
                self.comparisons[acc_1][acc_2][2] += n_residues

        if self.buffer_size and self.num_items >= self.buffer_size:
            self.sync()

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
    def __init__(self, dir: Optional[str]=None, buffer_size: int=0):
        fd, self.filepath = mkstemp(dir=dir)
        os.close(fd)
        os.remove(self.filepath)
        self.con = self.cur = None
        self._ensure_open()
        self.cur.execute(
            """
            CREATE TABLE data (
                id TEXT PRIMARY KEY NOT NULL,
                val TEXT NOT NULL
            )
            """
        )
        self.buffer = {}
        # < 0: unlimited, 0: disabled, > 0: normal (auto-scyn)
        self.buffer_size = buffer_size

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __del__(self):
        self.close()

    def __iter__(self):
        with sqlite3.connect(self.filepath) as con:
            for row in con.execute("SELECT id, val FROM data ORDER BY id"):
                yield row[0], pickle.loads(row[1])

    def __contains__(self, key: str) -> bool:
        found, value = self._getitem(key)
        return found

    def __setitem__(self, key: str, value: Any):
        if self.buffer_size:
            self.buffer[key] = value
            if len(self.buffer) == self.buffer_size:
                self.sync()
        else:
            self._ensure_open()
            self.cur.execute(
                "INSERT OR REPLACE INTO data (id, val) VALUES (?, ?)",
                (key, pickle.dumps(value))
            )
            self.con.commit()

    def __getitem__(self, key: str) -> Any:
        found, value = self._getitem(key)
        if found:
            return value
        raise KeyError(key)

    def __len__(self) -> int:
        if self.buffer_size:
            return len(self.buffer)
        else:
            self._ensure_open()
            self.cur.execute("SELECT count(*) FROM data")
            return cur.fetchone()[0]

    def _ensure_open(self):
        if not self.con:
            self.con = sqlite3.connect(self.filepath)
            self.cur = self.con.cursor()

    def _getitem(self, key: str):
        if key in self.buffer:
            return True, self.buffer[key]

        self._ensure_open()
        self.cur.execute("SELECT val FROM data WHERE id=?", (key,))
        row = self.cur.fetchone()
        if row:
            value = pickle.loads(row[0])
            if self.buffer_size:
                self.buffer[key] = value
            return True, value
        return False, None

    @property
    def size(self) -> int:
        return os.path.getsize(self.filepath)

    def sync(self):
        if not self.buffer:
            return

        self._ensure_open()
        self.cur.executemany(
            "INSERT OR REPLACE INTO data (id, val) VALUES (?, ?)",
            ((k, pickle.dumps(v)) for k, v in self.buffer.items())
        )
        self.con.commit()
        self.buffer = {}

    def close(self):
        if self.con is not None:
            self.sync()
            self.cur.close()
            self.con.close()
            self.con = self.cur = None

    def remove(self):
        self.close()
        os.remove(self.filepath)


def merge_comparators(comparators: Iterable[MatchComparator], remove: bool=False):
    signatures = {}
    comparisons = {}
    for c in comparators:
        for acc_1, val, _comparisons in c:
            if acc_1 in signatures:
                for i, v in enumerate(val):
                    signatures[acc_1][i] += v

                for acc_2, val in _comparisons.items():
                    if acc_2 in comparisons[acc_1]:
                        for i, v in enumerate(val):
                            comparisons[acc_1][acc_2][i] += v
                    else:
                        comparisons[acc_1][acc_2] = val
            else:
                signatures[acc_1] = val
                comparisons[acc_1] = _comparisons

        if remove:
            c.remove()

    return signatures, comparisons


def merge_kvdbs(iterable: Iterable[Kvdb], remove: bool=False):
    items = []
    _key = None
    for key, value in heapq.merge(*iterable, key=lambda x: x[0]):
        if key != _key:
            if _key is not None:
                yield _key, items
            _key = key
            items = []

        items.append(value)

    if _key is not None:
        yield _key, items

    if remove:
        for kvdb in iterable:
            kvdb.remove()


def merge_organizers(iterable: Iterable[Organizer], remove: bool=False):
    _key = None
    items = []
    for key, values in heapq.merge(*iterable):
        if key != _key:
            if _key is not None:
                yield _key, items
            _key = key
            items = []

        items += values

    if _key is not None:
        yield _key, items

    if remove:
        for organizer in iterable:
            organizer.remove()
