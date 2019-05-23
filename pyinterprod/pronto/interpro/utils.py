# -*- coding: utf-8 -*-

import bisect
import gzip
import heapq
import os
import pickle
from multiprocessing import Pool
from tempfile import mkdtemp, mkstemp
from typing import Any, Dict, Generator, Iterable, List, Optional, Union


MIN_OVERLAP = 0.5   # at least 50% of overlap
COMPRESS_LVL = 6


class Organizer(object):
    def __init__(self, keys: List, dir: Optional[str]=None, exists: bool=False):
        if exists:
            self.keys = os.listdir(dir)
            self.dir = dir
        else:
            self.keys = keys
            self.dir = mkdtemp(dir=dir)

        self.buckets = [
            {
                "path": os.path.join(self.dir, str(i+1)),
                "data": {}
            }
            for i in range(len(self.keys))
        ]

        # for the get() method only
        self.key = None
        self.val = None
        self.eof = False

    def __iter__(self):
        for b in self.buckets:
            if os.path.isfile(b["path"]):
                with gzip.open(b["path"], "rb") as fh:
                    while True:
                        try:
                            key, values = pickle.load(fh)
                        except EOFError:
                            break
                        else:
                            yield key, values

    def get(self, key: Union[int, str], default: Optional[Any]=None):
        if self.eof:
            return default
        elif key == self.key:
            return self.val
        elif self.key is None:
            self.key, self.val = next(self)

        while self.key < key:
            try:
                self.key, self.val = next(self)
            except StopIteration:
                self.eof = True
                break

        return self.val if key == self.key else default

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

    def write(self, key: Union[int, str], value: Any):
        i = bisect.bisect_right(self.keys, key)
        if i:
            bucket = self.buckets[i-1]
            with gzip.open(bucket["path"], "ab", COMPRESS_LVL) as fh:
                pickle.dump((key, value), fh)
        else:
            raise KeyError(key)

    def flush(self):
        for b in self.buckets:
            if b["data"]:
                with gzip.open(b["path"], "ab", COMPRESS_LVL) as fh:
                    pickle.dump(b["data"], fh)
                b["data"] = {}

    def merge(self, processes: int=1) -> int:
        size_before = self.size
        if processes > 1:
            with Pool(processes-1) as pool:
                pool.map(self._merge, [b["path"] for b in self.buckets])
        else:
            for b in self.buckets:
                self._merge(b["path"])

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
        for b in self.buckets:
            try:
                os.remove(b["path"])
            except FileNotFoundError:
                pass

        os.rmdir(self.dir)


def merge_organizers(organizers: Iterable[Organizer]):
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


class SignatureComparator(object):
    def __init__(self, dir: Optional[str]=None):
        self.signatures = {}
        fd, self.path = mkstemp(dir=dir)
        os.close(fd)
        os.remove(self.path)

    @property
    def size(self) -> int:
        return os.path.getsize(self.path)

    def sync(self):
        with gzip.open(self.path, "ab", COMPRESS_LVL) as fh:
            pickle.dump(self.signatures, fh)
        self.signatures = {}

    def remove(self):
        os.remove(self.path)

    def __iter__(self):
        with gzip.open(self.path, "rb") as fh:
            while True:
                try:
                    signatures = pickle.load(fh)
                except EOFError:
                    break
                else:
                    for acc_1, s in signatures.items():
                        yield acc_1, s

    def update(self, matches: List[tuple]) -> List[str]:
        signatures = self.prepare(matches)

        # Evaluate overlap between all signatures
        for acc_1 in signatures:
            # number of residues covered by the signature matches
            residues_1 = sum([e - s + 1 for s, e in signatures[acc_1]])

            if acc_1 in self.signatures:
                s = self.signatures[acc_1]
                s["proteins"] += 1
                s["residues"] += residues_1
            else:
                s = self.signatures[acc_1] = {
                    "proteins": 1,
                    "residues": residues_1,
                    "signatures": {}
                }

            for acc_2 in signatures:
                if acc_1 >= acc_2:
                    continue
                elif acc_2 not in s["signatures"]:
                    """
                    * number of collocations
                    * number of overlapping proteins 
                        (at least 1 overlap >= shortest * MIN_OVERLAP)
                    * number of overlapping proteins
                        (residue overlap / signature with least residues >= MIN_OVERLAP)
                    * number of overlapping residues 
                    """
                    s["signatures"][acc_2] = [0, 0, 0, 0]

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
                s["signatures"][acc_2][0] += 1

                # overlapping proteins (based on matches)
                if has_match_overlap:
                    s["signatures"][acc_2][1] += 1

                # overlapping proteins (based on residues)
                if n_residues / min(residues_1, residues_2) >= MIN_OVERLAP:
                    s["signatures"][acc_2][2] += 1

                # overlapping residues
                s["signatures"][acc_2][3] += n_residues

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
