# -*- coding: utf-8 -*-

import bisect
import heapq
import os
import pickle
from multiprocessing import Pool
from tempfile import mkdtemp, mkstemp
from typing import Any, Dict, Generator, Iterable, List, Optional, Union


MIN_OVERLAP = 0.5   # at least 50% of overlap


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
            with open(path, "rb") as fh:
                while True:
                    try:
                        chunk = pickle.load(fh)
                    except EOFError:
                        break
                    else:
                        for key, value in chunk.items():
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

    def from_organizers(self, organizers: Iterable):
        _key = None
        items = []
        for key, value in heapq.merge(*organizers):
            if key != _key:
                for item in items:
                    self.add(_key, item)

                _key = key
                items = []
                self.flush()

            items.append(value)

        for item in items:
            self.add(_key, item)
        self.flush()


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
        with open(self.path, "ab") as fh:
            pickle.dump(self.signatures, fh)
        self.signatures = {}

    def remove(self):
        os.remove(self.path)

    def __iter__(self):
        with open(self.path, "rb") as fh:
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
            residues = sum([e - s + 1 for s, e in signatures[acc_1]])

            if acc_1 in self.signatures:
                s = self.signatures[acc_1]
                s["proteins"] += 1
                s["residues"] += residues
            else:
                s = self.signatures[acc_1] = {
                    "proteins": 1,
                    "residues": residues,
                    "signatures": {}
                }

            for acc_2 in signatures:
                if acc_1 >= acc_2:
                    continue
                elif acc_2 not in s["signatures"]:
                    # collocations, protein overlaps, residue overlaps
                    s["signatures"][acc_2] = [0, 0, 0]

                s["signatures"][acc_2][0] += 1
                num_match_overlaps = 0
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
                        s["signatures"][acc_2][2] += o

                        # Shorted match
                        shortest = min(end_1 - start_1, end_2 - start_2) + 1
                        if o >= shortest * MIN_OVERLAP:
                            num_match_overlaps += 1

                if num_match_overlaps:
                    """
                    acc_1 and acc_2 overlaps
                    (by at least 50% of the shortest match)
                    at least once in this protein
                    """
                    s["signatures"][acc_2][1] += 1

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
