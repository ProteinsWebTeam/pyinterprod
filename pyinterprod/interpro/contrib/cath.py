# -*- coding: utf-8 -*-

import re
from typing import List

from .model import Method

_PREFIX = "G3DSA:"
_TYPE_SUPFAM = 'H'
_TYPE_FUNFAM = 'D'


def parse_superfamilies(filepath: str) -> List[Method]:
    """
    Parse the CathNames.txt file distributed with CATH-Gene3D releases

    :param filepath:
    :return:
    """
    signatures = []
    reg = re.compile(r"^(\d\.\d+\.\d+\.\d+)\s+([a-zA-Z0-9]+)\s+:(.+)?$")
    with open(filepath, "rt") as fh:
        for line in fh:
            if line[0] == '#':
                continue

            m = reg.match(line)
            if m is None:
                continue

            supfam, model, name = m.groups()
            accession = f"{_PREFIX}{supfam}"

            m = Method(accession, _TYPE_SUPFAM, description=name)
            signatures.append(m)

    return signatures


def parse_functional_families(filepath: str) -> List[Method]:
    """
    Parse the funfams HMM file distributed with CATH-Gene3D releases

    :param filepath:
    :return:
    """

    signatures = []
    reg_name = re.compile(r"^NAME\s+(\d+\.\d+\.\d+\.\d+)-FF-(\d+)$")
    dsc_name = re.compile(r"^DESC\s+(.+)$")
    with open(filepath, "rt") as fh:
        cnt = 0
        supfam = funfam = desc = None
        for line in fh:
            if line[:2] == "//":
                if cnt:
                    accession = f"{_PREFIX}{supfam}:{funfam}"
                    m = Method(accession, _TYPE_FUNFAM, description=desc)
                    signatures.append(m)

                supfam = funfam = desc = None
                cnt += 1
                continue

            m = reg_name.match(line)
            if m:
                supfam = m.group(1)
                funfam = int(m.group(2))
                continue

            m = dsc_name.match(line)
            if m:
                desc = m.group(1)

    if supfam is not None:
        accession = f"{_PREFIX}{supfam}:{funfam}"
        m = Method(accession, _TYPE_FUNFAM, description=desc)
        signatures.append(m)

    return signatures
