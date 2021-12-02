# -*- coding: utf-8 -*-

import re
from typing import List

from .common import Method

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
    reg = re.compile(r"^(\d\.\d+\.\d+\.\d+)\s+([a-zA-Z0-9]+)\s+:(.*)$")
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
    Parse the FunFam HMM file distributed with CATH-Gene3D releases.
    Version 4.3.0: http://download.cathdb.info/cath/releases/latest-release/sequence-data/funfam-hmm3.lib.gz

    :param filepath:
    :return:
    """

    signatures = []
    reg_name = re.compile(r"^NAME\s+(.+)$")
    with open(filepath, "rt") as fh:
        cnt = 0
        funfam = None
        for line in fh:
            if line[:2] == "//":
                signatures.append(Method(funfam, _TYPE_FUNFAM))
                funfam = None
                continue

            m = reg_name.match(line)
            if m:
                funfam = m.group(1)

    return signatures
