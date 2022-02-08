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
    Version 4.3.0:
        - http://download.cathdb.info/cath/releases/latest-release/sequence-data/funfam-hmm3.lib.gz
        - ftp://orengoftp.biochem.ucl.ac.uk//cath/releases/latest-release/sequence-data/funfam-hmm3.lib.gz

    :param filepath:
    :return:
    """

    signatures = []
    with open(filepath, "rt") as fh:
        supfam = funfam = None
        for line in fh:
            if line[:2] == "//":
                accession = f"{_PREFIX}{supfam}:FF:{funfam}"
                signatures.append(Method(accession, _TYPE_FUNFAM))
                supfam = funfam = None
                continue

            m = re.search(r"^NAME\s+(\d+\.\d+\.\d+\.\d+)-FF-(\d+)$", line)
            if m:
                supfam, funfam = m.groups()

    return signatures
