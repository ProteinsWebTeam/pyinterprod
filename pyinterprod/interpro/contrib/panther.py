# -*- coding: utf-8 -*-

import re
from typing import List

from .model import Method

_TYPE = 'F'


def parse_signatures(filepath: str) -> List[Method]:
    """
    Parse the names.tab file distributed with PANTHER releases

    :param filepath:
    :return:
    """
    signatures = []
    prog = re.compile(r"(PTHR\d+)\.(SF\d+)?")
    with open(filepath, "rt") as fh:
        for line in fh:
            accession, name = line.rstrip().split('\t')
            fam_id, sub_id = prog.search(accession).groups()
            if sub_id:
                accession = f"{fam_id}:{sub_id}"
            else:
                accession = fam_id

            if name.upper().endswith("NOT NAMED"):
                descr = None
            else:
                descr = name

            m = Method(accession, _TYPE, description=descr)
            signatures.append(m)

    return signatures
