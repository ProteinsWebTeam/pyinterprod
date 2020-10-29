# -*- coding: utf-8 -*-

import re
from typing import List

from .model import Method

_TYPE = 'D'


def parse_signatures(filepath: str) -> List[Method]:
    """
    Parse the cddid.tbl file distributed with CDD releases

    :param filepath:
    :return:
    """
    signatures = []
    with open(filepath, "rt") as fh:
        for line in fh:
            values = line.rstrip().split('\t')
            pssm_id, acc, name, descr, pssm_length = values
            if re.match(r"cd\d+$", acc):
                signatures.append(Method(accession=acc,
                                         sig_type=_TYPE,
                                         name=name.replace(',', '_'),
                                         description=name,
                                         abstract=descr))

    return signatures
