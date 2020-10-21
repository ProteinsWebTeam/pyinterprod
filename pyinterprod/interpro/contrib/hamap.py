# -*- coding: utf-8 -*-

import re
from datetime import datetime
from typing import List

from .model import Method

_TYPE = 'F'


def parse_signatures(filepath: str) -> List[Method]:
    """
    Parse the hamap.prf file distributed with HAMAP releases

    :param filepath:
    :return:
    """
    signatures = []
    reg_id = re.compile(r"^ID\s+(.+?);", flags=re.M)
    reg_acc = re.compile(r"^AC\s+(.+?);", flags=re.M)
    reg_dt = re.compile(r"^DT\s+(\d\d-[A-Z]-\d{4} CREATED)", flags=re.M)
    reg_de = re.compile(r"^DE\s+(.+?).$", flags=re.M)

    with open(filepath, "rt") as fh:
        profile = ""
        for line in fh:
            profile += line

            if line[:2] == "//":
                try:
                    date_string = reg_dt.search(profile).group(1)
                    date = datetime.strptime(date_string, "%d-%b-%Y")
                except (AttributeError, ValueError):
                    date = None

                signatures.append(Method(
                    accession=reg_acc.search(profile).group(1),
                    sig_type=_TYPE,
                    name=reg_id.search(profile).group(1),
                    description=reg_de.search(profile).group(1),
                    date=date
                ))

                profile = ""

    return signatures
