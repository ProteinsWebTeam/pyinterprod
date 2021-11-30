import re
from typing import List

from .common import Method


def parse_models(filepath: str) -> List[Method]:
    """
    Parse the AntiFam HMM file,
    available at https://ftp.ebi.ac.uk/pub/databases/Pfam/AntiFam/current/Antifam.tar.gz

    :param filepath:
    :return:
    """

    models = []
    reg_name = re.compile(r"^NAME\s+(.+)$")
    reg_acc = re.compile(r"^ACC\s+(.+)$")
    reg_desc = re.compile(r"^DESC\s+(.+)$")
    with open(filepath, "rt") as fh:
        name = acc = desc = None
        for line in fh:
            if line[:2] == "//":
                models.append(Method(acc, None, name, desc, None, None))
                continue

            m = reg_name.match(line)
            if m:
                name = m.group(1)
                continue

            m = reg_acc.match(line)
            if m:
                acc = m.group(1)
                continue

            m = reg_desc.match(line)
            if m:
                desc = m.group(1)

    return models
