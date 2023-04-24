from datetime import datetime

from .common import Method


def parse_instances(db_sources: dict) -> list[Method]:
    """
    Parse the ELM instances.tsv file, available at http://elm.eu.org/

    :param db_sources:
    :return:
    """
    instances = []
    with open(db_sources["sig_source"], "rt") as fh:
        date = None
        for line in fh:
            if line[0] == '#':
                if line[1:28] == 'ELM_Instance_Download_Date:':
                    date = line[29:-1]
                    date = datetime.strptime(date, "%Y-%m-%d %H:%M:%S.%f")
                else:
                    continue
            else:
                cols = line.split("\t")
                if cols[10].strip('"') == 'true positive':
                    instances.append(Method(cols[0].strip('"'), None, cols[9].strip('"'), None, None, date))

    return instances
