# -*- coding: utf-8 -*-

from typing import List
import xml.etree.ElementTree as ET

from .model import Method

_TYPE = 'D'


def parse_signatures(filepath: str) -> List[Method]:
    """
    Parse the cdd_interpro.xml file provided by CDD
    As of CDD 3.18, the file is available here:
        ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/cdd_interpro.xml.gz

    :param filepath:
    :return:
    """
    signatures = []

    tree = ET.parse(filepath)
    root = tree.getroot()
    namespace = "{http://www.ebi.ac.uk/schema/interpro}"
    for sig in root.findall(f"{namespace}signature"):
        attrib = sig.attrib
        abstract = sig.find(f"{namespace}abstract").text.strip()
        signatures.append(Method(accession=attrib["ac"],
                                 sig_type=_TYPE,
                                 name=attrib["name"].replace(',', '_'),
                                 description=attrib["desc"],
                                 abstract=abstract))

    return signatures
