from dataclasses import dataclass, field
from datetime import datetime
import xml.etree.ElementTree as ET
from typing import List, Optional


@dataclass
class Method:
    accession: str
    sig_type: Optional[str]
    name: Optional[str] = None
    description: Optional[str] = None
    abstract: Optional[str] = None
    date: Optional[datetime] = None


@dataclass
class Clan:
    accession: str
    name: str = None
    description: str = None
    members: list = field(default_factory=list)


def parse_xml(filepath: str, sig_type: str) -> List[Method]:
    """
    Parse the interpro XML file provided by member databases

    :param filepath: path the XML file
    :param sig_type: signature type
    :return:
    """

    with open(filepath, "rt", errors="replace") as fh:
        tree = ET.parse(fh)

    root = tree.getroot()
    namespace = "{http://www.ebi.ac.uk/schema/interpro}"
    signatures = []
    for sig in root.findall(f"{namespace}signature"):
        attrib = sig.attrib

        try:
            abstract = sig.find(f"{namespace}abstract").text.strip()
        except AttributeError:
            abstract = None

        signatures.append(Method(accession=attrib["ac"],
                                 sig_type=sig_type,
                                 name=attrib["name"].replace(',', '_'),
                                 description=attrib["desc"],
                                 abstract=abstract))

    return signatures
