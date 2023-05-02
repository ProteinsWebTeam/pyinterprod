import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from datetime import datetime
from typing import Optional


@dataclass
class Method:
    accession: str
    sig_type: Optional[str]
    name: Optional[str] = None
    description: Optional[str] = None
    abstract: Optional[str] = None
    date: Optional[datetime] = None
    references: list[int] = field(default_factory=list)
    model: Optional[str] = None


@dataclass
class Clan:
    accession: str
    name: str = None
    description: str = None
    members: list = field(default_factory=list)


def parse_hmm(filepath: str):
    reg_name = re.compile(r"^NAME\s+(.+)$", flags=re.M)
    reg_acc = re.compile(r"^ACC\s+(.+)$", flags=re.M)
    reg_desc = re.compile(r"^DESC\s+(.+)$", flags=re.M)
    reg_date = re.compile(r"^DATE\s+(.+)$", flags=re.M)

    with open(filepath, "rt") as fh:
        buffer = ""
        for line in fh:
            buffer += line

            if line[:2] == "//":
                # Mandatory field
                name = reg_name.search(buffer).group(1)
                acc = descr = dt = None

                try:
                    # Optional field
                    acc = reg_acc.search(buffer).group(1)
                except AttributeError:
                    pass

                try:
                    # Optional
                    descr = reg_desc.search(buffer).group(1)
                except AttributeError:
                    pass

                try:
                    date_string = reg_date.search(buffer).group(1)
                except AttributeError:
                    pass
                else:
                    # Example: Wed Sep  1 14:33:46 2021
                    parts = date_string.split()
                    if len(parts[2]) == 1:
                        parts[2] = f"0{parts[2]}"

                    date_string = " ".join(parts)
                    dt = datetime.strptime(date_string, "%a %b %d %H:%M:%S %Y")

                yield acc, name, descr, dt
                buffer = ""


def parse_xml(filepath: str, sig_type: str) -> list[Method]:
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
