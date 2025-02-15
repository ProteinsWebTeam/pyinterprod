import gzip
import re

from .common import Clan, Method, parse_xml

_TYPE = 'D'


def parse_signatures(filepath: str) -> list[Method]:
    """
    Parse the cdd_interpro.xml file provided by CDD
    As of CDD 3.20, the file is available here:
        https://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/cdd_interpro.xml.gz

    :param filepath:
    :return:
    """
    signatures = parse_xml(filepath, _TYPE)

    for s in signatures:
        if s.description is not None and s.description.endswith("..."):
            s.description = None

    return signatures


def get_clans(cddid: str, fam2supfam: str) -> list[Clan]:
    """
    Return CDD superfamilies (clans)
    :param cddid: path to file containing summary information about
        the CD models. Available on CDD's FTP:
        ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/cddid.tbl.gz
    :param fam2supfam: File containing the list of NCBI-curated and
        imported domain models that are members of CDD superfamilies.
        Available on CDD's FTP:
        ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/family_superfamily_links
    :return: a list of clans
    """
    superfamilies = {}
    families = set()
    if cddid.lower().endswith(".gz"):
        fh = gzip.open(cddid, "rt")
    else:
        fh = open(cddid, "rt")

    for line in fh:
        fields = line.rstrip().split("\t")
        accession = fields[1]
        name = fields[2]
        descr = fields[3].lstrip("N/A. ")

        if re.match(r"cl\d+", accession):
            superfamilies[accession] = Clan(accession, name, descr)
        elif re.match(r"[cs]d\d+", accession):
            families.add(accession)

    fh.close()

    with open(fam2supfam, "rt") as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue

            fields = line.split("\t")
            family_acc = fields[0]
            supfam_acc = fields[2]

            if supfam_acc in superfamilies and family_acc in families:
                superfamilies[supfam_acc].members.append({
                    "accession": family_acc,
                    "score": 1
                })

    return list(superfamilies.values())
