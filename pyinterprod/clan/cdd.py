# -*- coding: utf-8 -*-

import gzip
import os
import re
from typing import List

from .utils import download


def get_clans() -> List[dict]:
    """
    Summary information about the CD models in this distribution that are
    part of the CD-Search tool's default "cdd" database 
    and are indexed in NCBI's Entrez CDD database 
    """
    filepath = download("ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/cddid.tbl.gz")
    superfamilies = {}
    families = {}
    with gzip.open(filepath, "rt") as fh:
        for line in fh:
            fields = line.rstrip().split("\t")
            accession = fields[1]
            name = fields[2]
            descr = fields[3].lstrip("N/A. ")

            if re.match("cl\d+", accession):
                superfamilies[accession] = {
                    "accession": accession,
                    "name": name,
                    "description": descr,
                    "members": []
                }
            elif re.match("cd\d+", accession):
                families[accession] = {
                    "accession": accession,
                    "name": name,
                    "description": descr
                }

    os.remove(filepath)

    """
    List of NCBI-curated and imported domain models that are members 
    of CDD superfamilies, along with the superfamily accession (cl*) 
    to which each domain model belongs
    """
    filepath = download("ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/family_superfamily_links")
    with open(filepath, "rt") as fh:
        for line in fh:
            fields = line.rstrip().split("\t")
            family_acc = fields[0]
            supfam_acc = fields[2]

            if supfam_acc in superfamilies and family_acc in families:
                superfamilies[supfam_acc]["members"].append({
                    "accession": family_acc,
                    "score": 1
                })

    os.remove(filepath)

    return list(superfamilies.values())
