# -*- coding: utf-8 -*-

import os
import pickle
from typing import Dict, Set

import psycopg2
from pyinterprod.utils.pg import url2dict

SIGNATURES_DESCR_FILE = "swissprot_descr.dat"


def get_swissprot_descriptions(url: str) -> Dict[str, Set[str]]:
    con = psycopg2.connect(**url2dict(url))
    cur = con.cursor("spdecur")
    cur.execute(
        """
        SELECT DISTINCT s2p.signature_acc, pn.text
        FROM interpro.signature2protein s2p
        INNER JOIN interpro.protein_name pn ON s2p.name_id = pn.name_id
        WHERE s2p.is_reviewed            
        """
    )
    signatures = {}
    for acc, text in cur:
        try:
            signatures[acc].add(text)
        except KeyError:
            signatures[acc] = {text}

    cur.close()
    con.close()
    return signatures


def export_swissprot_description(url: str, data_dir: str):
    signatures = get_swissprot_descriptions(url)

    with open(os.path.join(data_dir, SIGNATURES_DESCR_FILE), "wb") as fh:
        pickle.dump(signatures, fh)
