# -*- coding: utf-8 -*-

import os
import pickle
from typing import Dict, Set
from zipfile import ZipFile, ZIP_DEFLATED

import cx_Oracle
import psycopg2
from pyinterprod.utils.pg import url2dict

from .match import ENTRIES_TSV
from pyinterprod.utils import email


_DESCRIPTIONS_DAT = "swissprot_descr.dat"


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

    with open(os.path.join(data_dir, _DESCRIPTIONS_DAT), "wb") as fh:
        pickle.dump(signatures, fh)


def send_report(ora_url: str, pg_url: str, data_dir: str, pronto_url: str,
                emails: dict):
    pronto_url = pronto_url.rstrip('/')

    con = cx_Oracle.connect(ora_url)
    cur = con.cursor()
    cur.execute("SELECT METHOD_AC, ENTRY_AC FROM INTERPRO.ENTRY2METHOD")
    integrated = dict(cur.fetchall())

    cur.execute(
        """
        SELECT ENTRY_AC, ENTRY_TYPE, NAME, CHECKED 
        FROM INTERPRO.ENTRY
        """
    )
    entries = {row[0]: row[1:] for row in cur}

    cur.execute("SELECT VERSION FROM INTERPRO.DB_VERSION WHERE DBCODE = 'u'")
    release, = cur.fetchone()
    cur.close()
    con.close()

    # Load entry -> descriptions BEFORE UniProt update
    entries_then = {}
    with open(os.path.join(data_dir, _DESCRIPTIONS_DAT), "rb") as fh:
        for signature_acc, descriptions in pickle.load(fh).items():
            try:
                entry_acc = integrated[signature_acc]
            except KeyError:
                continue

            try:
                entries_then[entry_acc] |= descriptions
            except KeyError:
                entries_then[entry_acc] = set(descriptions)

    # Load entry -> descriptions AFTER UniProt update
    signatures_now = get_swissprot_descriptions(pg_url)
    entries_now = {}
    for signature_acc, descriptions in signatures_now.items():
        try:
            entry_acc = integrated[signature_acc]
        except KeyError:
            continue

        try:
            entries_now[entry_acc] |= descriptions
        except KeyError:
            entries_now[entry_acc] = set(descriptions)

    changes = {}  # key: entry accession, value: (gained, lost)
    for entry_acc, descs_now in entries_now.items():
        try:
            descs_then = entries_then.pop(entry_acc)
        except KeyError:
            # This entry did not have descriptions before the update
            changes[entry_acc] = (descs_now, [])
        else:
            changes[entry_acc] = (
                descs_now - descs_then,
                descs_then - descs_now
            )

    # Entries that no longer have descriptions
    for entry_acc, descs_then in entries_then.items():
        changes[entry_acc] = ([], descs_then)

    # Write entries with changes to two files (family/others)
    files = {}
    header = ("Accession\tLink\tName\tType\tChecked\t# Lost\t# Gained\t"
              "Lost\tGained\n")
    for entry_acc in sorted(changes):
        gained, lost = changes[entry_acc]
        if gained or lost:
            name, type_code, checked_flag = entries[entry_acc]
            entry_type = "families" if type_code == 'F' else "others"

            try:
                fh, path = files[entry_type]
            except KeyError:
                path = os.path.join(data_dir, f"swiss_de_{entry_type}.tsv")
                fh = open(path, "wt")
                fh.write(header)
                files[entry_type] = (fh, path)
            finally:
                fh.write(f"{entry_acc}\t{pronto_url}/entry/{entry_acc}/\t"
                         f"{name}\t{type_code}\t{checked_flag}\t{len(lost)}\t"
                         f"{len(gained)}\t{' | '.join(lost)}\t"
                         f"{' | '.join(gained)}\n")

    filename = os.path.join(data_dir, f"protein_update_{release}.zip")
    with ZipFile(filename, 'w', compression=ZIP_DEFLATED) as fh:
        path = os.path.join(data_dir, ENTRIES_TSV)
        fh.write(path, arcname=os.path.basename(path))

        for ofh, path in files.values():
            ofh.close()
            fh.write(path, arcname=os.path.basename(path))

    email.send(
        emails,
        subject=f"Protein update report: UniProt {release}",
        content="""\
Dear curators,

Pronto has been refreshed. Please find attached a ZIP archive containing \
the report files for this month's protein update.

The InterPro Production Team
""",
        attachments=[filename]
    )
