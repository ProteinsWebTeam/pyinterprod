# -*- coding: utf-8 -*-

import os
import pickle
import shutil
from datetime import datetime
from tempfile import mkdtemp
from typing import Sequence
from zipfile import ZipFile, ZIP_DEFLATED

import cx_Oracle

from pyinterprod.utils import email
from pyinterprod.pronto.signature import get_swissprot_descriptions
from .database import Database
from .match import track_entry_changes, track_sig_changes
from .signature import FILE_DB_SIG_DIFF, FILE_SIG_DESCR


def send_db_update_report(ora_url: str, pg_url: str, dbs: Sequence[Database],
                          data_dir: str, pronto_link: str, emails: dict):
    pronto_link = pronto_link.rstrip('/')

    tmpdir = mkdtemp()
    id2dst = {}
    for db in dbs:
        name = db.name.replace(' ', '_')
        id2dst[db.identifier] = os.path.join(tmpdir, name)
        os.mkdir(id2dst[db.identifier])

    con = cx_Oracle.connect(ora_url)
    cur = con.cursor()
    cur.execute(
        """
        SELECT EM.METHOD_AC, E.ENTRY_AC, E.ENTRY_TYPE, E.NAME
        FROM INTERPRO.ENTRY2METHOD EM
        INNER JOIN INTERPRO.ENTRY E ON EM.ENTRY_AC = E.ENTRY_AC
        """
    )
    integrated = {}
    for sig_acc, entry_acc, entry_type, entry_name in cur:
        integrated[sig_acc] = (entry_acc, entry_type, entry_name)

    cur.execute("SELECT METHOD_AC, DBCODE FROM INTERPRO.METHOD")
    acc2dbid = dict(cur.fetchall())

    # Protein counts changes
    for db_id, db_changes in track_sig_changes(cur, dbs, data_dir).items():
        dst = id2dst[db_id]

        with open(os.path.join(dst, "protein_counts.tsv"), "wt") as fh:
            fh.write("Signature\tLink\tEntry\tPrevious count"
                     "\tNew count\tChange (%)\n")

            for sig_acc, old_cnt, new_cnt, change in db_changes:
                try:
                    entry_acc = integrated[sig_acc]
                except KeyError:
                    continue

                link = f"{pronto_link}/signature/{sig_acc}/"
                fh.write(f"{sig_acc}\t{link}\t{entry_acc}\t{old_cnt}"
                         f"\t{new_cnt}\t{change * 100:.0f}\n")

    cur.close()
    con.close()

    # Annotation changes
    with open(os.path.join(data_dir, FILE_DB_SIG_DIFF), "rb") as fh:
        for db_id, db_changes in pickle.load(fh).items():
            dst = id2dst[db_id]

            with open(os.path.join(dst, "name_changes.tsv"), "wt") as fh2:
                fh2.write("Signature\tEntry\tLink\tPrevious name\tNew name\n")

                changes = db_changes["changes"]["names"]
                for sig_acc, old_val, new_val in sorted(changes):
                    try:
                        entry_acc = integrated[sig_acc][0]
                    except KeyError:
                        continue

                    link = f"{pronto_link}/entry/{entry_acc}/"
                    fh2.write(f"{sig_acc}\t{entry_acc}\t{link}\t{old_val}"
                              f"\t{new_val}\n")

            with open(os.path.join(dst, "descr_changes.tsv"), "wt") as fh2:
                fh2.write("Signature\tEntry\tLink\tPrevious description"
                          "\tNew description\n")

                changes = db_changes["changes"]["descriptions"]
                for sig_acc, old_val, new_val in sorted(changes):
                    try:
                        entry_acc = integrated[sig_acc][0]
                    except KeyError:
                        continue

                    link = f"{pronto_link}/entry/{entry_acc}/"
                    fh2.write(f"{sig_acc}\t{entry_acc}\t{link}\t{old_val}"
                              f"\t{new_val}\n")

            with open(os.path.join(dst, "type_changes.tsv"), "wt") as fh2:
                fh2.write("Signature\tEntry\tLink\tPrevious type\tNew type\n")

                changes = db_changes["changes"]["types"]
                for sig_acc, old_val, new_val in sorted(changes):
                    try:
                        entry_acc = integrated[sig_acc][0]
                    except KeyError:
                        continue

                    link = f"{pronto_link}/entry/{entry_acc}/"
                    fh2.write(f"{sig_acc}\t{entry_acc}\t{link}\t{old_val}"
                              f"\t{new_val}\n")

    # Swiss-Prot DE
    new_sigs = get_swissprot_descriptions(pg_url)
    changes = []
    with open(os.path.join(data_dir, FILE_SIG_DESCR), "rb") as fh:
        for sig_acc, old_descrs in pickle.load(fh).items():
            new_descrs = new_sigs.pop(sig_acc, set())

            try:
                entry_acc, entry_type, entry_name = integrated[sig_acc]
                db_id = acc2dbid[sig_acc]
            except KeyError:
                # signature is unintegrated, or deleted
                continue

            changes.append((sig_acc, db_id, entry_acc, entry_name,
                            entry_type, old_descrs - new_descrs,
                            new_descrs - old_descrs))

    for sig_acc, new_descrs in new_sigs:
        try:
            entry_acc, entry_type, entry_name = integrated[sig_acc]
            db_id = acc2dbid[sig_acc]
        except KeyError:
            # signature is unintegrated, or deleted
            continue

        changes.append((sig_acc, db_id, entry_acc, entry_name,
                        entry_type, [], new_descrs))

    files_descr = {}
    changes.sort(key=lambda x: x[0])
    for sig_acc, db_id, e_acc, e_name, e_type, lost, gained in changes:
        try:
            dst = id2dst[db_id]
        except KeyError:
            continue  # signature not from updated member database

        if e_type == 'F':
            label = "families"
        elif e_type == 'D':
            label = "domains"
        else:
            label = "others"

        file = os.path.join(dst, f"swiss_de_{label}.tsv")
        try:
            fh = files_descr[file]
        except KeyError:
            fh = files_descr[file] = open(file, "wt")
            fh.write(f"Signature\tLink\tEntry\tName\tType\t# Lost"
                     f"\t# Gained\tLost\tGained\n")

        link = f"{pronto_link}/signatures/{sig_acc}/descriptions/?reviewed"
        fh.write(f"{sig_acc}\t{link}\t{e_acc}\t{e_name}\t{e_type}\t{len(lost)}"
                 f"\t{len(gained)}\t{' | '.join(sorted(lost))}"
                 f"\t{' | '.join(sorted(gained))}\n")

    for fh in files_descr.values():
        fh.close()

    date = datetime.today().strftime("%Y%m%d")
    filename = os.path.join(data_dir, f"member_database_update_{date}.zip")
    with ZipFile(filename, 'w', compression=ZIP_DEFLATED) as fh:
        for root, dirs, files in os.walk(tmpdir):
            for file in files:
                path = os.path.join(root, file)
                fh.write(path, arcname=os.path.relpath(path, tmpdir))

    shutil.rmtree(tmpdir)


def send_prot_update_report(ora_url: str, pg_url: str, data_dir: str,
                            pronto_link: str, emails: dict):
    pronto_link = pronto_link.rstrip('/')

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
    with open(os.path.join(data_dir, FILE_SIG_DESCR), "rb") as fh:
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

    # Write entries with changes (two different files: families/others)
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
                fh.write(f"{entry_acc}\t{pronto_link}/entry/{entry_acc}/\t"
                         f"{name}\t{type_code}\t{checked_flag}\t{len(lost)}\t"
                         f"{len(gained)}\t{' | '.join(lost)}\t"
                         f"{' | '.join(gained)}\n")

    # Write entries with protein count changes
    path = os.path.join(data_dir, "entries_count_changes.tsv")
    with open(path, "wt") as ofh:
        ofh.write("# Accession\tLink\tName\tType\tChecked\t"
                  "Previous count\tNew count\tChange (%)\n")

        changes = track_entry_changes(cur, data_dir)
        for entry_acc, prev_count, count, change in changes:
            name, type_code, checked_flag = entries[entry_acc]
            ofh.write(f"{entry_acc}\t{pronto_link}/entry/{entry_acc}/\t{name}\t"
                      f"{type_code}\t{checked_flag}\t"
                      f"{prev_count}\t{count}\t{change*100:.0f}\n")

    filename = os.path.join(data_dir, f"protein_update_{release}.zip")
    with ZipFile(filename, 'w', compression=ZIP_DEFLATED) as fh:
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
the report files for this protein update.

The InterPro Production Team
""",
        attachments=[filename]
    )
