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
from .match import track_entry_changes, get_sig_proteins_count
from .signature import FILE_DB_SIG, FILE_SIG_DESCR


def send_db_update_report(ora_url: str, pg_url: str, dbs: Sequence[Database],
                          data_dir: str, pronto_link: str, emails: dict):
    # Get Swiss-Prot descriptions (after the update)
    all_sig2descs = get_swissprot_descriptions(pg_url)

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
        SELECT M.METHOD_AC, E.ENTRY_AC, E.ENTRY_TYPE, E.NAME, M.DBCODE
        FROM INTERPRO.METHOD M
        INNER JOIN INTERPRO.ENTRY2METHOD EM ON M.METHOD_AC = EM.METHOD_AC
        INNER JOIN INTERPRO.ENTRY E ON EM.ENTRY_AC = E.ENTRY_AC
        """
    )
    integrated = {row[0]: row[1:] for row in cur}

    with open(os.path.join(data_dir, FILE_DB_SIG), "rb") as fh:
        databases = pickle.load(fh)

    for db_id, data in databases.items():
        dst = id2dst[db_id]

        # Deleted signatures
        with open(os.path.join(dst, "deleted.tsv"), "wt") as fh:
            fh.write("Signature\tName\tDescription\tEntry\n")
            for acc, name, descr, entry_acc in sorted(data["deleted"]):
                if entry_acc:
                    # Only report signatures that were integrated
                    fh.write(f"{acc}\t{name}\t{descr}\t{entry_acc}\n")

        # Annotation changes
        with open(os.path.join(dst, "name_changes.tsv"), "wt") as fh:
            fh.write("Signature\tEntry\tLink\tPrevious name\tNew name\n")

            changes = data["changes"]["names"]
            for acc, old_val, new_val in sorted(changes):
                try:
                    entry_acc = integrated[acc][0]
                except KeyError:
                    continue

                link = f"{pronto_link}/entry/{entry_acc}/"
                fh.write(f"{acc}\t{entry_acc}\t{link}\t"
                         f"{old_val or 'N/A'}\t{new_val or 'N/A'}\n")

        with open(os.path.join(dst, "description_changes.tsv"), "wt") as fh:
            fh.write("Signature\tEntry\tLink\tPrevious description"
                     "\tNew description\n")

            changes = data["changes"]["descriptions"]
            for acc, old_val, new_val in sorted(changes):
                try:
                    entry_acc = integrated[acc][0]
                except KeyError:
                    continue

                link = f"{pronto_link}/entry/{entry_acc}/"
                fh.write(f"{acc}\t{entry_acc}\t{link}\t"
                         f"{old_val or 'N/A'}\t{new_val or 'N/A'}\n")

        with open(os.path.join(dst, "type_changes.tsv"), "wt") as fh:
            fh.write("Signature\tEntry\tLink\tPrevious type\tNew type\n")

            changes = data["changes"]["types"]
            for acc, old_val, new_val in sorted(changes):
                try:
                    entry_acc = integrated[acc][0]
                except KeyError:
                    continue

                link = f"{pronto_link}/entry/{entry_acc}/"
                fh.write(f"{acc}\t{entry_acc}\t{link}\t{old_val}"
                         f"\t{new_val}\n")

        # Swiss-Prot descriptions
        old_sigs = data["descriptions"]
        new_sigs = {
            acc: all_sig2descs[acc]
            for acc in all_sig2descs
            if acc in integrated and integrated[acc][3] == db_id
        }
        changes = {}
        for acc, old_descrs in old_sigs.items():
            new_descrs = new_sigs.pop(acc, set())

            try:
                entry_acc, entry_type, entry_name, _ = integrated[acc]
            except KeyError:
                continue

            if old_descrs != new_descrs:
                changes[acc] = (entry_acc, entry_name, entry_type,
                                old_descrs - new_descrs,
                                new_descrs - old_descrs)

        for acc, new_descrs in new_sigs.items():
            entry_acc, entry_type, entry_name, _ = integrated[acc]
            changes[acc] = (entry_acc, entry_name, entry_type, [], new_descrs)

        files = {}  # file objects
        for acc in sorted(changes):
            entry_acc, entry_name, entry_type, lost, gained = changes[acc]
            if entry_type == 'F':
                entry_type = "families"
            elif entry_type == 'D':
                entry_type = "domains"
            else:
                entry_type = "others"

            file = os.path.join(dst, f"swiss_de_{entry_type}.tsv")
            try:
                fh = files[file]
            except KeyError:
                fh = files[file] = open(file, "wt")
                fh.write(f"Signature\tLink\tEntry\tName\t# Lost"
                         f"\t# Gained\tLost\tGained\n")

            link = f"{pronto_link}/signatures/{acc}/descriptions/?reviewed"
            fh.write(f"{acc}\t{link}\t{entry_acc}\t{entry_name}"
                     f"\t{len(lost)}\t{len(gained)}"
                     f"\t{' | '.join(sorted(lost))}"
                     f"\t{' | '.join(sorted(gained))}\n")

        for fh in files.values():
            fh.close()

        # Keep track of signatures with Swiss-Prot description changes
        sig_changes = set(changes.keys())

        # Protein count changes (total + per superkingdom)
        old_counts = data["proteins"]
        new_counts = get_sig_proteins_count(cur, db_id)
        changes = []
        superkingdoms = set()
        for acc in sorted(old_counts):  # sort by accession
            sig_old_cnts = old_counts[acc]
            sig_new_cnts = new_counts.get(acc, {})

            try:
                entry_acc, entry_type, entry_name, _ = integrated[acc]
            except KeyError:
                continue

            # Total number of proteins matched
            sig_old_tot = sum(sig_old_cnts.values())
            sig_new_tot = sum(sig_new_cnts.values())

            change = (sig_new_tot - sig_old_tot) / sig_old_tot

            # If the signature does not have any matches anymore,
            # we want to report it (it is integrated in InterPro)
            if sig_new_tot != 0 and abs(change) < 0.1:
                continue

            sig_superkingdoms = {}
            for superkingdom, old_cnt in sig_old_cnts.items():
                new_cnt = sig_new_cnts.pop(superkingdom, 0)
                sig_superkingdoms[superkingdom] = (old_cnt, new_cnt)
                superkingdoms.add(superkingdom)

            # superkingdoms with proteins only matches in new version of DB
            for superkingdom, new_cnt in sig_new_cnts.items():
                sig_superkingdoms[superkingdom] = (0, new_cnt)
                superkingdoms.add(superkingdom)

            changes.append((
                acc,
                entry_acc,
                entry_name,
                entry_type,
                sig_old_tot,
                sig_new_tot,
                change,
                sig_superkingdoms
            ))

        superkingdoms = sorted(superkingdoms)
        with open(os.path.join(dst, "protein_counts.tsv"), "wt") as fh:
            # Header
            line = ["Signature", "Link", "DE changes", "Entry", "Type",
                    "Name", "Previous count", "New count", "Change (%)"]
            line += superkingdoms
            fh.write('\t'.join(line) + '\n')

            line = [''] * 9
            line += ["Previous count", "New count"] * len(superkingdoms)
            fh.write('\t'.join(line) + '\n')

            # Body
            for obj in changes:
                acc = obj[0]
                entry_acc = obj[1]
                entry_name = obj[2]
                entry_type = obj[3]
                sig_old_tot = obj[4]
                sig_new_tot = obj[5]
                change = obj[6]
                sig_superkingdoms = obj[7]

                if entry_type == 'F':
                    entry_type = "Family"
                elif entry_type == 'D':
                    entry_type = "Domain"
                else:
                    entry_type = "Other"

                line = [
                    acc,
                    f"{pronto_link}/signature/{acc}/",
                    "Yes" if acc in sig_changes else "No",
                    entry_acc,
                    entry_type,
                    entry_name,
                    str(sig_old_tot),
                    str(sig_new_tot),
                    f"{change * 100:.0f}"
                ]
                for superkingdom in superkingdoms:
                    try:
                        old_cnt, new_cnt = sig_superkingdoms[superkingdom]
                    except KeyError:
                        line += ['0', '0']
                    else:
                        line += [str(old_cnt), str(new_cnt)]

                fh.write('\t'.join(line) + '\n')

    cur.close()
    con.close()

    date = datetime.today().strftime("%Y%m%d")
    filename = os.path.join(data_dir, f"member_database_update_{date}.zip")
    with ZipFile(filename, 'w', compression=ZIP_DEFLATED) as fh:
        for root, dirs, files in os.walk(tmpdir):
            for file in files:
                path = os.path.join(root, file)
                fh.write(path, arcname=os.path.relpath(path, tmpdir))

    shutil.rmtree(tmpdir)

    names = [f"{db.name} {db.version}" for db in dbs]
    email.send(
        emails,
        subject=f"Member database update report: {', '.join(names)}",
        content="""\
Dear curators,

Pronto has been refreshed. Please find attached a ZIP archive containing \
the report files for this member database update.

The InterPro Production Team
""",
        attachments=[filename]
    )


def send_prot_update_report(ora_url: str, pg_url: str, data_dir: str,
                            pronto_link: str, emails: dict):
    pronto_link = pronto_link.rstrip('/')

    con = cx_Oracle.connect(ora_url)
    cur = con.cursor()
    cur.execute("SELECT METHOD_AC, ENTRY_AC FROM INTERPRO.ENTRY2METHOD")
    integrated = dict(cur.fetchall())

    cur.execute(
        """
        SELECT ENTRY_AC, NAME, ENTRY_TYPE, CHECKED 
        FROM INTERPRO.ENTRY
        """
    )
    entries = {row[0]: row[1:] for row in cur}

    cur.execute("SELECT VERSION FROM INTERPRO.DB_VERSION WHERE DBCODE = 'u'")
    release, = cur.fetchone()

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
        descs_then = entries_then.pop(entry_acc, set())

        if descs_now != descs_then:
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
        type_code, name, checked_flag = entries[entry_acc]
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
                     f"{len(gained)}\t{' | '.join(sorted(lost))}\t"
                     f"{' | '.join(sorted(gained))}\n")

    # Write entries with protein count changes
    path = os.path.join(data_dir, "entries_count_changes.tsv")
    with open(path, "wt") as ofh:
        ofh.write("# Accession\tLink\tName\tType\tChecked\t"
                  "Previous count\tNew count\tChange (%)\n")

        changes = track_entry_changes(cur, data_dir)
        for entry_acc, prev_count, count, change in changes:
            name, type_code, checked_flag = entries[entry_acc]
            ofh.write(f"{entry_acc}\t{pronto_link}/entry/{entry_acc}/\t"
                      f"{name}\t{type_code}\t{checked_flag}\t"
                      f"{prev_count}\t{count}\t{change*100:.0f}\n")

    cur.close()
    con.close()

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
