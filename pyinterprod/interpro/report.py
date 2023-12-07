import os
import pickle
import re
import shutil
from datetime import datetime
from tempfile import mkdtemp
from zipfile import ZipFile, ZIP_DEFLATED

import oracledb

from pyinterprod.utils import email
from pyinterprod.pronto.signature import get_swissprot_descriptions
from .database import Database
from .match import track_entry_changes, get_sig_protein_counts
from .signature import FILE_DB_SIG, FILE_SIG_DESCR

MIN_ENTRY_CHANGE = 0.5
MIN_SIGNATURE_CHANGE = 0.1


def send_db_update_report(ora_url: str, pg_url: str, dbs: list[Database],
                          data_dir: str, pronto_link: str, emails: dict):
    # Get Swiss-Prot descriptions (after the update)
    all_sig2info = get_swissprot_descriptions(pg_url)

    pronto_link = pronto_link.rstrip('/')

    tmpdir = mkdtemp()
    id2dst = {}
    for db in dbs:
        name = db.name.replace(' ', '_')
        id2dst[db.identifier] = os.path.join(tmpdir, name)
        os.mkdir(id2dst[db.identifier])

    con = oracledb.connect(ora_url)
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

    cur.execute(
        """
        SELECT CODE, REPLACE(ABBREV, '_', ' ')
        FROM INTERPRO.CV_ENTRY_TYPE
        """
    )
    types = dict(cur.fetchall())

    with open(os.path.join(data_dir, FILE_DB_SIG), "rb") as fh:
        databases = pickle.load(fh)

    for db_id, data in databases.items():
        dst = id2dst[db_id]

        # Deleted signatures
        with open(os.path.join(dst, "deleted.tsv"), "wt") as fh:
            fh.write("Signature\tName\tDescription\tEntry\n")
            for acc, name, descr, entry_acc in data["deleted"]:
                if entry_acc:
                    # Only report signatures that were integrated
                    fh.write(f"{acc}\t{name or 'N/A'}\t{descr or 'N/A'}\t"
                             f"{entry_acc}\n")

        # Signatures with a different name
        with open(os.path.join(dst, "name_changes.tsv"), "wt") as fh:
            fh.write("Signature\tEntry\tLink\tPrevious name\tNew name\n")

            for acc, old_val, new_val in data["changes"]["names"]:
                try:
                    entry_acc = integrated[acc][0]
                except KeyError:
                    continue

                link = f"{pronto_link}/entry/{entry_acc}/"
                fh.write(f"{acc}\t{entry_acc}\t{link}\t"
                         f"{old_val or 'N/A'}\t{new_val or 'N/A'}\n")

        # Signatures with a different descriptions
        with open(os.path.join(dst, "description_changes.tsv"), "wt") as fh:
            fh.write("Signature\tEntry\tLink\tPrevious description"
                     "\tNew description\n")

            for acc, old_val, new_val in data["changes"]["descriptions"]:
                try:
                    entry_acc = integrated[acc][0]
                except KeyError:
                    continue

                link = f"{pronto_link}/entry/{entry_acc}/"
                fh.write(f"{acc}\t{entry_acc}\t{link}\t"
                         f"{old_val or 'N/A'}\t{new_val or 'N/A'}\n")

        # Signatures with a different type (if provided by member DB)
        with open(os.path.join(dst, "type_changes.tsv"), "wt") as fh:
            fh.write("Signature\tEntry\tLink\tPrevious type\tNew type\n")

            for acc, old_val, new_val in data["changes"]["types"]:
                try:
                    entry_acc = integrated[acc][0]
                except KeyError:
                    continue

                link = f"{pronto_link}/entry/{entry_acc}/"
                fh.write(f"{acc}\t{entry_acc}\t{link}\t{old_val}"
                         f"\t{new_val}\n")

        # Swiss-Prot descriptions
        old_sigs = data["descriptions"]
        new_sigs = {}
        descr2prots = {}
        for acc in all_sig2info:
            if acc in integrated and integrated[acc][3] == db_id:
                new_sigs[acc] = all_sig2info[acc]
            for descr, proteins in all_sig2info[acc].items():
                try:
                    descr2prots[descr] |= proteins
                except KeyError:
                    descr2prots[descr] = proteins

        changes = {}
        for acc, old_info in old_sigs.items():
            new_info = new_sigs.pop(acc, {})
            old_descrs = set(old_info.keys())
            new_descrs = set(new_info.keys())
            for descr, proteins in old_info.items():
                try:
                    descr2prots[descr] |= set(proteins)
                except KeyError:
                    descr2prots[descr] = set(proteins)
            try:
                entry_acc, entry_type, entry_name, _ = integrated[acc]
            except KeyError:
                continue

            if old_descrs != new_descrs:
                changes[acc] = (entry_acc, entry_name, entry_type,
                                old_descrs - new_descrs,
                                new_descrs - old_descrs)

        sig2prots = {}
        for acc, new_info in new_sigs.items():
            entry_acc, entry_type, entry_name, _ = integrated[acc]
            for descr, proteins in new_info.items():
                try:
                    sig2prots[acc] |= proteins
                except KeyError:
                    sig2prots[acc] = proteins

            new_descrs = list(new_info.keys())
            changes[acc] = (entry_acc, entry_name, entry_type, [], new_descrs)

        files = {}  # file objects
        for acc in sorted(changes):
            entry_acc, entry_name, entry_type, lost, gained = changes[acc]
            if entry_type == 'F':
                filename = "swiss_de_families.tsv"
            elif entry_type == 'D':
                filename = "swiss_de_domains.tsv"
            else:
                filename = "swiss_de_others.tsv"

            try:
                fh = files[filename]
            except KeyError:
                filepath = os.path.join(dst, filename)
                fh = files[filename] = open(filepath, "wt")
                fh.write(f"Signature\tLink\tEntry\tType\t# Swiss-Prot\tName\t# Lost"
                         f"\t# Gained\tLost\tGained\n")

            link = f"{pronto_link}/signatures/{acc}/descriptions/?reviewed"

            lost_descs = [
                f"{desc} ({list(descr2prots[desc])[0]})" for desc in sorted(lost)
            ]
            gained_descs = [
                f"{desc} ({list(descr2prots[desc])[0]})" for desc in sorted(gained)
            ]
            fh.write(f"{acc}\t{link}\t{entry_acc}\t{types[entry_type]}"
                     f"{len(sig2prots[acc])}\t"
                     f"\t{entry_name}\t{len(lost)}\t{len(gained)}"
                     f"\t{' | '.join(lost_descs)}"
                     f"\t{' | '.join(gained_descs)}\n")

        for fh in files.values():
            fh.close()

        # Keep track of signatures with Swiss-Prot description changes
        sig_changes = set(changes.keys())

        # Protein count changes (total + per superkingdom)
        old_counts = data["proteins"]
        new_counts = get_sig_protein_counts(cur, db_id)
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
            if sig_new_tot != 0 and abs(change) < MIN_SIGNATURE_CHANGE:
                continue

            sig_superkingdoms = {}
            for superkingdom, old_cnt in sig_old_cnts.items():
                new_cnt = sig_new_cnts.pop(superkingdom, 0)
                sig_superkingdoms[superkingdom] = (old_cnt, new_cnt)
                superkingdoms.add(superkingdom)

            # superkingdoms with proteins only matched in new version of DB
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
            for sk in superkingdoms:
                line += [sk, '']
            fh.write('\t'.join(line) + '\n')

            line = [''] * 9
            line += ["Previous count", "New count"] * len(superkingdoms)
            fh.write('\t'.join(line) + '\n')

            # Body
            for obj in changes:
                acc = obj[0]
                entry_acc = obj[1]
                entry_name = obj[2]
                entry_type = types[obj[3]]
                sig_old_tot = obj[4]
                sig_new_tot = obj[5]
                change = obj[6]
                sig_superkingdoms = obj[7]

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

        # New signatures
        with open(os.path.join(dst, "new.tsv"), "wt") as fh:
            fh.write("Signature\tName\tDescription\n")
            for acc, name, descr, _type in data["new"]:
                # Ignore PANTHER subfamilies (won't be integrated)
                if not re.fullmatch(r"PTHR\d+:SF\d+", acc):
                    fh.write(f"{acc}\t{name or 'N/A'}\t{descr or 'N/A'}\n")

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
        info=emails,
        to=["interpro"],
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

    con = oracledb.connect(ora_url)
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

    cur.execute(
        """
        SELECT CODE, REPLACE(ABBREV, '_', ' ')
        FROM INTERPRO.CV_ENTRY_TYPE
        """
    )
    types = dict(cur.fetchall())

    cur.execute("SELECT VERSION FROM INTERPRO.DB_VERSION WHERE DBCODE = 'u'")
    release, = cur.fetchone()

    descr2prots = {}
    # Load entry -> descriptions BEFORE UniProt update
    entries_then = {}
    with open(os.path.join(data_dir, FILE_SIG_DESCR), "rb") as fh:
        for signature_acc, old_info in pickle.load(fh).items():
            try:
                entry_acc = integrated[signature_acc]
            except KeyError:
                continue

            try:
                entry_descrs = entries_then[entry_acc]
            except KeyError:
                entry_descrs = entries_then[entry_acc] = set()

            for description, proteins in old_info.items():
                entry_descrs.add(description)
                try:
                    descr2prots[description] |= proteins
                except KeyError:
                    descr2prots[description] = proteins

    # Load entry -> descriptions AFTER UniProt update
    signatures_now = get_swissprot_descriptions(pg_url)
    entries_now = {}
    entry2prots = {}
    for signature_acc, info in signatures_now.items():
        try:
            entry_acc = integrated[signature_acc]
        except KeyError:
            continue

        for description, proteins in info.items():
            try:
                descr2prots[description] |= proteins
            except KeyError:
                descr2prots[description] = proteins

            try:
                entries_now[entry_acc].add(description)
            except KeyError:
                entries_now[entry_acc] = {description}
                entry2prots[entry_acc] = proteins
            else:
                entry2prots[entry_acc] |= proteins

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

    # Write entries with changes (by entry type: families, domains, others)
    tmpdir = mkdtemp()
    files = {}
    header = ("Accession\tLink\tName\tChecked\t# Swiss-Prot"
              "\t# Lost\t# Gained\tLost\tGained\n")
    for entry_acc in sorted(changes):
        gained, lost = changes[entry_acc]
        try:
            entry_type, name, checked_flag = entries[entry_acc]
        except KeyError:
            # Entry does not exist anymore
            continue

        if entry_type == 'F':
            filename = "swiss_de_families.tsv"
        elif entry_type == 'D':
            filename = "swiss_de_domains.tsv"
        else:
            filename = "swiss_de_others.tsv"

        try:
            fh = files[filename]
        except KeyError:
            filepath = os.path.join(tmpdir, filename)
            fh = files[filename] = open(filepath, "wt")
            fh.write(header)
        finally:
            lost_descs = [
                f"{desc} ({list(descr2prots[desc])[0]})" for desc in sorted(lost)
            ]
            gained_descs = [
                f"{desc} ({list(descr2prots[desc])[0]})" for desc in sorted(gained)
            ]
            fh.write(f"{entry_acc}\t{pronto_link}/entry/{entry_acc}/\t"
                     f"{name}\t{'Yes' if checked_flag == 'Y' else 'No'}\t"
                     f"{len(entry2prots[entry_acc])}\t"
                     f"{len(lost)}\t{len(gained)}\t"
                     f"{' | '.join(lost_descs)}\t"
                     f"{' | '.join(gained_descs)}\n")

    for fh in files.values():
        fh.close()

    # Keep track of entries with Swiss-Prot description changes
    entries_changes = set(changes.keys())

    # Write entries with protein count changes (total + per superkingdom)
    with open(os.path.join(tmpdir, "entries_count_changes.tsv"), "wt") as fh:
        changes = track_entry_changes(cur, data_dir, MIN_ENTRY_CHANGE)
        superkingdoms = sorted({sk for e in changes for sk in e[4]})

        # Header
        line = ["Accession", "Link", "Type", "Name", "Checked",
                "DE changes", "Previous count", "New count", "Change (%)"]
        for sk in superkingdoms:
            line += [sk, '']

        fh.write('\t'.join(line) + '\n')

        line = [''] * 9
        line += ["Previous count", "New count"] * len(superkingdoms)

        # Body
        for obj in changes:
            entry_acc = obj[0]
            old_total = obj[1]
            new_total = obj[2]
            change = obj[3]
            entry_superkingdoms = obj[4]
            try:
                entry_type, name, checked_flag = entries[entry_acc]
            except KeyError:
                # Entry does not exist anymore
                continue

            line = [
                entry_acc,
                f"{pronto_link}/entry/{entry_acc}/",
                types[entry_type],
                name,
                "Yes" if checked_flag == 'Y' else "No",
                "Yes" if entry_acc in entries_changes else "No",
                str(old_total),
                str(new_total),
                f"{change*100:.0f}"
            ]

            for sk in superkingdoms:
                try:
                    old_cnt, new_cnt = entry_superkingdoms[sk]
                except KeyError:
                    line += ['0', '0']
                else:
                    line += [str(old_cnt), str(new_cnt)]

            fh.write('\t'.join(line) + '\n')

    cur.close()
    con.close()

    filename = os.path.join(data_dir, f"protein_update_{release}.zip")
    with ZipFile(filename, 'w', compression=ZIP_DEFLATED) as fh:
        for root, dirs, files in os.walk(tmpdir):
            for file in files:
                path = os.path.join(root, file)
                fh.write(path, arcname=os.path.relpath(path, tmpdir))

    shutil.rmtree(tmpdir)

    email.send(
        info=emails,
        to=["interpro"],
        subject=f"Protein update report: UniProt {release}",
        content="""\
Dear curators,

Pronto has been refreshed. Please find attached a ZIP archive containing \
the report files for this protein update.

The InterPro Production Team
""",
        attachments=[filename]
    )
