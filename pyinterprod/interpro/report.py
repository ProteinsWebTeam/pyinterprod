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
    all_sig2descrs = get_swissprot_descriptions(pg_url)

    pronto_link = pronto_link.rstrip('/')

    tmpdir = mkdtemp()
    id2dst = {}
    for db in dbs:
        name = db.name.replace(' ', '_')
        id2dst[db.identifier] = os.path.join(tmpdir, name)
        os.mkdir(id2dst[db.identifier])

    con = oracledb.connect(ora_url)
    cur = con.cursor()
    integrated = {}
    cur.execute(
        """
        SELECT M.METHOD_AC, M.DBCODE, E.ENTRY_AC, E.ENTRY_TYPE, E.NAME, 
               E.LLM, E.LLM_CHECKED
        FROM INTERPRO.METHOD M
        INNER JOIN INTERPRO.ENTRY2METHOD EM ON M.METHOD_AC = EM.METHOD_AC
        INNER JOIN INTERPRO.ENTRY E ON EM.ENTRY_AC = E.ENTRY_AC
        """
    )
    for row in cur.fetchall():
        if row[5] == "N":
            origin = "Curator"
        elif row[6] == "Y":
            origin = "AI, reviewed"
        else:
            origin = "AI, unreviewed"

        # value: dbcode, entry_acc, type_code, entry_name, origin
        integrated[row[0]] = (row[1], row[2], row[3], row[4], origin)

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
                    _, entry_acc, _, _, _ = integrated[acc]
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
                    _, entry_acc, _, _, _ = integrated[acc]
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
                    _, entry_acc, _, _, _ = integrated[acc]
                except KeyError:
                    continue

                link = f"{pronto_link}/entry/{entry_acc}/"
                fh.write(f"{acc}\t{entry_acc}\t{link}\t{old_val}"
                         f"\t{new_val}\n")

        # Swiss-Prot descriptions
        sigs_then = data["descriptions"]
        sigs_now = {}
        sig2prots = {}
        for acc in all_sig2descrs:
            try:
                dbcode, _, _, _, _ = integrated[acc]
            except KeyError:
                continue

            if db_id == dbcode:
                sigs_now[acc] = all_sig2descrs[acc]
                sig2prots[acc] = set()

                for proteins in sigs_now[acc].values():
                    sig2prots[acc] |= proteins

        changes = {}
        for acc in (set(sigs_then.keys()) | set(sigs_now.keys())):
            try:
                _, entry_acc, type_code, entry_name, origin = integrated[acc]
            except KeyError:
                continue

            descrs_then = set(sigs_then.get(acc, {}).keys())
            descrs_now = set(sigs_now.get(acc, {}).keys())
            if descrs_then != descrs_now:
                changes[entry_acc] = (
                    entry_acc, entry_name, type_code, origin,
                    descrs_then - descrs_now,
                    descrs_now - descrs_then
                )

        files = {}  # file objects
        for acc, obj in sorted(changes.items(), key=lambda x: x[0]):
            entry_acc, entry_name, type_code, origin, lost, gained = obj
            if type_code == 'F':
                filename = "swiss_de_families.tsv"
            elif type_code == 'D':
                filename = "swiss_de_domains.tsv"
            else:
                filename = "swiss_de_others.tsv"

            try:
                fh = files[filename]
            except KeyError:
                filepath = os.path.join(dst, filename)
                fh = files[filename] = open(filepath, "wt")
                fh.write(f"Signature\tLink\tEntry\tType\t"
                         f"Source origin\t# Swiss-Prot\tName\t"
                         f"# Lost\t# Gained\tLost\tGained\n")

            link = f"{pronto_link}/signatures/{acc}/descriptions/?reviewed"

            lost_descrs = []
            for descr in sorted(lost):
                ex_acc = next(iter(sigs_then[acc][descr]))
                lost_descrs.append(f"{descr} ({ex_acc})")

            gained_descrs = []
            for descr in sorted(gained):
                ex_acc = next(iter(sigs_now[acc][descr]))
                gained_descrs.append(f"{descr} ({ex_acc})")

            fh.write(f"{acc}\t{link}\t{entry_acc}\t{types[type_code]}\t"
                     f"{origin}\t{len(sig2prots.get(acc, []))}\t{entry_name}\t"
                     f"{len(lost)}\t{len(gained)}\t"
                     f"{' | '.join(lost_descrs)}\t"
                     f"{' | '.join(gained_descrs)}\n")

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
                _, entry_acc, type_code, entry_name, origin = integrated[acc]
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
                type_code,
                origin,
                sig_old_tot,
                sig_new_tot,
                change,
                sig_superkingdoms
            ))

        superkingdoms = sorted(superkingdoms)
        with open(os.path.join(dst, "protein_counts.tsv"), "wt") as fh:
            # Header
            line = ["Signature", "Link", "DE changes", "Entry", "Type",
                    "Source origin", "Name", "Previous count", "New count",
                    "Change (%)"]
            for sk in superkingdoms:
                line += [sk, '']
            fh.write('\t'.join(line) + '\n')

            line = [''] * 10
            line += ["Previous count", "New count"] * len(superkingdoms)
            fh.write('\t'.join(line) + '\n')

            # Body
            for obj in changes:
                acc = obj[0]
                entry_acc = obj[1]
                entry_name = obj[2]
                entry_type = types[obj[3]]
                origin = obj[4],
                sig_old_tot = obj[5]
                sig_new_tot = obj[6]
                change = obj[7]
                sig_superkingdoms = obj[8]

                line = [
                    acc,
                    f"{pronto_link}/signature/{acc}/",
                    "Yes" if acc in sig_changes else "No",
                    entry_acc,
                    entry_type,
                    origin,
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

    entries = {}
    cur.execute(
        """
        SELECT ENTRY_AC, ENTRY_TYPE, NAME, CHECKED, LLM, LLM_CHECKED
        FROM INTERPRO.ENTRY
        """
    )
    for entry_acc, type_code, name, checked, llm, llm_checked in cur.fetchall():
        if llm == "N":
            origin = "Curator"
        elif llm_checked == "Y":
            origin = "AI, reviewed"
        else:
            origin = "AI, unreviewed"

        entries[entry_acc] = (type_code, name, checked == "Y", origin)

    cur.execute(
        """
        SELECT CODE, REPLACE(ABBREV, '_', ' ')
        FROM INTERPRO.CV_ENTRY_TYPE
        """
    )
    types = dict(cur.fetchall())

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
                entry_descrs = entries_then[entry_acc]
            except KeyError:
                entry_descrs = entries_then[entry_acc] = {}

            for description, proteins in descriptions.items():
                try:
                    entry_descrs[description] |= proteins
                except KeyError:
                    entry_descrs[description] = set(proteins)

    # Load entry -> descriptions AFTER UniProt update
    signatures_now = get_swissprot_descriptions(pg_url)
    entries_now = {}
    entry2prots = {}
    for signature_acc, descriptions in signatures_now.items():
        try:
            entry_acc = integrated[signature_acc]
        except KeyError:
            continue

        try:
            entry_descrs = entries_now[entry_acc]
        except KeyError:
            entry_descrs = entries_now[entry_acc] = {}
            entry2prots[entry_acc] = set()

        for description, proteins in descriptions.items():
            entry2prots[entry_acc] |= proteins

            try:
                entry_descrs[description] |= proteins
            except KeyError:
                entry_descrs[description] = set(proteins)

    changes = {}  # key: entry accession, value: (lost, gained)
    for entry_acc in (set(entries_then.keys()) | set(entries_now.keys())):
        descrs_then = set(entries_then.get(entry_acc, {}).keys())
        descrs_now = set(entries_now.get(entry_acc, {}).keys())
        if descrs_then != descrs_now:
            changes[entry_acc] = (
                descrs_then - descrs_now,
                descrs_now - descrs_then
            )

    # Write entries with changes (by entry type: families, domains, others)
    tmpdir = mkdtemp()
    files = {}
    header = ("Accession\tLink\tName\tChecked\tSource origin\t# Swiss-Prot"
              "\t# Lost\t# Gained\tLost\tGained\n")
    for entry_acc in sorted(changes):
        try:
            type_code, name, is_checked, origin = entries[entry_acc]
        except KeyError:
            # Entry does not exist anymore
            continue

        if type_code == 'F':
            filename = "swiss_de_families.tsv"
        elif type_code == 'D':
            filename = "swiss_de_domains.tsv"
        else:
            filename = "swiss_de_others.tsv"

        try:
            fh = files[filename]
        except KeyError:
            filepath = os.path.join(tmpdir, filename)
            fh = files[filename] = open(filepath, "wt")
            fh.write(header)

        lost, gained = changes[entry_acc]
        lost_descrs = []
        for descr in sorted(lost):
            ex_acc = next(iter(entries_then[entry_acc][descr]))
            lost_descrs.append(f"{descr} ({ex_acc})")

        gained_descrs = []
        for descr in sorted(gained):
            ex_acc = next(iter(entries_now[entry_acc][descr]))
            gained_descrs.append(f"{descr} ({ex_acc})")

        fh.write(f"{entry_acc}\t{pronto_link}/entry/{entry_acc}/\t"
                 f"{name}\t{'Yes' if is_checked else 'No'}\t{origin}\t"
                 f"{len(entry2prots[entry_acc])}\t"
                 f"{len(lost)}\t{len(gained)}\t"
                 f"{' | '.join(lost_descrs)}\t"
                 f"{' | '.join(gained_descrs)}\n")

    for fh in files.values():
        fh.close()

    # Keep track of entries with Swiss-Prot description changes
    entries_changes = set(changes.keys())

    # Write entries with protein count changes (total + per superkingdom)
    with open(os.path.join(tmpdir, "entries_count_changes.tsv"), "wt") as fh:
        changes = track_entry_changes(cur, data_dir, MIN_ENTRY_CHANGE)
        superkingdoms = sorted({sk for e in changes for sk in e[4]})

        # Header
        line = ["Accession", "Link", "Type", "Name", "Checked", "Source origin",
                "DE changes", "Previous count", "New count", "Change (%)"]
        for sk in superkingdoms:
            line += [sk, '']

        fh.write('\t'.join(line) + '\n')

        line = [''] * 10
        line += ["Previous count", "New count"] * len(superkingdoms)
        fh.write('\t'.join(line) + '\n')

        # Body
        for obj in changes:
            entry_acc = obj[0]
            old_total = obj[1]
            new_total = obj[2]
            change = obj[3]
            entry_superkingdoms = obj[4]
            try:
                type_code, name, is_checked, origin = entries[entry_acc]
            except KeyError:
                # Entry does not exist anymore
                continue

            line = [
                entry_acc,
                f"{pronto_link}/entry/{entry_acc}/",
                types[type_code],
                name,
                "Yes" if is_checked else "No",
                origin,
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
