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
    all_sig2swiss_new = get_swissprot_descriptions(pg_url)

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

        # Swiss-Prot descriptions before the update
        sig2swiss = {}
        for acc, proteins in data["descriptions"].items():
            sig2swiss[acc] = {}

            for protein_acc, description in proteins:
                # List: descr before, descr after
                sig2swiss[acc][protein_acc] = [description, None]

        # Add Swiss-Prot descriptions after the update
        for acc in all_sig2swiss_new:
            try:
                dbcode, _, _, _, _ = integrated[acc]
            except KeyError:
                continue

            if db_id == dbcode:
                proteins = all_sig2swiss_new[acc]
                try:
                    sig_swiss = sig2swiss[acc]
                except KeyError:
                    # First time we have matches for this signature
                    sig_swiss = sig2swiss[acc] = {}

                for protein_acc, description in proteins:
                    try:
                        sig_swiss[protein_acc][1] = description
                    except KeyError:
                        sig_swiss[protein_acc] = [None, description]

        files = {}  # file objects
        sig_changes = set()
        for acc in sorted(sig2swiss):
            try:
                _, entry_acc, type_code, entry_name, origin = integrated[acc]
            except KeyError:
                continue

            proteins = sig2swiss[acc]
            lost, gained = compare_descriptions(proteins)
            if not lost and not gained:
                continue
            elif type_code == "F":
                filename = "swiss_de_families.tsv"
            elif type_code == "D":
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

            cnt_protein = sum(1 for _, descr in proteins.values() if descr)
            link = f"{pronto_link}/signatures/{acc}/descriptions/?reviewed"
            fh.write(f"{acc}\t{link}\t{entry_acc}\t{types[type_code]}\t"
                     f"{origin}\t{cnt_protein}\t{entry_name}\t"
                     f"{len(lost)}\t{len(gained)}\t"
                     f"{' | '.join(lost)}\t"
                     f"{' | '.join(gained)}\n")

            # Keep track of signatures with Swiss-Prot description changes
            sig_changes.add(acc)

        for fh in files.values():
            fh.close()

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
                origin = obj[4]
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
            fh.write("Signature\tName\tDescription\t# Proteins\n")
            for acc, name, descr, _type in data["new"]:
                # Ignore PANTHER subfamilies (won't be integrated)
                if not re.fullmatch(r"PTHR\d+:SF\d+", acc):
                    sig_cnts = sum(new_counts.get(acc, {}).values())
                    fh.write(f"{acc}\t{name or 'N/A'}\t{descr or 'N/A'}\t"
                             f"{sig_cnts}\n")

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

    # Swiss-Prot descriptions before the update
    entry2swiss = {}
    with open(os.path.join(data_dir, FILE_SIG_DESCR), "rb") as fh:
        for signature_acc, proteins in pickle.load(fh).items():
            try:
                entry_acc = integrated[signature_acc]
            except KeyError:
                continue

            try:
                entry_swiss = entry2swiss[entry_acc]
            except KeyError:
                entry_swiss = entry2swiss[entry_acc] = {}

            for protein_acc, description in proteins:
                # List: descr before, descr after
                entry_swiss[protein_acc] = [description, None]

    # Add Swiss-Prot descriptions after the update
    for signature_acc, proteins in get_swissprot_descriptions(pg_url).items():
        try:
            entry_acc = integrated[signature_acc]
        except KeyError:
            continue

        try:
            entry_swiss = entry2swiss[entry_acc]
        except KeyError:
            # First time we have matches for this entry
            # (created during the updated?)
            entry_swiss = entry2swiss[entry_acc] = {}

        for protein_acc, description in proteins:
            try:
                entry_swiss[protein_acc][1] = description
            except KeyError:
                entry_swiss[protein_acc] = [None, description]

    # Write entries with changes (by entry type: families, domains, others)
    tmpdir = mkdtemp()
    files = {}
    entries_changes = set()
    header = ("Accession\tLink\tName\tChecked\tSource origin\t# Swiss-Prot"
              "\t# Lost\t# Gained\tLost\tGained\n")
    for entry_acc in sorted(entry2swiss):
        try:
            type_code, name, is_checked, origin = entries[entry_acc]
        except KeyError:
            # Entry does not exist anymore
            # (deleted while this function is running)
            continue

        proteins = entry2swiss[entry_acc]
        lost, gained = compare_descriptions(proteins)
        if not lost and not gained:
            continue
        elif type_code == "F":
            filename = "swiss_de_families.tsv"
        elif type_code == "D":
            filename = "swiss_de_domains.tsv"
        else:
            filename = "swiss_de_others.tsv"

        try:
            fh = files[filename]
        except KeyError:
            filepath = os.path.join(tmpdir, filename)
            fh = files[filename] = open(filepath, "wt")
            fh.write(header)

        cnt_protein = sum(1 for _, descr in proteins.values() if descr)
        fh.write(f"{entry_acc}\t{pronto_link}/entry/{entry_acc}/\t"
                 f"{name}\t{'Yes' if is_checked else 'No'}\t{origin}\t"
                 f"{cnt_protein}\t"
                 f"{len(lost)}\t{len(gained)}\t"
                 f"{' | '.join(lost)}\t"
                 f"{' | '.join(gained)}\n")

        # Write a "slim" version, with renamed proteins ignored
        slim_filename = f"slim_{filename}"

        try:
            fh = files[slim_filename]
        except KeyError:
            filepath = os.path.join(tmpdir, slim_filename)
            fh = files[slim_filename] = open(filepath, "wt")
            fh.write(header)

        lost, gained = compare_descriptions(proteins, ignore_renamed=True)
        fh.write(f"{entry_acc}\t{pronto_link}/entry/{entry_acc}/\t"
                 f"{name}\t{'Yes' if is_checked else 'No'}\t{origin}\t"
                 f"{cnt_protein}\t"
                 f"{len(lost)}\t{len(gained)}\t"
                 f"{' | '.join(lost)}\t"
                 f"{' | '.join(gained)}\n")

        entries_changes.add(entry_acc)

    for fh in files.values():
        fh.close()

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


def compare_descriptions(
        proteins: dict[str, list[str | None]],
        ignore_renamed: bool = False
) -> tuple[list[str], list[str]]:
    descrs_old = {}
    descrs_new = {}
    # Never ignore proteins whose name contained or contains these strings:
    never_ignore = [
        "unknown",
        "putative",
        "uncharacterized",
        "hypothetical",
        "fragment",
        "predicted"
    ]

    for protein_acc in sorted(proteins):
        descr_old, descr_new = proteins[protein_acc]

        if (ignore_renamed and
                descr_old and
                descr_new and
                descr_old != descr_new):
            for s in never_ignore:
                if s in descr_old.lower() or s in descr_new.lower():
                    break
            else:
                # We can safely ignore protein being renamed
                continue

        if descr_old and descr_old not in descrs_old:
            descrs_old[descr_old] = protein_acc

        if descr_new and descr_new not in descrs_new:
            descrs_new[descr_new] = protein_acc

    lost = []
    for descr in sorted(set(descrs_old.keys()) - set(descrs_new.keys())):
        protein_acc = descrs_old[descr]
        lost.append(f"{descr} ({protein_acc})")

    gained = []
    for descr in sorted(set(descrs_new.keys()) - set(descrs_old.keys())):
        protein_acc = descrs_new[descr]
        gained.append(f"{descr} ({protein_acc})")

    return lost, gained
