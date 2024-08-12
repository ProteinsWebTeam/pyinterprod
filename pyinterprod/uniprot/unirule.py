import oracledb

from pyinterprod.utils import email, oracle


def report_integration_changes(uri: str, emails: dict):
    """Sends a list of integration changes.
    If a signature is integrated in an unchecked entry, we consider
    the signature as unintegrated.

    :param uri: Oracle connection string
    :param emails: email info (SMTP server/port, sender, recipients, etc.)
    """
    con = oracledb.connect(uri)
    cur = con.cursor()
    cur.execute(
        """
        SELECT DISTINCT METHOD_AC, ENTRY_AC
        FROM INTERPRO.XREF_SUMMARY
        WHERE ENTRY_AC IS NOT NULL 
        """
    )
    previous_signatures = dict(cur.fetchall())
    cur.execute(
        """
        SELECT M.METHOD_AC, E.ENTRY_AC
        FROM INTERPRO.METHOD M
        LEFT OUTER JOIN INTERPRO.ENTRY2METHOD EM
          ON M.METHOD_AC = EM.METHOD_AC
        LEFT OUTER JOIN INTERPRO.ENTRY E
          ON EM.ENTRY_AC = E.ENTRY_AC
          AND E.CHECKED = 'Y'
        """
    )
    current_signatures = dict(cur.fetchall())
    cur.execute("SELECT VERSION FROM INTERPRO.DB_VERSION WHERE DBCODE = 'u'")
    release = cur.fetchone()[0]
    cur.close()
    con.close()

    changes = {}
    for signature, entry_then in previous_signatures.items():
        try:
            entry_now = current_signatures.pop(signature)
        except KeyError:
            # Signature does not exist any more (member database update)
            changes[signature] = ("deleted", entry_then, '')
            continue

        if entry_now is None:
            # Signature unintegrated
            changes[signature] = ("unintegrated", entry_then, '')
        elif entry_now != entry_then:
            # Integrated in a different entry
            changes[signature] = ("integrated", entry_then, entry_now)

    # All new integrations (signature was not integrated before)
    for signature, entry_now in current_signatures.items():
        if entry_now:
            changes[signature] = ("integrated", '', entry_now)

    content = """\
Dear UniProt team,

Please find below the list of recent integration changes.

Description of states:
  - integrated:       the signature has been integrated in an InterPro entry
  - unintegrated:     the signature has been removed from an InterPro entry
  - deleted:          the signature does not exist anymore

"""
    content += (f"{'Signature':<30}{'Status':<30}{'Previous entry':^30}"
                f"{'Current entry':^30}\n")
    content += '-' * 120 + '\n'

    for acc in sorted(changes, key=lambda k: k.lower()):
        status, prev, curr = changes[acc]
        content += f"{acc:<30}{status:<30}{prev:^30}{curr:^30}\n"

    content += "\nKind regards,\nThe InterPro Production Team\n"

    email.send(
        info=emails,
        to=["aa_dev"],
        cc=["unirule"],
        bcc=["sender"],
        subject=f"Protein update {release}: integration changes",
        content=content
    )


def ask_to_snapshot(uri: str, emails: dict):
    con = oracledb.connect(uri)
    cur = con.cursor()
    cur.execute("SELECT VERSION FROM INTERPRO.DB_VERSION WHERE DBCODE = 'u'")
    release, = cur.fetchone()
    cur.execute("SELECT MAX(UPI) FROM UNIPARC.PROTEIN")
    max_upi, = cur.fetchone()
    cur.execute(
        """
        SELECT I.IPRSCAN_SIG_LIB_REL_ID, D.DBNAME, V.VERSION
        FROM INTERPRO.IPRSCAN2DBCODE I
        INNER JOIN INTERPRO.CV_DATABASE D ON I.DBCODE = D.DBCODE
        LEFT OUTER JOIN INTERPRO.DB_VERSION V ON I.DBCODE = V.DBCODE
        WHERE I.IPRSCAN_SIG_LIB_REL_ID != I.PREV_ID
        ORDER BY D.DBNAME
        """
    )
    updates = []
    for analysis_id, name, version in cur:
        if version:
            name += " (" + version + ")"
        updates.append((name, analysis_id))

    content = f"""\
All tables for UniProt are ready.

Please, take a snapshot of IPPRO, and inform UniProt they can refresh IPREADU.

Recipients
----------
To:  {emails['uniprot_db']}
Cc:  {emails['aa']}, {emails['uniprot_prod']}, {emails['interpro']}

Subject
-------
Protein update {release} completed: please refresh IPREADU

Body (replace [DATE-TIME] by snapshot date/time)
------------------------------------------------
Dear UniProt team,

You may refresh IPREADU with the snapshot of IPPRO from [DATE-TIME]

The max UPI we processed up to in this update is {max_upi}.
"""

    if updates:
        content += """
We have updated the following databases since the previous protein update:
"""

        for name, analysis_id in updates:
            content += f"  - {name:<30}new analysis ID: {analysis_id}\n"
    else:
        content += """
There are no changes to the analysis IDs for this protein update.
"""

    content += """
Kind regards,
The InterPro Production Team
"""

    try:
        email.send(
            info=emails,
            to=["interpro"],
            subject=f"Protein update {release}: please snapshot IPPRO",
            content=content
        )

        # Update PREV_ID only if mail sent
        cur.execute(
            """
            UPDATE INTERPRO.IPRSCAN2DBCODE
            SET PREV_ID = IPRSCAN_SIG_LIB_REL_ID
            """
        )
        con.commit()
    finally:
        cur.close()
        con.close()


def update_signatures(filepath: str, uri: str):
    # Load accessions from file generated by UniRule
    accessions = []
    with open(filepath, "rt") as fh:
        for line in fh:
            database, accession = line.rstrip().split()
            if database == "GDAC":
                accessions.append(f"G3DSA:{accession}")
            elif database != "IPRO":
                # ignore InterPro entries
                accessions.append(accession)

    con = oracledb.connect(uri)
    cur = con.cursor()

    # Refresh data
    oracle.truncate_table(cur, "INTERPRO.METHOD_UNIRULE")

    cur.executemany(
        """
        INSERT INTO INTERPRO.METHOD_UNIRULE VALUES (:1)
        """, [(acc,) for acc in set(accessions)]
    )

    con.commit()
    cur.close()
    con.close()
