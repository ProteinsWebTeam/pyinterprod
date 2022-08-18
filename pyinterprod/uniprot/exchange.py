import os

import cx_Oracle

from pyinterprod import logger
from pyinterprod.utils import email


def export_sib(url: str, emails: dict):
    logger.info("exporting data")
    con = cx_Oracle.connect(url)
    cur = con.cursor()

    """
    The SIB_EXPORT package belongs to INTERPRODP and exports the tables:
    * INTERPRO.DB_VERSION
    * INTERPRO.ENTRY
    * INTERPRO.ENTRY2ENTRY
    * INTERPRO.ENTRY2METHOD
    * INTERPRO.MATCH
    * INTERPRO.METHOD
    * INTERPRO.PROTEIN
    * INTERPRO.XREF_CONDENSED
    """
    cur.callproc("SIB_EXPORT.EXP_SIB_DATA")

    logger.info("data successfully exported")

    cur.execute("SELECT VERSION FROM INTERPRO.DB_VERSION WHERE DBCODE = 'u'")
    version, = cur.fetchone()
    cur.close()
    con.close()

    email.send(
        info=emails,
        to=["interpro"],
        subject=f"Protein update {version}: data for SIB ready",
        content=f"""\
The data required by SIB was successfully exported.

Please, archive the dump on the FTP, and inform SIB that the archive is ready.

Recipients
----------
To: {emails['sib']}

Subject
-------
InterPro data for UniProt private release available

Body 
----
Dear Swiss-Prot team,

The interpro.tar.gz archive for UniProt private release {version} \
is available at ftp://ftp-private-beta.ebi.ac.uk/interpro/

Kind regards,
The InterPro Production Team
"""
    )


def export_xrefs(url: str, outdir: str, emails: dict):
    """
    Format for Uniprot dat files:
      <protein>    DR   <database>; <signature/entry>; <name>; <count>.
    Exceptions:
        - Gene3D: do not include prefix before accession (G3DSA:)
                  replace signature's name by hyphen (-)
        - PRINTS: do not include match count
        - InterPro: do not include match count
    """
    logger.info("exporting dat files")
    os.makedirs(outdir, 0o775, exist_ok=True)

    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.execute("SELECT VERSION FROM INTERPRO.DB_VERSION WHERE DBCODE = 'u'")
    release = cur.fetchone()[0]
    cur.execute(
        """
        SELECT
          PROTEIN_AC,
          METHOD_AC,
          MIN(METHOD_NAME),
          MIN(DBCODE),
          MIN(ENTRY_AC),
          MIN(SHORT_NAME),
          COUNT(*)
        FROM INTERPRO.XREF_SUMMARY
        GROUP BY PROTEIN_AC, METHOD_AC
        ORDER BY PROTEIN_AC
        """
    )

    dbcodes = {
        'B': ("SFLD", "SFLD"),
        # 'D': ("ProDom", "PD"),  # ProDom removed from InterPro
        'F': ("PRINTS", "PP"),
        'H': ("Pfam", "PF"),
        'J': ("CDD", "CDD"),
        'M': ("PROSITE", "PR"),
        'N': ("TIGRFAMs", "TF"),
        'P': ("PROSITE", "PR"),
        'Q': ("HAMAP", "HP"),
        'R': ("SMART", "SM"),
        'U': ("PIRSF", "PI"),
        'V': ("PANTHER", "PTHR"),
        'X': ("Gene3D", "G3D"),
        'Y': ("SUPFAM", "SF"),
    }

    files = []
    handlers = {}
    for dbcode in dbcodes:
        if dbcode != 'P':
            # P (PROSITE patterns) -> same file than M (PROSITE profiles)
            dbname, dbkey = dbcodes[dbcode]

            filepath = os.path.join(outdir, dbname + ".dat")
            files.append(filepath)
            handlers[dbcode] = open(filepath, "wt")

    handlers['P'] = handlers['M']

    filepath = os.path.join(outdir, "InterPro.dat")
    files.append(filepath)
    ifh = open(filepath, "wt")

    entries = {}
    prev_acc = None
    for row in cur:
        protein_acc = row[0]
        signature_acc = row[1]
        signature_name = row[2]
        dbcode = row[3]
        entry_acc = row[4]
        entry_name = row[5]
        num_matches = int(row[6])

        dbname, dbkey = dbcodes[dbcode]
        fh = handlers[dbcode]

        if dbcode == 'X':
            """
            Gene3D
              - accession: G3DSA:3.90.1580.10 -> 3.90.1580.10
              - do not print signature's name
            """
            fh.write(f"{protein_acc}    DR   {dbname}; {signature_acc[6:]}; "
                     f"-; {num_matches}.\n")
        elif dbcode == 'F':
            # PRINTS: do not print match count
            fh.write(f"{protein_acc}    DR   {dbname}; {signature_acc}; "
                     f"{signature_name}.\n")
        elif dbcode == "V" and "NOT NAMED" in signature_name:
            fh.write(f"{protein_acc}    DR   {dbname}; {signature_acc}; "
                     f"-; {num_matches}.\n")
        else:
            fh.write(f"{protein_acc}    DR   {dbname}; {signature_acc}; "
                     f"{signature_name}; {num_matches}.\n")

        if protein_acc != prev_acc:
            for entry in sorted(entries):
                name = entries[entry]
                ifh.write(f"{prev_acc}    DR   InterPro; {entry}; {name}.\n")

            entries = {}
            prev_acc = protein_acc

        entries[entry_acc] = entry_name

    # Last protein
    if prev_acc:
        for entry in sorted(entries):
            name = entries[entry]
            ifh.write(f"{prev_acc}    DR   InterPro; {entry}; {name}.\n")

    ifh.close()

    for fh in handlers.values():
        fh.close()

    for path in files:
        os.chmod(path, 0o775)

    email.send(
        info=emails,
        to=["uniprot_db"],
        cc=["uniprot_prod"],
        bcc=["sender"],
        subject="InterPro XREF files are ready",
        content=f"""\
Dear UniProt team,

The InterPro cross-references files for {release} are available \
in the following directory:
  {outdir}
  
Kind regards,
The InterPro Production Team
"""
        )

    logger.info("complete")
