import os

import oracledb

from pyinterprod import logger
from pyinterprod.utils import email


def export(url: str, outdir: str, emails: dict):
    """
    Format for Uniprot dat files:
      <protein>    DR   <database>; <signature/entry>; <name>; <count>.

    (https://web.expasy.org/docs/userman.html#DR_line)

    Exceptions:
        - Gene3D: do not include prefix before accession (G3DSA:)
        - PRINTS: do not include match count
        - InterPro: do not include match count
        - PANTHER: add a record for the subfamily
    """
    logger.info("exporting dat files")
    os.makedirs(outdir, 0o775, exist_ok=True)

    con = oracledb.connect(url)
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
        "a": "AntiFam",
        "B": "SFLD",
        "F": "PRINTS",
        "f": "FunFam",
        "H": "Pfam",
        "J": "CDD",
        "M": "PROSITE",
        "N": "NCBIfam",
        "P": "PROSITE",
        "Q": "HAMAP",
        "R": "SMART",
        "U": "PIRSF",
        "V": "PANTHER",
        "X": "Gene3D",
        "Y": "SUPFAM",
    }

    files = []
    handlers = {}
    for dbcode in dbcodes:
        if dbcode != 'P':
            # P (PROSITE patterns) -> same file than M (PROSITE profiles)
            dbname = dbcodes[dbcode]

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

        dbname = dbcodes[dbcode]
        fh = handlers[dbcode]

        if dbcode in ("X", "f"):
            # CATH-Gene3D/FunFam: G3DSA:3.50.70.10 -> 3.50.70.10
            identifier = signature_acc[6:]
        else:
            identifier = signature_acc

        optional_1 = ((signature_name or "-")
                      .replace("\"", "'")
                      .replace(";", ","))

        if dbcode in ("a", "F"):
            # AntiFam and PRINTS: no number of matches (who knows why)
            optional_2 = ""
        else:
            optional_2 = f"; {num_matches}"

        fh.write(f"{protein_acc}    DR   {dbname}; {identifier}; "
                 f"{optional_1}{optional_2}.\n")

        if protein_acc != prev_acc:
            for entry in sorted(entries):
                name = entries[entry]
                ifh.write(f"{prev_acc}    DR   InterPro; {entry}; {name}.\n")

            entries.clear()
            prev_acc = protein_acc

        if entry_acc:
            entries[entry_acc] = entry_name

    # Last protein
    for entry in sorted(entries):
        name = entries[entry]
        ifh.write(f"{prev_acc}    DR   InterPro; {entry}; {name}.\n")

    ifh.close()

    for fh in handlers.values():
        fh.close()

    for path in files:
        os.chmod(path, 0o664)

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
