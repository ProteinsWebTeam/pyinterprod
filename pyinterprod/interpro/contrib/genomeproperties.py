import os
import re
import sys

import oracledb

_DBCODE = "h"


def update_xrefs(uri: str, file_path: str):
    con = oracledb.connect(uri)
    cur = con.cursor()

    cur.execute("SELECT ENTRY_AC FROM INTERPRO.ENTRY")
    entries = {acc for acc, in cur.fetchall()}

    gp = set()
    with open(file_path, "rt") as fh:
        for line in map(str.rstrip, fh):
            if line.startswith("AC"):
                entry_ac = line.split()[1]
            elif line.startswith("DE"):
                description = line.split(maxsplit=1)[1]
            elif line.startswith("EV"):
                _, evidence = line.split(maxsplit=1)
                m = re.match(r"IPR\d+", evidence)
                if m and m.group(0) in entries:
                    gp.add((m.group(0), _DBCODE, description, entry_ac))

    cur.execute(
        """
        DELETE FROM INTERPRO.ENTRY_XREF
        WHERE DBCODE = :1
        """,
        [_DBCODE]
    )

    cur.executemany(
        """
        INSERT INTO INTERPRO.ENTRY_XREF (
        ENTRY_AC, DBCODE, NAME, AC
        ) VALUES (:1, :2, :3, :4)
        """,
        list(gp),
    )

    cur.close()
    con.commit()
    con.close()


if __name__ == "__main__":
    # Latest version: https://github.com/ebi-pf-team/genome-properties/blob/master/flatfiles/genomeProperties.txt
    gp_file = sys.argv[1]
    update_xrefs(os.environ["INTERPRO_CONNECTION_URL"], gp_file)
