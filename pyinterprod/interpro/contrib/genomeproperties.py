import os
import re
import sys

import cx_Oracle

_DBCODE = "h"


def update_xrefs(uri: str, file_path: str):
    con = cx_Oracle.connect(uri)
    cur = con.cursor()

    cur.execute("SELECT ENTRY_AC FROM INTERPRO.ENTRY")
    ip_entry_ac = [row[0] for row in cur.fetchall()]

    gp = []
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
                    gp.append((m.group(0), _DBCODE, description, entry_ac))

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
        ENTRY_AC, DBCODE, AC, NAME
        ) VALUES (:1, :2, :3, :4)
        """,
        gp,
    )

    cur.close()
    con.commit()
    con.close()


if __name__ == "__main__":
    gp_file = sys.argv[1]
    update_xrefs(os.environ["INTERPRO_CONNECTION_URL"], gp_file)
