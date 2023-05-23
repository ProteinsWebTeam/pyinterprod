import cx_Oracle


def import_gp(cur: cx_Oracle.Cursor, file_path: str):
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
                if m:
                    gp.append((m.group(0), "h", description, entry_ac)(

    cur.execute(
        """
        DELETE FROM INTERPRO.ENTRY_XREF
        WHERE DBCODE = 'h'
        """
    )

    cur.executemany(
        """
        INSERT INTO INTERPRO.ENTRY_XREF (
        ENTRY_AC, DBCODE, AC, NAME
        ) VALUES (:1, :2, :3, :4)
        """,
        gp,
    )


if __name__ == "__main__":
    con = cx_Oracle.connect("connection string")
    cur = con.cursor()
    import_gp(cur, "path/to/genomeProperties.txt")
    cur.close()
    con.commit()
    con.close()
