import cx_Oracle


def import_gp(cur: cx_Oracle.Cursor, file_path: str):
    gp = []
    with open(file_path, "rt") as fh:
        for line in fh:
            if line.startswith('AC'):
                entry_ac = line.split()[1]
            elif line.startswith('DE'):
                description = line.split(" ", maxsplit=1)[1].strip()
            elif line.startswith('EV'):
                ev = line.split()[1].split(';')[0]
                if ev.startswith('IPR'):
                    gp.append((ev, 'h', description, entry_ac))

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
        """, gp
    )


if __name__ == '__main__':
    con = cx_Oracle.connect('connection string')
    cur = con.cursor()
    import_gp(cur, 'path/to/genomeProperties.txt')
    cur.close()
    con.close()
