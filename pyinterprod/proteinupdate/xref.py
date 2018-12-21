import cx_Oracle

def delete_crc64_mismatches(url: str) -> int:
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.execute(
        """
        DELETE FROM INTERPRO.XREF_SUMMARY
        WHERE PROTEIN_AC IN (
            SELECT IP.PROTEIN_AC
            FROM UNIPARC.XREF UX
            INNER JOIN INTERPRO.PROTEIN IP
                ON UX.AC = IP.PROTEIN_AC
            INNER JOIN UNIPARC.PROTEIN UP
                ON UX.UPI = UP.UPI
            WHERE UX.DELETED = 'N'
            AND IP.CRC64 != UP.RC64
        )
        """
    )
    n = cur.rowcount
    con.commit()
    cur.close()
    con.close()
    return n
