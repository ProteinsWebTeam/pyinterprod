import cx_Oracle


def cross_pdb_uniprot_interpro(url: str, dst: str):
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.execute(
        """
        SELECT 
          UPPER(UX_PDBE.AC),
          UPPER(ASYM.ENTRY_ID),
          UPPER(ASYM.AUTH_ASYM_ID),
          SRC.TAX_ID,
          IP.ENTRY_AC,
          IP.GO_ID,
          UPPER(IP.AC)
        FROM (
          SELECT UPI, AC
          FROM UNIPARC.XREF
          WHERE DBID = 21
          AND DELETED = 'N' 
        ) UX_PDBE
        INNER JOIN PDBE.STRUCT_ASYM@PDBE_LIVE ASYM
          ON (UX_PDBE.AC = ASYM.ENTRY_ID || '_' || ASYM.AUTH_ASYM_ID)
        INNER JOIN PDBE.ENTITY_SRC@PDBE_LIVE SRC
          ON (ASYM.ENTRY_ID = SRC.ENTRY_ID AND ASYM.ENTITY_ID = SRC.ENTITY_ID)
        INNER JOIN (
          SELECT DISTINCT 
            E.ENTRY_AC, IG.GO_ID, X.AC, X.UPI
          FROM INTERPRO.ENTRY E
          INNER JOIN INTERPRO.INTERPRO2GO IG 
            ON E.ENTRY_AC = IG.ENTRY_AC
          INNER JOIN INTERPRO.ENTRY2METHOD EM 
            ON E.ENTRY_AC = EM.ENTRY_AC
          INNER JOIN INTERPRO.MATCH M
            ON EM.METHOD_AC = M.METHOD_AC
          INNER JOIN UNIPARC.XREF X
            ON (M.PROTEIN_AC = X.AC AND X.DBID IN (2, 3) AND X.DELETED = 'N')
          WHERE E.CHECKED = 'Y'
        ) IP
          ON UX_PDBE.UPI = IP.UPI
        """
    )

    with open(dst, "wt") as fh:
        for row in cur:
            fh.write('\t'.join(map(str, row)) + '\n')

    cur.close()
    con.close()
