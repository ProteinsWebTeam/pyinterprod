import oracledb


def get_sifts_mapping(uri: str) -> dict[str, set[str]]:
    sifts = {}

    con = oracledb.connect(uri)
    cur = con.cursor()

    cur.execute(
        """
        SELECT DISTINCT U.ACCESSION, E.ID, U.AUTH_ASYM_ID
        FROM PDBE.ENTRY E
        INNER JOIN SIFTS_ADMIN.SIFTS_XREF_SEGMENT U ON (
          E.ID = U.ENTRY_ID AND
          E.METHOD_CLASS IN ('nmr', 'x-ray', 'em')
        )
        """
    )

    for uniprot_acc, pdb_id, chain in cur.fetchall():
        pdb_chain = f"{pdb_id}_{chain}"
        try:
            sifts[pdb_chain].add(uniprot_acc)
        except KeyError:
            sifts[pdb_chain] = {uniprot_acc}

    cur.close()
    con.close()

    return sifts
