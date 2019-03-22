import cx_Oracle


from .. import logger
from ..orautils import TablePopulator


def _condense(matches: dict):
    for entry_acc in matches:
        fragments = []
        start = end = None
        for s, e in sorted(matches[entry_acc]):
            if start is None:
                # Leftmost fragment
                start = s
                end = e
            elif s > (end + 1):
                """
                        end
                    ----] 
                          [----
                          s
                -> new fragment
                
                end + 1: the `end` residue is included
                    so if end = 34 and start = 35, there is not gap
                """
                fragments.append((start, end))
                start = s
                end = e
            elif e > end:
                """
                        end
                    ----] 
                      ------]
                            e
                -> extend
                """
                end = e

        fragments.append((start, end))
        matches[entry_acc] = fragments


def condense_matches(url: str):
    logger.info("condensing matches")
    con = cx_Oracle.connect(url)
    cur = con.cursor()

    try:
        cur.execute("DROP TABLE INTERPRO.XREF_CONDENSED")
    except cx_Oracle.DatabaseError as e:
        _error = e.args[0]
        if _error.code != 942:
            # ORA-00942 (table or view does not exist)
            cur.close()
            con.close()
            raise e

    cur.execute(
        """
        CREATE TABLE INTERPRO.XREF_CONDENSED
        (
            PROTEIN_AC VARCHAR2(15) NOT NULL,
            ENTRY_AC VARCHAR2(9) NOT NULL,
            ENTRY_TYPE CHAR(1) NOT NULL,
            ENTRY_NAME VARCHAR2(100) NOT NULL,
            POS_FROM NUMBER(5) NOT NULL,
            POS_TO NUMBER(5) NOT NULL
        )
        PARTITION BY LIST (ENTRY_TYPE) (
          PARTITION XC_TYPE_A VALUES ('A'),
          PARTITION XC_TYPE_B VALUES ('B'),
          PARTITION XC_TYPE_C VALUES ('C'),
          PARTITION XC_TYPE_D VALUES ('D'),
          PARTITION XC_TYPE_F VALUES ('F'),
          PARTITION XC_TYPE_H VALUES ('H'),
          PARTITION XC_TYPE_P VALUES ('P'),
          PARTITION XC_TYPE_R VALUES ('R'),
          PARTITION XC_TYPE_U VALUES ('U')
        )
        NOLOGGING 
        """
    )

    cur.execute(
        """
        SELECT E.ENTRY_AC, E.NAME, E.ENTRY_TYPE, EM.METHOD_AC
        FROM INTERPRO.ENTRY E
        INNER JOIN INTERPRO.ENTRY2METHOD EM
          ON E.ENTRY_AC = EM.ENTRY_AC
        WHERE E.CHECKED = 'Y'
        """
    )

    entries = {}
    signatures = {}
    for entry_acc, name, entry_type, method_acc in cur:
        signatures[method_acc] = entry_acc
        entries[entry_acc] = (entry_type, name)

    cur.execute(
        """
        SELECT PROTEIN_AC, METHOD_AC, POS_FROM, POS_TO, FRAGMENTS
        FROM INTERPRO.MATCH
        ORDER BY PROTEIN_AC        
        """
    )

    matches = {}
    _protein_acc = None
    num_proteins = 0
    query = """
      INSERT /*+ APPEND */ INTO INTERPRO.XREF_CONDENSED
      VALUES (:1, :2, :3, :4, :5, :6)
    """
    table = TablePopulator(con, query, autocommit=True)

    for protein_acc, method_acc, pos_from, pos_to, fragments in cur:
        if protein_acc != _protein_acc:
            if matches:
                _condense(matches)
                for entry_acc, frags in matches.items():
                    entry_type, name = entries[entry_acc]
                    for start, end in frags:
                        table.insert((_protein_acc, entry_acc, entry_type,
                                      name, start, end))

                num_proteins += 1
                if not num_proteins % 10000000:
                    logger.info("proteins processed: "
                                "{:>15}".format(num_proteins))

            _protein_acc = protein_acc
            matches = {}

        try:
            entry_acc = signatures[method_acc]
        except KeyError:
            # Not integrated or integrated in an unchecked entry
            continue

        if entry_acc in matches:
            entry = matches[entry_acc]
        else:
            entry = matches[entry_acc] = []

        if fragments is not None:
            for frag in fragments.split(','):
                start, end, _ = frag.split('-')
                start = int(start)
                end = int(end)

                if start < end:
                    entry.append((start, end))

        if not entry:
            entry.append((pos_from, pos_to))

    if matches:
        _condense(matches)
        for entry_acc, frags in matches.items():
            entry_type, name = entries[entry_acc]
            for start, end in frags:
                table.insert((_protein_acc, entry_acc, entry_type,
                              name, start, end))

        num_proteins += 1

    table.close()

    logger.info("creating indexes")
    cur.execute(
        """
        CREATE INDEX I_XREF_CONDENSED$P
        ON INTERPRO.XREF_CONDENSED (PROTEIN_AC) NOLOGGING
        """
    )

    cur.execute(
        """
        CREATE INDEX I_XREF_CONDENSED$E
        ON INTERPRO.XREF_CONDENSED (ENTRY_AC) NOLOGGING
        """
    )

    cur.close()
    con.close()
    logger.info("proteins processed: {:>15}".format(num_proteins))
