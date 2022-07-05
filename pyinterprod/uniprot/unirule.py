from typing import MutableMapping, Sequence, Tuple

import cx_Oracle

from pyinterprod import logger
from pyinterprod.interpro import iprscan
from pyinterprod.utils import email, oracle, Table


def report_integration_changes(uri: str, emails: dict):
    """Sends a list of integration changes.
    If a signature is integrated in an unchecked entry, we consider
    the signature as unintegrated.

    :param uri: Oracle connection string
    :param emails: email info (SMTP server/port, sender, recipients, etc.)
    """
    con = cx_Oracle.connect(uri)
    cur = con.cursor()
    cur.execute(
        """
        SELECT DISTINCT METHOD_AC, ENTRY_AC
        FROM INTERPRO.XREF_SUMMARY
        """
    )
    previous_signatures = dict(cur.fetchall())
    cur.execute(
        """
        SELECT M.METHOD_AC, E.ENTRY_AC
        FROM INTERPRO.METHOD M
        LEFT OUTER JOIN INTERPRO.ENTRY2METHOD EM
          ON M.METHOD_AC = EM.METHOD_AC
        LEFT OUTER JOIN INTERPRO.ENTRY E
          ON EM.ENTRY_AC = E.ENTRY_AC
          AND E.CHECKED = 'Y'
        """
    )
    current_signatures = dict(cur.fetchall())
    cur.execute("SELECT VERSION FROM INTERPRO.DB_VERSION WHERE DBCODE = 'u'")
    release = cur.fetchone()[0]
    cur.close()
    con.close()

    changes = {}
    for signature, entry_then in previous_signatures.items():
        try:
            entry_now = current_signatures.pop(signature)
        except KeyError:
            # Signature does not exist any more (member database update)
            changes[signature] = ("deleted", entry_then, '')
            continue

        if entry_now is None:
            # Signature unintegrated
            changes[signature] = ("unintegrated", entry_then, '')
        elif entry_now != entry_then:
            # Integrated in a different entry
            changes[signature] = ("integrated", entry_then, entry_now)

    # All new integrations (signature was not integrated before)
    for signature, entry_now in current_signatures.items():
        if entry_now:
            changes[signature] = ("integrated", '', entry_now)

    content = """\
Dear UniProt team,

Please find below the list of recent integration changes.

Description of states:
  - integrated:       the signature has been integrated in an InterPro entry
  - unintegrated:     the signature has been removed from an InterPro entry
  - deleted:          the signature does not exist anymore

"""
    content += (f"{'Signature':<30}{'Status':<30}{'Previous entry':^30}"
                f"{'Current entry':^30}\n")
    content += '-' * 120 + '\n'

    for acc in sorted(changes, key=lambda k: k.lower()):
        status, prev, curr = changes[acc]
        content += f"{acc:<30}{status:<30}{prev:^30}{curr:^30}\n"

    content += "\nKind regards,\nThe InterPro Production Team\n"

    email.send(
        info=emails,
        to=["aa_dev"],
        cc=["unirule"],
        bcc=["sender"],
        subject=f"Protein update {release}: integration changes",
        content=content
    )


def _condense(matches: MutableMapping[str, Sequence[Tuple[int, int]]]):
    for entry_acc in matches:
        condensed_matches = []
        start = end = None
        for s, e in sorted(matches[entry_acc]):
            if start is None:
                # Leftmost match
                start = s
                end = e
            elif s > end:
                """
                      end
                   [----] [----]
                          s
                -> new match
                """
                condensed_matches.append((start, end))
                start = s
                end = e
            elif e > end:
                """
                        end
                   [----]
                     [------]
                            e
                -> extend
                """
                end = e

        condensed_matches.append((start, end))
        matches[entry_acc] = condensed_matches


def build_aa_alignment(uri: str):
    logger.info("building AA_ALIGNMENT")

    con = cx_Oracle.connect(uri)
    cur = con.cursor()

    analyses = {}
    for analysis in iprscan.get_analyses(cur, type="matches",
                                         status="production"):
        analyses[analysis.name] = (analysis.id, analysis.table)

    oracle.drop_table(cur, "IPRSCAN.AA_ALIGNMENT", purge=True)
    cur.execute(
        """
        CREATE TABLE IPRSCAN.AA_ALIGNMENT
        (
            UPI VARCHAR2(13) NOT NULL,
            LIBRARY VARCHAR2(25) NOT NULL,
            SIGNATURE VARCHAR2(255) NOT NULL,
            SEQ_START NUMBER(10) NOT NULL,
            SEQ_END NUMBER(10) NOT NULL,
            ALIGNMENT VARCHAR2(4000)
        ) COMPRESS NOLOGGING
        """
    )

    # TODO: add FunFam
    for name in ["HAMAP", "PROSITE patterns", "PROSITE profiles"]:
        logger.info(f"inserting data from {name}")
        analysis_id, table = analyses[name]

        cur.execute(
            f"""
           SELECT UPI, METHOD_AC, SEQ_START, SEQ_END, ALIGNMENT
           FROM IPRSCAN.{iprscan.PREFIX}{table}
           WHERE ANALYSIS_ID = :1
           """,
            [analysis_id]
        )

        rows = []
        library = name.replace(" ", "_").upper()
        for row in cur:
            rows.append((row[0], library, row[1], row[2], row[3], row[4]))

            if len(rows) == 1000:
                cur.executemany(
                    f"""
                    INSERT /*+ APPEND */ INTO IPRSCAN.AA_ALIGNMENT
                    VALUES (:1, :2, :3, :4, :5, :6)
                    """,
                    rows
                )
                con.commit()
                rows.clear()

        if rows:
            cur.executemany(
                f"""
                INSERT /*+ APPEND */ INTO IPRSCAN.AA_ALIGNMENT
                VALUES (:1, :2, :3, :4, :5, :6)
                """,
                rows
            )
            con.commit()
            rows.clear()

    logger.info("indexing")
    for col in ("UPI", "METHOD_AC"):
        cur.execute(
            f"""
            CREATE INDEX I_AA_ALIGNMENT${col}
            ON IPRSCAN.AA_ALIGNMENT ({col}) NOLOGGING
            """
        )

    cur.execute("GRANT SELECT ON IPRSCAN.AA_ALIGNMENT TO KRAKEN")
    cur.close()
    con.close()

    logger.info("AA_ALIGNMENT ready")


def build_aa_iprscan(uri: str):
    logger.info("building AA_IPRSCAN")

    con = cx_Oracle.connect(uri)
    cur = con.cursor()
    oracle.drop_table(cur, "IPRSCAN.AA_IPRSCAN", purge=True)
    cur.execute(
        """
        CREATE TABLE IPRSCAN.AA_IPRSCAN
        (
            UPI VARCHAR2(13) NOT NULL,
            -- LIBRARY_ID NUMBER(5) NOT NULL,
            LIBRARY VARCHAR2(25) NOT NULL,
            SIGNATURE VARCHAR2(255) NOT NULL,
            SEQ_START NUMBER(10) NOT NULL,
            SEQ_END NUMBER(10) NOT NULL,
            SEQ_FEATURE VARCHAR2(4000)
        ) COMPRESS NOLOGGING
        """
    )

    for db in ["COILS", "MobiDB Lite", "Phobius", "PROSITE patterns",
               "PROSITE profiles", "SignalP_Euk", "SignalP_Gram_positive",
               "SignalP_Gram_negative", "TMHMM"]:
        logger.info(f"inserting data from {db}")
        partition = iprscan.MATCH_PARTITIONS[db]["partition"]

        sql = f"""
            SELECT UPI, ANALYSIS_ID, METHOD_AC, SEQ_START, SEQ_END, SEQ_FEATURE
            FROM IPRSCAN.MV_IPRSCAN PARTITION ({partition})        
        """

        if db == "Phobius":
            sql += "WHERE METHOD_AC IN ('SIGNAL_PEPTIDE','TRANSMEMBRANE')"

        cur.execute(sql)

        rows = []
        library = db.replace(" ", "_").upper()
        for row in cur:
            # rows.append((row[0], row[1], library, row[2], row[3], row[4],
            #              row[5]))
            rows.append((row[0], library, row[2], row[3], row[4], row[5]))

            if len(rows) == 1000:
                cur.executemany(
                    f"""
                    INSERT /*+ APPEND */ INTO IPRSCAN.AA_IPRSCAN
                    VALUES (:1, :2, :3, :4, :5, :6)
                    """,
                    rows
                )
                con.commit()
                rows.clear()

        if rows:
            cur.executemany(
                f"""
                INSERT /*+ APPEND */ INTO IPRSCAN.AA_IPRSCAN
                VALUES (:1, :2, :3, :4, :5, :6)
                """,
                rows
            )
            con.commit()
            rows.clear()

    logger.info("indexing")
    for col in ("UPI", "SIGNATURE"):
        cur.execute(
            f"""
            CREATE INDEX I_AA_IPRSCAN${col}
            ON IPRSCAN.AA_IPRSCAN ({col}) NOLOGGING
            """
        )

    cur.execute("GRANT SELECT ON IPRSCAN.AA_IPRSCAN TO KRAKEN")
    cur.close()
    con.close()

    logger.info("AA_IPRSCAN ready")


def build_xref_condensed(uri: str):
    logger.info("building XREF_CONDENSED")
    con = cx_Oracle.connect(uri)
    cur = con.cursor()
    oracle.drop_table(cur, "INTERPRO.XREF_CONDENSED", purge=True)
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
          PARTITION PART_A VALUES ('A'),  -- Active_site
          PARTITION PART_B VALUES ('B'),  -- Binding_site
          PARTITION PART_C VALUES ('C'),  -- Conserved_site
          PARTITION PART_D VALUES ('D'),  -- Domain
          PARTITION PART_F VALUES ('F'),  -- Family
          PARTITION PART_H VALUES ('H'),  -- Homologous_superfamily
          PARTITION PART_P VALUES ('P'),  -- PTM
          PARTITION PART_R VALUES ('R')   -- Repeat
        ) COMPRESS NOLOGGING      
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

    sql = """
        INSERT /*+ APPEND */ 
        INTO INTERPRO.XREF_CONDENSED 
        VALUES (:1, :2, :3, :4, :5, :6)
    """
    with Table(con, sql, autocommit=True) as table:
        cur.execute(
            """
            SELECT PROTEIN_AC, METHOD_AC, POS_FROM, POS_TO, FRAGMENTS
            FROM INTERPRO.MATCH
            ORDER BY PROTEIN_AC
            """
        )

        prev_acc = None
        matches = {}
        for row in cur:
            protein_acc = row[0]

            if protein_acc != prev_acc:
                # New protein: condense previous protein's matches
                if matches:
                    # Condense in-place
                    _condense(matches)

                    for entry_acc, entry_matches in matches.items():
                        entry_type, entry_name = entries[entry_acc]
                        for pos_from, pos_end in entry_matches:
                            table.insert((
                                prev_acc,
                                entry_acc,
                                entry_type,
                                entry_name,
                                pos_from,
                                pos_end
                            ))

                matches = {}
                prev_acc = protein_acc

            method_acc = row[1]
            pos_from = row[2]
            pos_to = row[3]
            # fragments = row[4]

            try:
                entry_acc = signatures[method_acc]
            except KeyError:
                # Signature not integrated or integrated in an unchecked entry
                continue

            try:
                obj = matches[entry_acc]
            except KeyError:
                obj = matches[entry_acc] = []

            """
            As of May 2019, UniProt does not use discontinuous domains
            because their collaborators need to be able to
            distinguish between repeated matches and fragmented matches
            """
            # if fragments is not None:
            #     for frag in fragments.split(','):
            #         start, end, _ = frag.split('-')
            #         e.append((int(start), int(end)))
            # else:
            #     e.append((pos_from, pos_to))
            obj.append((pos_from, pos_to))

        # Last protein
        _condense(matches)
        for entry_acc, entry_matches in matches.items():
            entry_type, entry_name = entries[entry_acc]

            for pos_from, pos_end in entry_matches:
                table.insert((
                    prev_acc,
                    entry_acc,
                    entry_type,
                    entry_name,
                    pos_from,
                    pos_end
                ))

    logger.info("indexing")
    for col in ("PROTEIN_AC", "ENTRY_AC"):
        cur.execute(
            f"""
            CREATE INDEX I_XREF_CONDENSED${col}
            ON INTERPRO.XREF_CONDENSED ({col}) NOLOGGING
            """
        )

    cur.execute("GRANT SELECT ON INTERPRO.XREF_CONDENSED TO KRAKEN")
    cur.close()
    con.close()

    logger.info("XREF_CONDENSED ready")


def build_xref_summary(uri: str):
    logger.info("building XREF_SUMMARY")
    con = cx_Oracle.connect(uri)
    cur = con.cursor()
    oracle.drop_table(cur, "INTERPRO.XREF_SUMMARY", purge=True)
    cur.execute(
        """
        CREATE TABLE INTERPRO.XREF_SUMMARY
        (
            DBCODE CHAR(1) NOT NULL,
            PROTEIN_AC VARCHAR2(15) NOT NULL,
            ENTRY_AC VARCHAR2(9) NOT NULL,
            SHORT_NAME VARCHAR2(30),
            METHOD_AC VARCHAR2(25) NOT NULL,
            METHOD_NAME VARCHAR2(100),
            POS_FROM NUMBER(5) NOT NULL,
            POS_TO NUMBER(5) NOT NULL,
            MATCH_STATUS CHAR(1) NOT NULL,
            SCORE FLOAT(*),
            FRAGMENTS VARCHAR2(400)
        )
        PARTITION BY LIST (DBCODE) (
          PARTITION PART_B VALUES ('B'),
          PARTITION PART_F VALUES ('F'),
          PARTITION PART_H VALUES ('H'),
          PARTITION PART_J VALUES ('J'),
          PARTITION PART_M VALUES ('M'),
          PARTITION PART_N VALUES ('N'),
          PARTITION PART_P VALUES ('P'),
          PARTITION PART_Q VALUES ('Q'),
          PARTITION PART_R VALUES ('R'),
          PARTITION PART_U VALUES ('U'),
          PARTITION PART_V VALUES ('V'),
          PARTITION PART_X VALUES ('X'),
          PARTITION PART_Y VALUES ('Y')
        ) COMPRESS NOLOGGING        
        """
    )

    cur.execute(
        """
        INSERT /*+ APPEND */ INTO INTERPRO.XREF_SUMMARY
        SELECT
          MA.DBCODE,
          MA.PROTEIN_AC,
          E.ENTRY_AC,
          E.SHORT_NAME,
          MA.METHOD_AC,
          ME.NAME,
          MA.POS_FROM,
          MA.POS_TO,
          MA.STATUS,
          MA.SCORE,
          MA.FRAGMENTS
        FROM INTERPRO.MATCH MA
        INNER JOIN INTERPRO.METHOD ME ON MA.METHOD_AC = ME.METHOD_AC
        INNER JOIN INTERPRO.ENTRY2METHOD EM ON ME.METHOD_AC = EM.METHOD_AC
        INNER JOIN INTERPRO.ENTRY E ON EM.ENTRY_AC = E.ENTRY_AC
        WHERE EM.EVIDENCE != 'EXC'
        AND E.CHECKED = 'Y'
        """
    )
    con.commit()

    logger.info("indexing")
    for col in ("PROTEIN_AC", "ENTRY_AC", "METHOD_AC"):
        cur.execute(
            f"""
            CREATE INDEX I_XREF_SUMMARY${col}
            ON INTERPRO.XREF_SUMMARY ({col}) NOLOGGING
            """
        )

    cur.execute("GRANT SELECT ON INTERPRO.XREF_SUMMARY TO KRAKEN")
    cur.close()
    con.close()

    logger.info("XREF_SUMMARY ready")


def ask_to_snapshot(uri: str, emails: dict):
    con = cx_Oracle.connect(uri)
    cur = con.cursor()
    cur.execute("SELECT VERSION FROM INTERPRO.DB_VERSION WHERE DBCODE = 'u'")
    release, = cur.fetchone()
    cur.execute("SELECT MAX(UPI) FROM UNIPARC.PROTEIN")
    max_upi, = cur.fetchone()
    cur.execute(
        """
        SELECT I.IPRSCAN_SIG_LIB_REL_ID, D.DBNAME, V.VERSION
        FROM INTERPRO.IPRSCAN2DBCODE I
        INNER JOIN INTERPRO.CV_DATABASE D ON I.DBCODE = D.DBCODE
        LEFT OUTER JOIN INTERPRO.DB_VERSION V ON I.DBCODE = V.DBCODE
        WHERE I.IPRSCAN_SIG_LIB_REL_ID != I.PREV_ID
        ORDER BY D.DBNAME
        """
    )
    updates = []
    for analysis_id, name, version in cur:
        if version:
            name += " (" + version + ")"
        updates.append((name, analysis_id))

    content = f"""\
INTERPRO.XREF_SUMMARY, INTERPRO.XREF_CONDENSED, and IPRSCAN.AA_IPRSCAN are ready.

Please, take a snapshot of IPPRO, and inform UniProt they can refresh IPREADU.

Recipients
----------
To:  {emails['uniprot_db']}
Cc:  {emails['aa']}, {emails['uniprot_prod']}

Subject
-------
Protein update {release} completed: please refresh IPREADU

Body (change text between square brackets)
------------------------------------------
Dear UniProt team,

You may refresh IPREADU with the snapshot of IPPRO from [date/time].

The max UPI we processed up to in this update is {max_upi}.
"""

    if updates:
        content += """
We have updated the following databases since the previous protein update:
"""

        for name, analysis_id in updates:
            content += f"  - {name:<30}new analysis ID: {analysis_id}\n"
    else:
        content += """
There are no changes to the analysis IDs for this protein update.
"""

    content += """
Kind regards,
The InterPro Production Team
"""

    try:
        email.send(
            info=emails,
            to=["interpro"],
            subject=f"Protein update {release}: please snapshot IPPRO",
            content=content
        )

        # Update PREV_ID only if mail sent
        cur.execute(
            """
            UPDATE INTERPRO.IPRSCAN2DBCODE
            SET PREV_ID = IPRSCAN_SIG_LIB_REL_ID
            """
        )
        con.commit()
    finally:
        cur.close()
        con.close()


def update_signatures(filepath: str, uri: str):
    # Load accessions from file generated by UniRule
    accessions = []
    with open(filepath, "rt") as fh:
        for line in fh:
            database, accession = line.rstrip().split()
            if database == "GDAC":
                accessions.append(f"G3DSA:{accession}")
            elif database != "IPRO":
                # ignore InterPro entries
                accessions.append(accession)

    con = cx_Oracle.connect(uri)
    cur = con.cursor()

    # Refresh data
    oracle.truncate_table(cur, "INTERPRO.METHOD_UNIRULE")

    cur.executemany(
        """
        INSERT INTO INTERPRO.METHOD_UNIRULE VALUES (:1)
        """, [(acc,) for acc in set(accessions)]
    )

    con.commit()
    cur.close()
    con.close()
