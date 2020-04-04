# -*- coding: utf-8 -*-

from typing import MutableMapping, Sequence, Tuple

import cx_Oracle

from pyinterprod import logger
from pyinterprod.utils import email, oracle, Table


def report_integration_changes(url: str, send_email: bool=True):
    con = cx_Oracle.connect(url)
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
        SELECT M.METHOD_AC, EM.ENTRY_AC
        FROM INTERPRO.METHOD M
        LEFT OUTER JOIN INTERPRO.ENTRY2METHOD EM
          ON M.METHOD_AC = EM.METHOD_AC
        """
    )
    current_signatures = dict(cur.fetchall())
    cur.execute(
        """
        SELECT ENTRY_AC, CHECKED
        FROM INTERPRO.ENTRY
        """
    )
    current_entries = {row[0]: row[1] == 'Y' for row in cur}
    cur.execute("SELECT VERSION FROM INTERPRO.DB_VERSION WHERE DBCODE = 'u'")
    release = cur.fetchone()[0]
    cur.close()
    con.close()

    changes = {}
    for s_acc, e_acc in previous_signatures.items():
        if s_acc in current_signatures:
            # Signature still exists
            ne_acc = current_signatures.pop(s_acc)
            if ne_acc:
                # Signature still integrated
                if e_acc != ne_acc:
                    # In another entry
                    if current_entries[ne_acc]:
                        # The other entry is checked
                        changes[s_acc] = ("moved", e_acc, ne_acc)
                    else:
                        # The other entry is unchecked
                        changes[s_acc] = ("moved to unchecked", e_acc, ne_acc)
                elif current_entries[e_acc]:
                    pass  # No changes
                else:
                    # In the same (now unchecked) entry
                    changes[s_acc] = ("entry unchecked", '', e_acc)
            else:
                # Signature unintegrated
                changes[s_acc] = ("unintegrated", e_acc, '')
        else:
            # Signature does not exist any more (member database update)
            changes[s_acc] = ("deleted", e_acc, '')

    for s_acc, ne_acc in current_signatures.items():
        if ne_acc and current_entries[ne_acc]:
            changes[s_acc] = ("integrated", '', ne_acc)

    content = """\
Dear UniProt team,

Please find below the list of recent integration changes.

Description of states:
  - deleted:                    the signature does not exist any more
  - entry unchecked:            the signature is integrated in an InterPro \
entry that has recently been flagged as not ready to be made public
  - integrated:                 the signature has been integrated in \
an InterPro entry ready to be made public.
  - moved:                      the signature has been integrated in \
an InterPro entry ready to be made public
  - moved to unchecked:         the signature has been integrated in \
an InterPro entry not ready to be made public
  - unintegrated:               the signature has been removed from \
an InterPro entry
  
"""
    content += (f"{'Signature':<30}{'Status':<30}{'Previous entry':<20}"
                f"{'Current entry':<20}\n")
    content += '-' * 100 + '\n'

    for acc in sorted(changes, key=lambda k: k.lower()):
        status, prev, curr = changes[acc]
        content += f"{acc:<30}{status:<30}{prev:<20}{curr:<20}\n"

    content += "\nKind regards,\nThe InterPro Production Team\n"

    if send_email:
        email.send(
            to=[email.AA_DEV],
            subject=f"Protein update {release}: integration changes",
            content=content,
            cc=[email.UNIRULE],
            bcc=[email.INTERPRO]
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


def build_aa_iprscan(url: str):
    logger.info("building AA_IPRSCAN")

    con = cx_Oracle.connect(url)
    cur = con.cursor()
    oracle.drop_table(cur, "IPRSCAN.AA_IPRSCAN", purge=True)

    select_stmt = """
        SELECT
          UPI,
          ANALYSIS_ID AS LIBRARY_ID,
          METHOD_AC AS SIGNATURE,
          SEQ_START,
          SEQ_END,
          SEQ_FEATURE
    """

    cur.execute(
        f"""
        CREATE TABLE IPRSCAN.AA_IPRSCAN COMPRESS NOLOGGING
        AS {select_stmt} FROM IPRSCAN.MV_IPRSCAN WHERE 1 = 0
        """
    )

    partitions = ("TMHMM", "SIGNALP_EUK", "SIGNALP_GRAM_POSITIVE",
                  "SIGNALP_GRAM_NEGATIVE", "COILS", "PROSITE_PATTERNS",
                  "PROSITE_PROFILES", "MOBIDBLITE")
    for p in partitions:
        logger.debug(f"inserting partition {p}")
        cur.execute(
            f"""
            INSERT /*+ APPEND */ INTO IPRSCAN.AA_IPRSCAN
            {select_stmt}
            FROM IPRSCAN.MV_IPRSCAN PARTITION({p})
            """
        )
        con.commit()

    logger.debug("inserting partition PHOBIUS")
    cur.execute(
        f"""
        INSERT /*+ APPEND */ INTO IPRSCAN.AA_IPRSCAN
        {select_stmt}
        FROM IPRSCAN.MV_IPRSCAN PARTITION(PHOBIUS)
        WHERE METHOD_AC IN ('SIGNAL_PEPTIDE','TRANSMEMBRANE')
        """
    )
    con.commit()

    logger.info("indexing")
    for col in ("UPI", "SIGNATURE"):
        cur.execute(
            f"""
            CREATE INDEX I_AA_IPRSCAN${col}
            ON IPRSCAN.AA_IPRSCAN ({col}) NOLOGGING
            """
        )

    logger.info("gathering statistics")
    oracle.gather_stats(cur, "IPRSCAN", "AA_IPRSCAN")

    cur.execute("GRANT SELECT ON IPRSCAN.AA_IPRSCAN TO KRAKEN")

    cur.close()
    con.close()

    logger.info("AA_IPRSCAN ready")


def build_xref_condensed(url: str):
    logger.info("building XREF_CONDENSED")
    con = cx_Oracle.connect(url)
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
          PARTITION PART_A VALUES ('A'),
          PARTITION PART_B VALUES ('B'),
          PARTITION PART_C VALUES ('C'),
          PARTITION PART_D VALUES ('D'),
          PARTITION PART_F VALUES ('F'),
          PARTITION PART_H VALUES ('H'),
          PARTITION PART_P VALUES ('P'),
          PARTITION PART_R VALUES ('R'),
          PARTITION PART_U VALUES ('U')
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

    logger.info("gathering statistics")
    oracle.gather_stats(cur, "INTERPRO", "XREF_CONDENSED")

    cur.execute("GRANT SELECT ON INTERPRO.XREF_CONDENSED TO KRAKEN")

    cur.close()
    con.close()

    logger.info("XREF_CONDENSED ready")


def build_xref_summary(url: str):
    logger.info("building XREF_SUMMARY")
    con = cx_Oracle.connect(url)
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

    logger.info("gathering statistics")
    oracle.gather_stats(cur, "INTERPRO", "XREF_SUMMARY")

    cur.execute("GRANT SELECT ON INTERPRO.XREF_SUMMARY TO KRAKEN")

    cur.close()
    con.close()

    logger.info("XREF_SUMMARY ready")


def ask_to_snapshot(url: str, send_email: bool=True):
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.execute("SELECT VERSION FROM INTERPRO.DB_VERSION WHERE DBCODE = 'u'")
    release, = cur.fetchone()
    cur.execute("SELECT MAX(UPI) FROM UNIPARC.PROTEIN@UAREAD")
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
To:  {email.UNIPROT_DB}
Cc:  {email.AA}, {email.UNIPROT_PROD}

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
        if send_email:
            email.send(
                to=[email.INTERPRO],
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
