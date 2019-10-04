# -*- coding: utf-8 -*-

import os

import cx_Oracle

from .sendmail import send_mail
from .. import logger, orautils


def _condense(matches: dict):
    for entry_acc in matches:
        entry_matches = []
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
                entry_matches.append((start, end))
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

        entry_matches.append((start, end))
        matches[entry_acc] = entry_matches


def build_xref_condensed(user: str, dsn: str):
    logger.info("building XREF_CONDENSED")
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()

    orautils.drop_table(cur, "INTERPRO", "XREF_CONDENSED", purge=True)
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

    table = orautils.TablePopulator(con,
                                    query="INSERT /*+ APPEND */ "
                                          "INTO INTERPRO.XREF_CONDENSED "
                                          "VALUES (:1, :2, :3, :4, :5, :6)",
                                    autocommit=True)

    cur.execute(
        """
        SELECT PROTEIN_AC, METHOD_AC, POS_FROM, POS_TO, FRAGMENTS
        FROM INTERPRO.MATCH
        ORDER BY PROTEIN_AC
        """
    )
    _protein_acc = None
    protein_matches = []
    for protein_acc, method_acc, pos_from, pos_to, fragments in cur:
        if protein_acc != _protein_acc:
            if protein_matches:
                # Condense in-place
                _condense(protein_matches)
                for _entry_acc, entry_matches in protein_matches.items():
                    entry_type, entry_name = entries[_entry_acc]
                    for _pos_from, _pos_end in entry_matches:
                        table.insert((_protein_acc, _entry_acc, entry_type,
                                      entry_name, _pos_from, _pos_end))

            protein_matches = {}
            _protein_acc = protein_acc

        try:
            entry_acc = signatures[method_acc]
        except KeyError:
            # Signature not integrated or integrated in an unchecked entry
            continue

        if entry_acc in protein_matches:
            e = protein_matches[entry_acc]
        else:
            e = protein_matches[entry_acc] = []

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
        e.append((pos_from, pos_to))

    # Last protein
    _condense(protein_matches)
    for _entry_acc, entry_matches in protein_matches.items():
        entry_type, entry_name = entries[_entry_acc]

        for _pos_from, _pos_end in entry_matches:
            table.insert((_protein_acc, _entry_acc, entry_type,
                          entry_name, _pos_from, _pos_end))

    table.close()

    # logger.info("gathering statistics")
    # orautils.gather_stats(cur, "INTERPRO", "XREF_CONDENSED")

    logger.info("creating indexes")
    for col in ("PROTEIN_AC", "ENTRY_AC"):
        cur.execute(
            f"""
            CREATE INDEX I_XREF_CONDENSED${col}
            ON INTERPRO.XREF_CONDENSED ({col}) NOLOGGING
            """
        )

    orautils.grant(cur, "INTERPRO", "XREF_CONDENSED", "SELECT", "KRAKEN")
    cur.close()
    con.close()


def build_xref_summary(user: str, dsn: str):
    logger.info("building XREF_SUMMARY")
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    orautils.drop_table(cur, "INTERPRO", "XREF_SUMMARY", purge=True)
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

    # logger.info("gathering statistics")
    # orautils.gather_stats(cur, "INTERPRO", "XREF_SUMMARY")

    logger.info("creating indexes")
    for col in ("PROTEIN_AC", "ENTRY_AC", "METHOD_AC"):
        cur.execute(
            f"""
            CREATE INDEX I_XREF_SUMMARY${col}
            ON INTERPRO.XREF_SUMMARY ({col}) NOLOGGING
            """
        )

    orautils.grant(cur, "INTERPRO", "XREF_SUMMARY", "SELECT", "KRAKEN")
    cur.close()
    con.close()


def export_databases(user: str, dsn: str, dst: str, notify: bool=True):
    """
    Format for Uniprot dat files:
      <protein>    DR   <database>; <signature/entry>; <name>; <count>.

    Exceptions:
        - Gene3D: do not include prefix before accession (G3DSA:)
                  replace signature#s name by hyphen (-)
        - PRINTS: do not include match count
        - InterPro: do not include match count
    """
    logger.info("exporting dat/tab files")
    os.makedirs(dst, 0o775, exist_ok=True)
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    cur.execute("SELECT VERSION FROM INTERPRO.DB_VERSION WHERE DBCODE = 'u'")
    release = cur.fetchone()[0]
    cur.execute(
        """
        SELECT
          PROTEIN_AC,
          METHOD_AC,
          MIN(METHOD_NAME),
          MIN(DBCODE),
          MIN(ENTRY_AC),
          MIN(SHORT_NAME),
          COUNT(*)
        FROM INTERPRO.XREF_SUMMARY
        GROUP BY PROTEIN_AC, METHOD_AC
        ORDER BY PROTEIN_AC
        """
    )

    dbcodes = {
        'B': ("SFLD", "SFLD"),
        # 'D': ("ProDom", "PD"),  # ProDom removed from InterPro
        'F': ("PRINTS", "PP"),
        'H': ("Pfam", "PF"),
        'J': ("CDD", "CDD"),
        'M': ("PROSITE", "PR"),
        'N': ("TIGRFAMs", "TF"),
        'P': ("PROSITE", "PR"),
        'Q': ("HAMAP", "HP"),
        'R': ("SMART", "SM"),
        'U': ("PIRSF", "PI"),
        'V': ("PANTHER", "PTHR"),
        'X': ("Gene3D", "G3D"),
        'Y': ("SUPFAM", "SF"),
    }

    files = []
    handlers = {}
    for dbcode in dbcodes:
        if dbcode != 'P':
            # P (PROSITE patterns) -> same file than M (PROSITE profiles)
            dbname, dbkey = dbcodes[dbcode]

            files.append(os.path.join(dst, '.' + dbname.lower() + ".tab"))
            fh1 = open(files[-1], "wt")

            files.append(os.path.join(dst, dbname + ".dat"))
            fh2 = open(files[-1], "wt")

            handlers[dbcode] = (fh1, fh2)

    handlers['P'] = handlers['M']

    files.append(os.path.join(dst, ".interpro.tab"))
    ifh1 = open(files[-1], "wt")
    files.append(os.path.join(dst, "InterPro.dat"))
    ifh2 = open(files[-1], "wt")

    # starts at -1 because on the first row, `protein_acc != _protein` is True
    cnt = -1
    entries = {}
    _protein = None
    for row in cur:
        protein_acc = row[0]
        signature_acc = row[1]
        signature_name = row[2]
        dbcode = row[3]
        entry_acc = row[4]
        entry_name = row[5]
        num_matches = int(row[6])

        dbname, dbkey = dbcodes[dbcode]
        fh1, fh2 = handlers[dbcode]

        fh1.write(f"{protein_acc}\t{dbkey}\t{signature_acc}\t"
                  f"{signature_name}\t{num_matches}\n")

        if dbcode == 'X':
            """
            Gene3D
              - accession: G3DSA:3.90.1580.10 -> 3.90.1580.10
              - do not print signature's name
            """
            fh2.write(f"{protein_acc}    DR   {dbname}; {signature_acc[6:]}; "
                      f"-; {num_matches}.\n")
        elif dbcode == 'F':
            # PRINTS: do not print match count
            fh2.write(f"{protein_acc}    DR   {dbname}; {signature_acc}; "
                      f"{signature_name}.\n")
        else:
            fh2.write(f"{protein_acc}    DR   {dbname}; {signature_acc}; "
                      f"{signature_name}; {num_matches}.\n")

        if protein_acc != _protein:
            for _entry in sorted(entries):
                _name = entries[_entry]
                ifh1.write(f"{_protein}\t{_entry}\t{_name}\n")
                ifh2.write(f"{_protein}    DR   InterPro; {_entry}; {_name}.\n")

            entries = {}
            _protein = protein_acc
            cnt += 1
            if not cnt % 10000000:
                logger.info(f"  {cnt:,}")

        entries[entry_acc] = entry_name

    # Last protein
    cnt += 1
    logger.info(f"  {cnt:,}")
    for _entry in sorted(entries):
        _name = entries[_entry]
        ifh1.write(f"{_protein}\t{_entry}\t{_name}\n")
        ifh2.write(f"{_protein}    DR   InterPro; {_entry}; {_name}.\n")

    ifh1.close()
    ifh2.close()

    for fh1, fh2 in handlers.values():
        fh1.close()
        fh2.close()

    for path in files:
        os.chmod(path, 0o775)

    if notify:
        send_mail(
            to_addrs=["uniprot-database@ebi.ac.uk"],
            cc_addrs=["uniprot-prod@ebi.ac.uk"],
            subject="InterPro XREF files are ready",
            content=f"""\
Dear UniProt team,

The InterPro cross-references files for {release} are available in the following directory:
  {dst}

Kind regards,
The InterPro Production Team
"""
        )

    logger.info("complete")


def build_aa_iprscan(user: str, dsn: str):
    logger.info("creating AA_IPRSCAN")
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    orautils.drop_mview(cur, "IPRSCAN", "AA_IPRSCAN")
    orautils.drop_table(cur, "IPRSCAN", "AA_IPRSCAN", purge=True)

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

    # orautils.gather_stats(cur, "IPRSCAN", "AA_IPRSCAN")
    orautils.grant(cur, "IPRSCAN", "AA_IPRSCAN", "SELECT", "KRAKEN")

    logger.info("creating indices")
    for col in ("UPI", "SIGNATURE"):
        logger.debug(f"index on {col}")
        cur.execute(
            f"""
            CREATE INDEX I_AA_IPRSCAN${col}
            ON IPRSCAN.AA_IPRSCAN ({col}) NOLOGGING
            """
        )

    cur.close()
    con.close()
    logger.info("AA_IPRSCAN is ready")


def ask_to_snapshot(user: str, dsn: str, notify: bool=True):
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    cur.execute("SELECT VERSION FROM INTERPRO.DB_VERSION WHERE DBCODE = 'u'")
    release = cur.fetchone()[0]
    cur.execute("SELECT MAX(UPI) FROM UNIPARC.PROTEIN@UAREAD")
    max_upi = cur.fetchone()[0]
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
To:  uniprot-database@ebi.ac.uk
Cc:  automated_annotation@ebi.ac.uk, uniprot-prod@ebi.ac.uk
Bcc: interpro-team@ebi.ac.uk

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
        if notify:
            send_mail(
                to_addrs=["interpro-team@ebi.ac.uk"],
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


def report_integration_changes(user: str, dsn: str, notify: bool=True):
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
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
                    changes[s_acc] = ("entry unchecked", e_acc, '')
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

"""
    content += (f"{'Signature':<30}{'Comment':<30}{'Previous entry':<20}"
                f"{'Current entry':<20}\n")
    content += '-' * 100 + '\n'

    for acc in sorted(changes, key=lambda k: k.lower()):
        status, prev, curr = changes[acc]
        content += f"{acc:<30}{status:<30}{prev:<20}{curr:<20}\n"

    content += "\nKind regards,\nThe InterPro Production Team\n"

    if notify:
        send_mail(
            to_addrs=["aa_dev@ebi.ac.uk"],
            cc_addrs=["unirule@ebi.ac.uk"],
            bcc_addrs=["interpro-team@ebi.ac.uk"],
            subject=f"Protein update {release}: integration changes",
            content=content
        )


def import_unirules(user: str, dsn: str, src: str):
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    orautils.drop_table(cur, "INTERPRO", "UNIRULE", purge=True)
    cur.execute(
        """
        CREATE TABLE INTERPRO.UNIRULE (
          ENTRY_AC VARCHAR2(9),
          METHOD_AC VARCHAR2(25),
          CONSTRAINT UQ_UNIRULE UNIQUE (ENTRY_AC, METHOD_AC),
          CONSTRAINT FK_UNIRULE$ENTRY FOREIGN KEY (ENTRY_AC)
            REFERENCES INTERPRO.ENTRY (ENTRY_AC) ON DELETE CASCADE,
          CONSTRAINT FK_UNIRULE$METHOD FOREIGN KEY (METHOD_AC)
            REFERENCES INTERPRO.METHOD (METHOD_AC) ON DELETE CASCADE,
          CONSTRAINT CK_UNIRULE
            CHECK ((ENTRY_AC IS NULL AND METHOD_AC IS NOT NULL)
              OR (ENTRY_AC IS NOT NULL AND METHOD_AC IS NULL))
        ) NOLOGGING
        """
    )

    cur.execute("SELECT ENTRY_AC, CHECKED FROM INTERPRO.ENTRY")
    entries = dict(cur.fetchall())

    cur.execute(
        """
        SELECT M.METHOD_AC, EM.ENTRY_AC
        FROM INTERPRO.METHOD M
        LEFT OUTER JOIN INTERPRO.ENTRY2METHOD EM
          ON M.METHOD_AC = EM.METHOD_AC
        """
    )
    signatures = dict(cur.fetchall())
    cur.close()

    table = orautils.TablePopulator(
        con=con,
        query="INSERT INTO INTERPRO.UNIRULE VALUES (:1, :2)"
    )

    with open(src, "rt") as fh:
        for line in fh:
            db, acc = line.rstrip().split()

            if db == "GDAC":
                acc = "G3DSA:" + acc

            if acc in entries:
                table.insert((acc, None))
            elif acc in signatures:
                table.insert((None, acc))

    table.close()
    con.commit()

    cur = con.cursor()
    orautils.grant(cur, "INTERPRO", "UNIRULE", "SELECT", "INTERPRO_SELECT")
    orautils.gather_stats(cur, "INTERPRO", "UNIRULE")
    cur.close()
    con.close()
