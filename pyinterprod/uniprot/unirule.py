import oracledb

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
    con = oracledb.connect(uri)
    cur = con.cursor()
    cur.execute(
        """
        SELECT DISTINCT METHOD_AC, ENTRY_AC
        FROM INTERPRO.XREF_SUMMARY
        WHERE ENTRY_AC IS NOT NULL 
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


def _condense(matches: dict[str, list[tuple[int, int]]]):
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

    con = oracledb.connect(uri)
    cur = con.cursor()

    analyses = {}
    for analysis in iprscan.get_analyses(cur, type="matches"):
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
            HMMER_SEQ_START NUMBER(10),
            HMMER_SEQ_END NUMBER(10),
            ALIGNMENT VARCHAR2(4000)
        ) COMPRESS NOLOGGING
        """
    )

    # Open second cursor for INSERT statements (first used for SELECT)
    cur2 = con.cursor()

    for name in ["FunFam", "HAMAP", "PROSITE patterns", "PROSITE profiles"]:
        logger.info(f"inserting data from {name}")

        columns = ["UPI", "METHOD_AC", "SEQ_START", "SEQ_END"]

        if name == "FunFam":
            columns += ["HMMER_SEQ_START", "HMMER_SEQ_END"]
        else:
            columns += ["NULL", "NULL"]

        columns.append("ALIGNMENT")

        analysis_id, table = analyses[name]
        cur.execute(
            f"""
            SELECT {', '.join(columns)}
            FROM IPRSCAN.{iprscan.PREFIX}{table}
            WHERE ANALYSIS_ID = :1
           """,
            [analysis_id]
        )

        rows = []
        library = name.replace(" ", "_").upper()
        for row in cur:
            rows.append((row[0], library, row[1], row[2], row[3], row[4],
                         row[5], row[6]))

            if len(rows) == 1000:
                cur2.executemany(
                    f"""
                    INSERT /*+ APPEND */ INTO IPRSCAN.AA_ALIGNMENT
                    VALUES (:1, :2, :3, :4, :5, :6, :7, :8)
                    """,
                    rows
                )
                con.commit()
                rows.clear()

        if rows:
            cur2.executemany(
                f"""
                INSERT /*+ APPEND */ INTO IPRSCAN.AA_ALIGNMENT
                VALUES (:1, :2, :3, :4, :5, :6, :7, :8)
                """,
                rows
            )
            con.commit()
            rows.clear()

    cur2.close()

    logger.info("indexing")
    for col in ("UPI", "SIGNATURE"):
        cur.execute(
            f"""
            CREATE INDEX I_AA_ALIGNMENT${col}
            ON IPRSCAN.AA_ALIGNMENT ({col}) 
            TABLESPACE IPRSCAN_IND
            NOLOGGING
            """
        )

    cur.execute("GRANT SELECT ON IPRSCAN.AA_ALIGNMENT TO KRAKEN")
    cur.close()
    con.close()

    logger.info("AA_ALIGNMENT ready")


def build_aa_iprscan(uri: str):
    logger.info("building AA_IPRSCAN")

    con = oracledb.connect(uri)
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

    # Open second cursor for INSERT statements (first used for SELECT)
    cur2 = con.cursor()

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
                cur2.executemany(
                    f"""
                    INSERT /*+ APPEND */ INTO IPRSCAN.AA_IPRSCAN
                    VALUES (:1, :2, :3, :4, :5, :6)
                    """,
                    rows
                )
                con.commit()
                rows.clear()

        if rows:
            cur2.executemany(
                f"""
                INSERT /*+ APPEND */ INTO IPRSCAN.AA_IPRSCAN
                VALUES (:1, :2, :3, :4, :5, :6)
                """,
                rows
            )
            con.commit()
            rows.clear()

    cur2.close()

    logger.info("indexing")
    for col in ("UPI", "SIGNATURE"):
        cur.execute(
            f"""
            CREATE INDEX I_AA_IPRSCAN${col}
            ON IPRSCAN.AA_IPRSCAN ({col}) 
            TABLESPACE IPRSCAN_IND
            NOLOGGING
            """
        )

    cur.execute("GRANT SELECT ON IPRSCAN.AA_IPRSCAN TO KRAKEN")
    cur.close()
    con.close()

    logger.info("AA_IPRSCAN ready")


def build_xref_condensed(uri: str):
    logger.info("building XREF_CONDENSED")
    con = oracledb.connect(uri)
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
            ON INTERPRO.XREF_CONDENSED ({col}) 
            TABLESPACE INTERPRO_IND
            NOLOGGING
            """
        )

    cur.execute("GRANT SELECT ON INTERPRO.XREF_CONDENSED TO KRAKEN")
    cur.close()
    con.close()

    logger.info("XREF_CONDENSED ready")


def build_xref_summary(uri: str):
    logger.info("building XREF_SUMMARY")
    con = oracledb.connect(uri)
    cur = con.cursor()
    oracle.drop_table(cur, "INTERPRO.XREF_SUMMARY", purge=True)
    cur.execute(
        """
        CREATE TABLE INTERPRO.XREF_SUMMARY
        (
            DBCODE CHAR(1) NOT NULL,
            PROTEIN_AC VARCHAR2(15) NOT NULL,
            ENTRY_AC VARCHAR2(9),
            SHORT_NAME VARCHAR2(30),
            METHOD_AC VARCHAR2(50) NOT NULL,
            METHOD_NAME VARCHAR2(400),
            POS_FROM NUMBER(5) NOT NULL,
            POS_TO NUMBER(5) NOT NULL,
            MATCH_STATUS CHAR(1) NOT NULL,
            SCORE FLOAT(*),
            FRAGMENTS VARCHAR2(400)
        )
        PARTITION BY LIST (DBCODE) (
          PARTITION PART_A VALUES ('a'),
          PARTITION PART_B VALUES ('B'),
          PARTITION PART_F1 VALUES ('F'),
          PARTITION PART_F2 VALUES ('f'),
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
        INSERT INTO INTERPRO.XREF_SUMMARY
        SELECT
            MA.DBCODE,
            MA.PROTEIN_AC,
            E.ENTRY_AC,
            E.SHORT_NAME,
            MA.METHOD_AC,
            CASE 
                WHEN ME.NAME IS NOT NULL AND MA.METHOD_AC != ME.NAME 
                    THEN ME.NAME
                WHEN ME.DESCRIPTION IS NOT NULL
                    THEN REGEXP_REPLACE(DESCRIPTION, '"', '''')
                ELSE NULL 
            END,
            MA.POS_FROM,
            MA.POS_TO,
            MA.STATUS,
            MA.SCORE,
            MA.FRAGMENTS
        FROM INTERPRO.MATCH MA
        INNER JOIN INTERPRO.METHOD ME 
            ON MA.METHOD_AC = ME.METHOD_AC
        LEFT OUTER JOIN INTERPRO.ENTRY2METHOD EM 
            ON ME.METHOD_AC = EM.METHOD_AC
        LEFT OUTER JOIN INTERPRO.ENTRY E 
            ON EM.ENTRY_AC = E.ENTRY_AC AND E.CHECKED = 'Y'
        UNION ALL
        -- PANTHER subfamilies
        SELECT
            MA.DBCODE,
            MA.PROTEIN_AC,
            NULL,
            NULL,
            MA.MODEL_AC,
            CASE 
                WHEN ME.NAME IS NOT NULL AND MA.MODEL_AC != ME.NAME 
                    THEN ME.NAME
                WHEN ME.DESCRIPTION IS NOT NULL
                    THEN REGEXP_REPLACE(DESCRIPTION, '"', '''')
                ELSE NULL 
            END,
            MA.POS_FROM,
            MA.POS_TO,
            MA.STATUS,
            MA.SCORE,
            MA.FRAGMENTS
        FROM INTERPRO.MATCH PARTITION (MATCH_DBCODE_V) MA
        INNER JOIN INTERPRO.METHOD ME 
            ON MA.MODEL_AC = ME.METHOD_AC
        WHERE MA.MODEL_AC IS NOT NULL 
          AND REGEXP_LIKE(MA.MODEL_AC, '^PTHR\d+:SF\d+$')
        UNION ALL
        -- AntiFam, FunFam
        SELECT
            FM.DBCODE,
            FM.PROTEIN_AC,
            NULL,
            NULL,
            FM.METHOD_AC,
            CASE WHEN ME.NAME IS NOT NULL AND FM.METHOD_AC != ME.NAME
                THEN ME.NAME ELSE ME.DESCRIPTION END,
            FM.POS_FROM,
            FM.POS_TO,
            'T',
            NULL,
            NULL
        FROM INTERPRO.FEATURE_MATCH FM
        INNER JOIN INTERPRO.FEATURE_METHOD ME
            ON FM.METHOD_AC = ME.METHOD_AC
        WHERE FM.DBCODE IN ('a', 'f')
        """
    )

    con.commit()

    logger.info("indexing")
    for col in ("PROTEIN_AC", "ENTRY_AC", "METHOD_AC"):
        cur.execute(
            f"""
            CREATE INDEX I_XREF_SUMMARY${col}
            ON INTERPRO.XREF_SUMMARY ({col}) 
            TABLESPACE INTERPRO_IND
            NOLOGGING
            """
        )

    cur.execute("GRANT SELECT ON INTERPRO.XREF_SUMMARY TO KRAKEN")
    cur.close()
    con.close()

    logger.info("XREF_SUMMARY ready")


def ask_to_snapshot(uri: str, emails: dict):
    con = oracledb.connect(uri)
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
All tables for UniProt are ready.

Please, take a snapshot of IPPRO, and inform UniProt they can refresh IPREADU.

Recipients
----------
To:  {emails['uniprot_db']}
Cc:  {emails['aa']}, {emails['uniprot_prod']}, {emails['interpro']}

Subject
-------
Protein update {release} completed: please refresh IPREADU

Body (replace [DATE-TIME] by snapshot date/time)
------------------------------------------------
Dear UniProt team,

You may refresh IPREADU with the snapshot of IPPRO from [DATE-TIME]

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

    con = oracledb.connect(uri)
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


MAX_DOM_BY_GROUP = 20
DOM_OVERLAP_THRESHOLD = 0.3
REPR_DOM_DATABASES = ["H", "J", "M", "R", "N"]
DC_STATUSES = {
    "S": "CONTINUOUS",
    "N": "N_TERMINAL_DISC",
    "C": "C_TERMINAL_DISC",
    "NC": "NC_TERMINAL_DISC"
}


def get_repr_domains(ora_url: str, output: str = "repr_domains.tsv"):
    con = oracledb.connect(ora_url)
    cur = con.cursor()
    cur.execute(
        """
        SELECT PROTEIN_AC, H.METHOD_AC, MODEL_AC,
            POS_FROM, POS_TO, FRAGMENTS
        FROM INTERPRO.MATCH H
        INNER JOIN INTERPRO.METHOD D
        ON H.METHOD_AC = D.METHOD_AC
        WHERE H.DBCODE in ('H', 'J', 'M', 'R', 'N')
        AND (D.SIG_TYPE = 'D' OR D.SIG_TYPE = 'R')
        ORDER BY PROTEIN_AC
        """
    )

    logger.info(f"Writing {output}")
    proteins_domains = {}
    current_protein = None
    with open(output, "w") as f:
        for protein_acc, signature_acc, model_acc, dbcode, pos_start, pos_end, frags in cur:
            fragments = _get_fragments(pos_start, pos_end, frags)
            match = {
                "signature": signature_acc,
                "model": model_acc or signature_acc,
                "fragments": fragments,
                "rank": REPR_DOM_DATABASES.index(dbcode)
            }
            try:
                proteins_domains[protein_acc].append(match)
            except KeyError:
                if proteins_domains:
                    repr_domains = _select_repr_domains(list(proteins_domains[current_protein]))
                    for domain in repr_domains:
                        if domain.get('representative', None):
                            f.write(
                                f"{current_protein}\t{domain['signature']}\t{domain['model']}\t{domain['fragments']}\t{domain['representative']}\n"
                            )
                    proteins_domains.clear()
                current_protein = protein_acc
                proteins_domains[protein_acc] = [match]
        if proteins_domains:
            repr_domains = _select_repr_domains(list(proteins_domains[current_protein]))
            for domain in repr_domains:
                if domain.get('representative', None):
                    f.write(
                        f"{current_protein}\t{domain['signature']}\t{domain['model']}\t{domain['fragments']}\t{domain['representative']}\n"
                    )

    logger.info("Done")
    cur.close()
    con.close()


def _select_repr_domains(domains: list[dict]):
    # Sort by boundaries
    domains.sort(key=lambda d: (d["fragments"][0]["start"],
                                d["fragments"][-1]["end"]))

    # Group overlapping domains together
    domain = domains[0]
    domain["residues"] = _calc_coverage(domain)
    stop = domain["fragments"][-1]["end"]
    group = [domain]
    groups = []

    for domain in domains[1:]:
        domain["residues"] = _calc_coverage(domain)
        start = domain["fragments"][0]["start"]

        if start <= stop:
            group.append(domain)
            stop = max(stop, domain["fragments"][-1]["end"])
        else:
            groups.append(group)
            group = [domain]
            stop = domain["fragments"][-1]["end"]

    groups.append(group)

    # Select representative domain in each group
    for group in groups:
        """
        Only consider the "best" N domains of the group, 
        otherwise the number of possible combinations/sets is too high 
        (if M domains, max number of combinations is `2 ^ M`)
        """
        group = sorted(group,
                       key=lambda d: (-len(d["residues"]), d["rank"])
                       )[:MAX_DOM_BY_GROUP]

        nodes = set(range(len(group)))
        graph = {i: nodes - {i} for i in nodes}

        for i, dom_a in enumerate(group):
            for j in range(i + 1, len(group)):
                dom_b = group[j]
                if _eval_overlap(dom_a, dom_b, DOM_OVERLAP_THRESHOLD):
                    graph[i].remove(j)
                    graph[j].remove(i)

        # Find possible domains combinations
        subgroups = _resolve_domains(graph)

        # Find the best combination
        max_coverage = 0
        max_pfams = 0
        best_subgroup = None
        for subgroup in subgroups:
            coverage = set()
            pfams = 0
            _subgroup = []

            for i in subgroup:
                domain = group[i]
                coverage |= domain["residues"]
                if domain["rank"] == 0:
                    pfams += 1

                _subgroup.append(domain)

            coverage = len(coverage)
            if coverage < max_coverage:
                continue
            elif coverage > max_coverage or pfams > max_pfams:
                max_coverage = coverage
                max_pfams = pfams
                best_subgroup = _subgroup

        # Flag selected representative domains
        for domain in best_subgroup:
            domain["representative"] = True
    return domains


def _resolve_domains(graph: dict[int, set[int]]) -> list[set[int]]:
    def is_valid(candidate: list[int]) -> bool:
        for node_a in candidate:
            for node_b in candidate:
                if node_a != node_b and node_a not in graph[node_b]:
                    return False

        return True

    def make_sets(current_set: list[int], remaining_nodes: list[int]):
        if is_valid(current_set):
            if not remaining_nodes:
                all_sets.append(set(current_set))
                return True
        else:
            return False

        current_node = remaining_nodes[0]
        remaining_nodes = remaining_nodes[1:]

        # Explore two possibilities at each step of the recursion
        # 1) current node is added to the set under consideration
        make_sets(current_set + [current_node], remaining_nodes)
        # 2) current node is not added to the set
        make_sets(current_set, remaining_nodes)

    all_sets = []
    make_sets([], list(graph.keys()))
    return all_sets


def _calc_coverage(domain: dict) -> set[int]:
    residues = set()
    for f in domain["fragments"]:
        residues |= set(range(f["start"], f["end"] + 1))

    return residues


def _eval_overlap(dom_a: dict, dom_b: dict, threshold: float) -> bool:
    overlap = dom_a["residues"] & dom_b["residues"]
    if overlap:
        len_a = len(dom_a["residues"])
        len_b = len(dom_b["residues"])
        return len(overlap) / min(len_a, len_b) >= threshold

    return False


def _get_fragments(pos_start: int, pos_end: int, fragments: str) -> list[dict]:
    if fragments:
        result = []
        for frag in fragments.split(','):
            # Format: START-END-STATUS
            s, e, t = frag.split('-')
            result.append({
                "start": int(s),
                "end": int(e),
                "dc-status": DC_STATUSES[t]
            })

        result.sort(key=lambda x: (x["start"], x["end"]))
    else:
        result = [{
            "start": pos_start,
            "end": pos_end,
            "dc-status": DC_STATUSES['S']  # Continuous
        }]

    return result
