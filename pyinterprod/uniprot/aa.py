import math
import os
from multiprocessing import Process, Queue
from tempfile import mkstemp

import oracledb

from pyinterprod import logger
from pyinterprod.interpro import iprscan
from pyinterprod.pronto.match import load_index, iter_matches
from pyinterprod.utils import email, oracle, Table


MAX_DOM_BY_GROUP = 20
DOM_OVERLAP_THRESHOLD = 0.3
# Pfam, CDD, PROSITE profiles, SMART, NCBIFAM, CATH-Gene3D, SUPERFAMILY
REPR_DOM_DATABASES = ["H", "J", "M", "R", "N", "X", "Y"]
# Domain, Repeat, Homologous superfamily
REPR_DOM_TYPES = ["D", "R", "H"]


def create_aa_alignment(uri: str):
    logger.info("creating AA_ALIGNMENT")

    con = oracledb.connect(uri)
    cur = con.cursor()
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

    cur.execute(
        """
        INSERT /*+ APPEND */ INTO IPRSCAN.AA_ALIGNMENT
        SELECT M.UPI, UPPER(TRANSLATE(DB.DBNAME, ' ', '_')), 
               M.METHOD_AC, M.SEQ_START, M.SEQ_END, M.SEQ_FEATURE
        FROM INTERPRO.CV_DATABASE DB
        INNER JOIN INTERPRO.IPRSCAN2DBCODE I2D 
            ON DB.DBCODE = I2D.DBCODE
        INNER JOIN IPRSCAN.MV_IPRSCAN M 
            ON I2D.IPRSCAN_SIG_LIB_REL_ID = M.ANALYSIS_ID
        WHERE DB.DBCODE IN ('Q', 'P', 'M')
        """
    )
    con.commit()

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


def create_aa_iprscan(uri: str):
    logger.info("creating AA_IPRSCAN")

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


def create_xref_condensed(uri: str):
    logger.info("creating XREF_CONDENSED")
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
            SELECT PROTEIN_AC, METHOD_AC, POS_FROM, POS_TO
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


def create_xref_summary(uri: str):
    logger.info("creating XREF_SUMMARY")
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

    logger.info("inserting signature matches")
    cur.execute(
        """
        INSERT /*+ APPEND */ INTO INTERPRO.XREF_SUMMARY
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
        """
    )
    con.commit()

    logger.info("inserting PANTHER subfamily matches")
    cur.execute(
        r"""
        INSERT /*+ APPEND */ INTO INTERPRO.XREF_SUMMARY
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
        """
    )
    con.commit()

    logger.info("inserting AntiFam and CATH-FunFam matches")
    cur.execute(
        """
        INSERT /*+ APPEND */ INTO INTERPRO.XREF_SUMMARY
        SELECT
            FM.DBCODE,
            FM.PROTEIN_AC,
            NULL,
            NULL,
            FM.METHOD_AC,
            ME.DESCRIPTION,
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


def _repr_domains_worker(matches_file: str,
                         inqueue: Queue,
                         outqueue: Queue,
                         domain_signatures: dict[str: str],
                         output: str):
    with open(output, "wt") as fh:
        for prot_acc, _, _, _, matches in iter_matches(matches_file, inqueue, outqueue):
            domains = []
            for signature_acc, models in matches.items():
                try:
                    dbcode = domain_signatures[signature_acc]
                except KeyError:
                    continue

                for _, locations in models.values():
                    for fragments_as_str in locations:
                        fragments = _parse_fragments(fragments_as_str)
                        pos_start = fragments[0]["start"]
                        pos_end = max(f["end"] for f in fragments)
                        domains.append({
                            "signature": signature_acc,
                            "start": pos_start,
                            "end": pos_end,
                            "frag": fragments_as_str,
                            "fragments": fragments,
                            "rank": REPR_DOM_DATABASES.index(dbcode)
                        })

            if domains:
                for domain in _select_repr_domains(domains):
                    fh.write(
                        f"{prot_acc}\t{domain['signature']}\t"
                        f"{domain['start']}\t{domain['end']}\t"
                        f"{domain['frag']}\n"
                    )


def export_repr_domains(ora_url: str, matches_file: str, output: str, emails: dict, processes: int = 8):
    logger.info("starting")

    con = oracledb.connect(ora_url)
    cur = con.cursor()
    params_dbcode = ",".join(":1" for _ in REPR_DOM_DATABASES)
    params_types = ",".join(":1" for _ in REPR_DOM_TYPES)
    cur.execute(
        f"""
        SELECT METHOD_AC, DBCODE
        FROM INTERPRO.METHOD
        WHERE DBCODE in ({params_dbcode})
          AND SIG_TYPE IN ({params_types})
        """,
        REPR_DOM_DATABASES + REPR_DOM_TYPES
    )
    domain_signatures = dict(cur.fetchall())
    cur.execute("SELECT VERSION FROM INTERPRO.DB_VERSION WHERE DBCODE = 'u'")
    release, = cur.fetchone()
    cur.close()
    con.close()

    num_workers = max(1, processes - 1)
    index = load_index(matches_file)
    tasks_per_worker = math.ceil(len(index) / num_workers)

    inqueue = Queue()
    outqueue = Queue()
    workers = []
    i = total = 0
    for _ in range(num_workers):
        fd, tmpfile = mkstemp()
        os.close(fd)
        p = Process(
            target=_repr_domains_worker,
            args=(matches_file, inqueue, outqueue, domain_signatures, tmpfile),
        )
        p.start()
        workers.append(tmpfile)

        # Enqueue all tasks for the worker, so they are continuous and proteins
        # will be sorted by accession
        tasks = []
        for _ in range(tasks_per_worker):
            try:
                offset, count = index[i]
            except IndexError:
                break
            else:
                tasks.append((offset, count))
                total += count
                i += 1

        inqueue.put(tasks)

    for _ in workers:
        inqueue.put(None)

    done = 0
    milestone = step = 10
    for _ in index:
        count = outqueue.get()
        done += count
        progress = math.floor(done / total * 100)
        if progress >= milestone:
            logger.info(f"{progress}%")
            milestone += step

    with open(output, "wt") as fh:
        for p, tmpfile in workers:
            p.join()

            with open(tmpfile, "rt") as fh2:
                while (block := fh2.read(1024)) != '':
                    fh.write(block)

    os.chmod(output, 0o664)

    email.send(
        info=emails,
        to=["aa_dev"],
        bcc=["sender"],
        subject="InterPro representative domains are ready",
        content=f"""\
Dear UniProt team,

The file containing the representative domains for UniProt {release} is \
available at the following path:
  {output}

Kind regards,
The InterPro Production Team
"""
    )
    logger.info("done")


def _select_repr_domains(domains: list[dict]) -> list[dict]:
    repr_domains = []

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
            repr_domains.append(domain)
    return repr_domains


def _calc_coverage(domain: dict) -> set[int]:
    residues = set()
    for f in domain["fragments"]:
        residues |= set(range(f["start"], f["end"] + 1))

    return residues


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


def _eval_overlap(dom_a: dict, dom_b: dict, threshold: float) -> bool:
    overlap = len(dom_a["residues"] & dom_b["residues"])
    return overlap and overlap / min(len(dom_a["residues"]),
                                     len(dom_b["residues"])) >= threshold


def _parse_fragments(fragments_as_str: str) -> list[dict]:
    fragments = []
    for frag in fragments_as_str.split(','):
        # Format: START-END-STATUS
        s, e, t = frag.split('-')
        fragments.append({
            "start": int(s),
            "end": int(e),
            "dc-status": t
        })

    return sorted(fragments, key=lambda x: (x["start"], x["end"]))
