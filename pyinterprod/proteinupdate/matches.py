import json
import os

import cx_Oracle

from .. import logger, orautils


def _get_match_counts(cur: cx_Oracle.Cursor) -> (str, int):
    logger.info("counting entries protein counts")
    cur.execute(
        """
        SELECT E.ENTRY_AC, COUNT(DISTINCT M.PROTEIN_AC)
        FROM INTERPRO.ENTRY2METHOD E
        INNER JOIN INTERPRO.MATCH M
          ON E.METHOD_AC = M.METHOD_AC
        GROUP BY E.ENTRY_AC
        """
    )
    entries = dict(cur.fetchall())

    logger.info("counting databases match counts")
    cur.execute(
        """
        SELECT DBCODE, COUNT(*)
        FROM INTERPRO.MATCH
        GROUP BY DBCODE
        """
    )
    databases = dict(cur.fetchall())

    return entries, databases


def prepare_matches(user: str, dsn: str):
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()

    logger.info("adding new matches")
    # cur.execute("TRUNCATE TABLE INTERPRO.MATCH_NEW")
    orautils.drop_table(cur, "INTERPRO", "MATCH_NEW")
    cur.execute(
        """
        CREATE TABLE INTERPRO.MATCH_NEW NOLOGGING
        AS
        SELECT * FROM INTERPRO.MATCH WHERE 1 = 0
        """
    )

    cur.execute(
        """
        INSERT /*+ APPEND */ INTO INTERPRO.MATCH_NEW (
          PROTEIN_AC, METHOD_AC, POS_FROM, POS_TO, STATUS,
          DBCODE, EVIDENCE,
          SEQ_DATE, MATCH_DATE, TIMESTAMP, USERSTAMP,
          SCORE, MODEL_AC, FRAGMENTS
        )
        SELECT
          P.PROTEIN_AC, M.METHOD_AC, M.SEQ_START, M.SEQ_END, 'T',
          D.DBCODE, D.EVIDENCE,
          SYSDATE, SYSDATE, SYSDATE, 'INTERPRO',
          M.EVALUE, M.MODEL_AC, M.FRAGMENTS
        FROM INTERPRO.PROTEIN_TO_SCAN P
        INNER JOIN IPRSCAN.MV_IPRSCAN M
          ON P.UPI = M.UPI
        INNER JOIN INTERPRO.IPRSCAN2DBCODE D
          ON M.ANALYSIS_ID = D.IPRSCAN_SIG_LIB_REL_ID
        WHERE D.DBCODE NOT IN ('g', 'j', 'n', 'q', 's', 'v', 'x')
        AND M.SEQ_START != M.SEQ_END
        """
    )
    con.commit()
    orautils.gather_stats(cur, "INTERPRO", "MATCH_NEW")

    logger.info("building indices")
    for col in ("DBCODE", "PROTEIN_AC"):
        cur.execute(
            """
            CREATE INDEX I_MATCH_NEW${0}
            ON INTERPRO.MATCH_NEW ({0}) NOLOGGING
            """.format(col)
        )

    logger.info("SUPERFAMILY: deleting duplicated matches")
    cur.execute(
        """
        DELETE FROM INTERPRO.MATCH_NEW M1
        WHERE EXISTS(
          SELECT 1
          FROM INTERPRO.MATCH_NEW M2
          WHERE M2.DBCODE = 'Y'
          AND M1.PROTEIN_AC = M2.PROTEIN_AC
          AND M1.METHOD_AC = M2.METHOD_AC
          AND M1.POS_FROM = M2.POS_FROM
          AND M1.POS_TO = M2.POS_TO
          AND M1.SCORE > M2.SCORE
        )
        """
    )
    logger.info("SUPERFAMILY: {} rows deleted".format(cur.rowcount))
    con.commit()

    cur.close()
    con.close()


def check_matches(user: str, dsn: str, outdir: str):
    os.makedirs(outdir, exist_ok=True)

    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()

    # Match outside of the protein
    logger.info("checking out-of-bound matches")
    cur.execute(
        """
        SELECT M.PROTEIN_AC, M.METHOD_AC, M.POS_TO, P.LEN
        FROM INTERPRO.MATCH_NEW M
        INNER JOIN INTERPRO.PROTEIN P
          ON M.PROTEIN_AC = P.PROTEIN_AC
        WHERE M.POS_TO > P.LEN
        """
    )
    n = 0
    for row in cur:
        logger.critical("out-of-bound: {}\t{}\t{}\t{}".format(*row))
        n += 1

    if n:
        cur.close()
        con.close()
        raise RuntimeError("{} out-of-bound matches".format(n))

    # Match with invalid start/end positions
    logger.info("checking invalid matches")
    cur.execute(
        """
        SELECT PROTEIN_AC, METHOD_AC, POS_FROM, POS_TO
        FROM INTERPRO.MATCH_NEW
        WHERE POS_FROM < 1 OR POS_FROM > POS_TO
        """
    )
    n = 0
    for row in cur:
        logger.critical("invalid: {}\t{}\t{}\t{}".format(*row))
        n += 1

    if n:
        cur.close()
        con.close()
        raise RuntimeError("{} invalid matches".format(n))

    entries, databases = _get_match_counts(cur)
    cur.close()
    con.close()

    with open(os.path.join(outdir, "counts.json"), "wt") as fh:
        json.dump(dict(entries=entries, databases=databases), fh)


def update_matches(user: str, dsn: str, drop_indices: bool=False):
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()

    logger.info("updating MATCH")
    if drop_indices:
        enforced = []
        for c in orautils.get_constraints(cur, "INTERPRO", "MATCH"):
            if c["index_name"]:
                # This constraint enforces an index, so it cannot be dropped
                enforced.append(c["index_name"])

        to_recreate = []
        for index in orautils.get_indices(cur, "INTERPRO", "MATCH"):
            if index["name"] not in enforced:
                orautils.drop_index(cur, index["owner"], index["name"])
                to_recreate.append(index)
    else:
        to_recreate = []

    cur.execute(
        """
        DELETE FROM INTERPRO.MATCH
        WHERE PROTEIN_AC IN (
          SELECT PROTEIN_AC
          FROM INTERPRO.PROTEIN_TO_SCAN
        )
        """
    )
    logger.info("{} rows deleted".format(cur.rowcount))
    con.commit()

    cur.execute(
        """
        INSERT INTO INTERPRO.MATCH
        SELECT * FROM INTERPRO.MATCH_NEW
        """
    )
    logger.info("{} rows inserted".format(cur.rowcount))
    con.commit()

    for index in to_recreate:
        orautils.recreate_index(cur, index)

    cur.close()
    con.close()


def track_count_changes(user: str, dsn: str, outdir: str):
    with open(os.path.join(outdir, "counts.json"), "rt") as fh:
        prev = json.load(fh)

    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    entries, databases = _get_match_counts(cur)
    cur.close()
    con.close()

    changes = []
    for entry_acc, new_count in entries.items():
        prev_count = prev["entries"].pop(entry_acc, 0)

        try:
            change = (new_count - prev_count) / prev_count * 100
        except ZeroDivisionError:
            changes.append((entry_acc, prev_count, new_count, "N/A"))
        else:
            if abs(change) >= 50:
                changes.append((entry_acc, prev_count, new_count, change))

    with open(os.path.join(outdir, "entries_changes.tsv"), "wt") as fh:
        def _sort_entries(e):
            return 1 if isinstance(e[3], str) else 0, e[3], e[0]

        fh.write("# Accession\tPrevious protein count\t"
                 "New protein count\tChange (%)\n")

        for ac, pc, nc, c in sorted(changes, key=_sort_entries):
            if isinstance(c, str):
                fh.write("{}\t{}\t{}\t{}\n".format(ac, pc, nc, c))
            else:
                fh.write("{}\t{}\t{}\t{:.0f}\n".format(ac, pc, nc, c))

    changes = []
    for dbcode, new_count in databases.items():
        prev_count = prev["databases"].pop(dbcode, 0)
        changes.append((dbcode, prev_count, new_count))

    for dbcode, prev_count in prev["databases"].items():
        changes.append((dbcode, prev_count, 0))

    with open(os.path.join(outdir, "databases_changes.tsv"), "wt") as fh:
        fh.write("# Code\tPrevious match count\tNew match count\n")

        for dbcode, pc, nc in sorted(changes):
            fh.write("{}\t{}\t{}\n".format(dbcode, pc, nc))


def update_feature_matches(user: str, dsn: str, drop_indices: bool=False):
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()

    logger.info("updating FEATURE_MATCH")
    if drop_indices:
        enforced = []
        for c in orautils.get_constraints(cur, "INTERPRO", "FEATURE_MATCH"):
            if c["index_name"]:
                # This constraint enforces an index, so it cannot be dropped
                enforced.append(c["index_name"])

        to_recreate = []
        for index in orautils.get_indices(cur, "INTERPRO", "FEATURE_MATCH"):
            if index["name"] not in enforced:
                orautils.drop_index(cur, index["owner"], index["name"])
                to_recreate.append(index)
    else:
        to_recreate = []

    cur.execute(
        """
        DELETE FROM INTERPRO.FEATURE_MATCH
        WHERE PROTEIN_AC IN (
          SELECT PROTEIN_AC
          FROM INTERPRO.PROTEIN_TO_SCAN
        )
        """
    )
    logger.info("{} rows deleted".format(cur.rowcount))
    con.commit()

    cur.execute(
        """
        INSERT INTO INTERPRO.FEATURE_MATCH
        SELECT
          P.PROTEIN_AC, M.METHOD_AC, M.SEQ_FEATURE, M.SEQ_START, M.SEQ_END,
          D.DBCODE
        FROM INTERPRO.PROTEIN_TO_SCAN P
        INNER JOIN IPRSCAN.MV_IPRSCAN M
          ON P.UPI = M.UPI
        INNER JOIN INTERPRO.IPRSCAN2DBCODE D
          ON M.ANALYSIS_ID = D.IPRSCAN_SIG_LIB_REL_ID
        WHERE D.DBCODE IN ('g', 'j', 'n', 'q', 's', 'v', 'x')
        """
    )
    logger.info("{} rows inserted".format(cur.rowcount))
    con.commit()

    for index in to_recreate:
        orautils.recreate_index(cur, index)

    cur.close()
    con.close()


def update_alt_splicing_matches(user: str, dsn: str):
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    for t in ("VARSPLIC_MASTER", "VARSPLIC_MATCH", "VARSPLIC_NEW"):
        orautils.drop_table(cur, "INTERPRO", t)

    logger.info("building table")
    cur.execute(
        """
        CREATE TABLE INTERPRO.VARSPLIC (
            PROTEIN_AC VARCHAR2(15) NOT NULL,
            VARIANT NUMBER(3) NOT NULL,
            LEN NUMBER(5) NOT NULL,
            IS_CANONICAL CHAR(1) NOT NULL,
            METHOD_AC VARCHAR2(25) NOT NULL,
            DBCODE CHAR(1) NOT NULL,
            POS_FROM NUMBER NOT NULL,
            POS_TO NUMBER NOT NULL,
            FRAGMENTS VARCHAR2(400),
            MODEL_AC VARCHAR2(255)
        ) NOLOGGING
        AS
        SELECT
          XV.PROTEIN_AC,
          XV.VARIANT,
          P.LEN,
          CASE WHEN XV.UPI = XP.UPI THEN 'Y' ELSE 'N' END,
          M.METHOD_AC,
          I2D.DBCODE,
          M.SEQ_START,
          M.SEQ_END,
          M.FRAGMENTS,
          M.MODEL_AC
        FROM (
            SELECT
                SUBSTR(AC, 1, INSTR(AC, '-') - 1) AS PROTEIN_AC,
                SUBSTR(AC, INSTR(AC, '-') + 1) AS VARIANT,
                UPI
            FROM UNIPARC.XREF
            WHERE DBID IN (24, 25)
            AND DELETED = 'N'
        ) XV
        INNER JOIN UNIPARC.PROTEIN P
            ON XV.UPI = P.UPI
        INNER JOIN UNIPARC.XREF XP
            ON XV.PROTEIN_AC = XP.AC
            AND XP.DBID IN (2, 3)
            AND XP.DELETED = 'N'
        INNER JOIN IPRSCAN.MV_IPRSCAN M
            ON XV.UPI = M.UPI
        INNER JOIN INTERPRO.IPRSCAN2DBCODE I2D
            ON M.ANALYSIS_ID = I2D.IPRSCAN_SIG_LIB_REL_ID
        WHERE I2D.DBCODE NOT IN ('g', 'j', 'n', 'q', 's', 'v', 'x')
        """
    )
    orautils.grant(cur, "INTERPRO", "VARSPLIC", "SELECT", "INTERPRO_SELECT")

    logger.info("analyzing table")
    orautils.gather_stats(cur, "INTERPRO", "VARSPLIC")

    logger.info("indexing table")
    cur.execute(
        """
        CREATE INDEX I_VARSPLIC
        ON INTERPRO.VARSPLIC (PROTEIN_AC, VARIANT) NOLOGGING
        """
    )

    cur.close()
    con.close()


def update_site_matches(user: str, dsn: str, drop_indices: bool=False):
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()

    logger.info("creating SITE_MATCH_NEW")
    orautils.drop_table(cur, "INTERPRO", "SITE_MATCH_NEW")
    cur.execute(
        """
        CREATE TABLE INTERPRO.SITE_MATCH_NEW
        AS
        SELECT
            P.PROTEIN_AC, S.METHOD_AC, S.LOC_START, S.LOC_END, S.DESCRIPTION,
            S.RESIDUE, S.RESIDUE_START, S.RESIDUE_END, S.NUM_SITES, D.DBCODE
        FROM INTERPRO.PROTEIN_TO_SCAN P
        INNER JOIN IPRSCAN.SITE S
          ON P.UPI = S.UPI
        INNER JOIN INTERPRO.IPRSCAN2DBCODE D
          ON S.ANALYSIS_ID = D.IPRSCAN_SIG_LIB_REL_ID
        """
    )

    logger.info("indexing SITE_MATCH_NEW")
    cur.execute(
        """
        CREATE INDEX I_SITE_NEW 
        ON INTERPRO.SITE_MATCH_NEW 
        (PROTEIN_AC, METHOD_AC, LOC_START, LOC_END) 
        NOLOGGING 
        """
    )

    logger.info("checking SITE_MATCH_NEW")
    cur.execute(
        """
        SELECT COUNT(*)
        FROM (
            SELECT DISTINCT PROTEIN_AC, METHOD_AC, LOC_START, LOC_END
            FROM INTERPRO.SITE_MATCH_NEW 
            MINUS (
              SELECT DISTINCT PROTEIN_AC, METHOD_AC, POS_FROM, POS_TO 
              FROM INTERPRO.MATCH PARTITION (MATCH_DBCODE_J)
              UNION ALL
              SELECT DISTINCT PROTEIN_AC, METHOD_AC, POS_FROM, POS_TO
              FROM INTERPRO.MATCH PARTITION (MATCH_DBCODE_B)
            )        
        )
        """
    )
    num_rows = cur.fetchone()[0]

    if num_rows:
        cur.close()
        con.close()
        raise RuntimeError("{} matches in SITE_MATCH_NEW "
                           "that are not in MATCH".format(num_rows))

    if drop_indices:
        enforced = []
        for c in orautils.get_constraints(cur, "INTERPRO", "SITE_MATCH"):
            if c["index_name"]:
                # This constraint enforces an index, so it cannot be dropped
                enforced.append(c["index_name"])

        to_recreate = []
        for index in orautils.get_indices(cur, "INTERPRO", "SITE_MATCH"):
            if index["name"] not in enforced:
                orautils.drop_index(cur, index["owner"], index["name"])
                to_recreate.append(index)
    else:
        to_recreate = []

    logger.info("updating SITE_MATCH")
    cur.execute(
        """
        DELETE FROM INTERPRO.SITE_MATCH
        WHERE PROTEIN_AC IN (
          SELECT PROTEIN_AC
          FROM INTERPRO.PROTEIN_TO_SCAN
        )
        """
    )
    logger.info("{} rows deleted".format(cur.rowcount))
    con.commit()

    cur.execute(
        """
        INSERT INTO INTERPRO.SITE_MATCH
        SELECT * FROM INTERPRO.SITE_MATCH_NEW
        """
    )
    logger.info("{} rows inserted".format(cur.rowcount))
    con.commit()

    for index in to_recreate:
        orautils.recreate_index(cur, index)

    cur.close()
    con.close()
