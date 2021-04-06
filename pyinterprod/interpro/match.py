# -*- coding: utf-8 -*-

import os
import pickle
from typing import Dict, Sequence

import cx_Oracle

from pyinterprod import logger
from pyinterprod.utils import oracle
from .database import Database


FILE_ENTRY_PROT_COUNTS = "entries.prot.counts.pickle"
MATCH_PARTITIONS = {
    'B': "MATCH_DBCODE_B",  # SFLD
    'F': "MATCH_DBCODE_F",  # PRINTS
    'H': "MATCH_DBCODE_H",  # Pfam
    'J': "MATCH_DBCODE_J",  # CDD
    'M': "MATCH_DBCODE_M",  # PROSITE profiles
    'N': "MATCH_DBCODE_N",  # TIGRFAMs
    'P': "MATCH_DBCODE_P",  # PROSITE patterns
    'Q': "MATCH_DBCODE_Q",  # HAMAP
    'R': "MATCH_DBCODE_R",  # SMART
    'U': "MATCH_DBCODE_U",  # PIRSF
    'V': "MATCH_DBCODE_V",  # PANTHER
    'X': "MATCH_DBCODE_X",  # CATH-Gene3D
    'Y': "MATCH_DBCODE_Y",  # SUPERFAMILY
}
SITE_PARTITIONS = {
    # DB identifier -> (str:partition, bool:check against MATCH table)
    'B': ("SFLD", True),
    'J': ("CDD", True),
    'z': ("PIRSR", False)
}


def update_database_matches(url: str, databases: Sequence[Database]):
    con = cx_Oracle.connect(url)
    cur = con.cursor()

    for database in databases:
        logger.info(f"{database.name}")
        oracle.drop_table(cur, "INTERPRO.MATCH_NEW", purge=True)
        logger.debug(f"\tpopulating MATCH_MEW "
                     f"(ANALYSIS_ID: {database.analysis_id})")
        cur.execute(
            """
            CREATE TABLE INTERPRO.MATCH_NEW NOLOGGING
            AS SELECT * FROM INTERPRO.MATCH WHERE 1 = 0
            """
        )
        cur.execute(
            """
            INSERT /*+ APPEND */ INTO INTERPRO.MATCH_NEW
            SELECT
              X.AC, M.METHOD_AC, M.SEQ_START, M.SEQ_END, 'T',
              D.DBCODE, D.EVIDENCE,
              SYSDATE, SYSDATE, SYSDATE, 'INTERPRO',
              M.EVALUE, M.MODEL_AC, M.FRAGMENTS
            FROM IPRSCAN.MV_IPRSCAN M
            INNER JOIN UNIPARC.XREF X 
              ON M.UPI = X.UPI
            INNER JOIN INTERPRO.IPRSCAN2DBCODE D
              ON M.ANALYSIS_ID = D.IPRSCAN_SIG_LIB_REL_ID
            WHERE M.ANALYSIS_ID = :1
            AND M.SEQ_START != M.SEQ_END 
            AND X.DBID IN (2, 3)  -- Swiss-Prot or TrEMBL
            AND X.DELETED = 'N'            
            """, (database.analysis_id,)
        )
        con.commit()

        # Add constraints/indexes to be able to exchange partition
        logger.debug("\tcreating indexes and constraints")
        for col in ("DBCODE", "EVIDENCE", "STATUS", "METHOD_AC"):
            cur.execute(
                f"""
                CREATE INDEX MATCH_NEW${col[0]}
                ON INTERPRO.MATCH_NEW ({col}) 
                TABLESPACE INTERPRO_IND
                NOLOGGING
                """
            )

        cur.execute(
            """
            ALTER TABLE INTERPRO.MATCH_NEW
            ADD CONSTRAINT MATCH_NEW$CK1 
            CHECK ( POS_FROM >= 1 ) 
            """
        )
        cur.execute(
            """
            ALTER TABLE INTERPRO.MATCH_NEW
            ADD CONSTRAINT MATCH_NEW$CK2
            CHECK ( POS_TO - POS_FROM > 0 ) 
            """
        )
        cur.execute(
            """
            ALTER TABLE INTERPRO.MATCH_NEW
            ADD CONSTRAINT MATCH_NEW$CK3
            CHECK ( STATUS != 'N' OR (STATUS = 'N' AND DBCODE IN ('P', 'M', 'Q')) )
            """
        )
        cur.execute(
            """
            ALTER TABLE INTERPRO.MATCH_NEW
            ADD CONSTRAINT MATCH_NEW$PK
            PRIMARY KEY (PROTEIN_AC, METHOD_AC, POS_FROM, POS_TO)
            """
        )
        cur.execute(
            f"""
            ALTER TABLE INTERPRO.MATCH_NEW
            ADD CONSTRAINT MATCH_NEW$FK1 FOREIGN KEY (DBCODE) 
            REFERENCES INTERPRO.CV_DATABASE (DBCODE)
            """
        )
        cur.execute(
            """
            ALTER TABLE INTERPRO.MATCH_NEW
            ADD CONSTRAINT MATCH_NEW$FK2 FOREIGN KEY (EVIDENCE) 
            REFERENCES INTERPRO.CV_EVIDENCE (CODE)
            """
        )
        cur.execute(
            f"""
            ALTER TABLE INTERPRO.MATCH_NEW
            ADD CONSTRAINT MATCH_NEW$FK3 FOREIGN KEY (METHOD_AC)
            REFERENCES INTERPRO.METHOD (METHOD_AC)
            """
        )
        cur.execute(
            f"""
            ALTER TABLE INTERPRO.MATCH_NEW
            ADD CONSTRAINT MATCH_NEW$FK4 FOREIGN KEY (PROTEIN_AC) 
            REFERENCES INTERPRO.PROTEIN (PROTEIN_AC)
            """
        )
        cur.execute(
            f"""
            ALTER TABLE INTERPRO.MATCH_NEW
            ADD CONSTRAINT MATCH_NEW$FK5 FOREIGN KEY (STATUS) 
            REFERENCES INTERPRO.CV_STATUS (CODE)
            """
        )
        cur.execute(
            f"""
            ALTER TABLE INTERPRO.MATCH_NEW
            ADD CONSTRAINT MATCH_NEW$CK4 
            CHECK (PROTEIN_AC IS NOT NULL )
            """
        )
        cur.execute(
            f"""
            ALTER TABLE INTERPRO.MATCH_NEW
            ADD CONSTRAINT MATCH_NEW$CK5 
            CHECK (METHOD_AC IS NOT NULL )
            """
        )

        cur.execute("SELECT COUNT(*) FROM INTERPRO.MATCH_NEW")
        cnt, = cur.fetchone()
        if not cnt:
            raise RuntimeError(f"no rows inserted "
                               f"for analysis ID {database.analysis_id}")

        logger.debug(f"\texchanging partition")
        partition = MATCH_PARTITIONS[database.identifier]
        cur.execute(
            f"""
            ALTER TABLE INTERPRO.MATCH
            EXCHANGE PARTITION ({partition}) 
            WITH TABLE INTERPRO.MATCH_NEW
            """
        )
        oracle.drop_table(cur, "INTERPRO.MATCH_NEW", purge=True)

        logger.info("\tgathering statistics")
        oracle.gather_stats(cur, "INTERPRO", "MATCH", partition)

    for index in oracle.get_indexes(cur, "INTERPRO", "MATCH"):
        if index["unusable"]:
            logger.info(f"rebuilding index {index['name']}")
            oracle.rebuild_index(cur, index["name"])

    cur.close()
    con.close()

    logger.info("complete")


def update_database_site_matches(url: str, databases: Sequence[Database]):
    con = cx_Oracle.connect(url)
    cur = con.cursor()

    for database in databases:
        logger.info(database.name)

        """
        Same partition names in:
            - IPRSCAN.SITE
            - INTERPRO.SITE_MATCH
        """
        site_partition, ck_matches = SITE_PARTITIONS[database.identifier]

        oracle.drop_table(cur, "INTERPRO.SITE_MATCH_NEW", purge=True)
        cur.execute(
            """
            CREATE TABLE INTERPRO.SITE_MATCH_NEW NOLOGGING
            AS SELECT * FROM INTERPRO.SITE_MATCH WHERE 1 = 0
            """
        )

        logger.debug(f"\tinserting site matches")
        cur.execute(
            f"""
            INSERT /*+ APPEND */ INTO INTERPRO.SITE_MATCH_NEW
            SELECT
                X.AC, S.METHOD_AC, S.LOC_START, S.LOC_END, S.DESCRIPTION,
                S.RESIDUE, S.RESIDUE_START, S.RESIDUE_END, S.NUM_SITES, 
                D.DBCODE
            FROM IPRSCAN.SITE PARTITION ({site_partition}) S
            INNER JOIN UNIPARC.XREF X 
              ON S.UPI = X.UPI
            INNER JOIN INTERPRO.IPRSCAN2DBCODE D
              ON S.ANALYSIS_ID = D.IPRSCAN_SIG_LIB_REL_ID
            WHERE S.ANALYSIS_ID = :1
            AND X.DBID IN (2, 3)  -- Swiss-Prot or TrEMBL
            AND X.DELETED = 'N'
            """, (database.analysis_id,)
        )
        con.commit()

        logger.debug(f"\tindexing")
        cur.execute(
            """
            CREATE INDEX I_SITE_MATCH_NEW
            ON INTERPRO.SITE_MATCH_NEW (
                PROTEIN_AC, METHOD_AC, LOC_START, LOC_END
            )
            TABLESPACE INTERPRO_IND 
            NOLOGGING
            """
        )

        if ck_matches:
            logger.debug(f"\tchecking matches")
            match_partition = MATCH_PARTITIONS[database.identifier]
            cur.execute(
                f"""
                SELECT COUNT(*)
                FROM (
                    SELECT DISTINCT PROTEIN_AC, METHOD_AC, LOC_START, LOC_END
                    FROM INTERPRO.SITE_MATCH_NEW
                    MINUS
                    SELECT DISTINCT PROTEIN_AC, METHOD_AC, POS_FROM, POS_TO
                    FROM INTERPRO.MATCH PARTITION ({match_partition})
                )
                """
            )

            cnt, = cur.fetchone()
            if cnt:
                cur.close()
                con.close()
                raise RuntimeError(f"{database.name}: {cnt} matches "
                                   f"in SITE_MATCH_NEW that are not in MATCH")

        logger.debug(f"\tadding constraint")
        cur.execute(
            """
            ALTER TABLE INTERPRO.SITE_MATCH_NEW 
            ADD CONSTRAINT FK_SITE_MATCH_NEW 
            FOREIGN KEY (PROTEIN_AC) REFERENCES PROTEIN
            """
        )

        logger.debug(f"\texchanging partition")
        cur.execute(
            f"""
            ALTER TABLE INTERPRO.SITE_MATCH 
            EXCHANGE PARTITION ({site_partition})
            WITH TABLE INTERPRO.SITE_MATCH_NEW
            """
        )

        oracle.drop_table(cur, "INTERPRO.SITE_MATCH_NEW", purge=True)

    for index in oracle.get_indexes(cur, "INTERPRO", "SITE_MATCH"):
        if index["unusable"]:
            logger.info(f"rebuilding index {index['name']}")
            oracle.rebuild_index(cur, index["name"])

    cur.close()
    con.close()

    logger.info("complete")


def update_matches(url: str, data_dir: str):
    """
    Add protein matches for recently added/modified sequences

    :param url: Oracle connection string
    :param data_dir: output directory for data files
    """
    con = cx_Oracle.connect(url)
    _prepare_matches(con)
    _check_matches(con, os.path.join(data_dir, FILE_ENTRY_PROT_COUNTS))
    _insert_matches(con)
    con.close()


def update_feature_matches(url: str):
    """
    Add protein feature matches for recently added/modified sequences

    :param url: Oracle connection string
    """
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    logger.info("updating FEATURE_MATCH")
    cur.execute(
        """
        DELETE FROM INTERPRO.FEATURE_MATCH
        WHERE PROTEIN_AC IN (
          SELECT PROTEIN_AC
          FROM INTERPRO.PROTEIN_TO_SCAN
        )
        """
    )
    logger.info(f"{cur.rowcount} rows deleted")

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
        -- Include MobiDB-Lite, Phobius, SignalP (Euk, Gram+, Gram-), TMHMM, COILS
        WHERE D.DBCODE IN ('g', 'j', 'n', 's', 'v', 'q', 'x')
        """
    )
    logger.info(f"{cur.rowcount} rows inserted")
    con.commit()

    cur.close()
    con.close()


def update_variant_matches(url: str):
    """
    Recreate splice-variants table with the most recent data

    :param url: Oracle connection string
    """
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    logger.info("updating VARSPLIC_MASTER")
    oracle.truncate_table(cur, "INTERPRO.VARSPLIC_MASTER", reuse_storage=True)
    cur.execute(
        """
        INSERT INTO INTERPRO.VARSPLIC_MASTER
        SELECT 
          SUBSTR(X.AC, 1, INSTR(X.AC, '-') - 1),
          SUBSTR(X.AC, INSTR(X.AC, '-') + 1),
          P.CRC64,
          P.LEN
        FROM UNIPARC.XREF X
        INNER JOIN UNIPARC.PROTEIN P ON X.UPI = P.UPI
        WHERE X.DBID IN (24, 25) -- SWISSPROT_VARSPLIC, TREMBL_VARSPLIC
        AND X.DELETED = 'N'
        """
    )
    logger.info(f"{cur.rowcount} rows inserted")

    logger.info("updating VARSPLIC_MATCH")
    oracle.truncate_table(cur, "INTERPRO.VARSPLIC_MATCH", reuse_storage=True)
    cur.execute(
        """
        INSERT INTO INTERPRO.VARSPLIC_MATCH
        SELECT 
          X.AC, MV.METHOD_AC, MV.SEQ_START, MV.SEQ_END, 'T' AS STATUS, 
          I2D.DBCODE, I2D.EVIDENCE, SYSDATE, SYSDATE, SYSDATE, 
          'INTERPRO', MV.EVALUE, MV.MODEL_AC, MV.FRAGMENTS
        FROM UNIPARC.XREF X
        INNER JOIN IPRSCAN.MV_IPRSCAN MV
          ON X.UPI = MV.UPI
        INNER JOIN INTERPRO.IPRSCAN2DBCODE I2D 
          ON MV.ANALYSIS_ID = I2D.IPRSCAN_SIG_LIB_REL_ID
        INNER JOIN INTERPRO.METHOD M
          ON MV.METHOD_AC = M.METHOD_AC
        WHERE X.DBID IN (24, 25)  -- SWISSPROT_VARSPLIC, TREMBL_VARSPLIC
        AND X.DELETED = 'N'
        -- Exclude MobiDB-Lite, Phobius, SignalP (Euk, Gram+, Gram-), TMHMM, COILS
        AND I2D.DBCODE NOT IN ('g', 'j', 'n', 'q', 's', 'v', 'x')
        AND M.SKIP_FLAG = 'N'
        """
    )
    logger.info(f"{cur.rowcount} rows inserted")
    con.commit()

    cur.close()
    con.close()


def update_site_matches(url: str):
    """
    Add protein site matches for recently added/modified sequences

    :param url: Oracle connection string
    """
    con = cx_Oracle.connect(url)
    cur = con.cursor()

    logger.info("populating SITE_MATCH_NEW")
    oracle.drop_table(cur, "INTERPRO.SITE_MATCH_NEW", purge=True)

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

    logger.info("indexing")
    cur.execute(
        """
        CREATE INDEX I_SITE_NEW
        ON INTERPRO.SITE_MATCH_NEW (PROTEIN_AC, METHOD_AC, LOC_START, LOC_END)
        TABLESPACE INTERPRO_IND
        """
    )

    logger.info("gathering statistics")
    oracle.gather_stats(cur, "INTERPRO", "SITE_MATCH_NEW")

    logger.info("checking")
    queries = []
    for identifier, ck_matches in SITE_PARTITIONS:
        if ck_matches:
            queries.append(
                f"""
                SELECT DISTINCT PROTEIN_AC, METHOD_AC, POS_FROM, POS_TO
                FROM INTERPRO.MATCH PARTITION ({MATCH_PARTITIONS[identifier]})
                """
            )

    if queries:
        cur.execute(
            f"""
            SELECT COUNT(*)
            FROM (
                SELECT DISTINCT PROTEIN_AC, METHOD_AC, LOC_START, LOC_END
                FROM INTERPRO.SITE_MATCH_NEW
                MINUS (
                  {' UNION ALL '.join(queries)}
                )
            )
            """
        )

        cnt, = cur.fetchone()
        if cnt:
            cur.close()
            con.close()
            raise RuntimeError(f"{cnt} matches in SITE_MATCH_NEW "
                               f"that are not in MATCH")

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
    logger.info(f"{cur.rowcount} rows deleted")

    cur.execute(
        """
        INSERT INTO INTERPRO.SITE_MATCH
        SELECT * FROM INTERPRO.SITE_MATCH_NEW
        """
    )
    logger.info(f"{cur.rowcount} rows inserted")
    con.commit()

    oracle.drop_table(cur, "INTERPRO.SITE_MATCH_NEW", purge=True)

    cur.close()
    con.close()


def _prepare_matches(con: cx_Oracle.Connection):
    """
    Import protein matches in a staging table

    :param con: Oracle connection object
    """
    cur = con.cursor()

    logger.info("populating MATCH_NEW")
    oracle.drop_table(cur, "INTERPRO.MATCH_NEW", purge=True)

    cur.execute(
        """
        CREATE TABLE INTERPRO.MATCH_NEW NOLOGGING
        AS
        SELECT * 
        FROM INTERPRO.MATCH WHERE 1 = 0
        """
    )

    cur.execute(
        """
        INSERT /*+ APPEND */ INTO INTERPRO.MATCH_NEW
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
        -- Exclude MobiDB-Lite, Phobius, SignalP (Euk, Gram+, Gram-), TMHMM, COILS
        WHERE D.DBCODE NOT IN ('g', 'j', 'n', 's', 'v', 'q', 'x')
        AND M.SEQ_START != M.SEQ_END
        """
    )
    con.commit()

    logger.info("indexing")
    for col in ("DBCODE", "PROTEIN_AC"):
        cur.execute(
            f"""
            CREATE INDEX I_MATCH_NEW${col}
            ON INTERPRO.MATCH_NEW ({col})
            TABLESPACE INTERPRO_IND
            NOLOGGING
            """
        )

    # logger.info("gathering statistics")
    # oracle.gather_stats(cur, "INTERPRO", "MATCH_NEW")

    logger.info("deleting SUPERFAMILY duplicated matches")
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
    logger.info(f"{cur.rowcount} SUPERFAMILY matches deleted")
    con.commit()
    cur.close()


def _check_matches(con: cx_Oracle.Connection, output: str):
    """
    Check there are not errors in imported matches

    :param con: Oracle connection object
    :param output: output file to store the number of protein matches per entry
    """
    cur = con.cursor()

    # Matches outside of the protein
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
    cnt = 0
    for row in cur:
        logger.critical("out-of-bound: {}\t{}\t{}\t{}".format(*row))
        cnt += 1

    if cnt:
        cur.close()
        con.close()
        raise RuntimeError(f"{cnt} out-of-bound matches")

    # Matches with invalid start/end positions
    logger.info("checking invalid matches")
    cur.execute(
        """
        SELECT PROTEIN_AC, METHOD_AC, POS_FROM, POS_TO
        FROM INTERPRO.MATCH_NEW
        WHERE POS_FROM < 1 OR POS_FROM > POS_TO
        """
    )
    cnt = 0
    for row in cur:
        logger.critical("invalid: {}\t{}\t{}\t{}".format(*row))
        cnt += 1

    if cnt:
        cur.close()
        con.close()
        raise RuntimeError(f"{cnt} invalid matches")

    with open(output, "wb") as fh:
        pickle.dump(_get_entries_proteins_count(cur), fh)

    cur.close()


def _insert_matches(con: cx_Oracle.Connection):
    """
    Update the MATCH table with data from the staging table

    :param con: Oracle connection object
    """
    cur = con.cursor()

    logger.info("updating MATCH")
    cur.execute(
        """
        DELETE FROM INTERPRO.MATCH
        WHERE PROTEIN_AC IN (
          SELECT PROTEIN_AC
          FROM INTERPRO.PROTEIN_TO_SCAN
        )
        """
    )
    logger.info(f"{cur.rowcount} rows deleted")

    cur.execute(
        """
        INSERT INTO INTERPRO.MATCH
        SELECT * FROM INTERPRO.MATCH_NEW
        """
    )
    logger.info(f"{cur.rowcount} rows inserted")
    con.commit()

    oracle.drop_table(cur, "INTERPRO.MATCH_NEW", purge=True)
    cur.close()


def track_entry_changes(cur: cx_Oracle.Cursor, data_dir: str) -> list:
    """
    Find entries with significant protein count changes

    :param cur: Oracle cursor object
    :param data_dir: directory containing the file for protein counts
                     before the update
    :return: list of entries with a significant changes in proteins count
    """

    with open(os.path.join(data_dir, FILE_ENTRY_PROT_COUNTS), "rb") as fh:
        old_counts = pickle.load(fh)

    new_counts = _get_entries_proteins_count(cur)
    changes = []
    for acc in sorted(old_counts):
        entry_old_counts = old_counts[acc]
        entry_new_counts = new_counts.pop(acc, {})

        # Total number of proteins matched
        entry_old_total = sum(entry_old_counts.values())
        entry_new_total = sum(entry_new_counts.values())

        change = (entry_new_total - entry_old_total) / entry_old_total

        # If the entry does not have any matches anymore,
        # we want to report it
        if entry_new_total != 0 and abs(change) < 0.5:
            continue

        entry_superkingdoms = {}
        for superkingdom, old_cnt in entry_old_counts.items():
            new_cnt = entry_new_counts.pop(superkingdom, 0)
            entry_superkingdoms[superkingdom] = (old_cnt, new_cnt)

        # superkingdoms with proteins only matched in new UniProt release
        for superkingdom, new_cnt in entry_new_counts.items():
            entry_superkingdoms[superkingdom] = (0, new_cnt)

        changes.append((
            acc,
            entry_old_total,
            entry_new_total,
            change,
            entry_superkingdoms
        ))

    return changes


def _get_taxon2superkingdom(cur: cx_Oracle.Cursor) -> Dict[int, str]:
    # Load all taxa
    cur.execute(
        """
        SELECT TAX_ID, SCIENTIFIC_NAME, RANK, PARENT_ID
        FROM INTERPRO.ETAXI
        """
    )
    taxa = {}
    for tax_id, name, rank, parent_id in cur:
        if tax_id in (1, 131567):
            """
            Skip root and meta-superkingdom (131567) which contains:
                * Bacteria (2)
                * Archaea (2157)
                * Eukaryota (2759)
            """
            continue
        elif parent_id in (1, 131567):
            rank = "superkingdom"
            parent_id = None

        taxa[tax_id] = (name, rank, parent_id)

    # For each taxon, find its root (i.e. superkingdom)
    taxon2superkingdom = {}
    for tax_id in taxa:
        name, rank, parent_id = taxa[tax_id]

        while parent_id is not None:
            name, rank, parent_id = taxa[parent_id]
            if rank == "superkingdom":
                taxon2superkingdom[tax_id] = name
                break

    return taxon2superkingdom


def _get_entries_proteins_count(cur: cx_Oracle.Cursor) -> Dict[str, Dict[str, int]]:
    """
    Return the number of protein matched by each InterPro entry.
    Only complete sequences are considered.

    :param cur: Oracle cursor object
    :return: dictionary
    """
    taxon2superkingdom = _get_taxon2superkingdom(cur)
    cur.execute(
        """
        SELECT EM.ENTRY_AC, P.TAX_ID, COUNT(DISTINCT P.PROTEIN_AC)
        FROM INTERPRO.PROTEIN P
        INNER JOIN INTERPRO.MATCH M ON P.PROTEIN_AC = M.PROTEIN_AC
        INNER JOIN INTERPRO.ENTRY2METHOD EM ON EM.METHOD_AC = M.METHOD_AC
        WHERE P.FRAGMENT = 'N'
        GROUP BY EM.ENTRY_AC, P.TAX_ID
        """
    )
    counts = {}
    for entry_acc, tax_id, n_proteins in cur:
        try:
            e = counts[entry_acc]
        except KeyError:
            e = counts[entry_acc] = {}

        superkingdom = taxon2superkingdom[tax_id]
        try:
            e[superkingdom] += n_proteins
        except KeyError:
            e[superkingdom] = n_proteins

    return dict(counts)


def get_sig_proteins_count(cur: cx_Oracle.Cursor, dbid: str) -> Dict[str, Dict[str, int]]:
    """
    Return the number of protein matches by each member database signature.
    Only complete sequences are considered

    :param cur: Oracle cursor object
    :param dbid: member database identifier
    :return: dictionary
    """
    partition = MATCH_PARTITIONS[dbid]
    taxon2superkingdom = _get_taxon2superkingdom(cur)
    cur.execute(
        f"""
        SELECT M.METHOD_AC, P.TAX_ID, COUNT(DISTINCT P.PROTEIN_AC)
        FROM INTERPRO.MATCH PARTITION ({partition}) M 
        INNER JOIN INTERPRO.PROTEIN P
            ON P.PROTEIN_AC = M.PROTEIN_AC
        WHERE P.FRAGMENT = 'N'
        GROUP BY M.METHOD_AC, P.TAX_ID
        """
    )
    counts = {}
    for sig_acc, tax_id, n_proteins in cur:
        try:
            sig = counts[sig_acc]
        except KeyError:
            sig = counts[sig_acc] = {}

        superkingdom = taxon2superkingdom[tax_id]
        try:
            sig[superkingdom] += n_proteins
        except KeyError:
            sig[superkingdom] = n_proteins

    return counts


# def _get_databases_matches_count(cur: cx_Oracle.Cursor) -> Dict[str, int]:
#     """
#     Return the number of matches per member database
#     :param cur: Oracle cursor object
#     :return: dictionary
#     """
#     cur.execute(
#         """
#         SELECT DBCODE, COUNT(*)
#         FROM INTERPRO.MATCH
#         GROUP BY DBCODE
#         """
#     )
#     return dict(cur.fetchall())
