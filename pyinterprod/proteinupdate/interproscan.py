# -*- coding: utf-8 -*-

import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Collection, List, Optional

import cx_Oracle

from .. import logger, orautils


def _get_max_upi(cur: cx_Oracle.Cursor, analysis_id: int) -> Optional[str]:
    # TODO: decide *which* method is the most accurate
    cur.execute(
        """
        SELECT MAX(JOB_END)
        FROM IPRSCAN.IPM_PERSISTED_JOBS@ISPRO
        WHERE ANALYSIS_ID = :1
        AND PERSISTED > 0
        """, (analysis_id, )
    )
    row = cur.fetchone()
    return row[0] if row else None

    upis = []
    cur.execute(
        """
        SELECT MAX(HWM_SUBMITTED)
        FROM IPRSCAN.HWM@ISPRO
        WHERE ANALYSIS_ID = :1
        """, (analysis_id,)
    )
    row = cur.fetchone()
    if row:
        upis.append(row[0])

    cur.execute(
        """
        SELECT MAX(JOB_END)
        FROM IPRSCAN.IPM_RUNNING_JOBS@ISPRO
        WHERE ANALYSIS_ID = :1
        """, (analysis_id,)
    )
    row = cur.fetchone()
    if row:
        upis.append(row[0])

    cur.execute(
        """
        SELECT MIN(JOB_END)
        FROM (
            SELECT MIN(JOB_END) JOB_END
            FROM IPRSCAN.IPM_COMPLETED_JOBS@ISPRO
            WHERE ANALYSIS_ID = :analysisid
            AND PERSISTED = 0
            UNION ALL
            SELECT MAX(JOB_END) JOB_END
            FROM IPRSCAN.IPM_COMPLETED_JOBS@ISPRO
            WHERE ANALYSIS_ID = :analysisid
            AND PERSISTED = 1
        )
        """, dict(analysisid=analysis_id)
    )
    row = cur.fetchone()
    if row:
        upis.append(row[0])

    cur.execute(
        """
        SELECT MIN(JOB_END)
        FROM (
            SELECT MIN(JOB_END) JOB_END
            FROM IPRSCAN.IPM_PERSISTED_JOBS@ISPRO
            WHERE ANALYSIS_ID = :analysisid
            AND PERSISTED = 0
            UNION ALL
            SELECT MAX(JOB_END) JOB_END
            FROM IPRSCAN.IPM_PERSISTED_JOBS@ISPRO
            WHERE ANALYSIS_ID = :analysisid
            AND PERSISTED = 1
        )
        """, dict(analysisid=analysis_id)
    )
    row = cur.fetchone()
    if row:
        upis.append(row[0])

    try:
        return max([upi for upi in upis if upi is not None])
    except ValueError:
        return None


def _get_ispro_upis(url: str) -> List[dict]:
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.execute(
        """
        SELECT *
        FROM (
          SELECT
            A.ANALYSIS_ID,
            A.ANALYSIS_NAME,
            A.ANALYSIS_TYPE,
            B.MATCH_TABLE,
            B.SITE_TABLE,
            ROW_NUMBER() OVER (
              PARTITION
              BY A.ANALYSIS_MATCH_TABLE_ID
              ORDER BY A.ANALYSIS_ID DESC
            ) CNT
          FROM IPM_ANALYSIS@ISPRO A
            INNER JOIN IPM_ANALYSIS_MATCH_TABLE@ISPRO B
              ON A.ANALYSIS_MATCH_TABLE_ID = B.ID
          WHERE A.ACTIVE = 1
        )
        ORDER BY CNT
        """
    )

    match_tables = set()
    analyses = []
    for row in cur:
        """
        SignalP has EUK/Gram+/Gram- matches in the same tables so we need
        several times the same time for different analyses.

        But we shouldn't accept any other duplicates
        (ACTIVE might not always be up-to-date)
        """
        match_table = row[3]
        if match_table in match_tables and match_table != "ipm_signalp_match":
            continue

        match_tables.add(match_table)
        analyses.append({
            "id": row[0],
            "full_name": row[1],
            "name": row[2],
            "match_table": match_table,
            "site_table": row[4] if row[4] else "",
            "upi": None
        })

    for e in analyses:
        e["upi"] = _get_max_upi(cur, e["id"])

    cur.close()
    con.close()

    return analyses


def _check_ispro(url: str, max_attempts: int=1, secs: int=3600,
                 exclude: Optional[Collection[str]]=None) -> List[dict]:
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.execute("SELECT MAX(UPI) FROM UNIPARC.PROTEIN@UAREAD")
    max_upi_read = cur.fetchone()[0]
    cur.close()
    con.close()

    num_attempts = 0
    while True:
        analyses = _get_ispro_upis(url)
        not_ready = 0
        for e in analyses:
            if exclude and e["name"] in exclude:
                status = "ready (excluded)"
            elif e["upi"] and e["upi"] >= max_upi_read:
                status = "ready"
            else:
                status = "not ready"
                not_ready += 1

            logger.debug("{id:<5}{full_name:<30}{match_table:<30}"
                         "{site_table:<30}{upi:<20}"
                         "{status}".format(**e, status=status))

        num_attempts += 1
        if not_ready == 0:
            break
        elif num_attempts == max_attempts:
            raise RuntimeError("ISPRO tables not ready "
                               "after {} attempts".format(num_attempts))
        else:
            logger.info("{} analyses not ready".format(not_ready))
            time.sleep(secs)

    return analyses


def _import_table(url: str, owner: str, name: str):
    index = "I_" + name
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.execute(
        """
        SELECT MAX(UPI)
        FROM {}.{}@ISPRO
        """.format(owner, name)
    )
    src_upi = cur.fetchone()[0]
    cur.execute(
        """
        SELECT MAX(UPI)
        FROM {}.{}
        """.format(owner, name)
    )
    dst_upi = cur.fetchone()[0]

    if src_upi != dst_upi:
        orautils.drop_mview(cur, owner, name)
        orautils.drop_table(cur, owner, name)
        orautils.drop_index(cur, owner, index)
        cur.execute(
            """
            CREATE TABLE {0}.{1} NOLOGGING
            AS
            SELECT *
            FROM {0}.{1}@ISPRO
            """.format(owner, name)
        )
        orautils.gather_stats(cur, owner, name)
        cur.execute(
            """
            CREATE INDEX {0}.{1}
            ON {0}.{2}(ANALYSIS_ID)
            NOLOGGING
            """.format(owner, index, name)
        )
    cur.close()
    con.close()


def _init_ipm_table(cur: cx_Oracle.Cursor, name: str):
    orautils.drop_table(cur, "IPRSCAN", name)
    cur.execute(
        """
        CREATE TABLE IPRSCAN.{} NOLOGGING
        AS
        SELECT *
        FROM IPRSCAN.MV_IPRSCAN
        WHERE 1=0
        """.format(name)
    )


def _update_cdd_matches(url: str, analysis_id: int):
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    _init_ipm_table(cur, "IPM_CDD_MATCH_TMP")
    cur.execute(
        """
        INSERT /*+ APPEND */ INTO IPRSCAN.IPM_CDD_MATCH_TMP
        SELECT
          ANALYSIS_ID, UPI, METHOD_AC, RELNO_MAJOR, RELNO_MINOR, SEQ_START,
          SEQ_END, 0, 0, 0, NULL, SEQSCORE, SEQSCORE, SEQEVALUE, SEQEVALUE,
          0, 0, MODEL_AC, NULL, FRAGMENTS
        FROM IPRSCAN.IPM_CDD_MATCH
        WHERE ANALYSIS_ID = :1
        """, (analysis_id,)
    )
    con.commit()
    orautils.exchange_partition(cur, "IPRSCAN", "IPM_CDD_MATCH_TMP",
                                "MV_IPRSCAN", "CDD")
    orautils.drop_table(cur, "IPRSCAN", "IPM_CDD_MATCH_TMP")
    cur.close()
    con.close()


def _update_cdd_sites(url: str, analysis_id: int):
    con = cx_Oracle.connect(url)
    cur = con.cursor()

    cur.execute(
        """
        CREATE TABLE IPRSCAN.IPM_CDD_SITE_TMP
        NOLOGGING
        AS
        SELECT *
        FROM IPRSCAN.IPM_CDD_SITE
        WHERE ANALYSIS_ID = :1
        """, (analysis_id,)
    )
    con.commit()
    orautils.exchange_partition(cur, "IPRSCAN", "IPM_CDD_MATCH_TMP",
                                "SITE", "CDD")
    orautils.drop_table(cur, "IPRSCAN", "IPM_CDD_SITE_TMP")
    cur.close()
    con.close()


def _update_coils(url: str, analysis_id: int):
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    _init_ipm_table(cur, "IPM_COILS_MATCH_TMP")
    cur.execute(
        """
        INSERT /*+ APPEND */ INTO IPRSCAN.IPM_COILS_MATCH_TMP
        SELECT
          ANALYSIS_ID, UPI, METHOD_AC, RELNO_MAJOR, RELNO_MINOR, SEQ_START,
          SEQ_END, 0, 0, 0, NULL, 0, 0, 0, 0, 0, 0, MODEL_AC, NULL, FRAGMENTS
        FROM IPRSCAN.IPM_COILS_MATCH
        WHERE ANALYSIS_ID = :1
        """, (analysis_id,)
    )
    con.commit()
    orautils.exchange_partition(cur, "IPRSCAN", "IPM_COILS_MATCH_TMP",
                                "MV_IPRSCAN", "COILS")
    orautils.drop_table(cur, "IPRSCAN", "IPM_COILS_MATCH_TMP")
    cur.close()
    con.close()


def _update_gene3d(url: str, analysis_id: int):
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    _init_ipm_table(cur, "IPM_GENE3D_MATCH_TMP")
    cur.execute(
        """
        INSERT /*+ APPEND */ INTO IPRSCAN.IPM_GENE3D_MATCH_TMP
        SELECT
          ANALYSIS_ID, UPI, METHOD_AC, RELNO_MAJOR, RELNO_MINOR, SEQ_START,
          SEQ_END, HMM_START, HMM_END, HMM_LENGTH, HMM_BOUNDS, SCORE,
          SEQSCORE, EVALUE, SEQEVALUE, ENV_START, ENV_END, MODEL_AC,
          NULL, FRAGMENTS
        FROM IPRSCAN.IPM_GENE3D_MATCH
        WHERE ANALYSIS_ID = :1
        """, (analysis_id,)
    )
    con.commit()
    orautils.exchange_partition(cur, "IPRSCAN", "IPM_GENE3D_MATCH_TMP",
                                "MV_IPRSCAN", "GENE3D")
    orautils.drop_table(cur, "IPRSCAN", "IPM_GENE3D_MATCH_TMP")
    cur.close()
    con.close()


def _update_hamap(url: str, analysis_id: int):
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    _init_ipm_table(cur, "IPM_HAMAP_MATCH_TMP")
    cur.execute(
        """
        INSERT /*+ APPEND */ INTO IPRSCAN.IPM_HAMAP_MATCH_TMP
        SELECT
          ANALYSIS_ID, UPI, METHOD_AC, SUBSTR(RELNO_MAJOR, 1, 4),
          SUBSTR(RELNO_MAJOR, 6, 7), SEQ_START, SEQ_END, 0, 0, 0, NULL , 0,
          SEQSCORE, 0, 0, 0, 0, MODEL_AC, NULL, FRAGMENTS
        FROM IPRSCAN.IPM_HAMAP_MATCH
        WHERE ANALYSIS_ID = :1
        """, (analysis_id,)
    )
    con.commit()
    orautils.exchange_partition(cur, "IPRSCAN", "IPM_HAMAP_MATCH_TMP",
                                "MV_IPRSCAN", "HAMAP")
    orautils.drop_table(cur, "IPRSCAN", "IPM_HAMAP_MATCH_TMP")
    cur.close()
    con.close()


def _update_mobidblite(url: str, analysis_id: int):
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    _init_ipm_table(cur, "IPM_MOBIDBLITE_MATCH_TMP")
    cur.execute(
        """
        INSERT /*+ APPEND */ INTO IPRSCAN.IPM_MOBIDBLITE_MATCH_TMP
        SELECT
          ANALYSIS_ID, UPI, METHOD_AC, RELNO_MAJOR, RELNO_MINOR, SEQ_START,
          SEQ_END, 0, 0, 0, NULL , 0, 0, 0, 0, 0, 0, MODEL_AC, SEQ_FEATURE,
          FRAGMENTS
        FROM IPRSCAN.IPM_MOBIDBLITE_MATCH
        WHERE ANALYSIS_ID = :1
        """, (analysis_id,)
    )
    con.commit()
    orautils.exchange_partition(cur, "IPRSCAN", "IPM_MOBIDBLITE_MATCH_TMP",
                                "MV_IPRSCAN", "MOBIDBLITE")
    orautils.drop_table(cur, "IPRSCAN", "IPM_MOBIDBLITE_MATCH_TMP")
    cur.close()
    con.close()


def _update_panther(url: str, analysis_id: int):
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    _init_ipm_table(cur, "IPM_PANTHER_MATCH_TMP")
    cur.execute(
        """
        INSERT /*+ APPEND */ INTO IPRSCAN.IPM_PANTHER_MATCH_TMP
        SELECT
          ANALYSIS_ID, UPI, METHOD_AC, RELNO_MAJOR, RELNO_MINOR, SEQ_START,
          SEQ_END, HMM_START, HMM_END, HMM_LENGTH, HMM_BOUNDS, SEQSCORE,
          SEQSCORE, SEQEVALUE, SEQEVALUE, ENV_START, ENV_END, MODEL_AC,
          NULL, FRAGMENTS
        FROM IPRSCAN.IPM_PANTHER_MATCH
        WHERE ANALYSIS_ID = :1
        """, (analysis_id,)
    )
    con.commit()
    orautils.exchange_partition(cur, "IPRSCAN", "IPM_PANTHER_MATCH_TMP",
                                "MV_IPRSCAN", "PANTHER")
    orautils.drop_table(cur, "IPRSCAN", "IPM_PANTHER_MATCH_TMP")
    cur.close()
    con.close()


def _update_pfam(url: str, analysis_id: int):
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    _init_ipm_table(cur, "IPM_PFAM_MATCH_TMP")
    cur.execute(
        """
        INSERT /*+ APPEND */ INTO IPRSCAN.IPM_PFAM_MATCH_TMP
        SELECT
          ANALYSIS_ID, UPI, METHOD_AC, RELNO_MAJOR, RELNO_MINOR, SEQ_START,
          SEQ_END, HMM_START, HMM_END, HMM_LENGTH, HMM_BOUNDS, SCORE,
          SEQSCORE, EVALUE, SEQEVALUE, ENV_START, ENV_END, MODEL_AC,
          NULL, FRAGMENTS
        FROM IPRSCAN.IPM_PFAM_MATCH
        WHERE ANALYSIS_ID = :1
        """, (analysis_id,)
    )
    con.commit()
    orautils.exchange_partition(cur, "IPRSCAN", "IPM_PFAM_MATCH_TMP",
                                "MV_IPRSCAN", "PFAM")
    orautils.drop_table(cur, "IPRSCAN", "IPM_PFAM_MATCH_TMP")
    cur.close()
    con.close()


def _update_phobius(url: str, analysis_id: int):
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    _init_ipm_table(cur, "IPM_PHOBIUS_MATCH_TMP")
    cur.execute(
        """
        INSERT /*+ APPEND */ INTO IPRSCAN.IPM_PHOBIUS_MATCH_TMP
        SELECT
          ANALYSIS_ID, UPI, METHOD_AC, RELNO_MAJOR, RELNO_MINOR, SEQ_START,
          SEQ_END, 0, 0, 0, NULL, 0, 0, 0, 0, 0, 0, MODEL_AC, NULL, FRAGMENTS
        FROM IPRSCAN.IPM_PHOBIUS_MATCH
        WHERE ANALYSIS_ID = :1
        """, (analysis_id,)
    )
    con.commit()
    orautils.exchange_partition(cur, "IPRSCAN", "IPM_PHOBIUS_MATCH_TMP",
                                "MV_IPRSCAN", "PHOBIUS")
    orautils.drop_table(cur, "IPRSCAN", "IPM_PHOBIUS_MATCH_TMP")
    cur.close()
    con.close()


def _update_pirsf(url: str, analysis_id: int):
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    _init_ipm_table(cur, "IPM_PIRSF_MATCH_TMP")
    cur.execute(
        """
        INSERT /*+ APPEND */ INTO IPRSCAN.IPM_PIRSF_MATCH_TMP
        SELECT
          ANALYSIS_ID, UPI, METHOD_AC, RELNO_MAJOR, RELNO_MINOR, SEQ_START,
          SEQ_END, HMM_START, HMM_END, HMM_LENGTH, HMM_BOUNDS, SCORE,
          SEQSCORE, EVALUE, SEQEVALUE, ENV_START, ENV_END, MODEL_AC,
          NULL, FRAGMENTS
        FROM IPRSCAN.IPM_PIRSF_MATCH
        WHERE ANALYSIS_ID = :1
        """, (analysis_id,)
    )
    con.commit()
    orautils.exchange_partition(cur, "IPRSCAN", "IPM_PIRSF_MATCH_TMP",
                                "MV_IPRSCAN", "PIRSF")
    orautils.drop_table(cur, "IPRSCAN", "IPM_PIRSF_MATCH_TMP")
    cur.close()
    con.close()


def _update_prints(url: str, analysis_id: int):
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    _init_ipm_table(cur, "IPM_PRINTS_MATCH_TMP")
    cur.execute(
        """
        INSERT /*+ APPEND */ INTO IPRSCAN.IPM_PRINTS_MATCH_TMP
        SELECT
          ANALYSIS_ID, UPI, METHOD_AC, RELNO_MAJOR, RELNO_MINOR, SEQ_START,
          SEQ_END, 0, 0, MOTIF_NUMBER, NULL, 0, SEQSCORE, PVALUE, SEQEVALUE,
          0, 0, MODEL_AC, GRAPHSCAN, FRAGMENTS
        FROM IPRSCAN.IPM_PRINTS_MATCH
        WHERE ANALYSIS_ID = :1
        """, (analysis_id,)
    )
    con.commit()
    orautils.exchange_partition(cur, "IPRSCAN", "IPM_PRINTS_MATCH_TMP",
                                "MV_IPRSCAN", "PRINTS")
    orautils.drop_table(cur, "IPRSCAN", "IPM_PRINTS_MATCH_TMP")
    cur.close()
    con.close()


def _update_prodom(url: str, analysis_id: int):
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    _init_ipm_table(cur, "IPM_PRODOM_MATCH_TMP")
    cur.execute(
        """
        INSERT /*+ APPEND */ INTO IPRSCAN.IPM_PRODOM_MATCH_TMP
        SELECT
          ANALYSIS_ID, UPI, METHOD_AC, RELNO_MAJOR, RELNO_MINOR, SEQ_START,
          SEQ_END, 0, 0, 0, NULL, SEQSCORE, SEQSCORE, SEQEVALUE, SEQEVALUE,
          0, 0, MODEL_AC, NULL, FRAGMENTS
        FROM IPRSCAN.IPM_PRODOM_MATCH
        WHERE ANALYSIS_ID = :1
        """, (analysis_id,)
    )
    con.commit()
    orautils.exchange_partition(cur, "IPRSCAN", "IPM_PRODOM_MATCH_TMP",
                                "MV_IPRSCAN", "PRODOM")
    orautils.drop_table(cur, "IPRSCAN", "IPM_PRODOM_MATCH_TMP")
    cur.close()
    con.close()


def _update_prosite_patterns(url: str, analysis_id: int):
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    _init_ipm_table(cur, "IPM_PROSITE_PATTERNS_MATCH_TMP")
    cur.execute(
        """
        INSERT /*+ APPEND */ INTO IPRSCAN.IPM_PROSITE_PATTERNS_MATCH_TMP
        SELECT
          ANALYSIS_ID, UPI, METHOD_AC, SUBSTR(RELNO_MAJOR, 1, 4),
          SUBSTR(RELNO_MAJOR, 6, 7), SEQ_START,
          SEQ_END, 0, 0, 0, LOCATION_LEVEL, 0, 0, 0, 0, 0, 0, MODEL_AC,
          ALIGNMENT, FRAGMENTS
        FROM IPRSCAN.IPM_PROSITE_PATTERNS_MATCH
        WHERE ANALYSIS_ID = :1
        """, (analysis_id,)
    )
    con.commit()
    orautils.exchange_partition(cur, "IPRSCAN", "IPM_PROSITE_PATTERNS_MATCH_TMP",
                                "MV_IPRSCAN", "PROSITE_PATTERNS")
    orautils.drop_table(cur, "IPRSCAN", "IPM_PROSITE_PATTERNS_MATCH_TMP")
    cur.close()
    con.close()


def _update_prosite_profiles(url: str, analysis_id: int):
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    _init_ipm_table(cur, "IPM_PROSITE_PROFILES_MATCH_TMP")
    cur.execute(
        """
        INSERT /*+ APPEND */ INTO IPRSCAN.IPM_PROSITE_PROFILES_MATCH_TMP
        SELECT
          ANALYSIS_ID, UPI, METHOD_AC, SUBSTR(RELNO_MAJOR, 1, 4),
          SUBSTR(RELNO_MAJOR, 6, 7), SEQ_START,
          SEQ_END, 0, 0, 0, NULL, 0, SEQSCORE, 0, 0, 0, 0, MODEL_AC,
          ALIGNMENT, FRAGMENTS
        FROM IPRSCAN.IPM_PROSITE_PROFILES_MATCH
        WHERE ANALYSIS_ID = :1
        """, (analysis_id,)
    )
    con.commit()
    orautils.exchange_partition(cur, "IPRSCAN", "IPM_PROSITE_PROFILES_MATCH_TMP",
                                "MV_IPRSCAN", "PROSITE_PROFILES")
    orautils.drop_table(cur, "IPRSCAN", "IPM_PROSITE_PROFILES_MATCH_TMP")
    cur.close()
    con.close()


def _update_sfld_matches(url: str, analysis_id: int):
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    _init_ipm_table(cur, "IPM_SFLD_MATCH_TMP")
    cur.execute(
        """
        INSERT /*+ APPEND */ INTO IPRSCAN.IPM_SFLD_MATCH_TMP
        SELECT
          ANALYSIS_ID, UPI, METHOD_AC, RELNO_MAJOR, RELNO_MINOR, SEQ_START,
          SEQ_END, HMM_START, HMM_END, HMM_LENGTH, HMM_BOUNDS, SCORE,
          SEQSCORE, EVALUE, SEQEVALUE, ENV_START, ENV_END, MODEL_AC, NULL,
          FRAGMENTS
        FROM IPRSCAN.IPM_SFLD_MATCH
        WHERE ANALYSIS_ID = :1
        """, (analysis_id,)
    )
    con.commit()
    orautils.exchange_partition(cur, "IPRSCAN", "IPM_SFLD_MATCH_TMP",
                                "MV_IPRSCAN", "SFLD")
    orautils.drop_table(cur, "IPRSCAN", "IPM_SFLD_MATCH_TMP")
    cur.close()
    con.close()


def _update_sfld_sites(url: str, analysis_id: int):
    con = cx_Oracle.connect(url)
    cur = con.cursor()

    cur.execute(
        """
        CREATE TABLE IPRSCAN.IPM_SFLD_SITE_TMP
        NOLOGGING
        AS
        SELECT *
        FROM IPRSCAN.IPM_SFLD_SITE
        WHERE ANALYSIS_ID = :1
        """, (analysis_id,)
    )
    con.commit()
    orautils.exchange_partition(cur, "IPRSCAN", "IPM_SFLD_MATCH_TMP",
                                "SITE", "SFLD")
    orautils.drop_table(cur, "IPRSCAN", "IPM_SFLD_SITE_TMP")
    cur.close()
    con.close()


def _update_signalp_euk(url: str, analysis_id: int):
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    _init_ipm_table(cur, "IPM_SIGNALP_EUK_MATCH_TMP")
    cur.execute(
        """
        INSERT /*+ APPEND */ INTO IPRSCAN.IPM_SIGNALP_EUK_MATCH_TMP
        SELECT
          ANALYSIS_ID, UPI, METHOD_AC, RELNO_MAJOR, RELNO_MINOR, SEQ_START,
          SEQ_END, 0, 0, 0, NULL, SEQSCORE, SEQSCORE, 0, 0, 0, 0, MODEL_AC,
          NULL, FRAGMENTS
        FROM IPRSCAN.IPM_SIGNALP_MATCH
        WHERE ANALYSIS_ID = :1
        """, (analysis_id,)
    )
    con.commit()
    orautils.exchange_partition(cur, "IPRSCAN", "IPM_SIGNALP_EUK_MATCH_TMP",
                                "MV_IPRSCAN", "SIGNALP_EUK")
    orautils.drop_table(cur, "IPRSCAN", "IPM_SIGNALP_EUK_MATCH_TMP")
    cur.close()
    con.close()


def _update_signalp_gram_neg(url: str, analysis_id: int):
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    _init_ipm_table(cur, "IPM_SIGNALP_GRAM_NEG_MATCH_TMP")
    cur.execute(
        """
        INSERT /*+ APPEND */ INTO IPRSCAN.IPM_SIGNALP_GRAM_NEG_MATCH_TMP
        SELECT
          ANALYSIS_ID, UPI, METHOD_AC, RELNO_MAJOR, RELNO_MINOR, SEQ_START,
          SEQ_END, 0, 0, 0, NULL, SEQSCORE, SEQSCORE, 0, 0, 0, 0, MODEL_AC,
          NULL, FRAGMENTS
        FROM IPRSCAN.IPM_SIGNALP_MATCH
        WHERE ANALYSIS_ID = :1
        """, (analysis_id,)
    )
    con.commit()
    orautils.exchange_partition(cur, "IPRSCAN", "IPM_SIGNALP_GRAM_NEG_MATCH_TMP",
                                "MV_IPRSCAN", "SIGNALP_GRAM_NEGATIVE")
    orautils.drop_table(cur, "IPRSCAN", "IPM_SIGNALP_GRAM_NEG_MATCH_TMP")
    cur.close()
    con.close()


def _update_signalp_gram_pos(url: str, analysis_id: int):
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    _init_ipm_table(cur, "IPM_SIGNALP_GRAM_POS_MATCH_TMP")
    cur.execute(
        """
        INSERT /*+ APPEND */ INTO IPRSCAN.IPM_SIGNALP_GRAM_POS_MATCH_TMP
        SELECT
          ANALYSIS_ID, UPI, METHOD_AC, RELNO_MAJOR, RELNO_MINOR, SEQ_START,
          SEQ_END, 0, 0, 0, NULL, SEQSCORE, SEQSCORE, 0, 0, 0, 0, MODEL_AC,
          NULL, FRAGMENTS
        FROM IPRSCAN.IPM_SIGNALP_MATCH
        WHERE ANALYSIS_ID = :1
        """, (analysis_id,)
    )
    con.commit()
    orautils.exchange_partition(cur, "IPRSCAN", "IPM_SIGNALP_GRAM_POS_MATCH_TMP",
                                "MV_IPRSCAN", "SIGNALP_GRAM_POSITIVE")
    orautils.drop_table(cur, "IPRSCAN", "IPM_SIGNALP_GRAM_POS_MATCH_TMP")
    cur.close()
    con.close()


def _update_smart(url: str, analysis_id: int):
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    _init_ipm_table(cur, "IPM_SMART_MATCH_TMP")
    cur.execute(
        """
        INSERT /*+ APPEND */ INTO IPRSCAN.IPM_SMART_MATCH_TMP
        SELECT
          ANALYSIS_ID, UPI, METHOD_AC, RELNO_MAJOR, RELNO_MINOR, SEQ_START,
          SEQ_END, HMM_START, HMM_END, HMM_LENGTH, HMM_BOUNDS, SCORE,
          SEQSCORE, EVALUE, SEQEVALUE, 0, 0, MODEL_AC, NULL, FRAGMENTS
        FROM IPRSCAN.IPM_SMART_MATCH
        WHERE ANALYSIS_ID = :1
        """, (analysis_id,)
    )
    con.commit()
    orautils.exchange_partition(cur, "IPRSCAN", "IPM_SMART_MATCH_TMP",
                                "MV_IPRSCAN", "SMART")
    orautils.drop_table(cur, "IPRSCAN", "IPM_SMART_MATCH_TMP")
    cur.close()
    con.close()


def _update_superfamily(url: str, analysis_id: int):
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    _init_ipm_table(cur, "IPM_SUPERFAMILY_MATCH_TMP")
    cur.execute(
        """
        INSERT /*+ APPEND */ INTO IPRSCAN.IPM_SUPERFAMILY_MATCH_TMP
        SELECT
          ANALYSIS_ID, UPI, METHOD_AC, RELNO_MAJOR, RELNO_MINOR, SEQ_START,
          SEQ_END, 0, 0, HMM_LENGTH, NULL, 0, 0, SEQEVALUE, SEQEVALUE, 0, 0,
          MODEL_AC, NULL, FRAGMENTS
        FROM IPRSCAN.IPM_SUPERFAMILY_MATCH
        WHERE ANALYSIS_ID = :1
        """, (analysis_id,)
    )
    con.commit()
    orautils.exchange_partition(cur, "IPRSCAN", "IPM_SUPERFAMILY_MATCH_TMP",
                                "MV_IPRSCAN", "SUPERFAMILY")
    orautils.drop_table(cur, "IPRSCAN", "IPM_SUPERFAMILY_MATCH_TMP")
    cur.close()
    con.close()


def _update_tigrfam(url: str, analysis_id: int):
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    _init_ipm_table(cur, "IPM_TIGRFAM_MATCH_TMP")
    cur.execute(
        """
        INSERT /*+ APPEND */ INTO IPRSCAN.IPM_TIGRFAM_MATCH_TMP
        SELECT
          ANALYSIS_ID, UPI, METHOD_AC, RELNO_MAJOR, RELNO_MINOR, SEQ_START,
          SEQ_END, HMM_START, HMM_END, HMM_LENGTH, HMM_BOUNDS, SCORE,
          SEQSCORE, EVALUE, SEQEVALUE, ENV_START, ENV_END, MODEL_AC,
          NULL, FRAGMENTS
        FROM IPRSCAN.IPM_TIGRFAM_MATCH
        WHERE ANALYSIS_ID = :1
        """, (analysis_id,)
    )
    con.commit()
    orautils.exchange_partition(cur, "IPRSCAN", "IPM_TIGRFAM_MATCH_TMP",
                                "MV_IPRSCAN", "TIGRFAM")
    orautils.drop_table(cur, "IPRSCAN", "IPM_TIGRFAM_MATCH_TMP")
    cur.close()
    con.close()


def _update_tmhmm(url: str, analysis_id: int):
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    _init_ipm_table(cur, "IPM_TMHMM_MATCH_TMP")
    cur.execute(
        """
        INSERT /*+ APPEND */ INTO IPRSCAN.IPM_TMHMM_MATCH_TMP
        SELECT
          ANALYSIS_ID, UPI, METHOD_AC, RELNO_MAJOR, RELNO_MINOR, SEQ_START,
          SEQ_END, 0, 0, 0, NULL, SEQSCORE, SEQSCORE, 0, 0, 0, 0, MODEL_AC,
          NULL, FRAGMENTS
        FROM IPRSCAN.IPM_TMHMM_MATCH
        WHERE ANALYSIS_ID = :1
        """, (analysis_id,)
    )
    con.commit()
    orautils.exchange_partition(cur, "IPRSCAN", "IPM_TMHMM_MATCH_TMP",
                                "MV_IPRSCAN", "TMHMM")
    orautils.drop_table(cur, "IPRSCAN", "IPM_TMHMM_MATCH_TMP")
    cur.close()
    con.close()


def update_mv_iprscan(user: str, dsn: str, **kwargs):
    ispro_tables = kwargs.pop("tables", True)
    exchange_partitions = kwargs.pop("partitions", True)
    rebuild_indices = kwargs.pop("indices", True)

    url = orautils.make_connect_string(user, dsn)
    analyses = _check_ispro(url, **kwargs)

    with ThreadPoolExecutor(max_workers=len(analyses)) as executor:
        if ispro_tables:
            logger.info("copying tables from ISPRO")
            fs = {}
            for analysis in analyses:
                table = analysis["match_table"].upper()

                if table not in fs.values():
                    # Some analyses are in the same table!
                    f = executor.submit(_import_table, url, "IPRSCAN", table)
                    fs[f] = table

            num_errors = 0
            for f in as_completed(fs):
                table = fs[f]
                exc = f.exception()
                if exc is None:
                    logger.debug("{}: imported".format(table))
                else:
                    logger.error("{}: {}".format(table, exc))
                    num_errors += 1

            if num_errors:
                raise RuntimeError("{} tables "
                                   "were not imported".format(num_errors))

        if exchange_partitions:
            logger.info("updating MV_IPRSCAN")
            functions = {
                "ipm_cdd_match": _update_cdd_matches,
                "ipm_coils_match": _update_coils,
                "ipm_gene3d_match": _update_gene3d,
                "ipm_hamap_match": _update_hamap,
                "ipm_mobidblite_match": _update_mobidblite,
                "ipm_panther_match": _update_panther,
                "ipm_pfam_match": _update_pfam,
                "ipm_phobius_match": _update_phobius,
                "ipm_pirsf_match": _update_pirsf,
                "ipm_prints_match": _update_prints,
                "ipm_prodom_match": _update_prodom,
                "ipm_prosite_patterns_match": _update_prosite_patterns,
                "ipm_prosite_profiles_match": _update_prosite_profiles,
                "ipm_sfld_match": _update_sfld_matches,
                "ipm_smart_match": _update_smart,
                "ipm_superfamily_match": _update_superfamily,
                "ipm_tigrfam_match": _update_tigrfam,
                "ipm_tmhmm_match": _update_tmhmm
            }

            signalp = {
                "signalp_euk": _update_signalp_euk,
                "signalp_gram_negative": _update_signalp_gram_neg,
                "signalp_gram_positive": _update_signalp_gram_pos
            }

            fs = {}
            for analysis in analyses:
                _id = analysis["id"]
                name = analysis["name"]
                table = analysis["match_table"].lower()

                if table == "ipm_signalp_match":
                    fn = signalp[name]
                else:
                    fn = functions[table]

                f = executor.submit(fn, url, _id)
                fs[f] = analysis["full_name"]

            num_errors = 0
            for f in as_completed(fs):
                full_name = fs[f]
                exc = f.exception()
                if exc is None:
                    logger.debug("{}: exchanged".format(full_name))
                else:
                    logger.error("{}: {}".format(full_name, exc))
                    num_errors += 1

            if num_errors:
                raise RuntimeError("{} partitions "
                                   "were not exchanged".format(num_errors))

    if rebuild_indices:
        logger.info("indexing table")
        con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
        cur = con.cursor()

        for idx in orautils.get_indices(cur, "IPRSCAN", "MV_IPRSCAN"):
            # logger.debug("rebuilding {}".format(idx))
            # cur.execute("ALTER INDEX {} REBUILD NOLOGGING".format(idx))
            orautils.drop_index(cur, "IPRSCAN", idx["name"])

        cur.execute(
            """
            CREATE INDEX I_MV_IPRSCAN$UPI
            ON IPRSCAN.MV_IPRSCAN (UPI) NOLOGGING
            """
        )

        cur.close()
        con.close()

    logger.info("MV_IPRSCAN is ready")


def update_iprscan_site(user: str, dsn: str, **kwargs):
    ispro_tables = kwargs.pop("tables", True)
    exchange_partitions = kwargs.pop("partitions", True)
    rebuild_indices = kwargs.pop("indices", True)

    url = orautils.make_connect_string(user, dsn)
    analyses = _check_ispro(url, **kwargs)

    with ThreadPoolExecutor(max_workers=len(analyses)) as executor:
        if ispro_tables:
            logger.info("copying tables from ISPRO")
            fs = {}
            for analysis in analyses:
                table = analysis["site_table"].upper()

                if table not in fs.values():
                    # Some analyses are in the same table!
                    f = executor.submit(_import_table, url, "IPRSCAN", table)
                    fs[f] = table

            num_errors = 0
            for f in as_completed(fs):
                table = fs[f]
                exc = f.exception()
                if exc is None:
                    logger.debug("{}: imported".format(table))
                else:
                    logger.error("{}: {}".format(table, exc))
                    num_errors += 1

            if num_errors:
                raise RuntimeError("{} tables "
                                   "were not imported".format(num_errors))

        if exchange_partitions:
            logger.info("updating MV_IPRSCAN")
            functions = {
                "ipm_cdd_match": update_cdd,
                "ipm_coils_match": update_coils,
                "ipm_gene3d_match": update_gene3d,
                "ipm_hamap_match": update_hamap,
                "ipm_mobidblite_match": update_mobidblite,
                "ipm_panther_match": update_panther,
                "ipm_pfam_match": update_pfam,
                "ipm_phobius_match": update_phobius,
                "ipm_pirsf_match": update_pirsf,
                "ipm_prints_match": update_prints,
                "ipm_prodom_match": update_prodom,
                "ipm_prosite_patterns_match": update_prosite_patterns,
                "ipm_prosite_profiles_match": update_prosite_profiles,
                "ipm_sfld_match": update_sfld,
                "ipm_smart_match": update_smart,
                "ipm_superfamily_match": update_superfamily,
                "ipm_tigrfam_match": update_tigrfam,
                "ipm_tmhmm_match": update_tmhmm
            }

            signalp = {
                "signalp_euk": update_signalp_euk,
                "signalp_gram_negative": update_signalp_gram_neg,
                "signalp_gram_positive": update_signalp_gram_pos
            }

            fs = {}
            for analysis in analyses:
                _id = analysis["id"]
                name = analysis["name"]
                table = analysis["match_table"].lower()

                if table == "ipm_signalp_match":
                    fn = signalp[name]
                else:
                    fn = functions[table]

                f = executor.submit(fn, url, _id)
                fs[f] = analysis["full_name"]

            num_errors = 0
            for f in as_completed(fs):
                full_name = fs[f]
                exc = f.exception()
                if exc is None:
                    logger.debug("{}: exchanged".format(full_name))
                else:
                    logger.error("{}: {}".format(full_name, exc))
                    num_errors += 1

            if num_errors:
                raise RuntimeError("{} partitions "
                                   "were not exchanged".format(num_errors))

    if rebuild_indices:
        logger.info("indexing table")
        con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
        cur = con.cursor()

        for idx in orautils.get_indices(cur, "IPRSCAN", "MV_IPRSCAN"):
            # logger.debug("rebuilding {}".format(idx))
            # cur.execute("ALTER INDEX {} REBUILD NOLOGGING".format(idx))
            orautils.drop_index(cur, "IPRSCAN", idx["name"])

        cur.execute(
            """
            CREATE INDEX I_MV_IPRSCAN$UPI
            ON IPRSCAN.MV_IPRSCAN (UPI) NOLOGGING
            """
        )

        cur.close()
        con.close()

    logger.info("MV_IPRSCAN is ready")

