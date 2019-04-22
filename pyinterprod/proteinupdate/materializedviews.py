import cx_Oracle

from .. import logger, orautils


def refresh_interpro_views(user: str, dsn: str):
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()

    logger.info("creating MV_METHOD2PROTEIN")
    orautils.drop_table(cur, "INTERPRO", "MV_METHOD2PROTEIN")
    cur.execute(
        """
        CREATE TABLE INTERPRO.MV_METHOD2PROTEIN (
            METHOD_AC VARCHAR2(25) NOT NULL ,
            PROTEIN_AC VARCHAR2(15) NOT NULL,
            MATCH_COUNT NUMBER(7) NOT NULL
        ) NOLOGGING
        """
    )
    cur.execute(
        """
        INSERT /*+ APPEND */ INTO INTERPRO.MV_METHOD2PROTEIN
        SELECT METHOD_AC, PROTEIN_AC, COUNT(*)
        FROM INTERPRO.MATCH
        GROUP BY METHOD_AC, PROTEIN_AC
        """
    )
    con.commit()

    logger.info("indexing MV_METHOD2PROTEIN")
    cur.execute(
        """
        CREATE INDEX I_MV_METHOD2PROTEIN
        ON INTERPRO.MV_METHOD2PROTEIN (METHOD_AC)
        """
    )

    logger.info("creating MV_METHOD_MATCH")
    orautils.drop_table(cur, "INTERPRO", "MV_METHOD_MATCH")
    cur.execute(
        """
        CREATE TABLE INTERPRO.MV_METHOD_MATCH (
            METHOD_AC VARCHAR2(25) NOT NULL,
            PROTEIN_COUNT NUMBER(8) NOT NULL,
            MATCH_COUNT NUMBER(8) NOT NULL
        ) NOLOGGING
        """
    )
    cur.execute(
        """
        INSERT /*+ APPEND */ INTO INTERPRO.MV_METHOD_MATCH
        SELECT METHOD_AC, COUNT(*), SUM(MATCH_COUNT)
        FROM INTERPRO.MV_METHOD2PROTEIN
        GROUP BY METHOD_AC
        """
    )
    con.commit()

    logger.info("indexing MV_METHOD_MATCH")
    cur.execute(
        """
        CREATE UNIQUE INDEX UI_MV_METHOD_MATCH
        ON INTERPRO.MV_METHOD_MATCH (METHOD_AC)
        """
    )

    logger.info("creating MV_ENTRY2PROTEIN_TRUE")
    orautils.drop_table(cur, "INTERPRO", "MV_ENTRY2PROTEIN_TRUE")
    cur.execute(
        """
        CREATE TABLE INTERPRO.MV_ENTRY2PROTEIN_TRUE (
            ENTRY_AC VARCHAR2(9) NOT NULL,
            PROTEIN_AC VARCHAR2(15) NOT NULL,
            MATCH_COUNT NUMBER(7) NOT NULL
        ) NOLOGGING
        """
    )
    cur.execute(
        """
        INSERT /*+ APPEND */ INTO INTERPRO.MV_ENTRY2PROTEIN_TRUE
        SELECT E.ENTRY_AC, M.PROTEIN_AC, COUNT(*)
        FROM INTERPRO.ENTRY2METHOD E
        INNER JOIN INTERPRO.MV_METHOD2PROTEIN M
            ON E.METHOD_AC = M.METHOD_AC
        GROUP BY E.ENTRY_AC, M.PROTEIN_AC
        """
    )
    con.commit()

    logger.info("indexing MV_ENTRY2PROTEIN_TRUE")
    cur.execute(
        """
        CREATE INDEX I_MV_ENTRY2PROTEIN_TRUE
        ON INTERPRO.MV_ENTRY2PROTEIN (ENTRY_AC)
        """
    )

    cur.close()
    con.close()


def refresh_iprscan_views():
    pass


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


def update_cdd(url: str, analysis_id: int):
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


def update_coils(url: str, analysis_id: int):
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


def update_gene3d(url: str, analysis_id: int):
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


def update_hamap(url: str, analysis_id: int):
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


def update_mobidblite(url: str, analysis_id: int):
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


def update_panther(url: str, analysis_id: int):
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


def update_pfam(url: str, analysis_id: int):
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


def update_phobius(url: str, analysis_id: int):
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


def update_pirsf(url: str, analysis_id: int):
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


def update_prints(url: str, analysis_id: int):
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


def update_prodom(url: str, analysis_id: int):
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


def update_prosite_patterns(url: str, analysis_id: int):
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


def update_prosite_profiles(url: str, analysis_id: int):
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


def update_sfld(url: str, analysis_id: int):
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


def update_signalp_euk(url: str, analysis_id: int):
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


def update_signalp_gram_neg(url: str, analysis_id: int):
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


def update_signalp_gram_pos(url: str, analysis_id: int):
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


def update_smart(url: str, analysis_id: int):
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


def update_superfamily(url: str, analysis_id: int):
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


def update_tigrfam(url: str, analysis_id: int):
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


def update_tmhmm(url: str, analysis_id: int):
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
