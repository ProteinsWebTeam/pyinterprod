# -*- coding: utf-8 -*-

from typing import Optional

import cx_Oracle

from ... import orautils


def load_databases(user: str, dsn: str):
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    orautils.drop_table(cur, owner, "CV_DATABASE")
    cur.execute(
        """
        CREATE TABLE {}.CV_DATABASE
        (
            DBCODE VARCHAR2(10) NOT NULL,
            DBNAME VARCHAR2(50) NOT NULL,
            DBSHORT VARCHAR2(10) NOT NULL,
            VERSION VARCHAR2(20),
            FILE_DATE DATE,
            IS_READY CHAR(1) DEFAULT 'N',
            CONSTRAINT PK_DATABASE PRIMARY KEY (DBCODE)
        ) NOLOGGING
        """.format(owner)
    )

    cur.execute(
        """
        INSERT /*+ APPEND */ INTO {}.CV_DATABASE (
            DBCODE, DBNAME, DBSHORT, VERSION, FILE_DATE
        )
        SELECT DB.DBCODE, DB.DBNAME, DB.DBSHORT, V.VERSION, V.FILE_DATE
        FROM INTERPRO.CV_DATABASE DB
        LEFT OUTER JOIN INTERPRO.DB_VERSION V
          ON DB.DBCODE = V.DBCODE
        """.format(owner)
    )
    con.commit()

    orautils.gather_stats(cur, owner, "CV_DATABASE")
    orautils.grant(cur, owner, "CV_DATABASE", "SELECT", "INTERPRO_SELECT")
    cur.close()
    con.close()


def load_matches(user: str, dsn: str):
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    orautils.drop_table(cur, owner, "MATCH")
    cur.execute(
        """
        CREATE TABLE {}.MATCH
        (
            PROTEIN_AC VARCHAR2(15) NOT NULL,
            METHOD_AC VARCHAR2(25) NOT NULL,
            DBCODE CHAR(1) NOT NULL,
            MODEL_AC VARCHAR2(25),
            POS_FROM NUMBER(5) NOT NULL,
            POS_TO NUMBER(5) NOT NULL,
            FRAGMENTS VARCHAR2(200) DEFAULT NULL
        ) NOLOGGING
        """.format(owner)
    )

    """
    Insert directly signature matches and MobiDB-lite matches
    For some signature matches, the HMM model accession
        is the signature accession itself: in such cases, use NULL
    """
    cur.execute(
        """
        INSERT /*+ APPEND */ INTO {}.MATCH
        SELECT
          PROTEIN_AC, METHOD_AC, DBCODE,
          CASE
            WHEN MODEL_AC IS NOT NULL AND MODEL_AC != METHOD_AC
            THEN MODEL_AC
            ELSE NULL
          END AS MODEL_AC,
          POS_FROM, POS_TO, FRAGMENTS
        FROM INTERPRO.MATCH
        UNION ALL
        SELECT
          PROTEIN_AC, METHOD_AC, DBCODE,
          NULL,
          POS_FROM, POS_TO, NULL
        FROM INTERPRO.FEATURE_MATCH
        WHERE DBCODE = 'g'
        """.format(owner)
    )
    con.commit()

    orautils.gather_stats(cur, owner, "MATCH")
    orautils.grant(cur, owner, "MATCH", "SELECT", "INTERPRO_SELECT")

    cur.execute(
        """
        CREATE INDEX I_MATCH$PROTEIN
        ON {}.MATCH (PROTEIN_AC) NOLOGGING
        """.format(owner)
    )
    con.execute(
        """
        CREATE INDEX I_MATCH$METHOD
        ON {}.MATCH (METHOD_AC) NOLOGGING
        """.format(owner)
    )
    con.execute(
        """
        CREATE INDEX I_MATCH$DBCODE
        ON {}.MATCH (DBCODE) NOLOGGING
        """.format(owner)
    )

    cur.execute(
        """
        SELECT METHOD_AC, COUNT(DISTINCT PROTEIN_AC)
        FROM {}.MATCH
        GROUP BY METHOD_AC
        """.format(owner)
    )
    signatures = cur.fetchall()

    for acc, num_proteins in signatures:
        cur.execute(
            """
            UPDATE {0}.METHOD
            SET PROTEIN_COUNT = :1
            WHERE METHOD_AC = :2
            """.format(owner), (num_proteins, acc)
        )
    con.commit()
    cur.close()
    con.close()


def load_signatures(user: str, dsn: str):
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    orautils.drop_table(cur, owner, "METHOD")
    cur.execute(
        """
        CREATE TABLE {}.METHOD
        (
            METHOD_AC VARCHAR2(25) NOT NULL,
            NAME VARCHAR2(100),
            DBCODE CHAR(1) NOT NULL,
            CANDIDATE CHAR(1) NOT NULL,
            DESCRIPTION VARCHAR2(220),
            SIG_TYPE CHAR(1),
            ABSTRACT VARCHAR2(4000),
            ABSTRACT_LONG CLOB,
            PROTEIN_COUNT NUMBER(8) DEFAULT 0 NOT NULL
        ) NOLOGGING
        """.format(owner)
    )

    cur.execute(
        """
        INSERT /*+ APPEND */ INTO {}.METHOD (
            METHOD_AC, NAME, DBCODE, CANDIDATE,
            DESCRIPTION, SIG_TYPE, ABSTRACT, ABSTRACT_LONG
        )
        SELECT
            METHOD_AC, NAME, DBCODE, CANDIDATE,
            DESCRIPTION, SIG_TYPE, ABSTRACT, ABSTRACT_LONG
        FROM INTERPRO.METHOD
        """.format(owner)
    )
    con.commit()

    cur.execute(
        """
        ALTER TABLE {}.METHOD
        ADD CONSTRAINT PK_METHOD PRIMARY KEY (METHOD_AC)
        """.format(owner)
    )
    cur.execute(
        """
        CREATE INDEX I_METHOD$DBCODE
        ON {}.METHOD (DBCODE) NOLOGGING
        """.format(owner)
    )

    orautils.gather_stats(cur, owner, "METHOD")
    orautils.grant(cur, owner, "METHOD", "SELECT", "INTERPRO_SELECT")
    cur.close()
    con.close()


def load_taxa(user: str, dsn: str):
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    orautils.drop_table(cur, owner, "ETAXI")
    cur.execute(
        """
        CREATE TABLE {}.ETAXI
        NOLOGGING
        AS
        SELECT 
          TAX_ID, PARENT_ID, SCIENTIFIC_NAME, RANK, 
          LEFT_NUMBER, RIGHT_NUMBER, FULL_NAME
        FROM INTERPRO.ETAXI
        """.format(owner)
    )
    orautils.gather_stats(cur, owner, "ETAXI")
    orautils.grant(cur, owner, "ETAXI", "SELECT", "INTERPRO_SELECT")
    cur.execute(
        """
        ALTER TABLE {}.ETAXI
        ADD CONSTRAINT PK_ETAXI PRIMARY KEY (TAX_ID)
        """.format(owner)
    )

    orautils.drop_table(cur, owner, "LINEAGE")
    cur.execute(
        """
        CREATE TABLE {}.LINEAGE
        (
            LEFT_NUMBER NUMBER NOT NULL,
            TAX_ID NUMBER(10) NOT NULL,
            RANK VARCHAR2(50)
        ) NOLOGGING
        """.format(owner)
    )

    taxons = {}
    cur.execute(
        """
        SELECT TAX_ID, PARENT_ID, LEFT_NUMBER, RANK
        FROM {}.ETAXI
        """.format(owner)
    )
    for tax_id, parent_id, left_num, rank in cur:
        if tax_id != 131567:
            """
            taxID 131567 (cellular organisms) contains three superkingdoms:
                * Bacteria (2)
                * Archaea (2157)
                * Eukaryota (2759)

            therefore it is not needed (we don't want a meta-superkingdom)
            """
            taxons[tax_id] = {
                "parent": parent_id,
                "left_num": left_num,
                "rank": "superkingdom" if parent_id == 1 else rank
            }

    lineage = []
    for tax_id in taxons:
        t = taxons[tax_id]
        left_num = t["left_num"]
        if not left_num:
            continue
        elif t["rank"] != "no rank":
            lineage.append((left_num, tax_id, t["rank"]))

        parent_id = t["parent"]
        while parent_id:
            if parent_id not in taxons:
                # taxID 131567 missing from dictionary
                break

            t = taxons[parent_id]
            if t["rank"] != "no rank":
                lineage.append((left_num, parent_id, t["rank"]))

            parent_id = t["parent"]

    table = orautils.TablePopulator(con,
                                    query=" INSERT /*+ APPEND */ "
                                          "INTO {}.LINEAGE "
                                          "VALUES (:1, :2, :3)".format(owner),
                                    autocommit=True)
    for record in lineage:
        table.insert(record)
    table.close()

    orautils.gather_stats(cur, owner, "LINEAGE")
    orautils.grant(cur, owner, "LINEAGE", "SELECT", "INTERPRO_SELECT")
    cur.execute(
        """
        CREATE INDEX I_LINEAGE$L$R
        ON {}.LINEAGE (LEFT_NUMBER, RANK)
        NOLOGGING
        """.format(owner)
    )

    cur.close()
    con.close()


def load_proteins(user: str, dsn: str):
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    orautils.drop_table(cur, owner, "PROTEIN")
    cur.execute(
        """
        CREATE TABLE {}.PROTEIN
        NOLOGGING
        AS
        SELECT PROTEIN_AC, NAME, DBCODE, LEN, FRAGMENT, TAX_ID
        FROM INTERPRO.PROTEIN
        """.format(owner)
    )
    orautils.gather_stats(cur, owner, "PROTEIN")
    orautils.grant(cur, owner, "PROTEIN", "SELECT", "INTERPRO_SELECT")
    cur.execute(
        """
        ALTER TABLE {}.PROTEIN
        ADD CONSTRAINT PK_PROTEIN PRIMARY KEY (PROTEIN_AC)
        """.format(owner)
    )
    cur.execute(
        """
        CREATE INDEX I_PROTEIN$DBCODE
        ON {}.PROTEIN (DBCODE) NOLOGGING
        """.format(owner)
    )
    cur.execute(
        """
        CREATE INDEX I_PROTEIN$NAME
        ON {}.PROTEIN (NAME) NOLOGGING
        """.format(owner)
    )

    cur.close()
    con.close()


def load_signature2protein(user: str, dsn: str, processes: int=1,
                           tmpdir: Optional[str]=None):
    pass


def copy_schema(user_src: str, user_dst: str, dsn: str):
    _enable_schema(user_src, dsn)
    orautils.clear_schema(user_dst, dsn)


def _enable_schema(user: str, dsn: str):
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    cur.execute(
        """
        UPDATE {}.CV_DATABASE
        SET IS_READY = 'Y'
        """.format(owner)
    )
    con.commit()
    cur.close()
    con.close()
