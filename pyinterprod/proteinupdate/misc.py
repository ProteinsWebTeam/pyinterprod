# -*- coding: utf-8 -*-

import os
import time
from datetime import datetime, timedelta
from typing import List
from tempfile import mkdtemp
from zipfile import ZipFile, ZIP_DEFLATED

import cx_Oracle

from .. import logger, orautils
from .sendmail import send_mail


def refresh_mviews(user: str, dsn: str, notice: int=3600, plsql: bool=True):
    date = datetime.now() + timedelta(seconds=notice)
    send_mail(
        to_addrs=["interpro-team@ebi.ac.uk"],
        subject="MV update scheduled",
        content="""\
Dear curators,

Legacy materialized views will start being updated \
on {0:%a %d %b %Y} at {0:%H:%M}.

You may continue to use Pronto but please do not perform \
any of the following actions:
    - integrate a member database signature;
    - unintegrate a member database signature;
    - delete an Interpro entry (unless it does not have any signature).

InterPro entries may be created (but do NOT integrate any entry), \
modified (name, short name, abstracts, GO terms, etc.), or checked/unchecked.

An email notification will automatically be sent at the end of the update.

The InterPro Production Team
""".format(date)
    )

    time.sleep(notice)
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()

    # Ensure that legacy tables are dropped
    orautils.drop_table(cur, "INTERPRO", "MATCH_STATS")
    orautils.drop_table(cur, "INTERPRO", "MATCH_STATS_OLD")

    if plsql:
        cur.callproc("INTERPRO.REFRESH_MATCH_COUNTS.REFRESH")
    else:
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
        orautils.drop_table(cur, "INTERPRO", "MV_ENTRY2PROTEIN")
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
            ON INTERPRO.MV_ENTRY2PROTEIN_TRUE (ENTRY_AC)
            """
        )

    cur.close()
    con.close()

    send_mail(
        to_addrs=["interpro-team@ebi.ac.uk"],
        subject="MV update complete",
        content="""\
Dear curators,

Legacy materialized views are now updated. \
You may resume integrating or unintegrating signatures. Have fun!

The InterPro Production Team
""".format(date)
    )


def refresh_taxonomy(user: str, dsn: str):
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()

    logger.info("creating TAXONOMY_LOAD")
    orautils.drop_table(cur, "INTERPRO", "TAXONOMY_LOAD")
    cur.execute(
        """
        CREATE TABLE INTERPRO.TAXONOMY_LOAD
        NOLOGGING
        AS
        SELECT
          TAX_ID,
          PARENT_ID,
          SPTR_SCIENTIFIC AS SCIENTIFIC_NAME,
          RANK,
          NVL(SPTR_COMMON, NCBI_COMMON) AS COMMON_NAME
        FROM TAXONOMY.V_PUBLIC_NODE@SWPREAD
        """
    )
    orautils.gather_stats(cur, "INTERPRO", "TAXONOMY_LOAD")

    logger.info("creating ETAXI")
    orautils.drop_table(cur, "INTERPRO", "ETAXI")
    cur.execute(
        """
        CREATE TABLE INTERPRO.ETAXI
        AS
          SELECT
            N.TAX_ID, N.PARENT_ID, N.SCIENTIFIC_NAME,
            'X' COMPLETE_GENOME_FLAG, N.RANK, 0 HIDDEN,
            LR.TREE_LEFT LEFT_NUMBER, LR.TREE_RIGHT RIGHT_NUMBER,
            'X' ANNOTATION_SOURCE,
            N.SCIENTIFIC_NAME || CASE WHEN N.COMMON_NAME IS NULL THEN '' ELSE ' (' || N.COMMON_NAME || ')' END FULL_NAME
        FROM INTERPRO.TAXONOMY_LOAD N
        INNER JOIN (
            SELECT
              TAX_ID,
              MIN(TREE_NUMBER) TREE_LEFT,
              MAX(TREE_NUMBER) TREE_RIGHT
            FROM (
                SELECT PARENT_ID AS TAX_ID, ROWNUM AS TREE_NUMBER
                FROM (
                    SELECT TAX_ID, PARENT_ID
                    FROM (
                        SELECT TAX_ID, PARENT_ID
                        FROM INTERPRO.TAXONOMY_LOAD
                        UNION ALL
                        SELECT 9999999 AS TAX_ID, TAX_ID AS PARENT_ID
                        FROM INTERPRO.TAXONOMY_LOAD
                        UNION ALL
                        SELECT 0 AS TAX_ID, TAX_ID AS PARENT_ID
                        FROM INTERPRO.TAXONOMY_LOAD
                    )
                    START WITH TAX_ID = 1
                    CONNECT BY PRIOR TAX_ID=PARENT_ID
                    ORDER SIBLINGS BY TAX_ID
                )
                WHERE TAX_ID IN (9999999, 0)
            )
            GROUP BY TAX_ID
        ) LR
        ON LR.TAX_ID = N.TAX_ID
        """
    )
    orautils.gather_stats(cur, "INTERPRO", "ETAXI")
    cur.execute(
        """
        CREATE INDEX I_ETAXI$L$R$T
        ON INTERPRO.ETAXI (LEFT_NUMBER, RIGHT_NUMBER, TAX_ID)
        """
    )
    cur.execute(
        """
        CREATE INDEX I_ETAXI$P$T$R
        ON INTERPRO.ETAXI (PARENT_ID, TAX_ID, RANK)
        """
    )
    cur.execute(
        """
        CREATE INDEX I_ETAXI$T$P$R
        ON INTERPRO.ETAXI (TAX_ID, PARENT_ID, RANK)
        """
    )

    # Dropping temporary table
    orautils.drop_table(cur, "INTERPRO", "TAXONOMY_LOAD")

    orautils.grant(cur, "INTERPRO", "ETAXI", "SELECT", "PUBLIC")
    cur.close()
    con.close()
    logger.info("complete")


def report_to_curators(user: str, dsn: str, files: List[str]):
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    cur.execute("SELECT VERSION FROM INTERPRO.DB_VERSION WHERE DBCODE = 'u'")
    release = cur.fetchone()[0]
    cur.close()
    con.close()

    dirname = mkdtemp()
    filename = os.path.join(dirname, f"protein_update_{release}.zip")

    with ZipFile(filename, 'w', compression=ZIP_DEFLATED) as fh:
        for path in files:
            fh.write(path, arcname=os.path.basename(path))

    send_mail(
        to_addrs=["interpro-team@ebi.ac.uk"],
        subject=f"Protein update report: UniProt {release}",
        content="""\
Dear curators,

Pronto has been rereshed. Please find attached a ZIP archive containing \
the report files for this month's protein update.

The InterPro Production Team
""",
        attachments=[filename]
    )

    os.remove(filename)
    os.rmdir(dirname)
