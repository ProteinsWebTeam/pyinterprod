import cx_Oracle

from .. import logger, orautils
from ..proteinupdate import sendmail
 
def notify_curators(mail_interpro):
    sendmail.send_mail(
        to_addrs=mail_interpro,
        subject="[CURATORS] MV tables update in progress",
        content="""\
Dear Curators,

MV tables are being updated at the moment.
Pronto may not function properly.
This will take approximately 20 hours and you will be notified when this is done

Best regards
The production team"""
    )


def refresh_mv_and_match_stats(user: str, dsn: str, plsql: bool=True):
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    
    logger.info("refresh mv tables")
    if plsql:
        cur.callproc('INTERPRO.REFRESH_MATCH_COUNTS.REFRESH') 
    #else:
    logger.info("DONE: refresh mv tables")
    logger.info("refresh match stats")

    # cur.execute('CREATE TABLE INTERPRO.MATCH_STATS_OLD AS SELECT * FROM INTERPRO.MATCH_STATS') 
    # orautils.drop_table(cur, "INTERPRO", "MATCH_STATS", purge=True)
    # cur.execute("ALTER SESSION FORCE parallel dml")
    # cur.execute("""INSERT /*+ APPEND PARALLEL */
    #             INTO INTERPRO.MATCH_STATS
    #           SELECT C.DBNAME
    #                 ,C.DBCODE
    #                 ,M1.STATUS
    #                 ,COUNT(M1.STATUS) AS COUNT
    #            FROM INTERPRO.CV_DATABASE C
    #             JOIN INTERPRO.MATCH M1 ON C.DBCODE = M1.DBCODE
    #           GROUP BY C.DBNAME, C.DBCODE, M1.STATUS
    #         """)
    orautils.drop_table(cur, "INTERPRO", "MATCH_STATS_OLD", purge=True)
    con.commit()
    logger.info("DONE: refresh match stats")

    cur.close()
    con.close()
    
def refresh_tax_tables(user: str, dsn: str):
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()

    logger.info("creating TAXONOMY_LOAD")
    orautils.drop_table(cur, "INTERPRO", "TAXONOMY_LOAD", purge=True)

    cur.execute("""CREATE TABLE INTERPRO.TAXONOMY_LOAD NOLOGGING AS 
                SELECT /*+ PARALLEL */ 
                  N.TAX_ID, 
                  N.PARENT_ID, 
                  N.SPTR_SCIENTIFIC SCIENTIFIC_NAME, 
                  N.RANK, 
                  NVL(N.SPTR_COMMON, N.NCBI_COMMON) COMMON_NAME 
                FROM TAXONOMY.V_PUBLIC_NODE@SWPREAD N""")
    con.commit()

    orautils.gather_stats(cur, "INTERPRO", "TAXONOMY_LOAD")

    # cur.execute('GRANT DELETE, INSERT, UPDATE ON INTERPRO.TAXONOMY_LOAD TO INTERPRO_PRODUCTION')
    # cur.execute('GRANT SELECT ON TAXONOMY_LOAD TO INTERPRO_SELECT')
    # cur.execute('GRANT INSERT ON TAXONOMY_LOAD TO INTERPRO_WEBSERVER')
    # con.commit()

    # cur.callproc('INTERPRO.IPRO_UTL_PKG.TABLE_STATS', ['TAXONOMY_LOAD',])
    # con.commit()

    logger.info('computing left and right numbers for taxonomy tree')

    logger.info("creating ETAXI")
    cur.execute('DROP TABLE INTERPRO.ETAXI')

    # Table of taxonomic classifications
    cur.execute("""CREATE TABLE INTERPRO.ETAXI NOLOGGING
                AS SELECT /*+ PARALLEL */ 
                  N.TAX_ID, N.PARENT_ID, N.SCIENTIFIC_NAME, 'X' COMPLETE_GENOME_FLAG, 
                  N.RANK, 0 HIDDEN, LR.TREE_LEFT LEFT_NUMBER, LR.TREE_RIGHT RIGHT_NUMBER, 'X' ANNOTATION_SOURCE, 
                  N.SCIENTIFIC_NAME || CASE WHEN N.COMMON_NAME IS NULL THEN '' ELSE ' (' || N.COMMON_NAME || ')' END FULL_NAME 
                FROM INTERPRO.TAXONOMY_LOAD N 
                JOIN (
                  SELECT TAX_ID, MIN(TREE_NUMBER) TREE_LEFT, MAX(TREE_NUMBER) TREE_RIGHT
                  FROM (
                    SELECT PARENT_ID AS TAX_ID, ROWNUM AS TREE_NUMBER 
                    FROM (
                       SELECT TAX_ID, PARENT_ID 
                       FROM (
                         SELECT TAX_ID, PARENT_ID FROM INTERPRO.TAXONOMY_LOAD 
                         UNION ALL 
                         SELECT 9999999 AS TAX_ID, TAX_ID AS PARENT_ID FROM INTERPRO.TAXONOMY_LOAD 
                         UNION ALL 
                         SELECT 0 AS TAX_ID, TAX_ID AS PARENT_ID FROM INTERPRO.TAXONOMY_LOAD
                       ) 
                       START WITH TAX_ID = 1 
                       CONNECT BY PRIOR TAX_ID=PARENT_ID 
                       ORDER SIBLINGS BY TAX_ID
                    ) 
                    WHERE TAX_ID IN (9999999, 0)
                  )
                  GROUP BY TAX_ID 
                ) LR 
                ON (LR.TAX_ID = N.TAX_ID)""")
    con.commit()

    cur.execute('CREATE INDEX ETAXI$L$R$T ON INTERPRO.ETAXI (LEFT_NUMBER, RIGHT_NUMBER, TAX_ID) NOLOGGING')
    cur.execute('CREATE INDEX ETAXI$P$T$R ON INTERPRO.ETAXI (PARENT_ID, TAX_ID, RANK) TABLESPACE NOLOGGING')
    cur.execute('CREATE INDEX ETAXI$T$P$R ON INTERPRO.ETAXI (TAX_ID, PARENT_ID, RANK) TABLESPACE NOLOGGING')

    orautils.grant(cur, "INTERPRO", "ETAXI", "SELECT", "PUBLIC")
    

    # cur.execute('GRANT SELECT ON INTERPRO.ETAXI TO INTERPRO_DEVELOPER')
    # cur.execute('GRANT ALTER, DELETE, INSERT, SELECT, UPDATE, ON COMMIT REFRESH, QUERY REWRITE, DEBUG, FLASHBACK ON INTERPRO.ETAXI TO INTERPRO_PRODUCTION')
    # cur.execute('GRANT SELECT ON INTERPRO.ETAXI TO INTERPRO_SELECT')
    # cur.execute('GRANT SELECT ON INTERPRO.ETAXI TO INTERPRO_WEBSERVER')
    # cur.execute('GRANT SELECT ON INTERPRO.ETAXI TO PUBLIC')
    # con.commit()

    # Refresh stats
    #cur.callproc('INTERPRO.IPRO_UTL_PKG.TABLE_STATS', ('ETAXI',))
    orautils.gather_stats(cur, "INTERPRO", "ETAXI")
    con.commit()

    # Dropping temporary table
    orautils.drop_table(cur, "INTERPRO", "TAXONOMY_LOAD", purge=True)

    logger.info("creating UNIPROT_TAXONOMY by combining protein taxonomy with left numbers")
    orautils.drop_table(cur, "INTERPRO", "UNIPROT_TAXONOMY", purge=True)

    cur.execute("""CREATE TABLE INTERPRO.UNIPROT_TAXONOMY NOLOGGING AS
                SELECT /*+ PARALLEL */
                  P.PROTEIN_AC,
                  P.TAX_ID,
                  NVL(ET.LEFT_NUMBER, 0) LEFT_NUMBER,
                  NVL(ET.RIGHT_NUMBER, 0) RIGHT_NUMBER
                FROM INTERPRO.PROTEIN P
                LEFT OUTER JOIN INTERPRO.ETAXI ET ON (P.TAX_ID=ET.TAX_ID)""")
    con.commit()

    cur.execute("CREATE INDEX UNIPROT_TAXONOMY$L$P ON INTERPRO.UNIPROT_TAXONOMY (LEFT_NUMBER, PROTEIN_AC) NOLOGGING")
    cur.execute("CREATE INDEX UNIPROT_TAXONOMY$P$L ON INTERPRO.UNIPROT_TAXONOMY (PROTEIN_AC, LEFT_NUMBER) NOLOGGING")
    con.commit()

    # Refresh stats
    # cur.callproc('INTERPRO.IPRO_UTL_PKG.TABLE_STATS', ('UNIPROT_TAXONOMY',))
    orautils.gather_stats(cur, "INTERPRO", "UNIPROT_TAXONOMY")
    # con.commit()
    # orautils.grant(cur, "INTERPRO", "UNIPROT_TAXONOMY", "SELECT", "PUBLIC")
    # cur.execute('GRANT SELECT ON UNIPROT_TAXONOMY TO INTERPRO_DEVELOPER')
    # cur.execute('GRANT ALTER, DELETE, INSERT, SELECT, UPDATE, ON COMMIT REFRESH, QUERY REWRITE, DEBUG, FLASHBACK ON UNIPROT_TAXONOMY TO INTERPRO_PRODUCTION')
    # cur.execute('GRANT SELECT ON UNIPROT_TAXONOMY TO INTERPRO_SELECT')
    # cur.execute('GRANT SELECT ON UNIPROT_TAXONOMY TO INTERPRO_WEBSERVER')
    # cur.execute('GRANT SELECT ON UNIPROT_TAXONOMY TO PUBLIC')
    con.commit()

    logger.info("creating MV_TAX_ENTRY_COUNT")
    orautils.drop_table(cur, "INTERPRO", "MV_TAX_ENTRY_COUNT", purge=True)

    """
    Count of proteins with true matches to InterPro entries
        ENTRY_AC                        InterPro entry
        TAX_ID                          Taxonomic ID of proteins matching entry
        COUNT                           Count of proteins for this entry and tax Id, also including any child tax Ids
        COUNT_SPECIFIED_TAX_ID          Count of proteins for this entry and this tax Id only
    """
    cur.execute("CREATE TABLE INTERPRO.MV_TAX_ENTRY_COUNT NOLOGGING AS "
                "WITH QUERY1 AS ("
                "  SELECT ENTRY_AC, ANC.PARENT AS TAX_ID, COUNT(1) AS COUNT"
                "  FROM INTERPRO.UNIPROT_TAXONOMY UT "
                "  JOIN INTERPRO.MV_ENTRY2PROTEIN_TRUE MVEP "
                "ON UT.PROTEIN_AC=MVEP.PROTEIN_AC "
                "  JOIN ("
                "    SELECT "
                "      NVL(SUBSTR(SYS_CONNECT_BY_PATH(TAX_ID, '.'), 2, INSTR(SYS_CONNECT_BY_PATH (TAX_ID,'.'),'.',2) - 2), TAX_ID) AS CHILD, "
                "      TAX_ID AS PARENT "
                "    FROM INTERPRO.ETAXI ET "
                "    CONNECT BY PRIOR PARENT_ID=TAX_ID"
                "  ) ANC "
                "  ON ANC.CHILD=UT.TAX_ID "
                "  GROUP BY ENTRY_AC, ANC.PARENT"
                "), QUERY2 AS ("
                "  SELECT ENTRY_AC, TAX_ID, COUNT(1) AS COUNT "
                "  FROM INTERPRO.UNIPROT_TAXONOMY UT "
                "  JOIN INTERPRO.MV_ENTRY2PROTEIN_TRUE MVEP "
                "  ON UT.PROTEIN_AC=MVEP.PROTEIN_AC "
                "  GROUP BY ENTRY_AC, TAX_ID"
                ")"
                "SELECT /*+ PARALLEL */ QUERY1.ENTRY_AC, QUERY1.TAX_ID, QUERY1.COUNT AS COUNT, QUERY2.COUNT AS COUNT_SPECIFIED_TAX_ID "
                "FROM QUERY1 "
                "LEFT OUTER JOIN QUERY2 "
                "ON QUERY1.ENTRY_AC = QUERY2.ENTRY_AC AND QUERY1.TAX_ID = QUERY2.TAX_ID")
    con.commit()

    cur.execute('ALTER TABLE INTERPRO.MV_TAX_ENTRY_COUNT '
                'ADD CONSTRAINT PK_MV_TAX_ENTRY_COUNT '
                'PRIMARY KEY (ENTRY_AC, TAX_ID) '
                'USING INDEX TABLESPACE INTERPRO_IND')

    cur.execute('CREATE UNIQUE INDEX TEC_PERF_IND1 '
                'ON INTERPRO.MV_TAX_ENTRY_COUNT (TAX_ID, ENTRY_AC) '
                'TABLESPACE INTERPRO_IND')

    orautils.grant(cur, "INTERPRO", "MV_TAX_ENTRY_COUNT", "SELECT",
                   "INTERPRO_SELECT")            
    # cur.execute('GRANT ALTER, DELETE, INSERT, SELECT, UPDATE, ON COMMIT REFRESH, QUERY REWRITE, DEBUG, FLASHBACK '
    #            'ON INTERPRO.MV_TAX_ENTRY_COUNT TO INTERPRO_PRODUCTION')
    # cur.execute('GRANT SELECT ON INTERPRO.MV_TAX_ENTRY_COUNT TO INTERPRO_SELECT')
    con.commit()

    # Refresh stats
    orautils.gather_stats(cur, "INTERPRO", "MV_TAX_ENTRY_COUNT")
    #cur.callproc('INTERPRO.IPRO_UTL_PKG.TABLE_STATS', ('MV_TAX_ENTRY_COUNT',))
    con.commit()

    cur.close()
    con.close()

    
def refresh_tables(user: str, dsn: str, memberdb: str, mail_interpro: list, wdir: str):

    notify_curators(mail_interpro)

    refresh_mv_and_match_stats(user, dsn)

    refresh_tax_tables(*user, dsn)
