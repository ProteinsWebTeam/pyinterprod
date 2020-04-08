#!/usr/bin/env python3

import time
import re
import cx_Oracle
import sys
import os
from .. import logger, orautils
from .. import proteinupdate

def prepare_matches(cursor, connection, dbcode_list):
    logger.info("adding matches_tmp")
    orautils.drop_table(cursor, "INTERPRO", "MATCH_TMP", purge=True)
    logger.info("success dropping")

    query_part="CREATE TABLE INTERPRO.MATCH_TMP PARTITION BY LIST (DBCODE) \
        ( \
            PARTITION MATCH_DBCODE_R VALUES ('R'), \
            PARTITION MATCH_DBCODE_U VALUES ('U'), \
            PARTITION MATCH_DBCODE_P VALUES ('P'), \
            PARTITION MATCH_DBCODE_Y VALUES ('Y'), \
            PARTITION MATCH_DBCODE_D VALUES ('D'), \
            PARTITION MATCH_DBCODE_V VALUES ('V'), \
            PARTITION MATCH_DBCODE_M VALUES ('M'), \
            PARTITION MATCH_DBCODE_Q VALUES ('Q'), \
            PARTITION MATCH_DBCODE_X VALUES ('X'), \
            PARTITION MATCH_DBCODE_N VALUES ('N'), \
            PARTITION MATCH_DBCODE_H VALUES ('H'), \
            PARTITION MATCH_DBCODE_F VALUES ('F'), \
            PARTITION MATCH_DBCODE_J VALUES ('J'), \
            PARTITION MATCH_DBCODE_B VALUES ('B') \
        )AS SELECT * FROM INTERPRO.MATCH WHERE 1=0"
    #PARTITION MATCH_DBCODE_z VALUES ('z') \
    query_alter="ALTER SESSION FORCE parallel dml"

    query_insert=f"INSERT INTO INTERPRO.MATCH_TMP \
        SELECT /*+ parallel */ PS.PROTEIN_AC PROTEIN_AC, \
                    IPR.METHOD_AC METHOD_AC, \
                    IPR.SEQ_START POS_FROM, \
                    IPR.SEQ_END POS_TO, \
                    'T' STATUS, \
                    I2D.DBCODE DBCODE, \
                    I2D.EVIDENCE EVIDENCE, \
                    SYSDATE SEQ_DATE, \
                    SYSDATE MATCH_DATE, \
                    SYSDATE TIMESTAMP, \
                    'INTERPRO' USERSTAMP, \
                    IPR.EVALUE SCORE, \
                    IPR.MODEL_AC MODEL_AC, \
                    IPR.FRAGMENTS FRAGMENTS \
        FROM INTERPRO.PROTEIN_TO_SCAN PS, IPRSCAN.MV_IPRSCAN_MINI IPR, INTERPRO.IPRSCAN2DBCODE I2D \
        WHERE PS.UPI = IPR.UPI AND I2D.IPRSCAN_SIG_LIB_REL_ID = IPR.ANALYSIS_ID \
        AND I2D.DBCODE NOT IN ('g', 'j', 'n', 'q', 's', 'v', 'x') AND I2D.DBCODE IN ('{','.join(dbcode_list)}') \
        AND IPR.SEQ_START != IPR.SEQ_END"

    cursor.execute(query_part)
    cursor.execute(query_alter)

    logger.info("Inserting data in MATCH_TMP")
    cursor.execute(query_insert)
    connection.commit()

def add_indexes(cursor, connection, table, shorttable):
    logger.info("building indices")
    for index, val in {
        "DBCODE": f"I_FK_{shorttable}$DBCODE",
        "EVIDENCE": f"I_FK_{shorttable}$EVIDENCE",
        "STATUS": f"I_{shorttable}$FK_STATUS",
        "METHOD_AC": f"I_{shorttable}$METHOD_AC",
    }.items():
        cursor.execute(
            f"CREATE INDEX INTERPRO.{val} ON INTERPRO.{table}({index}) NOLOGGING"
        )
    connection.commit()

def add_constraints(cursor, connection, table, shorttable):
    cursor.execute(f"ALTER TABLE INTERPRO.{table} ADD CONSTRAINT CK_{shorttable}$FROM CHECK (pos_from >= 1) PARALLEL")
    cursor.execute(f"ALTER TABLE INTERPRO.{table} ADD CONSTRAINT CK_{shorttable}$NEG CHECK (pos_to - pos_from > 0) PARALLEL")
    cursor.execute(f"ALTER TABLE INTERPRO.{table} ADD CONSTRAINT CK_{shorttable}$STATUS_N_PROSITE CHECK (status<>'N' OR (status ='N' AND dbcode IN ('P','M','Q'))) PARALLEL")
    cursor.execute(f"ALTER TABLE INTERPRO.{table} ADD CONSTRAINT PK_{shorttable}$P$M$TO$FM PRIMARY KEY (PROTEIN_AC, METHOD_AC, POS_FROM, POS_TO) PARALLEL")
    cursor.execute(f"ALTER TABLE INTERPRO.{table} ADD CONSTRAINT FK_{shorttable}$DBCODE FOREIGN KEY (DBCODE) REFERENCES INTERPRO.CV_DATABASE (DBCODE) PARALLEL")
    cursor.execute(f"ALTER TABLE INTERPRO.{table} ADD CONSTRAINT FK_{shorttable}$EVIDENCE FOREIGN KEY (EVIDENCE) REFERENCES INTERPRO.CV_EVIDENCE (CODE) PARALLEL")
    cursor.execute(f"ALTER TABLE INTERPRO.{table} ADD CONSTRAINT FK_{shorttable}$METHOD FOREIGN KEY (METHOD_AC) REFERENCES INTERPRO.METHOD (METHOD_AC) ON DELETE CASCADE PARALLEL")
    cursor.execute(f"ALTER TABLE INTERPRO.{table} ADD CONSTRAINT FK_{shorttable}$PROTEIN FOREIGN KEY (PROTEIN_AC) REFERENCES INTERPRO.PROTEIN (PROTEIN_AC) ON  DELETE CASCADE PARALLEL")
    cursor.execute(f"ALTER TABLE INTERPRO.{table} ADD CONSTRAINT FK_{shorttable}$STATUS FOREIGN KEY (STATUS) REFERENCES INTERPRO.CV_STATUS (CODE) PARALLEL")
    cursor.execute(f"ALTER TABLE INTERPRO.{table} ADD CONSTRAINT CK_{shorttable}$PROTEIN_AC CHECK (PROTEIN_AC IS NOT NULL) PARALLEL")
    cursor.execute(f"ALTER TABLE INTERPRO.{table} ADD CONSTRAINT CK_{shorttable}$METHOD_AC CHECK (METHOD_AC IS NOT NULL) PARALLEL")
    cursor.execute(f"ALTER INDEX PK_{shorttable}$P$M$TO$FM RENAME TO UQ_{shorttable}")

    connection.commit()

def delete_duplicate_match(cursor, connection):
    logger.info("SUPERFAMILY: deleting duplicated matches")
    cursor.execute(
        """DELETE FROM INTERPRO.MATCH_TMP M1
    WHERE EXISTS(
        SELECT 1
        FROM INTERPRO.MATCH_TMP M2
        WHERE M2.DBCODE = 'Y'
        AND M1.PROTEIN_AC = M2.PROTEIN_AC
        AND M1.METHOD_AC = M2.METHOD_AC
        AND M1.POS_FROM = M2.POS_FROM
        AND M1.POS_TO = M2.POS_TO
        AND M1.SCORE > M2.SCORE
    )
    """
    )
    logger.info(f"SUPERFAMILY: {cursor.rowcount} rows deleted")
    connection.commit()

def check_match_length(cursor, mail_interpro):
    sqlLine = """SELECT M.METHOD_AC, M.PROTEIN_AC FROM INTERPRO.PROTEIN P, INTERPRO.MATCH_TMP M
                WHERE P.PROTEIN_AC = M.PROTEIN_AC
                AND M.POS_TO > P.LEN"""
    try:
        cursor.execute(sqlLine)
        results = []
        for row in cursor:
            method_ac, protein_ac = row
            result = ""
            result = ",".join([method_ac, protein_ac])
            results.append(result)
        if len(results) > 0:
            proteinupdate.sendmail.send_mail(
            to_addrs=mail_interpro,
            subject="FAILED Check match length",
            content="""\nMatches that are longer than the protein length: {} \
                \nPlease check if matches at iprscan.mv_iprscan_mini are correct.""".format("\n".join(results))
            )
            sys.exit(1)
        else:
            logger.info("All matches are within the protein lengths.")
    except cx_Oracle.DatabaseError as exception:
        logger.info("Failed to execute " + sqlLine)
        logger.info(exception)
        exit(1)

def get_count(cursor, sqlLine):
    try:
        cursor.execute(sqlLine)
        results = []
        dbcodes = []
        for row in cursor:
            dbcode, count = row
            result = ""
            result = ": ".join([dbcode, str(count)])
            results.append(result)
            dbcodes.append(dbcode)
        # logger.info("dbcode: counts")

        return results
    except cx_Oracle.DatabaseError as exception:
        logger.info("Failed to execute " + sqlLine)
        logger.info(exception)
        exit(1)

# def get_protein_count(cursor, table, dbcode):
#     query_count_match=f"select count(distinct mn.protein_ac) from {table} mn, protein p \
#                     where mn.protein_ac = p.protein_ac \
#                     and p.dbcode = 'S' \
#                     and p.fragment = 'N' \
#                     and mn.method_ac =:dbcode
#     cursor.execute(query_count_match, dbcode=dbcode)


def produce_initial_report(cursor, dbcodes, outputdir):

    query = """SELECT METHOD, NEW_COUNT, OLD_COUNT, IPR, PCENT_CHANGE
                FROM (
                SELECT PC.MAC METHOD
                    ,CASE WHEN PC.NEW != 0 THEN  TO_CHAR(PC.NEW) ELSE 'DELETED' END NEW_COUNT
                    ,CASE WHEN PC.OLD != 0 THEN TO_CHAR(PC.OLD) ELSE 'NEW' END OLD_COUNT
                    ,(SELECT NVL(MAX(E2M.ENTRY_AC), 'NI')
                        FROM INTERPRO.ENTRY2METHOD E2M
                        WHERE E2M.METHOD_AC = PC.MAC) IPR
                    ,CASE WHEN PC.OLD != 0
                    THEN TO_CHAR(PC.NEW*100/PC.OLD, '999.99')||'%'
                    ELSE 'NOT DEFINED' END PCENT_CHANGE
                FROM (SELECT NVL(MN.METHOD_AC,ME.METHOD_AC) MAC, NVL(MN.CNT,'0') NEW, NVL(ME.CNT,'0') OLD
                        FROM (SELECT METHOD_AC,COUNT(DISTINCT PROTEIN_AC) CNT
                                FROM INTERPRO.MATCH_TMP
                            WHERE DBCODE =:DBCODE
                            GROUP BY METHOD_AC) MN
                    FULL OUTER JOIN (SELECT METHOD_AC,COUNT(DISTINCT PROTEIN_AC) CNT
                                        FROM INTERPRO.MATCH
                                        WHERE DBCODE =:DBCODE
                                        GROUP BY METHOD_AC) ME ON MN.METHOD_AC = ME.METHOD_AC) PC
                )
            """
    for dbcode in dbcodes:
        try:
            cursor.execute(query, dbcode=dbcode)
            lines = []
            for row in cursor:
                method, new_count, old_count, entry, p_change = row
                line = "\t".join([method, new_count, old_count, entry, p_change])
                lines.append(line)
            write_report(dbcode, lines, outputdir)
        except cx_Oracle.DatabaseError as exception:
            logger.info("Failed to execute " + query)
            logger.info(exception)
            exit(1)

def write_report(dbcode, lines, outputdir):
    """
    Check MATCH_TMP against MATCH table
    """
    file_name = os.path.join(outputdir, f"match_counts_new_{dbcode}.csv")
    header = "\t".join(["Method", "New_Count", "Old_Count", "Entry", "%_change\n"])
    lines = "\n".join(lines)
    with open(file_name, "w") as output:
        output.write(header)
        output.write(lines)
        output.write("\n")


def refresh_site_matches(user: str, dsn: str, memberdb: list):

    dbcodes = []
    for member in memberdb:
        if member["dbcode"] == "J" or member["dbcode"] == "B":
            dbcodes.append(member["dbcode"])

    if len(dbcodes) > 0:
        con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
        cur = con.cursor()

        logger.info("Refreshing SITE_MATCH_NEW")
        orautils.drop_table(cur, "INTERPRO", "SITE_MATCH_NEW", purge=True)

        cur.execute(
            f"CREATE TABLE INTERPRO.SITE_MATCH_NEW \
                    AS \
                    SELECT \
                    P.PROTEIN_AC, S.METHOD_AC, S.LOC_START, S.LOC_END, S.DESCRIPTION, \
                    S.RESIDUE, S.RESIDUE_START, S.RESIDUE_END, S.NUM_SITES, D.DBCODE \
                    FROM INTERPRO.PROTEIN_TO_SCAN P \
                    INNER JOIN IPRSCAN.SITE S \
                    ON P.UPI=S.UPI \
                    INNER JOIN INTERPRO.IPRSCAN2DBCODE D \
                    ON S.ANALYSIS_ID=D.IPRSCAN_SIG_LIB_REL_ID \
                    WHERE D.DBCODE in ('{', '.join(dbcodes)}')"
        )

        logger.info("indexing SITE_MATCH_NEW")
        cur.execute(
            """CREATE INDEX I_SITE_NEW
            ON INTERPRO.SITE_MATCH_NEW
            (PROTEIN_AC, METHOD_AC, LOC_START, LOC_END)
            NOLOGGING
            """
        )

        logger.info("adding partition to SITE_MATCH_NEW")
        cur.execute(
            """ALTER TABLE SITE_MATCH_NEW ADD PARTITION BY LIST ("DBCODE")
            (
                PARTITION "CDD" VALUES ('J'),
                PARTITION "SFLD" VALUES ('B')
            )
            """
        )
        con.commit()

        cur.execute("SELECT DBCODE, COUNT(*) FROM INTERPRO.SITE_MATCH_NEW GROUP BY DBCODE")
        rows = cur.fetchall()
        logger.info("Count in SITE_MATCH_NEW:")
        for row in rows:
            logger.info(f"{row[0]}: row[1]")

        cur.execute("SELECT DBCODE, COUNT(*) FROM INTERPRO.SITE_MATCH GROUP BY DBCODE")
        rows = cur.fetchall()
        logger.info("Count in SITE_MATCH before refresh:")
        for row in rows:
            logger.info(f"{row[0]}: row[1]")

        # exchange partition between site_match_new and site_match
        logger.info("Exchanging partition between SITE_MATCH_NEW and SITE_MATCH")
        for dbcode in dbcodes:
            cur.execute(
                f"ALTER TABLE SITE_MATCH EXCHANGE PARTITION {dbcode} WITH TABLE SITE_MATCH_NEW"
            )
        con.commit()

        cur.execute("SELECT DBCODE, COUNT(*) FROM INTERPRO.SITE_MATCH GROUP BY DBCODE")
        rows = cur.fetchall()
        logger.info("Count in SITE_MATCH after refresh:")
        for row in rows:
            logger.info(f"{row[0]}: row[1]")

        #orautils.drop_table(cur, "INTERPRO", "SITE_MATCH_NEW", purge=True)

        cur.close()
        con.close()
    else:
        logger.info("No need to update SITE_MATCH table")

def drop_match_tmp(user: str, dsn: str):
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()

    logger.info("dropping MATCH_TMP")
    orautils.drop_table(cur, "INTERPRO", "MATCH_TMP", purge=True)
    con.commit()

    cur.close()
    con.close()


def create_match_tmp(user: str, dsn: str, memberdb: list, outputdir: str, mail_interpro):
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()

    dbcodes = []
    results_string=""
    for member in memberdb:
        dbcodes.append(member["dbcode"])
    logger.info(dbcodes)

    #not working
    logger.info("Creating match_tmp")
    start_time = time.time()
    prepare_matches(cur, con, dbcodes)
    logger.info(f"Elapsed time in seconds: {time.time() - start_time}")

    if 'Y' in dbcodes:
        delete_duplicate_match(cur, con)
    try:
        add_indexes(cur, con, "MATCH_TMP", "MATCHT")
    except cx_Oracle.DatabaseError as e:
        logger.info(f"Error adding indexes to MATCH_TMP {e}")
        sys.exit(1)

    try:
        add_constraints(cur, con, "MATCH_TMP", "MATCHT")
    except cx_Oracle.DatabaseError as e:
        logger.info(f"Error adding constraints to MATCH_TMP {e}")
        sys.exit(1)

    check_match_length(cur, mail_interpro)

    get_distinct_method_count = "SELECT DBCODE, COUNT(DISTINCT METHOD_AC) FROM INTERPRO.MATCH_TMP GROUP BY DBCODE"
    results_string+="Please check these numbers correspond to the method counts in method table:\n"
    
    results = get_count(cur, get_distinct_method_count)
    results_string+='\n'.join(results)

    produce_initial_report(cur, dbcodes, outputdir)
    results_string+=f"\n\nPlease check the {outputdir}/match_counts_new_*.csv files.\n"

    get_match_tmp_count = (
        "SELECT DBCODE, COUNT(*) FROM INTERPRO.MATCH_TMP GROUP BY DBCODE"
    )
    results_string+="\nThe total match_tmp count for each member database is:\n"
    results = get_count(cur, get_match_tmp_count)
    results_string+='\n'.join(results)

    get_match_count = f"SELECT DBCODE, COUNT(*) FROM INTERPRO.MATCH WHERE DBCODE IN ('{','.join(dbcodes)}') GROUP BY DBCODE"
    results_string+="\n\nThe total match count for each member database is:\n"
    results = get_count(cur, get_match_count)
    results_string+='\n'.join(results)

    cur.close()
    con.close()

    logger.info(
        "MATCH_TMP is ready. Please check the above to make sure you have the correct data."
    )
    proteinupdate.sendmail.send_mail(
        to_addrs=mail_interpro,
        subject="MATCH TMP",
        content=f"Counts for Create match_tmp\n {results_string}"
    )

