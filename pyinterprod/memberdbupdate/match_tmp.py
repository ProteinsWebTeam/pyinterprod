#!/usr/bin/env python3

import time
import re
import cx_Oracle
import sys
import os
from .. import logger, orautils


class Recreate_match(object):

    def __init__(self, user, dsn):
        self.user, self.dsn = user, dsn
        url = orautils.make_connect_string(user, dsn)
        self.connection = cx_Oracle.connect(url)
        self.cursor = self.connection.cursor()

    def prepare_matches(self, dbcode_list):
        logger.info("adding matches_tmp")
        orautils.drop_table(self.cursor, "INTERPRO", "MATCH_TMP", purge=True)
        self.cursor.execute(
            """CREATE TABLE INTERPRO.MATCH_TMP NOLOGGING
            AS
            SELECT * FROM INTERPRO.MATCH WHERE 1 = 0
            """
        )

        query_insert = f"INSERT /*+ APPEND */ INTO INTERPRO.MATCH_TMP ( \
          PROTEIN_AC, METHOD_AC, POS_FROM, POS_TO, STATUS, \
          DBCODE, EVIDENCE, \
          SEQ_DATE, MATCH_DATE, TIMESTAMP, USERSTAMP, \
          SCORE, MODEL_AC, FRAGMENTS \
        ) \
        SELECT \
          P.PROTEIN_AC, M.METHOD_AC, M.SEQ_START, M.SEQ_END, 'T', \
          D.DBCODE, D.EVIDENCE, \
          SYSDATE, SYSDATE, SYSDATE, 'INTERPRO', \
          M.EVALUE, M.MODEL_AC, M.FRAGMENTS \
        FROM INTERPRO.PROTEIN_TO_SCAN P \
        INNER JOIN IPRSCAN.MV_IPRSCAN_MINI M \
          ON P.UPI = M.UPI \
        INNER JOIN INTERPRO.IPRSCAN2DBCODE D \
          ON M.ANALYSIS_ID = D.IPRSCAN_SIG_LIB_REL_ID \
        WHERE D.DBCODE NOT IN ('g', 'j', 'n', 'q', 's', 'v', 'x') AND D.DBCODE IN ('{','.join(dbcode_list)}') \
        AND M.SEQ_START != M.SEQ_END "

        self.cursor.execute(query_insert)
        self.connection.commit()

    def add_constraint(self):
        logger.info("building indices")
        for index, val in {"DBCODE": "I_FK_MATCHT$DBCODE", "EVIDENCE": "I_FK_MATCHT$EVIDENCE", "STATUS": "I_MATCHT$FK_STATUS", "METHOD_AC": "I_MATCHT$METHOD_AC"}.items():
            self.cursor.execute(
                f" CREATE INDEX INTERPRO.{val} ON INTERPRO.MATCH_TMP({index}) NOLOGGING")

        self.connection.commit()

    def delete_duplicate_match(self):
        logger.info("SUPERFAMILY: deleting duplicated matches")
        self.cursor.execute(
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
        logger.info(f"SUPERFAMILY: {self.cursor.rowcount} rows deleted")
        self.connection.commit()

    def check_match_length(self):
        sqlLine = """SELECT M.METHOD_AC, M.PROTEIN_AC FROM INTERPRO.PROTEIN P, INTERPRO.MATCH_TMP M
                  WHERE P.PROTEIN_AC = M.PROTEIN_AC
                  AND M.POS_TO > P.LEN"""
        try:
            self.cursor.execute(sqlLine)
            results = []
            for row in self.cursor:
                method_ac, protein_ac = row
                result = ""
                result = ','.join([method_ac, protein_ac])
                results.append(result)
            if len(results) > 0:
                logger.info(
                    "\nMatches that are longer than the protein length: ")
                logger.info("\n".join(results))
                logger.info(
                    "Please check if matches at iprscan.mv_iprscan_mini are correct.")
                sys.exit(1)
            else:
                logger.info("\nAll matches are within the protein lengths.")
        except cx_Oracle.DatabaseError as exception:
            logger.info("Failed to execute " + sqlLine)
            logger.info(exception)
            exit(1)

    def create_match_tmp(self, file_name, id_list):
        analysis_id = ','.join(id_list)
        pattern = '--'
        with open(file_name, 'r') as sqlFile:
            sqlLine = ''
            for line in sqlFile:
                if re.search(pattern, line):
                    continue
                sqlLine += line
                if line.strip().endswith(';'):
                    try:
                        sqlLine = sqlLine.replace(";", "")
                        sqlLine = sqlLine.format(analysis_id)
                        self.cursor.execute(sqlLine)
                        self.connection.commit()
                        sqlLine = ''
                    except cx_Oracle.DatabaseError as exception:
                        logger.info('Failed to execute ' + sqlLine)
                        logger.info(exception)
                        exit(1)

    def execute_sql_stm(self, file_name):
        pattern = '--'
        with open(file_name, 'r') as sqlFile:
            sqlLine = ''
            for line in sqlFile:
                if re.search(pattern, line):
                    continue
                sqlLine += line
                if line.strip().endswith(';'):
                    try:
                        sqlLine = sqlLine.replace(";", "")
                        self.cursor.execute(sqlLine)
                        sqlLine = ''
                    except cx_Oracle.DatabaseError as exception:
                        logger.info("Failed to execute " + sqlLine)
                        logger.info(exception)
                        exit(1)

    def get_count(self, sqlLine):
        try:
            self.cursor.execute(sqlLine)
            results = []
            dbcodes = []
            for row in self.cursor:
                dbcode, count = row
                result = ''
                result = ': '.join([dbcode, str(count)])
                results.append(result)
                dbcodes.append(dbcode)
            # logger.info("dbcode: counts")
            logger.info('\n'.join(results))
            return dbcodes
        except cx_Oracle.DatabaseError as exception:
            logger.info("Failed to execute " + sqlLine)
            logger.info(exception)
            exit(1)

    def produce_initial_report(self, dbcodes, outputdir):

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
                self.cursor.execute(query, dbcode=dbcode)
                lines = []
                for row in self.cursor:
                    method, new_count, old_count, entry, p_change = row
                    line = '\t'.join(
                        [method, new_count, old_count, entry, p_change])
                    lines.append(line)
                self.write_report(dbcode, lines)
            except cx_Oracle.DatabaseError as exception:
                logger.info("Failed to execute " + query)
                logger.info(exception)
                exit(1)

    def write_report(self, dbcode, lines):
        """
        Check MATCH_TMP against MATCH table
        """
        file_name = "match_counts_new_" + dbcode + ".csv"
        header = '\t'.join(
            ['Method', 'New_Count', 'Old_Count', 'Entry', '%_change\n'])
        lines = '\n'.join(lines)
        with open(file_name, 'w') as output:
            output.write(header)
            output.write(lines)
            output.write('\n')

    def close(self):
        self.cursor.close()
        self.connection.close()


def refresh_site_matches(user: str, dsn: str, memberdb: list):
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()

    dbcodes = []
    for member in memberdb:
        if member["dbcode"] == 'J' or member["dbcode"] == 'B':
            dbcodes.append(member["dbcode"])

    logger.info("Refreshing SITE_MATCH_NEW")
    orautils.drop_table(cur, "INTERPRO", "SITE_MATCH_NEW", purge=True)

    cur.execute(f" CREATE TABLE INTERPRO.SITE_MATCH_NEW \
                AS \
                SELECT \
                P.PROTEIN_AC, S.METHOD_AC, S.LOC_START, S.LOC_END, S.DESCRIPTION, \
                S.RESIDUE, S.RESIDUE_START, S.RESIDUE_END, S.NUM_SITES, D.DBCODE \
                FROM INTERPRO.PROTEIN_TO_SCAN P \
                INNER JOIN IPRSCAN.SITE S \
                ON P.UPI=S.UPI \
                INNER JOIN INTERPRO.IPRSCAN2DBCODE D \
                ON S.ANALYSIS_ID=D.IPRSCAN_SIG_LIB_REL_ID \
                WHERE D.DBCODE in ('{', '.join(dbcodes)}') \
                ")

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

    cur.execute(
        "SELECT DBCODE, COUNT(*) FROM INTERPRO.SITE_MATCH_NEW GROUP BY DBCODE"
    )
    rows = cur.fetchall()
    logger.info("Count in SITE_MATCH_NEW:")
    for row in rows:
        logger.info(f"{row[0]}: row[1]")

    cur.execute(
        "SELECT DBCODE, COUNT(*) FROM INTERPRO.SITE_MATCH GROUP BY DBCODE"
    )
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

    cur.execute(
        "SELECT DBCODE, COUNT(*) FROM INTERPRO.SITE_MATCH GROUP BY DBCODE")
    rows = cur.fetchall()
    logger.info("Count in SITE_MATCH after refresh:")
    for row in rows:
        logger.info(f"{row[0]}: row[1]")

    cur.close()
    con.close()


def create_match_tmp(user: str, dsn: str, memberdb: list, outputdir: str):

    recreate_match = Recreate_match(user, dsn)
    dbcodes = []
    for member in memberdb:
        dbcodes.append(member["dbcode"])

    logger.info("Creating match_tmp")
    start_time = time.time()
    recreate_match.prepare_matches(dbcodes)
    logger.info(f"Elapsed time in seconds: {time.time() - start_time}")

    recreate_match.delete_duplicate_match()
    recreate_match.add_constraint()

    recreate_match.check_match_length()

    get_distinct_method_count = "SELECT DBCODE, COUNT(DISTINCT METHOD_AC) FROM INTERPRO.MATCH_TMP GROUP BY DBCODE"
    logger.info(
        "Please check these numbers correspond to the method counts in method table:")
    dbcode_list = recreate_match.get_count(get_distinct_method_count)

    recreate_match.produce_initial_report(dbcodes, outputdir)
    logger.info(
        f"\nPlease check the {outputdir}/match_counts_new_*.csv files.\n")

    get_match_tmp_count = "SELECT DBCODE, COUNT(*) FROM INTERPRO.MATCH_TMP GROUP BY DBCODE"
    logger.info("\nThe total match_tmp count for each member database is:")
    recreate_match.get_count(get_match_tmp_count)

    get_match_count = f"SELECT DBCODE, COUNT(*) FROM INTERPRO.MATCH WHERE DBCODE IN ('{','.join(dbcode_list)}') GROUP BY DBCODE"
    logger.info("\nThe total match count for each member database is:")
    recreate_match.get_count(get_match_count)

    recreate_match.close()
    logger.info(
        "MATCH_TMP is ready. Please check the above to make sure you have the correct data.")


def drop_match_tmp(user: str, dsn: str):
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()

    logger.info("dropping MATCH_TMP")
    cur.execute("DROP TABLE INTERPRO.MATCH_TMP PURGE")
    con.commit()

    cur.close()
    con.close()
