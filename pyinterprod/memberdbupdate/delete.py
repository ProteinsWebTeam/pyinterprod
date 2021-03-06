from .. import orautils, logger
import cx_Oracle


def create_temp_table(cur, con):
    try:
        query = """CREATE GLOBAL TEMPORARY TABLE METHODS_TO_DELETE_TEMP_TABLE (
                METHOD_AC  VARCHAR2(25 BYTE)
                ) ON COMMIT PRESERVE ROWS"""

        cur.execute(query)
        con.commit()
    except Exception as e:
        logger.info(e)

def delete_from_tables(cur, con, dbcode:str):

    query_copy = """INSERT INTO INTERPRO.METHODS_TO_DELETE_TEMP_TABLE
            (SELECT METHOD_AC
            FROM INTERPRO.METHOD
            WHERE DBCODE=:dbcode
            MINUS
            SELECT METHOD_AC
            FROM INTERPRO.METHOD_STG
            WHERE DBCODE=:dbcode)
            """

    cur.execute(query_copy, dbcode=dbcode)
    logger.info(f"\ninserted rows: + {cur.rowcount}")
    con.commit()

    query = "SELECT * FROM INTERPRO.METHODS_TO_DELETE_TEMP_TABLE"
    cur.execute(query)
    results = [row[0] for row in cur]
    logger.info(f"rows recovered: {results}")

    query_partition = f"DELETE FROM MATCH PARTITION (MATCH_DBCODE_{dbcode}) WHERE METHOD_AC IN (SELECT METHOD_AC FROM METHODS_TO_DELETE_TEMP_TABLE)"
    cur.execute(query_partition)
    con.commit()

    tables_list = [
        "ENTRY2METHOD",
        "MATCH_NEW",
        "METHOD2PUB",
        "MV_METHOD2PROTEIN",
        "MV_METHOD_MATCH",
        "VARSPLIC_MATCH",
        "METHOD",
    ]
    for table in tables_list:
        logger.info(table)
        query = f"DELETE FROM {table} WHERE METHOD_AC IN (SELECT METHOD_AC FROM METHODS_TO_DELETE_TEMP_TABLE)"
        cur.execute(query)
        con.commit()
        logger.info(cur.rowcount)


def delete_dead_signatures(user: str, dsn: str, memberdb: list):
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()

    create_temp_table(cur, con)

    for member in memberdb:
        delete_from_tables(cur, con, member["dbcode"])

    cur.close()
    con.close()
