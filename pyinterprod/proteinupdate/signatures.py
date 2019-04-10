import cx_Oracle

from ..orautils import make_connect_string


def update_method2descriptions(user: str, dsn: str):
    con = cx_Oracle.connect(make_connect_string(user, dsn))
    cur = con.cursor()
    cur.execute("TRUNCATE TABLE INTERPRO.METHOD2SWISS_DE")
    cur.execute(
        """
        INSERT /*+APPEND*/ INTO INTERPRO.METHOD2SWISS_DE
          SELECT M.PROTEIN_AC, M.METHOD_AC, D.TEXT
          FROM INTERPRO_ANALYSIS_LOAD.METHOD2PROTEIN M
          INNER JOIN INTERPRO_ANALYSIS_LOAD.DESC_VALUE D
            ON M.DESC_ID = D.DESC_ID
          WHERE M.DBCODE = 'S'
        """
    )
    con.commit()
    cur.close()
    con.close()
