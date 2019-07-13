# -*- coding: utf-8 -*-

import cx_Oracle

from ..orautils import drop_table, make_connect_string


def update_method2descriptions(user: str, dsn: str):
    con = cx_Oracle.connect(make_connect_string(user, dsn))
    cur = con.cursor()
    drop_table(cur, "INTERPRO", "METHOD2SWISS_DE", purge=True)
    cur.execute(
        """
        CREATE TABLE INTERPRO.METHOD2SWISS_DE
        NOLOGGING
        AS
        SELECT M.PROTEIN_AC, M.METHOD_AC, D.TEXT
        FROM INTERPRO_ANALYSIS_LOAD.METHOD2PROTEIN PARTITION(M2P_SWISSP) M
        INNER JOIN INTERPRO_ANALYSIS_LOAD.DESC_VALUE D
          ON M.DESC_ID = D.DESC_ID
        """
    )

    cur.execute(
        """
        ALTER TABLE INTERPRO.METHOD2SWISS_DE
        ADD CONSTRAINT PK_METHOD2SWISS
        PRIMARY KEY (PROTEIN_AC, METHOD_AC)
        """
    )

    cur.close()
    con.close()
