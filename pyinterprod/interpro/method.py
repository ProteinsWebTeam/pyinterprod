# -*- coding: utf-8 -*-

import cx_Oracle

from pyinterprod.utils import oracle


def update_method2descriptions(url: str):
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    oracle.truncate_table(cur, "INTERPRO.METHOD2SWISS_DE", reuse_storage=True)
    cur.execute(
        """
        INSERT INTO INTERPRO.METHOD2SWISS_DE
        SELECT M.PROTEIN_AC, M.METHOD_AC, D.TEXT AS DESCRIPTION
        FROM INTERPRO_ANALYSIS_LOAD.METHOD2PROTEIN PARTITION(M2P_SWISSP) M
        INNER JOIN INTERPRO_ANALYSIS_LOAD.DESC_VALUE D
          ON M.DESC_ID = D.DESC_ID
        """
    )
    cur.execute(
        """
        CREATE TABLE INTERPRO.METHOD2SWISS_DE
        NOLOGGING
        AS
        SELECT M.PROTEIN_AC, M.METHOD_AC, D.TEXT AS DESCRIPTION
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