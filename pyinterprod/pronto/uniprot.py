# -*- coding: utf-8 -*-

import math

import cx_Oracle

from .. import orautils


def load_comments(user: str, dsn: str):
    tables = ("CV_COMMENT_TOPIC", "COMMENT_VALUE", "PROTEIN_COMMENT")
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    for table in tables:
        orautils.drop_table(cur, owner, table)

    cur.execute(
        """
        CREATE TABLE {}.CV_COMMENT_TOPIC
        (
            TOPIC_ID NUMBER(2) NOT NULL,
            TOPIC VARCHAR2(30) NOT NULL
        ) NOLOGGING
        """.format(owner)
    )

    cur.execute(
        """
        CREATE TABLE {}.COMMENT_VALUE
        (
            TOPIC_ID NUMBER(2) NOT NULL,
            COMMENT_ID NUMBER(6) NOT NULL,
            TEXT VARCHAR2(4000) NOT NULL,
            CONSTRAINT PK_COMMENT_VALUE PRIMARY KEY (TOPIC_ID, COMMENT_ID)
        ) NOLOGGING
        """.format(owner)
    )

    cur.execute(
        """
        CREATE TABLE {}.PROTEIN_COMMENT
        (
            PROTEIN_AC VARCHAR2(15) NOT NULL,
            TOPIC_ID NUMBER(2) NOT NULL,
            COMMENT_ID NUMBER(6) NOT NULL,
            CONSTRAINT PK_PROTEIN_COMMENT PRIMARY KEY (
                PROTEIN_AC, TOPIC_ID, COMMENT_ID
            )
        ) NOLOGGING
        """.format(owner)
    )

    cur.execute(
        """
        SELECT E.ACCESSION, B.COMMENT_TOPICS_ID, T.TOPIC, SS.TEXT
        FROM SPTR.DBENTRY@SWPREAD E
        INNER JOIN SPTR.COMMENT_BLOCK@SWPREAD B
            ON E.DBENTRY_ID = B.DBENTRY_ID
        INNER JOIN SPTR.CV_COMMENT_TOPICS@SWPREAD T
            ON B.COMMENT_TOPICS_ID = T.COMMENT_TOPICS_ID
        INNER JOIN SPTR.COMMENT_STRUCTURE@SWPREAD S
            ON B.COMMENT_BLOCK_ID = S.COMMENT_BLOCK_ID
        INNER JOIN SPTR.COMMENT_SUBSTRUCTURE@SWPREAD SS
            ON S.COMMENT_STRUCTURE_ID = SS.COMMENT_STRUCTURE_ID
        WHERE E.ENTRY_TYPE = 0
        AND E.MERGE_STATUS != 'R'
        AND E.DELETED = 'N'
        AND E.FIRST_PUBLIC IS NOT NULL
        """
    )

    topics = {}
    protein2comment = set()
    for row in cur:
        accession = row[0]
        topic_id = math.floor(row[1])
        comment = row[3]

        if topic_id not in topics:
            topics[topic_id] = {
                "name": row[2],
                "comments": {}
            }

        topic_comments = topics[topic_id]["comments"]
        if comment in topic_comments:
            comment_id = topic_comments[comment]
        else:
            comment_id = topic_comments[comment] = len(topic_comments) + 1

        protein2comment.add((accession, topic_id, comment_id))

    comments = []
    for topic_id, topic in topics.items():
        for comment, comment_id in topic["comments"].items():
            comments.append((topic_id, comment_id, comment))

    topics = [
        (topic_id, topic["name"])
        for topic_id, topic in topics.items()
    ]

    cur.executemany(
        """
        INSERT /*+ APPEND */ INTO {}.CV_COMMENT_TOPIC
        VALUES (:1, :2)
        """.format(owner), topics
    )
    con.commit()

    table = orautils.TablePopulator(con,
                                    query="INSERT /*+ APPEND */ "
                                          "INTO {}.COMMENT_VALUE "
                                          "VALUES (:1, :2, :3)".format(owner),
                                    autocommit=True)
    for record in comments:
        table.insert(record)
    table.close()

    table = orautils.TablePopulator(con,
                                    query="INSERT /*+ APPEND */ "
                                          "INTO {}.PROTEIN_COMMENT "
                                          "VALUES (:1, :2, :3)".format(owner),
                                    autocommit=True)
    for record in protein2comment:
        table.insert(record)
    table.close()

    for table in tables:
        orautils.gather_stats(cur, owner, table)
        orautils.grant(cur, owner, table, "SELECT", "INTERPRO_SELECT")
    cur.close()
    con.close()


def load_descriptions(user: str, dsn: str):
    tables = ["DESC_VALUE", "PROTEIN_DESC"]
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    for table in tables:
        orautils.drop_table(cur, owner, table)

    cur.execute(
        """
        CREATE TABLE {}.DESC_VALUE
        (
            DESC_ID NUMBER(10) NOT NULL,
            TEXT VARCHAR2(4000) NOT NULL
        ) NOLOGGING
        """.format(owner)
    )

    cur.execute(
        """
        CREATE TABLE {}.PROTEIN_DESC
        (
            PROTEIN_AC VARCHAR2(15) NOT NULL,
            DESC_ID NUMBER(10) NOT NULL
        ) NOLOGGING
        """.format(owner)
    )

    cur.execute(
        """
        SELECT DESCR, ACCESSION
        FROM (
          SELECT
            E.ACCESSION,
            D.DESCR,
            ROW_NUMBER() OVER (
                PARTITION BY E.ACCESSION
                ORDER BY D.SECTION_GROUP_ID, D.DESC_ID
            ) R
          FROM SPTR.DBENTRY@SWPREAD E
          LEFT OUTER JOIN SPTR.DBENTRY_2_DESC@SWPREAD D
            ON E.DBENTRY_ID = D.DBENTRY_ID
          WHERE E.ENTRY_TYPE IN (0, 1)
                AND E.MERGE_STATUS != 'R'
                AND E.DELETED = 'N'
                AND E.FIRST_PUBLIC IS NOT NULL
        )
        WHERE R = 1        
        """
    )

    table1 = orautils.TablePopulator(con,
                                     query="INSERT INTO {}.DESC_VALUE "
                                           "VALUES (:1, :2)".format(owner))

    table2 = orautils.TablePopulator(con,
                                     query="INSERT /*+ APPEND */ "
                                           "INTO {}.PROTEIN_DESC "
                                           "VALUES (:1, :2)".format(owner),
                                     autocommit=True)
    descriptions = {}
    for text, accession in cur:
        if text in descriptions:
            desc_id = descriptions[text]
        else:
            desc_id = descriptions[text] = len(descriptions) + 1
            table1.insert((desc_id, text))
            table1.flush(commit=True)

        table2.insert((accession, desc_id))

    table1.close()
    table2.close()

    cur.execute(
        """
        ALTER TABLE {}.DESC_VALUE
        ADD CONSTRAINT PK_DESC_VALUE PRIMARY KEY (DESC_ID)
        """.format(owner)
    )

    cur.execute(
        """
        ALTER TABLE {}.PROTEIN_DESC
        ADD CONSTRAINT PK_PROTEIN_DESC PRIMARY KEY (PROTEIN_AC, DESC_ID)
        """.format(owner)
    )

    for table in tables:
        orautils.gather_stats(cur, owner, table)
        orautils.grant(cur, owner, table, "SELECT", "INTERPRO_SELECT")

    cur.close()
    con.close()


def load_enzymes(user: str, dsn: str):
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    orautils.drop_table(cur, owner, "ENZYME")
    cur.execute(
        """
        CREATE TABLE {}.ENZYME
        (
            PROTEIN_AC VARCHAR2(15) NOT NULL,
            ECNO VARCHAR2(15) NOT NULL
        ) NOLOGGING
        """.format(owner)
    )

    cur.execute(
        """
        INSERT /*+APPEND*/ INTO {}.ENZYME (PROTEIN_AC, ECNO)
        SELECT DISTINCT E.ACCESSION, D.DESCR
        FROM SPTR.DBENTRY@SWPREAD E
        LEFT OUTER JOIN SPTR.DBENTRY_2_DESC@SWPREAD D
            ON E.DBENTRY_ID = D.DBENTRY_ID
        LEFT OUTER JOIN SPTR.CV_DESC@SWPREAD C
            ON D.DESC_ID = C.DESC_ID
        WHERE E.ENTRY_TYPE IN (0, 1)
        AND E.MERGE_STATUS != 'R'
        AND E.DELETED = 'N'
        AND E.FIRST_PUBLIC IS NOT NULL
        AND C.SUBCATG_TYPE='EC'
        """.format(owner)
    )
    con.commit()

    cur.execute(
        """
        CREATE INDEX I_ENZYME$PROTEIN
        ON {}.ENZYME (PROTEIN_AC) NOLOGGING
        """.format(owner)
    )

    cur.execute(
        """
        CREATE INDEX I_ENZYME$EC
        ON {}.ENZYME (ECNO) NOLOGGING
        """.format(owner)
    )

    orautils.gather_stats(cur, owner, "ENZYME")
    orautils.grant(cur, owner, "ENZYME", "SELECT", "INTERPRO_SELECT")
    cur.close()
    con.close()
