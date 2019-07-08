# -*- coding: utf-8 -*-

import hashlib
import math
from multiprocessing import Queue
from typing import List, Optional

import cx_Oracle

from .. import orautils
from . import utils


_MAX_GAP = 20        # at least 20 residues between positions


def load_comments(user: str, dsn: str):
    tables = ("CV_COMMENT_TOPIC", "COMMENT_VALUE", "PROTEIN_COMMENT")
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    for table in tables:
        orautils.drop_table(cur, owner, table, purge=True)

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


def load_descriptions_ctas(user: str, dsn: str):
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    orautils.drop_table(cur, owner, "PROTEIN_DESC_STG", purge=True)
    cur.execute(
        """
        CREATE TABLE {}.PROTEIN_DESC_STG
        (
          TEXT VARCHAR2(4000) NOT NULL,
          PROTEIN_AC VARCHAR2(15) NOT NULL,
          DESC_ID NUMBER(10) NOT NULL
        ) NOLOGGING
        AS
        SELECT DESCR, ACCESSION, DENSE_RANK() OVER (ORDER BY DESCR)
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
          WHERE E.ENTRY_TYPE = 0
            AND E.MERGE_STATUS != 'R'
            AND E.DELETED = 'N'
            AND E.FIRST_PUBLIC IS NOT NULL
        )
        WHERE R = 1
        """
    )

    for table in ("DESC_VALUE", "PROTEIN_DESC"):
        orautils.drop_table(cur, owner, table, purge=True)

    cur.execute(
        """
        CREATE TABLE {0}.DESC_VALUE
        NOLOGGING
        AS
        SELECT DESC_ID, TEXT
        FROM {0}.PROTEIN_DESC_STG
        """.format(owner)
    )

    cur.execute(
        """
        CREATE TABLE {0}.PROTEIN_DESC
        NOLOGGING
        AS
        SELECT PROTEIN_AC, DESC_ID
        FROM {0}.PROTEIN_DESC_STG
        """.format(owner)
    )

    orautils.drop_table(cur, owner, "PROTEIN_DESC_STG", purge=True)

    cur.execute(
        """
        CREATE UNIQUE INDEX UI_DESC_VALUE
        ON {}.DESC_VALUE (DESC_ID) NOLOGGING
        """.format(owner)
    )

    cur.execute(
        """
        CREATE UNIQUE INDEX UI_PROTEIN_DESC
        ON {}.PROTEIN_DESC (PROTEIN_AC) NOLOGGING
        """.format(owner)
    )
    cur.execute(
        """
        CREATE INDEX I_PROTEIN_DESC
        ON {}.PROTEIN_DESC (DESC_ID) NOLOGGING
        """.format(owner)
    )

    for table in ("DESC_VALUE", "PROTEIN_DESC"):
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
        orautils.drop_table(cur, owner, table, purge=True)

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
                                     query="INSERT /*+ APPEND */ "
                                           "INTO {}.DESC_VALUE "
                                           "VALUES (:1, :2)".format(owner),
                                     autocommit=True)

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

        table2.insert((accession, desc_id))

    table1.close()
    table2.close()

    cur.execute(
        """
        CREATE UNIQUE INDEX UI_DESC_VALUE
        ON {}.DESC_VALUE (DESC_ID) NOLOGGING
        """.format(owner)
    )

    cur.execute(
        """
        CREATE UNIQUE INDEX UI_PROTEIN_DESC
        ON {}.PROTEIN_DESC (PROTEIN_AC) NOLOGGING
        """.format(owner)
    )
    cur.execute(
        """
        CREATE INDEX I_PROTEIN_DESC
        ON {}.PROTEIN_DESC (DESC_ID) NOLOGGING
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
    orautils.drop_table(cur, owner, "ENZYME", purge=True)
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


def load_proteins(user: str, dsn: str):
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    orautils.drop_table(cur, owner, "PROTEIN", purge=True)
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


def consume_proteins(user: str, dsn: str, task_queue: Queue,
                     done_queue: Queue, tmpdir: Optional[str]=None,
                     bucket_size: int=100):
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    cur.execute("SELECT METHOD_AC FROM {}.METHOD".format(owner))
    keys = sorted([row[0] for row in cur])
    keys = [keys[i] for i in range(0, len(keys), bucket_size)]

    proteins2go = {}
    cur.execute("SELECT PROTEIN_AC, GO_ID FROM {}.PROTEIN2GO".format(owner))
    for acc, go_id in cur:
        if acc in proteins2go:
            proteins2go[acc].add(go_id)
        else:
            proteins2go[acc] = {go_id}
    cur.close()

    names = utils.Organizer(keys, dir=tmpdir)
    taxa = utils.Organizer(keys, dir=tmpdir)
    terms = utils.Organizer(keys, dir=tmpdir)
    comparator = utils.MatchComparator(dir=tmpdir)
    table = orautils.TablePopulator(
        con=con,
        query="INSERT /*+ APPEND */ "
              "INTO {}.METHOD2PROTEIN "
              "VALUES (:1, :2, :3, :4, :5, :6, :7)".format(owner),
        autocommit=True
    )
    for chunk in iter(task_queue.get, None):
        for acc, dbcode, length, tax_id, desc_id, matches in chunk:
            md5 = _hash_protein(matches)
            protein_terms = proteins2go.get(acc, [])

            signatures = comparator.update(matches)

            # Update organizers and populate table
            for signature_acc in signatures:
                # UniProt descriptions
                names.add(signature_acc, (desc_id, dbcode))
                # Taxonomic origins
                taxa.add(signature_acc, tax_id)
                # GO terms
                for go_id in protein_terms:
                    terms.add(signature_acc, go_id)

                table.insert((signature_acc, acc, dbcode, md5, length,
                              tax_id, desc_id))

        names.flush()
        taxa.flush()
        terms.flush()

    table.close()
    con.close()

    comparator.sync()
    size = comparator.size
    organizers = (names, taxa, terms)
    for o in organizers:
        size += o.merge()

    done_queue.put((comparator, *organizers, size))


def _hash_protein(matches: List[tuple]) -> str:
    # flatten matches
    locations = []
    for acc, pos_start, pos_end in matches:
        locations.append((pos_start, acc))
        locations.append((pos_end, acc))

    """
    Evaluate the protein's match structure,
        i.e. how signatures match the proteins

    -----------------------------   Protein
     <    >                         Signature 1
       <    >                       Signature 2
                  < >               Signature 3

    Flattened:
    -----------------------------   Protein
     < <  > >     < >
     1 2  1 2     3 3

    Structure, with '-' representing a "gap"
        (more than N bp between two positions):
    1212-33
    """

    # Sort locations by position
    locations.sort()

    """
    Do not set the offset to 0, but to the first position:
    if two proteins have the same structure,
    but the first position of one protein is > max_gap
    while the first position of the other protein is <= max_gap,
    a gap will be used for the first protein and not for the other,
    which will results in two different structures
    """
    offset = locations[0][0]

    # overall match structure
    structure = []
    # close signatures (less than max_gap between two positions)
    signatures = []

    for pos, acc in locations:
        if pos > offset + _MAX_GAP:
            for _acc in signatures:
                structure.append(_acc)

            signatures = []
            structure.append('')  # add a gap

        offset = pos
        signatures.append(acc)

    for _acc in signatures:
        structure.append(_acc)

    return hashlib.md5(
        '/'.join(structure).encode("utf-8")
    ).hexdigest()
