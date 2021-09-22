# -*- coding: utf-8 -*-

import pickle

import cx_Oracle
import psycopg2
from psycopg2.extras import execute_values

from pyinterprod import logger
from pyinterprod.utils.pg import url2dict


def import_signatures(ora_url: str, pg_url: str, allseqs: str, compseqs: str):
    logger.info("populating")
    with open(allseqs, "rb") as fh:
        allseqs = pickle.load(fh)

    with open(compseqs, "rb") as fh:
        compseqs = pickle.load(fh)

    pg_con = psycopg2.connect(**url2dict(pg_url))
    with pg_con.cursor() as pg_cur:
        pg_cur.execute("SELECT name, id FROM database")
        databases = dict(pg_cur.fetchall())

    ora_con = cx_Oracle.connect(ora_url)
    ora_cur = ora_con.cursor()
    ora_cur.execute(
        """
        SELECT
            M.METHOD_AC, LOWER(D.DBSHORT), M.NAME, M.DESCRIPTION, T.ABBREV, 
            M.ABSTRACT, M.ABSTRACT_LONG
        FROM INTERPRO.METHOD M
        INNER JOIN INTERPRO.CV_DATABASE D ON M.DBCODE = D.DBCODE
        INNER JOIN INTERPRO.CV_ENTRY_TYPE T ON M.SIG_TYPE = T.CODE
        """
    )
    values = []
    for row in ora_cur:
        try:
            db_id = databases[row[1]]
        except KeyError:
            continue

        acc = row[0]
        try:
            n_rev_seqs, n_unrev_seqs = allseqs[acc]
        except KeyError:
            n_rev_seqs = n_unrev_seqs = 0

        try:
            num_compl_seqs, num_compl_rev_seqs, num_residues = compseqs[acc]
        except KeyError:
            num_compl_seqs = num_compl_rev_seqs = num_residues = 0

        values.append((
            acc,                        # accession
            db_id,                      # database_id
            row[2],                     # name
            row[3],                     # description
            row[4],                     # type
            row[6].read() if row[6] is not None else row[5],    # abstract
            n_rev_seqs + n_unrev_seqs,  # num_sequences
            n_rev_seqs,                 # num_reviewed_sequences
            num_compl_seqs,             # num_complete_sequences
            num_compl_rev_seqs,         # num_complete_reviewed_sequences
            num_residues                # num_residues
        ))
    ora_cur.close()
    ora_con.close()

    with pg_con.cursor() as pg_cur:
        pg_cur.execute("DROP TABLE IF EXISTS signature")
        pg_cur.execute(
            """
            CREATE TABLE signature (
                accession VARCHAR(25) NOT NULL 
                    CONSTRAINT signature_pkey PRIMARY KEY,
                database_id INTEGER NOT NULL,
                name VARCHAR(100) NOT NULL,
                description VARCHAR(400),
                type VARCHAR(25) NOT NULL,
                abstract TEXT,
                num_sequences INTEGER NOT NULL,
                num_reviewed_sequences INTEGER NOT NULL,
                num_complete_sequences INTEGER NOT NULL,
                num_complete_reviewed_sequences INTEGER NOT NULL,
                num_residues BIGINT NOT NULL
            )
            """
        )
        execute_values(pg_cur, "INSERT INTO signature VALUES %s", values,
                       page_size=1000)

        pg_cur.execute(
            """
            CREATE INDEX signature_database_idx
            ON signature (database_id)
            """
        )
        pg_con.commit()

    pg_con.close()
    logger.info("complete")


def get_swissprot_descriptions(pg_url: str) -> dict:
    con = psycopg2.connect(**url2dict(pg_url))
    with con.cursor() as cur:
        cur.execute(
            """
            SELECT DISTINCT s2p.signature_acc, pn.text
            FROM interpro.signature2protein s2p
            INNER JOIN interpro.protein_name pn ON s2p.name_id = pn.name_id
            WHERE s2p.is_reviewed            
            """
        )

        signatures = {}
        for signature_acc, text in cur:
            try:
                signatures[signature_acc].add(text)
            except KeyError:
                signatures[signature_acc] = {text}

    con.close()

    return signatures
