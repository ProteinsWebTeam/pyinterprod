from decimal import Decimal
from typing import Callable

import oracledb
from oracledb.exceptions import IntegrityError

from pyinterprod import logger


_INSERT_SIZE = 10000


def get_analyses(obj: str | oracledb.Cursor) -> dict:
    if isinstance(obj, str):
        con = oracledb.connect(obj)
        cur = con.cursor()
    else:
        con = None
        cur = obj

    cur.execute(
        """
        SELECT A.ID, A.NAME, A.VERSION, A.MAX_UPI, A.I5_DIR,
               T.MATCH_TABLE, T.SITE_TABLE
        FROM IPRSCAN.ANALYSIS A
        INNER JOIN IPRSCAN.ANALYSIS_TABLES T
            ON LOWER(A.NAME) = LOWER(T.NAME)
        WHERE A.ACTIVE = 'Y'
            AND I5_DIR IS NOT NULL
        """
    )

    analyses = {}
    for row in cur:
        analyses[row[0]] = {
            "name": row[1],
            "version": row[2],
            "max_upi": row[3],
            "i5_dir": row[4],
            "tables": {
                "matches": row[5],
                "sites": row[6],
            }
        }

    if con is not None:
        cur.close()
        con.close()

    return analyses


def persist_results(cur: oracledb.Cursor,
                    analysis_id: int,
                    matches_fn: Callable,
                    matches_file: str,
                    matches_table: str,
                    sites_fn: Callable | None,
                    sites_file: str | None,
                    sites_table: str | None) -> bool:
    try:
        matches_fn(cur, matches_file, analysis_id, matches_table)

        if sites_fn is not None:
            sites_fn(cur, sites_file, analysis_id, sites_table)
    except IntegrityError:
        cur.connection.rollback()
        return False
    else:
        cur.connection.commit()
        return True


def cdd_matches(cur: oracledb.Cursor, file: str, analysis_id: int, table: str):
    sql = f"""
        INSERT INTO {table} (
            ANALYSIS_ID, ANALYSIS_NAME, RELNO_MAJOR, RELNO_MINOR,
            UPI, METHOD_AC, MODEL_AC, SEQ_START, SEQ_END, FRAGMENTS,
            SEQSCORE, SEQEVALUE
        )
        VALUES (:analysis_id, :analysis_name, :relno_major, :relno_minor,
                :upi, :method_ac, :model_ac, :seq_start, :seq_end, :fragments,
                :seq_score, :seq_evalue)
    """

    cur.setinputsizes(seq_evalue=oracledb.DB_TYPE_BINARY_DOUBLE)

    values = []
    num_parsed = num_inserted = 0
    with open(file, "rt") as fh:
        for line in fh:
            num_parsed += 1
            cols = line.rstrip().split('\t')
            values.append({
                "analysis_id": analysis_id,
                "analysis_name": cols[0],
                "relno_major": cols[1],
                "relno_minor": cols[2],
                "upi": cols[3],
                "method_ac": cols[4],
                "model_ac": cols[5],
                "seq_start": int(cols[6]),
                "seq_end": int(cols[7]),
                "fragments": cols[8],
                "seq_score": Decimal(cols[9]),
                "seq_evalue": Decimal(cols[10])
            })

            if len(values) == _INSERT_SIZE:
                cur.executemany(sql, values)
                num_inserted += cur.rowcount
                values.clear()

    if values:
        cur.executemany(sql, values)
        num_inserted += cur.rowcount

    logger.debug(f"parsed: {num_parsed}; inserted: {num_inserted}")


def sites(cur: oracledb.Cursor, file: str, analysis_id: int, table: str):
    sql = f"""
        INSERT INTO {table} (
            ANALYSIS_ID, UPI, MD5, SEQ_LENGTH, ANALYSIS_NAME, 
            METHOD_AC, LOC_START, LOC_END, NUM_SITES, RESIDUE, 
            RES_START, RES_END, DESCRIPTION
        )
        VALUES (:analysis_id, :upi, :md5, :seq_length, :analysis_name, 
                :method_ac, :loc_start, :loc_end, :num_sites, :residue,
                :res_start, :res_end, :description)
    """

    values = []
    num_parsed = num_inserted = 0
    with open(file, "rt") as fh:
        for line in fh:
            num_parsed += 1
            cols = line.rstrip().split('\t')
            values.append({
                "analysis_id": analysis_id,
                "upi": cols[0],
                "md5": cols[1],
                "seq_length": int(cols[2]),
                "analysis_name": cols[3],
                "method_ac": cols[4],
                "loc_start": int(cols[5]),
                "loc_end": int(cols[6]),
                "num_sites": int(cols[7]),
                "residue": cols[8],
                "res_start": int(cols[9]),
                "res_end": int(cols[10]),
                "description": cols[11]
            })

            if len(values) == _INSERT_SIZE:
                cur.executemany(sql, values)
                num_inserted += cur.rowcount
                values.clear()

    if values:
        cur.executemany(sql, values)
        num_inserted += cur.rowcount

    logger.debug(f"parsed: {num_parsed}; inserted: {num_inserted}")


def coils_phobius_matches(cur: oracledb.Cursor, file: str, analysis_id: int,
                          table: str):
    sql = f"""
        INSERT INTO {table} (
            ANALYSIS_ID, ANALYSIS_NAME, RELNO_MAJOR, RELNO_MINOR,
            UPI, METHOD_AC, MODEL_AC, SEQ_START, SEQ_END, FRAGMENTS
        )
        VALUES (:analysis_id, :analysis_name, :relno_major, :relno_minor,
                :upi, :method_ac, :model_ac, :seq_start, :seq_end, :fragments)
    """

    values = []
    num_parsed = num_inserted = 0
    with open(file, "rt") as fh:
        for line in fh:
            num_parsed += 1
            cols = line.rstrip().split('\t')
            values.append({
                "analysis_id": analysis_id,
                "analysis_name": cols[0],
                "relno_major": cols[1],
                "relno_minor": cols[2],
                "upi": cols[3],
                "method_ac": cols[4],
                "model_ac": cols[5],
                "seq_start": int(cols[6]),
                "seq_end": int(cols[7]),
                "fragments": cols[8]
            })

            if len(values) == _INSERT_SIZE:
                cur.executemany(sql, values)
                num_inserted += cur.rowcount
                values.clear()

    if values:
        cur.executemany(sql, values)
        num_inserted += cur.rowcount

    logger.debug(f"parsed: {num_parsed}; inserted: {num_inserted}")


def hamap_matches(cur: oracledb.Cursor, file: str, analysis_id: int,
                  table: str):
    sql = f"""
        INSERT INTO {table} (
            ANALYSIS_ID, ANALYSIS_NAME, RELNO_MAJOR, RELNO_MINOR,
            UPI, METHOD_AC, MODEL_AC, SEQ_START, SEQ_END, FRAGMENTS,
            SEQSCORE, ALIGNMENT
        )
        VALUES (:analysis_id, :analysis_name, :relno_major, :relno_minor,
                :upi, :method_ac, :model_ac, :seq_start, :seq_end, :fragments,
                :seq_score, :alignment)
    """

    values = []
    num_parsed = num_inserted = 0
    with open(file, "rt") as fh:
        for line in fh:
            num_parsed += 1
            cols = line.rstrip().split('\t')
            values.append({
                "analysis_id": analysis_id,
                "analysis_name": cols[0],
                "relno_major": cols[1],
                "relno_minor": cols[2],
                "upi": cols[3],
                "method_ac": cols[4],
                "model_ac": cols[5],
                "seq_start": int(cols[6]),
                "seq_end": int(cols[7]),
                "fragments": cols[8],
                "seq_score": Decimal(cols[9]),
                "alignment": cols[10]
            })

            if len(values) == _INSERT_SIZE:
                cur.executemany(sql, values)
                num_inserted += cur.rowcount
                values.clear()

    if values:
        cur.executemany(sql, values)
        num_inserted += cur.rowcount

    logger.debug(f"parsed: {num_parsed}; inserted: {num_inserted}")


def _hmmer3_matches(cur: oracledb.Cursor, file: str, analysis_id: int,
                    table: str, relno_maj_as_int: bool):
    sql = f"""
        INSERT INTO {table} (
            ANALYSIS_ID, ANALYSIS_NAME, RELNO_MAJOR, RELNO_MINOR,
            UPI, METHOD_AC, MODEL_AC, SEQ_START, SEQ_END, FRAGMENTS,
            SEQSCORE, SEQEVALUE, HMM_BOUNDS, HMM_START, HMM_END,
            HMM_LENGTH, ENV_START, ENV_END, SCORE, EVALUE
        )
        VALUES (:analysis_id, :analysis_name, :relno_major, :relno_minor,
                :upi, :method_ac, :model_ac, :seq_start, :seq_end, :fragments,
                :seq_score, :seq_evalue, :hmm_bounds, :hmm_start, :hmm_end,
                :hmm_length, :env_start, :env_end, :score, :evalue)
    """

    cur.setinputsizes(seq_evalue=oracledb.DB_TYPE_BINARY_DOUBLE,
                      evalue=oracledb.DB_TYPE_BINARY_DOUBLE)

    values = []
    num_parsed = num_inserted = 0
    with open(file, "rt") as fh:
        for line in fh:
            num_parsed += 1
            cols = line.rstrip().split('\t')
            values.append({
                "analysis_id": analysis_id,
                "analysis_name": cols[0],
                "relno_major": int(cols[1]) if relno_maj_as_int else cols[1],
                "relno_minor": cols[2],
                "upi": cols[3],
                "method_ac": cols[4],
                "model_ac": cols[5],
                "seq_start": int(cols[6]),
                "seq_end": int(cols[7]),
                "fragments": cols[8],
                "seq_score": Decimal(cols[9]),
                "seq_evalue": Decimal(cols[10]),
                "hmm_bounds": cols[11],
                "hmm_start": int(cols[12]),
                "hmm_end": int(cols[13]),
                "hmm_length": int(cols[14]),
                "env_start": int(cols[15]),
                "env_end": int(cols[16]),
                "score": Decimal(cols[17]),
                "evalue": Decimal(cols[18])
            })

            if len(values) == _INSERT_SIZE:
                cur.executemany(sql, values)
                num_inserted += cur.rowcount
                values.clear()

    if values:
        cur.executemany(sql, values)
        num_inserted += cur.rowcount

    logger.debug(f"parsed: {num_parsed}; inserted: {num_inserted}")


def funfam_matches(cur: oracledb.Cursor, file: str, analysis_id: int,
                   table: str):
    sql = f"""
        INSERT INTO {table} (
            ANALYSIS_ID, ANALYSIS_NAME, RELNO_MAJOR, RELNO_MINOR,
            UPI, METHOD_AC, MODEL_AC, SEQ_START, SEQ_END, FRAGMENTS,
            SEQSCORE, SEQEVALUE, HMM_BOUNDS, HMM_START, HMM_END,
            HMM_LENGTH, ENV_START, ENV_END, SCORE, EVALUE, HMMER_SEQ_START,
            HMMER_SEQ_END, ALIGNMENT
        )
        VALUES (:analysis_id, :analysis_name, :relno_major, :relno_minor,
                :upi, :method_ac, :model_ac, :seq_start, :seq_end, :fragments,
                :seq_score, :seq_evalue, :hmm_bounds, :hmm_start, :hmm_end,
                :hmm_length, :env_start, :env_end, :score, :evalue, 
                :hmmer_seq_start, :hmmer_seq_end, :alignment)
        """

    cur.setinputsizes(seq_evalue=oracledb.DB_TYPE_BINARY_DOUBLE,
                      evalue=oracledb.DB_TYPE_BINARY_DOUBLE)

    values = []
    num_parsed = num_inserted = 0
    with open(file, "rt") as fh:
        for line in fh:
            num_parsed += 1
            cols = line.rstrip().split('\t')
            values.append({
                "analysis_id": analysis_id,
                "analysis_name": cols[0],
                "relno_major": int(cols[1]),
                "relno_minor": cols[2],
                "upi": cols[3],
                "method_ac": cols[4],
                "model_ac": cols[5],
                "seq_start": int(cols[6]),
                "seq_end": int(cols[7]),
                "fragments": cols[8],
                "seq_score": Decimal(cols[9]),
                "seq_evalue": Decimal(cols[10]),
                "hmm_bounds": cols[11],
                "hmm_start": int(cols[12]),
                "hmm_end": int(cols[13]),
                "hmm_length": int(cols[14]),
                "env_start": int(cols[15]),
                "env_end": int(cols[16]),
                "score": Decimal(cols[17]),
                "evalue": Decimal(cols[18]),
                "hmmer_seq_start": int(cols[19]),
                "hmmer_seq_end": int(cols[20]),
                "alignment": cols[21],
            })

            if len(values) == _INSERT_SIZE:
                cur.executemany(sql, values)
                num_inserted += cur.rowcount
                values.clear()

    if values:
        cur.executemany(sql, values)
        num_inserted += cur.rowcount

    logger.debug(f"parsed: {num_parsed}; inserted: {num_inserted}")


def hmmer3_matches(cur: oracledb.Cursor, file: str, analysis_id: int,
                   table: str):
    _hmmer3_matches(cur, file, analysis_id, table, relno_maj_as_int=True)


def mobidb_lite_matches(cur: oracledb.Cursor, file: str, analysis_id: int,
                        table: str):
    sql = f"""
        INSERT INTO {table} (
            ANALYSIS_ID, ANALYSIS_NAME, RELNO_MAJOR, RELNO_MINOR,
            UPI, METHOD_AC, MODEL_AC, SEQ_START, SEQ_END, FRAGMENTS,
            SEQ_FEATURE
        )
        VALUES (:analysis_id, :analysis_name, :relno_major, :relno_minor,
                :upi, :method_ac, :model_ac, :seq_start, :seq_end, :fragments,
                :seq_feature)
    """

    values = []
    num_parsed = num_inserted = 0
    with open(file, "rt") as fh:
        for line in fh:
            num_parsed += 1
            cols = line.rstrip().split('\t')

            try:
                seq_feature = cols[9].strip()
            except IndexError:
                seq_feature = None

            values.append({
                "analysis_id": analysis_id,
                "analysis_name": cols[0],
                "relno_major": int(cols[1]),
                "relno_minor": cols[2],
                "upi": cols[3],
                "method_ac": cols[4],
                "model_ac": cols[5],
                "seq_start": int(cols[6]),
                "seq_end": int(cols[7]),
                "fragments": cols[8],
                "seq_feature": seq_feature
            })

            if len(values) == _INSERT_SIZE:
                cur.executemany(sql, values)
                num_inserted += cur.rowcount
                values.clear()

    if values:
        cur.executemany(sql, values)
        num_inserted += cur.rowcount

    logger.debug(f"parsed: {num_parsed}; inserted: {num_inserted}")


def panther_matches(cur: oracledb.Cursor, file: str, analysis_id: int,
                    table: str):
    sql = f"""
        INSERT INTO {table} (
            ANALYSIS_ID, ANALYSIS_NAME, RELNO_MAJOR, RELNO_MINOR,
            UPI, METHOD_AC, MODEL_AC, SEQ_START, SEQ_END, FRAGMENTS,
            SEQSCORE, SEQEVALUE, HMM_BOUNDS, HMM_START, HMM_END,
            HMM_LENGTH, ENV_START, ENV_END, AN_NODE_ID
        )
        VALUES (:analysis_id, :analysis_name, :relno_major, :relno_minor,
                :upi, :method_ac, :model_ac, :seq_start, :seq_end, :fragments,
                :seq_score, :seq_evalue, :hmm_bounds, :hmm_start, :hmm_end,
                :hmm_length, :env_start, :env_end, :an_node_id)
    """

    cur.setinputsizes(seq_evalue=oracledb.DB_TYPE_BINARY_DOUBLE)

    values = []
    num_parsed = num_inserted = 0
    with open(file, "rt") as fh:
        for line in fh:
            num_parsed += 1
            cols = line.rstrip().split('\t')

            try:
                an_node_id = cols[17].strip()
            except IndexError:
                an_node_id = None
            else:
                if an_node_id == "-":
                    an_node_id = None

            values.append({
                "analysis_id": analysis_id,
                "analysis_name": cols[0],
                "relno_major": cols[1],
                "relno_minor": cols[2],
                "upi": cols[3],
                "method_ac": cols[4],
                "model_ac": cols[5],
                "seq_start": int(cols[6]),
                "seq_end": int(cols[7]),
                "fragments": cols[8],
                "seq_score": Decimal(cols[9]),
                "seq_evalue": Decimal(cols[10]),
                "hmm_bounds": cols[11],
                "hmm_start": int(cols[12]),
                "hmm_end": int(cols[13]),
                "hmm_length": int(cols[14]),
                "env_start": int(cols[15]),
                "env_end": int(cols[16]),
                "an_node_id": an_node_id
            })

            if len(values) == _INSERT_SIZE:
                cur.executemany(sql, values)
                num_inserted += cur.rowcount
                values.clear()

    if values:
        cur.executemany(sql, values)
        num_inserted += cur.rowcount

    logger.debug(f"parsed: {num_parsed}; inserted: {num_inserted}")


def pirsr_matches(cur: oracledb.Cursor, file: str, analysis_id: int,
                  table: str):
    _hmmer3_matches(cur, file, analysis_id, table, relno_maj_as_int=False)


def prints_matches(cur: oracledb.Cursor, file: str, analysis_id: int, table: str):
    sql = f"""
        INSERT INTO {table} (
            ANALYSIS_ID, ANALYSIS_NAME, RELNO_MAJOR, RELNO_MINOR,
            UPI, METHOD_AC, MODEL_AC, SEQ_START, SEQ_END, FRAGMENTS,
            SEQSCORE, SEQEVALUE, MOTIF_NUMBER, PVALUE, GRAPHSCAN
        )
        VALUES (:analysis_id, :analysis_name, :relno_major, :relno_minor,
                :upi, :method_ac, :model_ac, :seq_start, :seq_end, :fragments,
                :seq_score, :seq_evalue, :motif_number, :pvalue, :graphscan)
    """

    cur.setinputsizes(seq_evalue=oracledb.DB_TYPE_BINARY_DOUBLE,
                      pvalue=oracledb.DB_TYPE_BINARY_DOUBLE)

    values = []
    num_parsed = num_inserted = 0
    with open(file, "rt") as fh:
        for line in fh:
            num_parsed += 1
            cols = line.rstrip().split('\t')
            values.append({
                "analysis_id": analysis_id,
                "analysis_name": cols[0],
                "relno_major": cols[1],
                "relno_minor": cols[2],
                "upi": cols[3],
                "method_ac": cols[4],
                "model_ac": cols[5],
                "seq_start": int(cols[6]),
                "seq_end": int(cols[7]),
                "fragments": cols[8],
                "seq_score": Decimal(cols[9]),
                "seq_evalue": Decimal(cols[10]),
                "motif_number": int(cols[11]),
                "pvalue": Decimal(cols[12]),
                "graphscan": cols[13]
            })

            if len(values) == _INSERT_SIZE:
                cur.executemany(sql, values)
                num_inserted += cur.rowcount
                values.clear()

    if values:
        cur.executemany(sql, values)
        num_inserted += cur.rowcount

    logger.debug(f"parsed: {num_parsed}; inserted: {num_inserted}")


prosite_profiles_matches = hamap_matches


def prosite_patterns_matches(cur: oracledb.Cursor, file: str, analysis_id: int,
                             table: str):
    sql = f"""
        INSERT INTO {table} (
            ANALYSIS_ID, ANALYSIS_NAME, RELNO_MAJOR, RELNO_MINOR,
            UPI, METHOD_AC, MODEL_AC, SEQ_START, SEQ_END, FRAGMENTS,
            LOCATION_LEVEL, ALIGNMENT
        )
        VALUES (:analysis_id, :analysis_name, :relno_major, :relno_minor,
                :upi, :method_ac, :model_ac, :seq_start, :seq_end, :fragments,
                :loc_level, :alignment)
    """

    values = []
    num_parsed = num_inserted = 0
    with open(file, "rt") as fh:
        for line in fh:
            num_parsed += 1
            cols = line.rstrip().split('\t')
            values.append({
                "analysis_id": analysis_id,
                "analysis_name": cols[0],
                "relno_major": cols[1],
                "relno_minor": cols[2],
                "upi": cols[3],
                "method_ac": cols[4],
                "model_ac": cols[5],
                "seq_start": int(cols[6]),
                "seq_end": int(cols[7]),
                "fragments": cols[8],
                "loc_level": int(cols[9]),
                "alignment": cols[10]
            })

            if len(values) == _INSERT_SIZE:
                cur.executemany(sql, values)
                num_inserted += cur.rowcount
                values.clear()

    if values:
        cur.executemany(sql, values)
        num_inserted += cur.rowcount

    logger.debug(f"parsed: {num_parsed}; inserted: {num_inserted}")


def signalp_tmhmm_matches(cur: oracledb.Cursor, file: str, analysis_id: int,
                          table: str, relno_maj_as_int: bool):
    sql = f"""
        INSERT INTO {table} (
            ANALYSIS_ID, ANALYSIS_NAME, RELNO_MAJOR, RELNO_MINOR,
            UPI, METHOD_AC, MODEL_AC, SEQ_START, SEQ_END, FRAGMENTS,
            SEQSCORE
        )
        VALUES (:analysis_id, :analysis_name, :relno_major, :relno_minor,
                :upi, :method_ac, :model_ac, :seq_start, :seq_end, :fragments,
                :seq_score)
    """

    values = []
    num_parsed = num_inserted = 0
    with open(file, "rt") as fh:
        for line in fh:
            num_parsed += 1
            cols = line.rstrip().split('\t')
            values.append({
                "analysis_id": analysis_id,
                "analysis_name": cols[0],
                "relno_major": int(cols[1]) if relno_maj_as_int else cols[1],
                "relno_minor": cols[2],
                "upi": cols[3],
                "method_ac": cols[4],
                "model_ac": cols[5],
                "seq_start": int(cols[6]),
                "seq_end": int(cols[7]),
                "fragments": cols[8],
                "seq_score": Decimal(cols[9])
            })

            if len(values) == _INSERT_SIZE:
                cur.executemany(sql, values)
                num_inserted += cur.rowcount
                values.clear()

    if values:
        cur.executemany(sql, values)
        num_inserted += cur.rowcount

    logger.debug(f"parsed: {num_parsed}; inserted: {num_inserted}")


def signalp_matches(cur: oracledb.Cursor, file: str, analysis_id: int,
                    table: str):
    signalp_tmhmm_matches(cur, file, analysis_id, table,
                          relno_maj_as_int=False)


def sfld_matches(cur: oracledb.Cursor, file: str, analysis_id: int, table: str):
    _hmmer3_matches(cur, file, analysis_id, table, relno_maj_as_int=True)


def smart_matches(cur: oracledb.Cursor, file: str, analysis_id: int,
                  table: str):
    sql = f"""
        INSERT INTO {table} (
            ANALYSIS_ID, ANALYSIS_NAME, RELNO_MAJOR, RELNO_MINOR,
            UPI, METHOD_AC, MODEL_AC, SEQ_START, SEQ_END, FRAGMENTS,
            SEQSCORE, SEQEVALUE, HMM_BOUNDS, HMM_START, HMM_END,
            HMM_LENGTH, SCORE, EVALUE
        )
        VALUES (:analysis_id, :analysis_name, :relno_major, :relno_minor,
                :upi, :method_ac, :model_ac, :seq_start, :seq_end, :fragments,
                :seq_score, :seq_evalue, :hmm_bounds, :hmm_start, :hmm_end,
                :hmm_length, :score, :evalue)
    """

    cur.setinputsizes(seq_evalue=oracledb.DB_TYPE_BINARY_DOUBLE,
                      evalue=oracledb.DB_TYPE_BINARY_DOUBLE)

    values = []
    num_parsed = num_inserted = 0
    with open(file, "rt") as fh:
        for line in fh:
            num_parsed += 1
            cols = line.rstrip().split('\t')
            values.append({
                "analysis_id": analysis_id,
                "analysis_name": cols[0],
                "relno_major": int(cols[1]),
                "relno_minor": cols[2],
                "upi": cols[3],
                "method_ac": cols[4],
                "model_ac": cols[5],
                "seq_start": int(cols[6]),
                "seq_end": int(cols[7]),
                "fragments": cols[8],
                "seq_score": Decimal(cols[9]),
                "seq_evalue": Decimal(cols[10]),
                "hmm_bounds": cols[11],
                "hmm_start": int(cols[12]),
                "hmm_end": int(cols[13]),
                "hmm_length": int(cols[14]),
                "score": Decimal(cols[15]),
                "evalue": Decimal(cols[16])
            })

            if len(values) == _INSERT_SIZE:
                cur.executemany(sql, values)
                num_inserted += cur.rowcount
                values.clear()

    if values:
        cur.executemany(sql, values)
        num_inserted += cur.rowcount

    logger.debug(f"parsed: {num_parsed}; inserted: {num_inserted}")


def superfamily_matches(cur: oracledb.Cursor, file: str, analysis_id: int,
                        table: str):
    sql = f"""
        INSERT INTO {table} (
            ANALYSIS_ID, ANALYSIS_NAME, RELNO_MAJOR, RELNO_MINOR,
            UPI, METHOD_AC, MODEL_AC, SEQ_START, SEQ_END, FRAGMENTS,
            SEQEVALUE, HMM_LENGTH
        )
        VALUES (:analysis_id, :analysis_name, :relno_major, :relno_minor,
                :upi, :method_ac, :model_ac, :seq_start, :seq_end, :fragments,
                :seq_evalue, :hmm_length)
    """

    cur.setinputsizes(seq_evalue=oracledb.DB_TYPE_BINARY_DOUBLE)

    values = []
    num_parsed = num_inserted = 0
    with open(file, "rt") as fh:
        for line in fh:
            num_parsed += 1
            cols = line.rstrip().split('\t')
            values.append({
                "analysis_id": analysis_id,
                "analysis_name": cols[0],
                "relno_major": int(cols[1]),
                "relno_minor": cols[2],
                "upi": cols[3],
                "method_ac": cols[4],
                "model_ac": cols[5],
                "seq_start": int(cols[6]),
                "seq_end": int(cols[7]),
                "fragments": cols[8],
                "seq_evalue": Decimal(cols[9]),
                "hmm_length": int(cols[10]),
            })

            if len(values) == _INSERT_SIZE:
                cur.executemany(sql, values)
                num_inserted += cur.rowcount
                values.clear()

    if values:
        cur.executemany(sql, values)
        num_inserted += cur.rowcount

    logger.debug(f"parsed: {num_parsed}; inserted: {num_inserted}")


def tmhmm_matches(cur: oracledb.Cursor, file: str, analysis_id: int,
                  table: str):
    signalp_tmhmm_matches(cur, file, analysis_id, table,
                          relno_maj_as_int=True)
