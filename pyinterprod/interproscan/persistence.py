from decimal import Decimal

import cx_Oracle


_COMMIT_SIZE = 10000


def cdd_matches(uri: str, file: str, analysis_id: int, table: str):
    sql = f"""
        INSERT /*+ APPEND */ INTO {table} (
            ANALYSIS_ID, ANALYSIS_NAME, RELNO_MAJOR, RELNO_MINOR,
            UPI, METHOD_AC, MODEL_AC, SEQ_START, SEQ_END, FRAGMENTS,
            SEQSCORE, SEQEVALUE
        )
        VALUES (:analysis_id, :analysis_name, :relno_major, :relno_minor,
                :upi, :method_ac, :model_ac, :seq_start, :seq_end, :fragments,
                :seq_score, :seq_evalue)
    """

    con = cx_Oracle.connect(uri)
    cur = con.cursor()
    cur.setinputsizes(seq_evalue=cx_Oracle.DB_TYPE_BINARY_DOUBLE)

    values = []
    with open(file, "rt") as fh:
        for line in fh:
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

            if len(values) == _COMMIT_SIZE:
                cur.executemany(sql, values)
                con.commit()
                values.clear()

    if values:
        cur.executemany(sql, values)
        con.commit()

    cur.close()
    con.close()


def sites(uri: str, file: str, analysis_id: int, table: str):
    sql = f"""
        INSERT /*+ APPEND */ INTO {table} (
            ANALYSIS_ID, UPI, MD5, SEQ_LENGTH, ANALYSIS_NAME, 
            METHOD_AC, LOC_START, LOC_END, NUM_SITES, RESIDUE, 
            RES_START, RES_END, DESCRIPTION
        )
        VALUES (:analysis_id, :upi, :md5, :seq_length, :analysis_name, 
                :method_ac, :loc_start, :loc_end, :num_sites, :residue,
                :res_start, :res_end, :description)
    """

    con = cx_Oracle.connect(uri)
    cur = con.cursor()
    values = []
    with open(file, "rt") as fh:
        for line in fh:
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

            if len(values) == _COMMIT_SIZE:
                cur.executemany(sql, values)
                con.commit()
                values.clear()

    if values:
        cur.executemany(sql, values)
        con.commit()

    cur.close()
    con.close()


def coils_phobius_matches(uri: str, file: str, analysis_id: int, table: str):
    sql = f"""
        INSERT /*+ APPEND */ INTO {table} (
            ANALYSIS_ID, ANALYSIS_NAME, RELNO_MAJOR, RELNO_MINOR,
            UPI, METHOD_AC, MODEL_AC, SEQ_START, SEQ_END, FRAGMENTS
        )
        VALUES (:analysis_id, :analysis_name, :relno_major, :relno_minor,
                :upi, :method_ac, :model_ac, :seq_start, :seq_end, :fragments)
    """

    con = cx_Oracle.connect(uri)
    cur = con.cursor()

    values = []
    with open(file, "rt") as fh:
        for line in fh:
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

            if len(values) == _COMMIT_SIZE:
                cur.executemany(sql, values)
                con.commit()
                values.clear()

    if values:
        cur.executemany(sql, values)
        con.commit()

    cur.close()
    con.close()


def hamap_matches(uri: str, file: str, analysis_id: int, table: str):
    sql = f"""
        INSERT /*+ APPEND */ INTO {table} (
            ANALYSIS_ID, ANALYSIS_NAME, RELNO_MAJOR, RELNO_MINOR,
            UPI, METHOD_AC, MODEL_AC, SEQ_START, SEQ_END, FRAGMENTS,
            SEQSCORE, ALIGNMENT
        )
        VALUES (:analysis_id, :analysis_name, :relno_major, :relno_minor,
                :upi, :method_ac, :model_ac, :seq_start, :seq_end, :fragments,
                :seq_score, :alignment)
    """

    con = cx_Oracle.connect(uri)
    cur = con.cursor()

    values = []
    with open(file, "rt") as fh:
        for line in fh:
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

            if len(values) == _COMMIT_SIZE:
                cur.executemany(sql, values)
                con.commit()
                values.clear()

    if values:
        cur.executemany(sql, values)
        con.commit()

    cur.close()
    con.close()


def _hmmer3_matches(uri: str, file: str, analysis_id: int, table: str,
                    relno_maj_as_int: bool):
    sql = f"""
        INSERT /*+ APPEND */ INTO {table} (
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

    con = cx_Oracle.connect(uri)
    cur = con.cursor()
    cur.setinputsizes(seq_evalue=cx_Oracle.DB_TYPE_BINARY_DOUBLE,
                      evalue=cx_Oracle.DB_TYPE_BINARY_DOUBLE)

    values = []
    with open(file, "rt") as fh:
        for line in fh:
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

            if len(values) == _COMMIT_SIZE:
                cur.executemany(sql, values)
                con.commit()
                values.clear()

    if values:
        cur.executemany(sql, values)
        con.commit()

    cur.close()
    con.close()


def hmmer3_matches(uri: str, file: str, analysis_id: int, table: str):
    _hmmer3_matches(uri, file, analysis_id, table, relno_maj_as_int=True)


def mobidb_lite_matches(uri: str, file: str, analysis_id: int, table: str):
    sql = f"""
        INSERT /*+ APPEND */ INTO {table} (
            ANALYSIS_ID, ANALYSIS_NAME, RELNO_MAJOR, RELNO_MINOR,
            UPI, METHOD_AC, MODEL_AC, SEQ_START, SEQ_END, FRAGMENTS,
            SEQ_FEATURE
        )
        VALUES (:analysis_id, :analysis_name, :relno_major, :relno_minor,
                :upi, :method_ac, :model_ac, :seq_start, :seq_end, :fragments,
                :seq_feature)
    """

    con = cx_Oracle.connect(uri)
    cur = con.cursor()

    values = []
    with open(file, "rt") as fh:
        for line in fh:
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

            if len(values) == _COMMIT_SIZE:
                cur.executemany(sql, values)
                con.commit()
                values.clear()

    if values:
        cur.executemany(sql, values)
        con.commit()

    cur.close()
    con.close()


def panther_matches(uri: str, file: str, analysis_id: int, table: str):
    sql = f"""
        INSERT /*+ APPEND */ INTO {table} (
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

    con = cx_Oracle.connect(uri)
    cur = con.cursor()
    cur.setinputsizes(seq_evalue=cx_Oracle.DB_TYPE_BINARY_DOUBLE)

    values = []
    with open(file, "rt") as fh:
        for line in fh:
            cols = line.rstrip().split('\t')

            try:
                an_node_id = cols[17].strip()
            except IndexError:
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

            if len(values) == _COMMIT_SIZE:
                cur.executemany(sql, values)
                con.commit()
                values.clear()

    if values:
        cur.executemany(sql, values)
        con.commit()

    cur.close()
    con.close()


def pirsr_matches(uri: str, file: str, analysis_id: int, table: str):
    _hmmer3_matches(uri, file, analysis_id, table, relno_maj_as_int=False)


def prints_matches(uri: str, file: str, analysis_id: int, table: str):
    sql = f"""
        INSERT /*+ APPEND */ INTO {table} (
            ANALYSIS_ID, ANALYSIS_NAME, RELNO_MAJOR, RELNO_MINOR,
            UPI, METHOD_AC, MODEL_AC, SEQ_START, SEQ_END, FRAGMENTS,
            SEQSCORE, SEQEVALUE, MOTIF_NUMBER, PVALUE, GRAPHSCAN
        )
        VALUES (:analysis_id, :analysis_name, :relno_major, :relno_minor,
                :upi, :method_ac, :model_ac, :seq_start, :seq_end, :fragments,
                :seq_score, :seq_evalue, :motif_number, :pvalue, :graphscan)
    """

    con = cx_Oracle.connect(uri)
    cur = con.cursor()
    cur.setinputsizes(seq_evalue=cx_Oracle.DB_TYPE_BINARY_DOUBLE,
                      pvalue=cx_Oracle.DB_TYPE_BINARY_DOUBLE)

    values = []
    with open(file, "rt") as fh:
        for line in fh:
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

            if len(values) == _COMMIT_SIZE:
                cur.executemany(sql, values)
                con.commit()
                values.clear()

    if values:
        cur.executemany(sql, values)
        con.commit()

    cur.close()
    con.close()


prosite_profiles_matches = hamap_matches


def prosite_patterns_matches(uri: str, file: str, analysis_id: int,
                             table: str):
    sql = f"""
        INSERT /*+ APPEND */ INTO {table} (
            ANALYSIS_ID, ANALYSIS_NAME, RELNO_MAJOR, RELNO_MINOR,
            UPI, METHOD_AC, MODEL_AC, SEQ_START, SEQ_END, FRAGMENTS,
            LOCATION_LEVEL, ALIGNMENT
        )
        VALUES (:analysis_id, :analysis_name, :relno_major, :relno_minor,
                :upi, :method_ac, :model_ac, :seq_start, :seq_end, :fragments,
                :loc_level, :alignment)
    """

    con = cx_Oracle.connect(uri)
    cur = con.cursor()

    values = []
    with open(file, "rt") as fh:
        for line in fh:
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

            if len(values) == _COMMIT_SIZE:
                cur.executemany(sql, values)
                con.commit()
                values.clear()

    if values:
        cur.executemany(sql, values)
        con.commit()

    cur.close()
    con.close()


def signalp_tmhmm_matches(uri: str, file: str, analysis_id: int, table: str,
                          relno_maj_as_int: bool):
    sql = f"""
        INSERT /*+ APPEND */ INTO {table} (
            ANALYSIS_ID, ANALYSIS_NAME, RELNO_MAJOR, RELNO_MINOR,
            UPI, METHOD_AC, MODEL_AC, SEQ_START, SEQ_END, FRAGMENTS,
            SEQSCORE
        )
        VALUES (:analysis_id, :analysis_name, :relno_major, :relno_minor,
                :upi, :method_ac, :model_ac, :seq_start, :seq_end, :fragments,
                :seq_score)
    """

    con = cx_Oracle.connect(uri)
    cur = con.cursor()

    values = []
    with open(file, "rt") as fh:
        for line in fh:
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

            if len(values) == _COMMIT_SIZE:
                cur.executemany(sql, values)
                con.commit()
                values.clear()

    if values:
        cur.executemany(sql, values)
        con.commit()

    cur.close()
    con.close()


def signalp_matches(uri: str, file: str, analysis_id: int, table: str):
    signalp_tmhmm_matches(uri, file, analysis_id, table,
                          relno_maj_as_int=False)


def smart_matches(uri: str, file: str, analysis_id: int, table: str):
    sql = f"""
        INSERT /*+ APPEND */ INTO {table} (
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

    con = cx_Oracle.connect(uri)
    cur = con.cursor()
    cur.setinputsizes(seq_evalue=cx_Oracle.DB_TYPE_BINARY_DOUBLE,
                      evalue=cx_Oracle.DB_TYPE_BINARY_DOUBLE)

    values = []
    with open(file, "rt") as fh:
        for line in fh:
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

            if len(values) == _COMMIT_SIZE:
                cur.executemany(sql, values)
                con.commit()
                values.clear()

    if values:
        cur.executemany(sql, values)
        con.commit()

    cur.close()
    con.close()


def superfamily_matches(uri: str, file: str, analysis_id: int, table: str):
    sql = f"""
        INSERT /*+ APPEND */ INTO {table} (
            ANALYSIS_ID, ANALYSIS_NAME, RELNO_MAJOR, RELNO_MINOR,
            UPI, METHOD_AC, MODEL_AC, SEQ_START, SEQ_END, FRAGMENTS,
            SEQEVALUE, HMM_LENGTH
        )
        VALUES (:analysis_id, :analysis_name, :relno_major, :relno_minor,
                :upi, :method_ac, :model_ac, :seq_start, :seq_end, :fragments,
                :seq_evalue, :hmm_length)
    """

    con = cx_Oracle.connect(uri)
    cur = con.cursor()
    cur.setinputsizes(seq_evalue=cx_Oracle.DB_TYPE_BINARY_DOUBLE)

    values = []
    with open(file, "rt") as fh:
        for line in fh:
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

            if len(values) == _COMMIT_SIZE:
                cur.executemany(sql, values)
                con.commit()
                values.clear()

    if values:
        cur.executemany(sql, values)
        con.commit()

    cur.close()
    con.close()


def tmhmm_matches(uri: str, file: str, analysis_id: int, table: str):
    signalp_tmhmm_matches(uri, file, analysis_id, table,
                          relno_maj_as_int=True)
