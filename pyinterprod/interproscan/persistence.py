import cx_Oracle


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
                "seq_score": float(cols[9]),
                "seq_evalue": float(cols[10])
            })

            if len(values) == 1000:
                cur.executemany(sql, values)
                con.commit()
                values.clear()

    if values:
        cur.executemany(sql, values)
        con.commit()

    cur.close()
    con.close()


def cdd_sites(uri: str, file: str, analysis_id: int, table: str):
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
    cur.setinputsizes(seq_evalue=cx_Oracle.DB_TYPE_BINARY_DOUBLE)

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

            if len(values) == 1000:
                cur.executemany(sql, values)
                con.commit()
                values.clear()

    if values:
        cur.executemany(sql, values)
        con.commit()

    cur.close()
    con.close()


def hmmer3_matches(uri: str, file: str, analysis_id: int, table: str):
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
                "relno_major": cols[1],
                "relno_minor": cols[2],
                "upi": cols[3],
                "method_ac": cols[4],
                "model_ac": cols[5],
                "seq_start": int(cols[6]),
                "seq_end": int(cols[7]),
                "fragments": cols[8],
                "seq_score": float(cols[9]),
                "seq_evalue": float(cols[10]),
                "hmm_bounds": cols[11],
                "hmm_start": int(cols[12]),
                "hmm_end": int(cols[13]),
                "hmm_length": int(cols[14]),
                "env_start": int(cols[15]),
                "env_end": int(cols[16]),
                "score": float(cols[17]),
                "evalue": float(cols[18])
            })

            if len(values) == 1000:
                cur.executemany(sql, values)
                con.commit()
                values.clear()

    if values:
        cur.executemany(sql, values)
        con.commit()

    cur.close()
    con.close()
