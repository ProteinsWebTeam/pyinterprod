from datetime import datetime

import oracledb


def get_incomplete_jobs(cur: oracledb.Cursor) -> dict[int, tuple]:
    cur.execute(
        """
        SELECT ANALYSIS_ID, UPI_FROM, UPI_TO, END_TIME, SEQUENCES
        FROM IPRSCAN.ANALYSIS_JOBS
        WHERE END_TIME IS NULL AND SUCCESS = 'N'
        UNION 
        SELECT ANALYSIS_ID, UPI_FROM, UPI_TO, END_TIME, SEQUENCES
        FROM (
            SELECT ANALYSIS_ID, UPI_FROM, UPI_TO, END_TIME, SEQUENCES, SUCCESS, 
                   ROW_NUMBER() OVER (
                       PARTITION BY ANALYSIS_ID, UPI_FROM, UPI_TO 
                       ORDER BY END_TIME DESC
                   ) RN
            FROM ANALYSIS_JOBS
        )
        WHERE SUCCESS = 'N' AND RN = 1
        """
    )

    incomplete_jobs = {}
    for analysis_id, upi_from, upi_to, end_time, sequences in cur:
        is_running = end_time is None
        try:
            analysis_jobs = incomplete_jobs[analysis_id]
        except KeyError:
            analysis_jobs = incomplete_jobs[analysis_id] = []
        finally:
            analysis_jobs.append((upi_from, upi_to, is_running, sequences))

    return incomplete_jobs


def add_job(cur: oracledb.Cursor, analysis_id: int, upi_from: str, upi_to: str):
    cur.execute(
        """
        UPDATE IPRSCAN.ANALYSIS
        SET MAX_UPI = :upi
        WHERE ID = :analysis_id AND (MAX_UPI IS NULL OR MAX_UPI < :upi)
        """,
        analysis_id=analysis_id, upi=upi_to
    )
    cur.execute(
        """
        INSERT INTO IPRSCAN.ANALYSIS_JOBS (ANALYSIS_ID, UPI_FROM, UPI_TO)
        VALUES (:1, :2, :3)
        """,
        [analysis_id, upi_from, upi_to]
    )
    cur.connection.commit()


def update_job(cur: oracledb.Cursor,
               analysis_id: int,
               upi_from: str,
               upi_to: str,
               submit_time: datetime | None = None,
               start_time: datetime | None = None,
               end_time: datetime | None = None,
               max_mem: int | None = None,
               lim_mem: int | None = None,
               cpu_time: int | None = None,
               success: bool | None = None,
               sequences: int | None = None):
    columns = []
    params = []
    if submit_time is not None:
        columns.append("SUBMIT_TIME = :sbmtm")
        params.append(submit_time)
    if start_time is not None:
        columns.append("START_TIME = :strtm")
        params.append(start_time)
    if end_time is not None:
        columns.append("END_TIME = :endtm")
        params.append(end_time)
    if max_mem is not None:
        columns.append("MAX_MEMORY = :maxmem")
        params.append(max_mem)
    if lim_mem is not None:
        columns.append("LIM_MEMORY = :limmem")
        params.append(lim_mem)
    if cpu_time is not None:
        columns.append("CPU_TIME = :cputm")
        params.append(cpu_time)
    if success is not None:
        columns.append("SUCCESS = :success")
        params.append("Y" if success else None)
    if sequences is not None:
        columns.append("SEQUENCES = :sequences")
        params.append(sequences)

    if columns:
        cur.execute(
            f"""
            UPDATE IPRSCAN.ANALYSIS_JOBS
            SET {','.join(columns)}
            WHERE ANALYSIS_ID = :analysisid
                AND UPI_FROM = :upifrom
                AND UPI_TO = :upito
                AND END_TIME IS NULL
            """,
            params + [analysis_id, upi_from, upi_to]
        )
        cur.connection.commit()
