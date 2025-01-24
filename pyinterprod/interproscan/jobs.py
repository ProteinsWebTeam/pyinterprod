import os
import subprocess
import sys
import time
from datetime import datetime

import oracledb

from pyinterprod import logger


# String printed by I5 on successful completion
_I5_SUCCESS = "100% done:  InterProScan analyses completed"


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


def add_job(cur: oracledb.Cursor,
            analysis_id: int,
            upi_from: str,
            upi_to: str,
            num_sequences: int):
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
        INSERT INTO IPRSCAN.ANALYSIS_JOBS 
            (ANALYSIS_ID, UPI_FROM, UPI_TO, SEQUENCES)
        VALUES (:1, :2, :3, :4)
        """,
        [analysis_id, upi_from, upi_to, num_sequences]
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
               success: bool | None = None):
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
        params.append("Y" if success else "N")

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


def get_unfinished_lsf_jobs() -> dict[str, int]:
    stdout = subprocess.run(["bjobs", "-o", "jobid name"],
                            capture_output=True,
                            encoding="utf-8").stdout
    jobs = {}
    if "No unfinished job found" in stdout:
        return jobs

    for line in stdout.split("\n")[1:]:
        if line:
            job_id, name = line.split(maxsplit=1)
            jobs[name] = int(job_id)

    return jobs


def get_unfinished_slurm_jobs() -> dict[str, int]:
    stdout = subprocess.run(["squeue", "-h" "-O", "JobId,Name",
                             "-t", "pending,running"],
                            capture_output=True,
                            encoding="utf-8").stdout

    jobs = {}
    for line in stdout.split("\n")[1:]:
        if line:
            job_id, name = line.split(maxsplit=1)
            jobs[name] = int(job_id)

    return jobs


def run_job(i5_dir: str,
            applications: str,
            fasta_file: str,
            matches_output: str,
            sites_output: str | None = None,
            cpu: int | None = None,
            timeout: int | None = None):
    workdir = os.path.dirname(matches_output)

    args = [
        os.path.join(i5_dir, "interproscan.sh"),
        "-i", fasta_file,
        "-appl", applications,
        "-dp",
        "-f", "tsv-pro",
        "-o", matches_output,
        "-T", workdir
    ]

    if cpu is not None:
        args += ["-cpu", str(cpu)]

    if isinstance(timeout, int) and timeout > 0:
        # timeout in hours, but subprocess.run takes in seconds
        _timeout = timeout * 3600
    else:
        _timeout = None

    logger.info(f"Command: {' '.join(args)}")
    ts = time.time()
    process = subprocess.run(args, capture_output=True, timeout=_timeout)
    stdout = process.stdout.decode("utf-8")
    stderr = process.stderr.decode("utf-8")
    code = process.returncode
    runtime = time.time() - ts

    # Write captured streams
    print(stdout, file=sys.stdout)
    print(stderr, file=sys.stderr)

    logger.info(f"Process exited with code {code} after {runtime:.0f} seconds.")

    if code != 0 or _I5_SUCCESS not in stdout:
        raise RuntimeError("InterProScan error")
    elif not os.path.isfile(matches_output):
        raise RuntimeError(f"Matches output file not found")
    elif sites_output is not None and not os.path.isfile(sites_output):
        raise RuntimeError(f"Sites output file not found")
