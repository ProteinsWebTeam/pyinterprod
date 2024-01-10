import os
import random
import shutil
import subprocess
import sys
from dataclasses import dataclass
from typing import Callable

import oracledb
from mundone import Pool, Task
from mundone.states import PENDING, RUNNING

from pyinterprod import logger
from pyinterprod.uniprot.uniparc import int_to_upi, upi_to_int, range_upi
from pyinterprod.utils import oracle
from . import database, persistence

import time


"""
Key   -> analysis name in the database
Value -> tuple of three items
            - analysis name in InterProScan
            - persistence function for matches (.tsv-pro)
            - persistence function for sites (.tsv-pro.sites)
"""
_DB_TO_I5 = {
    "AntiFam": ("AntiFam", persistence.hmmer3_matches, None),
    "CATH-Gene3D": ("Gene3D", persistence.hmmer3_matches, None),
    "CDD": ("CDD", persistence.cdd_matches, persistence.sites),
    "COILS": ("Coils", persistence.coils_phobius_matches, None),
    "FunFam": ("FunFam", persistence.funfam_matches, None),
    "HAMAP": ("Hamap", persistence.hamap_matches, None),
    "MobiDB Lite": ("MobiDBLite", persistence.mobidb_lite_matches, None),
    "NCBIfam": ("NCBIfam", persistence.hmmer3_matches, None),
    "PANTHER": ("PANTHER", persistence.panther_matches, None),
    "Pfam": ("Pfam", persistence.hmmer3_matches, None),
    "Phobius": ("Phobius", persistence.coils_phobius_matches, None),
    "PIRSF": ("PIRSF", persistence.hmmer3_matches, None),
    "PIRSR": ("PIRSR", persistence.pirsr_matches, persistence.sites),
    "PRINTS": ("PRINTS", persistence.prints_matches, None),
    "PROSITE patterns": ("ProSitePatterns",
                         persistence.prosite_patterns_matches, None),
    "PROSITE profiles": ("ProSiteProfiles",
                         persistence.prosite_profiles_matches, None),
    "SFLD": ("SFLD", persistence.sfld_matches, persistence.sites),
    "SignalP_Euk": ("SignalP_EUK", persistence.signalp_matches, None),
    "SignalP_Gram_positive": ("SignalP_GRAM_POSITIVE",
                              persistence.signalp_matches, None),
    "SignalP_Gram_negative": ("SignalP_GRAM_NEGATIVE",
                              persistence.signalp_matches, None),
    "SMART": ("SMART", persistence.smart_matches, None),
    "SUPERFAMILY": ("SUPERFAMILY", persistence.superfamily_matches, None),
    "TIGRFAMs": ("TIGRFAM", persistence.hmmer3_matches, None),
    "TMHMM": ("TMHMM", persistence.tmhmm_matches, None),
}

# String printed by I5 on successful completion
_I5_SUCCESS = "100% done:  InterProScan analyses completed"

_JOB_PREFIX = "IPM_"


def sanitize_name(string: str) -> str:
    for c in [" ", "-", "_"]:
        string = string.replace(c, "")

    return string.lower()


@dataclass
class TaskFactory:
    uri: str
    i5_dir: str
    appl: str
    version: str
    work_dir: str
    analysis_id: int
    config: dict
    match_table: str
    persist_matches: Callable
    site_table: str | None = None
    persist_sites: Callable | None = None
    keep_files: str | None = None
    scheduler: str | None = None
    queue: str | None = None

    def make(self, upi_from: str, upi_to: str) -> Task:
        return Task(
            fn=run_job,
            args=(
                self.uri,
                upi_from,
                upi_to,
                self.i5_dir,
                self.appl,
                os.path.join(self.work_dir, f"{upi_from}_{upi_to}"),
                self.analysis_id,
                self.match_table,
                self.persist_matches,
                self.site_table,
                self.persist_sites
            ),
            kwargs=dict(cpu=self.config["job_cpu"],
                        keep_files=self.keep_files,
                        timeout=self.config["job_timeout"]),
            name=self.make_name(upi_from, upi_to),
            scheduler=dict(type=self.scheduler,
                           queue=self.queue,
                           cpu=self.config["job_cpu"],
                           mem=self.config["job_mem"],
                           hours=self.config["job_timeout"]),
            random_suffix=False
        )

    def make_name(self, upi_from: str, upi_to: str) -> str:
        return f"{_JOB_PREFIX}{self.appl}_{self.version}_{upi_from}_{upi_to}"


def run(uri: str, work_dir: str, temp_dir: str, **kwargs):
    base_config = {
        "job_cpu": kwargs.get("job_cpu", 8),
        "job_mem": kwargs.get("job_mem", 8 * 1024),
        "job_size": kwargs.get("job_size", 32000),
        "job_timeout": kwargs.get("job_timeout"),
    }
    custom_configs = kwargs.get("config", {})
    dry_run = kwargs.get("dry_run", False)
    infinite_mem = kwargs.get("infinite_mem", False)
    keep_files = kwargs.get("keep_files", None)
    max_retries = kwargs.get("max_retries", 0)
    max_timeout = kwargs.get("max_timeout", 0.016)  # 120
    max_running_jobs = kwargs.get("max_running_jobs", 1000)
    max_jobs_per_analysis = kwargs.get("max_jobs_per_analysis", -1)
    pool_threads = kwargs.get("pool_threads", 4)
    scheduler = kwargs.get("scheduler")
    queue = kwargs.get("queue")
    to_run = kwargs.get("analyses", [])
    to_exclude = kwargs.get("exclude", [])

    con = oracledb.connect(uri)
    cur = con.cursor()

    # Find analyses to run
    analyses = {}
    configs = {}
    name2id = {}

    for analysis_id, analysis in database.get_analyses(cur).items():
        name = sanitize_name(analysis["name"])
        try:
            name2id[name].append(analysis_id)
        except KeyError:
            name2id[name] = [analysis_id]

        if analysis_id in to_exclude:
            # Analysis flagged NOT to run
            continue
        elif to_run and analysis_id not in to_run:
            # Analysis not in the list of analyses to run
            continue

        # Either not list of analyses (i.e. run all), or part of this list
        analyses[analysis_id] = analysis

        # Apply default config
        configs[analysis_id] = base_config.copy()

    if not analyses:
        cur.close()
        con.close()
        logger.error("No analyses to process: exit")
        return

    # Override default config with custom configs
    for name in custom_configs:
        key = sanitize_name(name)
        if key not in name2id:
            cur.close()
            con.close()
            raise ValueError(f"Invalid analysis name '{name}'")

        for analysis_id in name2id[key]:
            if analysis_id in configs:
                for key, value in custom_configs[name].items():
                    configs[analysis_id][key] = value

    incomplete_jobs = database.get_incomplete_jobs(cur)
    cur.execute("SELECT MAX(UPI) FROM UNIPARC.PROTEIN")
    max_upi, = cur.fetchone()

    if scheduler.lower() == "lsf":
        name2id = get_unfinished_lsf_jobs()
    elif scheduler.lower() == "slurm":
        name2id = get_unfinished_slurm_jobs()
    else:
        raise ValueError(scheduler)

    with Pool(temp_dir, max_running_jobs,
              kill_on_exit=False,
              threads=pool_threads) as pool:
        running_jobs = []
        pending_jobs = []
        for analysis_id, analysis in analyses.items():
            name = analysis["name"]
            version = analysis["version"]

            appl, persist_matches, persist_sites = _DB_TO_I5[name]

            analysis_work_dir = os.path.join(work_dir, appl, version)
            os.makedirs(analysis_work_dir, exist_ok=True)

            config = configs[analysis_id]
            factory = TaskFactory(uri=uri, i5_dir=analysis["i5_dir"],
                                  appl=appl, version=version,
                                  work_dir=analysis_work_dir,
                                  analysis_id=analysis_id,
                                  config=config,
                                  match_table=analysis["tables"]["matches"],
                                  persist_matches=persist_matches,
                                  site_table=analysis["tables"]["sites"],
                                  persist_sites=persist_sites,
                                  keep_files=keep_files,
                                  scheduler=scheduler,
                                  queue=queue)

            n_tasks_analysis = 0
            jobs = incomplete_jobs.pop(analysis_id, [])
            for upi_from, upi_to, is_running in jobs:
                if dry_run:
                    """
                    We're not monitoring/submitting tasks, so we don't need
                    a real task object
                    """
                    if 0 <= max_jobs_per_analysis <= n_tasks_analysis:
                        break
                    else:
                        running_jobs.append(None)
                        n_tasks_analysis += 1
                        continue

                task = factory.make(upi_from, upi_to)
                task.workdir = os.path.join(temp_dir, task.name)

                if is_running:
                    # Flagged as running in the database

                    # Assumes task is running
                    task.status = RUNNING

                    # Checks if the associated job is running
                    if task.name in name2id:
                        # It is!
                        task.executor.id = name2id.pop(task.name)
                        running_jobs.append(task)
                        n_tasks_analysis += 1
                        continue

                    task.poll()  # Checks if output exists
                    if task.is_done():
                        """
                        Completed or failed. Will be submitted to the pool
                        which will send it back without restarting it.
                        """
                        running_jobs.append(task)
                    elif 0 <= max_jobs_per_analysis <= n_tasks_analysis:
                        break
                    else:
                        task.status = PENDING
                        running_jobs.append(task)
                        n_tasks_analysis += 1
                elif 0 <= max_jobs_per_analysis <= n_tasks_analysis:
                    break
                else:
                    # Flagged as failed in the database

                    # Create new job in the database
                    database.add_job(cur, analysis_id, upi_from, upi_to)

                    # Add task to queue
                    task.status = PENDING
                    running_jobs.append(task)
                    n_tasks_analysis += 1

            if analysis["max_upi"]:
                next_upi = int_to_upi(upi_to_int(analysis["max_upi"]) + 1)
            else:
                next_upi = int_to_upi(1)

            job_size = config["job_size"]
            for upi_from, upi_to in range_upi(next_upi, max_upi, job_size):
                if 0 <= max_jobs_per_analysis <= n_tasks_analysis:
                    break
                elif dry_run:
                    pending_jobs.append(None)
                    n_tasks_analysis += 1
                    continue

                database.add_job(cur, analysis_id, upi_from, upi_to)

                if count_sequences(cur, upi_from, upi_to) > 0:
                    pending_jobs.append(factory.make(upi_from, upi_to))
                    n_tasks_analysis += 1
                else:
                    database.set_job_done(cur, analysis_id, upi_from, upi_to, 0)
                    con.commit()

            logger.debug(f"{name} {version}: {n_tasks_analysis} tasks")

        cur.close()
        con.close()

        n_tasks = len(running_jobs) + len(pending_jobs)
        logger.info(f"Tasks: {n_tasks:,}")

        if dry_run:
            return

        """
        Tasks in the list are grouped by analysis, so if we submit them 
        in order, tasks from the "first" analysis will start to run first, 
        and tasks from other analyses will be pending for a while.
        By randomly selecting tasks to submit, we hope that each analysis 
        will rapidly have some tasks running.
        We could use random.shuffle, but we do not need the entire list 
        re-arranged, so we use random.randint to destructively iterate 
        over the list.
        """
        for obj in [running_jobs, pending_jobs]:
            # Take care of already submitted jobs first, then new ones
            while obj:
                i = random.randrange(len(obj))  # 0 <= i < len(obj)
                task = obj.pop(i)
                pool.submit(task)

        con = oracle.try_connect(uri)
        cur = con.cursor()

        n_completed = n_failed = 0
        milestone = step = 5
        retries = {}
        while (n_completed + n_failed) < n_tasks:
            for task in pool.as_completed(wait=True):
                upi_from = task.args[1]
                upi_to = task.args[2]
                analysis_id = task.args[6]

                # Check if job updated as completed in Oracle
                ok = database.is_job_done(cur, analysis_id, upi_from, upi_to)

                # Update job metadata
                database.update_job(cur, analysis_id, upi_from, upi_to, task)

                logfile = os.path.join(temp_dir, f"{task.name}.log")
                if ok:
                    # Remove the log file (exists if a previous run failed)
                    try:
                        os.unlink(logfile)
                    except FileNotFoundError:
                        pass

                    n_completed += 1
                else:
                    # Write to log file
                    with open(logfile, "wt") as fh:
                        fh.write(task.stdout)
                        fh.write(task.stderr)

                    # Number of times the task was re-submitted
                    num_retries = retries.get(task.name, 0)

                    # Did the job reached the memory usage limit?
                    maxmem = task.maxmem
                    if (maxmem is not None and
                            task.executor.memory is not None and
                            maxmem >= task.executor.memory):
                        mem_err = True
                    else:
                        mem_err = False

                    # Did the job reached the timeout limit?
                    start_time, end_time = task.executor.get_times()
                    runtime = (end_time - start_time).total_seconds() / 3600
                    if runtime >= task.executor.limit:
                        time_err = True
                    else:
                        time_err = False

                    if num_retries < max_retries or time_err or (mem_err and infinite_mem):
                        # Task allowed to be re-submitted

                        # Increase hours if time limit reached
                        if time_err and (task.executor.limit * 1.25 < max_timeout):
                            task.executor.limit *= 1.25

                        try:
                            # Increase memory requirement if needed
                            while task.executor.memory < maxmem:
                                task.executor.memory *= 1.5
                        except TypeError:
                            pass

                        # Resubmit task
                        task.status = PENDING
                        pool.submit(task)

                        # Add new job
                        database.add_job(cur, analysis_id, upi_from, upi_to)

                        # Increment retries counter
                        retries[task.name] = num_retries + 1
                    else:
                        # Max number of retries reached
                        n_failed += 1

                progress = (n_completed + n_failed) * 100 / n_tasks
                if progress >= milestone:
                    """
                    Skip intermediate milestones  (e.g. progress went 
                    from 0 to 10%): we don't want to print 5%
                    """
                    while progress >= milestone:
                        milestone += step

                    logger.info(f"Progress: {milestone - step:>3}%")

        if n_failed:
            logger.error(f"{n_failed} task(s) failed")
        else:
            logger.info("complete")


def count_sequences(cur: oracledb.Cursor, upi_from: str, upi_to: str) -> int:
    cur.execute(
        """
        SELECT COUNT(*)
        FROM UNIPARC.PROTEIN
        WHERE UPI BETWEEN :1 AND :2
        """,
        [upi_from, upi_to]
    )
    cnt, = cur.fetchone()
    return cnt


def export_fasta(uri: str, fasta_file: str, upi_from: str, upi_to: str) -> int:
    num_sequences = 0
    with open(fasta_file, "wt") as fh:
        con = oracle.try_connect(uri)
        cur = con.cursor()

        cur.execute(
            """
            SELECT UPI, SEQ_SHORT, SEQ_LONG
            FROM UNIPARC.PROTEIN
            WHERE UPI BETWEEN :1 AND :2
            """,
            [upi_from, upi_to]
        )

        for upi, seq_short, seq_long in cur:
            sequence = seq_short or seq_long.read()

            fh.write(f">{upi}\n")
            for i in range(0, len(sequence), 60):
                fh.write(f"{sequence[i:i+60]}\n")

            num_sequences += 1

        cur.close()
        con.close()

    return num_sequences


def run_i5(i5_dir: str, fasta_file: str, analysis_name: str, output: str,
           cpu: int | None = None, temp_dir: str | None = None,
           timeout: int | None = None) -> tuple[bool, str, str]:
    args = [
        os.path.join(i5_dir, "interproscan.sh"),
        "-i", fasta_file,
        "-appl", analysis_name,
        "-dp",
        "-f", "tsv-pro",
        "-o", output
    ]

    if cpu is not None:
        args += ["-cpu", str(cpu)]

    if temp_dir is not None:
        args += ["-T", temp_dir]

    if isinstance(timeout, int) and timeout > 0:
        # timeout in hours, but subprocess.run takes in seconds
        _timeout = timeout * 3600
    else:
        _timeout = None

    process = subprocess.run(args, capture_output=True, timeout=_timeout)
    return (
        process.returncode == 0,
        process.stdout.decode("utf-8"),
        process.stderr.decode("utf-8")
    )


def run_job(*args):
    time.sleep(0.005)


def _run_job(uri: str, upi_from: str, upi_to: str, i5_dir: str, appl: str,
            outdir: str, analysis_id: int, match_table: str,
            persist_matches: Callable, site_table: str | None,
            persist_sites: Callable | None, cpu: int | None = None,
            keep_files: str | None = None, timeout: int | None = None):
    try:
        shutil.rmtree(outdir)
    except FileNotFoundError:
        pass

    os.makedirs(outdir)
    fasta_file = os.path.join(outdir, "input.fasta")
    matches_output = os.path.join(outdir, "output.tsv-pro")
    sites_output = matches_output + ".sites"

    ok = False
    try:
        num_sequences = export_fasta(uri, fasta_file, upi_from, upi_to)

        if num_sequences > 0:
            i5_ok, stdout, stderr = run_i5(i5_dir, fasta_file, appl,
                                           matches_output, cpu=cpu,
                                           temp_dir=outdir, timeout=timeout)

            # Write captured streams
            sys.stdout.write(stdout)
            sys.stderr.write(stderr)

            if not i5_ok or _I5_SUCCESS not in stdout:
                raise RuntimeError("InterProScan error")
            elif not os.path.isfile(matches_output):
                raise RuntimeError(f"Cannot access output matches tsv-pro")
            elif site_table and not os.path.isfile(sites_output):
                raise RuntimeError(f"Cannot access output sites tsv-pro")

            con = oracle.try_connect(uri)

            """
            Use a different cursor for persist functions
            as they call cursor.setinputsizes()
            """
            cur = con.cursor()
            persist_matches(cur, matches_output, analysis_id, match_table)
            if site_table:
                persist_sites(cur, sites_output, analysis_id, site_table)
            cur.close()

            cur = con.cursor()
            database.set_job_done(cur, analysis_id, upi_from, upi_to,
                                  num_sequences)

            con.commit()
            cur.close()
            con.close()
        else:
            con = oracle.try_connect(uri)
            cur = con.cursor()
            database.set_job_done(cur, analysis_id, upi_from, upi_to,
                                  num_sequences)
            con.commit()
            cur.close()
            con.close()
    except Exception:
        raise
    else:
        ok = True
    finally:
        if keep_files != "all" and (keep_files != "failed" or ok):
            shutil.rmtree(outdir)


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
