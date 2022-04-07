import os
import random
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from typing import Callable, Optional

import cx_Oracle
from mundone import Pool, Task
from mundone.task import STATUS_PENDING, STATUS_RUNNING

from pyinterprod import logger
from . import database, persistence


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
    "FunFam": ("FunFam", persistence.hmmer3_matches, None),
    "HAMAP": ("Hamap", persistence.hamap_matches, None),
    "MobiDB Lite": ("MobiDBLite", persistence.mobidb_lite_matches, None),
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
    "SFLD": ("SFLD", persistence.hmmer3_matches, persistence.sites),
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
    parse_matches: Callable
    site_table: Optional[str] = None
    parse_sites: Optional[Callable] = None
    keep_files: Optional[str] = None
    lsf_queue: Optional[str] = None
    timeout: Optional[int] = None

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
                self.parse_matches,
                self.site_table,
                self.parse_sites
            ),
            kwargs=dict(cpu=self.config["job_cpu"],
                        keep_files=self.keep_files,
                        timeout=self.timeout),
            name=self.make_name(upi_from, upi_to),
            scheduler=dict(cpu=self.config["job_cpu"],
                           mem=self.config["job_mem"],
                           queue=self.lsf_queue),
            random_suffix=False
        )

    def make_name(self, upi_from: str, upi_to: str) -> str:
        return f"{_JOB_PREFIX}{self.appl}_{self.version}_{upi_from}_{upi_to}"


def int_to_upi(i):
    return f"UPI{i:010x}".upper()


def upi_to_int(upi):
    return int(upi[3:], 16)


def range_jobs(from_upi: str, to_upi: str, step: int):
    start = upi_to_int(from_upi)
    stop = upi_to_int(to_upi) + 1
    for i in range(start, stop, step):
        yield int_to_upi(i), int_to_upi(i + step)


def run(uri: str, work_dir: str, temp_dir: str, **kwargs):
    base_config = {
        "job_cpu": kwargs.get("job_cpu", 8),
        "job_mem": kwargs.get("job_mem", 8 * 1024),
        "job_size": kwargs.get("job_size", 100000)
    }
    custom_configs = kwargs.get("config", {})
    infinite_mem = kwargs.get("infinite_mem", False)
    keep_files = kwargs.get("keep_files", None)
    lsf_queue = kwargs.get("lsf_queue")
    max_retries = kwargs.get("max_retries", 0)
    max_running_jobs = kwargs.get("max_running_jobs", 1000)
    max_jobs_per_analysis = kwargs.get("max_jobs_per_analysis", -1)
    pool_threads = kwargs.get("pool_threads", 4)
    timeout = kwargs.get("timeout")
    to_run = kwargs.get("analyses", [])
    to_exclude = kwargs.get("exclude", [])

    con = cx_Oracle.connect(uri)
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

    name2id = get_lsf_unfinished_jobs()

    with Pool(temp_dir, max_running_jobs, pool_threads) as pool:
        to_run = []
        n_tasks = 0
        for analysis_id, analysis in analyses.items():
            name = analysis["name"]
            version = analysis["version"]

            if analysis["max_upi"]:
                next_upi = int_to_upi(upi_to_int(analysis["max_upi"]) + 1)
            else:
                next_upi = int_to_upi(1)

            appl, parse_matches, parse_sites = _DB_TO_I5[name]

            analysis_work_dir = os.path.join(work_dir, appl, version)
            os.makedirs(analysis_work_dir, exist_ok=True)

            config = configs[analysis_id]
            factory = TaskFactory(uri=uri, i5_dir=analysis["i5_dir"],
                                  appl=appl, version=version,
                                  work_dir=analysis_work_dir,
                                  analysis_id=analysis_id,
                                  config=config,
                                  match_table=analysis["tables"]["matches"],
                                  parse_matches=parse_matches,
                                  site_table=analysis["tables"]["sites"],
                                  parse_sites=parse_sites,
                                  keep_files=keep_files,
                                  lsf_queue=lsf_queue,
                                  timeout=timeout)

            for upi_from, upi_to in incomplete_jobs.pop(analysis_id, []):
                task = factory.make(upi_from, upi_to)
                task.status = STATUS_RUNNING  # assumes it's running
                task.workdir = os.path.join(temp_dir, task.name)

                if task.name in name2id:
                    task.jobid = name2id.pop(task.name)
                else:
                    task.poll()  # checks if output exists
                    if not task.done():
                        task.status = STATUS_PENDING

                pool.submit(task)
                n_tasks += 1

            n_tasks_analysis = 0
            for upi_from, upi_to in range_jobs(next_upi, max_upi,
                                               config["job_size"]):
                if n_tasks_analysis == max_jobs_per_analysis:
                    break

                to_run.append(factory.make(upi_from, upi_to))
                database.add_job(cur, analysis_id, upi_from, upi_to)
                n_tasks += 1
                n_tasks_analysis += 1

            logger.debug(f"{name} {version}: "
                         f"{n_tasks_analysis} tasks submitted")

        cur.close()
        con.close()

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
        while to_run:
            i = random.randrange(len(to_run))  # 0 <= i < len(to_run)
            task = to_run.pop(i)
            pool.submit(task)

        logger.info(f"Running tasks: {n_tasks:,}")

        n_completed = n_failed = 0
        milestone = step = 5
        retries = {}
        while (n_completed + n_failed) < n_tasks:
            for task in pool.as_completed(wait=True):
                upi_from = task.args[1]
                upi_to = task.args[2]
                analysis_id = task.args[6]

                logfile = os.path.join(temp_dir, f"{task.name}.log")

                # Check how much memory was used
                max_mem = get_lsf_max_memory(task.stdout)

                if task.completed():
                    # Remove the log file (exists if a previous run failed)
                    try:
                        os.unlink(logfile)
                    except FileNotFoundError:
                        pass

                    # Flag the job as completed in the database
                    cpu_time = get_lsf_cpu_time(task.stdout)
                    database.update_job(uri, analysis_id, upi_from, upi_to,
                                        task, max_mem, cpu_time)

                    n_completed += 1
                else:
                    # Write to log file
                    with open(logfile, "wt") as fh:
                        fh.write(task.stdout)
                        fh.write(task.stderr)

                    # Number of times the task was re-submitted
                    num_retries = retries.get(task.name, 0)

                    # Did the job reached the memory usage limit?
                    mem_err = max_mem >= task.scheduler["mem"]

                    if num_retries < max_retries or (mem_err and infinite_mem):
                        # Task allowed to be re-submitted

                        # Increase memory requirement if needed
                        while task.scheduler["mem"] < max_mem:
                            task.scheduler["mem"] *= 1.5

                        # Resubmit task
                        pool.submit(task)

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


def export_fasta(uri: str, fasta_file: str, upi_from: str, upi_to: str) -> int:
    num_sequences = 0
    with open(fasta_file, "wt") as fh:
        con = cx_Oracle.connect(uri)
        cur = con.cursor()

        cur.execute(
            """
            SELECT UPI, SEQ_SHORT, SEQ_LONG
            FROM UNIPARC.PROTEIN
            WHERE UPI BETWEEN :1 AND :2
            """,
            (upi_from, upi_to)
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
           cpu: Optional[int] = None, temp_dir: Optional[str] = None,
           timeout: Optional[int] = None) -> tuple[bool, str, str]:
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

    process = subprocess.run(args, capture_output=True, timeout=timeout)
    return (
        process.returncode == 0,
        process.stdout.decode("utf-8"),
        process.stderr.decode("utf-8")
    )


def run_job(uri: str, upi_from: str, upi_to: str, i5_dir: str, appl: str,
            outdir: str, analysis_id: int, match_table: str,
            parse_matches: Callable, site_table: Optional[str],
            parse_sites: Optional[Callable], cpu: Optional[int] = None,
            keep_files: Optional[str] = None, timeout: Optional[int] = None):
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
            ok, stdout, stderr = run_i5(i5_dir, fasta_file, appl,
                                        matches_output, cpu=cpu,
                                        temp_dir=outdir, timeout=timeout)

            # Write captured streams
            sys.stdout.write(stdout)
            sys.stderr.write(stderr)

            if not ok or _I5_SUCCESS not in stdout:
                raise RuntimeError("InterProScan error")
            elif not os.path.isfile(matches_output):
                raise RuntimeError(f"Cannot access output matches tsv-pro")
            elif site_table and not os.path.isfile(sites_output):
                raise RuntimeError(f"Cannot access output sites tsv-pro")

            parse_matches(uri, matches_output, analysis_id, match_table)
            if site_table:
                parse_sites(uri, sites_output, analysis_id, site_table)
    except Exception:
        raise
    else:
        ok = True
    finally:
        if keep_files != "all" and (keep_files != "failed" or ok):
            shutil.rmtree(outdir)


def get_lsf_max_memory(stdout: str) -> int:
    match = re.search(r"^\s*Max Memory :\s+(\d+\sMB|-)$", stdout, re.M)
    group = match.group(1)
    return 0 if group == "-" else int(group.split()[0])


def get_lsf_cpu_time(stdout: str) -> int:
    match = re.search(r"^\s*CPU time :\s+(\d+)\.\d+ sec.$", stdout, re.M)
    return int(match.group(1))


def get_lsf_unfinished_jobs() -> dict[str, int]:
    stdout = subprocess.run(["bjobs", "-o", "jobid name"],
                            capture_output=True,
                            encoding="utf-8").stdout
    jobs = {}
    if "No unfinished job found" in stdout:
        return jobs

    for line in stdout.split("\n")[1:]:
        if line:
            job_id, name = line.split()
            jobs[name] = int(job_id)

    return jobs


# def kill_lsf_job(name: str) -> bool:
#     process = subprocess.run(["bkill", "-r", "-J", name],
#                              stdout=subprocess.DEVNULL,
#                              stderr=subprocess.DEVNULL)
#     return process.returncode == 0
