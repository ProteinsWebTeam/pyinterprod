import os
import re
import shutil
import subprocess
import time
from dataclasses import dataclass
from typing import Callable, Optional

import cx_Oracle
from mundone import Pool, Task

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

    def make(self, upi_from: str, upi_to: str) -> Task:
        name = f"{self.appl}_{self.version}_{upi_from}_{upi_to}"

        return Task(
            fn=run_job,
            args=(
                self.uri,
                upi_from,
                upi_to,
                self.i5_dir,
                self.appl,
                os.path.join(self.work_dir, name),
                self.analysis_id,
                self.match_table,
                self.parse_matches,
                self.site_table,
                self.parse_sites
            ),
            kwargs=dict(timeout=self.config["timeout"]),
            name=f"IPM_{name}",
            scheduler=dict(cpu=self.config["job_cpu"],
                           mem=self.config["job_mem"],
                           queue=self.config["lsf_queue"]),
            random_suffix=False
        )


def run(uri: str, work_dir: str, temp_dir: str, **kwargs):
    base_config = {
        "job_cpu": kwargs.get("job_cpu", 8),
        "job_mem": kwargs.get("job_mem", 8 * 1024),
        "job_size": kwargs.get("job_size", 100000),
        "lsf_queue": kwargs.get("lsf_queue"),
        "timeout": kwargs.get("timeout")
    }
    custom_configs = kwargs.get("config", {})
    max_retries = kwargs.get("max_retries", 0)
    max_running_jobs = kwargs.get("max_running_jobs", 1000)
    max_jobs_per_analysis = kwargs.get("max_jobs_per_analysis", 0)
    pool_threads = kwargs.get("pool_threads", 4)
    to_run = kwargs.get("analyses", [])

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

        if to_run and analysis_id not in to_run:
            continue

        analyses[analysis_id] = analysis

        # Default config
        configs[analysis_id] = base_config.copy()

    if not analyses:
        cur.close()
        con.close()
        logger.error("No analyses to process: exit")

    # Override with custom configs
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

    with Pool(temp_dir, max_running_jobs, pool_threads) as pool:
        n_tasks = 0
        for analysis_id, analysis in analyses.items():
            name = analysis["name"]
            version = analysis["version"]
            max_upi = analysis["max_upi"]
            appl, parse_matches, parse_sites = _DB_TO_I5[name]

            to_restart = incomplete_jobs.get(analysis_id, [])
            killed_jobs = 0

            for upi_from, upi_to in to_restart:
                job_name = f"IPM_{appl}_{version}_{upi_from}_{upi_to}"
                if kill_lsf_job(job_name):
                    killed_jobs += 1

            if killed_jobs:
                time.sleep(5)

            config = configs[analysis_id]
            factory = TaskFactory(uri=uri, i5_dir=analysis["i5_dir"],
                                  appl=appl, version=version,
                                  work_dir=work_dir, analysis_id=analysis_id,
                                  config=config,
                                  match_table=analysis["tables"]["matches"],
                                  parse_matches=parse_matches,
                                  site_table=analysis["tables"]["sites"],
                                  parse_sites=parse_sites)

            n_tasks_analysis = 0
            for upi_from, upi_to in to_restart:
                pool.submit(factory.make(upi_from, upi_to))
                n_tasks += 1
                n_tasks_analysis += 1

                if n_tasks_analysis == max_jobs_per_analysis:
                    break

            if (max_jobs_per_analysis
                    and n_tasks_analysis == max_jobs_per_analysis):
                logger.debug(f"{name} {version}: "
                             f"{n_tasks_analysis} tasks submitted")
                continue

            for upi_from, upi_to in database.get_jobs(cur, config["job_size"],
                                                      max_upi):
                pool.submit(factory.make(upi_from, upi_to))
                database.add_job(cur, analysis_id, upi_from, upi_to)
                n_tasks += 1
                n_tasks_analysis += 1

                if n_tasks_analysis == max_jobs_per_analysis:
                    break

            logger.debug(f"{name} {version}: "
                         f"{n_tasks_analysis} tasks submitted")

        cur.close()
        con.close()

        logger.info(f"Tasks submitted: {n_tasks:,}")

        n_completed = n_failed = progress = 0
        step = 5
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
                    database.update_job(uri, analysis_id, upi_from, upi_to,
                                        task, max_mem)

                    n_completed += 1
                else:
                    # Write to log file (with this run's output/error)
                    with open(logfile, "at") as fh:
                        fh.write(task.stdout)
                        fh.write(task.stderr)

                    if retries.get(task.name, 0) == max_retries:
                        # Max number of retries reached
                        n_failed += 1
                        continue
                    elif max_mem >= task.scheduler["mem"]:
                        # Increase memory requirement
                        while task.scheduler["mem"] < max_mem:
                            task.scheduler["mem"] *= 1.5

                    # Resubmit task
                    pool.submit(task)
                    retries[task.name] = retries.get(task.name, 0) + 1

                pc = (n_completed + n_failed) * 100 / n_tasks
                if pc >= (progress + step):
                    logger.info(f"Progress: {progress + step:>3}%")
                    progress += step

        if n_failed:
            logger.error(f"{n_failed} task(s) failed")
        else:
            logger.info("complete")


def export_fasta(uri: str, fasta_file: str, upi_from: str, upi_to: str):
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

        cur.close()
        con.close()


def run_i5(i5_dir: str, fasta_file: str, analysis_name: str, output: str,
           temp_dir: Optional[str] = None,
           timeout: Optional[int] = None) -> bool:
    if temp_dir is None:
        temp_dir = os.path.dirname(output)

    args = [
        os.path.join(i5_dir, "interproscan.sh"),
        "-i", fasta_file,
        "-T", temp_dir,
        "-appl", analysis_name,
        "-dp",
        "-f", "tsv-pro",
        "-o", output
    ]

    process = subprocess.run(args, timeout=timeout)
    return process.returncode == 0 and os.path.isfile(output)


def run_job(uri: str, upi_from: str, upi_to: str, i5_dir: str, appl: str,
            outdir: str, analysis_id: int, match_table: str,
            parse_matches: Callable, site_table: Optional[str],
            parse_sites: Optional[Callable], timeout: Optional[int] = None):

    try:
        shutil.rmtree(outdir)
    except FileNotFoundError:
        pass

    os.makedirs(outdir)
    fasta_file = os.path.join(outdir, "input.fasta")
    matches_output = os.path.join(outdir, "output.tsv-pro")
    sites_output = matches_output + ".sites"

    try:
        export_fasta(uri, fasta_file, upi_from, upi_to)

        if not run_i5(i5_dir, fasta_file, appl, matches_output,
                      timeout=timeout):
            raise RuntimeError()

        if site_table and not os.path.isfile(sites_output):
            raise RuntimeError()

        parse_matches(uri, matches_output, analysis_id, match_table)
        if site_table:
            parse_sites(uri, sites_output, analysis_id, site_table)
    except Exception:
        raise
    finally:
        shutil.rmtree(outdir)


def get_lsf_max_memory(stdout: str) -> int:
    match = re.search(r"^\s*Max Memory :\s+(\d+\sMB|-)$", stdout, re.M)
    group = match.group(1)
    return 0 if group == "-" else int(group.split()[0])


def kill_lsf_job(name: str) -> bool:
    process = subprocess.run(["bkill", "-J", name],
                             stdout=subprocess.DEVNULL,
                             stderr=subprocess.DEVNULL)
    return process.returncode == 0
