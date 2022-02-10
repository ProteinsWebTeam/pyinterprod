import os
import re
import shutil
import subprocess
import time
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


def run(uri: str, work_dir: str, temp_dir: str, lsf_queue: str, **kwargs):
    analysis_ids = kwargs.get("analyses", [])
    job_cpu = kwargs.get("job_cpu", 8)
    job_mem = kwargs.get("job_mem", 8000)
    max_running_jobs = kwargs.get("max_running_jobs", 1000)
    max_jobs = kwargs.get("max_jobs", 0)

    con = cx_Oracle.connect(uri)
    cur = con.cursor()

    analyses = {}
    for analysis_id, analysis in database.get_analyses(cur).items():
        if not analysis_ids or analysis_id in analysis_ids:
            analyses[analysis_id] = analysis

    if not analyses:
        cur.close()
        con.close()
        logger.error("No analyses to process: exit")

    incomplete_jobs = database.get_incomplete_jobs(cur)

    with Pool(temp_dir, max_running_jobs) as pool:
        n_tasks = 0
        for analysis_id, analysis in analyses.items():
            name = analysis["name"]
            version = analysis["version"]
            max_upi = analysis["max_upi"]
            i5_dir = analysis["i5_dir"]
            match_table = analysis["tables"]["matches"]
            site_table = analysis["tables"]["sites"]
            appl, parse_matches, parse_sites = _DB_TO_I5[name]

            to_restart = incomplete_jobs.get(analysis_id, [])
            killed_jobs = 0

            for upi_from, upi_to in to_restart:
                job_name = f"IPM_{appl}_{version}_{upi_from}_{upi_to}"
                if kill_lsf_job(job_name):
                    killed_jobs += 1

            if killed_jobs:
                time.sleep(5)

            n_tasks_analysis = 0
            for upi_from, upi_to in to_restart:
                task = Task(
                    fn=run_job,
                    args=(
                        uri,
                        upi_from,
                        upi_to,
                        i5_dir,
                        appl,
                        os.path.join(work_dir, appl, version),
                        analysis_id,
                        match_table,
                        parse_matches,
                        site_table,
                        parse_sites
                    ),
                    name=f"IPM_{appl}_{version}_{upi_from}_{upi_to}",
                    scheduler=dict(cpu=job_cpu, mem=job_mem, queue=lsf_queue),
                    random_suffix=False
                )
                pool.submit(task)
                n_tasks += 1
                n_tasks_analysis += 1

                if n_tasks == max_jobs:
                    break

            if max_jobs and n_tasks == max_jobs:
                logger.debug(f"{name} {version}: "
                             f"{n_tasks_analysis} tasks submitted")
                continue

            for upi_from, upi_to in database.get_runnable_jobs(cur, max_upi):
                task = Task(
                    fn=run_job,
                    args=(
                        uri,
                        upi_from,
                        upi_to,
                        i5_dir,
                        appl,
                        os.path.join(work_dir, appl, version),
                        analysis_id,
                        match_table,
                        parse_matches,
                        site_table,
                        parse_sites
                    ),
                    name=f"IPM_{appl}_{version}_{upi_from}_{upi_to}",
                    scheduler=dict(cpu=job_cpu, mem=job_mem, queue=lsf_queue),
                    random_suffix=False
                )
                pool.submit(task)
                database.add_job(cur, analysis_id, upi_from, upi_to)
                n_tasks += 1
                n_tasks_analysis += 1

                if n_tasks == max_jobs:
                    break

            logger.debug(f"{name} {version}: "
                         f"{n_tasks_analysis} tasks submitted")

        cur.close()
        con.close()

        logger.info(f"Tasks submitted: {n_tasks:,}")

        n_completed = 0
        progress = 0
        step = 5
        while n_completed < n_tasks:
            for task in pool.as_completed(wait=True):
                logger.debug(f"{task.name}: {task.state}")

                upi_from = task.args[1]
                upi_to = task.args[2]
                analysis_id = task.args[6]

                logfile = os.path.join(temp_dir, f"{task.name}.log")
                try:
                    os.unlink(logfile)
                except FileNotFoundError:
                    pass

                max_mem = get_lsf_max_memory(task.stdout)

                if task.completed():
                    n_completed += 1

                    database.update_job(uri, analysis_id, upi_from, upi_to,
                                        task, max_mem)
                else:
                    with open(logfile, "wt") as fh:
                        fh.write(task.stdout)
                        fh.write(task.stderr)

                    if max_mem >= task.scheduler["mem"]:
                        # Resubmit task with more memory
                        while task.scheduler["mem"] < max_mem:
                            task.scheduler["mem"] *= 1.5

                        pool.submit(task)
                    else:
                        # TODO: decide what to do
                        n_completed += 1

                if (n_completed * 100 / n_tasks) >= (progress + step):
                    logger.info(f"Progress: {progress + step:>3}%")
                    progress += step

        logger.info("complete")


def run_job(uri: str, upi_from: str, upi_to: str, i5_dir: str,
            analysis_name: str, work_dir: str, analysis_id: int,
            match_table: str, parse_matches: Callable,
            site_table: Optional[str], parse_sites: Optional[Callable]):
    temp_dir = os.path.join(work_dir, f"{upi_from}_{upi_to}")
    try:
        shutil.rmtree(temp_dir)
    except FileNotFoundError:
        pass

    os.makedirs(temp_dir)

    fasta_file = os.path.join(temp_dir, "input.fasta")

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

    matches_output = os.path.join(temp_dir, "output.tsv-pro")
    sites_output = matches_output + ".sites"
    args = [
        os.path.join(i5_dir, "interproscan.sh"),
        "-i", fasta_file,
        "-T", temp_dir,
        "-appl", analysis_name,
        "-dp",
        "-f", "tsv-pro",
        "-o", matches_output
    ]

    process = subprocess.run(args)

    ok = True
    if process.returncode != 0:
        ok = False
    elif not os.path.isfile(matches_output):
        ok = False
    elif site_table is not None and not os.path.isfile(sites_output):
        ok = False
    else:
        parse_matches(uri, matches_output, analysis_id, match_table)

        if site_table is not None:
            parse_sites(uri, sites_output, analysis_id, site_table)

    shutil.rmtree(temp_dir)

    if not ok:
        raise RuntimeError()


def int_to_upi(i):
    return f"UPI{i:010x}".upper()


def upi_to_int(upi):
    return int(upi[3:], 16)


def get_lsf_max_memory(stdout: str) -> int:
    match = re.search(r"^\s*Max Memory :\s+(\d+\sMB|-)$", stdout, re.M)
    group = match.group(1)
    return 0 if group == "-" else int(group.split()[0])


def kill_lsf_job(name: str) -> bool:
    process = subprocess.run(["bkill", "-J", name],
                             stdout=subprocess.DEVNULL,
                             stderr=subprocess.DEVNULL)
    return process.returncode == 0
