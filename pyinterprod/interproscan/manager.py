import os
import shutil
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from logging import DEBUG

import oracledb
from mundone import Pool, Task
from mundone.states import PENDING, RUNNING

from pyinterprod import logger
from pyinterprod.uniprot.uniparc import int_to_upi, upi_to_int, range_upi
from pyinterprod.utils import oracle
from . import analyses, jobs


"""
Key   -> analysis name in the database
Value -> tuple of three items
            - analysis name in InterProScan
            - persistence function for matches (.tsv-pro)
            - persistence function for sites (.tsv-pro.sites)
"""
_DB_TO_I5 = {
    "AntiFam": ("AntiFam", analyses.hmmer3_matches, None),
    "CATH-Gene3D": ("Gene3D", analyses.hmmer3_matches, None),
    "CDD": ("CDD", analyses.cdd_matches, analyses.sites),
    "COILS": ("Coils", analyses.coils_phobius_matches, None),
    "FunFam": ("FunFam", analyses.funfam_matches, None),
    "HAMAP": ("Hamap", analyses.hamap_matches, None),
    "MobiDB Lite": ("MobiDBLite", analyses.mobidb_lite_matches, None),
    "NCBIfam": ("NCBIfam", analyses.hmmer3_matches, None),
    "PANTHER": ("PANTHER", analyses.panther_matches, None),
    "Pfam": ("Pfam", analyses.hmmer3_matches, None),
    "Phobius": ("Phobius", analyses.coils_phobius_matches, None),
    "PIRSF": ("PIRSF", analyses.hmmer3_matches, None),
    "PIRSR": ("PIRSR", analyses.pirsr_matches, analyses.sites),
    "PRINTS": ("PRINTS", analyses.prints_matches, None),
    "PROSITE patterns": ("ProSitePatterns",
                         analyses.prosite_patterns_matches, None),
    "PROSITE profiles": ("ProSiteProfiles",
                         analyses.prosite_profiles_matches, None),
    "SFLD": ("SFLD", analyses.sfld_matches, analyses.sites),
    "SignalP_Euk": ("SignalP_EUK", analyses.signalp_matches, None),
    "SignalP_Gram_positive": ("SignalP_GRAM_POSITIVE",
                              analyses.signalp_matches, None),
    "SignalP_Gram_negative": ("SignalP_GRAM_NEGATIVE",
                              analyses.signalp_matches, None),
    "SMART": ("SMART", analyses.smart_matches, None),
    "SUPERFAMILY": ("SUPERFAMILY", analyses.superfamily_matches, None),
    "TIGRFAMs": ("TIGRFAM", analyses.hmmer3_matches, None),
    "TMHMM": ("TMHMM", analyses.tmhmm_matches, None),
}

# Input/output file names
_INPUT_FASTA = "input.fa"
_OUTPUT_MATCHES = "output.tsv-pro"
_OUTPUT_SITES = f"{_OUTPUT_MATCHES}.sites"


def sanitize_name(string: str) -> str:
    for c in [" ", "-", "_"]:
        string = string.replace(c, "")

    return string.lower()


@dataclass
class TaskFactory:
    i5_dir: str
    work_dir: str
    appl: str
    version: str
    config: dict
    has_sites: bool
    scheduler: str | None = None
    queue: str | None = None

    def make(self, upi_from: str, upi_to: str) -> Task:
        run_dir = self.get_run_dir(upi_from, upi_to)

        if self.has_sites:
            sites_path = self.get_sites_path(run_dir)
        else:
            sites_path = None

        return Task(
            fn=jobs.run_job,
            args=(
                self.i5_dir,
                self.appl,
                self.get_fasta_path(upi_from, upi_to),
                self.get_matches_path(run_dir),
                sites_path
            ),
            kwargs=dict(cpu=self.config["job_cpu"]),
            name="_".join([
                "IPM",
                self.appl,
                self.version,
                upi_from,
                upi_to
            ]),
            scheduler=dict(type=self.scheduler,
                           queue=self.queue,
                           cpu=self.config["job_cpu"],
                           mem=self.config["job_mem"],
                           hours=self.config["job_timeout"]),
            random_suffix=False
        )

    def get_run_dir(self, upi_from: str, upi_to: str,
                    mkdir: bool = False) -> str:
        name = f"{upi_from}_{upi_to}"
        path = os.path.join(self.work_dir, self.appl, self.version, name)

        if mkdir:
            try:
                shutil.rmtree(path)
            except FileNotFoundError:
                pass

            os.makedirs(path)

        return path

    def get_fasta_path(self, upi_from: str, upi_to: str) -> str:
        return os.path.join(self.get_run_dir(upi_from, upi_to), _INPUT_FASTA)

    @staticmethod
    def get_matches_path(run_dir: str) -> str:
        return os.path.join(run_dir, _OUTPUT_MATCHES)

    @staticmethod
    def get_sites_path(run_dir: str) -> str:
        return os.path.join(run_dir, _OUTPUT_SITES)


def run(uri: str, work_dir: str, temp_dir: str, **kwargs):
    base_config = {
        "job_cpu": kwargs.get("job_cpu", 8),
        "job_mem": kwargs.get("job_mem", 8 * 1024),  # MB
        "job_size": kwargs.get("job_size", 32000),
        "job_timeout": kwargs.get("job_timeout", 12),  # Hours
    }
    custom_configs = kwargs.get("config", {})
    debug = kwargs.get("debug", False)
    dry_run = kwargs.get("dry_run", False)
    auto_retry = kwargs.get("auto_retry", False)
    keep_files = kwargs.get("keep_files", None)
    max_retries = kwargs.get("max_retries", 0)
    max_timeout = kwargs.get("max_timeout", 120)
    max_running_jobs = kwargs.get("max_running_jobs", 1000)
    max_jobs_per_analysis = kwargs.get("max_jobs_per_analysis", -1)
    pool_threads = kwargs.get("pool_threads", 4)
    scheduler = kwargs.get("scheduler")
    queue = kwargs.get("queue")
    to_run = kwargs.get("analyses", [])
    to_exclude = kwargs.get("exclude", [])

    if debug:
        logger.setLevel(DEBUG)

    logger.info("starting")
    con = oracledb.connect(uri)
    cur = con.cursor()

    # Find analyses to run
    analyses_info = {}
    configs = {}
    name2id = {}

    for analysis_id, analysis in analyses.get_analyses(cur).items():
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
        analyses_info[analysis_id] = analysis

        # Apply default config
        configs[analysis_id] = base_config.copy()

    if not analyses_info:
        cur.close()
        con.close()
        logger.error("no analyses to process: exit")
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

    incomplete_jobs = jobs.get_incomplete_jobs(cur)
    cur.execute("SELECT MAX(UPI) FROM UNIPARC.PROTEIN")
    max_upi, = cur.fetchone()

    if scheduler.lower() == "lsf":
        name2id = jobs.get_unfinished_lsf_jobs()
    elif scheduler.lower() == "slurm":
        name2id = jobs.get_unfinished_slurm_jobs()
    else:
        raise ValueError(scheduler)

    with Pool(temp_dir, max_running_jobs,
              kill_on_exit=False,
              threads=pool_threads) as pool:
        tasks_info = {}

        with ThreadPoolExecutor(max_workers=pool_threads) as executor:
            fs = {}

            for analysis_id, analysis in analyses_info.items():
                analysis_name = analysis["name"]
                analysis_version = analysis["version"]
                has_sites = analysis["tables"]["sites"] is not None

                appl, _, _ = _DB_TO_I5[analysis_name]
                config = configs[analysis_id]
                factory = TaskFactory(i5_dir=analysis["i5_dir"],
                                      work_dir=work_dir,
                                      appl=appl, version=analysis_version,
                                      config=config,
                                      has_sites=has_sites,
                                      scheduler=scheduler, queue=queue)

                n_tasks_analysis = 0

                # First, look at existing jobs marked as incomplete in the DB
                for info in incomplete_jobs.pop(analysis_id, []):
                    upi_from, upi_to, is_running, num_sequences = info

                    if dry_run:
                        if 0 <= max_jobs_per_analysis <= n_tasks_analysis:
                            break
                        else:
                            # Doesn't matter what is in tasks_info (dry run)
                            tasks_info[(analysis_id, upi_from, upi_to)] = None
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
                            pool.submit(task)
                            tasks_info[task.name] = (
                                analysis_id, upi_from, upi_to, num_sequences,
                                factory.get_run_dir(upi_from, upi_to)
                            )
                            n_tasks_analysis += 1
                            continue

                        task.poll()  # Checks if output exists
                        if task.is_done():
                            """
                            Completed or failed. Submitted to the pool
                            which will send it back without restarting it.
                            """
                            pool.submit(task)
                            tasks_info[task.name] = (
                                analysis_id, upi_from, upi_to, num_sequences,
                                factory.get_run_dir(upi_from, upi_to)
                            )
                        elif 0 <= max_jobs_per_analysis <= n_tasks_analysis:
                            break
                        else:
                            task.status = PENDING
                            pool.submit(task)
                            tasks_info[task.name] = (
                                analysis_id, upi_from, upi_to, num_sequences,
                                factory.get_run_dir(upi_from, upi_to)
                            )
                            n_tasks_analysis += 1
                    elif 0 <= max_jobs_per_analysis <= n_tasks_analysis:
                        break
                    else:
                        # Flagged as failed in the database: re-submit
                        # Change status to pending
                        task.status = PENDING

                        # Re-export sequences
                        run_dir = factory.get_run_dir(upi_from, upi_to,
                                                      mkdir=True)
                        fasta_path = factory.get_fasta_path(upi_from, upi_to)
                        f = executor.submit(export_sequences, uri, upi_from,
                                            upi_to, fasta_path)
                        fs[f] = (analysis_id, upi_from, upi_to, run_dir, task)

                        n_tasks_analysis += 1

                # Now submit new jobs
                if analysis["max_upi"]:
                    next_upi = int_to_upi(upi_to_int(analysis["max_upi"]) + 1)
                else:
                    next_upi = int_to_upi(1)

                job_size = config["job_size"]
                for upi_from, upi_to in range_upi(next_upi, max_upi, job_size):
                    if 0 <= max_jobs_per_analysis <= n_tasks_analysis:
                        break
                    elif dry_run:
                        # Doesn't matter what is in tasks_info (dry run)
                        tasks_info[(analysis_id, upi_from, upi_to)] = None
                        n_tasks_analysis += 1
                        continue

                    # Export sequences
                    run_dir = factory.get_run_dir(upi_from, upi_to, mkdir=True)
                    fasta_path = factory.get_fasta_path(upi_from, upi_to)
                    task = factory.make(upi_from, upi_to)
                    f = executor.submit(export_sequences, uri, upi_from,
                                        upi_to, fasta_path)
                    fs[f] = (analysis_id, upi_from, upi_to, run_dir, task)
                    n_tasks_analysis += 1

                logger.debug(f"{analysis_name} {analysis_version}: "
                             f"{n_tasks_analysis} tasks")

            for f in as_completed(fs):
                analysis_id, upi_from, upi_to, run_dir, task = fs[f]

                try:
                    num_sequences = f.result()
                except Exception as exc:
                    logger.error(f"{analysis_id} ({upi_from}-{upi_to}): {exc}")
                    continue

                # Add a (placeholder/inactive) job
                jobs.add_job(cur, analysis_id, upi_from, upi_to, num_sequences)

                if num_sequences > 0:
                    # Submit the task
                    pool.submit(task)
                    tasks_info[task.name] = (analysis_id, upi_from, upi_to,
                                             num_sequences, run_dir)
                else:
                    # Empty job: flag it as successful
                    jobs.update_job(cur, analysis_id, upi_from, upi_to,
                                    success=True)
                    shutil.rmtree(run_dir, ignore_errors=True)

            del fs

        cur.close()
        con.close()

        n_tasks = len(tasks_info)
        logger.info(f"tasks: {n_tasks:,}")

        if dry_run:
            return

        logger.info(f"monitoring")
        con = oracle.try_connect(uri)
        cur = con.cursor()

        n_completed = n_failed = 0
        milestone = step = 5
        retries = {}
        while (n_completed + n_failed) < n_tasks:
            for task in pool.as_completed(wait=True):
                info = tasks_info[task.name]
                analysis_id, upi_from, upi_to, num_sequences, run_dir = info
                analysis = analyses_info[analysis_id]
                analysis_name = analysis["name"]
                matches_table = analysis["tables"]["matches"]
                sites_table = analysis["tables"]["sites"]

                _, fn_matches, fn_sites = _DB_TO_I5[analysis_name]

                logfile = os.path.join(temp_dir, f"{task.name}.log")
                if task.is_successful():
                    logger.debug(f"Completed: {analysis_name} ({analysis_id})"
                                 f" {upi_from}-{upi_to}")

                    # Persist data
                    cur2 = con.cursor()
                    analyses.persist_results(
                        cur2,
                        analysis_id,
                        fn_matches,
                        factory.get_matches_path(run_dir),
                        matches_table,
                        fn_sites,
                        factory.get_sites_path(run_dir),
                        sites_table
                    )
                    cur2.close()

                    # Update job in database
                    jobs.update_job(
                        cur, analysis_id, upi_from, upi_to,
                        task.submit_time, task.start_time, task.end_time,
                        task.maxmem, task.executor.memory,
                        task.cputime, True
                    )

                    if keep_files == "all":
                        with open(logfile, "wt") as fh:
                            fh.write(task.stdout)
                            fh.write(task.stderr)
                    else:
                        # Remove the log file (exists if a previous run failed)
                        try:
                            os.unlink(logfile)
                        except FileNotFoundError:
                            pass

                    shutil.rmtree(run_dir, ignore_errors=True)
                    n_completed += 1
                else:
                    logger.debug(f"Failed: {analysis_name} ({analysis_id})"
                                 f" {upi_from}-{upi_to}")

                    if keep_files in ("all", "failed"):
                        with open(logfile, "wt") as fh:
                            fh.write(task.stdout)
                            fh.write(task.stderr)

                    # Number of times the task was re-submitted
                    num_retries = retries.get(task.name, 0)

                    # Did the job reached the memory usage limit?
                    mem_err = task.is_oom()

                    # Did the job reached the timeout limit?
                    time_err = False
                    task_limit = None
                    starttime, endtime = task.executor.get_times(task.stdout)
                    if (starttime is not None and
                            endtime is not None and
                            task.executor.limit is not None):
                        runtime = (endtime - starttime).total_seconds() / 3600
                        task_limit = task.executor.limit.total_seconds() / 3600
                        time_err = runtime >= task_limit

                    if ((auto_retry and (mem_err or time_err)) or
                            num_retries < max_retries):
                        # Task allowed to be re-submitted

                        # Increase hours if time limit reached
                        if time_err and (task_limit * 1.25 < max_timeout):
                            task.executor.limit *= 1.25

                        if mem_err:
                            # Increase memory requirement
                            maxmem = task.maxmem
                            try:
                                while True:
                                    task.executor.memory *= 1.5
                                    if task.executor.memory > maxmem:
                                        break
                            except TypeError:
                                pass

                        # Resubmit task
                        task.status = PENDING
                        pool.submit(task)

                        # Add new job
                        jobs.add_job(cur, analysis_id, upi_from, upi_to,
                                     num_sequences)

                        # Increment retries counter
                        retries[task.name] = num_retries + 1
                    else:
                        # Max number of retries reached
                        n_failed += 1
                        shutil.rmtree(run_dir, ignore_errors=True)

                progress = (n_completed + n_failed) * 100 / n_tasks
                if progress >= milestone:
                    while progress >= milestone:
                        milestone += step

                    logger.info(f"progress: {progress:>3.0f}%")

        if n_failed:
            logger.error(f"{n_failed} task(s) failed")
        else:
            logger.info("complete")


def export_sequences(uri: str, upi_from: str, upi_to: str, output: str) -> int:
    con = oracle.connect(uri)
    cur = con.cursor()
    cur.execute(
        """
        SELECT UPI, SEQ_SHORT, SEQ_LONG
        FROM UNIPARC.PROTEIN
        WHERE UPI BETWEEN :1 AND :2
        """,
        [upi_from, upi_to]
    )
    rows = cur.fetchall()
    num_sequences = len(rows)
    if num_sequences > 0:
        with open(output, "wt") as fh:
            for upi, seq_short, seq_long in rows:
                sequence = seq_short or seq_long.read()

                fh.write(f">{upi}\n")
                for i in range(0, len(sequence), 60):
                    fh.write(f"{sequence[i:i + 60]}\n")

    cur.close()
    con.close()
    return num_sequences
