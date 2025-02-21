import math
import os
import shutil
import time
from dataclasses import dataclass
from logging import DEBUG
from queue import Queue, Empty
from threading import Thread

import oracledb
from mundone import Task

from pyinterprod import logger
from pyinterprod.uniprot.uniparc import int_to_upi, upi_to_int, range_upi
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
    "CATH-FunFam": ("FunFam", analyses.funfam_matches, None),
    "CATH-Gene3D": ("Gene3D", analyses.hmmer3_matches, None),
    "CDD": ("CDD", analyses.cdd_matches, analyses.sites),
    "COILS": ("Coils", analyses.coils_phobius_matches, None),
    "HAMAP": ("Hamap", analyses.hamap_matches, None),
    "MobiDB Lite": ("MobiDBLite", analyses.mobidb_lite_matches, None),
    "NCBIFAM": ("NCBIFAM", analyses.hmmer3_matches, None),
    "PANTHER": ("PANTHER", analyses.panther_matches, None),
    "Pfam": ("Pfam", analyses.hmmer3_matches, None),
    "Phobius": ("Phobius", analyses.coils_phobius_matches, None),
    "PIRSF": ("PIRSF", analyses.hmmer3_matches, None),
    "PIRSR": ("PIRSR", analyses.pirsr_matches, analyses.sites),
    "PRINTS": ("PRINTS", analyses.prints_matches, None),
    "PROSITE patterns": ("ProSitePatterns", analyses.prosite_patterns_matches, None),
    "PROSITE profiles": ("ProSiteProfiles", analyses.prosite_profiles_matches, None),
    "SFLD": ("SFLD", analyses.sfld_matches, analyses.sites),
    "SignalP_Euk": ("SignalP_EUK", analyses.signalp_matches, None),
    "SignalP_Gram_positive": ("SignalP_GRAM_POSITIVE", analyses.signalp_matches, None),
    "SignalP_Gram_negative": ("SignalP_GRAM_NEGATIVE", analyses.signalp_matches, None),
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
class InterProScanTask(Task):
    analysis_id: int
    upi_from: str
    upi_to: str
    i5_dir: str
    work_dir: str
    appl: str
    version: str
    config: dict
    has_sites: bool
    scheduler: str | None = None
    queue: str | None = None

    def __post_init__(self):
        # Call the parent class's initializer
        super().__init__(
            fn=jobs.run_job,
            args=(
                self.i5_dir,
                self.appl,
                self.get_fasta_path(),
                self.get_matches_path(),
                self.get_sites_path(),
            ),
            kwargs=dict(cpu=self.config["job_cpu"]),
            name="_".join(["IPM", self.appl, self.version, self.upi_from, self.upi_to]),
            scheduler=dict(
                type=self.scheduler,
                queue=self.queue,
                cpu=self.config["job_cpu"],
                mem=self.config["job_mem"],
                hours=self.config["job_timeout"],
            ),
            random_suffix=False,
        )

    def mkdir(self):
        path = self.get_run_dir()

        try:
            shutil.rmtree(path)
        except FileNotFoundError:
            pass

        os.makedirs(path)

    def rmdir(self):
        path = self.get_run_dir()

        try:
            shutil.rmtree(path)
        except FileNotFoundError:
            pass

    def get_run_dir(self) -> str:
        return os.path.join(
            self.work_dir, self.appl, self.version, f"{self.upi_from}_{self.upi_to}"
        )

    def get_fasta_path(self) -> str:
        return os.path.join(self.get_run_dir(), _INPUT_FASTA)

    def get_matches_path(self) -> str:
        return os.path.join(self.get_run_dir(), _OUTPUT_MATCHES)

    def get_sites_path(self) -> str | None:
        if self.has_sites:
            return os.path.join(self.get_run_dir(), _OUTPUT_SITES)
        else:
            return None


def run(uri: str, work_dir: str, temp_dir: str, **kwargs):
    base_config = {
        "job_cpu": kwargs.get("job_cpu", 8),
        "job_mem": kwargs.get("job_mem", 8 * 1024),  # MB
        "job_size": kwargs.get("job_size", 32000),
        "job_timeout": kwargs.get("job_timeout", 12),  # Hours
    }
    custom_configs = kwargs.get("config", {})
    debug = kwargs.get("debug", False)
    auto_retry = kwargs.get("auto_retry", False)
    keep_files = kwargs.get("keep_files", None)
    max_retries = kwargs.get("max_retries", 0)
    max_timeout = kwargs.get("max_timeout", 120)
    max_running_jobs = kwargs.get("max_running_jobs", 1000)
    max_jobs_per_analysis = kwargs.get("max_jobs_per_analysis", -1)
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

        # Either we want to run all analyses (default),
        # or we explicitly specified the current analysis
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

    # Find unfinished (pending + running) jobs
    if scheduler.lower() == "lsf":
        name2id = jobs.get_unfinished_lsf_jobs()
    elif scheduler.lower() == "slurm":
        name2id = jobs.get_unfinished_slurm_jobs()
    else:
        raise ValueError(scheduler)

    """
    Find (in the database) jobs that are either:
        - flagged as failed 
        - flagged as incomplete
    """
    incomplete_jobs = jobs.get_incomplete_jobs(cur)

    cur.execute("SELECT MAX(UPI) FROM UNIPARC.PROTEIN")
    (max_upi,) = cur.fetchone()
    cur.close()
    con.close()

    tasks = []  # list of tuple (task, is_new, num_sequences)
    num_jobs_per_analysis = {}

    for analysis_id, analysis in analyses_info.items():
        analysis_name = analysis["name"]
        analysis_version = analysis["version"]
        has_sites = analysis["tables"]["sites"] is not None
        num_jobs_per_analysis[analysis_id] = 0

        appl, _, _ = _DB_TO_I5[analysis_name]
        config = configs[analysis_id]

        # Check jobs that were already submitted
        for upi_from, upi_to, num_sequences in incomplete_jobs.pop(analysis_id, []):
            task = InterProScanTask(
                analysis_id=analysis_id,
                upi_from=upi_from,
                upi_to=upi_to,
                i5_dir=analysis["i5_dir"],
                work_dir=work_dir,
                appl=appl,
                version=analysis_version,
                config=config,
                has_sites=has_sites,
                scheduler=scheduler,
                queue=queue,
            )

            task.workdir = os.path.join(temp_dir, task.name)

            # Assumes task is running
            task.set_running()

            if task.name in name2id:
                # There is an associated job running on the cluster
                task.executor.id = name2id.pop(task.name)
                tasks.append((task, False, num_sequences))
                continue

            # Checks if output exists
            task.poll()
            if task.is_done():
                # Completed, failed, or cancelled, but at least it ran
                tasks.append((task, False, num_sequences))
            else:
                # Unknown fate: resubmit
                task.set_pending()
                tasks.append((task, True, num_sequences))

        # Submit new jobs
        if analysis["max_upi"]:
            next_upi = int_to_upi(upi_to_int(analysis["max_upi"]) + 1)
        else:
            next_upi = int_to_upi(1)

        job_size = config["job_size"]
        for upi_from, upi_to in range_upi(next_upi, max_upi, job_size):
            task = InterProScanTask(
                analysis_id=analysis_id,
                upi_from=upi_from,
                upi_to=upi_to,
                i5_dir=analysis["i5_dir"],
                work_dir=work_dir,
                appl=appl,
                version=analysis_version,
                config=config,
                has_sites=has_sites,
                scheduler=scheduler,
                queue=queue,
            )
            tasks.append((task, True, 0))

    if max_jobs_per_analysis >= 0:
        tmp_tasks = []
        for task, is_new, num_sequences in tasks:
            analysis_id = task.analysis_id
            if num_jobs_per_analysis[analysis_id] < max_jobs_per_analysis:
                tmp_tasks.append((task, is_new, num_sequences))
                num_jobs_per_analysis[analysis_id] += 1

        tasks = tmp_tasks

    num_tasks = len(tasks)
    logger.info(f"tasks: {num_tasks}")
    if max_running_jobs == 0 or num_tasks == 0:
        return

    # Create pool of FASTA export workers
    fasta_workers = []
    fasta_queue = Queue()
    submit_queue = Queue()
    for _ in range(8):
        t = Thread(
            target=export_sequences_worker, args=(uri, fasta_queue, submit_queue)
        )
        t.start()
        fasta_workers.append(t)

    task_sequences = {}
    while tasks:
        task, is_new, num_sequences = tasks.pop(0)
        fasta_queue.put((task, is_new, num_sequences))
        task_sequences[task.name] = num_sequences

    # Add sentinel value to terminate pool when all FASTA have been exported
    for _ in fasta_workers:
        fasta_queue.put(None)

    # Start pool of task monitor workers
    task_workers = []
    collect_queue = Queue()
    for _ in range(8):
        t = Thread(
            target=task_worker,
            args=(
                submit_queue,
                math.floor(max_running_jobs / 8),
                temp_dir,
                collect_queue,
            ),
        )
        t.start()
        task_workers.append(t)

    # Persist successful tasks
    con = oracledb.connect(uri)
    cur = con.cursor()
    num_completed = num_failed = 0
    milestone = step = 5
    retries = {}

    while (num_completed + num_failed) < num_tasks:
        task, num_sequences = collect_queue.get()
        collect_queue.task_done()

        logfile = os.path.join(temp_dir, f"{task.name}.log")
        failed = False

        if task.is_successful():
            logger.debug(f"{task.name}: completed")
            analysis = analyses_info[task.analysis_id]
            analysis_name = analysis["name"]
            matches_table = analysis["tables"]["matches"]
            sites_table = analysis["tables"]["sites"]
            _, fn_matches, fn_sites = _DB_TO_I5[analysis_name]

            # Persist data
            cur2 = con.cursor()
            ok = analyses.persist_results(
                cur2,
                task.analysis_id,
                fn_matches,
                task.get_matches_path(),
                matches_table,
                fn_sites,
                task.get_sites_path(),
                sites_table,
            )
            cur2.close()

            if ok:
                # Data persisted successfully
                logger.debug(f"{task.name}: persisted")
                jobs.update_job(
                    cur,
                    task.analysis_id,
                    task.upi_from,
                    task.upi_to,
                    task.start_time,
                    task.end_time,
                    task.maxmem,
                    task.executor.memory,
                    task.cputime,
                    success=True,
                )

                if keep_files == "all":
                    with open(logfile, "wt") as fh:
                        fh.write(task.stdout)
                        fh.write(task.stderr)
                else:
                    # Remove the log file
                    try:
                        os.unlink(logfile)
                    except FileNotFoundError:
                        pass

                    try_rmtree(task.get_run_dir())

                num_completed += 1
            else:
                # Persistence error (duplicated matches in I5 output)
                logger.debug(f"{task.name}: not persisted")
                failed = True
        else:
            logger.debug(f"{task.name}: failed")
            failed = True

        if failed:
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
            if (
                starttime is not None
                and endtime is not None
                and task.executor.limit is not None
            ):
                runtime = (endtime - starttime).total_seconds() / 3600
                task_limit = task.executor.limit.total_seconds() / 3600
                time_err = runtime >= task_limit

            if (auto_retry and (mem_err or time_err)) or num_retries < max_retries:
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
                task.set_pending()
                submit_queue.put((task, num_sequences))

                # Add new job
                num_sequences = task_sequences[task.name]
                jobs.add_job(
                    cur,
                    task.analysis_id,
                    task.upi_from,
                    task.upi_to,
                    num_sequences,
                )

                # Increment retries counter
                retries[task.name] = num_retries + 1
                logger.debug(f"{task.name}: re-submitted")
            else:
                # Max number of retries reached
                num_failed += 1

                if keep_files not in ("all", "failed"):
                    try_rmtree(task.get_run_dir())

        progress = (num_completed + num_failed) * 100 / num_tasks
        if progress >= milestone:
            while progress >= milestone:
                milestone += step

            logger.info(f"progress: {progress:>3.0f}%")

    cur.close()
    con.close()

    for _ in task_workers:
        submit_queue.put(None)

    fasta_queue.join()
    submit_queue.join()
    collect_queue.join()

    if num_failed:
        logger.error(f"{num_failed} task(s) failed")
    else:
        logger.info("complete")


def export_sequences_worker(uri: str, inqueue: Queue, outqueue: Queue):
    con = oracledb.connect(uri)
    cur = con.cursor()

    try:
        while True:
            item = inqueue.get()
            if item is None:
                # Sentinel
                inqueue.task_done()
                break
            else:
                task, is_new, num_sequences = item

                if is_new:
                    # New task: we need to export the input FASTA file
                    task.mkdir()
                    num_sequences = export_sequences(
                        cur, task.upi_from, task.upi_to, task.get_fasta_path()
                    )

                    # Add a (placeholder/inactive) job
                    jobs.add_job(
                        cur, task.analysis_id, task.upi_from, task.upi_to, num_sequences
                    )

                    if num_sequences > 0:
                        outqueue.put((task, num_sequences))
                    else:
                        # No sequences: it won't be necessary to run the task
                        task.rmdir()
                        task.set_successful()
                        jobs.update_job(
                            cur,
                            task.analysis_id,
                            task.upi_from,
                            task.upi_to,
                            success=True,
                        )

                else:
                    # Not new task: we assume the input file already exists
                    outqueue.put((task, num_sequences))

                inqueue.task_done()
    finally:
        cur.close()
        con.close()


def export_sequences(
    cur: oracledb.Cursor, upi_from: str, upi_to: str, output: str
) -> int:
    cur.execute(
        """
        SELECT UPI, SEQ_SHORT, SEQ_LONG
        FROM UNIPARC.PROTEIN
        WHERE UPI BETWEEN :1 AND :2
        """,
        [upi_from, upi_to],
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

    return num_sequences


def task_worker(inqueue: Queue, max_running: int, workdir: str, outqueue: Queue):
    pending = []
    running = []
    while True:
        try:
            item = inqueue.get(block=True, timeout=1)
        except Empty:
            pass  # Queue is empty
        else:
            if item is None:
                # Sentinel: pending and running should be empty
                inqueue.task_done()
                if pending or running:
                    logger.warning(f"{len(pending)} and {len(running)} tasks")
                break
            else:
                # Task to submit
                task, num_sequences = item
                pending.append((task, num_sequences))
                inqueue.task_done()

        # Monitor running tasks
        tmp_running = []
        for task, num_sequences in running:
            task.poll()

            if task.is_running():
                # Still running
                tmp_running.append((task, num_sequences))
            else:
                # Done (successfully or not)
                outqueue.put((task, num_sequences))

                if pending:
                    # Free slot: replace with pending task
                    task, num_sequences = pending.pop(0)
                    if task.start(dir=workdir):
                        tmp_running.append((task, num_sequences))
                    else:
                        pending.append((task, num_sequences))

        # Fill free slots with pending tasks
        running = tmp_running
        while pending and len(running) < max_running:
            task, num_sequences = pending.pop(0)
            if task.start(dir=workdir):
                tmp_running.append((task, num_sequences))
            else:
                pending.append((task, num_sequences))
                break


def try_rmtree(path: str, max_attempts: int = 3):
    num_attempts = 1
    while num_attempts <= max_attempts:
        try:
            shutil.rmtree(path)
        except Exception as exc:
            time.sleep(1)
            num_attempts += 1
            if num_attempts == max_attempts:
                logger.error(exc)
        else:
            break
