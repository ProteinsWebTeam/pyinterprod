#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import json
import os
import shutil
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from tempfile import mkdtemp, mkstemp
from typing import Optional, Sequence

import cx_Oracle

from pyinterprod import logger, orautils
from pyinterprod.clan import cdd, pfam, panther, pirsf, utils


def create_tables(cur: cx_Oracle.Cursor):
    # Ensure these tables are deleted (legacy)
    orautils.drop_table(cur, "INTERPRO", "METHOD_SET", purge=True)
    orautils.drop_table(cur, "INTERPRO", "METHOD_SCAN", purge=True)

    orautils.drop_table(cur, "INTERPRO", "CLAN_MEMBER", purge=True)
    orautils.drop_table(cur, "INTERPRO", "CLAN", purge=True)
    orautils.drop_table(cur, "INTERPRO", "CLAN_MEMBER_ALN", purge=True)

    cur.execute(
        """
        CREATE TABLE INTERPRO.CLAN
        (
            CLAN_AC VARCHAR2(25) NOT NULL,
            DBCODE CHAR(1) NOT NULL,
            NAME VARCHAR2(100) DEFAULT NULL,
            DESCRIPTION VARCHAR2(4000) DEFAULT NULL,
            CONSTRAINT PK_CLAN 
              PRIMARY KEY (CLAN_AC),
            CONSTRAINT FK_CLAN$DBCODE
              FOREIGN KEY (DBCODE)
              REFERENCES INTERPRO.CV_DATABASE (DBCODE)
              ON DELETE CASCADE 
        )
        """
    )

    cur.execute(
        """
        CREATE TABLE INTERPRO.CLAN_MEMBER
        (
            CLAN_AC VARCHAR2(25) NOT NULL,
            METHOD_AC VARCHAR2(25) NOT NULL,
            SCORE BINARY_DOUBLE NOT NULL,
            SEQ CLOB DEFAULT NULL,
            CONSTRAINT PK_CLAN_MEMBER
              PRIMARY KEY (CLAN_AC, METHOD_AC),
            CONSTRAINT UQ_CLAN_MEMBER$METHOD_AC 
              UNIQUE (METHOD_AC),
            CONSTRAINT FK_CLAN_MEMBER$CLAN_AC
              FOREIGN KEY (CLAN_AC)
              REFERENCES INTERPRO.CLAN (CLAN_AC)
              ON DELETE CASCADE,
            CONSTRAINT FK_CLAN_MEMBER$METHOD_AC
              FOREIGN KEY (METHOD_AC)
              REFERENCES INTERPRO.METHOD (METHOD_AC)
              ON DELETE CASCADE
        )
        """
    )

    cur.execute(
        """
        CREATE TABLE INTERPRO.CLAN_MEMBER_ALN
        (
            QUERY_AC VARCHAR2(25) NOT NULL,
            TARGET_AC VARCHAR2(25) NOT NULL,
            EVALUE FLOAT NOT NULL,
            DOMAINS CLOB NOT NULL,
            CONSTRAINT PK_CLAN_MEMBER_ALN
                PRIMARY KEY (QUERY_AC, TARGET_AC),
            CONSTRAINT FK_CLAN_MEMBER_ALN
              FOREIGN KEY (QUERY_AC)
              REFERENCES INTERPRO.METHOD (METHOD_AC)
              ON DELETE CASCADE
        )
        """
    )


def _insert_clans(con: cx_Oracle.Connection, clans: Sequence[dict], dbcode: str):
    cur = con.cursor()
    cur.execute("SELECT METHOD_AC FROM INTERPRO.METHOD "
                "WHERE DBCODE=:1", (dbcode,))
    members = {row[0] for row in cur}
    cur.execute("DELETE FROM INTERPRO.CLAN WHERE DBCODE = :1", (dbcode,))
    con.commit()
    cur.close()

    t1 = orautils.TablePopulator(con, "INSERT INTO INTERPRO.CLAN "
                                      "VALUES (:1, :2, :3, :4)")
    t2 = orautils.TablePopulator(con, "INSERT INTO INTERPRO.CLAN_MEMBER "
                                      "(CLAN_AC, METHOD_AC ,SCORE) "
                                      "VALUES (:1, :2, :3)", depends_on=t1)

    for c in clans:
        clan_members = []
        for member in c["members"]:
            if member["accession"] in members:
                clan_members.append(member)

        if not clan_members:
            continue

        t1.insert((c["accession"], dbcode, c["name"], c["description"]))
        for member in clan_members:
            t2.insert((c["accession"], member["accession"], member["score"]))

    t1.close()
    t2.close()


def update_clans(user: str, dsn: str, databases: Sequence[str], pfam_url: str):
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    create_tables(cur)
    cur.close()

    for mem_db in databases:
        logger.info(f"updating clans for {mem_db.upper()}")
        if mem_db == "cdd":
            _insert_clans(con, cdd.get_clans(), cdd.DBCODE)
        elif mem_db == "panther":
            _insert_clans(con, panther.get_clans(con), panther.DBCODE)
        elif mem_db == "pfam":
            _insert_clans(con, pfam.get_clans(pfam_url), pfam.DBCODE)
        elif mem_db == "pirsf":
            _insert_clans(con, pirsf.get_clans(), pirsf.DBCODE)
        else:
            raise ValueError(mem_db)

    con.commit()
    con.close()
    logger.info("done")


def align_hmm(user: str, dsn: str, mem_db: str, hmm_db: str, threads: int=1,
              tmpdir: Optional[str]=None, progress: bool=False):
    if mem_db == "panther":
        dbcode = panther.DBCODE
        x = 7
    elif mem_db == "pfam":
        dbcode = pfam.DBCODE
        x = 5
    elif mem_db == "pirsf":
        dbcode = pirsf.DBCODE
        x = 8
    else:
        raise ValueError(mem_db)

    if tmpdir:
        os.makedirs(tmpdir, exist_ok=True)

    workdir = mkdtemp(dir=tmpdir)

    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    cur.execute(
        """
        SELECT METHOD_AC FROM INTERPRO.CLAN_MEMBER
        WHERE CLAN_AC IN (SELECT CLAN_AC FROM INTERPRO.CLAN WHERE DBCODE = :1)
        """, (dbcode,)
    )
    members = {row[0] for row in cur}
    cur.execute(
        """
        DELETE 
        FROM INTERPRO.CLAN_MEMBER_ALN 
        WHERE QUERY_AC IN (
          SELECT METHOD_AC 
          FROM INTERPRO.METHOD 
          WHERE DBCODE = :1
        )
        """, (dbcode,)
    )
    con.commit()
    cur.close()
    con.close()

    logger.info("writing consensus sequences")
    entries = []
    for acc, hmm in utils.iter_hmm_database(hmm_db):
        if acc in members:
            subdir = os.path.join(workdir, acc[:x])
            try:
                os.mkdir(subdir)
            except FileExistsError:
                pass

            hmmfile = os.path.join(subdir, acc) + ".hmm"
            with open(hmmfile, "wt") as fh:
                fh.write(hmm)

            seqfile = os.path.join(subdir, acc) + ".fa"
            utils.hmmemit(hmmfile, seqfile)
            os.remove(hmmfile)

            entries.append((acc, seqfile))

    logger.info("querying sequences")
    with ThreadPoolExecutor(max_workers=threads) as executor:
        fs = {}
        for acc, seqfile in entries:
            f = executor.submit(utils.hmmscan, seqfile, hmm_db)
            fs[f] = (acc, seqfile)

        to_persist = []
        cnt_done = 0
        for f in as_completed(fs):
            acc, seqfile = fs[f]

            try:
                outfile, tabfile = f.result()
            except Exception as exc:
                logger.error(f"{acc}: {exc}")
            else:
                to_persist.append((acc, seqfile, outfile, tabfile))
            finally:
                cnt_done += 1
                if progress:
                    sys.stderr.write(f"progress: "
                                     f"{cnt_done/len(fs)*100:>3.0f}%\r")

    logger.info("persisting data")
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    t1 = orautils.TablePopulator(con, "UPDATE INTERPRO.CLAN_MEMBER "
                                      "SET SEQ = :1 WHERE METHOD_AC = :2")
    t2 = orautils.TablePopulator(con, "INSERT INTO INTERPRO.CLAN_MEMBER_ALN "
                                      "VALUES (:1, :2, :3, :4)",
                                 buffer_size=1000)
    t2.cur.setinputsizes(25, 25, cx_Oracle.NATIVE_FLOAT, cx_Oracle.CLOB)

    cnt_done = 0
    for acc, seqfile, outfile, tabfile in to_persist:
        t1.update((utils.load_sequence(seqfile), acc))

        for t in utils.load_hmmscan_results(outfile, tabfile):
            if acc == t["accession"]:
                continue

            domains = []
            for dom in t["domains"]:
                domains.append({
                    "query": dom["sequences"]["query"],
                    "target": dom["sequences"]["target"],
                    "ievalue": dom["ievalue"],
                    "start": dom["coordinates"]["ali"]["start"],
                    "end": dom["coordinates"]["ali"]["end"],
                })

            t2.insert((
                acc,
                t["accession"],
                t["evalue"],
                json.dumps(domains)
            ))

        cnt_done += 1
        if progress:
            sys.stderr.write(f"progress: "
                             f"{cnt_done/len(to_persist)*100:>3.0f}%\r")

    t1.close()
    t2.close()
    con.commit()
    con.close()

    shutil.rmtree(workdir)
    logger.info("done")


def align_pssm(user: str, dsn: str, sequences_file: str, threads: int=1,
               tmpdir: Optional[str]=None, progress: bool=False):
    dbcode = cdd.DBCODE

    if tmpdir:
        os.makedirs(tmpdir, exist_ok=True)

    workdir = mkdtemp(dir=tmpdir)

    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    cur.execute(
        """
        SELECT METHOD_AC FROM INTERPRO.CLAN_MEMBER
        WHERE CLAN_AC IN (SELECT CLAN_AC FROM INTERPRO.CLAN WHERE DBCODE = :1)
        """, (dbcode,)
    )
    members = {row[0] for row in cur}
    cur.execute(
        """
        DELETE 
        FROM INTERPRO.CLAN_MEMBER_ALN 
        WHERE QUERY_AC IN (
          SELECT METHOD_AC 
          FROM INTERPRO.METHOD 
          WHERE DBCODE = :1
        )
        """, (dbcode,)
    )
    con.commit()
    cur.close()
    con.close()

    logger.info("parsing representative sequences")
    fd, files_list = mkstemp(dir=workdir)
    os.close(fd)

    seqfiles = {}
    id2acc = {}
    with open(files_list, "wt") as fh:
        for identifier, accession, seq in utils.iter_fasta(sequences_file):
            if accession not in members or identifier in id2acc:
                continue

            subdir = os.path.join(workdir, accession[:5])
            try:
                os.mkdir(subdir)
            except FileExistsError:
                pass

            seqfile = os.path.join(subdir, accession) + ".fa"
            with open(seqfile, "wt") as fh2:
                fh2.write(seq)

            fh.write(f"{seqfile}\n")
            id2acc[identifier] = accession
            seqfiles[identifier] = seqfile

    logger.info("building profile database")
    fd, profile_database = mkstemp(dir=workdir)
    os.close(fd)
    os.remove(profile_database)
    utils.mk_compass_db(files_list, profile_database)

    logger.info("querying sequences")
    with ThreadPoolExecutor(max_workers=threads) as executor:
        fs = {}
        for identifier, seqfile in seqfiles.items():
            f = executor.submit(utils.compass_vs_db, seqfile, profile_database)
            fs[f] = identifier

        to_persist = []
        cnt_done = 0
        for f in as_completed(fs):
            identifier = fs[f]

            try:
                outfile = f.result()
            except Exception as exc:
                logger.error(f"{identifier}: {exc}")
            else:
                to_persist.append((identifier, outfile))
            finally:
                cnt_done += 1
                if progress:
                    sys.stderr.write(f"progress: "
                                     f"{cnt_done/len(fs)*100:>3.0f}%\r")

    logger.info("persisting data")
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    t1 = orautils.TablePopulator(con, "UPDATE INTERPRO.CLAN_MEMBER "
                                      "SET SEQ = :1 WHERE METHOD_AC = :2")
    t2 = orautils.TablePopulator(con, "INSERT INTO INTERPRO.CLAN_MEMBER_ALN "
                                      "VALUES (:1, :2, :3, :4)")
    t2.cur.setinputsizes(25, 25, cx_Oracle.NATIVE_FLOAT, cx_Oracle.CLOB)

    cnt_done = 0
    for identifier, outfile in to_persist:
        query_acc = id2acc[identifier]
        seqfile = seqfiles[identifier]
        t1.update((utils.load_sequence(seqfile), query_acc))

        for t in utils.load_compass_results(outfile):
            target_acc = id2acc[t["id"]]

            if target_acc == query_acc:
                continue

            t2.insert((
                query_acc,
                target_acc,
                t["evalue"],
                json.dumps([{
                    "query": t["sequences"]["query"],
                    "target": t["sequences"]["target"],
                    "ievalue": None,
                    "start": t["start"],
                    "end": t["end"],
                }])
            ))

        cnt_done += 1
        if progress:
            sys.stderr.write(f"progress: "
                             f"{cnt_done/len(to_persist)*100:>3.0f}%\r")

    t1.close()
    t2.close()
    con.commit()
    con.close()

    shutil.rmtree(workdir)
    logger.info("done")


def main():
    parser = argparse.ArgumentParser(description="Clans/sets update")
    subparsers = parser.add_subparsers(dest="command")

    subparser = subparsers.add_parser("update", help="update clans")
    subparser.add_argument("-c", "--config",
                           help="config file", required=True)
    subparser.add_argument("--databases", nargs="+", required=True,
                           choices=["cdd", "panther", "pfam", "pirsf"],
                           help="member databases to update")
    subparser.add_argument("--pfam-mysql",
                           help="Pfam MySQL connection string "
                                "(format: user/password@host:port/database)")

    subparser = subparsers.add_parser("align", help="perform profile-profile "
                                                    "alignments")
    subparser.add_argument("-c", "--config",
                           help="config file", required=True)
    subparser.add_argument("--database", required=True,
                           choices=["cdd", "panther", "pfam", "pirsf"],
                           help="member databases to update")
    subparser.add_argument("--data", required=True,
                           help="HMM database (PANTHER, Pfam, PIRSF) or "
                                "FASTA file containing representative "
                                "sequences (CDD)")
    subparser.add_argument("-T", dest="temporary", help="temporary directory")
    subparser.add_argument("-t", dest="threads", help="number of threads",
                           type=int, default=1)
    subparser.add_argument("--progress", action="store_true",
                           help="show progress")

    args = parser.parse_args()
    if args.command is None:
        parser.print_usage()
        parser.exit(2)

    if not os.path.isfile(args.config):
        parser.error(f"No such file or directory: '{args.config}'")

    with open(args.config) as fh:
        config = json.load(fh)

    user = config["database"]["users"]["interpro"]
    dsn = config["database"]["dsn"]

    if args.command == "update":
        if "pfam" in args.databases and not args.pfam_mysql:
            parser.error("--pfam-mysql is required "
                         "when passing 'pfam' to --databases")

        update_clans(user, dsn, list(set(args.databases)), args.pfam_mysql)
    elif not os.path.isfile(args.data):
        parser.error(f"No such file or directory: '{args.data}'")
    elif args.database == "cdd":
        align_pssm(user, dsn, args.data, threads=args.threads,
                   tmpdir=args.temporary, progress=args.progress)
    else:
        align_hmm(user, dsn, args.database, args.data, threads=args.threads,
                  tmpdir=args.temporary, progress=args.progress)
