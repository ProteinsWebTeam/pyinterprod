# -*- coding: utf-8 -*-

import json
import os
import re
import shutil
import subprocess as sp
from concurrent import futures
from tempfile import mkdtemp, mkstemp
from typing import List, Tuple

import cx_Oracle

from pyinterprod import logger
from pyinterprod.utils import oracle, Table
from . import contrib

HMM_SUFFIX = ".hmm"
SEQ_SUFFIX = ".fa"
DOM_SUFFIX = ".tab"
OUT_SUFFIX = ".out"
DATABASES = {
    "cdd": 'J',
    "panther": 'V',
    "pfam": 'H',
    "pirsf": 'U'
}


def calc_dir_size(dirpath: str) -> int:
    size = 0
    for root, dirs, files in os.walk(dirpath):
        for f in files:
            size += os.path.getsize(os.path.join(root, f))

    return size


def create_tables(url: str):
    con = cx_Oracle.connect(url)
    cur = con.cursor()

    for table in ("CLAN_MEMBER_ALN", "CLAN_MATCH", "CLAN_MEMBER", "CLAN"):
        oracle.drop_table(cur, table, purge=True)

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
            MEMBER_AC VARCHAR2(25) NOT NULL,
            LEN NUMBER NOT NULL,
            SCORE FLOAT NOT NULL,
            CONSTRAINT PK_CLAN_MEMBER
              PRIMARY KEY (CLAN_AC, MEMBER_AC),
            CONSTRAINT UQ_CLAN_MEMBER$MEMBER_AC
              UNIQUE (MEMBER_AC),
            CONSTRAINT FK_CLAN_MEMBER$CLAN_AC
              FOREIGN KEY (CLAN_AC)
              REFERENCES INTERPRO.CLAN (CLAN_AC)
              ON DELETE CASCADE,
            CONSTRAINT FK_CLAN_MEMBER$MEMBER_AC
              FOREIGN KEY (MEMBER_AC)
              REFERENCES INTERPRO.METHOD (METHOD_AC)
              ON DELETE CASCADE
        )
        """
    )

    cur.execute(
        """
        CREATE TABLE INTERPRO.CLAN_MATCH
        (
            QUERY_AC VARCHAR2(25) NOT NULL,
            TARGET_AC VARCHAR2(25) NOT NULL,
            EVALUE FLOAT NOT NULL,
            DOMAINS CLOB NOT NULL,
            CONSTRAINT PK_CLAN_MATCH
                PRIMARY KEY (QUERY_AC, TARGET_AC),
            CONSTRAINT FK_CLAN_MATCH
              FOREIGN KEY (QUERY_AC)
              REFERENCES INTERPRO.CLAN_MEMBER (MEMBER_AC)
              ON DELETE CASCADE
        )
        """
    )

    cur.close()
    con.close()


def load_sequence(seqfile: str) -> str:
    seq = ""
    with open(seqfile, "rt") as fh:
        next(fh)
        for line in fh:
            seq += line.rstrip()

    return seq


def prepare_insert(con):
    sql = "INSERT INTO INTERPRO.CLAN VALUES (:1,:2,:3,:4)"
    t1 = Table(con, sql)

    sql = "INSERT INTO INTERPRO.CLAN_MEMBER VALUES (:1, :2, :3, :4)"
    t2 = Table(con, sql, depends_on=t1)

    sql = "INSERT INTO INTERPRO.CLAN_MATCH VALUES (:1, :2, :3, :4)"
    t3 = Table(con, sql, depends_on=t2)
    t3.cur.setinputsizes(25, 25, cx_Oracle.DB_TYPE_BINARY_DOUBLE,
                         cx_Oracle.DB_TYPE_CLOB)

    return t1, t2, t3


def update_hmm_clans(url: str, dbkey: str, hmmdb: str, **kwargs):
    clan_source = kwargs.get("source")
    threads = kwargs.get("threads")
    tmpdir = kwargs.get("tmpdir")

    dbcode = DATABASES[dbkey]

    logger.info("deleting old clans")
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.execute("DELETE FROM INTERPRO.CLAN WHERE DBCODE = :1", (dbcode,))
    con.commit()
    cur.close()
    con.close()

    logger.info("loading new clans")
    if dbkey == "panther":
        clans = contrib.panther.get_clans(url)

        def getsubdir(x): return x[:7]
    elif dbkey == "pfam":
        clans = contrib.pfam.get_clans(clan_source)

        def getsubdir(x): return x[:5]
    elif dbkey == "pirsf":
        clans = contrib.pirsf.get_clans(clan_source)

        def getsubdir(x): return x[:8]
    else:
        raise NotImplementedError()

    clans_to_insert = {}
    mem2clan = {}
    for c in clans:
        clans_to_insert[c.accession] = c
        for m in c.members:
            mem2clan[m["accession"]] = (c.accession, m["score"])

    workdir = mkdtemp(dir=tmpdir)
    num_duplicates = 0
    with futures.ThreadPoolExecutor(max_workers=threads) as executor:
        logger.info("emitting consensus sequences")
        fs = {}
        models = set()
        for model_acc, hmm in iter_models(hmmdb):
            if model_acc not in mem2clan:
                # Ignore models not belonging to a clan
                continue
            elif model_acc in models:
                num_duplicates += 1
                continue

            subdir = os.path.join(workdir, getsubdir(model_acc))
            try:
                os.mkdir(subdir)
            except FileExistsError:
                pass

            prefix = os.path.join(subdir, model_acc)
            hmmfile = prefix + HMM_SUFFIX
            with open(hmmfile, "wt") as fh:
                fh.write(hmm)

            seqfile = prefix + SEQ_SUFFIX
            f = executor.submit(hmmemit, hmmfile, seqfile)
            fs[f] = model_acc
            models.add(model_acc)

        done, not_done = futures.wait(fs)
        if not_done:
            shutil.rmtree(workdir)
            raise RuntimeError(f"{len(not_done)} error(s)")
        elif num_duplicates:
            shutil.rmtree(workdir)
            raise RuntimeError(f"HMM database {hmmdb} contains "
                               f"{num_duplicates} duplicated models.")

        logger.info("searching consensus sequences")
        fs = {}
        for model_acc in models:
            prefix = os.path.join(workdir, getsubdir(model_acc), model_acc)
            seqfile = prefix + SEQ_SUFFIX
            outfile = prefix + OUT_SUFFIX
            domfile = prefix + DOM_SUFFIX
            f = executor.submit(hmmscan, hmmdb, seqfile, domfile, outfile)
            fs[f] = model_acc

        con = cx_Oracle.connect(url)
        t1, t2, t3 = prepare_insert(con)

        completed = errors = progress = 0
        for f in futures.as_completed(fs):
            model_acc = fs[f]
            completed += 1

            if not f.result():
                logger.error(f"{model_acc}")
                errors += 1
                continue

            prefix = os.path.join(workdir, getsubdir(model_acc), model_acc)
            outfile = prefix + OUT_SUFFIX
            domfile = prefix + DOM_SUFFIX

            clan_acc, score = mem2clan[model_acc]
            sequence = load_sequence(prefix + SEQ_SUFFIX)

            try:
                clan = clans_to_insert.pop(clan_acc)
            except KeyError:
                # Clan already inserted
                pass
            else:
                t1.insert((
                    clan.accession,
                    dbcode,
                    clan.name,
                    clan.description
                ))

            t2.insert((clan_acc, model_acc, len(sequence), score))

            for target in load_hmmscan_results(outfile, domfile):
                if target["accession"] == model_acc:
                    continue

                domains = []
                for dom in target["domains"]:
                    domains.append((
                        dom["coordinates"]["ali"]["start"],
                        dom["coordinates"]["ali"]["end"]
                    ))

                t3.insert((
                    model_acc,
                    target["accession"],
                    target["evalue"],
                    json.dumps(domains)
                ))

            pc = completed * 100 // len(fs)
            if pc > progress:
                progress = pc
                logger.debug(f"{progress:>10}%")

        for t in (t1, t2, t3):
            t.close()

        con.commit()
        con.close()

        size = calc_dir_size(workdir)
        shutil.rmtree(workdir)
        if errors:
            raise RuntimeError(f"{errors} error(s)")

    logger.info(f"complete (disk usage: {size / 1024 ** 2:,.0f} MB)")


def iter_models(hmmdb: str):
    with open(hmmdb, "rt") as fh:
        reg_acc = re.compile(r"ACC\s+(\w+)", flags=re.M)
        reg_name = re.compile(r"^NAME\s+(PTHR\d+)\.(SF\d+)?", flags=re.M)
        hmm = ""
        for line in fh:
            hmm += line

            if line[:2] == "//":
                m = reg_acc.search(hmm)
                if m:
                    accession = m.group(1)
                else:
                    # PANTHER: accessions in the NAME field
                    m = reg_name.search(hmm)
                    accession, prefix = m.groups()
                    if prefix is not None:
                        accession += ':' + prefix

                yield accession, hmm
                hmm = ""


def hmmemit(hmmdb: str, seqfile: str):
    sp.run(args=["hmmemit", "-c", "-o", seqfile, hmmdb],
           stderr=sp.DEVNULL, stdout=sp.DEVNULL, check=True)


def hmmscan(hmmdb: str, seqfile: str, domfile: str, outfile: str) -> bool:
    args = ["hmmscan", "-o", outfile, "--domtblout", domfile, hmmdb, seqfile]
    process = sp.run(args=args, stderr=sp.DEVNULL, stdout=sp.DEVNULL)

    if process.returncode == 0:
        return True

    for f in (domfile, outfile):
        try:
            os.remove(f)
        except FileNotFoundError:
            pass

    return False


def load_hmmscan_results(outfile: str, tabfile: str) -> List[dict]:
    alignments = load_domain_alignments(outfile)
    targets = {}

    with open(tabfile, "rt") as fh:
        i = 0
        for line in fh:
            if line[0] == "#":
                continue

            cols = re.split(r"\s+", line.rstrip(), maxsplit=22)

            name = cols[0]

            # Pfam entries end with a mark followed by a number
            acc = cols[1].split(".")[0]

            if acc == "-":
                # Panther accessions are under the `target_name` column
                acc = name

            if acc in targets:
                t = targets[acc]
            else:
                t = targets[acc] = {
                    "name": name,
                    "accession": acc,
                    "tlen": int(cols[2]),
                    "qlen": int(cols[5]),

                    # full sequence
                    "evalue": float(cols[6]),
                    "evaluestr": cols[6],
                    "score": float(cols[7]),
                    "bias": float(cols[8]),

                    "domains": []
                }

            t["domains"].append({
                # this domain

                # conditional E-value
                "cevalue": float(cols[11]),
                "cevaluestr": cols[11],
                # independent E-value
                "ievalue": float(cols[12]),
                "ievaluestr": cols[12],
                "score": float(cols[13]),
                "bias": float(cols[14]),

                "coordinates": {
                    # target (as we scan an HMM DB)
                    "hmm": {
                        "start": int(cols[15]),
                        "end": int(cols[16])
                    },
                    # query
                    "ali": {
                        "start": int(cols[17]),
                        "end": int(cols[18])
                    },
                    "env": {
                        "start": int(cols[19]),
                        "end": int(cols[20])
                    },
                },
                "sequences": alignments[i]
            })
            i += 1

    return list(targets.values())


def load_domain_alignments(file: str) -> List[Tuple[str, str]]:
    """
    Parse the output file of hmmscan and load domain alignments.

    Example of alignments:
  == domain 1  score: 25.3 bits;  conditional E-value: 5.2e-09
                   Cytochrome_c4  11 llalaalal.alaaaadaeagaaklaea......gaaavkaCaaCHGadGnsaaaaayPrLAgqsaaYlakqLkdfrsg 82
                                     l++l+a+++ ++ a+++ e++a+k++ea       +   ++C +CHG+d ++a+    P+L    ++Y +++++++ ++
  Cytochrome_Bsub_c550-consensus  18 LVVLLAVNGgSKDAEEEKEEEAEKSEEAeaeaegEEIFKQKCISCHGKDLEGAVG---PNLEKVGSKYSEEEIAKIIEN 93
                                     34444444442223333333333333336666856777899***********766...***************999887 PP

                   Cytochrome_c4  83 errknpMaplakaLsdqdiedlaaYfaaq 111
                                        k +M   a+  sd++ +++a+++a++
  Cytochrome_Bsub_c550-consensus  94 G--KGAM--PAAIVSDDEAKAVAKWLAEK 118
                                     3..3344..46678999999999999986 PP

    Since the sequence is too long, the domain is represented with two "blocks".
    The "== domain" line might be followed by a consensus structure annotation line (not the case here).
    Each block has four lines:
        1. consensus of the target profile
        2. matches between the query sequence and target profile (**can be empty**)
        3. query sequence
        4. posterior probability of each aligned residue

    :param file: hmmscan output file
    :return: a list of alignments, represented by a tuple of two sequences (query, target)
    """
    alignments = []
    query_seq = target_seq = ""
    with open(file, "rt") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">> "):
                # New model
                # target_name = line[3:]
                if query_seq:
                    alignments.append((query_seq, target_seq))
                    query_seq = target_seq = ""
            elif line.startswith("== domain"):
                # New domain
                if query_seq:
                    alignments.append((query_seq, target_seq))
                    query_seq = target_seq = ""

                line = next(fh).strip()
                block = []
                while line or len(block) < 4:
                    block.append(line)
                    line = next(fh).strip()

                del block[:-4]
                target_seq += block[0].split()[2]
                query_seq += block[2].split()[2]
            elif line == "Internal pipeline statistics summary:":
                alignments.append((query_seq, target_seq))
                query_seq = target_seq = ""
            elif query_seq:
                # New block of domain
                block = []
                while line or len(block) < 4:
                    block.append(line)
                    line = next(fh).strip()

                del block[:-4]
                target_seq += block[0].split()[2]
                query_seq += block[2].split()[2]

    return alignments


def update_cdd_clans(url: str, cddmasters: str, cddid: str,
                     fam2supfam: str, **kwargs):
    threads = kwargs.get("threads")
    tmpdir = kwargs.get("tmpdir")

    dbcode = DATABASES["cdd"]

    logger.info("deleting old clans")
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.execute("DELETE FROM INTERPRO.CLAN WHERE DBCODE = :1", (dbcode,))
    con.commit()
    cur.close()
    con.close()

    clans = contrib.cdd.get_clans(cddid, fam2supfam)
    clans_to_insert = {}
    mem2clan = {}
    for c in clans:
        clans_to_insert[c.accession] = c
        for m in c.members:
            mem2clan[m["accession"]] = (c.accession, m["score"])

    logger.info("parsing representative sequences")
    workdir = mkdtemp(dir=tmpdir)
    fd, files_list = mkstemp(dir=workdir)

    seqfiles = {}
    with open(fd, "wt") as fh:
        for model_acc, sequence in iter_sequences(cddmasters):
            if model_acc not in mem2clan or model_acc in seqfiles:
                continue

            subdir = os.path.join(workdir, model_acc[:5])
            try:
                os.mkdir(subdir)
            except FileExistsError:
                pass

            prefix = os.path.join(subdir, model_acc)
            seqfile = prefix + SEQ_SUFFIX
            with open(seqfile, "wt") as fh2:
                fh2.write(sequence)

            fh.write(f"{seqfile}\n")
            seqfiles[model_acc] = prefix

    logger.info("building profile database")
    fd, database = mkstemp(dir=workdir)
    os.close(fd)
    os.remove(database)
    sp.run(["mk_compass_db", "-i", files_list, "-o", database],
           stderr=sp.DEVNULL, stdout=sp.DEVNULL, check=True)

    with futures.ThreadPoolExecutor(max_workers=threads) as executor:
        logger.info("querying sequences")
        fs = {}
        for model_acc, prefix in seqfiles.items():
            seqfile = prefix + SEQ_SUFFIX
            outfile = prefix + OUT_SUFFIX
            f = executor.submit(compass_vs_db, seqfile, database, outfile)
            fs[f] = prefix

        # con = cx_Oracle.connect(url)
        # t1, t2, t3 = prepare_insert(con)

        completed = errors = progress = 0
        for f in futures.as_completed(fs):
            prefix = fs[f]
            completed += 1

            if not f.result():
                logger.error(f"{model_acc}")
                errors += 1
                continue

            seqfile = prefix + SEQ_SUFFIX
            outfile = prefix + OUT_SUFFIX

            clan_acc, score = mem2clan[model_acc]

            pc = completed * 100 // len(fs)
            if pc > progress:
                progress = pc
                logger.debug(f"{progress:>10}%")

        # for t in (t1, t2, t3):
        #     t.close()
        #
        # con.commit()
        # con.close()

        size = calc_dir_size(workdir)
        shutil.rmtree(workdir)
        if errors:
            raise RuntimeError(f"{errors} error(s)")

    logger.info(f"complete (disk usage: {size / 1024 ** 2:,.0f} MB)")


def iter_sequences(seqfile: str):
    with open(seqfile, "rt") as fh:
        buffer = ""
        accession = identifier = None
        for line in fh:
            if line[0] == ">":
                if buffer and identifier:
                    yield accession, buffer

                m = re.match(r">(gnl\|CDD\|\d+)\s+(cd\d+),", line)
                if m:
                    identifier, accession = m.groups()
                else:
                    accession = identifier = None

                buffer = ""

            buffer += line

    if buffer and identifier:
        yield accession, buffer


def compass_vs_db(seqfile: str, database: str, outfile: str):
    args = ["compass_vs_db", "-i", seqfile, "-d", database, "-o", outfile]
    process = sp.run(args=args, stderr=sp.DEVNULL, stdout=sp.DEVNULL)

    if process.returncode == 0:
        return True

    try:
        os.remove(outfile)
    except FileNotFoundError:
        pass

    return False


def load_compass_results(outfile) -> List[dict]:
    # p1 = re.compile(r"length\s*=\s*(\d+)")
    p2 = re.compile(r"Evalue\s*=\s*([\d.e\-]+)")

    targets = {}
    block = 0
    query_id = None
    query_seq = ""
    target_id = None
    target_seq = ""
    length = None
    evalue = None
    evalue_str = None
    pos_start = None

    with open(outfile, "rt") as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith("Subject="):
                """
                Format:
                Subject= cd154/cd15468.fa
                length=413	filtered_length=413	Neff=1.000
                Smith-Waterman score = 254	Evalue = 3.36e-16

                (the path after "Subject=" might be truncated)
                """
                if target_id:
                    targets[target_id] = {
                        "id": target_id,
                        "evalue": evalue,
                        "evaluestr": evalue_str,
                        "length": length,
                        "start": pos_start,
                        "end": pos_start + len(query_seq.replace('=', '')) - 1,
                        "sequences": {
                            "query": query_seq,
                            "target": target_seq
                        }
                    }

                query_id = None
                query_seq = None
                target_id = None
                target_seq = None

                line = next(fh)
                # length = int(p1.match(line).group(1))

                line = next(fh)
                evalue_str = p2.search(line).group(1)
                try:
                    evalue = float(evalue_str)
                except ValueError:
                    evalue = 0

                block = 1
            elif line.startswith("Parameters:"):
                # Footer: end of results
                break
            elif not block:
                continue
            elif line:
                """
                First block:
                gnl|CDD|271233   1      PSFIPGPT==TPKGCTRIPSFSLSDTHWCYTHNVILSGCQDHSKSNQYLSLGVIKTNSDG
                CONSENSUS_1      1      PSFIPGPT==TPKGCTRIPSFSLSDTHWCYTHNVILSGCQDHSKSNQYLSLGVIKTNSDG
                                        P++IP+ T      C+R PSF++S+  + YT+ V  ++CQDH +  +Y+++GVI+ ++ G
                CONSENSUS_2      1      PNLIPADTGLLSGECVRQPSFAISSGIYAYTYLVRKGSCQDHRSLYRYFEVGVIRDDGLG
                gnl|CDD|271230   1      PNLIPADTGLLSGECVRQPSFAISSGIYAYTYLVRKGSCQDHRSLYRYFEVGVIRDDGLG

                (following blocks do not have the start position between the ID and the sequence)
                """
                query = line.split()
                next(fh)
                next(fh)
                next(fh)
                target = next(fh).split()

                if block == 1:
                    query_id = query[0]
                    pos_start = int(query[1])
                    query_seq = query[2]
                    target_id = target[0]
                    target_seq = target[2]
                else:
                    query_seq += query[1]
                    target_seq += target[1]

                block += 1

    targets[target_id] = {
        "id": target_id,
        "evalue": evalue,
        "evaluestr": evalue_str,
        "length": length,
        "start": pos_start,
        "end": pos_start + len(query_seq.replace('=', '')) - 1,
        "sequences": {
            "query": query_seq,
            "target": target_seq
        }
    }

    return list(targets.values())
