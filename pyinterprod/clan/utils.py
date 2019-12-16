# -*- coding: utf-8 -*-

import os
import re
from subprocess import Popen, PIPE, DEVNULL
from tempfile import mkstemp
from typing import Generator, List
from urllib.request import urlopen


def download(url: str) -> str:
    fd, dst = mkstemp()
    os.close(fd)

    with urlopen(url) as res, open(dst, "wb") as fh:
        fh.write(res.read())

    return dst


def load_sequence(seqfile: str) -> str:
    seq = ""
    with open(seqfile, "rt") as fh:
        next(fh)
        for line in fh:
            seq += line.rstrip()

    return seq


def parse_block(fh, line):
    block = []
    while line or len(block) < 4:
        block.append(line)
        line = next(fh).strip()

    i = 0
    for i, line in enumerate(block):
        if len(line.split()) == 4:
            break

    target = block[i].split()[2]
    query = block[i+2].split()[2]

    return target, query


def parse_hmmscan_alignments(filepath):
    domains = []
    target = ""
    query = ""
    # p_dom = re.compile(
    #     r"== domain (\d+)\s+score: ([\d.]+) bits;\s*conditional E-value: ([\d.\-e]+)"
    # )
    n_blank = 0
    with open(filepath, "rt") as fh:
        for line in fh:
            line = line.strip()
            if line:
                n_blank = 0
            else:
                n_blank += 1

            if line.startswith("== domain"):
                # new domain: flush previous one
                if target:
                    domains.append({"target": target, "query": query})
                    target = ""
                    query = ""

                t, q = parse_block(fh, next(fh).strip())
                target += t
                query += q

            elif line.startswith(">>"):
                # new complete sequence: flush previous domain
                if target:
                    domains.append({"target": target, "query": query})
                    target = ""
                    query = ""
            elif target:
                if line:
                    t, q = parse_block(fh, line)
                    target += t
                    query += q
                elif n_blank == 2:
                    domains.append({"target": target, "query": query})
                    target = ""
                    query = ""

    return domains


def load_hmmscan_results(outfile, tabfile) -> List[dict]:
    alignments = parse_hmmscan_alignments(outfile)

    targets = {}
    i = 0

    with open(tabfile, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue

            cols = re.split(r"\s+", line.rstrip(), maxsplit=22)

            # Pfam entries end with a mark followed by a number
            acc = cols[1].split(".")[0]

            if acc == "-":
                # Panther accessions are under the `target_name` column
                acc = cols[0]

            if acc in targets:
                t = targets[acc]
            else:
                t = targets[acc] = {
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


def hmmemit(hmmfile: str, seqfile: str):
    with open(seqfile, "wt") as fh:
        p = Popen(["hmmemit", "-c", hmmfile], stdout=fh, stderr=PIPE)
        p.wait()


def hmmscan(seqfile: str, hmmdb: str):
    tabfile = seqfile + ".tab"
    outfile = seqfile + ".out"

    with open(outfile, "wt") as fh:
        # (?) option for --cut_ga and -E
        p = Popen(["hmmscan", "--domtblout", tabfile, hmmdb, seqfile],
                  stdout=fh, stderr=DEVNULL)
        p.wait()

    return outfile, tabfile


def iter_hmm_database(hmmfile: str) -> Generator[tuple, None, None]:
    entries = set()
    with open(hmmfile, "rt") as fh:
        hmm = ""
        for line in fh:
            hmm += line

            if line[:2] == "//":
                m = re.search("^ACC\s+(\w+)", hmm, re.M)
                if m:
                    accession = m.group(1)
                else:
                    # PANTHER: accessions in the NAME field
                    m = re.search("^NAME\s+(PTHR\d+)\.(SF\d+)?", hmm, re.M)
                    accession, prefix = m.groups()
                    if prefix is not None:
                        accession += ':' + prefix

                if accession not in entries:
                    entries.add(accession)
                    yield accession, hmm

                hmm = ""


def iter_fasta(seqfile: str) -> Generator[tuple, None, None]:
    with open(seqfile, "rt") as fh:
        buffer = ""
        accession = identifier = None
        for line in fh:
            if line[0] == ">":
                if buffer and identifier:
                    yield identifier, accession, buffer

                m = re.match(">(gnl\|CDD\|\d+)\s+(cd\d+),", line)
                if m:
                    identifier, accession = m.groups()
                else:
                    accession = identifier = None

                buffer = ""

            buffer += line

    if buffer and identifier:
        yield identifier, accession, buffer


def mk_compass_db(files_list: str, profile_database: str):
    p = Popen(["mk_compass_db", "-i", files_list, "-o", profile_database],
          stderr=DEVNULL, stdout=DEVNULL)
    p.wait()


def compass_vs_db(seqfile: str, database: str) -> str:
    out_file = seqfile + ".out"
    p = Popen(["compass_vs_db", "-i", seqfile, "-d", database, "-o", out_file],
          stderr=DEVNULL, stdout=DEVNULL)
    p.wait()
    return out_file
