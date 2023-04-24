import os
import re
import shutil
import sys
from typing import List

from .common import Method, parse_xml

_TYPE_PATTERNS = 'C'  # conserved sites
_TYPE_PROFILES = 'D'


def parse_patterns(db_sources: dict) -> List[Method]:
    """
    Parse the pattern_model.cgi file provided by PROSITE (XML format)
    http://prosite.expasy.org/cgi-bin/prosite/pattern_model.cgi

    :param db_sources:
    :return:
    """
    return parse_xml(db_sources["sig_source"], _TYPE_PATTERNS)


def parse_profiles(db_sources: dict) -> List[Method]:
    """
    Parse the profile_model.cgi file provided by PROSITE (XML format)
    http://prosite.expasy.org/cgi-bin/prosite/profile_model.cgi

    :param db_sources:
    :return:
    """
    return parse_xml(db_sources["sig_source"], _TYPE_PROFILES)


def parse():
    """
    Parse the prosite.dat file provided by PROSITE at
    ftp://ftp.expasy.org/databases/prosite/prosite.dat
    and create two files:
        - prosite_patterns.dat: patterns
        - prosite_profiles.dat: profiles/matrices
    """
    try:
        file = sys.argv[1]
    except IndexError:
        sys.stderr.write("Usage: python -m pyinterprod.interpro.contrib.prosite "
                         "PROSITE-DAT-FILE\n")
        sys.exit(1)

    outdir = os.path.dirname(os.path.realpath(file))

    patterns_dir = os.path.join(outdir, "prosite_patterns")
    profiles_dir = os.path.join(outdir, "prosite_profiles")

    for dst in [patterns_dir, profiles_dir]:
        try:
            shutil.rmtree(dst)
        except FileNotFoundError:
            pass

        os.mkdir(dst)

    patterns_cnt = profiles_cnt = 0

    with open(file, "rt") as fh:
        buffer = ""
        for line in fh:
            buffer += line

            if line[:2] == "//":
                m = re.search(r"^ID\s+\w+;\s*(\w+)\.$", buffer, re.M)
                if m:
                    profile_type = m.group(1)

                    m = re.search(r"^AC\s+(PS\d+)\;$", buffer, re.M)
                    accession = m.group(1)

                    if profile_type == "MATRIX":
                        dst = os.path.join(profiles_dir, f"{accession}.prf")
                        profiles_cnt += 1
                    elif profile_type == "PATTERN":
                        dst = os.path.join(patterns_dir, f"{accession}.prf")
                        patterns_cnt += 1
                    else:
                        raise ValueError(f"Unknown type {profile_type}")

                    with open(dst, "wt") as fh2:
                        fh2.write(buffer)

                buffer = ""

    print(f"PROSITE Patterns: {patterns_cnt:>10,}")
    print(f"PROSITE Profiles: {profiles_cnt:>10,}")


if __name__ == '__main__':
    parse()
