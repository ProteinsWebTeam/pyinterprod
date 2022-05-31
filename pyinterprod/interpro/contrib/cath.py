import json
import re
from urllib.request import urlopen

from .common import Method


_PREFIX = "G3DSA:"
_TYPE_SUPFAM = 'H'
_TYPE_FUNFAM = 'D'


def parse_superfamilies(filepath: str) -> list[Method]:
    """
    Parse the CathNames.txt file distributed with CATH-Gene3D releases

    :param filepath:
    :return:
    """
    signatures = []
    reg = re.compile(r"^(\d\.\d+\.\d+\.\d+)\s+([a-zA-Z0-9]+)\s+:(.*)$")
    with open(filepath, "rt") as fh:
        for line in fh:
            if line[0] == '#':
                continue

            m = reg.match(line)
            if m is None:
                continue

            supfam, model, name = m.groups()
            accession = f"{_PREFIX}{supfam}"

            m = Method(accession, _TYPE_SUPFAM, description=name)
            signatures.append(m)

    return signatures


def parse_functional_families(file: str) -> list[Method]:
    """
    :param file: TSV file of FunFam names.
                 Can be generated using `get_funfam_names()`.
    :return: A list of FunFam signatures
    """

    signatures = []
    with open(file, "rt") as fh:
        for line in fh:
            accession, name = line.rstrip().split("\t")

            # accession format: 1.10.10.10-FF-000001
            supfam, _, funfam = accession.split("-")

            signatures.append(Method(
                accession=f"{_PREFIX}{supfam}:FF:{funfam}",
                sig_type=_TYPE_FUNFAM,
                name=None if name == "-" else name
            ))

    return signatures


def fetch(url):
    with urlopen(url) as f:
        response = f.read()

    return json.loads(response.decode("utf-8"))


def get_funfam_names(version: str = "v4_3_0") -> dict[str, str]:
    """
    Fetch FunFam names from the CATH REST API
    :param version: CATH version
    :return: dictionary of FunFam ID -> name
    """
    api_url = f"http://www.cathdb.info/version/{version}/api/rest"

    url = f"{api_url}/superfamily/"
    superfamilies = fetch(url)["data"]

    funfams = {}
    for superfam in superfamilies:
        supfam_id = superfam["superfamily_id"]

        url = f"{api_url}/superfamily/{supfam_id}/funfam"
        funfams = fetch(url)["data"]

        for funfam in funfams:
            funfam_id = int(funfam["funfam_number"])
            funfam_name = funfam["name"] or "-"

            model_id = f"{supfam_id}-FF-{funfam_id:06}"

            funfams[model_id] = funfam_name

    return funfams
