import json

from .common import Method, parse_hmm


_KNOWN_SOURCES = {
    "JCVI",
    "NCBIFAM",
    "NCBI Protein Cluster (PRK)"
}


def get_signatures(hmm_file: str, info_file: str):
    """
    NCBIFam HMM file: https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.LIB
    """

    info = {}
    with open(info_file, "rt") as fh:
        for e in json.load(fh):
            acc = e["accession"]
            info[acc] = e

    signatures = []
    for acc, name, descr, date in parse_hmm(hmm_file):
        if descr:
            parts = descr.split(":", 1)

            if parts[0] not in _KNOWN_SOURCES:
                raise ValueError(f"{name}: invalid DESC field {descr}")

            descr = descr[1].strip()

        try:
            obj = info[acc]
        except KeyError:
            abstract = None
            references = None
            _type = 'family'
        else:
            abstract = obj["abstract"]
            references = obj["pubmed"]
            _type = obj["type"]

        if _type == "repeat":
            _type = "R"
        elif _type == "domain":
            _type = "D"
        else:
            # Defaults to family
            _type = "F"

        signatures.append(Method(acc, _type, name, descr, abstract, date,
                                 references))

    return signatures
