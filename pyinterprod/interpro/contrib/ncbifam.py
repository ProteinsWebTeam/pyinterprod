from .common import Method, parse_hmm


_KNOWN_SOURCES = {
    "JCVI",
    "NCBIFAM",
    "NCBI Protein Cluster (PRK)"
}


def get_signatures(filepath: str):
    """
    Parse the NCBIFam HMM file,
    available at https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.LIB

    :param filepath:
    :return:
    """

    models = []
    for acc, name, descr, date in parse_hmm(filepath):
        if descr:
            parts = descr.split(":", 1)

            if parts[0] not in _KNOWN_SOURCES:
                raise ValueError(f"{name}: invalid DESC field {descr}")

            descr = descr[1].strip()

        models.append(Method(acc, "TODO", name, descr, "TODO", date))

    return models
