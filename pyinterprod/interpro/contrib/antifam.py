from .common import Method, parse_hmm


def parse_models(db_sources: dict) -> list[Method]:
    """
    Parse the AntiFam HMM file,
    available at https://ftp.ebi.ac.uk/pub/databases/Pfam/AntiFam/current/Antifam.tar.gz

    :param db_sources:
    :return:
    """
    models = []
    for acc, name, descr, date in parse_hmm(db_sources["sig_source"]):
        models.append(Method(acc, None, name, descr, None, date))

    return models
