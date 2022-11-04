from .common import Method, parse_hmm


def parse_models(filepath: str) -> list[Method]:
    """
    Parse the AntiFam HMM file,
    available at https://ftp.ebi.ac.uk/pub/databases/Pfam/AntiFam/current/Antifam.tar.gz

    :param filepath:
    :return:
    """

    models = []
    for acc, name, descr, date in parse_hmm(filepath):
        models.append(Method(acc, None, name, descr, None, date))

    return models
