import json
import xml.etree.ElementTree as xmlET
from urllib import parse, request
from concurrent.futures import as_completed, ThreadPoolExecutor

from .common import Method, parse_hmm
from pyinterprod import logger


_KNOWN_SOURCES = {
    "JCVI",
    "NCBIFAM",
    "NCBI Protein Cluster (PRK)"
}

EUTILS = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'
ESEARCH = f'{EUTILS}/esearch.fcgi'
ESUMMARY = f'{EUTILS}/esummary.fcgi'
NCBI_API = 'https://www.ncbi.nlm.nih.gov/genome/annotation_prok/evidence/api/data/'
INFO_FILTER_LIST = ['public_comment', 'product_name', 'short_name', 'go_terms', 'pubmed', 'family_type']


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


def get_ids() -> set:
    ids = set()
    params = {"db": "protfam", "term": "hmm", "field": "method", "retmax": 10000}
    r = _fetch_url_xml(ESEARCH, params)
    total_ids = r.find('Count').text
    for retstart in range(params["retmax"], int(total_ids)+params["retmax"], params["retmax"]):
        for i in r.findall('./IdList/Id'):
            ids.add(i.text)
        params["retstart"] = retstart
        r = _fetch_url_xml(ESEARCH, params)
    return ids


def get_all_accessions_list(ids_list: list) -> set:
    ncbi_accessions_list = []
    step = 500
    for i in range(0, len(ids_list), step):
        r = _fetch_url_xml(f"{BASE_ACCESSIONS_URL}&id={','.join(ids_list[i:i+step])}")
        accession = r.findall('./DocumentSummarySet/DocumentSummary/DispFamilyAcc')
        for a in accession:
            ncbi_accessions_list.append(a.text.split(".")[0])
    return set(ncbi_accessions_list)


def get_ncbifam_info(accessions: set) -> dict:
    result_dict = {}
    total_accessions_count = len(accessions)
    total_complete_requests = 0
    ncbi_infos_query = f'{NCBI_API}?collection=hmm_info&match=accession_._'
    with ThreadPoolExecutor(max_workers=8) as executor:
        futures_to_data = {executor.submit(request.urlopen, f"{ncbi_infos_query}{i}"): i for i in accessions}
        while total_complete_requests < total_accessions_count:
            for future in as_completed(futures_to_data):
                if future.exception():
                    _retry(future, futures_to_data, executor)
                else:
                    result = future.result().read()
                    filtered_info = _filter_ncbifam_info(result.decode('utf-8'))
                    result_dict.update(filtered_info)
                    total_complete_requests += 1
                futures_to_data.pop(future)
    info = result_dict
    return info


def _retry(future, futures_to_data, executor):
    data = futures_to_data[future]
    retry = executor.submit(request.urlopen, data)
    futures_to_data[retry] = data
    return data


def _filter_ncbifam_info(info: str) -> dict:
    dict_info = {}
    json_info = json.loads(info)
    for i in range(json_info['totalCount']):
        for filter_key in FILTER_LIST:
            try:
                value = json_info['data'][i][filter_key]
            except KeyError:
                value = None
            dict_info[filter_key+str(i)] = value
    return dict_info


def _fetch_url_xml(url: str, params: dict, post_request: bool = False) -> xmlET.Element:
    data = parse.urlencode(params)
    url = request.Request(url, data=data) if post_request else f'{url}?{data}'
    with request.urlopen(url) as f:
        response = f.read()
    return xmlET.fromstring(response.decode('utf-8'))


if __name__ == "__main__":
    ids = get_ids()
    logger.info(f"returned {len(ids)} ids")
    accessions = get_accessions(ids)
    logger.info(f"returned {len(accessions)} accessions")
    ncbifam_info = get_ncbifam_info(accessions)
    info_json_parsed = json.dumps(ncbifam_info, indent=4)
    print(info_json_parsed)
