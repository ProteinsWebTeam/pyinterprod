import json
import xml.etree.ElementTree as xmlET
from urllib import parse, request, error
from concurrent.futures import as_completed, ThreadPoolExecutor

from .common import Method, parse_hmm


_KNOWN_SOURCES = {
    "JCVI",
    "NCBIFAM",
    "NCBI Protein Cluster (PRK)"
}

EUTILS = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'
ESEARCH = f'{EUTILS}/esearch.fcgi'
ESUMMARY = f'{EUTILS}/esummary.fcgi'
NCBI_API = 'https://www.ncbi.nlm.nih.gov/genome/annotation_prok/evidence/api/data/'
INFO_FILTER_LIST = ['accession', 'public_comment', 'product_name', 'short_name',
                    'go_terms', 'pubmed', 'family_type']


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


def get_ids() -> set[str]:
    ids = set()
    params = {"db": "protfam", "term": "hmm", "field": "method", "retstart": 0, "retmax": 10000}
    while True:
        r = _fetch_url_xml(ESEARCH, params)
        cnt = 0
        for i in r.findall('./IdList/Id'):
            ids.add(i.text)
            cnt += 1
        if cnt == 0:
            break
        params["retstart"] += params["retmax"]
    return ids


def get_accessions(ids: set[str]) -> set[str]:
    ncbi_accessions = set()
    step = 500
    params = {"db": "protfam"}
    for i in range(0, len(ids), step):
        params["id"] = ','.join(list(ids)[i:i+step])
        r = _fetch_url_xml(ESUMMARY, params, post_request=True)
        elements = r.findall('./DocumentSummarySet/DocumentSummary/DispFamilyAcc')
        for a in elements:
            ncbi_accessions.add(a.text.split(".")[0])
    return ncbi_accessions


def get_ncbifam_info(accessions: set) -> list:
    results = []
    with ThreadPoolExecutor(max_workers=8) as executor:
        fs = {executor.submit(_request_ncbi_info, i): i for i in accessions}
        while fs:
            for future in as_completed(fs):
                try:
                    future.result()
                except error.HTTPError as e:
                    if e.code != 429:
                        raise
                    else:
                        data = fs[future]
                        retry = executor.submit(_request_ncbi_info, data)
                        fs[retry] = data
                else:
                    results.append(future.result())
                fs.pop(future)
    return results


def _request_ncbi_info(accession: str) -> dict:
    url = f'{NCBI_API}?collection=hmm_info&match=accession_._{accession}'
    result = request.urlopen(url)
    parsed_result = result.read().decode('utf-8')
    filtered_info = _filter_ncbifam_info(parsed_result)
    return filtered_info


def _filter_ncbifam_info(info: str) -> dict:
    infos = {}
    json_info = json.loads(info)
    if json_info['totalCount'] > 1:
        raise Exception(f"Returned more than one info version for the same accession: {json_info['data'][0]['accession']}.")
    for filter_key in INFO_FILTER_LIST:
        try:
            value = json_info['data'][0][filter_key]
            if filter_key in ['go_terms', "pubmed"]:
                value = value.split(';')
                if filter_key == "pubmed":
                    value = list(map(int, value))
        except KeyError:
            value = None
        infos[filter_key] = value
    return infos


def _fetch_url_xml(url: str, params: dict, post_request: bool = False) -> xmlET.Element:
    data = parse.urlencode(params)
    if post_request:
        response = request.urlopen(url, data=data.encode('ascii'))
    else:
        response = request.urlopen(f'{url}?{data}')
    return xmlET.fromstring(response.read().decode('utf-8'))


def main():
    ids = get_ids()
    accessions = get_accessions(ids)
    ncbifam_info = get_ncbifam_info(accessions)
    info_json_parsed = json.dumps(ncbifam_info, indent=4)
    print(info_json_parsed)


if __name__ == "__main__":
    main()
