import json
import xml.etree.ElementTree as xmlET
from urllib import parse, request, error
from concurrent.futures import as_completed, ThreadPoolExecutor

import oracledb

from pyinterprod.utils.oracle import drop_table
from .common import Method, parse_hmm


_KNOWN_SOURCES = {
    "JCVI",
    "NCBIFAM",
    "NCBI Protein Cluster (PRK)"
}

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ESEARCH = f"{EUTILS}/esearch.fcgi"
ESUMMARY = f"{EUTILS}/esummary.fcgi"
NCBI_API = "https://www.ncbi.nlm.nih.gov/genome/annotation_prok/evidence/api/data/"


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
        acc, _ = acc.split('.')

        if descr:
            parts = descr.split(":", 1)

            if parts[0] not in _KNOWN_SOURCES:
                raise ValueError(f"{name}: invalid DESC field {descr}")

            descr = parts[1].strip()

        try:
            obj = info[acc]
        except KeyError:
            abstract = None
            references = []
            _type = "family"
        else:
            abstract = obj["public_comment"]
            references = obj["pubmed"]
            _type = obj["family_type"]

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
    params = {
        "db": "protfam", 
        "term": "hmm", 
        "field": "method", 
        "retstart": 0, 
        "retmax": 10000
    }
    while True:
        r = _fetch_url_xml(ESEARCH, params)
        cnt = 0
        for i in r.findall("./IdList/Id"):
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
        params["id"] = ",".join(list(ids)[i:i+step])
        r = _fetch_url_xml(ESUMMARY, params, post_request=True)
        elements = r.findall("./DocumentSummarySet/DocumentSummary/DispFamilyAcc")
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
    url = f"{NCBI_API}?collection=hmm_info&match=accession_._{accession}"
    response = request.urlopen(url)
    payload = json.loads(response.read().decode("utf-8"))
    if payload["totalCount"] > 1:
        raise Exception(f"{accession}: more than one entry")

    entry = payload["data"][0]
    go_terms = set()
    for term in entry.get("go_terms", "").split(";"):
        if term.strip():
            go_terms.add(term.strip())

    references = set()
    for pmid in entry.get("pubmed", "").split(";"):
        if pmid.strip():
            references.add(int(pmid.strip()))

    return {
        "accession": entry["accession"],
        "family_type": entry.get("family_type"),
        "go_terms": list(go_terms),
        "computed_hmm_name": entry.get("computed_hmm_name"),
        "product_name": entry.get("product_name"),
        "public_comment": entry.get("public_comment"),
        "pubmed": list(references),
        "short_name": entry.get("short_name"),
    }


def _fetch_url_xml(url: str, params: dict,
                   post_request: bool = False) -> xmlET.Element:
    data = parse.urlencode(params)
    if post_request:
        response = request.urlopen(url, data=data.encode("ascii"))
    else:
        response = request.urlopen(f"{url}?{data}")
    return xmlET.fromstring(response.read().decode("utf-8"))


def update_go_terms(uri: str, file_path: str):
    con = oracledb.connect(uri)
    cur = con.cursor()
    cur.execute("SELECT METHOD_AC FROM INTERPRO.METHOD WHERE DBCODE = 'N'")
    signatures = {acc for acc, in cur.fetchall()}

    drop_table(cur, "INTERPRO.NCBIFAM2GO", purge=True)
    cur.execute(
        """
        CREATE TABLE INTERPRO.NCBIFAM2GO
        (
            METHOD_AC VARCHAR2(25) NOT NULL
                CONSTRAINT FK_NCBIFAM2GO
                REFERENCES INTERPRO.METHOD (METHOD_AC) ON DELETE CASCADE,
            GO_ID VARCHAR2(10) NOT NULL,
            CONSTRAINT PK_NCBIFAM2GO
            PRIMARY KEY (METHOD_AC, GO_ID)
        ) NOLOGGING
        """
    )

    sql = """
        INSERT /*+ APPEND */ 
        INTO INTERPRO.NCBIFAM2GO
        VALUES (:1, :2)
    """

    records = []

    with open(file_path, "rt") as fh:
        data = json.load(fh)
        for signature in data:
            accession = signature["accession"]
            if accession in signatures:
                for go_id in set(signature["go_terms"]):
                    records.append((accession, go_id))

    if records:
        cur.executemany(sql, records)
        con.commit()

    cur.close()
    con.close()


def main():
    ids = get_ids()
    accessions = get_accessions(ids)
    info = get_ncbifam_info(accessions)
    print(json.dumps(info, indent=4))


if __name__ == "__main__":
    main()
