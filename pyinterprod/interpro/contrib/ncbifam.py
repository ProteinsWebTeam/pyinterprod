import re
import sys

import cx_Oracle

from pyinterprod.utils.oracle import drop_table
from .common import Method, parse_hmm


_SOURCES_IN_HMM = {
    "JCVI",
    "NCBIFAM",
    "NCBI Protein Cluster (PRK)"
}
_SOURCES_IN_TSV = {
    "EBI-EMBL",
    "JCVI",
    "NCBIFAM",
    "NCBI Protein Cluster (PRK)"
}

# EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
# ESEARCH = f"{EUTILS}/esearch.fcgi"
# ESUMMARY = f"{EUTILS}/esummary.fcgi"
# NCBI_API = "https://www.ncbi.nlm.nih.gov/genome/annotation_prok/evidence/api/data/"


def parse_tsv(tsvfile: str):
    with open(tsvfile, "rt") as fh:
        for line in fh:
            if line[0] != "#":
                yield line.split("\t")


def get_signatures(tsvfile: str, txtfile: str) -> list[Method]:
    """Parse NCBIfam TSV file
    https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.tsv

    :param tsvfile: Path to TSV file
    :param txtfile: Path to TXT file of NCBIFAMs families to import
        (all TIGRFAMs are imported)
    :return:
    """
    ncbifams_to_import = set()
    with open(txtfile, "rt") as fh:
        for accession in map(str.rstrip, fh):
            ncbifams_to_import.add(accession)

    methods = []
    for values in parse_tsv(tsvfile):
        model_acc = values[0]
        signature_acc = model_acc.split(".")[0]

        name = values[2]
        _type = values[6]
        # for_AMRFinder = values[9] == "Y"

        if values[14]:
            references = list(map(int, set(values[14].split(","))))
        else:
            references = []

        source = values[19]
        if source not in _SOURCES_IN_TSV:
            raise ValueError(f"Unknown source in {tsvfile}: {source}")

        hmm_name = values[21]
        comment = values[22].strip() or None

        if source == "JCVI":
            pass
        elif source == "NCBIFAM" and signature_acc in ncbifams_to_import:
            pass
        else:
            continue

        if source != "JCVI" and (
                source != "NCBIFAM" or
                signature_acc not in ncbifams_to_import
        ):
            continue

        if _type == "repeat":
            _type = "R"
        elif _type == "domain":
            _type = "D"
        else:
            # Defaults to family
            _type = "F"

        m = Method(
            accession=signature_acc,
            sig_type=_type,
            name=name,
            description=hmm_name,
            abstract=comment,
            references=references,
            model=model_acc
        )
        methods.append(m)

    return methods


def create_custom_hmm(tsvfile: str, txtfile: str, hmmfile: str):
    models = {m.model: m for m in get_signatures(tsvfile, txtfile)}

    reg_acc = re.compile(r"^ACC\s+(\w+\.\d+)$", flags=re.M)
    reg_desc = re.compile(r"^(DESC\s+)(.+)$", flags=re.M)

    with open(hmmfile, "rt") as fh:
        buffer = ""
        for line in fh:
            buffer += line

            if line[:2] == "//":
                model_acc = reg_acc.search(buffer).group(1)
                if model_acc in models:
                    # desc = models[model_acc].description
                    # buffer = reg_desc.sub(replace_desc(desc), buffer)
                    print(buffer)

                buffer = ""

        if buffer:
            model_acc = reg_acc.search(buffer).group(1)
            if model_acc in models:
                # desc = models[model_acc].description
                # buffer = reg_desc.sub(replace_desc(desc), buffer)
                print(buffer)


def replace_desc(new_desc: str):
    def repl(match: re.Match):
        return f"{match.group(1)}{new_desc}"

    return repl


# def get_ids() -> set[str]:
#     ids = set()
#     params = {
#         "db": "protfam",
#         "term": "hmm",
#         "field": "method",
#         "retstart": 0,
#         "retmax": 10000
#     }
#     while True:
#         r = _fetch_url_xml(ESEARCH, params)
#         cnt = 0
#         for i in r.findall("./IdList/Id"):
#             ids.add(i.text)
#             cnt += 1
#         if cnt == 0:
#             break
#         params["retstart"] += params["retmax"]
#     return ids
#
#
# def get_accessions(ids: set[str]) -> set[str]:
#     ncbi_accessions = set()
#     step = 500
#     params = {"db": "protfam"}
#     for i in range(0, len(ids), step):
#         params["id"] = ",".join(list(ids)[i:i+step])
#         r = _fetch_url_xml(ESUMMARY, params, post_request=True)
#         elements = r.findall("./DocumentSummarySet/DocumentSummary/DispFamilyAcc")
#         for a in elements:
#             ncbi_accessions.add(a.text.split(".")[0])
#     return ncbi_accessions
#
#
# def get_ncbifam_info(accessions: set) -> list:
#     results = []
#     with ThreadPoolExecutor(max_workers=8) as executor:
#         fs = {executor.submit(_request_ncbi_info, i): i for i in accessions}
#         while fs:
#             for future in as_completed(fs):
#                 try:
#                     future.result()
#                 except error.HTTPError as e:
#                     if e.code != 429:
#                         raise
#                     else:
#                         data = fs[future]
#                         retry = executor.submit(_request_ncbi_info, data)
#                         fs[retry] = data
#                 else:
#                     results.append(future.result())
#                 fs.pop(future)
#     return results
#
#
# def _request_ncbi_info(accession: str) -> dict:
#     url = f"{NCBI_API}?collection=hmm_info&match=accession_._{accession}"
#     response = request.urlopen(url)
#     payload = json.loads(response.read().decode("utf-8"))
#     if payload["totalCount"] > 1:
#         raise Exception(f"{accession}: more than one entry")
#
#     entry = payload["data"][0]
#     go_terms = set()
#     for term in entry.get("go_terms", "").split(";"):
#         if term.strip():
#             go_terms.add(term.strip())
#
#     references = set()
#     for pmid in entry.get("pubmed", "").split(";"):
#         if pmid.strip():
#             references.add(int(pmid.strip()))
#
#     return {
#         "accession": entry["accession"],
#         "family_type": entry.get("family_type"),
#         "go_terms": list(go_terms),
#         "computed_hmm_name": entry.get("computed_hmm_name"),
#         "product_name": entry.get("product_name"),
#         "public_comment": entry.get("public_comment"),
#         "pubmed": list(references),
#         "short_name": entry.get("short_name"),
#     }
#
#
# def _fetch_url_xml(url: str, params: dict,
#                    post_request: bool = False) -> xmlET.Element:
#     data = parse.urlencode(params)
#     if post_request:
#         response = request.urlopen(url, data=data.encode("ascii"))
#     else:
#         response = request.urlopen(f"{url}?{data}")
#     return xmlET.fromstring(response.read().decode("utf-8"))


def update_go_terms(uri: str, tsvfile: str):
    con = cx_Oracle.connect(uri)
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

    records = []
    for values in parse_tsv(tsvfile):
        accession = values[0].split(".")[0]

        if accession in signatures and values[13]:
            for go_id in set(values[13].split(",")):
                records.append((accession, go_id))

    for i in range(0, len(records), 1000):
        cur.executemany("INSERT INTO INTERPRO.NCBIFAM2GO VALUES (:1, :2)",
                        records[i:i+1000])

    con.commit()
    cur.close()
    con.close()


def main():
    tsvfile = sys.argv[1]
    txtfile = sys.argv[2]
    hmmfile = sys.argv[3]
    create_custom_hmm(tsvfile, txtfile, hmmfile)


if __name__ == "__main__":
    main()
