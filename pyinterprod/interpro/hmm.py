import gzip
import re
from typing import Optional

import cx_Oracle

from pyinterprod import logger
from .database import Database


"""
CREATE TABLE INTERPRO.METHOD_HMM
(
    METHOD_AC VARCHAR2(25) NOT NULL,
    MODEL_AC VARCHAR2(255),
    HMM BLOB NOT NULL,
    CONSTRAINT UQ_METHOD_HMM
      UNIQUE (METHOD_AC, MODEL_AC),
    CONSTRAINT FK_METHOD_HMM$METHOD_AC
      FOREIGN KEY (METHOD_AC)
      REFERENCES INTERPRO.METHOD (METHOD_AC)
      ON DELETE CASCADE
)
"""


class _Mapper:
    def __init__(self, database: str, mapfile: Optional[str]):
        self.model2fam = {}
        if mapfile:
            with open(mapfile, "rt") as fh:
                for line in fh:
                    model, fam, _ = line.rstrip().split('\t')
                    self.model2fam[model] = fam

        if database.lower() in ("antifam", "pirsf", "sfld"):
            self.map = self.basic
        elif database.lower() == "ncbifam":
            self.map = self.ncbifam
        elif database.lower() == "cath-gene3d":
            self.map = self.cathgene3d
        elif database.lower() == "panther":
            self.map = self.panther
        elif database.lower() == "pfam":
            self.map = self.pfam
        elif database.lower() == "superfamily":
            self.map = self.superfamily
        else:
            raise ValueError(f"unsupported member database '{database}'")

        self.prog_cathgene3d = re.compile(r"([a-zA-Z0-9]+)-i\d")
        self.prog_panther = re.compile(r"(PTHR\d+)\.(SF\d+)?")

    def __call__(self, *args):
        return self.map(*args)

    @staticmethod
    def basic(acc, name):
        return acc, None

    def cathgene3d(self, acc, name):
        try:
            model = self.prog_cathgene3d.match(name).group(1)
            fam = self.model2fam[model]
        except AttributeError:
            return None, None
        else:
            return "G3DSA:" + fam, model

    def panther(self, acc, name):
        accession, prefix = self.prog_panther.match(name).groups()
        if prefix:
            accession += ':' + prefix

        return accession, None

    @staticmethod
    def pfam(acc, name):
        return re.match(r"PF\d+", acc).group(), None

    @staticmethod
    def ncbifam(acc, name):
        return acc.split(".")[0], None

    @staticmethod
    def superfamily(acc, name):
        return "SSF" + acc, name


def update(url: str, database: Database, hmmfile: str, mapfile: Optional[str]):
    con = cx_Oracle.connect(url)
    cur = con.cursor()

    logger.info(f"{database.name}: deleting HMMs")
    cur.execute(
        """
        DELETE FROM INTERPRO.METHOD_HMM
        WHERE METHOD_AC IN (
            SELECT METHOD_AC
            FROM INTERPRO.METHOD
            WHERE DBCODE = :1
        )
        """,
        [database.identifier]
    )

    logger.info(f"{database.name}: inserting HMMs")
    with open(hmmfile, "rt") as fh:
        mapper = _Mapper(database.name, mapfile)
        prog_acc = re.compile(r"^ACC\s+(.+)$", re.M)
        prog_name = re.compile(r"^NAME\s+(.+)$", re.M)

        buffer = ""
        cnt = 0
        for line in fh:
            buffer += line
            if line[:2] == "//":
                cnt += 1

                try:
                    # Optional
                    model_acc = prog_acc.search(buffer).group(1)
                except AttributeError:
                    model_acc = None

                # Mandatory
                model_name = prog_name.search(buffer).group(1)

                signature_acc, model_acc = mapper(model_acc, model_name)
                if signature_acc is None:
                    continue

                blob = gzip.compress(buffer.encode("utf-8"))
                cur.execute(
                    """
                    INSERT INTO INTERPRO.METHOD_HMM
                    VALUES (:1, :2, :3)
                    """, (signature_acc, model_acc, blob)
                )

                buffer = ""

    con.commit()
    cur.close()
    con.close()

    logger.info(f"{database.name}: complete")
