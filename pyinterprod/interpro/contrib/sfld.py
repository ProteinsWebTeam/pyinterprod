import re

from oracledb import Cursor

from .common import Clan, Method


def get_signatures(cur: Cursor) -> list[Method]:
    """
    We don't expect to ever do an actual SFLD update,
    as the PI retired and the website isn't even available anymore.
    If we need to do a data (matches) updates, we simply get all existing
    signatures.
    """
    method2pubmed = {}
    cur.execute(
        """
        SELECT M2P.METHOD_AC, C.PUBMED_ID
        FROM INTERPRO.METHOD2PUB M2P 
        INNER JOIN INTERPRO.CITATION C on M2P.PUB_ID = C.PUB_ID
        WHERE M2P.METHOD_AC IN (
            SELECT METHOD_AC 
            FROM INTERPRO.METHOD 
            WHERE DBCODE = 'B'
        )
        """
    )
    for method_acc, pubmed_id in cur.fetchall():
        try:
            method2pubmed[method_acc].append(pubmed_id)
        except KeyError:
            method2pubmed[method_acc] = [pubmed_id]

    cur.execute(
        """
        SELECT METHOD_AC, NAME, DESCRIPTION, SIG_TYPE, ABSTRACT, ABSTRACT_LONG
        FROM INTERPRO.METHOD
        WHERE DBCODE = 'B'
        """
    )

    methods = []
    for row in cur:
        m = Method(
            accession=row[0],
            sig_type=row[3],
            name=row[1],
            description=row[2],
            abstract=row[5].read() if row[5] else row[4],
            references=method2pubmed.get(row[0], [])
        )
        methods.append(m)

    return methods
