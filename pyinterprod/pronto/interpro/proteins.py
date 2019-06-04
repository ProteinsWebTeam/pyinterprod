# -*- coding: utf-8 -*-

import hashlib
from multiprocessing import Queue
from typing import List, Optional

import cx_Oracle

from .utils import Kvdb, MatchComparator, Organizer, TaxonomyComparator
from ... import orautils

MAX_GAP = 20        # at least 20 residues between positions


def consume_proteins(user: str, dsn: str, kvdb: Kvdb, task_queue: Queue,
                     done_queue: Queue, tmpdir: Optional[str]=None,
                     bucket_size: int=100):
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    cur.execute(
        """
        SELECT METHOD_AC, DBCODE, SIG_TYPE
        FROM {}.METHOD
        ORDER BY METHOD_AC
        """.format(owner)
    )
    keys = []
    signatures = {}
    for i, (acc, dbcode, _type) in enumerate(cur):
        signatures[acc] = (dbcode, _type)
        if not i % bucket_size:
            keys.append(acc)

    proteins2go = {}
    cur.execute("SELECT PROTEIN_AC, GO_ID FROM {}.PROTEIN2GO".format(owner))
    for acc, go_id in cur:
        if acc in proteins2go:
            proteins2go[acc].add(go_id)
        else:
            proteins2go[acc] = {go_id}
    cur.close()

    kvdb.open()
    names = Organizer(keys, dir=tmpdir)
    taxa = Organizer(keys, dir=tmpdir)
    terms = Organizer(keys, dir=tmpdir)
    m_comparator = MatchComparator(tmpdir)
    t_comparator = TaxonomyComparator(tmpdir)
    table = orautils.TablePopulator(
        con=con,
        query="INSERT /*+ APPEND */ "
              "INTO {}.METHOD2PROTEIN "
              "VALUES (:1, :2, :3, :4, :5, :6, :7)".format(owner),
        autocommit=True
    )
    for chunk in iter(task_queue.get, None):
        for acc, dbcode, length, taxid, descid, matches in chunk:
            md5 = hash_protein(matches)
            protein_terms = proteins2go.get(acc, [])
            signatures = m_comparator.update(matches)

            try:
                tax_ranks = kvdb[taxid]
            except KeyError:
                tax_ranks = {}
            else:
                t_comparator.update(signatures, tax_ranks.keys())

            for signature_acc in signatures:
                # UniProt descriptions
                names.add(signature_acc, (descid, dbcode))
                # Taxonomic origins
                for rank, rank_tax_id in tax_ranks.items():
                    taxa.add(signature_acc, (rank, rank_tax_id))
                # GO terms
                for go_id in protein_terms:
                    terms.add(signature_acc, go_id)

                table.insert((signature_acc, acc, dbcode, md5, length,
                              taxid, descid))

        names.flush()
        taxa.flush()
        terms.flush()

    table.close()
    con.close()
    kvdb.close()
    m_comparator.sync()
    t_comparator.sync()
    sizes = sum((names.merge(), taxa.merge(), terms.merge(), m_comparator.size, t_comparator.size))
    done_queue.put((names, taxa, terms, m_comparator, t_comparator, sizes))


def hash_protein(matches: List[tuple]) -> str:
    # flatten matches
    locations = []
    for acc, pos_start, pos_end in matches:
        locations.append((pos_start, acc))
        locations.append((pos_end, acc))

    """
    Evaluate the protein's match structure,
        i.e. how signatures match the proteins

    -----------------------------   Protein
     <    >                         Signature 1
       <    >                       Signature 2
                  < >               Signature 3

    Flattened:
    -----------------------------   Protein
     < <  > >     < >
     1 2  1 2     3 3

    Structure, with '-' representing a "gap"
        (more than N bp between two positions):
    1212-33
    """

    # Sort locations by position
    locations.sort()

    """
    Do not set the offset to 0, but to the first position:
    if two proteins have the same structure,
    but the first position of one protein is > max_gap
    while the first position of the other protein is <= max_gap,
    a gap will be used for the first protein and not for the other,
    which will results in two different structures
    """
    offset = locations[0][0]

    # overall match structure
    structure = []
    # close signatures (less than max_gap between two positions)
    signatures = []

    for pos, acc in locations:
        if pos > offset + MAX_GAP:
            for _acc in signatures:
                structure.append(_acc)

            signatures = []
            structure.append('')  # add a gap

        offset = pos
        signatures.append(acc)

    for _acc in signatures:
        structure.append(_acc)

    return hashlib.md5(
        '/'.join(structure).encode("utf-8")
    ).hexdigest()
