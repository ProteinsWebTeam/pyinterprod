# -*- coding: utf-8 -*-

import hashlib
from multiprocessing import Queue
from typing import List, Optional

import cx_Oracle

from . import utils
from ... import orautils

MAX_GAP = 20        # at least 20 residues between positions


def consume_proteins(user: str, dsn: str, task_queue: Queue,
                     done_queue: Queue, tmpdir: Optional[str]=None,
                     bucket_size: int=100):
    owner = user.split('/')[0]
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()
    cur.execute("SELECT METHOD_AC FROM {}.METHOD".format(owner))
    keys = sorted([row[0] for row in cur])
    keys = [keys[i] for i in range(0, len(keys), bucket_size)]

    proteins2go = {}
    cur.execute("SELECT PROTEIN_AC, GO_ID FROM {}.PROTEIN2GO".format(owner))
    for acc, go_id in cur:
        if acc in proteins2go:
            proteins2go[acc].add(go_id)
        else:
            proteins2go[acc] = {go_id}

    cur.execute(
        """
        SELECT DESC_ID
        FROM {0}.DESC_VALUE
        WHERE TEXT LIKE 'Predicted protein%'
          OR TEXT LIKE 'Uncharacterized protein%'
        """.format(owner)
    )
    excluded_descr = {row[0] for row in cur}
    cur.close()

    names = utils.Organizer(keys, dir=tmpdir)
    taxa = utils.Organizer(keys, dir=tmpdir)
    terms = utils.Organizer(keys, dir=tmpdir)
    m_comparator = utils.MatchComparator(dir=tmpdir)
    n_comparator = utils.Comparator(dir=tmpdir)
    ta_comparator = utils.TaxonomyComparator(dir=tmpdir)
    te_comparator = utils.Comparator(dir=tmpdir)
    table = orautils.TablePopulator(
        con=con,
        query="INSERT /*+ APPEND */ "
              "INTO {}.METHOD2PROTEIN "
              "VALUES (:1, :2, :3, :4, :5, :6, :7)".format(owner),
        autocommit=True
    )
    for chunk in iter(task_queue.get, None):
        for acc, dbcode, length, tax_id, ranks, desc_id, matches in chunk:
            md5 = hash_protein(matches)
            protein_terms = proteins2go.get(acc, [])

            # Update comparators
            signatures = m_comparator.update(matches)

            if desc_id not in excluded_descr:
                n_comparator.update(signatures)

            if ranks:
                ta_comparator.update(signatures, ranks.keys())

            te_comparator.update(signatures, incr=len(protein_terms))

            # Update organizers and populate table
            for signature_acc in signatures:
                # UniProt descriptions
                names.add(signature_acc, (desc_id, dbcode))
                # Taxonomic origins
                for rank, rank_tax_id in ranks.items():
                    taxa.add(signature_acc, (rank, rank_tax_id))
                # GO terms
                for go_id in protein_terms:
                    terms.add(signature_acc, go_id)

                table.insert((signature_acc, acc, dbcode, md5, length,
                              tax_id, desc_id))

        names.flush()
        taxa.flush()
        terms.flush()

    table.close()
    con.close()
    size = 0
    comparators = (m_comparator, n_comparator, ta_comparator, te_comparator)
    for c in comparators:
        c.sync()
        size += c.size

    organizers = (names, taxa, terms)
    for o in organizers:
        size += o.merge()

    done_queue.put((*comparators, *organizers, size))


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
