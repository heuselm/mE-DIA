#!/usr/bin/env python3
"""
Enzymatic specificity of peptides, computed from protein sequence.

Same definition as ``run_pipeline._digest_semi`` (fixed 2026-09-19): an
occurrence at [start, end) of a protein is fully specific for an enzyme if both
start and end are boundaries of that enzyme's digest
(``_find_cleavage_sites``; protein termini are boundaries). A peptide is
``fully_specific`` if ANY occurrence is fully specific for ANY enzyme of the
experiment, ``semi_specific`` if it occurs but never fully specifically, ``not_found`` if
it does not occur in the given proteins, and ``not_in_fasta`` if none of its
reported proteins (``candidates``) is in the FASTA, e.g. contaminants.

``met_excision=True`` additionally treats start == 1 after an initiator Met as
a boundary. Off by default: the pipeline does not pass ``--met-excision`` to
DIA-NN, so such peptides are not produced by the fully-specific digest.
"""
from __future__ import annotations

import os
import sys
from typing import Dict, Iterable, List, Optional, Sequence

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from run_pipeline import ENZYME_REGISTRY, _find_cleavage_sites, _parse_fasta  # noqa: E402


def classify(peptides: Iterable[str], proteins: Dict[str, str],
             enzyme_codes: Sequence[str],
             candidates: Optional[Dict[str, List[str]]] = None,
             met_excision: bool = False) -> Dict[str, str]:
    """
    peptides    : stripped sequences
    proteins    : accession -> sequence (``load_proteins``)
    candidates  : optional peptide -> accessions to search (e.g. a report's
                  Protein.Ids split on ';'); default every protein. Peptides
                  not found among their candidates are retried against all.
    """
    enzs = [ENZYME_REGISTRY[c] for c in enzyme_codes]
    site_cache: Dict[tuple, set] = {}

    def sites(acc: str, k: int) -> set:
        key = (acc, k)
        if key not in site_cache:
            site_cache[key] = set(_find_cleavage_sites(proteins[acc], enzs[k]))
        return site_cache[key]

    def one(pep: str, accs) -> str:
        found = False
        for acc in accs:
            seq = proteins.get(acc)
            if seq is None:
                continue
            i = seq.find(pep)
            while i != -1:
                found = True
                j = i + len(pep)
                met = met_excision and i == 1 and seq[0] == "M"
                for k in range(len(enzs)):
                    s = sites(acc, k)
                    if (i in s or met) and j in s:
                        return "fully_specific"
                i = seq.find(pep, i + 1)
        return "semi_specific" if found else "not_found"

    out: Dict[str, str] = {}
    for pep in set(peptides):
        if candidates is None:
            out[pep] = one(pep, list(proteins))
            continue
        accs = candidates.get(pep, [])
        if accs and not any(a in proteins for a in accs):
            # e.g. contaminant entries DIA-NN reports that are not in this FASTA:
            # searching other proteins would find the sequence by accident
            out[pep] = "not_in_fasta"
            continue
        lab = one(pep, accs)
        if lab == "not_found":
            lab = one(pep, list(proteins))
        out[pep] = lab
    return out


def load_proteins(fasta: str) -> Dict[str, str]:
    return _parse_fasta(fasta)
