#!/usr/bin/env python3
"""
proteoform_analysis.py — Detect endogenous N/C-terminal proteoform events.

A proteoform event at residue ``p`` of a protein is called when *at least*
``min_support`` distinct peptides share the SAME free terminus at ``p``:

  * N-terminal event at p : >= min_support peptides start at p AND have a
                            non-specific N-terminus under their observed enzyme.
  * C-terminal event at p : >= min_support peptides end   at p AND have a
                            non-specific C-terminus under their observed enzyme.

Read-through peptides (those that strictly span p but do NOT share its
terminus) are NOT counted as supporting evidence — they only show that the
protein is detected in that region, not that a free terminus exists at p.

Natural protein N-/C-termini (residue 1 and the last residue) are never
called as events.

In-source decay (ISD) filter (on by default, issue #42). Fragmentation in the ion
source, typically N-terminal to proline, yields y-type ions that are complete
peptides with a non-specific N-terminus. DIA-NN can identify them as genuine
shorter peptides; they sit at a different precursor m/z, in another isolation
window, so DIA-NN's no-shared-spectra default (disabled only by --int-removal 0)
cannot prevent them. A supporting precursor is
flagged (``isd_flag``) when an identified longer peptide that shares its
event-distal terminus (same C-terminus for N-terminal events, y-type; same
N-terminus for C-terminal events, b-type, checked for symmetry) has a precursor
whose median apex RT lies within ``isd_rt_window`` minutes and whose median
Precursor.Quantity is higher. A sequence supports an event only through
unflagged precursors; events left below ``min_support`` are dropped.

Outputs:
  - <out>/report.pf_matrix.tsv  precursor x run quantity matrix, one block per
                                proteoform event, with Precursor.Id +
                                Stripped.Sequence + Modified.Sequence and
                                event metadata columns.
  - <out>/report.pf_plot.pdf    multiplicity histogram + per-protein event map.
  - <out>/report.pf_isd.tsv     every ISD-flagged supporting precursor, with its
                                co-eluting parent, including those of dropped
                                events (written when the filter is on).

Usage:
  python proteoform_analysis.py [report_annotated.tsv] [fasta] [out_dir]
                                [--min-support N] [--no-isd-filter]
                                [--isd-rt-window MIN]

Drop-in replacement for proteoform_analysis.R — designed to be moved into the
SemiExtension Python pipeline (e.g. as 12_proteoform_analysis.py) without
modification; the ``run()`` function is the public entry point.
"""
from __future__ import annotations

import argparse
import os
import sys
from typing import Iterable

import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


# ── Enzyme cleavage rules ─────────────────────────────────────────────────────
# sense = 1  => cleavage C-terminal to a residue in ``cut`` (blocked if the
#               residue immediately after is in ``block``).
# sense = 0  => cleavage N-terminal to a residue in ``cut``.
ENZYME_RULES = {
    "No_enzyme":    {"sense": 0, "cut": "",     "block": ""},
    "Trypsin":      {"sense": 1, "cut": "KR",   "block": "P"},
    "Trypsin_P":    {"sense": 1, "cut": "KR",   "block": ""},
    "Lys_C":        {"sense": 1, "cut": "K",    "block": "P"},
    "Lys_N":        {"sense": 0, "cut": "K",    "block": ""},
    "Arg_C":        {"sense": 1, "cut": "R",    "block": "P"},
    "Asp_N":        {"sense": 0, "cut": "D",    "block": ""},
    "CNBr":         {"sense": 1, "cut": "M",    "block": ""},
    "Glu_C":        {"sense": 1, "cut": "DE",   "block": "P"},
    "GluC":         {"sense": 1, "cut": "DE",   "block": "P"},
    "PepsinA":      {"sense": 1, "cut": "FL",   "block": "P"},
    "Chymotrypsin": {"sense": 1, "cut": "FWYL", "block": "P"},
    "Elastase":     {"sense": 1, "cut": "AV",   "block": ""},
    # The pipeline's Enzyme column holds ENZYME_REGISTRY names (run_pipeline.py),
    # which are not the underscore spellings above. A name missing here (and from
    # the registry, see rule_for) counts every terminus as specific, so its
    # peptides could never call an event.
    "Elastase (extended)": {"sense": 1, "cut": "AVSGLI", "block": ""},
    "Trypsin/P":    {"sense": 1, "cut": "KR",   "block": ""},
    "LysC":         {"sense": 1, "cut": "K",    "block": "P"},
    "LysN":         {"sense": 0, "cut": "K",    "block": ""},
    "ArgC":         {"sense": 1, "cut": "R",    "block": "P"},
    "AspN":         {"sense": 0, "cut": "D",    "block": ""},
    # Registry codes, used in some annotated reports
    "T":  {"sense": 1, "cut": "KR",   "block": "P"},
    "T2": {"sense": 1, "cut": "KR",   "block": ""},
    "C":  {"sense": 1, "cut": "FWYL", "block": "P"},
    "G":  {"sense": 1, "cut": "DE",   "block": "P"},
    "E":  {"sense": 1, "cut": "AV",   "block": ""},
    "EX": {"sense": 1, "cut": "AVSGLI", "block": ""},
    "L":  {"sense": 1, "cut": "K",    "block": "P"},
    "LN": {"sense": 0, "cut": "K",    "block": ""},
    "R":  {"sense": 1, "cut": "R",    "block": "P"},
    "A":  {"sense": 0, "cut": "D",    "block": ""},
    "M":  {"sense": 1, "cut": "M",    "block": ""},
    "F":  {"sense": 1, "cut": "FL",   "block": "P"},
}


def rule_for(enzyme):
    """Cleavage rule for an Enzyme-column value: the table above, else the
    pipeline's ENZYME_REGISTRY (by name or code) when run_pipeline is loaded in
    this process, which covers enzymes added with --enzyme-config. None if
    neither knows it."""
    rule = ENZYME_RULES.get(enzyme)
    if rule is not None:
        return rule
    rp = sys.modules.get("run_pipeline")
    reg = getattr(rp, "ENZYME_REGISTRY", None) or {}
    for e in reg.values():
        if enzyme in (e.name, e.code):
            return {"sense": e.sense, "cut": e.cut_res, "block": e.noclv}
    return None


def nterm_specific(seq: str, start: int, enzyme: str) -> bool:
    """True if peptide N-terminus at 1-based ``start`` is enzyme-consistent
    (or is the natural protein N-terminus)."""
    if start <= 1:
        return True
    rule = rule_for(enzyme)
    if rule is None or not rule["cut"]:
        return True
    if rule["sense"] == 1:
        prev_aa = seq[start - 2]
        this_aa = seq[start - 1]
        return prev_aa in rule["cut"] and not (rule["block"] and this_aa in rule["block"])
    return seq[start - 1] in rule["cut"]


def cterm_specific(seq: str, stop: int, enzyme: str, length: int) -> bool:
    if stop >= length:
        return True
    rule = rule_for(enzyme)
    if rule is None or not rule["cut"]:
        return True
    if rule["sense"] == 1:
        this_aa = seq[stop - 1]
        next_aa = seq[stop]
        return this_aa in rule["cut"] and not (rule["block"] and next_aa in rule["block"])
    return seq[stop] in rule["cut"]


# ── I/O helpers ───────────────────────────────────────────────────────────────
def read_fasta(path: str) -> pd.DataFrame:
    accs, seqs = [], []
    cur_acc, cur = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if cur_acc is not None:
                    accs.append(cur_acc); seqs.append("".join(cur).upper())
                header = line[1:]
                if "|" in header:
                    parts = header.split("|")
                    cur_acc = parts[1] if len(parts) >= 2 else parts[0]
                else:
                    cur_acc = header.split()[0]
                cur = []
            else:
                cur.append(line.strip())
        if cur_acc is not None:
            accs.append(cur_acc); seqs.append("".join(cur).upper())
    return pd.DataFrame({"Accession": accs, "Sequence": seqs,
                         "Length": [len(s) for s in seqs]})


_LOG_FH = None


def log(msg: str) -> None:
    print(msg, flush=True)
    if _LOG_FH is not None:
        _LOG_FH.write(msg + "\n")
        _LOG_FH.flush()


# ── Core ──────────────────────────────────────────────────────────────────────
def map_peptides_to_proteins(rep: pd.DataFrame, fa: pd.DataFrame) -> pd.DataFrame:
    """Return rep augmented with Accession, Start, Stop, Length."""
    rep = rep.copy()
    rep["Accession"] = rep["Protein.Ids"].astype(str).str.split(";").str[0]
    pep = rep[["Accession", "Stripped.Sequence"]].drop_duplicates()
    pep = pep.merge(fa, on="Accession", how="inner")
    starts = [s.find(p) for s, p in zip(pep["Sequence"], pep["Stripped.Sequence"])]
    pep["Start"] = [s + 1 if s >= 0 else 0 for s in starts]
    pep = pep[pep["Start"] > 0].copy()
    pep["Stop"] = pep["Start"] + pep["Stripped.Sequence"].str.len() - 1
    out = rep.merge(
        pep[["Accession", "Stripped.Sequence", "Sequence", "Length", "Start", "Stop"]],
        on=["Accession", "Stripped.Sequence"], how="inner",
    )
    return out


def annotate_terminus_specificity(rep: pd.DataFrame) -> pd.DataFrame:
    """Add per-row N_specific, C_specific (under the observed enzyme)."""
    if "Enzyme" not in rep.columns:
        raise ValueError("Report has no 'Enzyme' column — run the annotation "
                         "step (11_annotate_matrices.py) first.")
    rep = rep.copy()
    unknown = sorted({str(e) for e in rep["Enzyme"].dropna().unique()
                      if rule_for(e) is None})
    if unknown:
        n = int(rep["Enzyme"].isin(unknown).sum())
        log(f"  WARNING: no cleavage rule for enzyme(s) {', '.join(unknown)} "
            f"({n:,} rows): their termini count as specific, so their peptides "
            f"cannot support an event.")
    rep["N_specific"] = [nterm_specific(s, st, e)
                        for s, st, e in zip(rep["Sequence"], rep["Start"], rep["Enzyme"])]
    rep["C_specific"] = [cterm_specific(s, sp, e, L)
                        for s, sp, e, L in zip(rep["Sequence"], rep["Stop"],
                                                rep["Enzyme"], rep["Length"])]
    return rep


def detect_events(rep: pd.DataFrame, min_support: int = 2) -> pd.DataFrame:
    """Return a long table of events; one row per supporting peptide.

    A peptide qualifies as evidence for an event at residue p if it SHARES the
    free terminus at p:
      * (N) Start == p AND at least one observation has a non-specific
            N-terminus under its observed enzyme (and p > 1).
      * (C) Stop  == p AND at least one observation has a non-specific
            C-terminus (and p < Length).
    Events with fewer than ``min_support`` distinct supporting peptides are
    discarded.  Read-through peptides spanning p are NOT counted.
    """
    pep = (rep.groupby(["Accession", "Stripped.Sequence", "Start", "Stop", "Length"],
                       as_index=False)
              .agg(N_nonspec=("N_specific", lambda x: (~x).any()),
                   C_nonspec=("C_specific", lambda x: (~x).any())))
    pep["N_evidence"] = pep["N_nonspec"] & (pep["Start"] > 1)
    pep["C_evidence"] = pep["C_nonspec"] & (pep["Stop"]  < pep["Length"])

    rows = []
    n_grp = pep[pep["N_evidence"]].groupby(["Accession", "Start"])
    for (acc, p), sub in n_grp:
        n = sub["Stripped.Sequence"].nunique()
        if n < min_support:
            continue
        for s in sub["Stripped.Sequence"].unique():
            rows.append({"Accession": acc, "aa_index": int(p), "terminus": "N",
                         "Stripped.Sequence": s, "n_support": int(n)})
    c_grp = pep[pep["C_evidence"]].groupby(["Accession", "Stop"])
    for (acc, p), sub in c_grp:
        n = sub["Stripped.Sequence"].nunique()
        if n < min_support:
            continue
        for s in sub["Stripped.Sequence"].unique():
            rows.append({"Accession": acc, "aa_index": int(p), "terminus": "C",
                         "Stripped.Sequence": s, "n_support": int(n)})

    if not rows:
        return pd.DataFrame(columns=["Accession", "aa_index", "terminus",
                                     "Stripped.Sequence", "n_support",
                                     "Proteoform_ID"])
    ev = pd.DataFrame(rows)
    # Proteoform_ID embeds the single-letter residue at aa_index, e.g. P69952_A110N
    acc2seq = dict(zip(rep["Accession"], rep["Sequence"]))
    ev["residue"] = [acc2seq[a][p - 1] for a, p in zip(ev["Accession"], ev["aa_index"])]
    ev["Proteoform_ID"] = (ev["Accession"] + "_" + ev["residue"]
                           + ev["aa_index"].astype(str) + ev["terminus"])
    return ev


ISD_COLUMNS = ["Proteoform_ID", "Accession", "aa_index", "terminus", "Precursor.Id",
               "Stripped.Sequence", "RT", "Precursor.Quantity", "isd_parent",
               "isd_parent_sequence", "isd_parent_RT", "isd_parent_quantity",
               "intensity_ratio", "event_removed"]


def flag_in_source_decay(rep: pd.DataFrame, events: pd.DataFrame,
                         min_support: int = 2, rt_window: float = 0.05):
    """Flag supporting precursors explained as in-source fragments of a longer,
    co-eluting, more intense identified peptide; drop events left below
    ``min_support``.

    Returns (events_kept, flags, stats). ``flags`` has one row per
    (Proteoform_ID, Precursor.Id) of the INPUT events with the boolean
    ``isd_flag`` and the parent precursor; ``events_kept`` carries the
    recomputed ``n_support`` (distinct sequences with an unflagged precursor).
    """
    stats = {"N": {}, "C": {}}
    need = {"Precursor.Id", "Stripped.Sequence", "RT", "Precursor.Quantity"}
    if events.empty or not need.issubset(rep.columns):
        if not events.empty:
            log(f"  WARNING: ISD filter skipped, report lacks {sorted(need - set(rep.columns))}")
        flags = pd.DataFrame(columns=["Proteoform_ID", "Precursor.Id", "isd_flag", "isd_parent"])
        return events, flags, stats
    prec = (rep.groupby("Precursor.Id")
               .agg(seq=("Stripped.Sequence", "first"), rt=("RT", "median"),
                    q=("Precursor.Quantity", "median")))
    by_seq = {}
    for pid, sq, rt, q in zip(prec.index, prec["seq"], prec["rt"], prec["q"]):
        by_seq.setdefault(sq, []).append((pid, rt, q))
    seqs = list(by_seq)
    # support precursors are those mapped to the event's protein (as in the pf
    # matrix); parents may come from any protein
    own = {}
    for acc, sq, pid in rep[["Accession", "Stripped.Sequence", "Precursor.Id"]].drop_duplicates().itertuples(index=False):
        own.setdefault((acc, sq), set()).add(pid)
    # parents of a peptide of length L: longer sequences whose last (N events) or
    # first (C events) L residues equal it; indexed once per needed length
    idx = {}

    def parents(sq, term):
        key = (len(sq), term)
        if key not in idx:
            L, d = len(sq), {}
            for t in seqs:
                if len(t) > L:
                    d.setdefault(t[-L:] if term == "N" else t[:L], []).append(t)
            idx[key] = d
        return idx[key].get(sq, [])

    rows = []
    for pfid, acc, p, term, sq in events[["Proteoform_ID", "Accession", "aa_index", "terminus",
                                          "Stripped.Sequence"]].itertuples(index=False):
        cand = [c for par in parents(sq, term) for c in [(par,) + x for x in by_seq[par]]]
        for pid, rt, q in [x for x in by_seq.get(sq, []) if x[0] in own.get((acc, sq), ())]:
            co = [c for c in cand if abs(c[2] - rt) <= rt_window and c[3] > q]
            best = max(co, key=lambda c: c[3]) if co else None
            rows.append({"Proteoform_ID": pfid, "Accession": acc, "aa_index": p, "terminus": term,
                         "Precursor.Id": pid, "Stripped.Sequence": sq, "RT": rt,
                         "Precursor.Quantity": q, "isd_flag": best is not None,
                         "isd_parent": best[1] if best else "",
                         "isd_parent_sequence": best[0] if best else "",
                         "isd_parent_RT": best[2] if best else float("nan"),
                         "isd_parent_quantity": best[3] if best else float("nan"),
                         "intensity_ratio": (best[3] / q if best and q else float("nan"))})
    flags = pd.DataFrame(rows)
    ok = flags[~flags["isd_flag"]].groupby("Proteoform_ID")["Stripped.Sequence"].nunique()
    ev = events.copy()
    ev["n_support"] = ev["Proteoform_ID"].map(ok).fillna(0).astype(int)
    keep = ev["n_support"] >= min_support
    removed = set(ev.loc[~keep, "Proteoform_ID"])
    flags["event_removed"] = flags["Proteoform_ID"].isin(removed)
    for term in ("N", "C"):
        f = flags[flags["terminus"] == term]
        e = ev[ev["terminus"] == term]
        stats[term] = {"events_in": int(e["Proteoform_ID"].nunique()),
                       "precursors_in": int(len(f)),
                       "precursors_flagged": int(f["isd_flag"].sum()),
                       "events_with_flag": int(f.loc[f["isd_flag"], "Proteoform_ID"].nunique()),
                       "events_removed": int(len({x for x in removed if x.endswith(term)}))}
    return ev[keep].copy(), flags, stats


def _per_protein_proteoform_count(rep: pd.DataFrame,
                                  events: pd.DataFrame) -> pd.DataFrame:
    """1 baseline proteoform + one per distinct event site (an aa_index seen
    at both N and C counts as one cleavage point)."""
    if events.empty:
        sites = pd.DataFrame(columns=["Accession", "aa_index"])
    else:
        sites = events[["Accession", "aa_index"]].drop_duplicates()
    n_pf = (sites.groupby("Accession").size().add(1)
                  .rename("n_proteoforms").reset_index()
            if not sites.empty
            else pd.DataFrame(columns=["Accession", "n_proteoforms"]))
    all_acc = rep["Accession"].unique()
    n_pf = (pd.DataFrame({"Accession": all_acc})
              .merge(n_pf, on="Accession", how="left")
              .fillna({"n_proteoforms": 1}))
    n_pf["n_proteoforms"] = n_pf["n_proteoforms"].astype(int)
    return n_pf


def build_pf_matrix(rep: pd.DataFrame, events: pd.DataFrame,
                    n_pf: pd.DataFrame, flags: pd.DataFrame | None = None) -> pd.DataFrame:
    """Precursor-level proteoform matrix with one column per run."""
    if events.empty:
        return pd.DataFrame()

    # Quantities for the supporting peptides only
    keep_pep = events[["Accession", "Stripped.Sequence"]].drop_duplicates()
    sub = rep.merge(keep_pep, on=["Accession", "Stripped.Sequence"], how="inner")

    prec_cols = ["Accession", "Stripped.Sequence", "Modified.Sequence",
                 "Precursor.Id", "Precursor.Charge", "Enzyme", "Run",
                 "Precursor.Quantity"]
    prec_cols = [c for c in prec_cols if c in sub.columns]
    sub = sub[prec_cols].copy()

    grp_keys = [c for c in ("Accession", "Stripped.Sequence", "Modified.Sequence",
                            "Precursor.Id", "Precursor.Charge", "Enzyme")
                if c in sub.columns]
    mat = (sub.groupby(grp_keys + ["Run"], as_index=False)
              ["Precursor.Quantity"].sum())
    mat = mat.pivot_table(index=grp_keys, columns="Run",
                          values="Precursor.Quantity", aggfunc="sum")
    mat.columns.name = None
    mat = mat.reset_index()

    pf = events.merge(mat, on=["Accession", "Stripped.Sequence"], how="left")
    pf = pf.merge(n_pf, on="Accession", how="left")
    if flags is not None and "Precursor.Id" in pf.columns:
        f = flags[["Proteoform_ID", "Precursor.Id", "isd_flag", "isd_parent"]]
        pf = pf.merge(f, on=["Proteoform_ID", "Precursor.Id"], how="left")
        pf["isd_flag"] = pf["isd_flag"].fillna(False).astype(bool)
        pf["isd_parent"] = pf["isd_parent"].fillna("")

    front_cols = ["Proteoform_ID", "Accession", "aa_index", "terminus",
                  "n_support", "n_proteoforms",
                  "Precursor.Id", "Precursor.Charge",
                  "Modified.Sequence", "Stripped.Sequence", "Enzyme",
                  "isd_flag", "isd_parent"]
    front_cols = [c for c in front_cols if c in pf.columns]
    run_cols = [c for c in pf.columns if c not in front_cols]
    pf = pf[front_cols + run_cols]
    pf = pf.sort_values(
        by=["n_proteoforms", "Accession", "aa_index", "terminus",
            "Stripped.Sequence", "Precursor.Charge"]
            if "Precursor.Charge" in pf.columns else
           ["n_proteoforms", "Accession", "aa_index", "terminus",
            "Stripped.Sequence"],
        ascending=[False, True, True, True, True, True]
            if "Precursor.Charge" in pf.columns else
           [False, True, True, True, True],
    )
    return pf


ENZ_COLORS: dict = {
    "Trypsin":      "#4C72B0",
    "Chymotrypsin": "#DD8452",
    "GluC":         "#55A868",
    "Glu_C":        "#55A868",
    "Elastase":     "#C44E52",
    "AspN":         "#8172B2",
    "Asp_N":        "#8172B2",
    "combined":     "#444444",
    "semi":         "#2ca02c",
    "theirs":       "#bcbd22",
    "shared":       "#cccccc",
    "other":        "#8c8c8c",
}

# Hatch styles by number of synergizing enzymes (n>=2). Singletons get no hatch.
COMBO_HATCH = {2: "//", 3: "xx", 4: "++", 5: "**", 6: "OO"}
SINGLE_ALPHA = 0.65   # transparency applied to single-enzyme fills
COMBO_ALPHA  = 0.85   # transparency for blended combination fills


def _hex_to_rgb(h: str):
    h = h.lstrip("#")
    return tuple(int(h[i:i+2], 16) / 255.0 for i in (0, 2, 4))


def _blend_enzyme_colors(enzymes) -> tuple:
    """Average RGB of the enzyme palette colors so each combination gets a
    distinct, semantically-meaningful fill (e.g. Trypsin+Chymotrypsin =
    blue-orange blend, Chymotrypsin+Elastase = orange-red blend, etc.)."""
    rgbs = [_hex_to_rgb(ENZ_COLORS.get(e, ENZ_COLORS["other"])) for e in enzymes]
    r = sum(c[0] for c in rgbs) / len(rgbs)
    g = sum(c[1] for c in rgbs) / len(rgbs)
    b = sum(c[2] for c in rgbs) / len(rgbs)
    return (r, g, b)


def _enzyme_combo_counts(events: pd.DataFrame,
                         rep: pd.DataFrame):
    """For each Proteoform_ID, find the set of enzymes among its supporting
    peptide observations and count Proteoform_IDs per enzyme-combination.

    Returns (counts: pd.Series indexed by label, sets: list[tuple]) where
    ``sets[i]`` is the sorted tuple of enzymes for the i-th category.
    """
    if events.empty or "Enzyme" not in rep.columns:
        return pd.Series(dtype=int), []
    enz_per_pep = (rep.groupby(["Accession", "Stripped.Sequence"])["Enzyme"]
                      .agg(lambda x: set(x)))
    sup = events[["Proteoform_ID", "Accession", "Stripped.Sequence"]].copy()
    sup["enz_set"] = [enz_per_pep.get((a, s), set())
                      for a, s in zip(sup["Accession"], sup["Stripped.Sequence"])]
    combos = (sup.groupby("Proteoform_ID")["enz_set"]
                 .apply(lambda lists: tuple(sorted(frozenset().union(*lists)))))
    labels = combos.apply(lambda t: " + ".join(t) if t else "unknown")
    counts = labels.value_counts()
    label_to_set = dict(zip(labels.values, combos.values))
    # Order by (number of enzymes ASC, count DESC) so all single-enzyme
    # categories come first, then 2-enzyme combinations, then 3-enzyme, ...
    ordered_labels = sorted(
        counts.index,
        key=lambda lbl: (len(label_to_set[lbl]), -counts[lbl], lbl),
    )
    counts = counts.loc[ordered_labels]
    sets = [label_to_set[lbl] for lbl in ordered_labels]
    return counts, sets


def _style_for_combo(enzyme_tuple) -> tuple:
    """Return (facecolor_rgba, hatch, edgecolor) for an enzyme combination.

    Single enzymes use their palette color at SINGLE_ALPHA transparency.
    Multi-enzyme combinations get a distinct color by averaging the RGB of
    the participating enzymes, plus a hatch pattern keyed to combo size.
    """
    n = len(enzyme_tuple)
    if n == 0:
        return (*_hex_to_rgb(ENZ_COLORS["other"]), SINGLE_ALPHA), "", "white"
    if n == 1:
        rgb = _hex_to_rgb(ENZ_COLORS.get(enzyme_tuple[0], ENZ_COLORS["other"]))
        return (*rgb, SINGLE_ALPHA), "", "white"
    rgb = _blend_enzyme_colors(enzyme_tuple)
    hatch = COMBO_HATCH.get(n, "..")
    return (*rgb, COMBO_ALPHA), hatch, "black"


def write_plot(events: pd.DataFrame, n_pf: pd.DataFrame,
               rep: pd.DataFrame, out_pdf: str) -> None:
    with PdfPages(out_pdf) as pdf:
        # 1) Multiplicity histogram (with inset pie of class shares)
        fig, ax = plt.subplots(figsize=(7, 5))
        counts = n_pf["n_proteoforms"].value_counts().sort_index()
        ax.bar(counts.index, counts.values, color="steelblue", edgecolor="black")
        ax.set_xlabel("Differentiable sequence proteoforms per protein")
        ax.set_ylabel("Number of proteins")
        ax.set_title("Proteoform multiplicity distribution")
        ax.spines[["top", "right"]].set_visible(False)

        # Inset pie: proteins grouped into 1 / 2 / 3-10 / 11-100 / 100+
        nvals = n_pf["n_proteoforms"].to_numpy()
        bins = [
            ("1",       int(((nvals == 1)).sum()),                  "#cccccc"),
            ("2",       int(((nvals == 2)).sum()),                  "#9ecae1"),
            ("3–10",    int(((nvals >= 3) & (nvals <= 10)).sum()),  "#4292c6"),
            ("11–100",  int(((nvals >= 11) & (nvals <= 100)).sum()),"#08519c"),
        ]
        over = int((nvals > 100).sum())
        if over:
            bins.append((">100", over, "#08306b"))
        bins = [b for b in bins if b[1] > 0]
        total = sum(b[1] for b in bins)
        #inset = ax.inset_axes([0.55, 0.45, 0.42, 0.50])
        inset = ax.inset_axes([0.42, 0.32, 0.56, 0.66])
        wedges, _ = inset.pie(
            [b[1] for b in bins],
            colors=[b[2] for b in bins],
            startangle=90, counterclock=False,
            wedgeprops=dict(edgecolor="white", linewidth=0.8),
        )
        inset.set_title("Proteins by proteoform class", fontsize=9, pad=2)
        inset.legend(
            wedges,
            [f"{lbl}: {n} ({n/total:.0%})" for lbl, n, _ in bins],
            loc="center left", bbox_to_anchor=(1.0, 0.5),
            fontsize=7, frameon=False, handlelength=1.2, borderaxespad=0,
        )
        fig.tight_layout(); pdf.savefig(fig); plt.close(fig)

        if events.empty:
            return

        # 2) Enzyme cooperativity per Proteoform_ID — pie chart
        combo_counts, combo_sets = _enzyme_combo_counts(events, rep)
        if not combo_counts.empty:
            styles = [_style_for_combo(t) for t in combo_sets]
            faces  = [s[0] for s in styles]
            hatches = [s[1] for s in styles]
            edges  = [s[2] for s in styles]
            total  = int(combo_counts.sum())

            fig, ax = plt.subplots(figsize=(7.5, 6))
            wedges, _ = ax.pie(
                combo_counts.values, colors=faces, startangle=90,
                wedgeprops=dict(linewidth=0.8),
            )
            for w, h, ec in zip(wedges, hatches, edges):
                w.set_hatch(h)
                w.set_edgecolor(ec)
            legend_labels = [f"{lbl}  ({n}, {n/total:.0%})"
                             for lbl, n in combo_counts.items()]
            ax.legend(wedges, legend_labels, title="Enzyme support",
                      loc="center left", bbox_to_anchor=(1.02, 0.5),
                      frameon=False, fontsize=8)
            ax.set_title(f"Enzyme cooperativity per Proteoform_ID  (N={total})")
            fig.tight_layout(); pdf.savefig(fig); plt.close(fig)

            # 3) Same data, horizontal bar chart
            fig, ax = plt.subplots(figsize=(8, max(3, 0.35 * len(combo_counts))))
            order = combo_counts.iloc[::-1]
            order_sets = list(reversed(combo_sets))
            order_styles = [_style_for_combo(t) for t in order_sets]
            bars = ax.barh(
                range(len(order)), order.values,
                color=[s[0] for s in order_styles],
                edgecolor=[s[2] for s in order_styles],
                linewidth=0.8,
            )
            for b, st in zip(bars, order_styles):
                b.set_hatch(st[1])
            ax.set_yticks(range(len(order)))
            ax.set_yticklabels(order.index, fontsize=8)
            for i, v in enumerate(order.values):
                ax.text(v, i, f" {v} ({v/total:.0%})",
                        va="center", fontsize=8)
            ax.set_xlabel("Number of Proteoform_IDs")
            ax.set_title("Enzyme cooperativity per Proteoform_ID")
            ax.spines[["top", "right"]].set_visible(False)

            # Synergy-size legend (hatch key)
            from matplotlib.patches import Patch
            sizes_present = sorted({len(t) for t in combo_sets if len(t) >= 2})
            handles = [
                Patch(facecolor=(0.7, 0.7, 0.7, COMBO_ALPHA),
                      hatch=COMBO_HATCH.get(n, ".."),
                      edgecolor="black",
                      label=f"{n}-enzyme synergy (hatch)")
                for n in sizes_present
            ]
            handles.insert(0, Patch(facecolor="white", edgecolor="white",
                                    label="combo color = blend of enzymes"))
            ax.legend(handles=handles, loc="lower right", frameon=False,
                      fontsize=8, title="Legend")
            fig.tight_layout(); pdf.savefig(fig); plt.close(fig)

        # 4) Top-25 proteoform-rich proteins: event map
        top = (n_pf[n_pf["n_proteoforms"] > 1]
               .sort_values("n_proteoforms", ascending=False).head(25))
        if top.empty:
            return
        ev_top = events[events["Accession"].isin(top["Accession"])]
        evt_uniq = ev_top.drop_duplicates(["Accession", "aa_index", "terminus"])
        order = (top.set_index("Accession")["n_proteoforms"]
                    .sort_values().index.tolist())
        fig, ax = plt.subplots(figsize=(8, max(4, 0.3 * len(order))))
        ypos = {a: i for i, a in enumerate(order)}
        for term, marker, color in (("N", "o", "tab:blue"),
                                    ("C", "s", "tab:red")):
            s = evt_uniq[evt_uniq["terminus"] == term]
            if s.empty: continue
            ax.scatter(s["aa_index"], [ypos[a] for a in s["Accession"]],
                       c=color, marker=marker, s=30 + 6 * s["n_support"],
                       alpha=0.75, label=f"{term}-term")
        ax.set_yticks(range(len(order)))
        ax.set_yticklabels(order, fontsize=7)
        ax.set_xlabel("Residue position")
        ax.set_title("Proteoform events in top proteoform-rich proteins")
        ax.legend(loc="upper right", frameon=False)
        ax.spines[["top", "right"]].set_visible(False)
        fig.tight_layout(); pdf.savefig(fig); plt.close(fig)


# ── Public entry point ────────────────────────────────────────────────────────
def run(report_path: str, fasta_path: str, out_dir: str,
        min_support: int = 2, isd_filter: bool = True,
        isd_rt_window: float = 0.05) -> dict:
    global _LOG_FH
    os.makedirs(out_dir, exist_ok=True)
    out_log = os.path.join(out_dir, "report.pf_log.txt")
    _LOG_FH = open(out_log, "w")
    try:
        log(f"[pf] report     : {report_path}")
        log(f"[pf] fasta      : {fasta_path}")
        log(f"[pf] out_dir    : {out_dir}")
        log(f"[pf] min_support: {min_support}")
        log(f"[pf] isd_filter : {'on' if isd_filter else 'off'}"
            + (f" (rt window {isd_rt_window} min)" if isd_filter else ""))

        fa = read_fasta(fasta_path)
        log(f"[pf] FASTA entries : {len(fa):,}")

        rep = pd.read_csv(report_path, sep="\t", low_memory=False)
        log(f"[pf] report rows   : {len(rep):,}")

        rep = map_peptides_to_proteins(rep, fa)
        log(f"[pf] mapped rows   : {len(rep):,}")
        rep = annotate_terminus_specificity(rep)

        events = detect_events(rep, min_support=min_support)
        flags, isd_stats = None, None
        if isd_filter:
            n_before = events[["Accession", "aa_index", "terminus"]].drop_duplicates().shape[0]
            events, flags, isd_stats = flag_in_source_decay(
                rep, events, min_support=min_support, rt_window=isd_rt_window)
            for term, st in isd_stats.items():
                if st:
                    log(f"[pf] ISD {term}-terminal: {st['precursors_flagged']} of "
                        f"{st['precursors_in']} supporting precursors flagged; "
                        f"{st['events_with_flag']} of {st['events_in']} events affected, "
                        f"{st['events_removed']} removed")
            log(f"[pf] ISD events: {n_before} before, "
                f"{events[['Accession', 'aa_index', 'terminus']].drop_duplicates().shape[0]} after")
        n_sites = events[["Accession", "aa_index"]].drop_duplicates().shape[0]
        n_proteins = events["Accession"].nunique() if not events.empty else 0
        log(f"[pf] events: {n_sites} sites in {n_proteins} proteins "
            f"(min_support={min_support})")

        n_detected = int(rep["Accession"].nunique())
        pct = round(100 * n_proteins / n_detected) if n_detected else 0
        log(f"[pf] proteins with sequence divergence events: "
            f"{n_proteins}/{n_detected} = {pct}% of detected proteins")

        n_pf = _per_protein_proteoform_count(rep, events)

        pf_matrix = build_pf_matrix(rep, events, n_pf, flags)
        out_isd = None
        if flags is not None:
            out_isd = os.path.join(out_dir, "report.pf_isd.tsv")
            fl = flags[flags["isd_flag"]] if len(flags) else flags
            fl.reindex(columns=ISD_COLUMNS).to_csv(out_isd, sep="\t", index=False)
            log(f"[pf] wrote {out_isd}  ({len(fl):,} flagged precursors)")
        out_matrix = os.path.join(out_dir, "report.pf_matrix.tsv")
        pf_matrix.to_csv(out_matrix, sep="\t", index=False)
        log(f"[pf] wrote {out_matrix}  ({len(pf_matrix):,} rows)")

        out_plot = os.path.join(out_dir, "report.pf_plot.pdf")
        write_plot(events, n_pf, rep, out_plot)
        log(f"[pf] wrote {out_plot}")
        log(f"[pf] wrote {out_log}")

        return {"matrix": out_matrix, "plot": out_plot, "log": out_log, "isd": out_isd,
                "isd_stats": isd_stats,
                "n_events": int(n_sites),
                "n_proteins_with_events": int((n_pf["n_proteoforms"] > 1).sum()),
                "n_proteins_detected": n_detected,
                "pct_proteins_with_events": pct}
    finally:
        if _LOG_FH is not None:
            _LOG_FH.close()
            _LOG_FH = None


# ── CLI ───────────────────────────────────────────────────────────────────────
def main(argv: Iterable[str] | None = None) -> int:
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("report", nargs="?",
                    default=os.path.join(here, "report_annotated.tsv"))
    ap.add_argument("fasta", nargs="?",
                    default=os.path.join(here,
                        "UP000000750_301447_1perGene_SpyoSF370_Date202601.fasta"))
    ap.add_argument("out_dir", nargs="?", default=here)
    ap.add_argument("--min-support", type=int, default=2,
                    help="minimum distinct peptides sharing a free terminus "
                         "(default: 2)")
    ap.add_argument("--no-isd-filter", action="store_true",
                    help="do not flag or remove in-source-decay products")
    ap.add_argument("--isd-rt-window", type=float, default=0.05,
                    help="max apex RT difference (min) to a longer parent for an "
                         "ISD flag (default: 0.05)")
    args = ap.parse_args(argv)
    run(args.report, args.fasta, args.out_dir, min_support=args.min_support,
        isd_filter=not args.no_isd_filter, isd_rt_window=args.isd_rt_window)
    return 0


if __name__ == "__main__":
    sys.exit(main())
