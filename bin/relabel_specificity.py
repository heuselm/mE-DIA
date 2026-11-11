#!/usr/bin/env python3
"""
Recompute the Specificity column of an existing me-DIA report.

Reports annotated before 2026-09-19 carry a Specificity label that marks most
fully-specific peptides ``semi_specific`` (run_pipeline._digest_semi bug). This
writes a NEW file next to the input, or into --out-dir (never overwrites):

    <input stem>.specificity_fixed.tsv

with ``Specificity`` recomputed from protein sequence (bin/specificity.py) and
the old label kept as ``Specificity_pre20260919``. A per-label summary is
printed and written to ``<input stem>.specificity_fixed.summary.tsv``.

Usage:
    python relabel_specificity.py REPORT.tsv --fasta FASTA --enzymes TCE
    (enzymes = every enzyme of the experiment; the multiplex rule applies)
"""
from __future__ import annotations

import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import specificity as sp  # noqa: E402


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("report")
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--enzymes", required=True, help="enzyme codes of the experiment, e.g. TCE")
    ap.add_argument("--met-excision", action="store_true")
    ap.add_argument("--out-dir", default=None,
                    help="write here instead of next to the input (use for read-only data)")
    ap.add_argument("--q", type=float, default=0.01, help="Q.Value cut for the summary only")
    a = ap.parse_args(argv)

    stem = a.report[:-4] if a.report.endswith(".tsv") else a.report
    if a.out_dir:
        os.makedirs(a.out_dir, exist_ok=True)
        stem = os.path.join(a.out_dir, os.path.basename(stem))
    out = stem + ".specificity_fixed.tsv"
    summ = stem + ".specificity_fixed.summary.tsv"
    if os.path.exists(out):
        print(f"exists, not overwriting: {out}")
        return 1

    df = pd.read_csv(a.report, sep="\t", low_memory=False)
    ids = "Protein.Ids" if "Protein.Ids" in df.columns else "Protein.Group"
    cand = (df.drop_duplicates("Stripped.Sequence").set_index("Stripped.Sequence")[ids]
              .astype(str).str.split(";").to_dict())
    lab = sp.classify(df["Stripped.Sequence"].unique(), sp.load_proteins(a.fasta),
                      list(a.enzymes), cand, met_excision=a.met_excision)
    if "Specificity" in df.columns:
        df["Specificity_pre20260919"] = df["Specificity"]
    df["Specificity"] = df["Stripped.Sequence"].map(lab)

    df.to_csv(out + ".part", sep="\t", index=False)
    os.replace(out + ".part", out)

    q = df[df["Q.Value"] <= a.q].drop_duplicates("Precursor.Id") if "Q.Value" in df.columns else df
    if "Specificity_pre20260919" in q.columns:
        t = pd.crosstab(q["Specificity_pre20260919"].fillna("NA"), q["Specificity"]).reset_index()
    else:
        t = q["Specificity"].value_counts().rename_axis("Specificity").reset_index()
    t.to_csv(summ + ".part", sep="\t", index=False)
    os.replace(summ + ".part", summ)
    print(f"wrote {out}\nprecursors at Q.Value <= {a.q} (old label x new label):")
    print(t.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
