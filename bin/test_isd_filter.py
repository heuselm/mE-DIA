#!/usr/bin/env python3
"""Checks for the in-source-decay (ISD) filter of proteoform_analysis (issue #42).
Run: python bin/test_isd_filter.py   (or pytest bin/)

Synthetic protein with an N-terminal event at P7, supported by two sequences:
  PGTIKPHTK   tryptic C-terminus, non-specific N-terminus (after K|P)
  PGTIKPHTKF  chymotryptic C-terminus, non-specific N-terminus
and a longer tryptic peptide AKPGTIKPHTK that ends like PGTIKPHTK (its y-type parent).
"""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import proteoform_analysis as pf    # noqa: E402

FAIL = []
SEQ = "MSGQVLAKPGTIKPHTKFESEVYILSKDEGGRHTPFFKGYRPQFYFRTTDVTGTIELPEGVEMVMPGDNIKMVVTLIHPIAMDDGLRFAIREGGRTVGAGVVAKVLG"


def check(name, cond, detail=""):
    print(f"  {'PASS' if cond else 'FAIL'}  {name}{('  ' + str(detail)) if detail else ''}")
    if not cond:
        FAIL.append(name)
        if "PYTEST_CURRENT_TEST" in os.environ:
            raise AssertionError(f"{name} {detail}")


def report(rows):
    """rows: (sequence, charge, enzyme, RT, quantity) -> annotated report, two runs."""
    out = []
    for run in ("r1", "r2"):
        for s, z, e, rt, q in rows:
            out.append({"Run": run, "Protein.Ids": "P1", "Stripped.Sequence": s, "Modified.Sequence": s,
                        "Precursor.Id": f"{s}{z}", "Precursor.Charge": z, "Enzyme": e, "RT": rt,
                        "Precursor.Quantity": q})
    fa = pd.DataFrame({"Accession": ["P1"], "Sequence": [SEQ], "Length": [len(SEQ)]})
    rep = pf.map_peptides_to_proteins(pd.DataFrame(out), fa)
    return pf.annotate_terminus_specificity(rep)


def call(rows, **kw):
    rep = report(rows)
    ev = pf.detect_events(rep, min_support=2)
    kept, flags, stats = pf.flag_in_source_decay(rep, ev, min_support=2, **kw)
    return ev, kept, flags, stats


BASE = [("AKPGTIKPHTK", 3, "Trypsin", 5.57, 1.0e7),     # parent of PGTIKPHTK
        ("PGTIKPHTK", 2, "Trypsin", 5.58, 5.0e5),       # co-elutes, 20x weaker
        ("PGTIKPHTKF", 2, "Chymotrypsin", 11.70, 6.0e5)]  # separate elution


def test_flag_and_remove():
    print("co-eluting, more intense y-type parent -> flagged; event drops below min_support")
    ev, kept, flags, st = call(BASE)
    check("event called before the filter", (ev.Proteoform_ID == "P1_P9N").any(), sorted(ev.Proteoform_ID.unique()))
    f = flags.set_index("Precursor.Id")
    check("PGTIKPHTK2 flagged", bool(f.loc["PGTIKPHTK2", "isd_flag"]))
    check("parent recorded", f.loc["PGTIKPHTK2", "isd_parent"] == "AKPGTIKPHTK3", f.loc["PGTIKPHTK2", "isd_parent"])
    check("PGTIKPHTKF2 not flagged (no co-eluting parent)", not bool(f.loc["PGTIKPHTKF2", "isd_flag"]))
    check("event removed", not (kept.Proteoform_ID == "P1_P9N").any())
    check("stats count it", st["N"]["events_removed"] == 1 and st["N"]["precursors_flagged"] == 1, st["N"])


def test_independent_support_keeps_event():
    print("a third, unflagged sequence keeps the event; n_support counts unflagged sequences")
    rows = BASE + [("PGTIKPHTKFE", 2, "GluC", 9.46, 3.0e5)]
    _, kept, flags, _ = call(rows)
    e = kept[kept.Proteoform_ID == "P1_P9N"]
    check("event kept", len(e) > 0)
    check("n_support = 2", set(e.n_support) == {2}, set(e.n_support))


def test_not_flagged_cases():
    print("weaker parent, distant parent, or filter window 0 -> not flagged")
    weak = [(s, z, e, rt, 1.0e4 if s == "AKPGTIKPHTK" else q) for s, z, e, rt, q in BASE]
    check("weaker parent", not call(weak)[2].isd_flag.any())
    far = [(s, z, e, 5.70 if s == "AKPGTIKPHTK" else rt, q) for s, z, e, rt, q in BASE]
    check("apex 0.13 min apart", not call(far)[2].isd_flag.any())
    check("window 0.005 min", not call(BASE, rt_window=0.005)[2].isd_flag.any())


def test_c_terminal_mirror():
    print("C-terminal event: parent shares the N-terminus (b-type) and is flagged the same way")
    rows = [("GQVLAKPGTIKPHTKF", 3, "Chymotrypsin", 12.0, 8.0e6),   # parent: same N-term, longer C-term
            ("GQVLAKPGTIK", 2, "Elastase", 12.01, 2.0e5),           # C-term at K? specific for trypsin only
            ("SGQVLAKPGTIK", 2, "GluC", 8.0, 3.0e5)]
    rep = report(rows)
    ev = pf.detect_events(rep, min_support=1)
    kept, flags, st = pf.flag_in_source_decay(rep, ev, min_support=1)
    fc = flags[(flags.terminus == "C") & (flags["Precursor.Id"] == "GQVLAKPGTIK2")]
    check("C-terminal support flagged against same-N-terminus parent", len(fc) and bool(fc.isd_flag.iloc[0]),
          fc[["Proteoform_ID", "isd_flag", "isd_parent"]].to_dict("records"))


def test_matrix_columns():
    print("pf matrix carries isd_flag and isd_parent per precursor")
    rows = BASE + [("PGTIKPHTKFE", 2, "GluC", 9.46, 3.0e5)]
    rep = report(rows)
    ev = pf.detect_events(rep, min_support=2)
    kept, flags, _ = pf.flag_in_source_decay(rep, ev, min_support=2)
    m = pf.build_pf_matrix(rep, kept, pf._per_protein_proteoform_count(rep, kept), flags)
    check("columns present", {"isd_flag", "isd_parent"} <= set(m.columns), list(m.columns)[:14])
    r = m[m["Precursor.Id"] == "PGTIKPHTK2"]
    check("flagged precursor kept in matrix, marked", len(r) == 1 and bool(r.isd_flag.iloc[0]))


if __name__ == "__main__":
    for t in (test_flag_and_remove, test_independent_support_keeps_event, test_not_flagged_cases,
              test_c_terminal_mirror, test_matrix_columns):
        t()
    print(f"\n{'FAILED: ' + ', '.join(FAIL) if FAIL else 'all checks passed'}")
    sys.exit(1 if FAIL else 0)
