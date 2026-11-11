#!/usr/bin/env python3
"""Checks that proteoform_analysis knows every pipeline enzyme's cleavage rule.
Run: python bin/test_proteoform_rules.py

The report's Enzyme column holds ENZYME_REGISTRY names (run_pipeline.py). An
enzyme proteoform_analysis has no rule for counts every terminus as specific,
so its peptides could never call an event; LysC, LysN, ArgC, AspN and
Trypsin/P were in that state until 2026-10.
"""
import contextlib
import io
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import run_pipeline as rp           # noqa: E402
import proteoform_analysis as pf    # noqa: E402

FAIL = []


def check(name, cond, detail=""):
    print(f"  {'PASS' if cond else 'FAIL'}  {name}{('  ' + str(detail)) if detail else ''}")
    if not cond:
        FAIL.append(name)
        # Without this a failing check only sets an exit code, which pytest never
        # sees, so the check would pass silently under `pytest bin/`.
        if "PYTEST_CURRENT_TEST" in os.environ:
            raise AssertionError(f"{name} {detail}")


def want(e):
    return {"sense": e.sense, "cut": e.cut_res, "block": e.noclv}


def test_table_matches_registry():
    print("every registry enzyme has its rule in the table, under name and code")
    for code, e in rp.ENZYME_REGISTRY.items():
        check(f"{code} {e.name!r}: by name", pf.ENZYME_RULES.get(e.name) == want(e),
              pf.ENZYME_RULES.get(e.name))
        check(f"{code}: by code", pf.ENZYME_RULES.get(code) == want(e),
              pf.ENZYME_RULES.get(code))


def test_terminus_calls():
    print("\nterminus specificity under the formerly missing enzymes")
    seq = "MAAAAKPEPTIDERWWWWDQQQQ"
    #        1-based: K at 6, P at 7 (K|P blocked), R at 14, D at 19
    check("LysC: N-term after K is specific, after another residue is not",
          pf.nterm_specific("MAAAAKSEPTIDE", 7, "LysC")
          and not pf.nterm_specific("MAAAAKSEPTIDE", 8, "LysC"))
    check("LysC: K|P is not a cleavage site", not pf.nterm_specific(seq, 7, "LysC"))
    check("ArgC: C-term at R is specific, at K is not",
          pf.cterm_specific(seq, 14, "ArgC", len(seq))
          and not pf.cterm_specific(seq, 6, "ArgC", len(seq)))
    check("AspN (N-terminal cutter): N-term at D is specific, elsewhere not",
          pf.nterm_specific(seq, 19, "AspN") and not pf.nterm_specific(seq, 18, "AspN"))
    check("LysN: N-term at K is specific", pf.nterm_specific(seq, 6, "LysN"))
    check("Trypsin/P: K|P is a site (no proline rule)", pf.nterm_specific(seq, 7, "Trypsin/P"))


def test_registry_fallback_and_warning():
    print("\nan --enzyme-config enzyme: registry fallback; an unknown one: warned")
    import pandas as pd
    saved = dict(rp.ENZYME_REGISTRY)
    rp.ENZYME_REGISTRY["Q"] = rp.EnzymeDef("Q", "Qase", "Q*", 1, "Q", "")
    try:
        check("not in the table, found through the registry by name and code",
              "Qase" not in pf.ENZYME_RULES and pf.rule_for("Qase") == pf.rule_for("Q")
              == {"sense": 1, "cut": "Q", "block": ""}, pf.rule_for("Qase"))
        check("... and used: N-term after Q specific, after A not",
              pf.nterm_specific("MQAAAA", 3, "Qase") and not pf.nterm_specific("MQAAAA", 4, "Qase"))
    finally:
        rp.ENZYME_REGISTRY.clear()
        rp.ENZYME_REGISTRY.update(saved)
    check("unknown name: no rule", pf.rule_for("Mysteryase") is None)
    rep = pd.DataFrame({"Sequence": ["MKAAAAKAAAA"] * 3, "Start": [3, 3, 3], "Stop": [7, 7, 7],
                        "Length": [11] * 3, "Enzyme": ["LysC", "Mysteryase", "Mysteryase"]})
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        out = pf.annotate_terminus_specificity(rep)
    check("an enzyme without a rule is named in a warning, with its row count",
          "Mysteryase" in buf.getvalue() and "2 rows" in buf.getvalue()
          and "LysC" not in buf.getvalue(), buf.getvalue().strip())
    check("LysC row: an N-terminus right after K is specific",
          bool(out.loc[0, "N_specific"]))


TESTS = (test_table_matches_registry, test_terminus_calls, test_registry_fallback_and_warning)


def main():
    for t in TESTS:
        t()
    print(f"\n{'all checks passed' if not FAIL else f'{len(FAIL)} FAILED: {FAIL}'}")
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
