#!/usr/bin/env python3
"""Checks for semi-specific digest labelling. Run: python bin/test_specificity.py"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import run_pipeline as rp  # noqa: E402
import specificity as sp  # noqa: E402

FAIL = []


def check(name, cond):
    print(f"  {'PASS' if cond else 'FAIL'}  {name}")
    if not cond:
        FAIL.append(name)
        if "PYTEST_CURRENT_TEST" in os.environ:     # otherwise pytest never sees a failure
            raise AssertionError(name)


def digest_labels(seq, code, missed=1):
    import tempfile
    enz = rp.ENZYME_REGISTRY[code]
    with tempfile.TemporaryDirectory() as td:
        fa = os.path.join(td, "p.fasta")
        with open(fa, "w") as fh:
            fh.write(f">sp|P1|X\n{seq}\n")
        vf, vm = os.path.join(td, "v.fasta"), os.path.join(td, "v.tsv")
        rp.generate_semi_vfasta(fa, code, enz, vf, vm, min_len=7, max_len=30, missed=missed)
        import pandas as pd
        m = pd.read_csv(vm, sep="\t")
    return dict(zip(m["PeptideSequence"], m["Specificity"]))


def main():
    #        trypsin sites after K4 (idx 4), R12 (idx 12), K20 (idx 20); protein ends 29
    prot = "MAAAKDDDDDDDRGGGGGGGKLLLLLLLE"
    lab = digest_labels(prot, "T")
    print("generate_semi_vfasta labels")
    check("tryptic peptide inside a missed-cleavage span is fully specific",
          lab.get("DDDDDDDR") == "fully_specific")
    check("C-terminal tryptic peptide with 1 missed cleavage is fully specific",
          lab.get("GGGGGGGKLLLLLLLE") == "fully_specific")
    check("N-truncated peptide is semi specific", lab.get("DDDDDDR") == "semi_specific")
    check("C-truncated peptide is semi specific", lab.get("DDDDDDDRGGGG") == "semi_specific")
    check("every fully-specific label has both termini on a cleavage site",
          all(sp.classify([p], {"P1": prot}, ["T"])[p] == v for p, v in lab.items()))

    print("specificity.classify")
    c = sp.classify(["AAAKDDDDDDDR"], {"P1": prot}, ["T"])
    check("Met-excised N-terminal peptide is semi by default", c["AAAKDDDDDDDR"] == "semi_specific")
    c = sp.classify(["AAAKDDDDDDDR"], {"P1": prot}, ["T"], met_excision=True)
    check("... and fully with met_excision=True", c["AAAKDDDDDDDR"] == "fully_specific")
    p2 = {"P2": "AAEGGGGGGGKGGGEAA"}
    check("semi-tryptic peptide", sp.classify(["GGGGGGGKGGGE"], p2, ["T"])["GGGGGGGKGGGE"]
          == "semi_specific")
    check("multiplex: fully specific for GluC -> fully specific",
          sp.classify(["GGGGGGGKGGGE"], p2, ["T", "G"])["GGGGGGGKGGGE"] == "fully_specific")
    check("absent peptide -> not_found", sp.classify(["WWWWWWW"], p2, ["T"])["WWWWWWW"] == "not_found")
    check("peptide whose proteins are absent from the FASTA -> not_in_fasta",
          sp.classify(["GGGGGGGKGGGE"], p2, ["T"], {"GGGGGGGKGGGE": ["cont_crap|X|"]})["GGGGGGGKGGGE"]
          == "not_in_fasta")
    aspn = digest_labels("GGGGGGGDAAAAAAADKKKKKKK", "A")
    check("N-terminal cutter (AspN): fully-specific peptide labelled fully",
          aspn.get("DAAAAAAA") == "fully_specific")

    print(f"\n{'all checks passed' if not FAIL else f'{len(FAIL)} FAILED: {FAIL}'}")
    return 1 if FAIL else 0


def test_specificity():
    """Entry point for `pytest bin/`; the checks themselves live in main()."""
    FAIL.clear()
    assert main() == 0, f"failed checks: {FAIL}"


if __name__ == "__main__":
    sys.exit(main())
