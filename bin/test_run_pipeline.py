#!/usr/bin/env python3
"""Checks for the digest, virtual-FASTA and modification code in run_pipeline.

Run: python bin/test_run_pipeline.py   (no DIA-NN, no data, about a second)
Also collected by pytest, if installed: a failed check fails its test.

The semi-specific digest is checked against an independent brute-force
enumeration: every substring with at least one terminus on a cleavage site,
length within [min_len, max_len] and at most `missed` cleavage sites strictly
inside it; fully specific when both termini are sites.
"""
import contextlib
import io
import os
import random
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import run_pipeline as rp  # noqa: E402

FAIL = []
REG = rp.ENZYME_REGISTRY


def check(name, cond, detail=""):
    print(f"  {'PASS' if cond else 'FAIL'}  {name}{('  ' + str(detail)) if detail else ''}")
    if not cond:
        FAIL.append(name)
        if "PYTEST_CURRENT_TEST" in os.environ:
            raise AssertionError(f"{name} {detail}")


def quiet(fn, *args, **kw):
    """Call fn with its progress prints swallowed."""
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*args, **kw)


# ── independent reference implementation ────────────────────────────────────
def ref_sites(seq, enz):
    """Peptide boundaries from the rule alone: protein termini + cleavage sites."""
    s = {0, len(seq)}
    for i, aa in enumerate(seq):
        if aa not in enz.cut_res:
            continue
        if enz.sense == 1:
            if i + 1 == len(seq) or seq[i + 1] not in enz.noclv:
                s.add(i + 1)
        elif i > 0:
            s.add(i)
    return s


def ref_digest(seq, enz, min_len, max_len, missed):
    """Set of (peptide, label) for every qualifying occurrence in seq."""
    s = ref_sites(seq, enz)
    out = set()
    for a in range(len(seq)):
        for b in range(a + min_len, min(a + max_len, len(seq)) + 1):
            if a not in s and b not in s:
                continue
            if sum(1 for x in s if a < x < b) > missed:
                continue
            out.add((seq[a:b], "fully_specific" if a in s and b in s else "semi_specific"))
    return out


def digest(seq, code, min_len=7, max_len=30, missed=1, pid="P1"):
    return rp._digest_semi(seq, pid, REG[code], min_len, max_len, missed)


def write(path, text):
    with open(path, "w", newline="") as fh:
        fh.write(text)


def read_vfasta(path):
    recs, hdr = [], None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                hdr = line[1:]
            elif hdr is not None:
                recs.append((hdr, line))
                hdr = None
    return recs


def vfasta(td, fasta_text, code="T", tag="v", **kw):
    """Run generate_semi_vfasta; return (count, [(hdr, seq)], map DataFrame)."""
    import pandas as pd
    fa = os.path.join(td, f"{tag}.fasta")
    write(fa, fasta_text)
    vf, vm = os.path.join(td, f"{tag}.vfasta"), os.path.join(td, f"{tag}.map.tsv")
    n = quiet(rp.generate_semi_vfasta, fa, code, REG[code], vf, vm, **kw)
    m = pd.read_csv(vm, sep="\t", dtype=str, keep_default_na=False)
    return n, read_vfasta(vf), m


# ── tests ────────────────────────────────────────────────────────────────────
def test_registry():
    print("enzyme registry")
    check("the enzymes the GUI and docs name are registered",
          {"T", "T2", "C", "G", "E", "L", "LN", "R", "A", "M", "F"} <= set(REG), sorted(REG))
    check("every entry's key is its code", all(k == e.code for k, e in REG.items()))
    check("sense is 0 or 1", all(e.sense in (0, 1) for e in REG.values()))

    def cut_matches(e):
        terms = set(e.cut.split(","))
        want = {f"{r}*" if e.sense == 1 else f"*{r}" for r in e.cut_res}
        want |= {f"!*{r}" for r in e.noclv}
        return terms == want
    bad = [k for k, e in REG.items() if not cut_matches(e)]
    check("DIA-NN --cut agrees with cut_res / noclv / sense", not bad, bad)
    check("trypsin: C-terminal K/R, not before P",
          (REG["T"].sense, REG["T"].cut_res, REG["T"].noclv) == (1, "KR", "P"))
    check("trypsin/P has no proline rule", REG["T2"].noclv == "")
    check("C-terminal cutter: cterm = cut_res, nterm empty",
          REG["T"].cterm == "KR" and REG["T"].nterm == "")
    check("N-terminal cutter: nterm = cut_res, cterm empty",
          REG["A"].nterm == "D" and REG["A"].cterm == "" and REG["LN"].nterm == "K")


def test_cleavage_sites():
    print("_find_cleavage_sites")
    cases = [
        ("trypsin, K/R sites and protein termini", "MAAAKDDDDDDDRGGGGGGGKLLLLLLLE", "T",
         [0, 5, 13, 21, 29]),
        ("trypsin: K before P is not a site", "AAKPAARAA", "T", [0, 7, 9]),
        ("trypsin/P: K before P is a site", "AAKPAARAA", "T2", [0, 3, 7, 9]),
        ("C-terminal K is not listed twice", "AAAK", "T", [0, 4]),
        ("no cleavage residue: protein termini only", "AAAAAA", "T", [0, 6]),
        ("chymotrypsin: F before P blocked, W and Y cut", "AFPAWAYA", "C", [0, 5, 7, 8]),
        ("LysN cuts before K", "AKAAKK", "LN", [0, 1, 4, 5, 6]),
        ("LysN: K at position 0 does not duplicate 0", "KAAK", "LN", [0, 3, 4]),
        ("AspN cuts before D", "DAADAA", "A", [0, 3, 6]),
    ]
    for name, seq, code, want in cases:
        got = rp._find_cleavage_sites(seq, REG[code])
        check(name, got == want, got)

    rng = random.Random(7)
    bad = []
    for code in REG:
        for _ in range(30):
            seq = "".join(rng.choice("ACDEFGKLMPRVWY") for _ in range(rng.randint(1, 60)))
            got = rp._find_cleavage_sites(seq, REG[code])
            if got != sorted(ref_sites(seq, REG[code])):
                bad.append((code, seq))
    check("sorted, unique, termini included: matches the rule on 330 random sequences",
          not bad, bad[:3])


def test_digest():
    print("_digest_semi")
    prot = "MAAAKDDDDDDDRGGGGGGGKLLLLLLLE"     # tryptic sites 0, 5, 13, 21, 29
    out = digest(prot, "T", missed=1)
    peps = {p: s for p, _, s in out}
    check("returns (peptide, protein id, label) with the given id",
          all(len(t) == 3 and t[1] == "P1" for t in out))
    check("labels are fully_specific or semi_specific only",
          {s for _, _, s in out} <= {"fully_specific", "semi_specific"})
    check("fully-specific peptide with 0 missed cleavages", peps.get("DDDDDDDR") == "fully_specific")
    check("fully-specific peptide with 1 missed cleavage",
          peps.get("DDDDDDDRGGGGGGGK") == "fully_specific")
    check("2 missed cleavages excluded at missed=1", "DDDDDDDRGGGGGGGKLLLLLLLE" not in peps)
    check("... included at missed=2", "DDDDDDDRGGGGGGGKLLLLLLLE" in
          {p for p, _, _ in digest(prot, "T", missed=2)})
    check("... and 1 missed cleavage excluded at missed=0",
          "DDDDDDDRGGGGGGGK" not in {p for p, _, _ in digest(prot, "T", missed=0)})
    check("protein N-terminal peptide is fully specific", peps.get("MAAAKDDDDDDDR") == "fully_specific")
    check("N-truncated peptide is semi specific", peps.get("DDDDDDR") == "semi_specific")
    check("C-truncated peptide is semi specific", peps.get("GGGGGGGKLLLL") == "semi_specific")
    check("peptide with neither terminus on a site is absent", "DDDDDDRGG" not in peps)
    check("peptides shorter than min_len are absent", min(len(p) for p in peps) >= 7)
    short = digest(prot, "T", min_len=7, max_len=10, missed=2)
    check("peptides longer than max_len are absent",
          short and max(len(p) for p, _, _ in short) <= 10)
    check("LysN: fully specific peptides start with K",
          all(p[0] == "K" or prot.startswith(p) for p, _, s in digest(prot, "LN") if s == "fully_specific"))

    # the full set, against the brute-force reference, over enzymes and bounds
    rng = random.Random(11)
    seqs = [prot, "AAKPAARAAAKAAAAAAAAKRAAAAAAAAPKA", "DAAAAAAADAAPDAAAAAAAAAAD"]
    seqs += ["".join(rng.choice("ACDEFGKLMPRVWY") for _ in range(rng.randint(5, 70)))
             for _ in range(25)]
    bad, n = [], 0
    for code in REG:
        for missed in (0, 1, 2):
            for min_len, max_len in ((7, 30), (4, 12), (1, 70)):
                for seq in seqs:
                    n += 1
                    got = set((p, s) for p, _, s in digest(seq, code, min_len, max_len, missed))
                    if got != ref_digest(seq, REG[code], min_len, max_len, missed):
                        bad.append((code, missed, min_len, max_len, seq))
    check(f"peptide set and labels equal the brute-force reference ({n} cases: "
          f"{len(REG)} enzymes x missed 0-2 x 3 length windows x {len(seqs)} sequences)",
          not bad, bad[:2])


FASTA = (
    ">sp|P1|ONE_TEST first protein\n"
    "MAAAKDDDDDDDRGG\n"
    "gggggKLLLLLLLE\n"                 # lower case and line wrap: same protein
    ">tr|P2|TWO_TEST second protein\n"
    "WWDDDDDDDRWWFFFFFFFK\n"
)


def test_parse_fasta():
    print("_parse_fasta / _load_fasta_dict")
    text = (FASTA + ">P3 plain header, no pipes\r\nACDEF\r\nGHIK\r\n"
            ">sp|P4|LAST\nMKW")       # no trailing newline on the last record
    with tempfile.TemporaryDirectory() as td:
        fa = os.path.join(td, "p.fasta")
        write(fa, "stray line before any header\n" + text)
        a = rp._parse_fasta(fa)
        b = rp._load_fasta_dict(fa)
    want = {"P1": "MAAAKDDDDDDDRGGGGGGGKLLLLLLLE", "P2": "WWDDDDDDDRWWFFFFFFFK",
            "P3": "ACDEFGHIK", "P4": "MKW"}
    check("_parse_fasta: sp/tr accession, plain first token, joined, upper case, CRLF",
          a == want, a)
    check("_load_fasta_dict agrees on UniProt and plain headers", b == want, b)


def test_vfasta():
    print("generate_semi_vfasta")
    p1, p2 = "MAAAKDDDDDDDRGGGGGGGKLLLLLLLE", "WWDDDDDDDRWWFFFFFFFK"
    with tempfile.TemporaryDirectory() as td:
        n, recs, m = vfasta(td, FASTA, "T", tag="a")

        check("return value is the number of peptides", n == len(m) == len(recs), (n, len(m), len(recs)))
        check("map columns", list(m.columns) ==
              ["VirtualHeader", "PeptideSequence", "ProteinGroup", "Enzyme", "Specificity"])
        check("virtual FASTA and map list the same (header, peptide) pairs in the same order",
              recs == list(zip(m["VirtualHeader"], m["PeptideSequence"])))
        check("headers are vpep_<code>_<7-digit index>, consecutive from 0",
              list(m["VirtualHeader"]) == [f"vpep_T_{i:07d}" for i in range(n)])
        check("one entry per peptide, sorted",
              list(m["PeptideSequence"]) == sorted(set(m["PeptideSequence"])))
        check("Enzyme column is the display name", set(m["Enzyme"]) == {"Trypsin"})

        exp = {}
        for pid, seq in (("P1", p1), ("P2", p2)):
            for pep, lab in ref_digest(seq, REG["T"], 7, 30, 1):
                prots, labs = exp.setdefault(pep, (set(), set()))
                prots.add(pid)
                labs.add(lab)
        check("peptide set = union of both proteins' digests",
              set(m["PeptideSequence"]) == set(exp))
        row = m.set_index("PeptideSequence")
        check("ProteinGroup lists every protein whose digest yields the peptide",
              all(row.at[p, "ProteinGroup"] == ";".join(sorted(exp[p][0])) for p in exp))
        check("Specificity: fully specific if fully specific in any protein",
              all(row.at[p, "Specificity"] ==
                  ("fully_specific" if "fully_specific" in exp[p][1] else "semi_specific")
                  for p in exp))
        check("shared peptide: fully in P1, semi in P2 -> fully, group P1;P2",
              (row.at["DDDDDDDR", "ProteinGroup"], row.at["DDDDDDDR", "Specificity"])
              == ("P1;P2", "fully_specific"))
        check("length bounds hold", m["PeptideSequence"].str.len().between(7, 30).all())

        rev = ">tr|P2|TWO_TEST\n" + p2 + "\n>sp|P1|ONE_TEST\n" + p1 + "\n"
        _, recs_r, m_r = vfasta(td, rev, "T", tag="r")
        check("output does not depend on FASTA entry order",
              recs_r == recs and m_r.equals(m))

        _, _, w = vfasta(td, FASTA, "T", tag="w", protein_whitelist={"P2"})
        check("protein whitelist restricts the digest to the listed proteins",
              set(w["PeptideSequence"]) == {p for p, _ in ref_digest(p2, REG["T"], 7, 30, 1)}
              and set(w["ProteinGroup"]) == {"P2"})
        check("... and a peptide's label then comes from those proteins only",
              w.set_index("PeptideSequence").at["DDDDDDDR", "Specificity"] == "semi_specific")

        n0, recs0, m0 = vfasta(td, FASTA, "T", tag="z", protein_whitelist={"NOPE"})
        check("no surviving protein: 0 peptides, empty FASTA, map with header only",
              n0 == 0 and recs0 == [] and m0.empty and "Specificity" in m0.columns)

        _, _, s = vfasta(td, FASTA, "T", tag="s", min_len=8, max_len=12, missed=0)
        check("min_len / max_len / missed are passed through",
              set(s["PeptideSequence"]) ==
              {p for q in (p1, p2) for p, _ in ref_digest(q, REG["T"], 8, 12, 0)})

        _, _, g = vfasta(td, FASTA, "A", tag="g")
        check("enzyme code names the headers (AspN)",
              len(g) > 0 and g["VirtualHeader"].str.startswith("vpep_A_").all())


# ── modifications ───────────────────────────────────────────────────────────
MOD_FLAG_TOKENS = {"--fixed-mod", "--var-mod", "--var-mods",
                   "--monitor-mod", "--no-cut-after-mod"}


def import_server():
    """bin/server.py with the web framework stubbed if it is not installed.

    build_pipeline_cmd needs none of fastapi/uvicorn, but server.py imports them at
    module level. The workstation that runs the pipeline has no fastapi, so without
    this the suite stops before the GUI-to-CLI checks -- on the machine where the
    pipeline actually runs.
    """
    try:
        import server
        return server, True
    except ModuleNotFoundError as exc:
        if exc.name not in ("fastapi", "uvicorn"):
            raise
    from unittest.mock import MagicMock
    for name in ("fastapi", "fastapi.responses", "fastapi.staticfiles",
                 "fastapi.middleware", "fastapi.middleware.cors", "uvicorn"):
        sys.modules.setdefault(name, MagicMock())
    import server
    return server, False


def mod_args_of(cmd):
    """The modification arguments of one DIA-NN command, in order."""
    out, i = [], 0
    while i < len(cmd):
        if cmd[i] in MOD_FLAG_TOKENS:
            out += cmd[i:i + 2]; i += 2
        elif cmd[i] == "--met-excision":
            out.append(cmd[i]); i += 1
        else:
            i += 1
    return out


class FakeDiann:
    """Stands in for run_cmd: records each DIA-NN command, writes its outputs."""

    def __init__(self):
        self.calls = []

    def __call__(self, cmd, log_path, label="", diann_exe="diann"):
        cmd = [str(c) for c in cmd]
        self.calls.append((label, cmd))
        val = lambda flag: cmd[cmd.index(flag) + 1] if flag in cmd else None
        out_lib, report = val("--out-lib"), val("--out")
        if out_lib and "--predictor" in cmd:
            write(out_lib.replace(".tsv", ".predicted.speclib"), "x" * 200)
        elif out_lib:
            write(out_lib, "transition_group_id\tPeptideSequence\tProteinGroup\n"
                           + "".join(f"AAAK{i}2\tAAAK\tP1\n" for i in range(20)))
        if report and "--predictor" not in cmd:
            write(report, "Run\tProtein.Group\n" + "r\tP1\n" * 20)
        return " ".join(cmd), 0.0


def run_all_passes(td, mods):
    """Run passes 1-4 against FakeDiann; returns [(label, cmd)] for all six calls."""
    fake = FakeDiann()
    saved = rp.run_cmd, rp.run_plotter
    rp.run_cmd, rp.run_plotter = fake, (lambda *a, **k: None)
    try:
        fasta = os.path.join(td, "db.fasta")
        write(fasta, ">sp|P1|A\nMSTEKLLAGRPEPTIDECKSAMPLEQQWR\n")
        out = os.path.join(td, "out")
        os.makedirs(out, exist_ok=True)
        log = os.path.join(out, "pipeline.log")
        p1 = quiet(rp.run_pass1, ["a.raw"], out, ["T"], fasta, {}, None, 1, "diann",
                   log, {}, mods=mods)
        comb = os.path.join(out, "combined.tsv")
        quiet(rp.combine_libraries, {"Trypsin": p1["T"][0]}, comb)
        quiet(rp.run_pass2, ["a.raw"], out, ["T"], comb, fasta, 1, "diann", log, {}, mods=mods)
        quiet(rp.run_pass3, ["a.raw"], out, ["T"], p1, fasta, 1, "diann", log, {},
              mods=mods)
        quiet(rp.run_pass4, ["a.raw"], out, ["T"], comb, {}, fasta, 1, "diann", log, {},
              mods=mods)
    finally:
        rp.run_cmd, rp.run_plotter = saved
    return fake.calls


def test_mod_flags_all_passes():
    print("\nmodifications reach all six DIA-NN calls")
    mods = rp.make_mod_config(["carbamidomethyl", "oxidation", "acetyl-nterm"],
                              met_excision=True)
    with tempfile.TemporaryDirectory() as td:
        calls = run_all_passes(td, mods)
    labels = [l.split(" —")[0].split(" [")[0] for l, _ in calls]
    check("six DIA-NN calls: pass1 predict+search, pass2, pass3 predict+search, pass4",
          labels == ["Pass1 predict", "Pass1 search", "Pass2", "Pass3 predict virtual FASTA",
                     "Pass3 semi search (inline reannotate)", "Pass4"], labels)
    cam = ["--fixed-mod", "UniMod:4,57.021464,C"]
    ox  = ["--var-mod", "UniMod:35,15.994915,M"]
    for label, cmd in calls:
        a = mod_args_of(cmd)
        check(f"{label[:40]}: carbamidomethyl C fixed", cam == a[0:2], a)
        check(f"{label[:40]}: oxidation M variable + --var-mods 1",
              ox == a[2:4] and a[a.index("--var-mods") + 1] == "1", a)
    full = mod_args_of(calls[0][1])
    check("protein-level run: acetyl *n and --met-excision emitted",
          "UniMod:1,42.010565,*n" in full and "--met-excision" in full, full)
    for i in (0, 2, 5):
        check(f"{labels[i]} uses the full configuration", mod_args_of(calls[i][1]) == full)
    check("pass1: library-generating and searching calls agree",
          mod_args_of(calls[0][1]) == mod_args_of(calls[1][1]))
    v = mod_args_of(calls[3][1])
    check("pass3 (virtual FASTA): no protein-N-term mod, no --met-excision",
          "UniMod:1,42.010565,*n" not in v and "--met-excision" not in v, v)
    check("pass3: library-generating and searching calls agree",
          mod_args_of(calls[3][1]) == mod_args_of(calls[4][1]))
    with tempfile.TemporaryDirectory() as td:
        bare = run_all_passes(td, rp.make_mod_config([]))
    check("no modifications configured: no modification flag anywhere",
          all(not mod_args_of(c) for _, c in bare))


def test_mod_settings_roundtrip():
    print("\nmodification settings round-trip through CLI parsing")
    cfg = rp.make_mod_config(["carbamidomethyl", "oxidation"],
                             variable=["UniMod:21,79.966331,STY"], max_var_mods=2,
                             met_excision=True, extra=["--peak-center"])
    base = ["--raw", "a.raw", "--out-dir", "o"]
    args = rp.build_arg_parser().parse_args(base + rp.mod_config_to_cli(cfg))
    check("pipeline CLI: settings → arguments → settings unchanged",
          rp.mod_config_from_args(args) == cfg, rp.mod_config_from_args(args))
    default = rp.mod_config_from_args(rp.build_arg_parser().parse_args(base))
    check("CLI default is carbamidomethyl (C) fixed and nothing else",
          default["modifications"] == [rp.MOD_PRESETS["carbamidomethyl"]]
          and default["metExcision"] is False and default["extraDiannArgs"] == [], default)
    server, _real = import_server()
    job = {"rawFiles": ["a.raw"], "outDir": "o", **{k: cfg[k] for k in
           ("modifications", "maxVarMods", "metExcision", "extraDiannArgs")}}
    cmd = server.build_pipeline_cmd({}, job)
    gui = rp.mod_config_from_args(
        rp.build_arg_parser().parse_args(cmd[cmd.index("--raw"):]))
    check("GUI job → server command → pipeline settings unchanged", gui == cfg, gui)
    legacy = server.build_pipeline_cmd({}, {"rawFiles": ["a.raw"], "outDir": "o"})
    check("GUI job without the key (older GUI): pipeline default applies",
          "--mods" not in legacy and "--fixed-mod" not in legacy)

    # The two user-facing ways of switching modifications off (#30).
    off = server.build_pipeline_cmd({}, {"rawFiles": ["a.raw"], "outDir": "o",
                                         "modifications": []})
    check("GUI, carbamidomethyl unchecked: --mods none, no --fixed-mod/--var-mod",
          "--mods" in off and off[off.index("--mods") + 1] == "none"
          and "--fixed-mod" not in off and "--var-mod" not in off, off)
    gui_off = rp.mod_config_from_args(
        rp.build_arg_parser().parse_args(off[off.index("--raw"):]))
    check("GUI, carbamidomethyl unchecked: pipeline settings have no modifications",
          gui_off["modifications"] == [], gui_off)
    cli_off = rp.mod_config_from_args(
        rp.build_arg_parser().parse_args(base + ["--mods", "none"]))
    check("CLI --mods none: no modifications", cli_off["modifications"] == [], cli_off)
    for name, cfg_off in (("GUI unchecked", gui_off), ("CLI --mods none", cli_off)):
        with tempfile.TemporaryDirectory() as td:
            calls = run_all_passes(td, cfg_off)
        check(f"{name}: no modification flag in any of the {len(calls)} DIA-NN calls",
              len(calls) == 6 and all(not mod_args_of(c) for _, c in calls),
              [mod_args_of(c) for _, c in calls])


def test_var_mods_flag():
    print("\n--var-mods only with variable modifications")
    f = rp.mod_flags
    check("fixed only: no --var-mods",
          "--var-mods" not in f(rp.make_mod_config(["carbamidomethyl"])))
    check("variable present: --var-mods N",
          f(rp.make_mod_config(["oxidation"], max_var_mods=3))[-2:] == ["--var-mods", "3"])
    nterm = rp.make_mod_config(["acetyl-nterm"])
    check("protein N-term only, virtual FASTA: no --var-mod, no --var-mods",
          f(nterm, protein_level=False) == [], f(nterm, protein_level=False))
    check("protein N-term only, protein level: emitted with its monitoring",
          f(nterm) == ["--var-mod", "UniMod:1,42.010565,*n", "--monitor-mod", "UniMod:1",
                       "--var-mods", "1"], f(nterm))


def test_prebuilt_lib_guard():
    print("\nprebuilt / reused library guard (configuration read from pipeline.settings)")
    import json
    cam_cfg = rp.make_mod_config(["carbamidomethyl"])
    cam = rp.mod_flags(cam_cfg, include_extra=False)

    def settings(run_dir, data):
        write(os.path.join(run_dir, "pipeline.settings"), json.dumps(data))

    with tempfile.TemporaryDirectory() as td:
        run = os.path.join(td, "otherrun")            # a previous mE-DIA run
        os.makedirs(os.path.join(run, "pass1_T"))
        lib = os.path.join(run, "pass1_T", "trypsin.predicted.speclib")
        write(lib, "x" * 200)

        settings(run, rp.make_mod_config([]))
        try:
            rp.check_lib_mods(lib, cam, "Pass1 Trypsin prebuilt speclib")
            msg = ""
        except rp.ModConfigError as e:
            msg = str(e)
        check("mismatch refused", bool(msg))
        check("message names both configurations",
              "(no modifications)" in msg and "--fixed-mod UniMod:4,57.021464,C" in msg, msg)

        settings(run, cam_cfg)
        try:
            rp.check_lib_mods(lib, cam, "x"); ok = True
        except rp.ModConfigError:
            ok = False
        check("matching configuration accepted", ok)

        settings(run, {"enzymes": ["T"]})               # written before modification support
        try:
            rp.check_lib_mods(lib, cam, "x"); refused_old = False
        except rp.ModConfigError:
            refused_old = True
        check("a run from before modification support counts as no modifications",
              refused_old)

        # pass 3 builds with protein_level=False flags; the guard must recompute the same way
        ptm = rp.make_mod_config(["carbamidomethyl", "acetyl-nterm"])
        settings(run, ptm)
        try:
            rp.check_lib_mods(lib, rp.mod_flags(ptm, protein_level=False, include_extra=False),
                              "Pass3", protein_level=False); ok3 = True
        except rp.ModConfigError:
            ok3 = False
        check("pass-3 library accepted when recomputed with pass-3 flags", ok3)

        os.remove(os.path.join(run, "pipeline.settings"))
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            rp.check_lib_mods(lib, cam, "Pass1 Trypsin prebuilt speclib")
        check("no pipeline.settings: loud warning, not silence", "WARNING" in buf.getvalue()
              and "Cannot verify" in buf.getvalue())

        # wired into pass 1: a supplied speclib from a run with other modifications stops the run
        settings(run, rp.make_mod_config([]))
        fasta = os.path.join(td, "db.fasta"); write(fasta, ">sp|P1|A\nMSTEKLLAGR\n")
        saved = rp.run_cmd
        rp.run_cmd = FakeDiann()
        try:
            quiet(rp.run_pass1, ["a.raw"], os.path.join(td, "o"), ["T"], fasta,
                  {"T": lib}, None, 1, "diann", os.path.join(td, "log"), {},
                  mods=cam_cfg)
            stopped = False
        except rp.ModConfigError:
            stopped = True
        finally:
            rp.run_cmd = saved
        check("run_pass1 refuses a mismatching --speclib before any DIA-NN call", stopped)

        # a full run leaves no per-library record files in the results
        with tempfile.TemporaryDirectory() as td2:
            run_all_passes(td2, rp.make_mod_config(["carbamidomethyl", "phospho"]))
            sidecars = [f for _, _, fs in os.walk(td2) for f in fs if f.endswith(".mods.json")]
        check("no .mods.json files written anywhere in the output", not sidecars, sidecars)

        # resuming into a folder made before modification support
        out = os.path.join(td, "res"); os.makedirs(os.path.join(out, "pass1_T"))
        write(os.path.join(out, "pipeline.settings"), '{"enzymes": ["T"]}')
        try:
            rp.check_resume_mods(out, rp.make_mod_config(["carbamidomethyl"])); refused = False
        except rp.ModConfigError as e:
            refused = "(no modifications)" in str(e)
        check("resume into older results with new modifications refused", refused)
        try:
            rp.check_resume_mods(out, rp.make_mod_config([])); same = True
        except rp.ModConfigError:
            same = False
        check("resume with the original (none) configuration allowed", same)


def test_malformed_mods():
    print("\nmalformed modifications rejected")
    bad = {
        "mass not a number":   ("fixed", "UniMod:4,abc,C"),
        "mass zero":           ("fixed", "UniMod:4,0,C"),
        "mass infinite":       ("fixed", "UniMod:4,inf,C"),
        "unknown site B":      ("variable", "UniMod:35,15.994915,B"),
        "lower-case residue":  ("variable", "UniMod:35,15.994915,m"),
        "empty sites":         ("variable", "UniMod:35,15.994915,"),
        "missing field":       ("fixed", "UniMod:4,57.021464"),
        "space in name":       ("fixed", "Carbamido methyl,57.021464,C"),
    }
    for what, (kind, spec) in bad.items():
        try:
            rp.parse_mod_spec(spec, kind); err = ""
        except rp.ModConfigError as e:
            err = str(e)
        check(f"rejected: {what}", bool(err), err)
    for preset in (["nonsense"],):
        try:
            rp.make_mod_config(preset); err = ""
        except rp.ModConfigError as e:
            err = str(e)
        check("unknown preset rejected with the known ones listed",
              "carbamidomethyl" in err, err)
    ok = rp.parse_mod_spec("UniMod:21,79.966331,STY", "variable")
    check("valid custom modification accepted",
          ok == {"name": "UniMod:21", "mass": 79.966331, "sites": "STY", "type": "variable"})
    check("peptide and protein N-term sites accepted",
          rp.parse_mod_spec("UniMod:1,42.010565,*nK", "variable")["sites"] == "*nK")
    server, _real = import_server()
    try:
        server.build_pipeline_cmd({}, {"rawFiles": ["a"], "outDir": "o",
                                       "modifications": [{"name": "X", "mass": 1,
                                                          "sites": "C", "type": "weird"}]})
        err = False
    except ValueError:
        err = True
    check("server rejects a modification type other than fixed/variable", err)


# The flags the DIA-NN 1.8.1 GUI itself emitted for carbamidomethyl + oxidation +
# acetyl(protein N-term) + phospho(STY) + GlyGly(K), read off report.log.txt of
# results/20261010_diannFlagCheck/181_oxM_acN_phospho_kgg. --unimod4 is the GUI's
# shorthand for carbamidomethyl: a probe run on the same FASTA showed that --unimod4 and
# --fixed-mod UniMod:4,57.021464,C give byte-identical predicted libraries, so the two are
# interchangeable and the pipeline emits the explicit form.
GUI_REFERENCE_PTM_FLAGS = [
    "--fixed-mod", "UniMod:4,57.021464,C",
    "--var-mod", "UniMod:35,15.994915,M",
    "--var-mod", "UniMod:1,42.010565,*n", "--monitor-mod", "UniMod:1",
    "--var-mod", "UniMod:21,79.966331,STY", "--monitor-mod", "UniMod:21",
    "--var-mod", "UniMod:121,114.042927,K", "--monitor-mod", "UniMod:121",
    "--no-cut-after-mod", "UniMod:121",
    "--var-mods", "1", "--met-excision",
]


def test_ptm_monitoring():
    print("\nPTM monitoring: --monitor-mod and --no-cut-after-mod")
    cfg = rp.make_mod_config(
        ["carbamidomethyl", "oxidation", "acetyl-nterm", "phospho", "gg"],
        max_var_mods=1, met_excision=True)
    flags = rp.mod_flags(cfg)
    check("the GUI's PTM flag set is reproduced exactly",
          flags == GUI_REFERENCE_PTM_FLAGS, "" if flags == GUI_REFERENCE_PTM_FLAGS else flags)
    check("--monitor-mod follows the modification it belongs to",
          flags[flags.index("--monitor-mod") - 2] == "--var-mod"
          and flags[flags.index("--monitor-mod") - 1].startswith("UniMod:1,"))
    check("oxidation is not monitored, matching the GUI",
          "--monitor-mod UniMod:35" not in " ".join(flags))
    check("GlyGly blocks cleavage after the modified residue",
          "--no-cut-after-mod" in flags
          and flags[flags.index("--no-cut-after-mod") + 1] == "UniMod:121")

    # pass 3: the protein-N-term modification is dropped, so its monitor flag must go too,
    # or DIA-NN is told to monitor a modification it was never given
    p3 = rp.mod_flags(cfg, protein_level=False)
    check("pass3 drops acetyl(*n) and its --monitor-mod together",
          "UniMod:1,42.010565,*n" not in p3 and "UniMod:1" not in p3, p3)
    check("pass3 keeps the residue PTMs and their monitoring",
          p3.count("--monitor-mod") == 2 and "UniMod:21" in p3 and "UniMod:121" in p3)

    def monitored_vs_declared(args):
        mon = [args[i + 1] for i, t in enumerate(args) if t == "--monitor-mod"]
        dec = [v.split(",")[0] for i, v in enumerate(args)
               if i and args[i - 1] in ("--var-mod", "--fixed-mod")]
        return mon, dec
    mon, dec = monitored_vs_declared(p3)
    check("pass3 monitors nothing it did not declare",
          all(m in dec for m in mon), (mon, dec))

    with tempfile.TemporaryDirectory() as td:
        calls = run_all_passes(td, cfg)
    check("six DIA-NN calls", len(calls) == 6, len(calls))
    for i, (label, cmd) in enumerate(calls):
        mon, dec = monitored_vs_declared(mod_args_of(cmd))
        check(f"call {i + 1} ({label}): monitoring present and self-consistent",
              bool(mon) and all(m in dec for m in mon), (mon, dec))
    check("pass3: monitoring identical between library generation and search",
          mod_args_of(calls[3][1]) == mod_args_of(calls[4][1]))

    args = rp.build_arg_parser().parse_args(
        ["--raw", "a.raw", "--out-dir", "o"] + rp.mod_config_to_cli(cfg))
    check("settings with PTM monitoring round-trip through the CLI unchanged",
          rp.mod_config_from_args(args) == cfg, rp.mod_config_from_args(args))
    server, _real = import_server()
    job = {"rawFiles": ["a.raw"], "outDir": "o", "modifications": cfg["modifications"],
           "maxVarMods": 1, "metExcision": True, "extraDiannArgs": []}
    cmd = server.build_pipeline_cmd({}, job)
    gui = rp.mod_config_from_args(rp.build_arg_parser().parse_args(cmd[cmd.index("--raw"):]))
    check("GUI job → server command → settings keep the monitoring", gui == cfg, gui)

    for bad, why in (
            (dict(presets=["oxidation"], monitor=["UniMod:99"]),
             "monitoring a modification that is not configured"),
            (dict(presets=["carbamidomethyl"], monitor=["UniMod:4"]),
             "monitoring a fixed modification"),
            (dict(presets=["oxidation"], no_cut_after=["nope"]),
             "blocking cleavage for a modification that is not configured")):
        try:
            rp.make_mod_config(**bad)
            check(f"rejected: {why}", False, "no error raised")
        except rp.ModConfigError as e:
            check(f"rejected: {why}", True, str(e)[:60])


def test_pass3_empty_vfasta():
    print("\npass 3 with an empty virtual FASTA (#19)")
    real_gen = rp.generate_semi_vfasta

    def gen(fasta, code, enz, vfasta, vmap, **kw):
        if code == "E":                      # Elastase: no candidates
            write(vfasta, "")
            write(vmap, "VirtualHeader\tPeptideSequence\tSpecificity\n")
            return 0
        return real_gen(fasta, code, enz, vfasta, vmap, **kw)

    saved = rp.run_cmd, rp.run_plotter, rp.generate_semi_vfasta
    rp.run_cmd, rp.run_plotter = FakeDiann(), (lambda *a, **k: None)
    rp.generate_semi_vfasta = gen
    try:
        for order in (["E", "T"], ["T", "E"]):
            with tempfile.TemporaryDirectory() as td:
                fasta = os.path.join(td, "db.fasta")
                write(fasta, ">sp|P1|A\nMSTEKLLAGRPEPTIDECKSAMPLEQQWR\n")
                out = os.path.join(td, "out"); os.makedirs(out)
                try:
                    res = quiet(rp.run_pass3, ["a.raw"], out, order, {}, fasta, 1,
                                "diann", os.path.join(out, "log"), {})
                    err = ""
                except Exception as e:     # noqa: BLE001 - the bug was an UnboundLocalError
                    res, err = {}, f"{type(e).__name__}: {e}"
                tag = "first" if order[0] == "E" else "later"
                check(f"empty virtual FASTA as {tag} enzyme: no crash", not err, err)
                if err:
                    continue
                e_lib, t_lib = res["E"][0], res["T"][0]
                check(f"{tag}: the empty enzyme reports its own library path, "
                      f"not the other enzyme's",
                      e_lib != t_lib and "pass3_semiExt_E" in e_lib, (e_lib, t_lib))
                check(f"{tag}: that library does not exist, so pass 4 skips it",
                      not os.path.exists(e_lib) and os.path.exists(t_lib))
    finally:
        rp.run_cmd, rp.run_plotter, rp.generate_semi_vfasta = saved


def main():
    test_registry()
    test_cleavage_sites()
    test_digest()
    test_parse_fasta()
    test_vfasta()
    test_mod_flags_all_passes()
    test_mod_settings_roundtrip()
    test_ptm_monitoring()
    test_var_mods_flag()
    test_prebuilt_lib_guard()
    test_malformed_mods()
    test_pass3_empty_vfasta()
    print(f"\n{'all checks passed' if not FAIL else f'{len(FAIL)} FAILED: {FAIL}'}")
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
