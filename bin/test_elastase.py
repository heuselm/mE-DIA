#!/usr/bin/env python3
"""Checks for the extended elastase EX and per-enzyme missed cleavages (#35).
Run: python bin/test_elastase.py

No DIA-NN needed: DIA-NN calls go to test_run_pipeline's FakeDiann, which
records each command and writes stand-in outputs.
"""
import contextlib
import io
import json
import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import run_pipeline as rp            # noqa: E402
import test_run_pipeline as trp      # noqa: E402  (FakeDiann, write, quiet, import_server)

REG = rp.ENZYME_REGISTRY
FAIL = []
CAM = rp.make_mod_config(["carbamidomethyl"])


def check(name, cond, detail=""):
    print(f"  {'PASS' if cond else 'FAIL'}  {name}{('  ' + str(detail)) if detail else ''}")
    if not cond:
        FAIL.append(name)
        # Without this a failing check only sets an exit code, which pytest never
        # sees, so the check would pass silently under `pytest bin/`.
        if "PYTEST_CURRENT_TEST" in os.environ:
            raise AssertionError(f"{name} {detail}")


def refused(fn, *a, **kw):
    """The ModConfigError message, or '' if fn accepted."""
    try:
        fn(*a, **kw)
        return ""
    except rp.ModConfigError as e:
        return str(e)


@contextlib.contextmanager
def fake_diann():
    fake = trp.FakeDiann()
    saved = rp.run_cmd, rp.run_plotter
    rp.run_cmd, rp.run_plotter = fake, (lambda *a, **k: None)
    try:
        yield fake
    finally:
        rp.run_cmd, rp.run_plotter = saved


@contextlib.contextmanager
def registry_restored():
    saved = dict(REG)
    try:
        yield
    finally:
        REG.clear()
        REG.update(saved)


def arg(cmd, flag):
    return cmd[cmd.index(flag) + 1] if flag in cmd else None


# ── tests ────────────────────────────────────────────────────────────────────
def test_registry_and_parsing():
    print("EX in the registry; E unchanged; code parsing")
    check("EX: A,V,S,G,L,I C-terminal, no proline rule, own mc 3",
          REG["EX"] == rp.EnzymeDef("EX", "Elastase (extended)", "A*,V*,S*,G*,L*,I*", 1,
                                    "AVSGLI", "", missed_cleavages=3), REG["EX"])
    check("E exactly as before: A*,V*, no own mc",
          REG["E"] == rp.EnzymeDef("E", "Elastase", "A*,V*", 1, "AV", ""), REG["E"])
    check("positional construction without the new field still works",
          rp.EnzymeDef("Q", "Q", "Q*", 1, "Q", "", True).missed_cleavages is None)
    for text, want in (("TCGEX", ["T", "C", "G", "EX"]), ("TCGE", ["T", "C", "G", "E"]),
                       ("T,C,G,EX", ["T", "C", "G", "EX"]), ("tcgex", ["T", "C", "G", "EX"]),
                       ("EEX", ["E", "EX"]), ("EXE", ["EX", "E"])):
        got = rp.parse_enzyme_codes(text)
        check(f"parse {text!r} -> {','.join(want)}", got == (want, []), got)
    check("file stems: elastase unchanged, EX filename-safe",
          REG["E"].stem == "elastase" and REG["EX"].stem == "elastase_extended",
          (REG["E"].stem, REG["EX"].stem))
    changed = [c for c, e in REG.items() if e.stem != e.name.lower() and c != "EX"]
    check("no existing enzyme's file names change except T2 ('trypsin/p' was a path bug)",
          changed == ["T2"], changed)


def test_effective_mc():
    print("\neffective missed cleavages: override > enzyme default > global")
    r = rp.resolve_enzyme_mc
    check("defaults: EX 3, others the global value",
          r(["T", "E", "EX"], 1) == {"T": 1, "E": 1, "EX": 3})
    check("global 2 reaches enzymes without their own, not EX",
          r(["T", "EX"], 2) == {"T": 2, "EX": 3})
    check("--enzyme-mc EX=2 overrides EX's default", r(["EX"], 1, {"EX": 2}) == {"EX": 2})
    check("--enzyme-mc T=0 overrides the global value", r(["T", "EX"], 1, {"T": 0})
          == {"T": 0, "EX": 3})
    check("parse: repeatable, case-insensitive",
          rp.parse_enzyme_mc(["ex=2", "T=0"]) == {"EX": 2, "T": 0})
    for bad in ("ZZ=1", "EX=-1", "EX=two", "EX", "=3"):
        try:
            rp.parse_enzyme_mc([bad]); err = ""
        except ValueError as e:
            err = str(e)
        check(f"rejected: --enzyme-mc {bad}", bool(err), err)
    with tempfile.TemporaryDirectory() as td:
        buf = io.StringIO()
        with contextlib.redirect_stderr(buf), contextlib.redirect_stdout(io.StringIO()):
            try:
                rp.main(["--raw", "a.raw", "--out-dir", td, "--enzymes", "EX",
                         "--skip-pass1", "--enzyme-mc", "EX=-1"]); code = 0
            except SystemExit as e:
                code = e.code
        check("main(): a negative override stops the run with a usage error",
              code == 2 and "0 or more" in buf.getvalue(), buf.getvalue().strip().splitlines()[-1:])


def test_pass1_commands():
    print("\npass 1: each enzyme's prediction gets its own rule and limit")
    with tempfile.TemporaryDirectory() as td, fake_diann() as fake:
        fasta = os.path.join(td, "db.fasta")
        trp.write(fasta, ">sp|P1|A\nMSTEKLLAGRPEPTIDECKSAMPLEQQWR\n")
        out = os.path.join(td, "out"); os.makedirs(out)
        mc = rp.resolve_enzyme_mc(["E", "EX"], 1)
        trp.quiet(rp.run_pass1, ["a.raw"], out, ["E", "EX"], fasta, {}, None, 1, "diann",
                  os.path.join(out, "pipeline.log"), {}, mods=CAM, enzyme_mc=mc)
        calls = {lab: cmd for lab, cmd in fake.calls}
    pred = {c: next(cmd for lab, cmd in calls.items() if lab.startswith("Pass1 predict")
                    and f"— {REG[c].name} [" in lab) for c in ("E", "EX")}
    check("EX predict: --cut A*,V*,S*,G*,L*,I* --missed-cleavages 3",
          (arg(pred["EX"], "--cut"), arg(pred["EX"], "--missed-cleavages"))
          == ("A*,V*,S*,G*,L*,I*", "3"), pred["EX"])
    p1 = os.path.join(out, "pass1_E")
    golden_e = ["diann", "--fasta", fasta, "--fasta-search", "--cut", "A*,V*",
                "--missed-cleavages", "1", "--min-pep-len", "7", "--max-pep-len", "30",
                "--fixed-mod", "UniMod:4,57.021464,C", "--predictor", "--gen-spec-lib",
                "--out-lib", os.path.join(p1, "elastase_predicted.tsv"),
                "--out", os.path.join(p1, "report_predict.tsv"), "--threads", "1"]
    check("E predict: byte for byte the command main produced", pred["E"] == golden_e,
          pred["E"])
    e_search = calls["Pass1 search — Elastase"]
    check("E search: library and outputs under the unchanged names",
          arg(e_search, "--lib").endswith(os.path.join("pass1_E",
                                                       "elastase_predicted.predicted.speclib"))
          and arg(e_search, "--out-lib").endswith(os.path.join("pass1_E",
                                                               "elastase_empirical.tsv")),
          e_search)
    ex_paths = [c for cmd in calls.values() for c in cmd if "pass1_EX" in c]
    check("EX paths are filename-safe (no space or parenthesis)",
          ex_paths and not any(ch in p for p in ex_paths for ch in " ()"), ex_paths)


# P1: EX segments ending S, G, L, I, A (no A/V inside): a span with 3 internal
#     S/G/L/I sites is in EX's limit-3 digest, one with 4 is not.
# P2: E segments ending V, A, V (no S/G/L/I): 1 internal A/V in, 2 out at mc 1.
SEG = ["DDDDS", "EEEEG", "KKKKL", "RRRRI", "TTTTA"]
P1 = "MTTTTA" + "".join(SEG) + "TTTT"
P2 = "MKKKKA" + "DDDDV" + "EEEEA" + "KKKKV" + "TTTT"


def test_pass3_vfasta_limit():
    print("\npass 3: the virtual FASTA is digested with each enzyme's limit")
    import pandas as pd
    with tempfile.TemporaryDirectory() as td, fake_diann():
        fasta = os.path.join(td, "db.fasta")
        trp.write(fasta, f">sp|P1|A\n{P1}\n>sp|P2|B\n{P2}\n")
        out = os.path.join(td, "out"); os.makedirs(out)
        mc = rp.resolve_enzyme_mc(["E", "EX"], 1)
        trp.quiet(rp.run_pass3, ["a.raw"], out, ["E", "EX"], {}, fasta, 1, "diann",
                  os.path.join(out, "pipeline.log"), {}, mods=CAM, enzyme_mc=mc)
        maps = {c: pd.read_csv(os.path.join(out, f"pass3_semiExt_{c}",
                                            f"{c}_semi_vfasta_map.tsv"), sep="\t")
                for c in ("E", "EX")}
    spec = {c: dict(zip(m["PeptideSequence"], m["Specificity"])) for c, m in maps.items()}
    three, four = "".join(SEG[:4]), "".join(SEG)
    check("EX: 3 internal S/G/L/I sites -> in, fully specific",
          spec["EX"].get(three) == "fully_specific", spec["EX"].get(three))
    check("EX: 4 internal sites -> absent", four not in spec["EX"])
    one, two = "DDDDVEEEEA", "DDDDVEEEEAKKKKV"
    check("E (global 1): 1 internal A/V -> in, fully specific",
          spec["E"].get(one) == "fully_specific", spec["E"].get(one))
    check("E (global 1): 2 internal A/V -> absent", two not in spec["E"])
    check("E did not get EX's limit: EX's 3-site span is no E peptide either",
          spec["E"].get(three) != "fully_specific")


def test_settings_provenance():
    print("\npipeline.settings records the definitions used")
    with tempfile.TemporaryDirectory() as td, registry_restored():
        cfg = os.path.join(td, "enz.json")
        trp.write(cfg, json.dumps({"Q": {"code": "Q", "name": "Qase", "cut": "Q*",
                                         "sense": 1, "cut_res": "Q", "noclv": ""}}))
        out = os.path.join(td, "o")
        trp.quiet(rp.main, ["--raw", "a.raw", "--out-dir", out, "--enzymes", "TEXQ",
                            "--enzyme-config", cfg, "--missed-cleavages", "2",
                            "--enzyme-mc", "T=0", "--skip-pass1"])
        with open(os.path.join(out, "pipeline.settings")) as f:
            s = json.load(f)
    defs = {d["code"]: d for d in s.get("enzymeDefinitions", [])}
    check("one record per enzyme of the run", sorted(defs) == ["EX", "Q", "T"], sorted(defs))
    check("EX: full definition, its own default 3",
          defs.get("EX") == {"code": "EX", "name": "Elastase (extended)",
                             "cut": "A*,V*,S*,G*,L*,I*", "sense": 1, "cutResidues": "AVSGLI",
                             "blocking": "", "missedCleavages": 3}, defs.get("EX"))
    check("T: the override (0), not the global value",
          defs.get("T", {}).get("missedCleavages") == 0, defs.get("T"))
    check("an --enzyme-config enzyme is recorded too, with the global value",
          defs.get("Q", {}).get("cut") == "Q*" and defs["Q"]["missedCleavages"] == 2,
          defs.get("Q"))
    check("overrides and the global value recorded as given",
          s.get("enzymeMissedCleavages") == {"T": 0} and s.get("missedCleavages") == 2)


def _folder(td, name, settings, passdir):
    d = os.path.join(td, name)
    os.makedirs(os.path.join(d, passdir))
    trp.write(os.path.join(d, "pipeline.settings"), json.dumps(settings))
    return d


def test_resume_guard():
    print("\nresume guard: same enzymes and definitions, or refused")
    e1, ex3 = rp.enzyme_spec("E", 1), rp.enzyme_spec("EX", 3)
    with tempfile.TemporaryDirectory() as td:
        old = _folder(td, "pre", {"enzymes": ["E"], "missedCleavages": 1, **CAM}, "pass1_E")
        check("pre-change E folder resumed as E: accepted",
              refused(rp.check_resume_enzymes, old, [e1]) == "")
        msg = refused(rp.check_resume_enzymes, old, [ex3])
        check("pre-change E folder resumed as EX: refused", bool(msg))
        check("... the message names both definitions",
              "A*,V*,S*,G*,L*,I*" in msg and "cut=A*,V*," in msg and "EX" in msg, msg)
        msg2 = refused(rp.check_resume_enzymes, old, [rp.enzyme_spec("E", 2)])
        check("pre-change E folder resumed with mc 2: refused",
              "missed cleavages=1" in msg2 and "missed cleavages=2" in msg2, msg2)
        new = _folder(td, "ex", {"enzymes": ["EX"], "enzymeDefinitions": [ex3]}, "pass1_EX")
        check("EX folder resumed as EX, mc 3: accepted",
              refused(rp.check_resume_enzymes, new, [ex3]) == "")
        msg3 = refused(rp.check_resume_enzymes, new, [rp.enzyme_spec("EX", 2)])
        check("EX folder resumed with mc 2: refused, both values named",
              "missed cleavages=3" in msg3 and "missed cleavages=2" in msg3, msg3)
        fresh = os.path.join(td, "fresh"); os.makedirs(fresh)
        trp.write(os.path.join(fresh, "pipeline.settings"), json.dumps({"enzymes": ["E"]}))
        check("a folder with settings but no results yet: nothing to protect",
              refused(rp.check_resume_enzymes, fresh, [ex3]) == "")
        buf = io.StringIO()
        with contextlib.redirect_stderr(buf), contextlib.redirect_stdout(io.StringIO()):
            try:
                rp.main(["--raw", "a.raw", "--out-dir", old, "--enzymes", "EX",
                         "--skip-pass1", "--mods", "carbamidomethyl"]); code = 0
            except SystemExit as e:
                code = e.code
        check("main(): resuming the E folder as EX stops before anything runs",
              code == 2 and "Elastase (extended)" in buf.getvalue(), buf.getvalue().strip().splitlines()[-1:])


def test_library_guard():
    print("\nlibrary guard: a library built under another enzyme definition is refused")
    lib_flags = rp.mod_flags(CAM, include_extra=False)
    e1, ex3 = rp.enzyme_spec("E", 1), rp.enzyme_spec("EX", 3)
    with tempfile.TemporaryDirectory() as td:
        run = _folder(td, "erun", {"enzymes": ["E"], "missedCleavages": 1, **CAM}, "pass1_E")
        lib = os.path.join(run, "pass1_E", "elastase_predicted.predicted.speclib")
        trp.write(lib, "x" * 200)
        msg = refused(rp.check_lib_mods, lib, lib_flags, "Pass1 EX prebuilt speclib",
                      enzymes=[ex3])
        check("E-built library supplied for EX: refused, both definitions named",
              "cut=A*,V*," in msg and "A*,V*,S*,G*,L*,I*" in msg, msg)
        check("... and accepted for E",
              refused(rp.check_lib_mods, lib, lib_flags, "x", enzymes=[e1]) == "")
        exrun = _folder(td, "exrun", {"enzymes": ["EX"], "enzymeDefinitions": [ex3], **CAM},
                        "pass1_EX")
        exlib = os.path.join(exrun, "pass1_EX", "elastase_extended_predicted.predicted.speclib")
        trp.write(exlib, "x" * 200)
        check("EX library with another mc: refused",
              "missed cleavages=3" in refused(rp.check_lib_mods, exlib, lib_flags, "x",
                                              enzymes=[rp.enzyme_spec("EX", 2)]))
        comb = _folder(td, "tcge", {"enzymes": ["T", "E"], "missedCleavages": 1, **CAM},
                       "pass2_TE")
        clib = os.path.join(comb, "pass2_TE", "combined_TE.tsv"); trp.write(clib, "x" * 200)
        check("combined library (--lib) of a T+E run, used for T+EX: refused",
              bool(refused(rp.check_lib_mods, clib, lib_flags, "--lib",
                           enzymes=[rp.enzyme_spec("T", 1), ex3])))
        check("... and accepted for T+E",
              refused(rp.check_lib_mods, clib, lib_flags, "--lib",
                      enzymes=[rp.enzyme_spec("T", 1), e1]) == "")
        fasta = os.path.join(td, "db.fasta"); trp.write(fasta, f">sp|P1|A\n{P1}\n")
        with fake_diann() as fake:
            msg = refused(trp.quiet, rp.run_pass1, ["a.raw"], os.path.join(td, "o"), ["EX"],
                          fasta, {"EX": lib}, None, 1, "diann", os.path.join(td, "log"), {},
                          mods=CAM, enzyme_mc={"EX": 3})
        check("run_pass1 refuses --speclib EX=<E library> before any DIA-NN call",
              bool(msg) and not fake.calls, msg)


def test_find_speclib():
    print("\n--speclib-dir: E does not take EX's library")
    with tempfile.TemporaryDirectory() as td:
        exf = os.path.join(td, "elastase_extended_predicted.predicted.speclib")
        trp.write(exf, "x" * 200)
        check("only an EX library present: none for E",
              rp.find_speclib({}, "E", "Elastase", td) is None)
        check("... and it is found for EX",
              rp.find_speclib({}, "EX", "Elastase (extended)", td) == exf)
        ef = os.path.join(td, "elastase_predicted.predicted.speclib"); trp.write(ef, "x" * 200)
        check("both present: E gets its own", rp.find_speclib({}, "E", "Elastase", td) == ef)


def test_proteoform_rules_and_viewer():
    print("\nproteoform rules and viewer colour know EX")
    try:
        import proteoform_analysis as pf
    except Exception as e:                       # matplotlib missing
        check("proteoform_analysis importable", False, e)
        return
    for c in ("E", "EX"):
        e = REG[c]
        want = {"sense": e.sense, "cut": e.cut_res, "block": e.noclv}
        check(f"{c}: rule under its name ({e.name!r}) matches the registry",
              pf.ENZYME_RULES.get(e.name) == want, pf.ENZYME_RULES.get(e.name))
        check(f"{c}: rule under its code matches the registry",
              pf.ENZYME_RULES.get(c) == want, pf.ENZYME_RULES.get(c))
    seq = "MKKKKSPEPTIDELGG"       # peptide starting after S (position 7): EX-specific only
    check("an N-terminus after S is specific for EX, not for E",
          pf.nterm_specific(seq, 7, "Elastase (extended)")
          and not pf.nterm_specific(seq, 7, "Elastase"))
    viewer = os.path.join(os.path.dirname(os.path.abspath(__file__)), "viewer", "viewer.html")
    with open(viewer, encoding="utf-8") as f:
        check("viewer has a colour for 'Elastase (extended)'",
              '"Elastase (extended)":' in f.read())


def test_gui_roundtrip():
    print("\nGUI job -> server command -> CLI -> pipeline.settings")
    server, _real = trp.import_server()
    job = {"rawFiles": ["a.raw"], "outDir": "o", "enzymes": ["T", "EX"],
           "missedCleavages": 1, "enzymeMissedCleavages": {"EX": 2, "G": 5}}
    cmd = server.build_pipeline_cmd({}, job)
    check("server: --enzyme-mc EX=2 for the selected enzyme only",
          [cmd[i + 1] for i, c in enumerate(cmd) if c == "--enzyme-mc"] == ["EX=2"], cmd)
    plain = server.build_pipeline_cmd({}, {**job, "enzymeMissedCleavages": {}})
    check("no overrides (or an older GUI): no --enzyme-mc", "--enzyme-mc" not in plain)
    legacy = {k: v for k, v in job.items() if k != "enzymeMissedCleavages"}
    check("job without the key: no --enzyme-mc",
          "--enzyme-mc" not in server.build_pipeline_cmd({}, legacy))
    try:
        server.build_pipeline_cmd({}, {**job, "enzymeMissedCleavages": {"EX": -1}}); err = ""
    except ValueError as e:
        err = str(e)
    check("server: a negative override is rejected", "0 or more" in err, err)
    with tempfile.TemporaryDirectory() as td:
        argv = cmd[cmd.index("--raw"):]
        argv[argv.index("--out-dir") + 1] = os.path.join(td, "o")
        trp.quiet(rp.main, argv + ["--skip-pass1"])
        with open(os.path.join(td, "o", "pipeline.settings")) as f:
            s = json.load(f)
    defs = {d["code"]: d["missedCleavages"] for d in s["enzymeDefinitions"]}
    check("pipeline: EX ran with the override 2, T with the global 1",
          defs == {"T": 1, "EX": 2}, defs)
    check("settings carry the override back for the GUI to load",
          s.get("enzymeMissedCleavages") == {"EX": 2}, s.get("enzymeMissedCleavages"))


TESTS = (test_registry_and_parsing, test_effective_mc, test_pass1_commands,
         test_pass3_vfasta_limit, test_settings_provenance, test_resume_guard,
         test_library_guard, test_find_speclib, test_proteoform_rules_and_viewer,
         test_gui_roundtrip)


def main():
    for t in TESTS:
        t()
    print(f"\n{'all checks passed' if not FAIL else f'{len(FAIL)} FAILED: {FAIL}'}")
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
