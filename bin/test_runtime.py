#!/usr/bin/env python3
"""Checks for the runtime record in pipeline.settings (#38).
Run: python bin/test_runtime.py

No DIA-NN needed: a Python one-liner stands in for it and prints DIA-NN 1.8.1's
banner (or nothing, for the unreadable case).
"""
import contextlib
import io
import json
import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import build_info as bi              # noqa: E402
import run_pipeline as rp            # noqa: E402

FAIL = []
BANNER = ["DIA-NN 1.8.1 (Data-Independent Acquisition by Neural Networks)",
          "Compiled on Apr 14 2022 15:31:19",
          "Current date and time: Sun Oct 12 09:00:00 2026"]
BANNER_ONE = BANNER[0] + " " + BANNER[1]


def check(name, cond, detail=""):
    print(f"  {'PASS' if cond else 'FAIL'}  {name}{('  ' + str(detail)) if detail else ''}")
    if not cond:
        FAIL.append(name)
        # Without this a failing check only sets an exit code, which pytest never
        # sees, so the check would pass silently under `pytest bin/`.
        if "PYTEST_CURRENT_TEST" in os.environ:
            raise AssertionError(f"{name} {detail}")


def fake_diann(lines, exit_code=0):
    """A command run_cmd treats as DIA-NN: 'diann' is replaced by this Python."""
    code = "import sys\n" + "".join(f"print({l!r})\n" for l in lines) + f"sys.exit({exit_code})"
    return ["diann", "-c", code]


def run_quiet(cmd, log):
    with contextlib.redirect_stdout(io.StringIO()) as out:
        try:
            rp.run_cmd(cmd, log, "test", sys.executable)
        except RuntimeError:
            pass
    return out.getvalue()


def test_runtime_block():
    print("the runtime block: present, well formed, what is loaded")
    rt = bi.runtime(r"C:\DIA-NN\1.8.1\diann.exe")
    check("keys", set(rt) == {"python", "python_impl", "platform", "executable", "frozen",
                              "packages", "diann"}, sorted(rt))
    check("python is this interpreter's X.Y.Z",
          rt["python"] == "%d.%d.%d" % tuple(sys.version_info[:3]), rt["python"])
    check("executable is sys.executable", rt["executable"] == sys.executable)
    import pandas, numpy
    check("pandas and numpy: the loaded versions",
          rt["packages"].get("pandas") == pandas.__version__
          and rt["packages"].get("numpy") == numpy.__version__, rt["packages"])
    check("every third-party package the pipeline imports is listed",
          set(rt["packages"]) >= {"pandas", "numpy", "matplotlib"}, rt["packages"])
    check("DIA-NN: the exe as given, banner unknown until DIA-NN runs",
          rt["diann"] == {"exe": r"C:\DIA-NN\1.8.1\diann.exe", "banner": "unknown"})
    check("frozen flag passed through", bi.runtime(frozen=True)["frozen"] is True)
    check("JSON-serialisable", bool(json.dumps(rt)))
    check("one-line summary names Python and pandas",
          "Python " in bi.runtime_one_line(rt) and "pandas " in bi.runtime_one_line(rt),
          bi.runtime_one_line(rt))


def test_unknowns_never_raise():
    print("\nwhat cannot be determined is 'unknown', never an exception")
    rt = bi.runtime(packages=("pandas", "no_such_package_xyz"))
    check("missing package -> unknown", rt["packages"]["no_such_package_xyz"] == "unknown")
    check("no DIA-NN exe -> unknown", bi.runtime()["diann"]["exe"] == "unknown")
    for lines, label in ((BANNER, "full banner"), (BANNER[:1], "first line only"),
                         ([], "nothing printed"), (["ERROR: cannot open"], "no banner"),
                         (None, "not a list")):
        try:
            got, err = bi.parse_diann_banner(lines), ""
        except Exception as e:      # noqa: BLE001
            got, err = None, repr(e)
        want = {"full banner": BANNER_ONE, "first line only": BANNER[0]}.get(label, "unknown")
        check(f"banner from {label}: {want!r}", got == want and not err, got or err)
    check("one-line summary of garbage -> unknown", bi.runtime_one_line(None) == "unknown")


def _settings(td):
    data = {"runtime": bi.runtime("diann")}
    with open(os.path.join(td, "pipeline.settings"), "w") as f:
        json.dump(data, f)


def _banner_in(td):
    with open(os.path.join(td, "pipeline.settings")) as f:
        return json.load(f)["runtime"]["diann"]["banner"]


def test_banner_from_first_diann_call():
    print("\nthe DIA-NN banner is taken from the first DIA-NN call's own output")
    with tempfile.TemporaryDirectory() as td:
        _settings(td)
        log = os.path.join(td, "pipeline.log")
        rp._DIANN_BANNER = None
        run_quiet([sys.executable, "-c", f"print({BANNER[0]!r})"], log)
        check("a non-DIA-NN command is not taken for DIA-NN", _banner_in(td) == "unknown")
        out = run_quiet(fake_diann(BANNER + ["Loading spectral library ..."]), log)
        check("pipeline.settings: runtime.diann.banner filled in", _banner_in(td) == BANNER_ONE,
              _banner_in(td))
        check("... and printed to the log", f"[runtime] DIA-NN: {BANNER_ONE}" in out, out[-300:])
        with open(log) as f:
            check("DIA-NN's own banner lines are still in pipeline.log",
                  BANNER[0] in f.read())
        run_quiet(fake_diann(["DIA-NN 9.9 (something else)"]), log)
        check("only the first DIA-NN call is recorded", _banner_in(td) == BANNER_ONE)
    with tempfile.TemporaryDirectory() as td:
        _settings(td)
        rp._DIANN_BANNER = None
        out = run_quiet(fake_diann(["garbled", "", "no version here"], exit_code=1),
                        os.path.join(td, "pipeline.log"))
        check("unreadable banner (and a failing call): stays unknown, no exception",
              _banner_in(td) == "unknown" and "Traceback" not in out, _banner_in(td))
    with tempfile.TemporaryDirectory() as td:
        with open(os.path.join(td, "pipeline.settings"), "w") as f:
            f.write("{ not json")
        rp._DIANN_BANNER = None
        out = run_quiet(fake_diann(BANNER), os.path.join(td, "pipeline.log"))
        check("a corrupt pipeline.settings: a warning, the call still completes",
              "could not record the DIA-NN banner" in out, out[-200:])
    rp._DIANN_BANNER = None


def test_settings_and_header():
    print("\na fresh run: runtime block in pipeline.settings, summary in the log header")
    with tempfile.TemporaryDirectory() as td:
        out_dir = os.path.join(td, "o")
        with contextlib.redirect_stdout(io.StringIO()) as out:
            rp.main(["--raw", "a.raw", "--out-dir", out_dir, "--enzymes", "T",
                     "--diann", r"C:\DIA-NN\1.8.1\diann.exe", "--skip-pass1"])
        with open(os.path.join(out_dir, "pipeline.settings")) as f:
            s = json.load(f)
    rt = s.get("runtime", {})
    check("runtime.python, packages.pandas, packages.numpy, diann.banner present",
          all(k for k in (rt.get("python"), rt.get("packages", {}).get("pandas"),
                          rt.get("packages", {}).get("numpy"), rt.get("diann", {}).get("banner"))),
          rt)
    # The exe the run actually used. On a WSL host normalise_cli_paths (#29) turns the
    # Windows-form --diann value into /mnt/c/..., so the literal argument is the wrong
    # expectation there; diannExe in the same settings file is what the run executed.
    check("runtime.diann.exe is the exe the run used (settings diannExe)",
          rt.get("diann", {}).get("exe") == s.get("diannExe")
          and str(s.get("diannExe", "")).replace("\\", "/").lower().endswith("dia-nn/1.8.1/diann.exe"),
          (rt.get("diann"), s.get("diannExe")))
    check("buildInfo still there, next to it", "buildInfo" in s)
    head = out.getvalue()
    check("log header: Runtime line with Python and pandas",
          f"Runtime : Python {rt.get('python')}" in head and "pandas " in head, head[-600:])
    check("log header: DIA-NN line with the exe the run used",
          f"DIA-NN  : {s.get('diannExe')}" in head, s.get("diannExe"))


TESTS = (test_runtime_block, test_unknowns_never_raise, test_banner_from_first_diann_call,
         test_settings_and_header)


def main():
    for t in TESTS:
        t()
    print(f"\n{'all checks passed' if not FAIL else f'{len(FAIL)} FAILED: {FAIL}'}")
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
