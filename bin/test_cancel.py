#!/usr/bin/env python3
"""Checks for cancelling a run (#31). Run: python bin/test_cancel.py

Cancel must stop the pipeline *and* the DIA-NN it started, wait for both, and
report `cancelled`. DIA-NN is replaced by surrogates:
  - POSIX: a real process tree (pipeline -> long-running grandchild), so the
    process-group kill is exercised for real;
  - WSL and Windows: diann.exe is a Windows process that a POSIX signal cannot
    reach, simulated behind recording fakes for powershell.exe and taskkill.
    Only taskkill makes it go away, so a fix that killed the Linux side alone
    fails here, as it did on the workstation.
Also covers the restart guard: an interrupted DIA-NN call leaves no output a
restart would accept as complete.
"""
import asyncio
import base64
import contextlib
import io
import json
import os
import re
import signal
import subprocess
import sys
import tempfile
import time
import types

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import run_pipeline as rp          # noqa: E402
import server as sv                # noqa: E402

FAIL = []


def check(name, cond, detail=""):
    print(f"  {'PASS' if cond else 'FAIL'}  {name}{('  ' + str(detail)) if detail else ''}")
    if not cond:
        FAIL.append(name)
        # Without this a failing check only sets an exit code, which pytest never
        # sees, so the check would pass silently under `pytest bin/`.
        if "PYTEST_CURRENT_TEST" in os.environ:
            raise AssertionError(f"{name} {detail}")


@contextlib.contextmanager
def patched(obj, **attrs):
    saved = {k: getattr(obj, k) for k in attrs}
    for k, v in attrs.items():
        setattr(obj, k, v)
    try:
        yield
    finally:
        for k, v in saved.items():
            setattr(obj, k, v)


def alive(pid):
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    # a zombie is dead for our purposes
    try:
        with open(f"/proc/{pid}/stat") as f:
            return f.read().split(") ")[1][0] != "Z"
    except OSError:
        return True


# The pipeline stand-in: starts a long-running "DIA-NN" grandchild, records its
# PID, then waits on it, as run_pipeline.run_cmd does.
PIPELINE = r"""
import subprocess, sys, time
child = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(120)"])
open(sys.argv[1], "w").write(str(child.pid))
print("##STATUS## {\"pass\": \"pipeline\", \"status\": \"running\"}", flush=True)
print("##STATUS## {\"pass\": \"pass1_T\", \"status\": \"running\"}", flush=True)
print("##STATUS## {\"pass\": \"pass0_T\", \"status\": \"done\"}", flush=True)
child.wait()
"""


def wait_for(cond, timeout=10):
    t0 = time.time()
    while time.time() - t0 < timeout:
        if cond():
            return True
        time.sleep(0.05)
    return False


def test_posix_tree_cancel():
    print("\nPOSIX: cancel stops the pipeline and its DIA-NN child, then reports")
    if sys.platform == "win32":
        print("  (skipped on Windows)")
        return

    sent = []

    class WS:                                   # a GUI client: records every message
        async def send_text(self, m):
            sent.append(json.loads(m))
    ws = WS()

    async def scenario(td):
        pidfile = os.path.join(td, "child.pid")
        run = sv.RunState("t1", [sys.executable, "-c", PIPELINE, pidfile])
        run.clients.add(ws)
        sv.active_runs[run.run_id] = run
        run.task = asyncio.create_task(sv._run_pipeline(run))
        for _ in range(200):
            if os.path.exists(pidfile) and open(pidfile).read():
                break
            await asyncio.sleep(0.05)
        child = int(open(pidfile).read())
        before = alive(child)
        with patched(sv, _windows_run_pids=lambda cmd: []):
            resp = await sv.cancel_run(run.run_id)
        return run, child, before, resp

    with tempfile.TemporaryDirectory() as td:
        run, child, before, resp = asyncio.run(scenario(td))
    check("the DIA-NN stand-in was running before the cancel", before)
    check("after the response, the DIA-NN stand-in has exited", not alive(child), child)
    check("... and so has the pipeline", run.process.returncode is not None)
    check("response: ok, status cancelled", resp == {"ok": True, "status": "cancelled",
                                                    "stillRunning": []}, resp)
    check("run status is cancelled, not failed", run.status == "cancelled", run.status)
    passes = {k: v["status"] for k, v in run.pass_statuses.items()}
    check("passes that were running are cancelled, finished ones kept",
          passes == {"pipeline": "cancelled", "pass1_T": "cancelled", "pass0_T": "done"},
          passes)
    lines = [m["data"] for m in sent if m["type"] == "line"]
    check("... and the GUI is sent the same ##STATUS## lines",
          '##STATUS## {"pass": "pass1_T", "status": "cancelled"}' in lines
          and '##STATUS## {"pass": "pipeline", "status": "cancelled"}' in lines, lines)
    statuses = [m["data"]["status"] for m in sent if m["type"] == "status"]
    check("the GUI is only ever told 'cancelled'", statuses and
          set(statuses) == {"cancelled"}, statuses)


def test_posix_finished_run_not_relabelled():
    print("\nA run that ends on its own keeps its own status")
    if sys.platform == "win32":
        print("  (skipped on Windows)")
        return

    async def scenario():
        out = []
        for code in (0, 3):
            run = sv.RunState(f"e{code}", [sys.executable, "-c", f"raise SystemExit({code})"])
            await sv._run_pipeline(run)
            out.append(run.status)
        return out

    check("exit 0 -> done, exit 3 -> failed", asyncio.run(scenario()) == ["done", "failed"])


class FakeWindows:
    """powershell.exe and taskkill as seen from WSL or Windows.

    `procs` maps PID -> command line of live Windows processes. The PowerShell
    query is answered by actually filtering them the way the script does; only
    taskkill removes a process (unless it is in `unkillable`).
    """

    def __init__(self, procs, unkillable=()):
        self.procs, self.unkillable, self.calls = dict(procs), set(unkillable), []

    def __call__(self, cmd, **kw):
        self.calls.append(list(cmd))
        exe = os.path.basename(cmd[0]).lower()
        if exe.startswith("powershell"):
            script = base64.b64decode(cmd[cmd.index("-EncodedCommand") + 1]).decode("utf-16-le")
            lit = re.findall(r"FromBase64String\('([^']*)'\)", script)
            folder = base64.b64decode(lit[0]).decode()
            names = base64.b64decode(lit[1]).decode().split("|")
            pids = [pid for pid, (name, cl) in self.procs.items()
                    if name.lower() in names and folder in cl.replace("/", "\\").lower()]
            return types.SimpleNamespace(returncode=0, stdout="\r\n".join(map(str, pids)),
                                         stderr="")
        if exe.startswith("taskkill"):
            pid = int(cmd[cmd.index("/PID") + 1])
            if pid not in self.unkillable:
                self.procs.pop(pid, None)
            return types.SimpleNamespace(returncode=0, stdout="", stderr="")
        raise AssertionError(f"unexpected command {cmd}")


DIANN_CL = (r"C:\DIA-NN\1.8.1\diann.exe --fasta E:\HighContentDIA\results\run_T\pass3_semiExt_T"
            r"\T_semi_vfasta.fasta --fasta-search --out-lib E:\HighContentDIA\results\run_T"
            r"\pass3_semiExt_T\T_semi.tsv --threads 14")
OTHER = {
    5001: ("diann.exe", r"C:\DIA-NN\1.8.1\diann.exe --out E:\HighContentDIA\results\run_T2"
                        r"\pass2_T\report.tsv"),            # another run: sibling folder
    5002: ("notepad.exe", r"notepad.exe E:\HighContentDIA\results\run_T\pipeline.log"),
}


class DeadProc:
    """A pipeline process that has already been stopped (its Linux side is gone)."""
    pid, returncode = 4000, -9

    async def wait(self):
        return self.returncode


def wsl_run(cmd_out_dir="/mnt/e/HighContentDIA/results/run_T"):
    run = sv.RunState("w1", ["python3", "run_pipeline.py", "--out-dir", cmd_out_dir,
                             "--diann", "/mnt/c/DIA-NN/1.8.1/diann.exe"])
    run.process = DeadProc()
    return run


def fake_to_windows(p):
    m = re.match(r"^/mnt/([a-z])(/.*)?$", p)
    return f"{m.group(1).upper()}:" + (m.group(2) or "").replace("/", "\\")


def test_wsl_windows_diann_killed():
    print("\nWSL: the Windows diann.exe is found and killed with taskkill")
    fake = FakeWindows({4032: ("diann.exe", DIANN_CL), **OTHER})
    run = wsl_run()
    with patched(sv, _is_wsl=lambda: True), patched(sv.sys, platform="linux"), \
         patched(sv.subprocess, run=fake), patched(sv.shutil, which=lambda n: "/x/" + n), \
         patched(sv._wslpaths, to_windows=fake_to_windows):
        remaining = asyncio.run(sv._cancel_tree(run))
    kills = [c for c in fake.calls if os.path.basename(c[0]).startswith("taskkill")]
    check("taskkill.exe /T /F on the run's diann.exe", kills ==
          [["taskkill.exe", "/T", "/F", "/PID", "4032"]], kills)
    check("diann.exe is gone, nothing reported left", 4032 not in fake.procs and remaining == [],
          remaining)
    check("another run's DIA-NN (sibling folder) and a non-DIA-NN program are untouched",
          5001 in fake.procs and 5002 in fake.procs)
    check("checked again after the kill, before answering",
          [os.path.basename(c[0]) for c in fake.calls][-1] == "powershell.exe", fake.calls)


def test_posix_kill_alone_is_not_enough():
    print("\nWSL: without taskkill the Windows process survives (the observed bug)")
    fake = FakeWindows({4032: ("diann.exe", DIANN_CL)}, unkillable={4032})
    run = wsl_run()
    with patched(sv, _is_wsl=lambda: True, CANCEL_WAIT_S=1), patched(sv.sys, platform="linux"), \
         patched(sv.subprocess, run=fake), patched(sv.shutil, which=lambda n: "/x/" + n), \
         patched(sv._wslpaths, to_windows=fake_to_windows):
        remaining = asyncio.run(sv._cancel_tree(run))
    check("a process that will not die is reported, not hidden", remaining == [4032], remaining)


def test_windows_native():
    print("\nWindows Python: taskkill /T on the pipeline, then the folder sweep")
    fake = FakeWindows({4032: ("diann.exe", DIANN_CL)})
    proc = types.SimpleNamespace(pid=777, returncode=None)
    with patched(sv.subprocess, run=fake), patched(sv.sys, platform="win32"):
        sv._stop_pipeline(proc)
    check("pipeline tree: taskkill /T /F /PID <pipeline>",
          fake.calls == [["taskkill", "/T", "/F", "/PID", "777"]], fake.calls)
    with patched(sv, _is_wsl=lambda: False), patched(sv.sys, platform="win32"):
        d = sv._windows_run_dir(["x", "--out-dir", "E:/HighContentDIA/results/run_T/"])
    check("output folder in DIA-NN's form", d == "e:\\highcontentdia\\results\\run_t\\", d)
    with patched(sv, _is_wsl=lambda: True):
        d = sv._windows_run_dir(["x", "--out-dir", "/home/u/res"])
    check("WSL, Linux-filesystem folder: no Windows sweep", d is None, d)
    names = sv._tool_names(["x", "--diann", r"D:\tools\DiaNN-1.8.1.exe"])
    check("a renamed DIA-NN exe is matched too", "diann-1.8.1.exe" in names
          and "diann.exe" in names, names)


def test_cancel_response_shapes():
    print("\ncancel_run: response and status")

    async def scenario(remaining):
        run = wsl_run()
        sv.active_runs[run.run_id] = run

        async def fake_tree(r):
            return remaining
        with patched(sv, _cancel_tree=fake_tree):
            return await sv.cancel_run(run.run_id), run.status

    resp, status = asyncio.run(scenario([4032]))
    check("leftover process: ok false, PID listed", resp == {
        "ok": False, "status": "cancelled", "stillRunning": [4032]}, resp)
    check("... and the run is still marked cancelled", status == "cancelled")


# ─── pipeline: restart guard ─────────────────────────────────────────────────

FAKE_DIANN = r"""
import sys
a = sys.argv
out = a[a.index("--out") + 1]
open(out, "w").write("Run\tProtein.Group\n" + "r\tP1\n" * 50)   # > 100 bytes
raise SystemExit(int(a[a.index("--exit") + 1]))
"""


def test_inflight_marker_and_recovery():
    print("\nrestart guard: an interrupted DIA-NN call's output is set aside")
    with tempfile.TemporaryDirectory() as td:
        log = os.path.join(td, "pipeline.log")
        report = os.path.join(td, "pass2_T", "report.tsv")
        os.makedirs(os.path.dirname(report))
        marker = os.path.join(td, rp.INFLIGHT)
        cmd = [sys.executable, "-c", FAKE_DIANN, "--out", report]
        with contextlib.redirect_stdout(io.StringIO()):
            rp.run_cmd(cmd + ["--exit", "0"], log, "ok")
        check("clean exit: no in-flight record left", not os.path.exists(marker))
        os.remove(report)
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                rp.run_cmd(cmd + ["--exit", "1"], log, "Pass2")
        except RuntimeError:
            pass
        check("interrupted call: in-flight record stays", os.path.exists(marker))
        check("... and its partial output passes the old restart check",
              os.path.getsize(report) > 100)
        with contextlib.redirect_stdout(io.StringIO()) as out:
            moved = rp.recover_interrupted(td)
        check("next start: the output is set aside, not deleted",
              not os.path.exists(report) and os.path.exists(report + ".interrupted"), moved)
        check("... the record is cleared, and the notice names the step",
              not os.path.exists(marker) and "Pass2" in out.getvalue(), out.getvalue())
        check("nothing to do on a clean folder", rp.recover_interrupted(td) == [])
    lib = "/x/pass1_T/trypsin.tsv"
    check("targets: --out, --out-lib and the predicted speclib",
          rp.diann_targets(["diann", "--out", "/x/r.tsv", "--out-lib", lib, "--predictor"])
          == ["/x/r.tsv", lib, "/x/pass1_T/trypsin.predicted.speclib"])


def test_write_tsv_atomic():
    print("\nwrite_tsv: never a partial file under the final name")
    import pandas as pd

    class Boom(pd.DataFrame):
        def to_csv(self, path, **kw):
            open(path, "w").write("x" * 500)
            raise KeyboardInterrupt

    with tempfile.TemporaryDirectory() as td:
        p = os.path.join(td, "combined.tsv")
        try:
            rp.write_tsv(Boom({"a": [1]}), p)
        except KeyboardInterrupt:
            pass
        check("interrupted write: no file under the final name", not os.path.exists(p))
        rp.write_tsv(pd.DataFrame({"a": [1, 2]}), p)
        check("complete write: the file, and no .part left",
              open(p).read() == "a\n1\n2\n" and not os.path.exists(p + ".part"))


TESTS = (test_posix_tree_cancel, test_posix_finished_run_not_relabelled,
         test_wsl_windows_diann_killed, test_posix_kill_alone_is_not_enough,
         test_windows_native, test_cancel_response_shapes,
         test_inflight_marker_and_recovery, test_write_tsv_atomic)


def main():
    for t in TESTS:
        t()
    print(f"\n{'all checks passed' if not FAIL else f'{len(FAIL)} FAILED: {FAIL}'}")
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
