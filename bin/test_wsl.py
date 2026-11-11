#!/usr/bin/env python3
"""Checks for running from WSL (#25). Run: python bin/test_wsl.py

No WSL, PowerShell, display or DIA-NN needed: wslpath, subprocess and the
platform checks are replaced by recording fakes. Covers the Windows file
dialog, the path boundary to Windows executables, path normalisation in the
server and the CLI (including the #18 library guard), and start.py's npm calls.
"""
import base64
import contextlib
import io
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import types

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wslpaths                    # noqa: E402
import run_pipeline as rp          # noqa: E402
import server as sv                # noqa: E402
import start                       # noqa: E402
from fastapi import HTTPException  # noqa: E402

FAIL = []


def check(name, cond, detail=""):
    print(f"  {'PASS' if cond else 'FAIL'}  {name}{('  ' + str(detail)) if detail else ''}")
    if not cond:
        FAIL.append(name)
        # Without this a failing check only sets an exit code, which pytest never
        # sees, so the check would pass silently under `pytest bin/`.
        if "PYTEST_CURRENT_TEST" in os.environ:
            raise AssertionError(f"{name} {detail}")


def fake_wslpath(flag, path):
    """wslpath for the default automount root: /mnt/e/x <-> E:\\x."""
    if flag == "-w":
        m = re.match(r"^/mnt/([a-z])(/.*)?$", path)
        assert m, path
        return f"{m.group(1).upper()}:" + (m.group(2) or "/").replace("/", "\\")
    m = re.match(r"^([A-Za-z]):[\\/](.*)$", path)
    assert m, path
    return f"/mnt/{m.group(1).lower()}/" + m.group(2).replace("\\", "/")


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


class FakeRun:
    """Stands in for subprocess.run: records argv, returns a canned result."""

    def __init__(self, stdout="", returncode=0, raise_exc=None):
        self.stdout, self.rc, self.raise_exc, self.calls = stdout, returncode, raise_exc, []

    def __call__(self, cmd, **kw):
        self.calls.append((cmd, kw))
        if self.raise_exc:
            raise self.raise_exc
        return subprocess.CompletedProcess(cmd, self.rc, self.stdout, "boom")


def dialog_on_wsl(fake, mode="file", title="Select file", filetypes=None,
                  powershell="/mnt/c/Windows/powershell.exe"):
    """Call the dialog choke point as if on WSL; returns (result, HTTPException)."""
    with patched(sv, _is_wsl=lambda: True), \
         patched(sv.subprocess, run=fake), \
         patched(wslpaths, wslpath=fake_wslpath), \
         patched(shutil, which=lambda n: powershell if n == "powershell.exe" else None):
        try:
            return sv._tk_dialog(mode, title, filetypes), None
        except HTTPException as e:
            return None, e


def decoded_script(argv):
    """The PowerShell script carried by -EncodedCommand."""
    return base64.b64decode(argv[argv.index("-EncodedCommand") + 1]).decode("utf-16-le")


# ── dialogs ──────────────────────────────────────────────────────────────────

def test_dialog_wsl_uses_powershell():
    print("\nfile dialog under WSL")
    fake = FakeRun(stdout="E:\\data\\a.raw\r\n")
    res, err = dialog_on_wsl(fake)
    argv = fake.calls[0][0] if fake.calls else []
    check("routes to powershell.exe, not tkinter", argv[:1] == ["powershell.exe"], argv[:1])
    check("-STA (WinForms needs a single-threaded apartment)", "-STA" in argv)
    check("-NoProfile -NonInteractive", "-NoProfile" in argv and "-NonInteractive" in argv)
    check("Windows path returned as POSIX via wslpath -u", res == "/mnt/e/data/a.raw", res)
    check("generous timeout (minutes)", fake.calls[0][1].get("timeout", 0) >= 300)
    script = decoded_script(argv)
    check("dialog shown with a TopMost owner form", "TopMost = $true" in script
          and "ShowDialog($owner)" in script)
    check("UTF-8 console output", "[Console]::OutputEncoding = [Text.Encoding]::UTF8" in script)
    check("stdout decoded as UTF-8", fake.calls[0][1].get("encoding") == "utf-8")


def test_dialog_modes_and_filter():
    print("\ndialog modes, filter, multi-select, cancel")
    check("filter from filetypes",
          sv._win_filter([("RAW files", "*.raw"), ("All files", "*.*")])
          == "RAW files (*.raw)|*.raw|All files (*.*)|*.*")
    check("several patterns become ';'-separated",
          sv._win_filter([("FASTA files", "*.fasta *.fa")])
          == "FASTA files (*.fasta;*.fa)|*.fasta;*.fa")
    check("no filetypes: all files", sv._win_filter(None) == "All files (*.*)|*.*")
    fake = FakeRun(stdout="E:\\r\\ä b.raw\r\nE:\\r\\2.raw\r\n")
    res, _ = dialog_on_wsl(fake, "files", "Select RAW files",
                           [("RAW files", "*.raw"), ("All files", "*.*")])
    check("multi-select: several POSIX paths, non-ASCII and space intact",
          res == ["/mnt/e/r/ä b.raw", "/mnt/e/r/2.raw"], res)
    script = decoded_script(fake.calls[0][0])
    b64 = lambda t: base64.b64encode(t.encode()).decode()
    check("filter and mode reach the script",
          b64("RAW files (*.raw)|*.raw|All files (*.*)|*.*") in script and b64("files") in script)
    fake = FakeRun(stdout="E:\\out\r\n")
    res, _ = dialog_on_wsl(fake, "dir", "Select folder")
    check("folder dialog: FolderBrowserDialog, POSIX result",
          res == "/mnt/e/out" and "FolderBrowserDialog" in decoded_script(fake.calls[0][0]), res)
    for mode, empty in (("file", ""), ("files", []), ("dir", "")):
        res, err = dialog_on_wsl(FakeRun(stdout=""), mode)
        check(f"cancel ({mode}): empty result, no error", res == empty and err is None, (res, err))


def test_dialog_errors():
    print("\ndialog failures are explicit, never silent")
    _, err = dialog_on_wsl(FakeRun(raise_exc=subprocess.TimeoutExpired("powershell.exe", 600)))
    check("timeout → HTTP error with a message, not an empty path",
          err is not None and err.status_code == 504 and "paste" in err.detail, err)
    _, err = dialog_on_wsl(FakeRun(returncode=1))
    check("non-zero exit → HTTP 500 naming the failure",
          err is not None and err.status_code == 500 and "exit 1" in err.detail, err)
    _, err = dialog_on_wsl(FakeRun(), powershell=None)
    check("WSL without interop → 501, type or paste",
          err is not None and err.status_code == 501 and "paste" in err.detail, err)
    env = {k: v for k, v in os.environ.items() if k not in ("DISPLAY", "WAYLAND_DISPLAY")}
    with patched(sv, _is_wsl=lambda: False), patched(sv.sys, platform="linux"), \
         patched(os, environ=env):
        try:
            sv._tk_dialog("file", "x"); err = None
        except HTTPException as e:
            err = e
    check("Linux without a display → 501", err is not None and err.status_code == 501, err)


def test_dialog_title_injection():
    print("\nhostile dialog title cannot inject")
    hostile = "x'; Remove-Item C:\\ -Recurse; $(calc) \" `whoami` ;"
    fake = FakeRun(stdout="")
    dialog_on_wsl(fake, title=hostile)
    argv = fake.calls[0][0]
    check("title text appears in no argument", all(hostile not in a and "Remove-Item" not in a
                                                   for a in argv))
    script = decoded_script(argv)
    check("title absent from the script text itself", "Remove-Item" not in script)
    lit = re.search(r"\$title\s*=\s*Dec '([A-Za-z0-9+/=]*)'", script)
    check("title travels as a base64 literal and decodes back unchanged",
          lit and base64.b64decode(lit.group(1)).decode() == hostile)


def test_dialog_non_wsl_tkinter():
    print("\nnon-WSL keeps tkinter")
    fd = types.SimpleNamespace(askopenfilename=lambda **k: "/home/u/a.raw",
                               askopenfilenames=lambda **k: ("/a", "/b"),
                               askdirectory=lambda **k: "/d")

    class Tk:
        def withdraw(self): pass
        def attributes(self, *a): pass
        def destroy(self): pass
    tk = types.ModuleType("tkinter"); tk.Tk = Tk; tk.filedialog = fd
    ran = FakeRun()
    with patched(sv, _is_wsl=lambda: False), patched(sv.sys, platform="linux"), \
         patched(sv.subprocess, run=ran), \
         patched(os, environ={**os.environ, "DISPLAY": ":0"}), \
         patched(sys, modules={**sys.modules, "tkinter": tk, "tkinter.filedialog": fd}):
        res = sv._tk_dialog("file", "x")
    check("tkinter result returned unchanged", res == "/home/u/a.raw", res)
    check("no PowerShell started", ran.calls == [])


# ── boundary to Windows executables ──────────────────────────────────────────

class FakePopen:
    calls = []

    def __init__(self, cmd, **kw):
        FakePopen.calls.append(cmd)
        self.stdout, self.returncode = iter(["ok\n"]), 0

    def wait(self):
        return 0


def run_cmd_as(cmd, wsl, diann="diann"):
    FakePopen.calls = []
    with tempfile.TemporaryDirectory() as td, \
         patched(rp.subprocess, Popen=FakePopen), \
         patched(wslpaths, is_wsl=lambda: wsl, wslpath=fake_wslpath), \
         contextlib.redirect_stdout(io.StringIO()):
        try:
            rp.run_cmd(list(cmd), os.path.join(td, "log"), "t", diann)
            return FakePopen.calls[-1] if FakePopen.calls else None, None
        except ValueError as e:
            return (FakePopen.calls[-1] if FakePopen.calls else None), e


def test_boundary_run_cmd():
    print("\npath boundary to diann.exe")
    cmd = ["diann", "--f", "/mnt/e/x/a.raw", "--lib", "/mnt/e/x/lib.speclib",
           "--fixed-mod", "UniMod:4,57.021464,C", "--cut", "K*,R*,!*P",
           "--threads", "4", "--matrix-qvalue", "0.01"]
    got, _ = run_cmd_as(cmd, True, "C:\\DIA-NN\\1.8.1\\diann.exe")
    check("Windows-form exe executed as /mnt/c/...",
          got and got[0] == "/mnt/c/DIA-NN/1.8.1/diann.exe", got and got[0])
    check("/mnt/e/x/a.raw becomes E:\\x\\a.raw", got and got[2] == "E:\\x\\a.raw", got)
    check("non-path arguments untouched",
          got and got[5:] == ["--fixed-mod", "UniMod:4,57.021464,C", "--cut", "K*,R*,!*P",
                              "--threads", "4", "--matrix-qvalue", "0.01"], got)
    got, err = run_cmd_as(["diann", "--f", "/home/u/a.raw"], True, "/mnt/c/DIA-NN/diann.exe")
    check("Linux-filesystem path refused before DIA-NN runs, path named",
          got is None and err is not None and "/home/u/a.raw" in str(err), (got, err))
    for wsl, exe in ((False, "C:\\DIA-NN\\diann.exe"), (True, "/usr/local/bin/diann"),
                     (False, "/usr/local/bin/diann")):
        got, _ = run_cmd_as(cmd, wsl, exe)
        check(f"unchanged byte for byte (wsl={wsl}, exe={exe})",
              got == [exe] + cmd[1:], got)


def test_boundary_plotter():
    print("\nplotter gets Windows paths too")
    ran = FakeRun()
    with tempfile.TemporaryDirectory() as td, \
         patched(rp, _find_plotter=lambda exe: "/mnt/c/DIA-NN/dia-nn-plotter.exe"), \
         patched(rp.subprocess, run=ran), \
         patched(wslpaths, is_wsl=lambda: True, wslpath=fake_wslpath), \
         contextlib.redirect_stdout(io.StringIO()):
        d = os.path.join(td, "p")
        os.makedirs(d)
        for f in ("report.stats.tsv", "report.tsv"):
            open(os.path.join(d, f), "w").write("x")
        # the report folder must be on a Windows drive to be convertible
        with patched(rp.os.path, join=lambda *a: "/mnt/e/res/" + a[-1]), \
             patched(rp.os.path, exists=lambda p: not p.endswith(".pdf")):
            rp.run_plotter("/mnt/e/res", "/mnt/c/DIA-NN/diann.exe")
    argv = ran.calls[0][0] if ran.calls else []
    check("plotter argv in Windows form",
          argv == ["/mnt/c/DIA-NN/dia-nn-plotter.exe", "E:\\res\\report.stats.tsv",
                   "E:\\res\\report.tsv", "E:\\res\\report.pdf"], argv)


# ── normalisation: CLI, library guard, server ────────────────────────────────

def test_cli_paths_and_library_guard():
    print("\nCLI: Windows-form paths normalised; #18 guard still refuses")
    with tempfile.TemporaryDirectory() as td:
        run1 = os.path.join(td, "run1")
        os.makedirs(os.path.join(run1, "pass2_T"))
        lib = os.path.join(run1, "pass2_T", "combined_T.tsv")
        open(lib, "w").write("x" * 200)
        # built before modification support: no modifications
        open(os.path.join(run1, "pipeline.settings"), "w").write('{"enzymes": ["T"]}')

        def to_posix(p):          # E:\ is mounted at td in this test
            return p if not wslpaths.is_windows_path(p) else \
                td + "/" + p[3:].replace("\\", "/")
        argv = ["--raw", "E:\\raw\\a.raw", "E:\\raw\\b.raw", "--out-dir", "E:\\out",
                "--fasta", "E:\\db.fasta", "--lib", "E:\\run1\\pass2_T\\combined_T.tsv",
                "--speclib", "T=E:\\s\\t.speclib", "--diann", "E:\\DIA-NN\\diann.exe"]
        args = rp.normalise_cli_paths(rp.build_arg_parser().parse_args(argv), to_posix)
        check("--raw / --fasta / --out-dir / --diann made POSIX",
              args.raw == [[td + "/raw/a.raw", td + "/raw/b.raw"]]
              and args.fasta == td + "/db.fasta" and args.out_dir == td + "/out"
              and args.diann == td + "/DIA-NN/diann.exe", vars(args))
        check("--speclib CODE=PATH keeps its code", args.speclib == [f"T={td}/s/t.speclib"],
              args.speclib)
        expected = rp.mod_flags(rp.make_mod_config(["carbamidomethyl"]), include_extra=False)
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                rp.check_lib_mods(args.lib, expected, "Provided combined library (--lib)")
            refused = False
        except rp.ModConfigError:
            refused = True
        check("normalised Windows-form --lib: mismatched library refused", refused)
        buf = io.StringIO()
        try:
            with contextlib.redirect_stdout(buf):
                rp.check_lib_mods("E:\\run1\\pass2_T\\combined_T.tsv", expected, "x")
            raw_refused = False
        except rp.ModConfigError:
            raw_refused = True
        check("(why: un-normalised, the guard could only warn)",
              not raw_refused and "WARNING" in buf.getvalue())


def test_server_normalisation():
    print("\nserver: pasted Windows paths read the POSIX file")
    with tempfile.TemporaryDirectory() as td:
        os.makedirs(os.path.join(td, "res"))
        open(os.path.join(td, "res", "pipeline.settings"), "w").write('{"enzymes": ["T"]}')
        open(os.path.join(td, "res", "pipeline.log"), "w").write("x")

        def wp(flag, p):
            assert flag == "-u"
            return td + "/" + p[3:].replace("\\", "/")
        with patched(sv, _is_wsl=lambda: True), patched(wslpaths, wslpath=wp):
            try:
                got = sv.load_settings("E:\\res\\pipeline.settings")
            except HTTPException as e:
                got = e
            check("load-settings with E:\\...\\pipeline.settings reads the file",
                  got == {"enzymes": ["T"]}, got)
            chk = sv.check_outdir("E:\\res")
            check("check-outdir sees the existing run",
                  chk["hasExisting"] and chk["settingsPath"].startswith(td), chk)
            sv.asyncio.run(sv.save_settings_endpoint(
                {"path": "E:\\new\\pipeline.settings", "data": {"a": 1}}))
            check("save-settings writes to the POSIX location",
                  json.load(open(os.path.join(td, "new", "pipeline.settings"))) == {"a": 1})
    with patched(sv, _is_wsl=lambda: True), patched(wslpaths, wslpath=fake_wslpath):
        cmd = sv.build_pipeline_cmd(
            {"diannExe": "C:\\DIA-NN\\1.8.1\\diann.exe"},
            {"rawFiles": ["E:\\r\\a.raw"], "fasta": "E:\\db.fasta", "outDir": "E:\\out",
             "speclibPaths": {"T": "E:\\s\\t.speclib"}, "combinedLib": "E:\\c.tsv"})
    want = ["--raw", "/mnt/e/r/a.raw", "--fasta", "/mnt/e/db.fasta",
            "--out-dir", "/mnt/e/out", "--diann", "/mnt/c/DIA-NN/1.8.1/diann.exe",
            "--speclib", "T=/mnt/e/s/t.speclib", "--lib", "/mnt/e/c.tsv"]
    check("job submission: every path POSIX",
          all(w in cmd for w in want) and not any("\\" in c for c in cmd[3:]), cmd)
    with patched(sv, _is_wsl=lambda: False):
        cmd = sv.build_pipeline_cmd({}, {"rawFiles": ["E:\\r\\a.raw"], "outDir": "E:\\out"})
    check("not WSL: paths passed as given", "E:\\r\\a.raw" in cmd and "E:\\out" in cmd)


# ── start.py ─────────────────────────────────────────────────────────────────

def test_start_npm():
    print("\nstart.py: npm invocation and the rollup binary")
    ran = FakeRun()
    with patched(start.subprocess, run=ran), patched(start.sys, platform="linux"), \
         contextlib.redirect_stdout(io.StringIO()):
        start.npm("run", "build", cwd="/ui")
    cmd, kw = ran.calls[0]
    check("POSIX: a list, without shell=True", cmd == ["npm", "run", "build"]
          and not kw.get("shell"), (cmd, kw))
    ran = FakeRun()
    with patched(start.subprocess, run=ran), patched(start.sys, platform="win32"), \
         contextlib.redirect_stdout(io.StringIO()):
        start.npm("run", "build", cwd="C:\\ui")
    cmd, kw = ran.calls[0]
    check("Windows: one command string through the shell (npm.cmd)",
          cmd == "npm run build" and kw.get("shell") is True, (cmd, kw))
    with tempfile.TemporaryDirectory() as ui, patched(start.sys, platform="linux"):
        r = os.path.join(ui, "node_modules", "@rollup")
        os.makedirs(os.path.join(r, "rollup-win32-x64-msvc"))
        os.makedirs(os.path.join(r, "pluginutils"))
        check("Windows-installed modules: Linux rollup binary reported missing",
              start.rollup_platform_missing(ui))
        os.makedirs(os.path.join(r, "rollup-linux-x64-gnu"))
        check("Linux binary present: nothing missing", not start.rollup_platform_missing(ui))

def test_to_windows_idempotent():
    print("\nwslpaths.to_windows: a Windows-form path passes through")
    def runner(flag, path):
        raise AssertionError(f"wslpath {flag} called on {path}")
    for w in (r"E:\HighContentDIA\x.raw", "C:/DIA-NN/1.8.1/diann.exe"):
        check(f"{w!r} unchanged, wslpath not called", wslpaths.to_windows(w, runner) == w)
    check("POSIX path still converted",
          wslpaths.to_windows("/mnt/e/x.raw", lambda f, p: "E:\\x.raw") == "E:\\x.raw")


def test_start_ensure_ui():
    print("\nstart.ensure_ui: Windows-installed modules, end to end")
    with tempfile.TemporaryDirectory() as tmp, patched(start.sys, platform="linux"), \
         patched(start.wslpaths, is_wsl=lambda: True):
        ui = start.Path(tmp)
        os.makedirs(ui / "node_modules" / "@rollup" / "rollup-win32-x64-msvc")
        out = io.StringIO()
        code = None
        with contextlib.redirect_stdout(out):
            try:
                start.ensure_ui(force_build=True, ui=ui)
            except SystemExit as e:
                code = e.code
        check("not on a Windows drive: exits 1", code == 1, code)
        check("... and prints both options", "cmd.exe /c" in out.getvalue()
              and "npm install" in out.getvalue())
        ran = []
        with patched(start, run=lambda cmd, **kw: ran.append((cmd, kw))), \
             patched(start.wslpaths, is_windows_drive_path=lambda p: True), \
             contextlib.redirect_stdout(io.StringIO()):
            start.ensure_ui(force_build=True, ui=ui)
        check("on a Windows drive: builds with Windows npm in the given ui",
              ran == [(["cmd.exe", "/c", "npm run build"], {"cwd": str(ui)})], ran)


def main():
    for t in (test_dialog_wsl_uses_powershell, test_dialog_modes_and_filter,
              test_dialog_errors, test_dialog_title_injection, test_dialog_non_wsl_tkinter,
              test_boundary_run_cmd, test_boundary_plotter,
              test_cli_paths_and_library_guard, test_server_normalisation, test_start_npm,
              test_to_windows_idempotent, test_start_ensure_ui):
        t()
    print(f"\n{'all checks passed' if not FAIL else f'{len(FAIL)} FAILED: {FAIL}'}")
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
