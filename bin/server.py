#!/usr/bin/env python3
"""
mE-DIA Pipeline UI Server
==========================
FastAPI backend for the React frontend.

  python bin/server.py [--port 8765] [--no-browser]

Endpoints:
  GET  /api/version         build provenance (branch, commit, dirty, built when)
  GET  /api/config          load saved config
  POST /api/config          save config
  POST /api/runs            start pipeline run → {runId}
  GET  /api/runs/{id}       get run state
  DELETE /api/runs/{id}     cancel run
  GET  /api/browse/file     open native file dialog → {path}
  GET  /api/browse/files    open native multi-file dialog → {paths}
  GET  /api/browse/dir      open native folder dialog → {path}
  WS   /ws/{id}             stream pipeline output
  GET  /                    serve React app
"""

import asyncio
import json
import ntpath
import os
import re
import shutil
import signal
import subprocess
import sys
import uuid
import webbrowser
import importlib.util
from pathlib import Path
from typing import Dict, Any, Optional

from fastapi import FastAPI, WebSocket, WebSocketDisconnect, HTTPException
from fastapi.responses import FileResponse, JSONResponse
from fastapi.staticfiles import StaticFiles
from fastapi.middleware.cors import CORSMiddleware
import uvicorn

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wslpaths as _wslpaths  # noqa: E402  (bin/: WSL detection + path boundary)

# ─── Frozen self-exec ─────────────────────────────────────────────────────────
# In a self-contained build the .exe doubles as its own Python interpreter:
#
#   mEDIA.exe --exec-script <path/to/run_pipeline.py> [pipeline args...]
#
# runs that script as __main__ inside the frozen environment, so pandas and
# matplotlib resolve from the bundle and no system Python is needed. This must
# happen before the server starts.
if (getattr(sys, "frozen", False)
        and len(sys.argv) >= 3 and sys.argv[1] == "--exec-script"):
    import runpy
    _script  = sys.argv[2]
    sys.argv = [_script] + sys.argv[3:]
    runpy.run_path(_script, run_name="__main__")
    sys.exit(0)

# ─── Paths ────────────────────────────────────────────────────────────────────
# When bundled with PyInstaller (sys.frozen=True):
#   - Read-only data (UI dist) lives in sys._MEIPASS
#   - Writable files (config.json) live next to the .exe
#   - The pipeline sources stay on disk in bin/ next to the .exe, so they can be
#     patched without rebuilding:
#         mEDIA.exe               ← repo root
#         bin/run_pipeline.py
#         bin/proteoform_analysis.py
#         bin/viewer/viewer.html
#   - sys.executable is the .exe itself, not Python → find a Python for subprocess
#
# Python resolution order (for spawning run_pipeline.py):
#   1. <exe_dir>/python_portable/python.exe   ← shipped portable Python, if present
#   2. the .exe itself via --exec-script       ← self-contained build (pandas bundled)
#   3. system python3 / python on PATH         ← lean build, needs pandas installed
def _resolve_python(exe_dir: Path, portable_only: bool = False) -> Optional[str]:
    import shutil as _shutil
    for portable in (exe_dir / "python_portable" / "python.exe",
                     exe_dir / "bin" / "python_portable" / "python.exe",
                     exe_dir.parent / "python_portable" / "python.exe"):
        if portable.exists():
            return str(portable)
    if portable_only:
        return None
    return (_shutil.which("python3") or _shutil.which("python") or "python")


def _resolve_pipeline(exe_dir: Path, bundle: Path) -> Path:
    """Locate run_pipeline.py for the frozen .exe.

    On-disk sources next to the .exe win, so the pipeline can be patched without
    rebuilding. A self-contained build also carries a copy inside the bundle, so
    a bare mEDIA.exe with no bin/ folder still runs.
    """
    candidates = (exe_dir / "bin" / "run_pipeline.py",          # canonical repo layout
                  exe_dir.parent / "bin" / "run_pipeline.py",    # .exe inside dist/
                  exe_dir / "run_pipeline.py",                   # flat layout
                  bundle / "run_pipeline.py")                    # bundled fallback
    for cand in candidates:
        if cand.exists():
            return cand
    return candidates[0]   # canonical path; error surfaces at launch time


if getattr(sys, "frozen", False):
    _bundle     = Path(sys._MEIPASS)          # type: ignore[attr-defined]
    _exe_dir    = Path(sys.executable).parent
    BIN_DIR     = _bundle
    UI_DIST     = _bundle / "ui" / "dist"
    CONFIG_FILE = _exe_dir / "config.json"
    PIPELINE_PY = _resolve_pipeline(_exe_dir, _bundle)
    # A self-contained build ships pandas inside the bundle: run the pipeline
    # through the .exe itself (--exec-script) instead of hunting for a system
    # Python. A shipped python_portable/ still wins if the user provides one.
    _portable   = _resolve_python(_exe_dir, portable_only=True)
    if _portable:
        _PYTHON_EXE, _EXEC_PREFIX = _portable, []
    elif importlib.util.find_spec("pandas") is not None:
        _PYTHON_EXE, _EXEC_PREFIX = sys.executable, ["--exec-script"]
    else:
        _PYTHON_EXE, _EXEC_PREFIX = _resolve_python(_exe_dir), []
else:
    BIN_DIR      = Path(__file__).parent
    UI_DIST      = BIN_DIR / "ui" / "dist"
    CONFIG_FILE  = BIN_DIR / "config.json"
    PIPELINE_PY  = BIN_DIR / "run_pipeline.py"
    # In dev: prefer a project-level python_portable if a developer has built one,
    # otherwise use the dev's own interpreter.
    _dev_portable = BIN_DIR.parent / "python_portable" / "python.exe"
    _PYTHON_EXE   = str(_dev_portable) if _dev_portable.exists() else sys.executable
    _EXEC_PREFIX  = []

DEFAULT_CONFIG = {
    "diannExe": "diann",
    "threads":  14,
}

# ─── App setup ────────────────────────────────────────────────────────────────

app = FastAPI(title="mE-DIA Pipeline UI")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)


# ─── Run state management ────────────────────────────────────────────────────

class RunState:
    def __init__(self, run_id: str, cmd: list):
        self.run_id       = run_id
        self.cmd          = cmd
        self.status       = "running"
        self.output_lines = []
        self.pass_statuses: Dict[str, Any] = {}
        self.metrics      = []
        self.chart_data   = []
        self.clients      = set()
        self.process      = None
        self.task         = None
        self.cancel_requested = False

    def to_dict(self):
        return {
            "runId":        self.run_id,
            "status":       self.status,
            "passStatuses": self.pass_statuses,
            "metrics":      self.metrics,
            "chartData":    self.chart_data,
            "outputLines":  self.output_lines[-500:],  # last 500 lines via REST
        }


active_runs: Dict[str, RunState] = {}


def _parse_structured_line(run: RunState, line: str):
    """Parse ##STATUS##, ##METRIC##, and ##CHARTDATA## lines emitted by run_pipeline.py."""
    for tag in ("STATUS", "METRIC", "CHARTDATA"):
        prefix = f"##{tag}## "
        if line.startswith(prefix):
            try:
                data = json.loads(line[len(prefix):])
                if tag == "STATUS":
                    pass_name = data.get("pass", "pipeline")
                    run.pass_statuses[pass_name] = data
                    if data.get("pass") == "pipeline":
                        run.status = data.get("status", "running")
                elif tag == "METRIC":
                    run.metrics.append(data)
                elif tag == "CHARTDATA":
                    run.chart_data.append(data)
            except json.JSONDecodeError:
                pass
            break


async def _broadcast(run: RunState, message: str):
    dead = set()
    for ws in run.clients:
        try:
            await ws.send_text(message)
        except Exception:
            dead.add(ws)
    run.clients -= dead


async def _run_pipeline(run: RunState):
    """Execute the pipeline subprocess and stream output."""
    env = os.environ.copy()
    env["PYTHONUNBUFFERED"] = "1"
    try:
        proc = await asyncio.create_subprocess_exec(
            *run.cmd,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.STDOUT,
            env=env,
            # POSIX: its own process group, so a cancel reaches DIA-NN too
            **({} if sys.platform == "win32" else {"start_new_session": True}),
        )
        run.process = proc
        async for raw_line in proc.stdout:
            line = raw_line.decode("utf-8", errors="replace").rstrip()
            run.output_lines.append(line)
            _parse_structured_line(run, line)
            await _broadcast(run, json.dumps({"type": "line", "data": line}))
        await proc.wait()
        run.status = ("cancelled" if run.cancel_requested
                      else "done" if proc.returncode == 0 else "failed")
    except Exception as e:
        run.status = "cancelled" if run.cancel_requested else "failed"
        err = f"[server error] {e}"
        run.output_lines.append(err)
        await _broadcast(run, json.dumps({"type": "line", "data": err}))
    finally:
        await _broadcast(run, json.dumps({
            "type": "status",
            "data": {"status": run.status},
        }))


# ─── Config endpoints ─────────────────────────────────────────────────────────

@app.get("/api/version")
def get_version():
    """Provenance of this build: branch, commit, dirty flag, build time."""
    try:
        import build_info
        info = build_info.get()
        return {**info, "line": build_info.one_line(info)}
    except Exception as e:
        return {"commit": "unknown", "branch": "unknown", "line": f"unknown ({e})"}


@app.get("/api/config")
def get_config():
    cfg = dict(DEFAULT_CONFIG)
    if CONFIG_FILE.exists():
        try:
            cfg.update(json.loads(CONFIG_FILE.read_text()))
        except Exception:
            pass
    return cfg


@app.post("/api/config")
async def save_config(body: dict):
    cfg = dict(DEFAULT_CONFIG)
    cfg.update(body)
    CONFIG_FILE.write_text(json.dumps(cfg, indent=2))
    return {"ok": True}


# ─── Run endpoints ────────────────────────────────────────────────────────────

def build_pipeline_cmd(cfg: dict, job: dict) -> list:
    """run_pipeline.py command line for a GUI job. Raises ValueError on bad input.

    Modifications: when the job carries a ``modifications`` list it is passed
    exactly (``--mods none`` plus one flag per entry), so what the GUI shows is
    what runs. A job without the key (an older GUI or settings file) gets the
    pipeline default, which is recorded in pipeline.settings either way.
    """
    # Under WSL, Windows-form paths typed or pasted into the GUI become POSIX
    # here, so the pipeline sees one convention (to_posix raises ValueError).
    wsl = _is_wsl()
    pp  = (lambda p: _wslpaths.to_posix(p)) if wsl else (lambda p: p)

    raw_files = [pp(r) for r in job.get("rawFiles", [])]
    if not raw_files:
        raise ValueError("No raw files specified")

    fasta         = pp(job.get("fasta", "").strip())
    out_dir       = pp(job.get("outDir", "").strip())
    enzymes       = "".join(job.get("enzymes", ["T", "C", "G", "E"]))
    diann_exe     = pp(cfg.get("diannExe") or os.environ.get("MEDIA_DIANN_EXE")
                       or "diann")
    threads       = int(cfg.get("threads", 14))
    speclib_paths: Dict[str, Any] = job.get("speclibPaths", {})

    if not out_dir:
        raise ValueError("No output directory specified")

    # NOTE: --raw uses nargs="+", so all files go after a single --raw flag
    cmd = [_PYTHON_EXE, *_EXEC_PREFIX, str(PIPELINE_PY)]
    cmd += ["--raw"] + raw_files
    if fasta:
        cmd += ["--fasta", fasta]
    cmd += ["--out-dir", out_dir]
    cmd += ["--enzymes", enzymes]
    cmd += ["--threads", str(threads)]
    cmd += ["--diann", diann_exe]
    # Per-enzyme speclib paths: --speclib CODE=PATH
    for code, path in speclib_paths.items():
        if path and str(path).strip():
            cmd += ["--speclib", f"{code}={pp(path.strip())}"]
    if job.get("semiExt"):
        cmd.append("--semi-ext")
    if job.get("combinedLib"):
        cmd += ["--lib", pp(job["combinedLib"])]
    if job.get("minPepLen"):
        cmd += ["--min-pep-len", str(job["minPepLen"])]
    if job.get("maxPepLen"):
        cmd += ["--max-pep-len", str(job["maxPepLen"])]
    if job.get("missedCleavages") is not None:
        cmd += ["--missed-cleavages", str(job["missedCleavages"])]
    # Per-enzyme overrides of the missed cleavages: --enzyme-mc CODE=N. Only for
    # enzymes of this job; an enzyme without one keeps its own default (EX: 3)
    # or the global value above.
    selected = set(job.get("enzymes", ["T", "C", "G", "E"]))
    for code, n in (job.get("enzymeMissedCleavages") or {}).items():
        if code not in selected or n is None or str(n).strip() == "":
            continue
        try:
            n = int(n)
        except (TypeError, ValueError):
            raise ValueError(f"missed cleavages for {code} must be a whole number: {n!r}")
        if n < 0:
            raise ValueError(f"missed cleavages for {code} must be 0 or more: {n}")
        cmd += ["--enzyme-mc", f"{code}={n}"]

    mods = job.get("modifications")
    if mods is not None:
        cmd += ["--mods", "none"]
        for m in mods:
            kind = m.get("type")
            if kind not in ("fixed", "variable"):
                raise ValueError(f"modification type must be fixed or variable: {m}")
            cmd += ["--fixed-mod" if kind == "fixed" else "--var-mod",
                    f"{m.get('name', '')},{m.get('mass', '')},{m.get('sites', '')}"]
            # PTM-of-interest scoring and cleavage blocking travel with the modification,
            # so what the GUI shows is still what runs
            if m.get("monitor"):
                cmd += ["--monitor-mod", str(m.get("name", ""))]
            if m.get("noCutAfter"):
                cmd += ["--no-cut-after-mod", str(m.get("name", ""))]
        cmd += ["--var-mods", str(int(job.get("maxVarMods", 1)))]
    if job.get("metExcision"):
        cmd.append("--met-excision")
    for a in job.get("extraDiannArgs") or []:
        if str(a).strip():
            cmd.append(f"--diann-arg={str(a).strip()}")
    return cmd


@app.post("/api/runs")
async def start_run(body: dict):
    """
    body = {
      config: {diannExe, threads},
      job:    {rawFiles, fasta, outDir, enzymes, speclibPaths, semiExt, skipPass1,
               combinedLib, minPepLen, maxPepLen, missedCleavages,
               enzymeMissedCleavages: {CODE: N} (overrides only),
               modifications: [{name, mass, sites, type: fixed|variable}],
               maxVarMods, metExcision, extraDiannArgs: [str]}
    }
    """
    try:
        cmd = build_pipeline_cmd(body.get("config", {}), body.get("job", {}))
    except ValueError as e:
        raise HTTPException(400, str(e))

    if not PIPELINE_PY.exists():
        raise HTTPException(
            500,
            f"Pipeline script not found: {PIPELINE_PY}. Keep mEDIA.exe in the "
            f"repository root, with the bin/ folder alongside it."
        )

    run_id = str(uuid.uuid4())[:8]
    run    = RunState(run_id, cmd)
    active_runs[run_id] = run
    run.task = asyncio.create_task(_run_pipeline(run))
    return {"runId": run_id}


@app.get("/api/runs/{run_id}")
def get_run(run_id: str):
    run = active_runs.get(run_id)
    if not run:
        raise HTTPException(404, "Run not found")
    return run.to_dict()


async def _cancel_running_passes(run: RunState):
    """Mark every pass still 'running' as cancelled, the way the pipeline itself
    reports a pass: a ##STATUS## line, so the REST state and the GUI's pass rows
    both update. A killed pipeline cannot send these itself."""
    for name, st in list(run.pass_statuses.items()):
        if isinstance(st, dict) and st.get("status") == "running":
            line = "##STATUS## " + json.dumps({"pass": name, "status": "cancelled"})
            run.output_lines.append(line)
            _parse_structured_line(run, line)
            await _broadcast(run, json.dumps({"type": "line", "data": line}))


@app.delete("/api/runs/{run_id}")
async def cancel_run(run_id: str):
    """Stop the run and everything it started; answers once they have exited.

    {"ok": true, "status": "cancelled"}, or "ok": false with "stillRunning"
    listing PIDs that did not exit within CANCEL_WAIT_S.
    """
    run = active_runs.get(run_id)
    if not run:
        raise HTTPException(404, "Run not found")
    if run.status in ("done", "failed", "cancelled"):
        return {"ok": True, "status": run.status, "stillRunning": []}
    run.cancel_requested = True
    remaining = await _cancel_tree(run)
    if run.task:
        try:
            await asyncio.wait_for(asyncio.shield(run.task), CANCEL_WAIT_S)
        except asyncio.TimeoutError:
            pass
    await _cancel_running_passes(run)
    run.status = "cancelled"
    if remaining:
        note = (f"[server] cancel: process(es) {', '.join(map(str, remaining))} did "
                f"not exit within {CANCEL_WAIT_S} s; stop them before reusing the "
                f"output folder.")
        run.output_lines.append(note)
        await _broadcast(run, json.dumps({"type": "line", "data": note}))
    await _broadcast(run, json.dumps({
        "type": "status", "data": {"status": "cancelled"},
    }))
    return {"ok": not remaining, "status": "cancelled", "stillRunning": remaining}


@app.post("/api/shutdown")
async def shutdown_server():
    """Terminate any in-flight pipeline subprocesses, then exit the server.

    The response returns *before* the process exits so the UI sees a clean
    success. uvicorn doesn't expose a portable graceful-stop hook from a
    request handler, so we schedule ``os._exit(0)`` on the loop. This is
    deliberately hard — it bypasses uvicorn's atexit cleanups, which is what
    we want when the user has clicked 'Shutdown' and expects the port freed
    immediately.
    """
    terminated = 0
    for run in active_runs.values():
        try:
            if run.process and run.process.returncode is None:
                run.cancel_requested = True
                await _cancel_tree(run)
                terminated += 1
        except Exception:
            pass
    loop = asyncio.get_event_loop()
    loop.call_later(0.4, lambda: os._exit(0))
    return {"ok": True, "terminatedRuns": terminated}


# ─── WebSocket ────────────────────────────────────────────────────────────────

@app.websocket("/ws/{run_id}")
async def ws_endpoint(websocket: WebSocket, run_id: str):
    await websocket.accept()
    run = active_runs.get(run_id)
    if not run:
        await websocket.send_text(json.dumps({
            "type": "error", "data": "Run not found",
        }))
        await websocket.close()
        return

    # Replay buffered output
    for line in run.output_lines:
        await websocket.send_text(json.dumps({"type": "line", "data": line}))

    # If already finished, send final status and close
    if run.status in ("done", "failed", "cancelled"):
        await websocket.send_text(json.dumps({
            "type": "status", "data": {"status": run.status},
        }))
        await websocket.close()
        return

    run.clients.add(websocket)
    try:
        while True:
            # Keep alive — client sends pings; we just wait
            try:
                await asyncio.wait_for(websocket.receive_text(), timeout=30)
            except asyncio.TimeoutError:
                pass
            if run.status in ("done", "failed", "cancelled"):
                break
    except WebSocketDisconnect:
        pass
    finally:
        run.clients.discard(websocket)


# ─── Native file/folder dialogs ───────────────────────────────────────────────
# One choke point for all five browse endpoints. Native Windows and Linux with a
# display use tkinter, as before. Under WSL a Tk dialog would be a Linux window
# (WSLg) that blocks the request unseen behind the Windows browser and returns
# paths diann.exe cannot read, so the dialog is opened on the Windows host via
# PowerShell instead, and its Windows paths are converted to POSIX.

WIN_DIALOG_TIMEOUT_S = 600     # a person is choosing a file: minutes, not seconds

# Title, filter and mode arrive as base64 literals (substituted below) and are
# decoded inside PowerShell, so text from the browser is never spliced into the
# command: base64 cannot contain a quote, ';' or '$('.
_PS_DIALOG = r"""
$ErrorActionPreference = 'Stop'
[Console]::OutputEncoding = [Text.Encoding]::UTF8
function Dec($b) { [Text.Encoding]::UTF8.GetString([Convert]::FromBase64String($b)) }
$mode   = Dec '__MODE__'
$title  = Dec '__TITLE__'
$filter = Dec '__FILTER__'
Add-Type -AssemblyName System.Windows.Forms
Add-Type -AssemblyName System.Drawing
# An invisible, top-most owner keeps the dialog in front of the browser.
$owner = New-Object System.Windows.Forms.Form
$owner.TopMost = $true
$owner.ShowInTaskbar = $false
$owner.Opacity = 0
$owner.StartPosition = 'CenterScreen'
$owner.Size = New-Object System.Drawing.Size(1, 1)
$owner.Show()
$owner.Activate()
if ($mode -eq 'dir') {
    $d = New-Object System.Windows.Forms.FolderBrowserDialog
    $d.Description = $title
    $d.ShowNewFolderButton = $true
    if ($d.ShowDialog($owner) -eq [System.Windows.Forms.DialogResult]::OK) {
        [Console]::Out.WriteLine($d.SelectedPath)
    }
} else {
    $d = New-Object System.Windows.Forms.OpenFileDialog
    $d.Title = $title
    $d.Filter = $filter
    $d.Multiselect = ($mode -eq 'files')
    if ($d.ShowDialog($owner) -eq [System.Windows.Forms.DialogResult]::OK) {
        foreach ($f in $d.FileNames) { [Console]::Out.WriteLine($f) }
    }
}
$owner.Close()
"""


def _b64(text: str) -> str:
    import base64
    return base64.b64encode(text.encode("utf-8")).decode("ascii")


def _win_filter(filetypes) -> str:
    """tkinter filetypes -> WinForms Filter, e.g. 'RAW files (*.raw)|*.raw|All files (*.*)|*.*'."""
    parts = []
    for label, patterns in (filetypes or [("All files", "*.*")]):
        pats = ";".join(str(patterns).split())
        parts.append(f"{label} ({pats})|{pats}")
    return "|".join(parts)


def _powershell_dialog(mode: str, title: str, filetypes=None):
    """Open the dialog on the Windows host; returns POSIX path(s), '' / [] on cancel."""
    import shutil
    if not shutil.which("powershell.exe"):
        raise HTTPException(501, "No file dialog available: WSL cannot reach "
                                 "powershell.exe (Windows interop is off). Type or "
                                 "paste the path instead.")
    script = (_PS_DIALOG.replace("__MODE__", _b64(mode))
                        .replace("__TITLE__", _b64(title))
                        .replace("__FILTER__", _b64(_win_filter(filetypes))))
    encoded = _b64_utf16(script)
    try:
        done = subprocess.run(
            ["powershell.exe", "-NoProfile", "-NonInteractive", "-STA",
             "-EncodedCommand", encoded],
            capture_output=True, text=True, encoding="utf-8", errors="replace",
            timeout=WIN_DIALOG_TIMEOUT_S,
            # cmd/PowerShell cannot start in a \\wsl$ directory; a drive path avoids
            # the warning and the fallback to C:\Windows.
            cwd="/mnt/c" if os.path.isdir("/mnt/c") else None)
    except subprocess.TimeoutExpired:
        raise HTTPException(504, f"The file dialog was not answered within "
                                 f"{WIN_DIALOG_TIMEOUT_S // 60} min. Try again, or "
                                 f"type or paste the path.")
    if done.returncode != 0:
        msg = (done.stderr or done.stdout or "").strip()[-400:]
        raise HTTPException(500, f"The Windows file dialog failed (exit "
                                 f"{done.returncode}): {msg}")
    picked = [ln.strip() for ln in done.stdout.splitlines() if ln.strip()]
    try:
        posix = [_wslpaths.to_posix(p) for p in picked]
    except ValueError as e:
        raise HTTPException(500, f"Could not convert the chosen path for WSL: {e}")
    if mode == "files":
        return posix
    return posix[0] if posix else ""


def _b64_utf16(script: str) -> str:
    """PowerShell -EncodedCommand wants base64 of UTF-16LE."""
    import base64
    return base64.b64encode(script.encode("utf-16-le")).decode("ascii")


# ─── Cancelling a run: the whole process tree ────────────────────────────────
#
# The pipeline starts DIA-NN as its child, so stopping only the pipeline leaves
# DIA-NN searching on (#31). On Windows, TerminateProcess does not reach
# children. Under WSL, diann.exe is a Windows process behind an interop relay,
# and a signal to the Linux process group does not stop it either. So:
#   1. stop the pipeline, so it starts nothing new: taskkill /T /F on Windows
#      (which takes DIA-NN with it), a process-group kill on POSIX;
#   2. on Windows and under WSL, find any DIA-NN process still working in this
#      run's output folder and kill its tree with taskkill;
#   3. wait until both have exited before answering.

CANCEL_WAIT_S = 30
WINDOWS_TOOLS = ("diann.exe", "dia-nn-plotter.exe")

# Prints the PIDs of processes named in __NAMES__ whose command line mentions
# the folder __DIR__ (lower case, backslashes, trailing backslash). Both are
# base64 literals, like the dialog's.
_PS_RUN_PIDS = r"""
$dir   = [Text.Encoding]::UTF8.GetString([Convert]::FromBase64String('__DIR__'))
$names = [Text.Encoding]::UTF8.GetString([Convert]::FromBase64String('__NAMES__')).Split('|')
Get-CimInstance Win32_Process | Where-Object {
    $_.ProcessId -ne $PID -and $_.CommandLine -and ($names -contains $_.Name.ToLower()) -and
    $_.CommandLine.Replace('/', '\').ToLower().Contains($dir)
} | ForEach-Object { $_.ProcessId }
"""


def _arg_after(cmd, flag) -> Optional[str]:
    cmd = [str(c) for c in cmd]
    return cmd[cmd.index(flag) + 1] if flag in cmd[:-1] else None


def _windows_run_dir(cmd) -> Optional[str]:
    """The run's output folder as DIA-NN sees it: lower case, Windows form, ending in \\."""
    d = _arg_after(cmd, "--out-dir")
    if not d:
        return None
    if _is_wsl():
        if not (_wslpaths.is_windows_drive_path(d) or _wslpaths.is_windows_path(d)):
            return None              # a Linux folder: no Windows program writes there
        try:
            d = _wslpaths.to_windows(d)
        except Exception:
            return None
    elif sys.platform == "win32":
        d = ntpath.abspath(d)          # also when run from a test on POSIX
    else:
        return None
    return d.replace("/", "\\").rstrip("\\").lower() + "\\"


def _tool_names(cmd) -> list:
    names = set(WINDOWS_TOOLS)
    exe = _arg_after(cmd, "--diann")
    if exe:
        names.add(re.split(r"[\\/]", exe)[-1].lower())
    return sorted(names)


def _win_cwd():
    # cmd/PowerShell cannot start in a \\wsl$ directory
    return "/mnt/c" if _is_wsl() and os.path.isdir("/mnt/c") else None


def _windows_run_pids(cmd, runner=None) -> list:
    """PIDs of DIA-NN processes working in this run's output folder (Windows/WSL)."""
    d = _windows_run_dir(cmd)
    if not d or not shutil.which("powershell.exe"):
        return []
    script = (_PS_RUN_PIDS.replace("__DIR__", _b64(d))
                          .replace("__NAMES__", _b64("|".join(_tool_names(cmd)))))
    try:
        done = (runner or subprocess.run)(["powershell.exe", "-NoProfile", "-NonInteractive",
                       "-EncodedCommand", _b64_utf16(script)],
                      capture_output=True, text=True, encoding="utf-8",
                      errors="replace", timeout=60, cwd=_win_cwd())
    except Exception:
        return []
    return [int(t) for t in (done.stdout or "").split() if t.isdigit()]


def _taskkill(pid: int, runner=None):
    """Kill a Windows process and its children (from Windows or WSL)."""
    exe = "taskkill" if sys.platform == "win32" else "taskkill.exe"
    try:
        (runner or subprocess.run)([exe, "/T", "/F", "/PID", str(pid)],
                                   capture_output=True, text=True, timeout=60,
                                   cwd=_win_cwd())
    except Exception:
        pass


def _stop_pipeline(proc, runner=None):
    """Stop the pipeline process and, as far as the OS allows, its children."""
    if proc is None or proc.returncode is not None:
        return
    if sys.platform == "win32":
        _taskkill(proc.pid, runner)
        return
    try:
        os.killpg(proc.pid, signal.SIGKILL)      # start_new_session: pgid == pid
    except (ProcessLookupError, PermissionError):
        pass


async def _cancel_tree(run: "RunState") -> list:
    """Stop the run's pipeline and every DIA-NN process working for it, and wait
    for them to exit. Returns the PIDs still alive after CANCEL_WAIT_S ([] = all gone)."""
    loop = asyncio.get_running_loop()
    deadline = loop.time() + CANCEL_WAIT_S
    proc = run.process
    await asyncio.to_thread(_stop_pipeline, proc)
    if proc is not None:
        try:
            await asyncio.wait_for(proc.wait(), CANCEL_WAIT_S)
        except asyncio.TimeoutError:
            pass
    while True:
        pids = await asyncio.to_thread(_windows_run_pids, run.cmd)
        if not pids or loop.time() > deadline:
            break
        for pid in pids:
            await asyncio.to_thread(_taskkill, pid)
        await asyncio.sleep(1)
    if proc is not None and proc.returncode is None:
        pids = [proc.pid] + list(pids)
    return pids


def _tk_dialog(mode: str, title: str, filetypes=None):
    """Open a file/folder dialog. Raises HTTPException (501) if none can be shown."""
    if _is_wsl():
        return _powershell_dialog(mode, title, filetypes)
    if (sys.platform.startswith("linux")
            and not (os.environ.get("DISPLAY") or os.environ.get("WAYLAND_DISPLAY"))):
        raise HTTPException(501, "No file dialog available: there is no display. "
                                 "Type or paste the path instead.")
    try:
        import tkinter as tk
        from tkinter import filedialog
        root = tk.Tk()
        root.withdraw()
        root.attributes("-topmost", True)
        if mode == "file":
            result = filedialog.askopenfilename(title=title, filetypes=filetypes or [])
        elif mode == "files":
            result = list(filedialog.askopenfilenames(title=title, filetypes=filetypes or []))
        else:  # dir
            result = filedialog.askdirectory(title=title)
        root.destroy()
        return result
    except Exception as e:
        raise HTTPException(501, f"No file dialog available ({e}). "
                                 f"Type or paste the path instead.")


def _posix_path(path: str) -> str:
    """Under WSL, a pasted Windows path (E:\\...) -> POSIX; otherwise unchanged."""
    if path and _is_wsl() and _wslpaths.is_windows_path(path):
        try:
            return _wslpaths.to_posix(path)
        except ValueError as e:
            raise HTTPException(400, str(e))
    return path


@app.get("/api/browse/file")
def browse_file(title: str = "Select file"):
    result = _tk_dialog("file", title)
    return {"path": result or ""}


@app.get("/api/browse/files")
def browse_files(title: str = "Select RAW files",
                 ext: str = ".raw"):
    ft = [(f"{ext.upper()} files", f"*{ext}"), ("All files", "*.*")]
    result = _tk_dialog("files", title, ft)
    return {"paths": result if isinstance(result, list) else []}


@app.get("/api/browse/dir")
def browse_dir(title: str = "Select folder"):
    result = _tk_dialog("dir", title)
    return {"path": result or ""}


@app.get("/api/check-outdir")
def check_outdir(path: str = ""):
    """Return whether pipeline.log already exists in the given output directory."""
    if not path:
        return {"hasExisting": False, "settingsPath": ""}
    d = Path(_posix_path(path))
    log_exists = (d / "pipeline.log").exists()
    settings_p = d / "pipeline.settings"
    return {
        "hasExisting":  log_exists,
        "settingsPath": str(settings_p) if settings_p.exists() else "",
    }


@app.get("/api/load-settings")
def load_settings(path: str = ""):
    """Read a pipeline.settings JSON file and return its contents."""
    if not path:
        raise HTTPException(400, "No path specified")
    p = Path(_posix_path(path))
    if not p.exists():
        raise HTTPException(404, "File not found")
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except Exception as e:
        raise HTTPException(400, f"Cannot parse settings: {e}")


@app.post("/api/save-settings")
async def save_settings_endpoint(body: dict):
    """Write a pipeline.settings JSON file to the requested path."""
    path = body.get("path", "").strip()
    data = body.get("data", {})
    if not path:
        raise HTTPException(400, "No path specified")
    p = Path(_posix_path(path))
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(data, indent=2), encoding="utf-8")
    return {"ok": True, "path": str(p)}


@app.get("/api/browse/fasta")
def browse_fasta():
    ft = [("FASTA files", "*.fasta *.fa"), ("All files", "*.*")]
    result = _tk_dialog("file", "Select FASTA database", ft)
    return {"path": result or ""}


@app.get("/api/browse/exe")
def browse_exe():
    ft = [("Executable", "*.exe diann*"), ("All files", "*.*")]
    result = _tk_dialog("file", "Select DIA-NN executable", ft)
    return {"path": result or ""}


# ─── Static React app ─────────────────────────────────────────────────────────

if UI_DIST.exists():
    app.mount("/assets", StaticFiles(directory=str(UI_DIST / "assets")), name="assets")

    @app.get("/{full_path:path}")
    def serve_spa(full_path: str = ""):
        index = UI_DIST / "index.html"
        return FileResponse(str(index))
else:
    @app.get("/")
    def no_ui():
        return JSONResponse({
            "message": "React UI not built. Run: cd bin/ui && npm install && npm run build",
        })


# ─── Opening the browser ─────────────────────────────────────────────────────

# One detector for the whole code base (bin/wslpaths.py); the alias keeps the
# name this module and its tests have always used.
_is_wsl = _wslpaths.is_wsl


def _wsl_open(url: str) -> bool:
    """Hand `url` to the Windows host's default browser, from inside WSL.

    A WSL distro normally has no browser of its own, so `webbrowser.open()`
    falls through to `xdg-open`, which walks its candidate list and gives up
    with "no method available for opening ...". The Windows side has to do it.
    Tried in order: wslview (wslu, honours the Windows default browser),
    PowerShell, cmd.exe, and finally explorer.exe -- which works but exits
    non-zero even then, so its status cannot be checked.
    """
    # cmd.exe cannot cd into a \\wsl$ UNC path and warns about it; starting it
    # from a Windows-visible directory keeps that off the console.
    cwd = "/mnt/c" if os.path.isdir("/mnt/c") else None
    attempts = [
        (["wslview", url], None),
        (["powershell.exe", "-NoProfile", "-NonInteractive",
          "-Command", f"Start-Process '{url}'"], None),
        (["cmd.exe", "/c", "start", "", url], cwd),
    ]
    for cmd, cmd_cwd in attempts:
        try:
            done = subprocess.run(cmd, cwd=cmd_cwd, timeout=20,
                                  stdout=subprocess.DEVNULL,
                                  stderr=subprocess.DEVNULL)
        except (OSError, subprocess.SubprocessError):
            continue
        if done.returncode == 0:
            return True
    try:
        subprocess.run(["explorer.exe", url], timeout=20,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return True          # exit status is meaningless here
    except (OSError, subprocess.SubprocessError):
        return False


def _platform() -> str:
    """`sys.platform`, behind a function so tests can fake it, as they fake _is_wsl.

    Without this, tests run under Windows python.exe could never reach the WSL
    or plain-Linux branch of open_browser(), and would instead open a real
    browser tab through os.startfile (#40).
    """
    return sys.platform


def _startfile(url: str) -> None:
    """`os.startfile`, behind a function so tests never launch a real browser."""
    os.startfile(url)  # type: ignore[attr-defined]  # Windows only


def open_browser(url: str) -> bool:
    """Open `url` in a browser. False if no launcher could be found."""
    if _platform() == "win32":
        try:
            _startfile(url)     # most reliable on Windows
            return True
        except Exception:
            pass
    elif _is_wsl():
        if _wsl_open(url):
            return True
        # fall through: an X/Wayland session with a real Linux browser is
        # possible under WSLg, and webbrowser would find it.
    try:
        return bool(webbrowser.open(url))
    except Exception:
        return False


# ─── Entry point ─────────────────────────────────────────────────────────────

def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--port",       type=int, default=8765)
    ap.add_argument("--no-browser", action="store_true")
    ap.add_argument("--version",    action="store_true",
                    help="print build provenance (branch, commit, build time) and exit")
    args = ap.parse_args()

    if args.version:
        try:
            import build_info
            print(build_info.one_line())
            print(build_info.as_json())
        except Exception as e:
            print(f"build provenance unavailable: {e}")
        return

    if not args.no_browser:
        import threading
        url = f"http://localhost:{args.port}"
        def _open():
            import time
            time.sleep(2.0)   # wait for uvicorn to be ready
            if not open_browser(url):
                print(f"\n  No browser could be opened automatically."
                      f"\n  Open {url} yourself"
                      + ("  (WSL: install wslu for `wslview`, "
                         "or just paste the URL into Windows)" if _is_wsl() else "")
                      + "\n")
        threading.Thread(target=_open, daemon=True).start()

    print(f"\n  mE-DIA Pipeline UI")
    print(f"  Open: http://localhost:{args.port}")
    if not UI_DIST.exists():
        print(f"\n  NOTE: React UI not built.")
        print(f"  Run: cd bin/ui && npm install && npm run build\n")
    uvicorn.run(app, host="0.0.0.0", port=args.port, log_level="warning")


if __name__ == "__main__":
    main()
