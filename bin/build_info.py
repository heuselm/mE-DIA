#!/usr/bin/env python3
"""
Where did this build come from?

A built `mEDIA.exe` is not tracked by git, so its file name is the only hint at
what is inside it — and a file name survives renaming, copying and downloading
from a Release. That is how a binary built from one branch ends up sitting next
to the sources of another. So provenance is recorded in three places:

1. inside the artifact: ``bin/build_exe.py`` writes ``bin/_build_info.py``
   (generated, gitignored) before PyInstaller runs, and it is bundled;
2. on the artifact: the file is named ``mEDIA-<branch>-<sha>[-dirty].exe``;
3. in every result folder: ``run_pipeline.py`` writes this dict into
   ``pipeline.settings`` and prints it at the top of its log, so a result can
   always be traced to the code that produced it.

Reading order:
- a frozen build reads the generated file -- there is no repository to ask;
- a source checkout reads live ``git``, and falls back to the generated file
  only if git cannot answer. The generated file describes the last exe *build*,
  which a checkout keeps (it is gitignored) and then moves away from, so
  preferring it there records the wrong commit (#27);
- otherwise ``unknown``.
A value read from the generated file is always labelled ``source: "stamped"``.
Nothing here raises.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
from typing import Dict, Optional

_HERE = os.path.dirname(os.path.abspath(__file__))
UNKNOWN = "unknown"


def _git(*args: str, cwd: Optional[str] = None) -> Optional[str]:
    try:
        out = subprocess.run(("git",) + args, cwd=cwd or _HERE, capture_output=True,
                             text=True, timeout=10)
    except (OSError, subprocess.SubprocessError):
        return None
    return out.stdout.strip() if out.returncode == 0 else None


def from_git(cwd: Optional[str] = None) -> Dict[str, object]:
    """Provenance of the working tree, or ``unknown`` outside a repository."""
    sha = _git("rev-parse", "--short=7", "HEAD", cwd=cwd)
    if not sha:
        return {"commit": UNKNOWN, "branch": UNKNOWN, "dirty": False,
                "describe": UNKNOWN, "source": "no-git"}
    branch = _git("rev-parse", "--abbrev-ref", "HEAD", cwd=cwd) or UNKNOWN
    status = _git("status", "--porcelain", cwd=cwd)
    return {"commit": sha,
            "branch": branch,
            "dirty": bool(status),
            "describe": _git("describe", "--tags", "--always", "--dirty", cwd=cwd) or sha,
            "source": "git"}


def _is_frozen() -> bool:
    return bool(getattr(sys, "frozen", False)) or hasattr(sys, "_MEIPASS")


def read_stamp() -> Optional[Dict[str, object]]:
    """The dict bin/build_exe.py wrote into bin/_build_info.py, or None."""
    try:
        import _build_info                      # written by build_exe.py
        return dict(_build_info.BUILD_INFO)
    except Exception:
        return None


def get(frozen: Optional[bool] = None, stamp=read_stamp, git=from_git) -> Dict[str, object]:
    """Provenance of the code that is running now.

    ``frozen``, ``stamp`` and ``git`` are injectable for tests; callers pass nothing.
    """
    frozen = _is_frozen() if frozen is None else frozen

    def stamped(s: Dict[str, object]) -> Dict[str, object]:
        # the stamp was built from from_git() and so says source "git"; say where
        # this value actually came from, or a stale stamp passes for live git
        s = dict(s)
        s["source"] = "stamped"
        s["frozen"] = frozen
        return s

    if frozen:
        s = stamp()
        if s:
            return stamped(s)
    info = git()
    if info.get("source") == "git":
        info["frozen"] = frozen
        return info
    s = stamp()
    if s:
        return stamped(s)
    info["frozen"] = frozen
    return info


def artifact_name(info: Dict[str, object], base: str = "mEDIA",
                  ext: str = ".exe") -> str:
    """
    File name that cannot lie: ``mEDIA-<branch>-<sha>[-dirty].exe``.

    The branch is sanitised to ``[A-Za-z0-9._-]`` so a name like
    ``feature/ms1-rt-constrained-semispecific`` cannot create a directory.
    """
    branch = str(info.get("branch") or UNKNOWN)
    branch = "".join(c if (c.isalnum() or c in "._-") else "-" for c in branch).strip("-")
    branch = branch or UNKNOWN
    parts = [base, branch, str(info.get("commit") or UNKNOWN)]
    if info.get("dirty"):
        parts.append("dirty")
    return "-".join(parts) + ext


def one_line(info: Optional[Dict[str, object]] = None) -> str:
    """Human-readable single line, for logs, ``--version`` and the GUI."""
    i = info if info is not None else get()
    bits = [f"{i.get('branch', UNKNOWN)}@{i.get('commit', UNKNOWN)}"]
    if i.get("dirty"):
        bits.append("DIRTY (uncommitted changes at build time)" if i.get("source") == "stamped"
                    else "DIRTY (uncommitted changes)")
    if i.get("built_utc"):
        bits.append(f"built {i['built_utc']}")
    if i.get("ui_built_utc"):
        bits.append(f"UI {i['ui_built_utc']}")
    return "mE-DIA " + ", ".join(bits)


# ── What ran it: interpreter, packages, DIA-NN (#38) ───────────────────────────
# Third-party packages the pipeline imports: run_pipeline.py (pandas) and
# proteoform_analysis.py (pandas, matplotlib); numpy carries pandas' arithmetic.
PIPELINE_PACKAGES = ("pandas", "numpy", "matplotlib")


def package_version(name: str) -> str:
    """Version of the package actually loaded in this process if it is, else the
    installed one; ``unknown`` if neither can be read (e.g. a frozen build
    without package metadata)."""
    try:
        mod = sys.modules.get(name)
        v = getattr(mod, "__version__", None) if mod is not None else None
        if v:
            return str(v)
        from importlib import metadata
        return metadata.version(name)
    except Exception:
        return UNKNOWN


def runtime(diann_exe: Optional[str] = None, packages=PIPELINE_PACKAGES,
            frozen: Optional[bool] = None) -> Dict[str, object]:
    """The interpreter, packages and DIA-NN a run uses. The DIA-NN banner is not
    known until DIA-NN has run once; it starts ``unknown`` and the pipeline fills
    it in from the first call's output. Nothing here raises."""
    def safe(fn):
        try:
            v = fn()
            return str(v) if v not in (None, "") else UNKNOWN
        except Exception:
            return UNKNOWN
    import platform
    return {
        "python":      safe(lambda: "%d.%d.%d" % tuple(sys.version_info[:3])),
        "python_impl": safe(platform.python_implementation),
        "platform":    safe(platform.platform),
        "executable":  safe(lambda: sys.executable),
        "frozen":      _is_frozen() if frozen is None else bool(frozen),
        "packages":    {n: package_version(n) for n in packages},
        "diann":       {"exe": str(diann_exe) if diann_exe else UNKNOWN, "banner": UNKNOWN},
    }


def parse_diann_banner(lines) -> str:
    """DIA-NN's own version line from the first lines it prints, e.g.
    'DIA-NN 1.8.1 (Data-Independent Acquisition by Neural Networks) Compiled on
    Apr 14 2022 15:31:19'; ``unknown`` if there is none."""
    try:
        text = [str(x).strip() for x in lines]
        for i, line in enumerate(text):
            if line.startswith("DIA-NN"):
                nxt = text[i + 1] if i + 1 < len(text) else ""
                return f"{line} {nxt}" if nxt.startswith("Compiled on") else line
    except Exception:
        pass
    return UNKNOWN


def runtime_one_line(rt: Optional[Dict[str, object]]) -> str:
    """'Python 3.14.0 (CPython, Windows-...) | pandas 2.2.3, numpy 2.1.0, ...'"""
    try:
        pk = ", ".join(f"{k} {v}" for k, v in (rt.get("packages") or {}).items())
        frozen = ", frozen" if rt.get("frozen") else ""
        return (f"Python {rt.get('python', UNKNOWN)} ({rt.get('python_impl', UNKNOWN)}, "
                f"{rt.get('platform', UNKNOWN)}{frozen}) | {pk}")
    except Exception:
        return UNKNOWN


def as_json(indent: int = 2) -> str:
    return json.dumps(get(), indent=indent, default=str)


if __name__ == "__main__":
    print(one_line())
    print(as_json())
