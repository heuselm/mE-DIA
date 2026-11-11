#!/usr/bin/env python3
"""Checks for build provenance. Run: python bin/test_build_info.py"""
import os
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import build_info as bi  # noqa: E402

FAIL = []


def check(name, cond, detail=""):
    print(f"  {'PASS' if cond else 'FAIL'}  {name}{('  ' + str(detail)) if detail else ''}")
    if not cond:
        FAIL.append(name)
        if "PYTEST_CURRENT_TEST" in os.environ:     # otherwise pytest never sees a failure
            raise AssertionError(f"{name} {detail}")


def git(*args, cwd):
    return subprocess.run(("git",) + args, cwd=cwd, capture_output=True, text=True)


def main():
    print("artifact names cannot lie")
    n = bi.artifact_name({"branch": "main", "commit": "abc1234", "dirty": False})
    check("clean build: mEDIA-main-abc1234.exe", n == "mEDIA-main-abc1234.exe", n)
    n = bi.artifact_name({"branch": "main", "commit": "abc1234", "dirty": True})
    check("dirty build is marked", n == "mEDIA-main-abc1234-dirty.exe", n)
    n = bi.artifact_name({"branch": "feature/ms1-rt/x", "commit": "abc1234"})
    check("a slash in the branch cannot create a directory",
          "/" not in n and n == "mEDIA-feature-ms1-rt-x-abc1234.exe", n)
    n = bi.artifact_name({"branch": "", "commit": ""})
    check("empty provenance still yields a usable name",
          n == "mEDIA-unknown-unknown.exe", n)
    n = bi.artifact_name({"branch": "main", "commit": "abc1234"}, ext="")
    check("extension is optional (linux build)", n == "mEDIA-main-abc1234", n)

    print("reading provenance from git")
    with tempfile.TemporaryDirectory() as td:
        check("outside a repository: unknown, never an exception",
              bi.from_git(cwd=td).get("commit") == bi.UNKNOWN)
        for cmd in (("init", "-q", "-b", "trunk"), ("config", "user.email", "t@t"),
                    ("config", "user.name", "t")):
            git(*cmd, cwd=td)
        open(os.path.join(td, "a.txt"), "w").write("1\n")
        git("add", "a.txt", cwd=td)
        git("commit", "-qm", "first", cwd=td)
        info = bi.from_git(cwd=td)
        head = git("rev-parse", "--short=7", "HEAD", cwd=td).stdout.strip()
        check("clean tree: commit and branch are read",
              info["commit"] == head and info["branch"] == "trunk", info)
        check("clean tree is not dirty", info["dirty"] is False)
        open(os.path.join(td, "a.txt"), "w").write("2\n")
        check("uncommitted change makes it dirty", bi.from_git(cwd=td)["dirty"] is True)
        check("name of a dirty build is marked",
              bi.artifact_name(bi.from_git(cwd=td)).endswith("-dirty.exe"))

    print("one-line summary")
    line = bi.one_line({"branch": "main", "commit": "abc1234", "dirty": False,
                        "built_utc": "2026-10-08T12:00:00Z"})
    check("names branch, commit and build time",
          "main@abc1234" in line and "2026-10-08" in line, line)
    line = bi.one_line({"branch": "main", "commit": "abc1234", "dirty": True})
    check("says DIRTY when it is", "DIRTY" in line, line)

    print("which source wins (#27)")
    STALE = {"commit": "626ff0f", "branch": "main", "dirty": False, "source": "git",
             "built_utc": "2026-10-08T20:19:50Z"}
    LIVE = {"commit": "abc1234", "branch": "feat/x", "dirty": True, "source": "git"}
    NOGIT = {"commit": bi.UNKNOWN, "branch": bi.UNKNOWN, "dirty": False,
             "describe": bi.UNKNOWN, "source": "no-git"}
    i = bi.get(frozen=False, stamp=lambda: dict(STALE), git=lambda: dict(LIVE))
    check("source checkout: live git beats a stale stamp",
          i["commit"] == "abc1234" and i["branch"] == "feat/x" and i["source"] == "git", i)
    i = bi.get(frozen=True, stamp=lambda: dict(STALE), git=lambda: dict(LIVE))
    check("frozen build: the stamp is used", i["commit"] == "626ff0f" and i["frozen"] is True, i)
    check("a stamp is labelled stamped, although it carries source 'git'",
          i["source"] == "stamped", i)
    i = bi.get(frozen=False, stamp=lambda: dict(STALE), git=lambda: dict(NOGIT))
    check("source tree without git: falls back to the stamp, labelled stamped",
          i["commit"] == "626ff0f" and i["source"] == "stamped", i)
    i = bi.get(frozen=False, stamp=lambda: None, git=lambda: dict(NOGIT))
    check("neither: unknown, never an exception", i["commit"] == bi.UNKNOWN, i)
    i = bi.get(frozen=True, stamp=lambda: None, git=lambda: dict(LIVE))
    check("frozen without a stamp: whatever git says", i["commit"] == "abc1234", i)
    check("dirty live tree is not described as 'at build time'",
          "at build time" not in bi.one_line(dict(LIVE)), bi.one_line(dict(LIVE)))
    check("dirty stamp still says 'at build time'",
          "at build time" in bi.one_line({**STALE, "dirty": True, "source": "stamped"}))

    print("this checkout")
    info = bi.get()
    head = git("rev-parse", "--short=7", "HEAD", cwd=os.path.dirname(os.path.abspath(__file__)))
    if head.returncode == 0:
        # the regression: with a bin/_build_info.py from an earlier exe build present,
        # get() used to return that build's commit instead of HEAD
        check("get() reports this working tree's HEAD, not an earlier build",
              info.get("commit") == head.stdout.strip() and info.get("source") == "git",
              bi.one_line(info))
    else:
        check("get() returns a commit for this working tree",
              isinstance(info.get("commit"), str) and info.get("commit") != "",
              bi.one_line(info))
    check("as_json is valid JSON", __import__("json").loads(bi.as_json()) is not None)

    print(f"\n{'all checks passed' if not FAIL else f'{len(FAIL)} FAILED: {FAIL}'}")
    return 1 if FAIL else 0


def test_build_info():
    """Entry point for `pytest bin/`; the checks themselves live in main()."""
    FAIL.clear()
    assert main() == 0, f"failed checks: {FAIL}"


if __name__ == "__main__":
    sys.exit(main())
