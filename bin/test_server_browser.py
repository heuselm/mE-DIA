#!/usr/bin/env python3
"""Checks for opening the GUI in a browser. Run: python bin/test_server_browser.py

The WSL case is the one that regressed: a distro has no browser, so
webbrowser/xdg-open fail and the URL has to be handed to Windows.
"""
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import server as sv  # noqa: E402

FAIL = []
URL = "http://localhost:8765"


def check(name, cond):
    print(f"  {'PASS' if cond else 'FAIL'}  {name}")
    if not cond:
        FAIL.append(name)
        # Without this a failing check only sets an exit code, which pytest never
        # sees, so the check would pass silently under `pytest bin/`.
        if "PYTEST_CURRENT_TEST" in os.environ:
            raise AssertionError(name)


class Fake:
    """Stands in for subprocess.run: records calls, succeeds for `ok_exe`."""

    def __init__(self, ok_exe=None, missing=()):
        self.ok_exe, self.missing, self.calls = ok_exe, set(missing), []

    def run(self, cmd, **kw):
        self.calls.append(cmd[0])
        if cmd[0] in self.missing:
            raise FileNotFoundError(cmd[0])
        rc = 0 if cmd[0] == self.ok_exe else 1
        return subprocess.CompletedProcess(cmd, rc)


STARTED = []   # every call that would have reached os.startfile


def as_if(platform, wsl, fake=None, webbrowser_open=None, startfile=None):
    """Call open_browser() with everything outside the process faked out.

    Platform, WSL detection, subprocess.run, os.startfile and webbrowser.open
    are all replaced, so the scenario is the same on every interpreter and no
    test ever opens a real browser (#40: under Windows python.exe the real
    sys.platform used to win, and os.startfile opened a tab).
    """
    saved = (sv._platform, sv._is_wsl, sv.subprocess.run, sv._startfile, sv.webbrowser.open)
    sv._platform = lambda: platform
    sv._is_wsl = lambda: wsl
    sv.subprocess.run = (fake or Fake(missing=["wslview", "powershell.exe", "cmd.exe",
                                               "explorer.exe"])).run
    sv._startfile = startfile or STARTED.append
    sv.webbrowser.open = webbrowser_open or (lambda u: False)
    try:
        return sv.open_browser(URL)
    finally:
        (sv._platform, sv._is_wsl, sv.subprocess.run,
         sv._startfile, sv.webbrowser.open) = saved


def wsl_open(fake, monkey_webbrowser=None):
    """Call open_browser() as if on WSL (a Linux interpreter with interop)."""
    return as_if("linux", True, fake=fake, webbrowser_open=monkey_webbrowser)


def main():
    STARTED.clear()
    print("\nWSL detection")
    on_wsl = "microsoft" in open("/proc/sys/kernel/osrelease").read().lower() \
        if os.path.exists("/proc/sys/kernel/osrelease") else False
    check(f"_is_wsl() agrees with /proc/sys/kernel/osrelease ({on_wsl})",
          sv._is_wsl() == (on_wsl or bool(os.environ.get("WSL_DISTRO_NAME"))))

    print("\nlauncher chain on WSL")
    f = Fake(ok_exe="wslview")
    check("wslview is tried first and short-circuits", wsl_open(f) and f.calls == ["wslview"])

    f = Fake(ok_exe="powershell.exe", missing=["wslview"])
    check("no wslview -> powershell.exe",
          wsl_open(f) and f.calls == ["wslview", "powershell.exe"])

    f = Fake(ok_exe="cmd.exe", missing=["wslview"])
    check("powershell failing -> cmd.exe",
          wsl_open(f) and f.calls == ["wslview", "powershell.exe", "cmd.exe"])

    f = Fake(ok_exe=None, missing=["wslview"])
    check("all non-zero -> explorer.exe, reported as success",
          wsl_open(f) and f.calls[-1] == "explorer.exe")

    f = Fake(ok_exe=None, missing=["wslview", "powershell.exe", "cmd.exe", "explorer.exe"])
    check("nothing available -> False, so main() can print the URL",
          wsl_open(f) is False)

    f = Fake(ok_exe=None, missing=["wslview", "powershell.exe", "cmd.exe", "explorer.exe"])
    check("WSLg fallback: webbrowser still consulted last",
          wsl_open(f, monkey_webbrowser=lambda u: True) is True)

    print("\nthe URL that is handed over")
    f = Fake(ok_exe="powershell.exe", missing=["wslview"])
    seen = []
    f.run = (lambda cmd, **kw: (seen.append(cmd), Fake.run(f, cmd, **kw))[1])
    wsl_open(f)
    ps = [c for c in seen if c[0] == "powershell.exe"][0]
    check("powershell gets Start-Process with the quoted URL",
          ps[-1] == f"Start-Process '{URL}'")

    print("\nnon-WSL Linux")
    hits = []
    ok = as_if("linux", False, webbrowser_open=lambda u: (hits.append(u), True)[1])
    check("plain Linux still goes through webbrowser", ok and hits == [URL])

    check("WSL and plain-Linux scenarios never reach os.startfile", STARTED == [])

    print("\nWindows")
    hits = []; f = Fake()
    ok = as_if("win32", False, fake=f, webbrowser_open=lambda u: (hits.append(u), True)[1])
    check("win32 -> os.startfile with the URL, nothing else",
          ok and STARTED == [URL] and f.calls == [] and hits == [])

    def refuse(u):
        raise OSError("no association")
    hits = []
    ok = as_if("win32", False, webbrowser_open=lambda u: (hits.append(u), True)[1],
               startfile=refuse)
    check("win32, os.startfile failing -> falls back to webbrowser", ok and hits == [URL])


    print(f"\n{'all checks passed' if not FAIL else f'{len(FAIL)} FAILED: {FAIL}'}")
    return 1 if FAIL else 0


def test_browser_launch():
    """Entry point for `pytest bin/`; the checks themselves live in main().

    Without a test_* function pytest collects nothing from this file and reports
    success, leaving the WSL browser-launch logic uncovered.
    """
    FAIL.clear()
    assert main() == 0, f"failed checks: {FAIL}"


if __name__ == "__main__":
    sys.exit(main())
