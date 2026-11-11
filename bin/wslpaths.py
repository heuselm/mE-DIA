"""WSL <-> Windows path handling, kept to one boundary.

Under WSL the server and the pipeline work in POSIX paths, but DIA-NN and
dia-nn-plotter are Windows executables that can only read Windows paths, and
users paste Windows paths (``E:\\data\\a.raw``) into the GUI and the CLI.
Paths are therefore made POSIX where they enter (CLI arguments, GUI input,
dialog results) and made Windows only where they leave, in the argument list
of a Windows executable. Nothing in between converts anything.

Every converter takes the ``wslpath`` runner as an argument so it can be
tested without WSL; by default it calls the real ``wslpath``.
"""
import os
import re
import subprocess
from typing import Callable, List, Optional

_WIN_PATH = re.compile(r"^[A-Za-z]:[\\/]")
_MNT_DRIVE = re.compile(r"^/mnt/[A-Za-z](/|$)")

Runner = Callable[[str, str], str]   # (flag, path) -> converted path


def is_wsl() -> bool:
    """True inside the Windows Subsystem for Linux."""
    if os.environ.get("WSL_DISTRO_NAME") or os.environ.get("WSL_INTEROP"):
        return True
    try:
        with open("/proc/sys/kernel/osrelease") as fh:
            return "microsoft" in fh.read().lower()
    except OSError:
        return False


def wslpath(flag: str, path: str) -> str:
    """The real converter: ``wslpath -u|-w PATH``."""
    out = subprocess.run(["wslpath", flag, path], capture_output=True,
                         text=True, encoding="utf-8", timeout=30)
    if out.returncode != 0:
        raise ValueError(f"wslpath {flag} {path!r} failed: "
                         f"{(out.stderr or out.stdout).strip()}")
    return out.stdout.strip()


def is_windows_path(path) -> bool:
    """``C:\\...`` or ``C:/...``, or a ``\\\\server\\share`` UNC path."""
    s = str(path)
    return bool(_WIN_PATH.match(s)) or s.startswith("\\\\")


def is_windows_drive_path(path) -> bool:
    """A POSIX path on a mounted Windows drive, ``/mnt/<letter>/...``."""
    return bool(_MNT_DRIVE.match(str(path)))


def is_windows_exe(path) -> bool:
    """A Windows executable: ``*.exe`` (any case), or anything on a Windows drive."""
    s = str(path)
    return s.lower().endswith(".exe") or is_windows_drive_path(s) or is_windows_path(s)


def to_posix(path, runner: Optional[Runner] = None):
    """Windows-form path -> POSIX; anything else is returned unchanged."""
    if path is None or not is_windows_path(path):
        return path
    return (runner or wslpath)("-u", str(path))


def to_windows(path, runner: Optional[Runner] = None) -> str:
    """POSIX path on a Windows drive -> Windows form; a Windows-form path is
    returned unchanged (wslpath -w would mangle it)."""
    if is_windows_path(path):
        return str(path)
    return (runner or wslpath)("-w", str(path))


def windows_argv(cmd: List, runner: Optional[Runner] = None) -> List[str]:
    """Argument list for a Windows executable started from WSL.

    The executable itself is made POSIX (WSL has to find it). Every argument
    that is an absolute POSIX path is converted to Windows form; a Windows
    program cannot read the Linux filesystem, so a path outside /mnt/<drive>/
    raises ValueError naming it, before anything runs. All other arguments
    (flags, numbers, cleavage rules, modification specs) are left alone.
    """
    exe = to_posix(str(cmd[0]), runner)
    out = [exe]
    for a in cmd[1:]:
        s = str(a)
        if s.startswith("/"):
            if not is_windows_drive_path(s):
                raise ValueError(
                    f"{s} is inside the Linux filesystem, which the Windows program "
                    f"{os.path.basename(exe)} cannot read. Keep raw files, FASTA, "
                    f"libraries and the output folder on a Windows drive (/mnt/c, "
                    f"/mnt/e, ...).")
            s = to_windows(s, runner)
        out.append(s)
    return out
