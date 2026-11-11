#!/usr/bin/env python3
"""
Assemble a portable Python folder with pandas + matplotlib pre-installed.

  python bin/build_portable_python.py [--out PATH] [--py-version 3.11.9]

Downloads the official Windows embeddable Python distribution (zip), unpacks
it, bootstraps pip, and installs pandas + matplotlib. The resulting folder is
self-contained: copy it next to mEDIA.exe and the FastAPI server will spawn
run_pipeline.py via <portable>/python.exe (see server.py:_resolve_python).

Idempotent: re-running with an existing target reuses python.exe if the
expected packages import correctly.

Notes
-----
* Embeddable Python ships without the standard library on sys.path. We patch
  ``pythonXX._pth`` to uncomment ``import site`` so pip + installed packages
  become importable.
* This is Windows-only; on non-Windows the script exits without action — the
  fall-through to system Python via PATH covers POSIX cases.
"""
import argparse
import os
import shutil
import subprocess
import sys
import urllib.request
import zipfile
from pathlib import Path

DEFAULT_VERSION = "3.11.9"
PACKAGES = ["pandas", "matplotlib"]


def download(url: str, dst: Path) -> None:
    print(f"  GET {url}")
    with urllib.request.urlopen(url) as r, open(dst, "wb") as f:
        shutil.copyfileobj(r, f)


def patch_pth(out: Path, py_minor: str) -> None:
    """Embeddable Python's pythonXX._pth disables site.py. Re-enable it so
    pip + installed packages become importable."""
    pth = next(out.glob(f"python{py_minor}._pth"), None)
    if not pth:
        return
    text = pth.read_text(encoding="utf-8")
    if "#import site" in text:
        pth.write_text(text.replace("#import site", "import site"), encoding="utf-8")
        print(f"  patched {pth.name}: enabled site")


def verify(out: Path) -> bool:
    py = out / "python.exe"
    if not py.exists():
        return False
    code = "; ".join([f"import {p}" for p in PACKAGES])
    try:
        subprocess.run([str(py), "-c", code], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return True
    except subprocess.CalledProcessError:
        return False


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(Path("dist") / "python_portable"),
                    help="Target folder for the portable Python")
    ap.add_argument("--py-version", default=DEFAULT_VERSION,
                    help=f"Python version to download (default: {DEFAULT_VERSION})")
    ap.add_argument("--force", action="store_true",
                    help="Rebuild even if the target verifies cleanly")
    args = ap.parse_args()

    if os.name != "nt":
        print("Not on Windows — portable Python is Windows-only. "
              "POSIX users get system Python via PATH.")
        return 0

    out = Path(args.out).resolve()

    if out.exists() and verify(out) and not args.force:
        print(f"[OK] {out} already has python.exe + pandas + matplotlib. Skipping.")
        return 0

    if out.exists() and args.force:
        print(f"[..] --force: removing existing {out}")
        shutil.rmtree(out)
    out.mkdir(parents=True, exist_ok=True)

    py_minor = "".join(args.py_version.split(".")[:2])  # "3.11.9" -> "311"
    url = (f"https://www.python.org/ftp/python/{args.py_version}/"
           f"python-{args.py_version}-embed-amd64.zip")
    zip_path = out / "embed.zip"

    print(f"[1/4] Download embeddable Python {args.py_version}")
    download(url, zip_path)

    print(f"[2/4] Unpack into {out}")
    with zipfile.ZipFile(zip_path) as zf:
        zf.extractall(out)
    zip_path.unlink()
    patch_pth(out, py_minor)

    print(f"[3/4] Bootstrap pip")
    getpip = out / "get-pip.py"
    download("https://bootstrap.pypa.io/get-pip.py", getpip)
    subprocess.run([str(out / "python.exe"), str(getpip), "--no-warn-script-location"],
                   check=True)
    getpip.unlink()

    print(f"[4/4] Install {', '.join(PACKAGES)}")
    subprocess.run([str(out / "python.exe"), "-m", "pip", "install",
                    "--no-warn-script-location", *PACKAGES],
                   check=True)

    if verify(out):
        print(f"\n[OK] Portable Python ready: {out}")
        return 0
    print("\n[ERR] Verification failed — pandas/matplotlib import did not succeed.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
