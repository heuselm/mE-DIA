#!/usr/bin/env python3
"""
mE-DIA Pipeline — single-click launcher
  python bin/start.py [--port 8765] [--no-browser] [--build-ui]

Installs Python deps, optionally builds the React UI, then starts the server.
"""
import argparse
import os
import subprocess
import sys
from pathlib import Path

BIN = Path(__file__).parent
UI  = BIN / "ui"
REQ = BIN / "requirements.txt"

sys.path.insert(0, str(BIN))
import wslpaths  # noqa: E402


def run(cmd, **kw):
    print(f"  $ {' '.join(cmd)}")
    subprocess.run(cmd, check=True, **kw)


def npm(*args, cwd):
    """Run npm. On Windows npm is npm.cmd and needs the shell, given a string;
    on POSIX a list with shell=True would run bare `npm` and drop the arguments."""
    if sys.platform == "win32":
        line = " ".join(("npm",) + args)
        print(f"  $ {line}")
        subprocess.run(line, cwd=cwd, check=True, shell=True)
    else:
        run(["npm", *args], cwd=cwd)


def rollup_platform_missing(ui=UI) -> bool:
    """node_modules installed by Windows npm has only rollup's win32 binaries;
    a Linux build then fails with "Cannot find module @rollup/rollup-linux-...".
    """
    if sys.platform == "win32":
        return False
    rollup = Path(ui) / "node_modules" / "@rollup"
    if not rollup.is_dir():
        return False
    names = [p.name for p in rollup.iterdir()]
    plat = "rollup-darwin-" if sys.platform == "darwin" else "rollup-linux-"
    return any(n.startswith("rollup-") for n in names) and \
        not any(n.startswith(plat) for n in names)


def build_ui_with_windows_npm(ui=UI):
    """Under WSL, build with Windows npm, which matches the installed modules."""
    if not wslpaths.is_windows_drive_path(str(Path(ui).resolve())):
        return False
    run(["cmd.exe", "/c", "npm run build"], cwd=str(ui))
    return True


ROLLUP_HELP = """
  The UI's node_modules were installed by Windows npm and hold no Linux rollup
  binary, so this build cannot run from here. Either
    - build with Windows npm: in a Windows terminal, cd bin\\ui && npm run build
      (from WSL: cmd.exe /c "npm run build", with the repository on a Windows drive), or
    - reinstall the modules from WSL: rm -rf bin/ui/node_modules && cd bin/ui && npm install
      (Windows npm then needs its own reinstall before building from Windows again).
"""


def ensure_python_deps():
    print("\n[1/3] Checking Python dependencies...")
    try:
        import fastapi, uvicorn
        print("  fastapi + uvicorn: OK")
    except ImportError:
        print("  Installing: fastapi uvicorn[standard] websockets pandas")
        run([sys.executable, "-m", "pip", "install",
             "fastapi", "uvicorn[standard]", "websockets", "pandas"])


def ensure_ui(force_build=False, ui=None):
    ui = ui or UI
    dist = ui / "dist" / "index.html"
    if dist.exists() and not force_build:
        print("\n[2/3] React UI already built — skipping  (use --build-ui to rebuild)")
        return
    print("\n[2/3] Building React UI...")
    if not (ui / "node_modules").exists():
        print("  npm install...")
        npm("install", cwd=str(ui))
    if rollup_platform_missing(ui):
        if wslpaths.is_wsl() and build_ui_with_windows_npm(ui):
            print("  React UI built (with Windows npm: node_modules are Windows').")
            return
        print(ROLLUP_HELP)
        sys.exit(1)
    print("  npm run build...")
    npm("run", "build", cwd=str(ui))
    print("  React UI built.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--port",       type=int, default=8765)
    ap.add_argument("--no-browser", action="store_true")
    ap.add_argument("--build-ui",   action="store_true",
                    help="Force rebuild of React UI")
    args = ap.parse_args()

    ensure_python_deps()
    ensure_ui(force_build=args.build_ui)

    print(f"\n[3/3] Starting server on port {args.port}...")
    server_cmd = [sys.executable, str(BIN / "server.py"),
                  "--port", str(args.port)]
    if args.no_browser:
        server_cmd.append("--no-browser")
    # execv does not flush, so anything still buffered (every message above,
    # whenever stdout is a pipe or a log file rather than a terminal) would be
    # lost with the old process image.
    sys.stdout.flush()
    sys.stderr.flush()
    os.execv(sys.executable, server_cmd)  # replace process


if __name__ == "__main__":
    main()
