"""Start AptaRank with one command.

Resolves an interpreter, checks the environment, and opens the dashboard in a
browser. Written to be run by double-clicking `start.bat` — so it explains
itself in plain language and never fails silently.

Three rules it keeps, deliberately:

* **It never uses `python` from PATH.** On this project's development machine
  that resolves to an MSYS2 build with none of the scientific packages, and the
  resulting import errors are impossible for a non-programmer to interpret.
* **It never installs without asking.** First launch does need to download and
  build packages; that is a reasonable thing to do, but not a reasonable thing
  to do quietly on someone's machine.
* **It never generates synthetic reference data.** Silently creating a
  placeholder corpus is the single easiest way for a non-expert to end up with
  convincing-looking results that mean nothing. If no real reference library is
  present, the dashboard says so and offers the demonstration as an explicit,
  clearly-labelled choice.
"""

from __future__ import annotations

import argparse
import os
import platform
import shutil
import subprocess
import sys
import time
import webbrowser
from datetime import datetime
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
VENV_DIR = REPO_ROOT / ".venv"
LOG_PATH = REPO_ROOT / "logs" / "setup.log"
DEFAULT_PORT = 8501
MIN_PYTHON = (3, 10)

#: Import name -> what the user loses without it.
REQUIREMENTS = {
    "aptarank": "the ranking pipeline itself",
    "RNA": "RNA folding (ViennaRNA)",
    "forgi": "structure element parsing",
    "ushuffle": "shuffled control sequences",
    "Bio": "protein structure parsing (Biopython)",
    "pandas": "tables",
    "scipy": "statistics",
    "streamlit": "the dashboard",
    "altair": "the charts",
    "matplotlib": "the evaluation figures",
}

BANNER = r"""
  ___        _        ___             _
 / _ \      | |      | _ \__ _ _ _ __| |__
| |_| |_ __ | |_ __ _|   / _` | ' \/ /  / /
|  _  | '_ \| __/ _` | |_\ \__,_|_||_\_\_\_\
|_| |_| .__/ \__\__,_|_|  Interpretable RNA aptamer ranking
      |_|
"""


def log(message: str, echo: bool = True) -> None:
    LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
    with LOG_PATH.open("a", encoding="utf-8") as handle:
        handle.write(f"{datetime.now().isoformat(timespec='seconds')}  {message}\n")
    if echo:
        print(message)


def venv_python() -> Path:
    return VENV_DIR / ("Scripts/python.exe" if os.name == "nt" else "bin/python")


def in_project_venv() -> bool:
    try:
        return Path(sys.executable).resolve() == venv_python().resolve()
    except OSError:
        return False


def ask(question: str, default_yes: bool = True) -> bool:
    if not sys.stdin or not sys.stdin.isatty():
        return default_yes         # launched without a console; proceed sanely
    suffix = "[Y/n]" if default_yes else "[y/N]"
    answer = input(f"\n{question} {suffix} ").strip().lower()
    if not answer:
        return default_yes
    return answer.startswith("y")


# -- interpreter resolution ---------------------------------------------


def find_base_python() -> Path | None:
    """A real CPython >= 3.10 to build the environment from.

    Never falls through to whatever `python` happens to be on PATH.
    """
    if os.name == "nt":
        launcher = shutil.which("py")
        if launcher:
            for version in ("-3.12", "-3.11", "-3.10"):
                try:
                    found = subprocess.run(
                        [launcher, version, "-c", "import sys; print(sys.executable)"],
                        capture_output=True, text=True, timeout=30,
                    )
                except (OSError, subprocess.SubprocessError):
                    continue
                if found.returncode == 0 and found.stdout.strip():
                    return Path(found.stdout.strip())

    current = Path(sys.executable)
    healthy = (
        sys.version_info >= MIN_PYTHON
        and "msys" not in str(current).lower()
        and "mingw" not in str(current).lower()
    )
    return current if healthy else None


def create_environment() -> Path:
    base = find_base_python()
    if base is None:
        print(
            "\n  AptaRank needs Python "
            f"{MIN_PYTHON[0]}.{MIN_PYTHON[1]} or newer, and could not find a "
            "suitable copy on this computer.\n\n"
            "  Install it from https://www.python.org/downloads/ — tick "
            '"Add python.exe to PATH" during setup — then run this again.\n'
        )
        raise SystemExit(2)

    log(f"Creating a local Python environment with {base}")
    print("  This takes a minute or two. Nothing outside this folder is changed.\n")
    subprocess.run([str(base), "-m", "venv", str(VENV_DIR)], check=True)
    return venv_python()


# -- environment doctor --------------------------------------------------


def missing_packages(python: Path) -> list[str]:
    """Which requirements cannot be imported by the given interpreter."""
    script = (
        "import importlib, json, sys\n"
        f"names = {list(REQUIREMENTS)!r}\n"
        "missing = []\n"
        "for name in names:\n"
        "    try:\n"
        "        importlib.import_module(name)\n"
        "    except Exception:\n"
        "        missing.append(name)\n"
        "print(json.dumps(missing))\n"
    )
    result = subprocess.run(
        [str(python), "-c", script], capture_output=True, text=True,
    )
    try:
        import json

        return json.loads(result.stdout.strip().splitlines()[-1])
    except (ValueError, IndexError):
        return list(REQUIREMENTS)


def install(python: Path) -> None:
    """Install AptaRank and its dependencies, including the ushuffle workaround."""
    log("Installing AptaRank and its dependencies")

    print("  Installing build tools…")
    result = subprocess.run(
        [str(python), "-m", "pip", "install", "--upgrade", "pip", "setuptools",
         "wheel", "cython"],
        capture_output=True, text=True,
    )
    log(f"build tools:\n{result.stdout[-2000:]}\n{result.stderr[-2000:]}", echo=False)
    if result.returncode != 0:
        print(f"\n  Could not install the build tools. Details in {LOG_PATH}")
        raise SystemExit(1)

    # ushuffle before the package that depends on it: installing AptaRank first
    # would make pip build ushuffle under build isolation, where it always
    # fails (see install_ushuffle).
    if "ushuffle" in missing_packages(python):
        install_ushuffle(python)

    print("  Installing AptaRank and its scientific packages…")
    result = subprocess.run(
        [str(python), "-m", "pip", "install", "-e", f"{REPO_ROOT}[dashboard]"],
        capture_output=True, text=True,
    )
    log(f"aptarank:\n{result.stdout[-4000:]}\n{result.stderr[-4000:]}", echo=False)
    if result.returncode != 0:
        print(f"\n  Could not install AptaRank. The details are in {LOG_PATH}")
        raise SystemExit(1)


def install_ushuffle(python: Path) -> None:
    """Build ushuffle from source.

    Its published package ships pre-generated Cython C that no longer compiles
    against modern CPython — `tp_print` on 3.9–3.11, `longintrepr.h` on 3.12 —
    so the normal install always fails. Building from the `.pyx` with a current
    Cython works, which means `--no-build-isolation` (build isolation hides the
    Cython we just installed). It needs a C compiler, the one prerequisite a
    biologist's machine may genuinely lack.
    """
    print("  Building the shuffling library from source…")
    workdir = REPO_ROOT / "logs" / "ushuffle_build"
    workdir.mkdir(parents=True, exist_ok=True)
    commands = [
        [str(python), "-m", "pip", "download", "ushuffle", "--no-binary", ":all:",
         "--no-deps", "-d", str(workdir)],
    ]
    for command in commands:
        result = subprocess.run(command, capture_output=True, text=True)
        log(f"$ {' '.join(command)}\n{result.stdout[-2000:]}\n{result.stderr[-2000:]}",
            echo=False)

    archives = list(workdir.glob("ushuffle-*.tar.gz"))
    if archives:
        import tarfile

        with tarfile.open(archives[0]) as tar:
            tar.extractall(workdir)
        sources = [p for p in workdir.glob("ushuffle-*") if p.is_dir()]
        if sources:
            result = subprocess.run(
                [str(python), "-m", "pip", "install", "--no-build-isolation",
                 str(sources[0])],
                capture_output=True, text=True,
            )
            log(f"ushuffle build:\n{result.stdout[-4000:]}\n{result.stderr[-4000:]}",
                echo=False)
            if result.returncode == 0:
                return

    print(
        "\n  The shuffling library could not be built. It needs a C compiler:\n"
        "    Windows — install 'Microsoft C++ Build Tools' from\n"
        "              https://visualstudio.microsoft.com/visual-cpp-build-tools/\n"
        "    macOS   — run: xcode-select --install\n"
        "    Linux   — install build-essential\n\n"
        "  AptaRank will still start, but the shuffled-control check — which "
        "shows a score reflects structure rather than nucleotide composition — "
        "will not be available.\n"
    )


# -- reference data ------------------------------------------------------


def report_reference_state() -> None:
    """Say what reference data exists. Never create any."""
    corpora = sorted((REPO_ROOT / "data" / "corpus").glob("*.csv"))
    real = [p for p in corpora if "placeholder" not in p.stem.lower()]
    targets = sorted((REPO_ROOT / "cache" / "targets").glob("*.bundle.json"))

    print("\n  Reference data on this computer:")
    if real:
        for path in real:
            print(f"    validated library   {path.name}")
    else:
        print("    validated library   none found")
    if corpora and not real:
        print("    example library     available (synthetic — testing only)")
    print(f"    prepared targets    {len(targets) if targets else 'none found'}")

    if not real:
        print(
            "\n  No validated reference library is installed yet. The dashboard "
            "will\n  ask you to upload one — or to start a clearly-labelled "
            "demonstration\n  using synthetic data, which cannot support any "
            "scientific claim.\n"
        )


# -- launching -----------------------------------------------------------


def launch(python: Path, port: int, open_browser: bool = True) -> int:
    app = REPO_ROOT / "dashboard" / "streamlit_app.py"
    command = [
        str(python), "-m", "streamlit", "run", str(app),
        "--server.port", str(port),
        "--server.headless", "true",
        "--browser.gatherUsageStats", "false",
    ]
    log(f"Launching: {' '.join(command)}", echo=False)
    print(f"\n  Starting AptaRank at http://localhost:{port}")
    print("  Leave this window open while you use it. Close it to stop.\n")

    process = subprocess.Popen(command, cwd=str(REPO_ROOT))
    if open_browser:
        time.sleep(4)
        webbrowser.open(f"http://localhost:{port}")
    try:
        return process.wait()
    except KeyboardInterrupt:
        process.terminate()
        return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Start the AptaRank dashboard.")
    parser.add_argument("--port", type=int, default=DEFAULT_PORT)
    parser.add_argument("--no-browser", action="store_true")
    parser.add_argument("--check-only", action="store_true",
                        help="run the environment check and exit")
    parser.add_argument("--assume-yes", action="store_true",
                        help="skip the install confirmation prompt")
    args = parser.parse_args(argv)

    print(BANNER)
    log(f"start.py on {platform.platform()} with {sys.executable}", echo=False)

    # Step into the project's own environment before doing anything else.
    if not in_project_venv():
        python = venv_python()
        if not python.exists():
            print("  First launch: AptaRank needs to set up a local Python")
            print("  environment inside this folder and install its scientific")
            print("  packages. This downloads software and takes a few minutes.")
            if not (args.assume_yes or ask("Set it up now?")):
                print("\n  Nothing was changed.")
                return 1
            python = create_environment()
        # Re-run this script inside the environment, so everything below runs
        # with the right interpreter.
        forwarded = [a for a in (argv or sys.argv[1:])]
        if args.assume_yes and "--assume-yes" not in forwarded:
            forwarded.append("--assume-yes")
        return subprocess.run([str(python), str(Path(__file__).resolve()), *forwarded]).returncode

    python = Path(sys.executable)
    missing = missing_packages(python)
    if missing:
        readable = ", ".join(f"{name} ({REQUIREMENTS[name]})" for name in missing)
        print(f"  Missing components: {readable}")
        if not (args.assume_yes or ask("Install them now?")):
            print("\n  Nothing was changed.")
            return 1
        install(python)
        missing = missing_packages(python)
        if [m for m in missing if m != "ushuffle"]:
            print(
                f"\n  Setup did not complete: {', '.join(missing)} still unavailable."
                f"\n  Details are in {LOG_PATH}. AptaRank will not start in this state."
            )
            return 1
    print("  Environment ready.")

    report_reference_state()
    if args.check_only:
        return 0
    return launch(python, args.port, open_browser=not args.no_browser)


if __name__ == "__main__":
    raise SystemExit(main())
