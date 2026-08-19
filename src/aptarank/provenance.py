"""Provenance: tool versions, content hashes and deterministic seed derivation.

Every claim the paper makes is traceable to a run artifact, and every run
artifact has to say exactly what produced it.
"""

from __future__ import annotations

import contextlib
import hashlib
import io
import json
import re
import platform
import struct
import subprocess
import sys
from pathlib import Path
from typing import Any

from . import ARTIFACT_SCHEMA_VERSION, __version__


def sha256_file(path: str | Path, chunk: int = 1 << 20) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as fh:
        while block := fh.read(chunk):
            digest.update(block)
    return digest.hexdigest()


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


#: Sidecar written beside a staged upload, recording what the user called it.
#: Uploads are stored under a content hash — two people uploading the same file
#: must not collide, and a client-supplied name must never be used as a path —
#: but showing that hash back to the user as "your input file" is a provenance
#: record of a name they have never seen. Both are kept: the hash for integrity,
#: the original name for the human.
ORIGIN_SUFFIX = ".origin.json"


def write_origin(path: str | Path, original_filename: str, **extra: Any) -> Path:
    target = Path(str(path) + ORIGIN_SUFFIX)
    payload = {"original_filename": original_filename, **extra}
    target.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return target


def original_filename(path: str | Path) -> str | None:
    """What the user called this file, if anything recorded it."""
    sidecar = Path(str(path) + ORIGIN_SUFFIX)
    if not sidecar.is_file():
        return None
    try:
        data = json.loads(sidecar.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return None
    name = data.get("original_filename") if isinstance(data, dict) else None
    return str(name) if name else None


def derive_seed(run_seed: int, *parts: Any) -> int:
    """A stable 32-bit seed from the run seed plus arbitrary identifying parts.

    Used so that per-candidate randomness (shuffling, stochastic sampling) does
    not depend on process-pool scheduling order. Calling this with the same
    arguments always returns the same seed, on any platform.
    """
    payload = "|".join(str(p) for p in (run_seed, *parts)).encode("utf-8")
    digest = hashlib.blake2b(payload, digest_size=4).digest()
    return struct.unpack("<I", digest)[0]


def tool_versions() -> dict[str, str | None]:
    """Versions of everything whose behaviour can change a number we report."""
    versions: dict[str, str | None] = {
        "aptarank": __version__,
        "artifact_schema": ARTIFACT_SCHEMA_VERSION,
        "python": sys.version.split()[0],
        "platform": f"{platform.system()} {platform.release()}",
    }
    versions["viennarna"] = _module_version("RNA")
    versions["forgi"] = _module_version("forgi")
    versions["ushuffle"] = _module_version("ushuffle")
    versions["biopython"] = _module_version("Bio")
    versions["numpy"] = _module_version("numpy")
    versions["scipy"] = _module_version("scipy")
    versions["fpocket"] = _cli_version(["fpocket", "--version"])
    versions["pdb2pqr"] = _cli_version([sys.executable, "-m", "pdb2pqr", "--version"])
    versions["apbs"] = _cli_version(["apbs", "--version"])
    return versions


def git_state(repo_root: str | Path) -> dict[str, Any]:
    """Commit and dirty flag, or nulls when this is not a git checkout."""
    root = Path(repo_root)
    try:
        commit = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=root, capture_output=True, text=True, timeout=10,
        )
        if commit.returncode != 0:
            return {"commit": None, "dirty": None}
        status = subprocess.run(
            ["git", "status", "--porcelain"],
            cwd=root, capture_output=True, text=True, timeout=10,
        )
        return {
            "commit": commit.stdout.strip(),
            "dirty": bool(status.stdout.strip()),
        }
    except (OSError, subprocess.SubprocessError):
        return {"commit": None, "dirty": None}


def _module_version(name: str) -> str | None:
    # forgi's __init__ pulls in an optional 3D module compiled against NumPy 1.x,
    # which prints a multi-line traceback to stderr before forgi swallows the
    # error. tier1/elements.py silences it at its own import; this path imports
    # forgi independently, and an unexplained traceback in a job log reads as a
    # crash to whoever is looking at it.
    try:
        with contextlib.redirect_stderr(io.StringIO()):
            module = __import__(name)
    except Exception:
        return None
    for attr in ("__version__", "version", "VERSION"):
        value = getattr(module, attr, None)
        if isinstance(value, str):
            return value
    return "unknown"


#: A version line has a number in it. fpocket answers `--version` with a
#: pocket-hunting banner and an error about the missing input file before
#: naming itself; APBS puts its banner on stdout and its version on stderr.
#: Taking "the first line of whichever stream spoke" recorded
#: "***** POCKET HUNTING BEGINS *****" as the fpocket version — in the field
#: whose whole job is to say what produced these numbers.
_VERSION_LINE = re.compile(r"\d+\.\d+")
_ANSI = re.compile(r"\x1b\[[0-9;]*m")


def _cli_version(cmd: list[str]) -> str | None:
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=15)
    except (OSError, subprocess.SubprocessError):
        return None
    for stream in (proc.stdout, proc.stderr):
        for line in (stream or "").splitlines():
            cleaned = _ANSI.sub("", line).strip().strip(":|").strip()
            if cleaned and _VERSION_LINE.search(cleaned):
                return cleaned
    output = ((proc.stdout or "") + (proc.stderr or "")).strip()
    return output.splitlines()[0] if output else None
