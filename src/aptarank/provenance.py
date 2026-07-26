"""Provenance: tool versions, content hashes and deterministic seed derivation.

Every claim the paper makes is traceable to a run artifact, and every run
artifact has to say exactly what produced it.
"""

from __future__ import annotations

import hashlib
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
    versions["pdb2pqr"] = _cli_version(["pdb2pqr30", "--version"])
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
    try:
        module = __import__(name)
    except Exception:
        return None
    for attr in ("__version__", "version", "VERSION"):
        value = getattr(module, attr, None)
        if isinstance(value, str):
            return value
    return "unknown"


def _cli_version(cmd: list[str]) -> str | None:
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=15)
    except (OSError, subprocess.SubprocessError):
        return None
    output = (proc.stdout or proc.stderr).strip()
    return output.splitlines()[0] if output else None
