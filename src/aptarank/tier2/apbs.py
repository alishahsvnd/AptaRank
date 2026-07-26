"""Electrostatic compatibility (spec §5.7). Stretch goal, target-level only.

RNA carries a uniformly negative backbone charge, so a positively charged
cavity is more electrostatically hospitable to *any* RNA. Because that depends
only on the target, it is a single target-level flag — folding it into a
per-candidate score would add the same constant to every candidate, changing
no ordering while making the score harder to explain.

Deliberately never allowed to block anything: PDB2PQR routinely struggles with
metals and non-standard residues, which is precisely the case for the demo
targets. Any failure is recorded with a reason and the bundle is still valid.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import Any

import numpy as np

from .fpocket import Pocket


def compute(structure_path: str | Path, pocket: Pocket, work_dir: str | Path) -> dict[str, Any]:
    """Run PDB2PQR + APBS and sample the grid at the pocket's alpha spheres."""
    path, work = Path(structure_path), Path(work_dir)
    work.mkdir(parents=True, exist_ok=True)

    for tool in ("pdb2pqr30", "apbs"):
        if shutil.which(tool) is None:
            return _skipped(f"{tool} not found on PATH")

    pqr = work / f"{path.stem}.pqr"
    apbs_in = work / f"{path.stem}.in"
    commands = [
        ["pdb2pqr30", "--ff=AMBER", f"--apbs-input={apbs_in}", str(path), str(pqr)],
        ["apbs", str(apbs_in)],
    ]
    executions = []
    for command in commands:
        proc = subprocess.run(command, capture_output=True, text=True, cwd=work)
        executions.append(
            {"command": command, "exit_code": proc.returncode,
             "status": "success" if proc.returncode == 0 else "failed"}
        )
        if proc.returncode != 0:
            return {
                **_failed(f"{command[0]} exited {proc.returncode}: "
                          f"{(proc.stderr or proc.stdout).strip()[:300]}"),
                "pdb2pqr": executions[0],
                "apbs": executions[1] if len(executions) > 1 else None,
            }

    dx_files = sorted(work.glob("*.dx"))
    if not dx_files:
        return _failed("APBS produced no OpenDX potential grid")

    try:
        grid = read_opendx(dx_files[0])
        centres = np.array([s.center_A for s in pocket.alpha_spheres], dtype=float)
        values = sample_grid(grid, centres)
    except Exception as exc:  # noqa: BLE001 - never let this break a bundle
        return _failed(f"could not sample the potential grid: {exc}")

    inside = values[np.isfinite(values)]
    if inside.size == 0:
        return _failed("no alpha-sphere centre fell inside the APBS grid")

    mean = float(inside.mean())
    return {
        "requested": True,
        "status": "success",
        "reason_code": None,
        "message": None,
        "pdb2pqr": executions[0],
        "apbs": executions[1],
        "grid": {
            "units": "kT/e",
            "origin_A": [float(v) for v in grid["origin"]],
            "spacing_A": [float(v) for v in grid["spacing"]],
            "counts": [int(v) for v in grid["counts"]],
            "file": dx_files[0].name,
        },
        "selected_pocket_sampling": {
            "n_points_requested": int(centres.shape[0]),
            "n_points_inside_grid": int(inside.size),
            "mean_potential_kT_per_e": mean,
            "median_potential_kT_per_e": float(np.median(inside)),
            # Positive mean potential = hospitable to a negatively charged
            # backbone. A target-level badge, never part of a candidate score.
            "electrostatic_compatible": bool(mean > 0),
        },
    }


def read_opendx(path: str | Path) -> dict[str, Any]:
    """Minimal OpenDX reader for the APBS potential grid."""
    counts = origin = None
    deltas: list[list[float]] = []
    values: list[float] = []
    with Path(path).open("r", encoding="utf-8") as fh:
        for line in fh:
            stripped = line.strip()
            if stripped.startswith("object 1 class gridpositions counts"):
                counts = [int(v) for v in stripped.split()[-3:]]
            elif stripped.startswith("origin"):
                origin = [float(v) for v in stripped.split()[1:4]]
            elif stripped.startswith("delta"):
                deltas.append([float(v) for v in stripped.split()[1:4]])
            elif stripped and not stripped[0].isalpha() and not stripped.startswith("#"):
                values.extend(float(v) for v in stripped.split())
    if counts is None or origin is None or len(deltas) < 3:
        raise ValueError("incomplete OpenDX header")
    data = np.array(values, dtype=float)
    if data.size < counts[0] * counts[1] * counts[2]:
        raise ValueError(
            f"OpenDX grid has {data.size} values, expected {np.prod(counts)}"
        )
    return {
        "counts": counts,
        "origin": origin,
        "spacing": [deltas[0][0], deltas[1][1], deltas[2][2]],
        "data": data[: int(np.prod(counts))].reshape(counts),
    }


def sample_grid(grid: dict[str, Any], points: np.ndarray) -> np.ndarray:
    """Trilinear interpolation of the potential at arbitrary coordinates."""
    origin = np.array(grid["origin"], dtype=float)
    spacing = np.array(grid["spacing"], dtype=float)
    counts = np.array(grid["counts"], dtype=int)
    data = grid["data"]

    fractional = (points - origin) / spacing
    out = np.full(points.shape[0], np.nan)
    for i, f in enumerate(fractional):
        base = np.floor(f).astype(int)
        if np.any(base < 0) or np.any(base + 1 >= counts):
            continue  # outside the grid; reported as n_points_inside_grid
        frac = f - base
        acc = 0.0
        for dx in (0, 1):
            for dy in (0, 1):
                for dz in (0, 1):
                    weight = (
                        (frac[0] if dx else 1 - frac[0])
                        * (frac[1] if dy else 1 - frac[1])
                        * (frac[2] if dz else 1 - frac[2])
                    )
                    acc += weight * data[base[0] + dx, base[1] + dy, base[2] + dz]
        out[i] = acc
    return out


def _skipped(message: str) -> dict[str, Any]:
    return {
        "requested": True, "status": "skipped", "reason_code": "tool_unavailable",
        "message": message, "pdb2pqr": None, "apbs": None, "grid": None,
        "selected_pocket_sampling": None,
    }


def _failed(message: str) -> dict[str, Any]:
    return {
        "requested": True, "status": "failed", "reason_code": "execution_failed",
        "message": message, "pdb2pqr": None, "apbs": None, "grid": None,
        "selected_pocket_sampling": None,
    }
