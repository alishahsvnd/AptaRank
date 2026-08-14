"""Electrostatics (spec §5.7, refinements §5.2).

RNA carries a uniformly negative backbone charge, so a positively charged region
is more electrostatically hospitable to *any* RNA. That is a property of the
target alone, which has one important consequence: the signal is identical for
every candidate. It can shift what a score means, but it can never move one
candidate past another, and it can never change a band. Surface mode therefore
carries it as a documented target-level term of the composite rather than
pretending it discriminates between candidates.

Deliberately never allowed to block anything: PDB2PQR routinely struggles with
metals and non-standard residues, which is precisely the case for the demo
targets. Any failure is recorded with a reason and the bundle is still valid.
"""

from __future__ import annotations

import math
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any, Sequence

import numpy as np

from .fpocket import Pocket

#: Points are sampled this far outside the protein atoms, along the outward
#: radial direction. Sampling at an atom centre reads the potential inside the
#: low-dielectric interior, which is not what a solvated aptamer approaching the
#: surface would experience.
SURFACE_PROBE_OFFSET_A = 3.0


def pdb2pqr_command() -> list[str] | None:
    """How to invoke PDB2PQR here, preferring the interpreter we are running in.

    `pdb2pqr30` lives in the virtualenv's bin directory, which is only on PATH
    when the environment has been activated — and the pipeline is deliberately
    launched by absolute interpreter path, not through an activated shell.
    """
    module = subprocess.run(
        [sys.executable, "-m", "pdb2pqr", "--version"],
        capture_output=True, text=True, check=False,
    )
    if module.returncode == 0:
        return [sys.executable, "-m", "pdb2pqr"]
    executable = shutil.which("pdb2pqr30") or shutil.which("pdb2pqr")
    return [executable] if executable else None


def compute(
    structure_path: str | Path,
    points: Pocket | Sequence[Sequence[float]],
    work_dir: str | Path,
    label: str = "selected_pocket",
) -> dict[str, Any]:
    """Run PDB2PQR + APBS and sample the potential grid at the given points.

    `points` may be a Pocket (sampled at its alpha-sphere centres, as in pocket
    mode) or an explicit array of coordinates (surface mode's patch).
    """
    # Both tools run with cwd set to the work directory (APBS writes its grid
    # relative to cwd), so a relative input path would resolve against the wrong
    # place and fail with a confusing "file not found" from inside pdb2pqr.
    path, work = Path(structure_path).resolve(), Path(work_dir).resolve()
    work.mkdir(parents=True, exist_ok=True)

    centres = (
        np.asarray([s.center_A for s in points.alpha_spheres], dtype=float)
        if isinstance(points, Pocket)
        else np.asarray(points, dtype=float)
    )
    if centres.ndim != 2 or centres.shape[1] != 3 or centres.shape[0] == 0:
        return _failed(f"no valid sampling points for {label}")

    pdb2pqr = pdb2pqr_command()
    if pdb2pqr is None:
        return _skipped("pdb2pqr not found (neither as a module nor on PATH)")
    if shutil.which("apbs") is None:
        return _skipped("apbs not found on PATH")

    pqr = work / f"{path.stem}.pqr"
    apbs_in = work / f"{path.stem}.in"
    commands = [
        [*pdb2pqr, "--ff=AMBER", f"--apbs-input={apbs_in}", str(path), str(pqr)],
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
        values = sample_grid(grid, centres)
    except Exception as exc:  # noqa: BLE001 - never let this break a bundle
        return _failed(f"could not sample the potential grid: {exc}")

    inside = values[np.isfinite(values)]
    if inside.size == 0:
        return _failed(f"no {label} sampling point fell inside the APBS grid")

    mean = float(inside.mean())
    sampling = {
        "label": label,
        "n_points_requested": int(centres.shape[0]),
        "n_points_inside_grid": int(inside.size),
        "mean_potential_kT_per_e": mean,
        "median_potential_kT_per_e": float(np.median(inside)),
        # Positive mean potential = hospitable to a negatively charged backbone.
        # Target-level: the same for every candidate.
        "electrostatic_compatible": bool(mean > 0),
    }
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
        "sampling": sampling,
        # Kept under its original name too: v1 bundles and their readers use it.
        "selected_pocket_sampling": sampling,
    }


def surface_sample_points(
    patch_coords: Sequence[Sequence[float]],
    reference_centroid: Sequence[float],
    offset_A: float = SURFACE_PROBE_OFFSET_A,
) -> np.ndarray:
    """Move patch atom positions outward, to where solvent (and RNA) actually is.

    The outward direction is approximated as radially away from the chain's
    centroid. That is exact for a convex protein and adequate for a patch on the
    outside of a globular one, which is the only place a surface-mode patch can
    be. Recorded in the bundle so the approximation is visible.
    """
    xyz = np.asarray(patch_coords, dtype=float)
    centre = np.asarray(reference_centroid, dtype=float)
    directions = xyz - centre
    norms = np.linalg.norm(directions, axis=1, keepdims=True)
    # A point sitting exactly on the centroid has no outward direction; leave it.
    safe = np.where(norms > 1e-6, norms, 1.0)
    return xyz + float(offset_A) * (directions / safe)


def charge_complementarity(
    mean_potential_kT_per_e: float | None, scale_kT_per_e: float = 1.0
) -> dict[str, Any]:
    """Map a patch potential to a complementarity factor in (0, 1).

        f = 1 / (1 + exp(-phi / scale))

    Negative potential (repels the backbone) tends to 0, positive tends to 1,
    and a neutral patch sits at exactly 0.5 rather than at a flattering default.
    A directional signal, never an affinity estimate.
    """
    if mean_potential_kT_per_e is None or not math.isfinite(float(mean_potential_kT_per_e)):
        return {"factor": None, "status": "not_computed",
                "mean_potential_kT_per_e": None, "scale_kT_per_e": float(scale_kT_per_e)}
    phi = float(mean_potential_kT_per_e) / float(scale_kT_per_e)
    # exp overflows for a strongly negative patch; the limit is 0 either way.
    factor = 1.0 / (1.0 + math.exp(-phi)) if -700 < phi < 700 else (1.0 if phi > 0 else 0.0)
    return {
        "factor": factor,
        "status": "computed",
        "mean_potential_kT_per_e": float(mean_potential_kT_per_e),
        "scale_kT_per_e": float(scale_kT_per_e),
        "hospitable_to_rna": bool(mean_potential_kT_per_e > 0),
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
