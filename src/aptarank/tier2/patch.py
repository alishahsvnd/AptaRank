"""Surface-patch measurement (refinements §4.1, surface mode).

A binding interface is not a cavity, so fpocket's alpha spheres are the wrong
instrument for it. What surface mode needs is the *exposed* surface of the
residues the biologist named, measured on the isolated target chain — that is,
after the binding partner has been stripped, which is exactly the surface an
aptamer would have to cover.

Three numbers, all taken at the configured binding-site residues:

    patch_area_A2   solvent-accessible area of the residue set (freeSASA)
    planarity_A     extent along the thinnest principal axis; small = flat
    elongation      longest / shortest principal spread; large = a groove

The area definition is deliberately explicit and recorded in the bundle. The
whole selected-residue set is the default; the alternative (only the subset
lining a detected cavity) is computed alongside it so the choice can be revisited
against real numbers rather than argued in the abstract.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np

from ..errors import ExternalToolError, TargetError
from .geometry import LOWER_Q, UPPER_Q, elongation, planarity

PATCH_ALGORITHM = "residue-patch-sasa-v1"
AREA_DEFINITIONS = ("selected_residues", "pocket_overlap")


@dataclass
class PatchGeometry:
    """Shape and area of one surface patch, measured at named residues."""

    algorithm: str
    definition: str
    n_residues: int
    n_atoms: int
    residue_numbers: list[int]
    area_A2: float
    per_residue_area_A2: dict[int, float]
    centroid_A: list[float]
    extents_A: list[float]
    planarity_A: float
    elongation: float
    buried_residue_numbers: list[int]
    shape_warning: bool
    total_chain_area_A2: float | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "algorithm": self.algorithm,
            "definition": self.definition,
            "n_residues": self.n_residues,
            "n_atoms": self.n_atoms,
            "residue_numbers": self.residue_numbers,
            "patch_area_A2": self.area_A2,
            "per_residue_area_A2": {str(k): v for k, v in self.per_residue_area_A2.items()},
            "centroid_A": self.centroid_A,
            "extents_A": self.extents_A,
            "planarity_A": self.planarity_A,
            "elongation": self.elongation,
            "buried_residue_numbers": self.buried_residue_numbers,
            "shape_warning": self.shape_warning,
            "total_chain_area_A2": self.total_chain_area_A2,
            "quantiles": {"lower": LOWER_Q, "upper": UPPER_Q, "method": "linear"},
        }


def residue_sasa(structure_path: str | Path) -> tuple[dict[int, float], float]:
    """Per-residue solvent-accessible area of a single-chain PDB, in A^2.

    Keyed by residue number, which is unambiguous here because the structure has
    already been reduced to one chain with no insertion codes (§3.6).
    """
    try:
        import freesasa
    except ImportError as exc:  # pragma: no cover - environment problem
        raise ExternalToolError(
            "freesasa is required for surface-mode targets. Install it with "
            "`pip install freesasa`."
        ) from exc

    path = Path(structure_path)
    if not path.exists():
        raise TargetError(f"prepared structure not found: {path}")

    # freesasa writes classifier complaints to stderr for anything non-standard;
    # they are not errors and would otherwise land in a user's job log.
    freesasa.setVerbosity(freesasa.silent)
    try:
        structure = freesasa.Structure(str(path))
        result = freesasa.calc(structure)
    except Exception as exc:  # noqa: BLE001 - freesasa raises bare exceptions
        raise ExternalToolError(f"freesasa failed on {path.name}: {exc}") from exc

    per_residue: dict[int, float] = {}
    for i in range(structure.nAtoms()):
        try:
            number = int(str(structure.residueNumber(i)).strip())
        except ValueError:
            continue
        per_residue[number] = per_residue.get(number, 0.0) + float(result.atomArea(i))
    return per_residue, float(result.totalArea())


def patch_geometry(
    structure_path: str | Path,
    residue_numbers: Sequence[int],
    definition: str = "selected_residues",
    buried_threshold_A2: float = 1.0,
    flat_threshold_A: float = 12.0,
) -> PatchGeometry:
    """Measure the named residues as one surface patch."""
    from Bio.PDB import PDBParser

    if definition not in AREA_DEFINITIONS:
        raise TargetError(
            f"unknown patch area definition {definition!r}; expected one of "
            f"{AREA_DEFINITIONS}"
        )
    wanted = [int(n) for n in residue_numbers]
    if not wanted:
        raise TargetError(
            "surface mode needs binding-site residues: there is no patch to "
            "measure without them. Set tier2.target.target_site_residues."
        )

    areas, total = residue_sasa(structure_path)
    structure = PDBParser(QUIET=True).get_structure("patch", str(structure_path))
    coords, seen = [], []
    for residue in structure.get_residues():
        number = int(residue.id[1])
        if number not in wanted or residue.id[0] != " ":
            continue
        seen.append(number)
        coords.extend(atom.coord for atom in residue)

    missing = sorted(set(wanted) - set(seen))
    if missing:
        raise TargetError(
            f"binding-site residue(s) {missing[:8]} are absent from the prepared "
            f"structure; the cleaning step may have removed them"
        )

    xyz = np.asarray(coords, dtype=float)
    centred = xyz - xyz.mean(axis=0)
    _u, sv, vt = np.linalg.svd(centred, full_matrices=False)
    projections = centred @ vt.T
    extents = [
        float(np.quantile(projections[:, j], UPPER_Q) - np.quantile(projections[:, j], LOWER_Q))
        for j in range(3)
    ]

    per_residue = {n: float(areas.get(n, 0.0)) for n in sorted(set(wanted))}
    # A configured residue with no exposed area is buried inside the fold. It
    # cannot be part of a binding surface, and silently summing a zero would
    # hide a residue-numbering mistake behind a plausible total.
    buried = [n for n, area in per_residue.items() if area < buried_threshold_A2]
    patch_planarity = planarity(sv, xyz.shape[0])

    return PatchGeometry(
        algorithm=PATCH_ALGORITHM,
        definition=definition,
        n_residues=len(per_residue),
        n_atoms=int(xyz.shape[0]),
        residue_numbers=sorted(per_residue),
        area_A2=float(sum(per_residue.values())),
        per_residue_area_A2=per_residue,
        centroid_A=[float(v) for v in xyz.mean(axis=0)],
        extents_A=extents,
        planarity_A=patch_planarity,
        elongation=elongation(sv),
        buried_residue_numbers=buried,
        # Surface mode assumes something an aptamer can lie against. A thick,
        # highly curved "patch" is flagged rather than quietly scored.
        shape_warning=bool(
            np.isfinite(patch_planarity) and patch_planarity > flat_threshold_A
        ),
        total_chain_area_A2=total,
    )


def overlapping_residue_numbers(
    pocket_lining: Iterable[Any], residue_numbers: Sequence[int]
) -> list[int]:
    """Configured residues that also line a detected cavity (the alternative
    area definition, computed for comparison)."""
    wanted = {int(n) for n in residue_numbers}
    return sorted({r.residue_number for r in pocket_lining if r.residue_number in wanted})
