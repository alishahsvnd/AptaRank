"""The target bundle: immutable, checksummed target evidence.

Everything that needs an external Linux tool happens once, on Linux, and lands
in a bundle. Scoring then consumes the bundle and needs no external tools at
all — which is what lets the pipeline run identically on a laptop, in CI, and
during the live demo, and lets the paper cite a repository that regenerates the
evidence rather than a machine that happened to have fpocket installed.

The bundle carries *all* pockets and the full selection evidence, not just the
chosen cavity, so an automatic pocket choice can be audited after the fact.
Alpha-sphere coordinates are stored inline: they are small, and every geometric
quantity must be recomputable without re-running fpocket.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence

from .. import __version__
from ..config import REPO_ROOT
from ..errors import TargetError
from ..provenance import git_state, tool_versions
from .fpocket import Pocket
from .geometry import PocketGeometry

BUNDLE_SCHEMA_VERSION = "target-bundle-v1"

#: Top-level keys excluded from the bundle id: they change between two builds
#: that are scientifically identical, which is exactly what the id must not do.
_VOLATILE_TOP_LEVEL = ("created_utc", "bundle_id", "provenance")

#: Provenance fields that DO belong in the hash. Excluding provenance wholesale
#: would let the fpocket version or the command line be edited without
#: detection, and those determine every number in the bundle. Excluded from
#: here: git state, aptarank version, OS and CI metadata — all of which differ
#: between two honest builds of the same evidence.
_HASHED_PROVENANCE = ("fpocket_version", "fpocket_command", "fpocket_exit_code",
                      "biopython_version")


def build(
    prepared,
    pockets: Sequence[Pocket],
    geometries: Mapping[int, PocketGeometry],
    selection: Mapping[str, Any],
    fpocket_provenance: Mapping[str, Any],
    electrostatics: Mapping[str, Any] | None = None,
    config_hash: str | None = None,
    synthetic: bool = False,
) -> dict[str, Any]:
    """Assemble a bundle from prepared target + parsed fpocket output.

    `synthetic=True` marks a bundle whose cavity was fabricated rather than
    detected by fpocket. Such bundles exist so the demo and the tests can run
    on machines without fpocket; they carry the flag everywhere it can be seen
    and make any run that uses them ineligible for publication.
    """
    bundle: dict[str, Any] = {
        "schema_version": BUNDLE_SCHEMA_VERSION,
        "bundle_id": None,
        "created_utc": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        # Part of the scientific payload, not provenance: a fabricated cavity
        # changes what every downstream number means, so it changes the id.
        "synthetic": bool(synthetic),
        "config_hash": config_hash,
        "target": {
            "pdb_id": prepared.pdb_id,
            "name": prepared.name,
            "model_index": prepared.model_index,
            "chain_id": prepared.chain_id,
            "assembly": "asymmetric_unit",
            "source": prepared.source,
        },
        "preparation": {
            "applied": prepared.applied,
            "cleaned_structure": {
                "filename": Path(prepared.path).name,
                "sha256": _sha256_file(prepared.path),
            },
        },
        "pockets": [
            {
                "index": pocket.index,
                "fpocket": {
                    "score": pocket.score,
                    "druggability_score": pocket.druggability,
                    "volume_A3": pocket.volume_A3,
                    "n_alpha_spheres_reported": pocket.n_alpha_spheres_reported,
                    "metrics": pocket.metrics,
                },
                "lining_residues": [r.to_dict() for r in pocket.lining_residues],
                "alpha_spheres": [
                    {
                        "index": s.index,
                        "center_A": list(s.center_A),
                        "radius_A": s.radius_A,
                        "kind": s.kind,
                    }
                    for s in pocket.alpha_spheres
                ],
                "geometry": geometries[pocket.index].to_dict(),
            }
            for pocket in pockets
        ],
        "selection": dict(selection),
        "electrostatics": dict(electrostatics) if electrostatics else _skipped_electrostatics(),
        "provenance": {
            "aptarank_version": __version__,
            "git": git_state(REPO_ROOT),
            "versions": tool_versions(),
            "fpocket": {
                k: v for k, v in fpocket_provenance.items()
                if k not in ("stdout", "stderr", "output_dir")
            },
        },
    }
    bundle["bundle_id"] = compute_bundle_id(bundle)
    validate(bundle)
    return bundle


def compute_bundle_id(bundle: Mapping[str, Any]) -> str:
    """SHA-256 over the scientific payload plus the tool identity that produced it.

    Two bundles built from the same structure with the same tool and settings
    get the same id even though their timestamps, git commits and CI run ids
    differ. That is what makes "this run used this target evidence" a checkable
    claim — and why the fpocket version and command line are inside the hash
    while the CI metadata is not.
    """
    payload = {k: v for k, v in bundle.items() if k not in _VOLATILE_TOP_LEVEL}
    provenance = bundle.get("provenance") or {}
    fpocket = provenance.get("fpocket") or {}
    payload["_tooling"] = {
        "fpocket_version": fpocket.get("version"),
        "fpocket_command": fpocket.get("command"),
        "fpocket_exit_code": fpocket.get("exit_code"),
        "biopython_version": (provenance.get("versions") or {}).get("biopython"),
    }
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha256(blob.encode()).hexdigest()[:32]


def validate(bundle: Mapping[str, Any]) -> None:
    """Reject a bundle that cannot support a Tier 2 score."""
    if bundle.get("schema_version") != BUNDLE_SCHEMA_VERSION:
        raise TargetError(
            f"unsupported target bundle schema {bundle.get('schema_version')!r}; "
            f"expected {BUNDLE_SCHEMA_VERSION}"
        )
    for key in ("target", "preparation", "pockets", "selection", "provenance"):
        if key not in bundle:
            raise TargetError(f"target bundle is missing required section {key!r}")
    if not bundle["pockets"]:
        raise TargetError("target bundle contains no pockets")

    selection = bundle["selection"]
    if selection.get("status") != "selected":
        raise TargetError(f"target bundle selection status is {selection.get('status')!r}")
    index = selection.get("selected_pocket_index")
    if index is None or not any(p["index"] == index for p in bundle["pockets"]):
        raise TargetError(f"selected pocket {index!r} is not present in the bundle")

    pocket = selected_pocket(bundle)
    d_pocket = pocket["geometry"]["d_pocket_A"]
    if not isinstance(d_pocket, (int, float)) or not d_pocket > 0:
        raise TargetError(f"selected pocket has a non-positive d_pocket_A: {d_pocket!r}")


def selected_pocket(bundle: Mapping[str, Any]) -> dict[str, Any]:
    index = bundle["selection"]["selected_pocket_index"]
    for pocket in bundle["pockets"]:
        if pocket["index"] == index:
            return pocket
    raise TargetError(f"selected pocket {index} not found in bundle")


def summary(bundle: Mapping[str, Any]) -> dict[str, Any]:
    """The target-level block embedded in a run artifact (spec §7.2)."""
    pocket = selected_pocket(bundle)
    geometry = pocket["geometry"]
    electro = bundle.get("electrostatics", {})
    return {
        "pdb_id": bundle["target"]["pdb_id"],
        "chain": bundle["target"]["chain_id"],
        "name": bundle["target"]["name"],
        "bundle_id": bundle["bundle_id"],
        "synthetic": bool(bundle.get("synthetic", False)),
        "n_pockets": len(bundle["pockets"]),
        "pocket_selection": bundle["selection"]["method"],
        "selection_warnings": bundle["selection"].get("warnings", []),
        "selected_pocket": {
            "index": pocket["index"],
            "volume_A3": pocket["fpocket"]["volume_A3"],
            "fpocket_score": pocket["fpocket"]["score"],
            "d_pocket_A": geometry["d_pocket_A"],
            "d_pocket_centres_A": geometry["d_pocket_centres_A"],
            "d_equiv_A": geometry["d_equiv_A"],
            "n_alpha_spheres": geometry["n_alpha_spheres"],
            "shape_warning": geometry["shape_warning"],
        },
        "electrostatics_status": electro.get("status"),
        "electrostatic_mean_potential": (
            (electro.get("selected_pocket_sampling") or {}).get("mean_potential_kT_per_e")
        ),
        "electrostatic_compatible": (
            (electro.get("selected_pocket_sampling") or {}).get("electrostatic_compatible")
        ),
        "retained_hetero": [
            r["residue_name"]
            for r in bundle["preparation"]["applied"].get("retained_hetero_residues", [])
        ],
    }


def write(bundle: Mapping[str, Any], directory: str | Path) -> Path:
    out_dir = Path(directory)
    out_dir.mkdir(parents=True, exist_ok=True)
    target = bundle["target"]
    name = f"{target['pdb_id']}_{target['chain_id']}_{bundle['bundle_id'][:12]}.bundle.json"
    path = out_dir / name
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(bundle, indent=2, allow_nan=False, default=str), encoding="utf-8")
    tmp.replace(path)
    return path


def load(path: str | Path) -> dict[str, Any]:
    p = Path(path)
    if not p.exists():
        raise TargetError(f"target bundle not found: {p}")
    bundle = json.loads(p.read_text(encoding="utf-8"))
    validate(bundle)

    recomputed = compute_bundle_id(bundle)
    if bundle["bundle_id"] != recomputed:
        raise TargetError(
            f"{p}: bundle_id does not match its contents "
            f"(stored {bundle['bundle_id']}, recomputed {recomputed}). "
            f"The bundle has been edited or truncated."
        )
    return bundle


def find(directory: str | Path, pdb_id: str, chain_id: str | None = None) -> Path:
    """Most recent bundle for a target in `directory`."""
    d = Path(directory)
    pattern = f"{pdb_id.upper()}_{chain_id}_*.bundle.json" if chain_id else f"{pdb_id.upper()}_*.bundle.json"
    matches = sorted(d.glob(pattern), key=lambda p: p.stat().st_mtime, reverse=True)
    if not matches:
        raise TargetError(
            f"no target bundle for {pdb_id}"
            f"{f' chain {chain_id}' if chain_id else ''} in {d}. "
            f"Build one on Linux with `aptarank target build --pdb-id {pdb_id}`, "
            f"or fetch it from the repository's target-bundles CI workflow."
        )
    return matches[0]


def _skipped_electrostatics() -> dict[str, Any]:
    return {
        "requested": False,
        "status": "skipped",
        "reason_code": "not_requested",
        "message": None,
        "pdb2pqr": None,
        "apbs": None,
        "grid": None,
        "selected_pocket_sampling": None,
    }


def _sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as fh:
        while block := fh.read(1 << 20):
            digest.update(block)
    return digest.hexdigest()
