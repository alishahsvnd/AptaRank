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

BUNDLE_SCHEMA_VERSION = "target-bundle-v2"
#: v1 bundles predate binding modes. They only ever described a cavity, so they
#: are read as pocket-mode targets rather than rejected.
SUPPORTED_SCHEMA_VERSIONS = ("target-bundle-v2", "target-bundle-v1")
LEGACY_SCHEMA_VERSION = "target-bundle-v1"

#: Top-level keys excluded from the bundle id: they change between two builds
#: that are scientifically identical, which is exactly what the id must not do.
#: `target_signature` is a cache tag over the *request*, not evidence about the
#: protein, so it stays outside the hash too.
_VOLATILE_TOP_LEVEL = ("created_utc", "bundle_id", "provenance", "target_signature")

#: Provenance fields that DO belong in the hash. Excluding provenance wholesale
#: would let the fpocket version or the command line be edited without
#: detection, and those determine every number in the bundle. Excluded from
#: here: git state, aptarank version, OS and CI metadata — all of which differ
#: between two honest builds of the same evidence.
_HASHED_PROVENANCE = ("fpocket_version", "fpocket_command", "fpocket_exit_code",
                      "biopython_version")

#: fpocket estimates cavity volume by Monte-Carlo integration (300 iterations by
#: default) and exposes no seed, so the same structure measured twice gives
#: volumes that differ by a few percent. Everything derived from that volume
#: inherits the wobble.
#:
#: These fields are therefore excluded from the bundle id — otherwise no two
#: builds of the same evidence would ever agree, which is exactly what the id
#: exists to certify. They stay in the bundle as measured, and the exclusion is
#: declared in the bundle itself rather than being a silent convention: a reader
#: must be able to see which numbers the id does *not* cover.
#:
#: Nothing scored depends on them. `d_pocket_A`, which drives every pocket-mode
#: band, comes from alpha-sphere coordinates and is exactly reproducible.
NONDETERMINISTIC_FIELDS = (
    "pockets[].fpocket.volume_A3",
    "pockets[].fpocket.metrics.Volume",
    "pockets[].geometry.d_equiv_A",
    "pockets[].geometry.envelope_to_equiv_ratio",
    "pockets[].geometry.shape_warning",
)

#: Recorded command lines name the tool *and* the directory it happened to run
#: in. The tool is evidence; the directory is not — a bundle built in
#: ~/aptarank-data and one built in a CI checkout describe the same protein. For
#: hashing, path arguments are reduced to their filenames, so the command still
#: cannot be edited undetectably while the location no longer matters.
PATH_NORMALISED_FIELDS = (
    "provenance.fpocket.command",
    "electrostatics.pdb2pqr.command",
    "electrostatics.apbs.command",
    # Where the downloaded structure was cached. Its sha256 is hashed in full,
    # so the file itself is still pinned exactly.
    "target.source.path",
)


def build(
    prepared,
    pockets: Sequence[Pocket],
    geometries: Mapping[int, PocketGeometry],
    selection: Mapping[str, Any],
    fpocket_provenance: Mapping[str, Any],
    electrostatics: Mapping[str, Any] | None = None,
    config_hash: str | None = None,
    synthetic: bool = False,
    binding_mode: str = "pocket",
    patch: Mapping[str, Any] | None = None,
    target_signature: str | None = None,
) -> dict[str, Any]:
    """Assemble a bundle from a prepared target and its measurements.

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
        # Which geometry this target was measured for. A bundle prepared for one
        # mode cannot be scored in the other, and says so rather than being
        # compared against whichever number happens to be present.
        "binding_mode": binding_mode,
        "config_hash": config_hash,
        # A tag over the request that produced this bundle, so a repeat run can
        # recognise its own work without re-measuring. Outside the id: it
        # describes what was asked for, not what was found.
        "target_signature": target_signature,
        # Declared, not assumed: these numbers are real measurements but are not
        # covered by bundle_id, because fpocket re-estimates them stochastically.
        "nondeterministic_fields": list(NONDETERMINISTIC_FIELDS),
        "target": {
            "pdb_id": prepared.pdb_id,
            "identifier": getattr(prepared, "identifier", prepared.pdb_id),
            "target_source": getattr(prepared, "source_kind", "pdb"),
            # "experimental" or "predicted" — a measurement taken on a model is
            # not a measurement, and downstream must be able to see which it is.
            "structure_kind": getattr(prepared, "structure_kind", "experimental"),
            "name": prepared.name,
            "model_index": prepared.model_index,
            "chain_id": prepared.chain_id,
            "assembly": "asymmetric_unit",
            "source": prepared.source,
        },
        "preparation": {
            "applied": prepared.applied,
            "site_residues": [
                r.to_dict() for r in getattr(prepared, "site_residues", []) or []
            ],
            "partner_evidence": getattr(prepared, "partner_evidence", {}) or {},
            "cleaned_structure": {
                "filename": Path(prepared.path).name,
                "sha256": _sha256_file(prepared.path),
            },
        },
        "patch": dict(patch) if patch else None,
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
    """SHA-256 over the reproducible payload plus the tool identity behind it.

    Two bundles built from the same structure with the same tool and settings
    get the same id even though their timestamps, git commits and CI run ids
    differ. That is what makes "this run used this target evidence" a checkable
    claim — and why the fpocket version and command line are inside the hash
    while the CI metadata is not.

    fpocket's Monte-Carlo volume estimate (and the handful of numbers derived
    from it) is excluded for the same reason, from the other direction: it
    differs between two honest builds of identical evidence, so including it
    would mean no bundle could ever certify anything.
    """
    payload = {
        k: v for k, v in _without_nondeterministic(bundle).items()
        if k not in _VOLATILE_TOP_LEVEL
    }
    provenance = bundle.get("provenance") or {}
    fpocket = provenance.get("fpocket") or {}
    command = fpocket.get("command")
    payload["_tooling"] = {
        "fpocket_version": fpocket.get("version"),
        # Same reason as PATH_NORMALISED_FIELDS: which directory fpocket ran in
        # is not evidence about the protein.
        "fpocket_command": (
            [_basename_argument(a) for a in command] if isinstance(command, list)
            else command
        ),
        "fpocket_exit_code": fpocket.get("exit_code"),
        "biopython_version": (provenance.get("versions") or {}).get("biopython"),
    }
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha256(blob.encode()).hexdigest()[:32]


def _without_nondeterministic(bundle: Mapping[str, Any]) -> dict[str, Any]:
    """A copy with fpocket's stochastic volume estimates blanked out.

    Blanked rather than deleted, so a bundle that simply omits the field cannot
    hash the same as one that carries it.
    """
    import copy as _copy

    out = _copy.deepcopy(dict(bundle))
    for pocket in out.get("pockets") or []:
        # Only ever overwrite a key that is already there. Adding one would make
        # a bundle with the field stripped out hash the same as one carrying it,
        # which would put a hole in exactly the guarantee this id provides.
        for node, keys in (
            (pocket.get("fpocket"), ("volume_A3",)),
            ((pocket.get("fpocket") or {}).get("metrics"), ("Volume",)),
            (pocket.get("geometry"),
             ("d_equiv_A", "envelope_to_equiv_ratio", "shape_warning")),
        ):
            if not isinstance(node, dict):
                continue
            for key in keys:
                if key in node:
                    node[key] = "<nondeterministic>"

    for dotted in PATH_NORMALISED_FIELDS:
        node: Any = out
        parts = dotted.split(".")
        for part in parts[:-1]:
            node = node.get(part) if isinstance(node, dict) else None
            if node is None:
                break
        if not isinstance(node, dict):
            continue
        value = node.get(parts[-1])
        if isinstance(value, list):
            node[parts[-1]] = [_basename_argument(a) for a in value]
        elif isinstance(value, str):
            node[parts[-1]] = _basename_argument(value)
    return out


def _basename_argument(argument: Any) -> Any:
    """`--apbs-input=/long/path/x.in` -> `--apbs-input=x.in`; leave the rest."""
    if not isinstance(argument, str) or "/" not in argument and "\\" not in argument:
        return argument
    prefix, sep, value = argument.partition("=")
    if sep:
        return f"{prefix}={Path(value).name}"
    return Path(argument).name


def binding_mode(bundle: Mapping[str, Any]) -> str:
    """The mode this bundle was measured for; v1 bundles only ever had cavities."""
    return bundle.get("binding_mode") or "pocket"


def validate(bundle: Mapping[str, Any]) -> None:
    """Reject a bundle that cannot support a Tier 2 score in its own mode."""
    version = bundle.get("schema_version")
    if version not in SUPPORTED_SCHEMA_VERSIONS:
        raise TargetError(
            f"unsupported target bundle schema {version!r}; expected one of "
            f"{SUPPORTED_SCHEMA_VERSIONS}"
        )
    for key in ("target", "preparation", "pockets", "selection", "provenance"):
        if key not in bundle:
            raise TargetError(f"target bundle is missing required section {key!r}")

    mode = binding_mode(bundle)
    if mode == "surface":
        patch = bundle.get("patch") or {}
        area = patch.get("patch_area_A2")
        if not isinstance(area, (int, float)) or not area > 0:
            raise TargetError(
                f"surface-mode bundle has no measured patch area (got {area!r}); "
                f"it cannot support a surface-mode score"
            )
        if not patch.get("residue_numbers"):
            raise TargetError("surface-mode bundle records no binding-site residues")
        return

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
    target = bundle["target"]
    applied = bundle["preparation"].get("applied", {})
    electro = bundle.get("electrostatics", {})
    sampling = electro.get("sampling") or electro.get("selected_pocket_sampling") or {}
    patch = bundle.get("patch") or {}
    mode = binding_mode(bundle)

    out: dict[str, Any] = {
        "pdb_id": target["pdb_id"],
        "identifier": target.get("identifier", target["pdb_id"]),
        "target_source": target.get("target_source", "pdb"),
        "structure_kind": target.get("structure_kind", "experimental"),
        "chain": target["chain_id"],
        "name": target["name"],
        "bundle_id": bundle["bundle_id"],
        "binding_mode": mode,
        "synthetic": bool(bundle.get("synthetic", False)),
        "n_pockets": len(bundle["pockets"]),
        "pocket_selection": bundle["selection"].get("method"),
        "selection_warnings": list(bundle["selection"].get("warnings", [])),
        "preparation_warnings": list(applied.get("warnings", [])),
        "was_multi_chain": bool(applied.get("was_multi_chain")),
        "chains_removed": applied.get("chains_removed", []),
        "partner_chains_removed": applied.get("partner_chains_removed", []),
        "site_residues": [
            r["residue_number"] for r in bundle["preparation"].get("site_residues", [])
        ],
        "electrostatics_status": electro.get("status"),
        "electrostatic_mean_potential": sampling.get("mean_potential_kT_per_e"),
        "electrostatic_compatible": sampling.get("electrostatic_compatible"),
        "retained_hetero": [
            r["residue_name"] for r in applied.get("retained_hetero_residues", [])
        ],
        "patch": None,
        "selected_pocket": None,
    }

    if patch:
        out["patch"] = {
            "patch_area_A2": patch.get("patch_area_A2"),
            "planarity_A": patch.get("planarity_A"),
            "elongation": patch.get("elongation"),
            "n_residues": patch.get("n_residues"),
            "definition": patch.get("definition"),
            "buried_residue_numbers": patch.get("buried_residue_numbers", []),
            "shape_warning": bool(patch.get("shape_warning")),
        }

    if bundle["pockets"] and bundle["selection"].get("selected_pocket_index") is not None:
        pocket = selected_pocket(bundle)
        geometry = pocket["geometry"]
        out["selected_pocket"] = {
            "index": pocket["index"],
            "volume_A3": pocket["fpocket"]["volume_A3"],
            "fpocket_score": pocket["fpocket"]["score"],
            "d_pocket_A": geometry["d_pocket_A"],
            "d_pocket_centres_A": geometry["d_pocket_centres_A"],
            "d_equiv_A": geometry["d_equiv_A"],
            "n_alpha_spheres": geometry["n_alpha_spheres"],
            "shape_warning": geometry["shape_warning"],
            "elongation": geometry.get("elongation"),
            "planarity_A": geometry.get("planarity_A"),
        }
    return out


def write(bundle: Mapping[str, Any], directory: str | Path) -> Path:
    out_dir = Path(directory)
    out_dir.mkdir(parents=True, exist_ok=True)
    target = bundle["target"]
    # The mode is in the filename so the same protein can hold a pocket bundle
    # and a surface bundle side by side, and so a lookup never has to open every
    # bundle on disk to find out which is which.
    name = (
        f"{target['pdb_id']}_{target['chain_id']}_{binding_mode(bundle)}_"
        f"{bundle['bundle_id'][:12]}.bundle.json"
    )
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


def find(
    directory: str | Path,
    identifier: str,
    chain_id: str | None = None,
    mode: str | None = None,
) -> Path:
    """Most recent prepared target for an identifier in `directory`.

    A mode-specific match wins; the mode-less pattern is the fallback that keeps
    v1 bundles (which had no mode in their name) findable.
    """
    d = Path(directory)
    stem = f"{identifier.upper()}_{chain_id}" if chain_id else identifier.upper()
    patterns = [f"{stem}_*.bundle.json"]
    if mode:
        patterns.insert(0, f"{stem}_{mode}_*.bundle.json")

    for pattern in patterns:
        matches = sorted(d.glob(pattern), key=lambda p: p.stat().st_mtime, reverse=True)
        if matches:
            return matches[0]
    raise TargetError(
        f"no prepared target for {identifier}"
        f"{f' chain {chain_id}' if chain_id else ''}"
        f"{f' in {mode} mode' if mode else ''} in {d}. "
        f"Prepare one with `aptarank target build --id {identifier}`."
    )


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
