"""Prepare a target: fetch it, measure it, write an immutable bundle.

This is the only module that needs fpocket (and, for surface mode, APBS). It
runs **server-side** now rather than on the user's machine: the refinements spec
moves the biology-literate preparation step off the biologist and into the
pipeline, so the dashboard asks for an identifier and a chain rather than for a
prepared JSON file.

Everything downstream still consumes the bundle and needs no external tools, so
scoring runs identically on a laptop, in CI, and during a live demo.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from ..config import Config
from ..errors import ExternalToolError, TargetError
from . import bundle as bundle_mod
from . import fpocket, modes, patch as patch_mod, selection, target
from .geometry import pocket_geometry


def target_signature(cfg: Config) -> str:
    """Everything about a target that changes what was measured.

    Used as the cache key, so re-running the same target is free while changing
    the chain, the site residues or the mode rebuilds rather than silently
    reusing the previous measurement.
    """
    import hashlib
    import json

    relevant = {
        "schema": bundle_mod.BUNDLE_SCHEMA_VERSION,
        "source": cfg.get("tier2.target.source"),
        "id": (cfg.get("tier2.target.id") or "").upper(),
        "chain": cfg.get("tier2.target.chain"),
        "partner_chains": cfg.get("tier2.target.partner_chains"),
        "site_residues": sorted(cfg.get("tier2.target.target_site_residues") or []),
        "binding_mode": cfg.get("tier2.binding_mode"),
        "model": cfg.get("tier2.target.model"),
        "strip_hetatm": cfg.get("tier2.target.strip_hetatm"),
        "hetero_default": cfg.get("tier2.target.hetero_default"),
        "retain_hetero": sorted(cfg.get("tier2.target.retain_hetero_resnames") or []),
        "remove_hetero": sorted(cfg.get("tier2.target.remove_hetero_resnames") or []),
        "patch_area_definition": cfg.get("tier2.surface.patch_area_definition"),
        "electrostatics": bool(_electrostatics_wanted(cfg)),
    }
    blob = json.dumps(relevant, sort_keys=True, default=str)
    return hashlib.sha256(blob.encode()).hexdigest()[:12]


def _electrostatics_wanted(cfg: Config) -> bool:
    """Surface mode asks for the patch potential; pocket mode keeps it opt-in.

    In surface mode the charge term is part of the documented comparison, so it
    is computed by default. It still never blocks a build.
    """
    if cfg.get("tier2.binding_mode") == modes.SURFACE:
        return bool(cfg.get("tier2.surface.charge.enabled", True))
    return bool(cfg.get("tier2.electrostatics.enabled", False))


def build_target_bundle(
    cfg: Config,
    identifier: str | None = None,
    chain_id: str | None = None,
    out_dir: str | Path | None = None,
    work_dir: str | Path | None = None,
    progress: Any = None,
) -> tuple[dict[str, Any], Path]:
    """Prepare the structure, measure it in the configured mode, write a bundle."""
    identifier = identifier or cfg.get("tier2.target.id", None)
    if not identifier:
        raise TargetError("no target configured (tier2.target.id)")
    chain_id = chain_id or cfg.get("tier2.target.chain", None)
    mode = modes.check_mode(cfg.get("tier2.binding_mode"))
    source = cfg.get("tier2.target.source")
    site_residues = list(cfg.get("tier2.target.target_site_residues") or [])

    if mode == modes.SURFACE and not site_residues:
        raise TargetError(
            "surface mode needs binding-site residues: the patch it measures is "
            "defined by them. Supply tier2.target.target_site_residues, or use "
            "pocket mode, where the cavity can be chosen automatically."
        )

    def report(message: str) -> None:
        if progress:
            progress(message)

    work = Path(work_dir or Path(cfg.get("tier2.bundle_dir")) / "work")
    report(f"Fetching {identifier} from {source}")
    prepared = target.prepare(
        identifier=identifier,
        source=source,
        cache_dir=cfg.get("tier2.structure_cache_dir"),
        work_dir=work,
        chain_id=chain_id,
        partner_chains=cfg.get("tier2.target.partner_chains", []),
        target_site_residues=site_residues,
        model_index=int(cfg.get("tier2.target.model", 0)),
        strip_hetatm=bool(cfg.get("tier2.target.strip_hetatm", False)),
        retain_hetero_resnames=cfg.get("tier2.target.retain_hetero_resnames", []),
        remove_hetero_resnames=cfg.get("tier2.target.remove_hetero_resnames", []),
        hetero_default=cfg.get("tier2.target.hetero_default", "retain"),
    )

    report("Detecting cavities with fpocket")
    pockets, provenance = _detect_pockets(prepared, required=mode == modes.POCKET)

    geometries = {
        pocket.index: pocket_geometry(
            centres=[s.center_A for s in pocket.alpha_spheres],
            radii=[s.radius_A for s in pocket.alpha_spheres],
            volume_A3=pocket.volume_A3,
        )
        for pocket in pockets
    }

    requested = selection.parse_residue_specs(site_residues, prepared.chain_id)
    chosen = selection.select_pocket(
        pockets,
        requested=requested,
        structure_residues=prepared.residues,
        allow_zero_overlap_fallback=bool(
            cfg.get("tier2.target.allow_zero_overlap_fallback", False)
        ),
        # In surface mode the patch is the evidence and any cavity is only a
        # cross-reference, so a binding site that lines no cavity is expected
        # rather than a build failure.
        require_overlap=mode == modes.POCKET,
    )

    measured_patch = None
    if mode == modes.SURFACE:
        report("Measuring the binding-site patch")
        geometry = patch_mod.patch_geometry(
            prepared.path,
            site_residues,
            definition=cfg.get("tier2.surface.patch_area_definition"),
        )
        measured_patch = geometry.to_dict()
        if pockets and chosen.get("selected_pocket_index") is not None:
            index = chosen["selected_pocket_index"]
            lining = next(p.lining_residues for p in pockets if p.index == index)
            measured_patch["pocket_overlap_residue_numbers"] = (
                patch_mod.overlapping_residue_numbers(lining, site_residues)
            )

    electrostatics = None
    if _electrostatics_wanted(cfg):
        report("Computing electrostatics (PDB2PQR + APBS)")
        electrostatics = _electrostatics(cfg, prepared, pockets, chosen, mode, work)

    result = bundle_mod.build(
        prepared=prepared,
        pockets=pockets,
        geometries=geometries,
        selection=chosen,
        fpocket_provenance=provenance,
        electrostatics=electrostatics,
        config_hash=cfg.scoring_signature(),
        binding_mode=mode,
        patch=measured_patch,
        target_signature=target_signature(cfg),
    )
    path = bundle_mod.write(result, out_dir or cfg.get("tier2.bundle_dir"))
    return result, path


def _detect_pockets(prepared, required: bool) -> tuple[list, dict[str, Any]]:
    """Run fpocket. In surface mode its absence is a note, not a failure."""
    try:
        provenance = fpocket.run(prepared.path)
        pockets = fpocket.load_pockets(provenance["output_dir"], prepared.path.stem)
        return pockets, provenance
    except ExternalToolError:
        if required:
            raise
        return [], {
            "status": "skipped",
            "reason": "fpocket unavailable or found no cavity; surface mode does "
                      "not require one",
            "version": None,
            "command": None,
            "exit_code": None,
        }


def _electrostatics(cfg, prepared, pockets, chosen, mode, work) -> dict[str, Any]:
    """Sample the potential where this mode needs it: patch or cavity."""
    import numpy as np
    from Bio.PDB import PDBParser

    from . import apbs

    if mode == modes.SURFACE:
        wanted = {int(n) for n in cfg.get("tier2.target.target_site_residues") or []}
        structure = PDBParser(QUIET=True).get_structure("t", str(prepared.path))
        atoms = np.array([a.coord for a in structure.get_atoms()], dtype=float)
        patch_atoms = np.array(
            [a.coord for r in structure.get_residues() if int(r.id[1]) in wanted for a in r],
            dtype=float,
        )
        if patch_atoms.size == 0:
            return {
                "requested": True, "status": "failed",
                "reason_code": "no_patch_atoms", "message":
                "no binding-site atoms were found in the prepared structure",
                "pdb2pqr": None, "apbs": None, "grid": None, "sampling": None,
                "selected_pocket_sampling": None,
            }
        points = apbs.surface_sample_points(patch_atoms, atoms.mean(axis=0))
        result = apbs.compute(prepared.path, points, work_dir=work,
                              label="binding_site_patch")
        result["sampling_geometry"] = {
            "method": "patch atoms displaced outward from the chain centroid",
            "offset_A": apbs.SURFACE_PROBE_OFFSET_A,
            "n_patch_atoms": int(patch_atoms.shape[0]),
        }
        return result

    index = chosen.get("selected_pocket_index")
    pocket = next((p for p in pockets if p.index == index), None)
    if pocket is None:
        return {
            "requested": True, "status": "skipped",
            "reason_code": "no_selected_pocket", "message":
            "no cavity was selected, so there was nothing to sample",
            "pdb2pqr": None, "apbs": None, "grid": None, "sampling": None,
            "selected_pocket_sampling": None,
        }
    return apbs.compute(prepared.path, pocket, work_dir=work, label="selected_pocket")
