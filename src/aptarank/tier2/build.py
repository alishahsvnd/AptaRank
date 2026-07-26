"""Build a target bundle. Runs the Linux-only tools; produces portable JSON.

This is the only module that requires fpocket. It is intended to run on Linux
— locally, in WSL, or in the pinned CI job — and its whole purpose is to make
everything downstream platform-independent.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from ..config import Config
from ..errors import TargetError
from . import bundle as bundle_mod
from . import fpocket, selection, target
from .geometry import pocket_geometry


def build_target_bundle(
    cfg: Config,
    pdb_id: str | None = None,
    chain_id: str | None = None,
    out_dir: str | Path | None = None,
    work_dir: str | Path | None = None,
) -> tuple[dict[str, Any], Path]:
    """Prepare the structure, detect cavities, select one, write the bundle."""
    pdb_id = pdb_id or cfg.get("tier2.target.pdb_id", None)
    if not pdb_id:
        raise TargetError("no target PDB ID configured (tier2.target.pdb_id)")
    chain_id = chain_id or cfg.get("tier2.target.chain", None)

    work = Path(work_dir or Path(cfg.get("tier2.bundle_dir")) / "work")
    prepared = target.prepare(
        pdb_id=pdb_id,
        cache_dir=cfg.get("tier2.structure_cache_dir"),
        work_dir=work,
        chain_id=chain_id,
        model_index=int(cfg.get("tier2.target.model", 0)),
        retain_hetero_resnames=cfg.get("tier2.target.retain_hetero_resnames", []),
        remove_hetero_resnames=cfg.get("tier2.target.remove_hetero_resnames", []),
        hetero_default=cfg.get("tier2.target.hetero_default", "retain"),
    )

    provenance = fpocket.run(prepared.path)
    pockets = fpocket.load_pockets(provenance["output_dir"], prepared.path.stem)

    geometries = {}
    for pocket in pockets:
        geometries[pocket.index] = pocket_geometry(
            centres=[s.center_A for s in pocket.alpha_spheres],
            radii=[s.radius_A for s in pocket.alpha_spheres],
            volume_A3=pocket.volume_A3,
        )

    requested = selection.parse_residue_specs(
        cfg.get("tier2.target.active_site_residues", []), prepared.chain_id
    )
    chosen = selection.select_pocket(
        pockets,
        requested=requested,
        structure_residues=prepared.residues,
        allow_zero_overlap_fallback=bool(
            cfg.get("tier2.target.allow_zero_overlap_fallback", False)
        ),
    )

    electrostatics = None
    if cfg.get("tier2.electrostatics.enabled", False):
        from . import apbs

        electrostatics = apbs.compute(
            prepared.path,
            pockets[[p.index for p in pockets].index(chosen["selected_pocket_index"])],
            work_dir=work,
        )

    result = bundle_mod.build(
        prepared=prepared,
        pockets=pockets,
        geometries=geometries,
        selection=chosen,
        fpocket_provenance=provenance,
        electrostatics=electrostatics,
        config_hash=cfg.scoring_signature(),
    )
    path = bundle_mod.write(result, out_dir or cfg.get("tier2.bundle_dir"))
    return result, path
