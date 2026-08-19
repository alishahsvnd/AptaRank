"""Verify a target bundle beyond schema validity.

Run in CI immediately after a bundle is built. The point is that a permissive
parser or a silent format drift must not be able to produce a bundle that looks
fine: every geometric quantity is recomputed from the inline alpha spheres and
compared against what was stored, and the build is repeated to confirm the
scientific payload is stable.

    python scripts/verify_target_bundle.py bundles/*.bundle.json [--rebuild]
"""

from __future__ import annotations

import argparse
import glob
import sys
from pathlib import Path

from aptarank.tier2 import bundle as bundle_mod
from aptarank.tier2.geometry import pocket_geometry

TOLERANCE = 1e-6


def verify(path: Path, rebuild: bool = False) -> list[str]:
    problems: list[str] = []
    bundle = bundle_mod.load(path)          # also re-checks the bundle id
    mode = bundle_mod.binding_mode(bundle)
    target = bundle["target"]
    print(f"  bundle_id       {bundle['bundle_id']}")
    print(f"  target          {target['pdb_id']} chain {target['chain_id']}")
    print(f"  binding mode    {mode}")
    print(f"  structure       {target.get('structure_kind', 'experimental')} "
          f"(from {target.get('target_source', 'pdb')})")
    print(f"  pockets         {len(bundle['pockets'])}")
    print(f"  selection       {bundle['selection']['method']} "
          f"-> pocket {bundle['selection']['selected_pocket_index']}")

    problems.extend(_check_patch(bundle))

    for pocket in bundle["pockets"]:
        index = pocket["index"]
        spheres = pocket["alpha_spheres"]
        stored = pocket["geometry"]

        if len(spheres) != pocket["fpocket"]["n_alpha_spheres_reported"]:
            problems.append(
                f"pocket {index}: {len(spheres)} inline spheres but fpocket reported "
                f"{pocket['fpocket']['n_alpha_spheres_reported']}"
            )
            continue

        recomputed = pocket_geometry(
            [s["center_A"] for s in spheres],
            [s["radius_A"] for s in spheres],
            pocket["fpocket"]["volume_A3"],
        ).to_dict()

        for key in ("d_pocket_A", "d_pocket_centres_A", "d_equiv_A"):
            a, b = stored[key], recomputed[key]
            if a != b and abs(a - b) > TOLERANCE * max(1.0, abs(a)):
                problems.append(
                    f"pocket {index}: stored {key}={a} but recomputed {b}"
                )

        if not pocket["lining_residues"]:
            problems.append(f"pocket {index}: no lining residues recorded")
        if not (pocket["fpocket"]["volume_A3"] > 0):
            problems.append(f"pocket {index}: non-positive volume")

    if bundle["selection"].get("selected_pocket_index") is not None:
        selected = bundle_mod.selected_pocket(bundle)
        if selected["geometry"]["shape_warning"]:
            print("  !! selected cavity is oddly shaped (envelope/equivalent-sphere "
                  "ratio out of range); the geometric comparison assumes a roughly "
                  "convex pocket")
        if mode == "pocket" and bundle["selection"]["method"] != "target_site_overlap":
            print(f"  !! cavity was NOT selected by binding-site overlap "
                  f"({bundle['selection']['method']}); the dashboard must caveat this")
    for warning in bundle["selection"].get("warnings", []):
        print(f"  !! {warning}")
    for warning in bundle["preparation"].get("applied", {}).get("warnings", []):
        print(f"  !! {warning}")

    if rebuild:
        problems.extend(_check_rebuild(bundle))
    return problems


def _check_patch(bundle: dict) -> list[str]:
    """A surface bundle's patch must be internally consistent and exposed."""
    patch = bundle.get("patch")
    if not patch:
        if bundle_mod.binding_mode(bundle) == "surface":
            return ["surface-mode bundle carries no measured patch"]
        return []

    problems = []
    per_residue = patch.get("per_residue_area_A2", {})
    total = sum(float(v) for v in per_residue.values())
    stored = float(patch["patch_area_A2"])
    if abs(total - stored) > TOLERANCE * max(1.0, abs(stored)):
        problems.append(
            f"patch area {stored} does not equal the sum of its per-residue "
            f"areas ({total})"
        )
    if len(per_residue) != patch["n_residues"]:
        problems.append(
            f"patch reports {patch['n_residues']} residues but stores "
            f"{len(per_residue)} per-residue areas"
        )
    print(f"  patch           {patch['n_residues']} residues, {stored:.0f} A^2, "
          f"planarity {patch['planarity_A']:.1f} A")
    if patch.get("buried_residue_numbers"):
        # Configured residues with no exposed surface cannot be part of a
        # binding face; almost always a numbering or chain mistake.
        print(f"  !! buried binding-site residues: {patch['buried_residue_numbers']}")
    if patch.get("shape_warning"):
        print("  !! patch is not flat; surface-mode agreement assumes a "
              "roughly planar face")
    return problems


def _check_rebuild(bundle: dict) -> list[str]:
    """Rebuild from the same inputs; the scientific payload must be identical."""
    from aptarank.config import load_config
    from aptarank.tier2.build import build_target_bundle

    import tempfile

    target = bundle["target"]
    applied = bundle["preparation"].get("applied", {})
    # A scratch directory, not the working tree: rebuilding is a check, and a
    # check should not leave a second copy of the evidence lying next to the
    # first one for someone to pick up by mistake.
    rebuild_dir = tempfile.mkdtemp(prefix="aptarank-rebuild-")
    cfg = load_config(
        overrides={
            "tier2": {
                "binding_mode": bundle_mod.binding_mode(bundle),
                "target": {
                    "id": target.get("identifier", target["pdb_id"]),
                    "source": target.get("target_source", "pdb"),
                    "chain": target["chain_id"],
                    "model": target["model_index"],
                    "partner_chains": applied.get("partner_chains_removed", []),
                    "strip_hetatm": bool(applied.get("strip_hetatm", False)),
                },
                "bundle_dir": rebuild_dir,
            }
        }
    )
    requested = [
        r["residue_number"]
        for r in bundle["preparation"].get("site_residues", [])
    ] or [
        r["residue_number"]
        for r in bundle["selection"].get("target_site", {}).get("requested_residues", [])
    ]
    if requested:
        cfg = cfg.with_overrides(
            {"tier2": {"target": {"target_site_residues": requested}}}
        )
    retained = sorted({
        r["residue_name"]
        for r in bundle["preparation"]["applied"].get("retained_hetero_residues", [])
    })
    if retained:
        cfg = cfg.with_overrides({"tier2": {"target": {"retain_hetero_resnames": retained}}})

    rebuilt, _path = build_target_bundle(cfg)
    if rebuilt["bundle_id"] != bundle["bundle_id"]:
        return [
            f"rebuild produced bundle_id {rebuilt['bundle_id']} but the original is "
            f"{bundle['bundle_id']}: the build is not reproducible"
        ]
    print("  rebuild         identical scientific payload")
    return _check_volume_drift(bundle, rebuilt)


def _check_volume_drift(original: dict, rebuilt: dict, tolerance: float = 0.10) -> list[str]:
    """How far fpocket's Monte-Carlo volume moved between two identical builds.

    The bundle id deliberately excludes this number (see bundle.py), so
    something has to keep an eye on it: a few percent is the tool, and a large
    jump means the two builds did not measure the same cavity after all.
    """
    problems, drifts = [], []
    rebuilt_by_index = {p["index"]: p for p in rebuilt["pockets"]}
    for pocket in original["pockets"]:
        other = rebuilt_by_index.get(pocket["index"])
        if other is None:
            problems.append(f"rebuild is missing pocket {pocket['index']}")
            continue
        before = float(pocket["fpocket"]["volume_A3"])
        after = float(other["fpocket"]["volume_A3"])
        if before <= 0:
            continue
        drift = abs(after - before) / before
        drifts.append(drift)
        if drift > tolerance:
            problems.append(
                f"pocket {pocket['index']}: volume moved {drift:.1%} between two "
                f"identical builds ({before:.1f} -> {after:.1f} A^3), beyond the "
                f"{tolerance:.0%} expected from fpocket's Monte-Carlo estimate"
            )
    if drifts:
        print(f"  volume drift    max {max(drifts):.1%} across {len(drifts)} pockets "
              f"(Monte-Carlo; excluded from bundle_id by design)")
    return problems


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="+")
    parser.add_argument("--rebuild", action="store_true",
                        help="rebuild the bundle and require an identical bundle_id")
    args = parser.parse_args(argv)

    failures = 0
    for pattern in args.paths:
        # glob.glob rather than Path.glob: the documented invocation passes an
        # absolute pattern (~/aptarank-data/cache/targets/*.bundle.json), which
        # Path.glob refuses outright.
        matches = [Path(p) for p in sorted(glob.glob(str(Path(pattern).expanduser())))]
        for path in matches or [Path(pattern).expanduser()]:
            print(f"\n{path}")
            problems = verify(path, rebuild=args.rebuild)
            for problem in problems:
                print(f"  FAIL  {problem}", file=sys.stderr)
            failures += len(problems)

    print(f"\n{'FAILED' if failures else 'OK'}: {failures} problem(s)")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
