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
import sys
from pathlib import Path

from aptarank.tier2 import bundle as bundle_mod
from aptarank.tier2.geometry import pocket_geometry

TOLERANCE = 1e-6


def verify(path: Path, rebuild: bool = False) -> list[str]:
    problems: list[str] = []
    bundle = bundle_mod.load(path)          # also re-checks the bundle id
    print(f"  bundle_id       {bundle['bundle_id']}")
    print(f"  target          {bundle['target']['pdb_id']} chain {bundle['target']['chain_id']}")
    print(f"  pockets         {len(bundle['pockets'])}")
    print(f"  selection       {bundle['selection']['method']} "
          f"-> pocket {bundle['selection']['selected_pocket_index']}")

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

    selected = bundle_mod.selected_pocket(bundle)
    if selected["geometry"]["shape_warning"]:
        print("  !! selected cavity is oddly shaped (envelope/equivalent-sphere "
              "ratio out of range); the geometric comparison assumes a roughly "
              "convex pocket")
    if bundle["selection"]["method"] != "active_site_overlap":
        print(f"  !! pocket was NOT selected by active-site overlap "
              f"({bundle['selection']['method']}); the dashboard must caveat this")
    for warning in bundle["selection"].get("warnings", []):
        print(f"  !! {warning}")

    if rebuild:
        problems.extend(_check_rebuild(bundle))
    return problems


def _check_rebuild(bundle: dict) -> list[str]:
    """Rebuild from the same inputs; the scientific payload must be identical."""
    from aptarank.config import load_config
    from aptarank.tier2.build import build_target_bundle

    target = bundle["target"]
    cfg = load_config(
        overrides={
            "tier2": {
                "target": {
                    "pdb_id": target["pdb_id"],
                    "chain": target["chain_id"],
                    "model": target["model_index"],
                },
                "bundle_dir": "bundles/_rebuild",
            }
        }
    )
    requested = bundle["selection"]["active_site"].get("requested_residues", [])
    if requested:
        cfg = cfg.with_overrides(
            {"tier2": {"target": {"active_site_residues": requested}}}
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
    return []


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="+")
    parser.add_argument("--rebuild", action="store_true",
                        help="rebuild the bundle and require an identical bundle_id")
    args = parser.parse_args(argv)

    failures = 0
    for pattern in args.paths:
        for path in sorted(Path().glob(pattern)) or [Path(pattern)]:
            print(f"\n{path}")
            problems = verify(path, rebuild=args.rebuild)
            for problem in problems:
                print(f"  FAIL  {problem}", file=sys.stderr)
            failures += len(problems)

    print(f"\n{'FAILED' if failures else 'OK'}: {failures} problem(s)")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
