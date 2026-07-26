"""Write a SYNTHETIC target bundle so Tier 2 can run without fpocket.

fpocket has no Windows build. This produces a bundle with a fabricated cavity
of a chosen size so the dashboard, the demo and the tests exercise the full
two-tier path on any machine.

The bundle is flagged `synthetic: true`, which propagates into the run artifact
and marks it `publication_eligible: false`. Real bundles come from the pinned
Linux CI workflow (.github/workflows/target-bundle.yml).

    python scripts/make_dev_target.py [--d-pocket 24.0] [--pdb-id DEMO]
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from aptarank.tier2 import bundle as bundle_mod
from aptarank.tier2.fpocket import AlphaSphere, Pocket, Residue
from aptarank.tier2.geometry import pocket_geometry

REPO_ROOT = Path(__file__).resolve().parent.parent


def synthetic_pockets(d_pocket_target: float) -> list[Pocket]:
    """Two cavities: one elongated at roughly the requested width, one small."""
    rng = np.random.default_rng(11)
    pockets = []
    for index, (scale, score, volume, n) in enumerate(
        [(1.0, 0.612, 412.0, 48), (0.45, 0.298, 96.0, 22)], start=1
    ):
        half = max(1.0, (d_pocket_target * scale - 6.0) / 2.0) / 0.9
        xs = np.linspace(-half, half, n)
        centres = np.column_stack(
            [xs, rng.normal(0, 1.6, n), rng.normal(0, 1.6, n)]
        )
        pocket = Pocket(
            index=index,
            metrics={
                "Score": score,
                "Druggability Score": 0.78 if index == 1 else 0.10,
                "Volume": volume,
                "Number of Alpha Spheres": n,
                "Mean alpha sphere radius": 3.5,
            },
            lining_residues=[
                Residue("A", r, "", name)
                for r, name in (
                    [(120, "HIS"), (122, "HIS"), (189, "ASP"), (208, "CYS"), (301, "ZN")]
                    if index == 1
                    else [(45, "LEU"), (46, "VAL"), (77, "SER")]
                )
            ],
        )
        pocket.alpha_spheres = [
            AlphaSphere(i, tuple(float(v) for v in c), 3.5, "apolar" if i % 3 else "polar")
            for i, c in enumerate(centres)
        ]
        pockets.append(pocket)
    return pockets


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdb-id", default="DEMO")
    parser.add_argument("--chain", default="A")
    parser.add_argument("--d-pocket", type=float, default=24.0,
                        help="approximate characteristic cavity width, Å")
    parser.add_argument("-o", "--output-dir", default=str(REPO_ROOT / "cache" / "targets"))
    args = parser.parse_args()

    pockets = synthetic_pockets(args.d_pocket)
    geometries = {
        p.index: pocket_geometry(
            [s.center_A for s in p.alpha_spheres],
            [s.radius_A for s in p.alpha_spheres],
            p.volume_A3,
        )
        for p in pockets
    }

    work = Path(args.output_dir) / "work"
    work.mkdir(parents=True, exist_ok=True)
    clean = work / f"{args.pdb_id}_{args.chain}_clean.pdb"
    clean.write_text(
        "REMARK  SYNTHETIC PLACEHOLDER STRUCTURE - NOT A REAL PDB ENTRY\n"
        "ATOM      1  CA  HIS A 120       0.000   0.000   0.000  1.00 20.00           C\n"
        "END\n",
        encoding="utf-8",
    )

    class _Prepared:
        pdb_id = args.pdb_id.upper()
        name = "SYNTHETIC development target — not a real structure"
        model_index = 0
        chain_id = args.chain
        path = clean
        source = {"url": "synthetic://development", "format": "pdb",
                  "sha256": "0" * 64, "size_bytes": clean.stat().st_size}
        applied = {
            "input_atom_count": 1, "output_atom_count": 1,
            "output_protein_residue_count": 1, "removed_water_residue_count": 0,
            "hetero_summary": [], "retained_hetero_residues": [
                {"chain_id": "A", "residue_number": 301, "insertion_code": "",
                 "residue_name": "ZN", "atom_count": 1}
            ],
            "altloc_policy": "first_altloc_only",
            "warnings": ["synthetic structure; no real atoms were parsed"],
        }

    selection = {
        "status": "selected",
        "method": "active_site_overlap",
        "selected_pocket_index": 1,
        "active_site": {
            "requested": True, "allow_zero_overlap_fallback": False,
            "requested_residues": [
                {"chain_id": "A", "residue_number": r, "insertion_code": "", "residue_name": None}
                for r in (120, 122, 189)
            ],
            "n_requested": 3, "total_overlap": 3,
        },
        "pocket_evidence": [
            {"pocket_index": p.index, "overlap_count": 3 if p.index == 1 else 0,
             "overlapping_residues": [], "fpocket_score": p.score,
             "druggability_score": p.druggability, "selected": p.index == 1}
            for p in pockets
        ],
        "tie_break_order": ["overlap_count_desc", "fpocket_score_desc",
                            "druggability_score_desc", "pocket_index_asc"],
        "warnings": ["SYNTHETIC bundle: this cavity was fabricated, not detected"],
    }

    result = bundle_mod.build(
        prepared=_Prepared(),
        pockets=pockets,
        geometries=geometries,
        selection=selection,
        fpocket_provenance={
            "status": "synthetic", "version": "none (no fpocket)",
            "command": [], "exit_code": None,
        },
        synthetic=True,
    )
    path = bundle_mod.write(result, args.output_dir)
    selected = bundle_mod.selected_pocket(result)
    print(
        f"wrote SYNTHETIC bundle {path}\n"
        f"  d_pocket = {selected['geometry']['d_pocket_A']:.2f} A "
        f"(requested ~{args.d_pocket})\n"
        f"  runs using it are stamped publication_eligible: false"
    )


if __name__ == "__main__":
    main()
