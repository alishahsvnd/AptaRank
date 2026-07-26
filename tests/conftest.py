from __future__ import annotations

from pathlib import Path

import pytest

FIXTURES = Path(__file__).parent / "fixtures"

# The worked example from spec §2.1, with the values verified against
# ViennaRNA 2.7.2. Used as a regression fixture throughout.
SPEC_SEQUENCE = "GUUCCAUGGGCCUUGACUUGCUGUGUCAUCACCAUGGGAC"
SPEC_STRUCTURE = "(((((((((....((((.......))))...)))))))))"
SPEC_MFE = -17.9
SPEC_MFE_NORM = -0.4475
SPEC_ENSEMBLE_DEFECT = 0.0092


@pytest.fixture(scope="session")
def mini_corpus_path() -> Path:
    return FIXTURES / "mini_corpus.csv"


@pytest.fixture(scope="session")
def mini_candidates_path() -> Path:
    return FIXTURES / "mini_candidates.csv"


def make_synthetic_bundle(tmp_path, d_pocket_target: float = 24.0):
    """A valid target bundle built without fpocket.

    fpocket has no Windows build, so everything downstream of the bundle is
    developed and tested against synthetic bundles like this one. The real
    tool's output is validated separately by the pinned Linux CI job.
    """
    import numpy as np

    from aptarank.tier2 import bundle as bundle_mod
    from aptarank.tier2.fpocket import Pocket, Residue
    from aptarank.tier2.geometry import pocket_geometry

    # A rod of alpha spheres whose robust envelope is close to d_pocket_target.
    n = 40
    half = (d_pocket_target - 6.0) / 2.0
    xs = np.linspace(-half, half, n) / 0.9  # widen for the 5-95% quantile trim
    centres = np.column_stack(
        [xs, np.sin(np.arange(n)) * 2.0, np.cos(np.arange(n)) * 2.0]
    )
    radii = np.full(n, 3.0)

    pockets = []
    for index, (score, scale) in enumerate([(0.62, 1.0), (0.31, 0.4)], start=1):
        pocket = Pocket(
            index=index,
            metrics={
                "Score": score,
                "Druggability Score": 0.5,
                "Volume": 400.0 * scale,
                "Number of Alpha Spheres": n,
            },
            lining_residues=[Residue("A", r, "", "HIS") for r in (120, 122, 124)],
        )
        from aptarank.tier2.fpocket import AlphaSphere

        pocket.alpha_spheres = [
            AlphaSphere(i, tuple(c * scale), float(r), "apolar")
            for i, (c, r) in enumerate(zip(centres, radii))
        ]
        pockets.append(pocket)

    geometries = {
        p.index: pocket_geometry(
            [s.center_A for s in p.alpha_spheres],
            [s.radius_A for s in p.alpha_spheres],
            p.volume_A3,
        )
        for p in pockets
    }

    clean = tmp_path / "TEST_A_clean.pdb"
    clean.write_text("ATOM      1  CA  HIS A 120       0.000   0.000   0.000\n", encoding="utf-8")

    class _Prepared:
        pdb_id = "TEST"
        name = "synthetic test target"
        model_index = 0
        chain_id = "A"
        path = clean
        source = {"url": "synthetic", "format": "pdb", "sha256": "0" * 64, "size_bytes": 1}
        applied = {
            "input_atom_count": 1, "output_atom_count": 1,
            "output_protein_residue_count": 1, "removed_water_residue_count": 0,
            "hetero_summary": [], "retained_hetero_residues": [],
            "altloc_policy": "first_altloc_only", "warnings": [],
        }

    selection_payload = {
        "status": "selected",
        "method": "active_site_overlap",
        "selected_pocket_index": 1,
        "active_site": {"requested": True, "allow_zero_overlap_fallback": False,
                        "requested_residues": [], "n_requested": 3, "total_overlap": 3},
        "pocket_evidence": [
            {"pocket_index": p.index, "overlap_count": 3 if p.index == 1 else 0,
             "overlapping_residues": [], "fpocket_score": p.score,
             "druggability_score": p.druggability, "selected": p.index == 1}
            for p in pockets
        ],
        "tie_break_order": [],
        "warnings": [],
    }

    return bundle_mod.build(
        prepared=_Prepared(),
        pockets=pockets,
        geometries=geometries,
        selection=selection_payload,
        fpocket_provenance={"status": "success", "version": "fpocket 4.1 (synthetic)",
                            "command": ["fpocket", "-f", "TEST_A_clean.pdb"], "exit_code": 0},
    )


@pytest.fixture
def synthetic_bundle(tmp_path):
    return make_synthetic_bundle(tmp_path)


@pytest.fixture
def dev_config(tmp_path, mini_corpus_path):
    """A fast, self-contained config pointed at the placeholder corpus."""
    from aptarank.config import load_config

    return load_config(
        overrides={
            "corpus": {
                "path": str(mini_corpus_path),
                "is_placeholder": True,
                "allow_placeholder": True,
                "cache_dir": str(tmp_path / "corpus_cache"),
            },
            "tier1": {
                "n_ensemble_samples": 20,
                "shuffle": {"n_shuffles": 20},
                "parallel": {"workers": 1},
            },
            "output": {"dir": str(tmp_path / "runs"), "n_diagrams": 2},
        }
    )
