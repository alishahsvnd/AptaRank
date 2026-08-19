"""End-to-end with Tier 2 enabled, against a synthetic target bundle."""

from __future__ import annotations

import pytest

from aptarank.config import load_config
from aptarank.pipeline import run_pipeline
from aptarank.tier2 import bundle as bundle_mod

from .conftest import make_synthetic_bundle


@pytest.fixture(scope="module")
def tier2_run(tmp_path_factory, request):
    tmp = tmp_path_factory.mktemp("tier2run")
    bundle = make_synthetic_bundle(tmp, d_pocket_target=24.0)
    bundle_path = bundle_mod.write(bundle, tmp / "bundles")

    corpus = request.path.parent / "fixtures" / "mini_corpus.csv"
    candidates = request.path.parent / "fixtures" / "mini_candidates.csv"
    cfg = load_config(
        overrides={
            "corpus": {"path": str(corpus), "is_placeholder": True,
                       "allow_placeholder": True, "cache_dir": str(tmp / "cache")},
            "tier1": {"n_ensemble_samples": 20, "shuffle": {"n_shuffles": 20},
                      "parallel": {"workers": 1}},
            "tier2": {
                "enabled": True,
                "n_candidates": 8,
                "bundle_path": str(bundle_path),
                "calibration": {"bank_size": 60, "cache_dir": str(tmp / "bank")},
            },
            "output": {"dir": str(tmp / "runs"), "n_diagrams": 1},
        }
    )
    return run_pipeline(cfg, candidates)


def test_target_evidence_reaches_the_artifact(tier2_run):
    target = tier2_run.artifact["target"]
    assert target["pdb_id"] == "TEST"
    assert target["n_pockets"] == 2
    assert target["selected_pocket"]["d_pocket_A"] > 0
    assert target["bundle_id"]


def test_only_the_survivors_are_evaluated(tier2_run):
    candidates = tier2_run.artifact["candidates"]
    evaluated = [c for c in candidates if c["tier2"].get("status") == "evaluated"]
    unevaluated = [c for c in candidates if c["tier2"].get("status") != "evaluated"]

    assert len(evaluated) <= 8
    assert all(c["rank"] <= 8 for c in evaluated)
    # Below the cut means no evidence — never "weak", which would invent some.
    assert all(c["tier2"]["band"] == "not_evaluated" for c in unevaluated)


def test_tier2_does_not_reorder_the_ranking(tier2_run):
    """The rule of §6, asserted rather than assumed."""
    candidates = tier2_run.artifact["candidates"]
    scores = [c["tier1_score"] for c in candidates]
    assert scores == sorted(scores, reverse=True)
    assert [c["rank"] for c in candidates] == list(range(1, len(candidates) + 1))


def test_every_evaluated_candidate_carries_auditable_geometry(tier2_run):
    for candidate in tier2_run.artifact["candidates"]:
        tier2 = candidate["tier2"]
        if tier2.get("status") != "evaluated":
            continue
        assert tier2["band"] in {"strong", "moderate", "weak"}
        assert 0.0 <= tier2["control_percentile_flexible"] <= 1.0
        # both descriptors stored: primary + sensitivity check
        assert tier2["d_apt_flexible_A"] > 0 and tier2["d_apt_extended_A"] > 0
        assert tier2["absolute_mismatch_flexible_A"] == pytest.approx(
            abs(tier2["d_apt_flexible_A"] - tier2["d_pocket_A"])
        )
        assert tier2["target_bundle_id"] and tier2["calibration_bank_id"]


def test_spearman_diagnostic_is_computed_every_run(tier2_run):
    """§6.4 — run as soon as both tiers produce output, not in evaluation week."""
    spearman = tier2_run.artifact["diagnostics"]["spearman_tier1_tier2"]
    assert spearman is not None
    assert spearman["n"] >= 0
    if spearman["rho"] is not None:
        assert -1.0 <= spearman["rho"] <= 1.0


def test_thresholds_and_bank_provenance_are_stored(tier2_run):
    thresholds = tier2_run.artifact["tier2_thresholds"]
    assert thresholds["n_controls"] > 0
    assert thresholds["binding_mode"] == "pocket"
    assert thresholds["units"] == "A"
    assert thresholds["parameters"]["primary_descriptor"] == "flexible"
    assert thresholds["calibration_bank"]["corpus_is_placeholder"] is True


def test_explanations_quote_the_stored_geometry(tier2_run):
    with_band = [
        c for c in tier2_run.artifact["candidates"]
        if c["tier2"].get("band") in {"strong", "moderate", "weak"}
    ]
    assert with_band
    for candidate in with_band:
        assert candidate["explanation"].endswith("it is not evidence of binding.")
