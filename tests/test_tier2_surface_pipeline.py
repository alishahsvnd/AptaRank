"""End to end in surface mode, and the guards around mode/bundle mismatch."""

from __future__ import annotations

import pytest

from aptarank.config import load_config
from aptarank.errors import TargetError
from aptarank.pipeline import run_pipeline
from aptarank.tier2 import bundle as bundle_mod

from .conftest import make_surface_bundle, make_synthetic_bundle


@pytest.fixture(scope="module")
def surface_run(tmp_path_factory, request):
    tmp = tmp_path_factory.mktemp("surfacerun")
    bundle = make_surface_bundle(tmp, patch_area_A2=1440.0)
    bundle_path = bundle_mod.write(bundle, tmp / "targets")

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
                "binding_mode": "surface",
                "n_candidates": 8,
                "bundle_path": str(bundle_path),
                "calibration": {"bank_size": 60, "cache_dir": str(tmp / "bank")},
            },
            "output": {"dir": str(tmp / "runs"), "n_diagrams": 1},
        }
    )
    return run_pipeline(cfg, candidates)


def test_the_artifact_records_which_geometry_was_compared(surface_run):
    artifact = surface_run.artifact
    assert artifact["binding_mode"] == "surface"
    assert artifact["target"]["binding_mode"] == "surface"
    assert artifact["tier2_thresholds"]["units"] == "A^2"


def test_surface_candidates_carry_area_evidence(surface_run):
    evaluated = [
        c for c in surface_run.artifact["candidates"]
        if (c["tier2"] or {}).get("status") == "evaluated"
    ]
    assert evaluated
    for record in (c["tier2"] for c in evaluated):
        assert record["binding_mode"] == "surface"
        assert record["footprint_area_A2"] > 0
        assert record["patch_area_A2"] == pytest.approx(1440.0)
        assert record["disagreement_units"] == "A^2"
        # The composite is reported, and says plainly that its charge term is a
        # property of the target rather than of the candidate.
        assert 0 < record["geometric_agreement_score"] <= 1
        assert record["charge_is_target_level"] is True


def test_surface_mode_still_never_reorders_the_ranking(surface_run):
    ranks = [c["rank"] for c in surface_run.artifact["candidates"]]
    scores = [c["tier1_score"] for c in surface_run.artifact["candidates"]]
    assert ranks == sorted(ranks)
    assert scores == sorted(scores, reverse=True)


def test_explanations_describe_the_mode_that_actually_ran(surface_run):
    banded = [
        c for c in surface_run.artifact["candidates"]
        if (c["tier2"] or {}).get("band") in ("strong", "moderate", "weak")
    ]
    assert banded
    text = " ".join(c["explanation"] for c in banded)
    # Surface mode talks about covering an area, never about a loop in a cavity.
    assert "Å²" in text or "A^2" in text
    assert "cavity detected on the target" not in text


def test_a_bundle_built_for_one_mode_is_refused_by_the_other(tmp_path, request):
    """Otherwise the run would compare against whatever number was present."""
    pocket_bundle = make_synthetic_bundle(tmp_path)
    path = bundle_mod.write(pocket_bundle, tmp_path / "targets")
    corpus = request.path.parent / "fixtures" / "mini_corpus.csv"
    candidates = request.path.parent / "fixtures" / "mini_candidates.csv"

    cfg = load_config(
        overrides={
            "corpus": {"path": str(corpus), "is_placeholder": True,
                       "allow_placeholder": True, "cache_dir": str(tmp_path / "cache")},
            "tier1": {"n_ensemble_samples": 5, "shuffle": {"enabled": False},
                      "parallel": {"workers": 1}},
            "tier2": {"enabled": True, "binding_mode": "surface",
                      "bundle_path": str(path),
                      "calibration": {"bank_size": 20,
                                      "cache_dir": str(tmp_path / "bank")}},
            "output": {"dir": str(tmp_path / "runs"), "embed_svg": False},
        }
    )
    with pytest.raises(TargetError, match="prepared for 'pocket' mode"):
        run_pipeline(cfg, candidates, write=False)


def test_a_surface_bundle_without_a_patch_is_invalid(tmp_path):
    bundle = make_surface_bundle(tmp_path)
    bundle["patch"] = None
    with pytest.raises(TargetError, match="no measured patch area"):
        bundle_mod.validate(bundle)


def test_v1_bundles_are_still_readable_as_pocket_targets(tmp_path):
    """Bundles built before binding modes existed only ever had cavities."""
    bundle = make_synthetic_bundle(tmp_path)
    bundle.pop("binding_mode")
    bundle["schema_version"] = bundle_mod.LEGACY_SCHEMA_VERSION
    bundle_mod.validate(bundle)
    assert bundle_mod.binding_mode(bundle) == "pocket"


def test_the_mode_is_in_the_filename_so_both_can_coexist(tmp_path):
    pocket_path = bundle_mod.write(make_synthetic_bundle(tmp_path), tmp_path / "t")
    surface_path = bundle_mod.write(make_surface_bundle(tmp_path), tmp_path / "t")
    assert "_pocket_" in pocket_path.name
    assert "_surface_" in surface_path.name
    assert bundle_mod.find(tmp_path / "t", "TEST", "A", mode="surface") == surface_path
    assert bundle_mod.find(tmp_path / "t", "TEST", "A", mode="pocket") == pocket_path


def test_a_prepared_target_carries_its_own_mode_into_the_run(tmp_path):
    """Picking a prepared target is the assertion; the config default must not
    override it.

    The bug: choosing the prepared 7WRQ surface target in the dashboard launched
    `--target-bundle <path>` with no mode, so the default (pocket) applied and
    the run died with "the run asks for 'pocket'" — a mode nobody had asked for.
    """
    from aptarank.cli import build_parser, _config_from_args

    surface_path = bundle_mod.write(make_surface_bundle(tmp_path), tmp_path / "t")
    pocket_path = bundle_mod.write(make_synthetic_bundle(tmp_path), tmp_path / "t")

    for path, expected in ((surface_path, "surface"), (pocket_path, "pocket")):
        args = build_parser().parse_args(
            ["run", "candidates.csv", "--target-bundle", str(path)]
        )
        assert _config_from_args(args).get("tier2.binding_mode") == expected


def test_an_explicit_mode_still_wins_and_still_has_to_agree(tmp_path):
    """Adopting the bundle's mode must not swallow a genuine conflict."""
    from aptarank.cli import build_parser, _config_from_args

    surface_path = bundle_mod.write(make_surface_bundle(tmp_path), tmp_path / "t")
    args = build_parser().parse_args(
        ["run", "candidates.csv", "--target-bundle", str(surface_path),
         "--binding-mode", "pocket"]
    )
    cfg = _config_from_args(args)
    assert cfg.get("tier2.binding_mode") == "pocket"

    from aptarank.tier2.service import _binding_mode

    with pytest.raises(TargetError, match="prepared for 'surface' mode"):
        _binding_mode(cfg, bundle_mod.load(surface_path))


def test_an_unreadable_bundle_is_left_to_the_pipeline_to_diagnose(tmp_path):
    from aptarank.cli import build_parser, _config_from_args

    broken = tmp_path / "broken.bundle.json"
    broken.write_text("{not json", encoding="utf-8")
    args = build_parser().parse_args(
        ["run", "candidates.csv", "--target-bundle", str(broken)]
    )
    # Falls back to the default rather than crashing in argument parsing.
    assert _config_from_args(args).get("tier2.binding_mode") == "pocket"


def test_the_dashboard_launches_a_prepared_target_with_its_mode(tmp_path):
    from dashboard import jobs
    from dashboard.inputs import build_target_request, inspect_target

    path = bundle_mod.write(make_surface_bundle(tmp_path), tmp_path / "targets")
    request = build_target_request("prepared", prepared=inspect_target(path))
    assert request.binding_mode == "surface"

    job = jobs.submit(
        tmp_path / "runs",
        candidates_path=tmp_path / "c.csv",
        corpus_path=tmp_path / "lib.csv",
        target=request,
        python_executable="/usr/bin/true",
    )
    command = job.request["command"]
    assert "--binding-mode" in command
    assert command[command.index("--binding-mode") + 1] == "surface"
