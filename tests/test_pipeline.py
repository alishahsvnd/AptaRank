"""End-to-end: candidate file in, valid run artifact out."""

from __future__ import annotations

import json

import pytest

from aptarank.artifacts import read_artifact
from aptarank.config import load_config
from aptarank.errors import CorpusError
from aptarank.pipeline import run_pipeline


@pytest.fixture(scope="module")
def run_output(tmp_path_factory, request):
    tmp = tmp_path_factory.mktemp("run")
    corpus = request.path.parent / "fixtures" / "mini_corpus.csv"
    candidates = request.path.parent / "fixtures" / "mini_candidates.csv"
    cfg = load_config(
        overrides={
            "corpus": {"path": str(corpus), "is_placeholder": True,
                       "allow_placeholder": True, "cache_dir": str(tmp / "cache")},
            "tier1": {"n_ensemble_samples": 20, "shuffle": {"n_shuffles": 20},
                      "parallel": {"workers": 1}},
            "output": {"dir": str(tmp / "runs"), "n_diagrams": 2},
        }
    )
    return run_pipeline(cfg, candidates)


def test_artifact_is_written_and_readable(run_output):
    assert run_output.path.exists()
    assert read_artifact(run_output.path)["run_id"] == run_output.artifact["run_id"]


def test_artifact_is_valid_json_without_nan(run_output):
    """NaN is not valid JSON; unavailable values must be null."""
    text = run_output.path.read_text(encoding="utf-8")
    assert "NaN" not in text and "Infinity" not in text
    json.loads(text)


def test_placeholder_corpus_marks_the_run_ineligible(run_output):
    assert run_output.artifact["run_mode"] == "development"
    assert run_output.artifact["publication_eligible"] is False
    assert run_output.artifact["corpus"]["is_placeholder"] is True


def test_ranks_are_dense_and_ordered_by_tier1_score(run_output):
    candidates = run_output.artifact["candidates"]
    assert [c["rank"] for c in candidates] == list(range(1, len(candidates) + 1))
    scores = [c["tier1_score"] for c in candidates]
    assert scores == sorted(scores, reverse=True)


def test_every_candidate_carries_decomposable_evidence(run_output):
    for candidate in run_output.artifact["candidates"]:
        assert set(candidate["criteria"]) == {
            "mfe_norm", "ensemble_defect", "positional_entropy_mean",
            "stem_fraction", "gc_fraction",
        }
        for entry in candidate["criteria"].values():
            assert entry["value"] is not None
            assert 0.0 <= entry["score"] <= 1.0
        assert candidate["structure"]["dot_bracket"]
        assert candidate["explanation"].endswith(
            "it is not evidence of binding."
        )


def test_tier2_absent_means_not_evaluated_never_weak(run_output):
    bands = {c["tier2"]["band"] for c in run_output.artifact["candidates"]}
    statuses = {c["tier2"]["status"] for c in run_output.artifact["candidates"]}
    assert bands == {"not_evaluated"}
    # Tier 2 was never run here; that is a different fact from "below the cut".
    assert statuses == {"not_run"}


def test_diagrams_are_embedded_for_the_top_n_only(run_output):
    with_svg = [c for c in run_output.artifact["candidates"] if c["structure"]["svg"]]
    assert len(with_svg) == 2
    assert all(c["rank"] <= 2 for c in with_svg)
    assert with_svg[0]["structure"]["svg"].lstrip().startswith("<?xml")


def test_provenance_is_recorded(run_output):
    artifact = run_output.artifact
    assert artifact["versions"]["viennarna"]
    assert artifact["input"]["sha256"]
    assert artifact["corpus"]["corpus_sha256"]
    assert artifact["scoring_signature"]


def test_running_without_a_corpus_is_blocked(tmp_path, request):
    cfg = load_config(overrides={"output": {"dir": str(tmp_path)}})
    candidates = request.path.parent / "fixtures" / "mini_candidates.csv"
    with pytest.raises(CorpusError, match="no reference corpus configured"):
        run_pipeline(cfg, candidates, write=False)


def test_placeholder_corpus_without_permission_is_blocked(tmp_path, request):
    corpus = request.path.parent / "fixtures" / "mini_corpus.csv"
    candidates = request.path.parent / "fixtures" / "mini_candidates.csv"
    cfg = load_config(
        overrides={
            "corpus": {"path": str(corpus), "is_placeholder": True,
                       "allow_placeholder": False, "cache_dir": str(tmp_path / "c")},
            "output": {"dir": str(tmp_path / "runs")},
        }
    )
    with pytest.raises(CorpusError, match="Refusing to score"):
        run_pipeline(cfg, candidates, write=False)
