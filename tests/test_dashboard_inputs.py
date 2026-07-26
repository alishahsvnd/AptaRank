"""Input discovery and the gate that decides whether a run may proceed."""

from __future__ import annotations

import json
import shutil

import pytest

from dashboard.inputs import (
    DEVELOPMENT,
    UNVERIFIED,
    VERIFIED,
    discover_libraries,
    discover_targets,
    inspect_library,
    review,
)

CORPUS_HEADER = "id,sequence,target_name,target_pdb_id,source_reference\n"


def write_corpus(path, n=200, name="lib"):
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = [
        f"{name}{i},{'ACGU' * 8},NDM-1,3SPU,ref\n" for i in range(n)
    ]
    path.write_text(CORPUS_HEADER + "".join(rows), encoding="utf-8")
    return path


def test_the_same_library_in_two_places_is_listed_once(tmp_path):
    """A corpus commonly exists in both the checkout and the data directory;
    two identical entries give the user no way to choose between them."""
    first = write_corpus(tmp_path / "data" / "corpus" / "validated.csv")
    second = tmp_path / "checkout" / "data" / "corpus" / "validated.csv"
    second.parent.mkdir(parents=True)
    shutil.copy(first, second)

    found = discover_libraries(first.parent, second.parent)
    assert len(found) == 1


def test_different_libraries_with_the_same_name_are_both_listed(tmp_path):
    a = write_corpus(tmp_path / "a" / "validated.csv", n=200, name="x")
    b = write_corpus(tmp_path / "b" / "validated.csv", n=300, name="y")
    found = discover_libraries(a.parent, b.parent)
    assert len(found) == 2


def test_provenance_state_is_read_from_a_manifest(tmp_path):
    path = write_corpus(tmp_path / "validated.csv")
    assert inspect_library(path).state == UNVERIFIED

    path.with_suffix(".manifest.json").write_text(
        json.dumps({"source": "SELEX literature", "curator": "Laura",
                    "curated_date": "2026-08-01"}),
        encoding="utf-8",
    )
    assert inspect_library(path).state == VERIFIED


def test_synthetic_data_is_recognised_and_sorted_last(tmp_path):
    write_corpus(tmp_path / "dev_placeholder_corpus.csv")
    write_corpus(tmp_path / "validated.csv", n=250, name="v")
    found = discover_libraries(tmp_path)
    assert [lib.state for lib in found] == [UNVERIFIED, DEVELOPMENT]


def test_a_candidate_file_is_not_offered_as_a_reference_library(tmp_path):
    """Calibrating a user's scores against their own unvalidated sequences
    would be silently meaningless."""
    path = tmp_path / "candidates.csv"
    path.write_text("id,sequence\nc1,ACGUACGUACGUACGUACGUACGU\n", encoding="utf-8")
    library = inspect_library(path)
    assert not library.usable
    assert "candidate sequences rather than a reference library" in library.problem


def test_a_tiny_library_is_rejected(tmp_path):
    library = inspect_library(write_corpus(tmp_path / "small.csv", n=10))
    assert not library.usable
    assert "distribution" in library.problem


def test_an_unreadable_target_bundle_is_reported_not_crashed(tmp_path):
    (tmp_path / "BROKEN_A_deadbeef.bundle.json").write_text("{not json", encoding="utf-8")
    found = discover_targets(tmp_path)
    assert len(found) == 1 and not found[0].usable


# -- the gate ------------------------------------------------------------


def ok_validation(**overrides):
    return {"ok": True, "n_submitted": 100, "n_valid": 100, "n_rejected": 0, **overrides}


def test_a_rigorous_run_on_synthetic_data_is_refused_not_downgraded(tmp_path):
    library = inspect_library(write_corpus(tmp_path / "dev_placeholder_corpus.csv"))
    verdict = review(ok_validation(), library, None, preset="evaluation")
    assert not verdict["can_run"]
    assert any("rigorous" in r for r in verdict["refusals"])

    # The same inputs at the standard setting are allowed, but flagged.
    verdict = review(ok_validation(), library, None, preset="standard")
    assert verdict["can_run"]
    assert verdict["expected_status"] == "Development only"
    assert any("synthetic" in w.lower() for w in verdict["warnings"])


def test_running_without_a_library_is_refused(tmp_path):
    verdict = review(ok_validation(), None, None, preset="standard")
    assert not verdict["can_run"]


def test_a_file_with_no_usable_sequences_is_refused(tmp_path):
    library = inspect_library(write_corpus(tmp_path / "validated.csv"))
    verdict = review(ok_validation(n_valid=0), library, None, preset="standard")
    assert not verdict["can_run"]


def test_a_verified_library_with_no_target_is_publication_eligible(tmp_path):
    path = write_corpus(tmp_path / "validated.csv")
    path.with_suffix(".manifest.json").write_text(
        json.dumps({"source": "s", "curator": "c", "curated_date": "d"}), encoding="utf-8"
    )
    verdict = review(ok_validation(), inspect_library(path), None, preset="evaluation")
    assert verdict["can_run"]
    assert verdict["expected_status"] == "Publication-eligible"
    # Tier 1 only is legitimate, but the user must know they chose it.
    assert any("No target selected" in w for w in verdict["warnings"])


def test_excluded_rows_are_surfaced_as_a_warning(tmp_path):
    library = inspect_library(write_corpus(tmp_path / "validated.csv"))
    verdict = review(
        ok_validation(n_valid=90, n_rejected=10), library, None, preset="standard"
    )
    assert verdict["can_run"]
    assert any("10 of 100" in w for w in verdict["warnings"])
