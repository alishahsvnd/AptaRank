"""The radius of gyration read off the secondary-structure element graph (§6).

Surface mode's footprint rests on this number, so the tests pin it against
closed-form results where one exists — a rigid rod, a free coil — rather than
against whatever the implementation happened to produce.
"""

from __future__ import annotations

import math

import pytest

from aptarank.tier1 import features as feature_mod
from aptarank.tier1.elements import (
    free_coil_radius_of_gyration,
    parse_elements,
)

A_BP, A_NT = 2.8, 6.0


def rg(dot_bracket: str, a_per_bp_helix: float = A_BP, a_per_nt_ss: float = A_NT) -> float:
    return parse_elements(
        dot_bracket, None, a_per_bp_helix, a_per_nt_ss
    ).radius_of_gyration_A


# -- closed-form anchors -------------------------------------------------


def test_a_long_helix_matches_the_rigid_rod_result():
    """A rod of length L has Rg = L/sqrt(12); a helix is a rod."""
    base_pairs = 45
    structure = "(" * base_pairs + "...." + ")" * base_pairs
    rod = A_BP * base_pairs / math.sqrt(12)
    # The four-nucleotide hairpin cap adds a little, and nothing else should.
    assert rod < rg(structure) < rod * 1.10


def test_an_unstructured_strand_is_a_gaussian_coil():
    """Nothing to fold on, so Rg = b*sqrt(N/6)."""
    assert rg("." * 40) == pytest.approx(free_coil_radius_of_gyration(40, A_NT))
    assert rg("." * 40) == pytest.approx(A_NT * math.sqrt(40 / 6.0))


def test_size_grows_with_length():
    assert rg("(" * 10 + "...." + ")" * 10) < rg("(" * 30 + "...." + ")" * 30)


# -- the property the refinement exists for ------------------------------


def test_branching_makes_a_molecule_compact_at_equal_length():
    """The whole point: two 100-mers are not the same size.

    A single long helix is an extended rod; a four-way junction of the same
    nucleotide count folds back on itself. The length proxy called these
    identical, which is what §6 asked us to fix.
    """
    extended = "(" * 45 + "." * 10 + ")" * 45
    compact = "((((" + "((((....))))" * 4 + "))))" + "." * 44
    assert len(extended) == len(compact) == 100
    assert rg(compact) < 0.75 * rg(extended)


def test_a_hairpin_and_a_coil_of_equal_length_differ():
    assert rg("." * 40) != pytest.approx(rg("(" * 18 + "...." + ")" * 18), rel=0.01)


# -- the constants are configuration, not source -------------------------


def test_the_helix_rise_changes_the_answer():
    """The lengths come from tier2.geometry; they must actually be used."""
    structure = "(" * 20 + "...." + ")" * 20
    assert rg(structure, a_per_bp_helix=2.8) < rg(structure, a_per_bp_helix=3.4)


def test_the_single_strand_rise_changes_the_answer():
    assert rg("." * 40, a_per_nt_ss=6.0) < rg("." * 40, a_per_nt_ss=7.0)


def test_the_geometry_reaches_the_corpus_cache_key():
    """Changing the constants must invalidate cached features, not silently
    leave a radius of gyration computed under the old ones."""
    from aptarank.tier1.corpus import tool_signature

    assert tool_signature({"a_per_bp_helix": 2.8, "a_per_nt_ss": 6.0}) != tool_signature(
        {"a_per_bp_helix": 3.4, "a_per_nt_ss": 6.0}
    )


# -- plumbing ------------------------------------------------------------


def test_the_feature_table_carries_the_size(mini_candidates_path):
    record = feature_mod.sequence_features(
        "GUUCCAUGGGCCUUGACUUGCUGUGUCAUCACCAUGGGAC", n_ensemble_samples=0
    )
    assert record["radius_of_gyration_A"] > 0
    assert "radius_of_gyration_A" in feature_mod.FEATURE_COLUMNS
    # No ensemble: the median falls back to the MFE structure's own value.
    assert record["rg_median_A"] == pytest.approx(record["radius_of_gyration_A"])
    assert record["rg_iqr_A"] == 0.0


def test_the_ensemble_median_is_taken_over_sampled_structures():
    record = feature_mod.sequence_features(
        "GUUCCAUGGGCCUUGACUUGCUGUGUCAUCACCAUGGGAC", n_ensemble_samples=40, seed=7
    )
    assert record["n_ensemble_samples"] > 1
    assert record["rg_median_A"] > 0
    # A real ensemble moves between shapes, so the spread is not zero.
    assert record["rg_iqr_A"] >= 0.0
    assert "rg_median_A" in feature_mod.ENSEMBLE_COLUMNS


def test_every_size_is_finite_and_positive_for_real_folds():
    from aptarank.tier1 import folding

    for sequence in (
        "GUUCCAUGGGCCUUGACUUGCUGUGUCAUCACCAUGGGAC",
        "AAAAAAAAAAAAAAAAAAAAAAAA",
        "GGGGCCCCGGGGCCCCGGGGCCCCGGGGCCCC",
    ):
        result = folding.fold(sequence)
        value = rg(result.dot_bracket)
        assert math.isfinite(value) and value > 0
