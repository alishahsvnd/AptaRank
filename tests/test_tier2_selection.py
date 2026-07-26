from __future__ import annotations

import pytest

from aptarank.errors import TargetError
from aptarank.tier2 import selection
from aptarank.tier2.fpocket import Pocket, Residue


def pocket(index, score, residues, druggability=0.5, volume=100.0, spheres=10):
    return Pocket(
        index=index,
        metrics={
            "Score": score,
            "Druggability Score": druggability,
            "Volume": volume,
            "Number of Alpha Spheres": spheres,
        },
        lining_residues=[Residue("A", n, "") for n in residues],
    )


ALL_RESIDUES = [Residue("A", n, "") for n in range(1, 400)]


def test_active_site_overlap_beats_a_higher_fpocket_score():
    """The functional cavity is not always fpocket's top-scoring one."""
    pockets = [pocket(1, 0.90, [10, 11, 12]), pocket(2, 0.30, [120, 124, 208])]
    result = selection.select_pocket(
        pockets,
        requested=selection.parse_residue_specs([120, 124, 208], "A"),
        structure_residues=ALL_RESIDUES,
    )
    assert result["selected_pocket_index"] == 2
    assert result["method"] == "active_site_overlap"
    assert [e["overlap_count"] for e in result["pocket_evidence"]] == [0, 3]


def test_ties_break_on_fpocket_score_then_index():
    pockets = [pocket(1, 0.20, [120]), pocket(2, 0.80, [124]), pocket(3, 0.80, [208])]
    result = selection.select_pocket(
        pockets,
        requested=selection.parse_residue_specs([120, 124, 208], "A"),
        structure_residues=ALL_RESIDUES,
    )
    assert result["selected_pocket_index"] == 2  # tied overlap 1, higher score, lower index


def test_zero_overlap_fails_the_build_by_default():
    """Almost always a numbering or chain mismatch, not a real finding."""
    pockets = [pocket(1, 0.9, [10, 11]), pocket(2, 0.4, [12, 13])]
    with pytest.raises(TargetError, match="overlapped no pocket|line any"):
        selection.select_pocket(
            pockets,
            requested=selection.parse_residue_specs([120, 124], "A"),
            structure_residues=ALL_RESIDUES,
        )


def test_zero_overlap_fallback_is_explicit_and_loudly_labelled():
    pockets = [pocket(1, 0.9, [10, 11]), pocket(2, 0.4, [12, 13])]
    result = selection.select_pocket(
        pockets,
        requested=selection.parse_residue_specs([120], "A"),
        structure_residues=ALL_RESIDUES,
        allow_zero_overlap_fallback=True,
    )
    assert result["method"] == "active_site_zero_overlap_fallback"
    assert result["selected_pocket_index"] == 1
    assert any("must NOT be described as active-site" in w for w in result["warnings"])


def test_residue_absent_from_the_structure_is_caught_before_selection():
    pockets = [pocket(1, 0.9, [10])]
    with pytest.raises(TargetError, match="absent from the prepared structure"):
        selection.select_pocket(
            pockets,
            requested=selection.parse_residue_specs([9999], "A"),
            structure_residues=ALL_RESIDUES,
        )


def test_automatic_selection_is_recorded_with_a_caveat():
    pockets = [pocket(1, 0.30, [10]), pocket(2, 0.95, [20])]
    result = selection.select_pocket(pockets)
    assert result["method"] == "automatic_fpocket_score"
    assert result["selected_pocket_index"] == 2
    assert result["warnings"]


def test_residue_specs_accept_numbers_and_explicit_chains():
    specs = selection.parse_residue_specs(
        [120, "124", {"chain_id": "B", "residue_number": 208, "insertion_code": "A"}], "A"
    )
    assert [s.key() for s in specs] == [("A", 120, ""), ("A", 124, ""), ("B", 208, "A")]
