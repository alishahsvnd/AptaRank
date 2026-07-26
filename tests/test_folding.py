"""Regression tests against the verified spec §4.5–4.7 values.

If any of these change, either ViennaRNA changed or the wrapper broke its call
ordering. Both are things we must find out about immediately, because every
score in the project is downstream of these numbers.
"""

from __future__ import annotations

import pytest

from aptarank.errors import FoldingError
from aptarank.tier1 import elements, folding

from .conftest import (
    SPEC_ENSEMBLE_DEFECT,
    SPEC_MFE,
    SPEC_MFE_NORM,
    SPEC_SEQUENCE,
    SPEC_STRUCTURE,
)


def test_fold_reproduces_spec_fixture():
    result = folding.fold(SPEC_SEQUENCE)
    assert result.dot_bracket == SPEC_STRUCTURE
    assert result.mfe == pytest.approx(SPEC_MFE, abs=1e-3)
    assert result.mfe_norm == pytest.approx(SPEC_MFE_NORM, abs=1e-4)
    assert result.ensemble_defect == pytest.approx(SPEC_ENSEMBLE_DEFECT, abs=1e-4)


def test_positional_entropy_has_one_value_per_position():
    """ViennaRNA returns a length n+1 list whose index 0 is unused."""
    result = folding.fold(SPEC_SEQUENCE)
    assert len(result.positional_entropy) == len(SPEC_SEQUENCE)
    assert result.positional_entropy_mean == pytest.approx(
        sum(result.positional_entropy) / len(SPEC_SEQUENCE)
    )


def test_sampling_returns_structures_and_is_seed_reproducible():
    _r1, s1 = folding.fold_and_sample(SPEC_SEQUENCE, 50, seed=123)
    _r2, s2 = folding.fold_and_sample(SPEC_SEQUENCE, 50, seed=123)
    assert len(s1) >= 45  # ViennaRNA may drop a couple of redundant samples
    assert s1 == s2
    assert all(len(db) == len(SPEC_SEQUENCE) for db in s1)


def test_empty_sequence_is_rejected():
    with pytest.raises(FoldingError):
        folding.fold("")


def test_forgi_elements_match_spec_fixture():
    feats = elements.parse_elements(SPEC_STRUCTURE, SPEC_SEQUENCE)
    assert feats.n_hairpins == 1
    assert feats.n_interior == 1
    assert feats.n_multiloop == 0
    assert feats.n_stems == 2
    assert feats.longest_stem_bp == 9
    assert feats.max_loop_nt == 7
    assert feats.total_unpaired == 14
    assert feats.element_string == "sssssssssiiiisssshhhhhhhssssiiisssssssss"
    # 9 bp + 4 bp stems, both strands: (9 + 4) * 2 = 26 of 40 positions
    assert feats.stem_fraction == pytest.approx(26 / 40)


def test_unfolded_structure_does_not_raise():
    feats = elements.parse_elements("." * 30)
    assert feats.n_stems == 0
    assert feats.stem_fraction == 0.0
    assert feats.total_unpaired == 30
