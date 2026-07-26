from __future__ import annotations

import pytest

from aptarank.tier1 import shuffles

from .conftest import SPEC_SEQUENCE


def test_shuffles_preserve_composition_and_length():
    out = shuffles.generate_shuffles(SPEC_SEQUENCE, 25, k=2, seed=1)
    assert len(out) == 25
    for s in out:
        assert sorted(s) == sorted(SPEC_SEQUENCE)


def test_shuffles_are_seed_reproducible():
    a = shuffles.generate_shuffles(SPEC_SEQUENCE, 10, k=2, seed=99)
    b = shuffles.generate_shuffles(SPEC_SEQUENCE, 10, k=2, seed=99)
    c = shuffles.generate_shuffles(SPEC_SEQUENCE, 10, k=2, seed=100)
    assert a == b
    assert a != c


def test_repeated_shuffles_are_not_memory_corrupted():
    """Regression: ushuffle.Shuffler keeps a raw pointer to its input bytes.

    Letting that temporary be collected produced silently corrupted controls
    (non-ACGU bytes) from the second shuffle onwards.
    """
    seq = "UGAUGGCGCCCUGUUAAGCCACCACCUGAA"
    out = shuffles.generate_shuffles(seq, 50, k=2, seed=7)
    assert all(set(s) <= set("ACGU") for s in out)


def test_dinucleotide_counts_are_preserved():
    """The property §4.8 actually claims, checked rather than assumed."""
    seq = SPEC_SEQUENCE
    expected = shuffles.kmer_counts(seq, 2)
    for s in shuffles.generate_shuffles(seq, 5, k=2, seed=3):
        assert shuffles.kmer_counts(s, 2) == expected


def test_kmer_counts_counts_overlapping_windows():
    assert shuffles.kmer_counts("AAGA", 2) == {"AA": 1, "AG": 1, "GA": 1}


def test_monte_carlo_p_value_convention():
    """p = (1 + #{control >= real}) / (M + 1); beating all 20 gives 1/21."""
    controls = [0.1] * 20
    beat_all = shuffles.evaluate(0.9, controls, alpha=0.05)
    assert beat_all.wins == 20
    assert beat_all.percentile == pytest.approx(1.0)
    assert beat_all.p_value == pytest.approx(1 / 21)
    assert beat_all.passed

    # Beating 19 of 20 is a 95% win rate but NOT significant at alpha=0.05.
    beat_19 = shuffles.evaluate(0.9, [0.1] * 19 + [1.0], alpha=0.05)
    assert beat_19.percentile == pytest.approx(0.95)
    assert beat_19.p_value == pytest.approx(2 / 21)
    assert not beat_19.passed


def test_margin_is_relative_to_the_control_median():
    outcome = shuffles.evaluate(0.8, [0.2, 0.3, 0.4], alpha=0.05)
    assert outcome.median_score == pytest.approx(0.3)
    assert outcome.margin == pytest.approx(0.5)
