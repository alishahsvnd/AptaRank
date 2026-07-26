from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from aptarank.tier1 import scoring
from aptarank.tier1.scoring import ReferenceDistributions


def test_empirical_cdf_uses_midranks_for_ties():
    ref = np.array([1.0, 2.0, 2.0, 3.0])
    # below everything / at a tied value / above everything
    got = scoring.empirical_cdf(ref, np.array([0.0, 2.0, 9.0]))
    assert got[0] == 0.0
    assert got[1] == pytest.approx((1 + 0.5 * 2) / 4)
    assert got[2] == 1.0


def test_two_sided_score_peaks_at_the_median():
    assert scoring.two_sided_score(np.array([0.5]))[0] == pytest.approx(1.0)
    assert scoring.two_sided_score(np.array([0.0]))[0] == pytest.approx(0.0)
    assert scoring.two_sided_score(np.array([1.0]))[0] == pytest.approx(0.0)
    assert scoring.two_sided_score(np.array([0.25]))[0] == pytest.approx(0.5)


def test_one_sided_score_flips_for_lower_is_better():
    p = np.array([0.2, 0.8])
    assert list(scoring.one_sided_score(p, higher_is_better=True)) == [0.2, 0.8]
    assert list(scoring.one_sided_score(p, higher_is_better=False)) == pytest.approx([0.8, 0.2])


def _refs() -> ReferenceDistributions:
    rng = np.random.default_rng(0)
    return ReferenceDistributions(
        values={
            "mfe_norm": rng.normal(-0.4, 0.1, 500),
            "ensemble_defect": rng.uniform(0, 0.3, 500),
            "positional_entropy_mean": rng.uniform(0, 0.5, 500),
            "stem_fraction": rng.uniform(0.2, 0.8, 500),
            "gc_fraction": rng.uniform(0.3, 0.7, 500),
        },
        n_sequences=500,
        corpus_id="test",
    )


def test_corpus_scores_do_not_depend_on_the_batch():
    """The property that makes runs comparable — and evaluation E1 valid."""
    refs = _refs()
    criteria = list(refs.values)
    one = pd.DataFrame([{"mfe_norm": -0.5, "ensemble_defect": 0.05,
                         "positional_entropy_mean": 0.1, "stem_fraction": 0.5,
                         "gc_fraction": 0.5}])
    many = pd.concat([one, pd.DataFrame([{"mfe_norm": -0.9, "ensemble_defect": 0.01,
                                          "positional_entropy_mean": 0.02,
                                          "stem_fraction": 0.5, "gc_fraction": 0.5}])],
                     ignore_index=True)
    weights = {c: 1.0 for c in criteria}

    alone = scoring.composite(
        scoring.criterion_scores(one, refs, criteria), weights, "corpus_weighted_mean"
    ).iloc[0]
    with_company = scoring.composite(
        scoring.criterion_scores(many, refs, criteria), weights, "corpus_weighted_mean"
    ).iloc[0]
    assert alone == pytest.approx(with_company)


def test_batch_rank_aggregation_does_depend_on_the_batch():
    """Kept available, but this is why it is not the default."""
    refs = _refs()
    criteria = list(refs.values)
    weights = {c: 1.0 for c in criteria}
    frame = pd.DataFrame(
        [
            {"mfe_norm": -0.5, "ensemble_defect": 0.05, "positional_entropy_mean": 0.1,
             "stem_fraction": 0.5, "gc_fraction": 0.5},
            {"mfe_norm": -0.9, "ensemble_defect": 0.01, "positional_entropy_mean": 0.02,
             "stem_fraction": 0.5, "gc_fraction": 0.5},
        ]
    )
    scores = scoring.criterion_scores(frame, refs, criteria)
    single = scoring.composite(scores.head(1), weights, "batch_rank_aggregation").iloc[0]
    assert single == pytest.approx(1.0)  # a batch of one scores perfectly


def test_non_finite_reference_values_are_rejected():
    from aptarank.errors import CorpusError

    with pytest.raises(CorpusError, match="non-finite"):
        scoring.empirical_cdf(np.array([1.0, np.nan, 3.0]), np.array([2.0]))


def test_non_finite_candidate_values_are_rejected():
    """NaN would otherwise sort to the top and score as a perfect candidate."""
    from aptarank.errors import CorpusError

    with pytest.raises(CorpusError, match="non-finite"):
        scoring.empirical_cdf(np.array([1.0, 2.0, 3.0]), np.array([np.nan]))


def test_missing_weight_is_an_error_not_an_implicit_zero():
    from aptarank.errors import CorpusError

    scores = pd.DataFrame({"a": [0.5], "b": [0.5]})
    with pytest.raises(CorpusError, match="no weight declared"):
        scoring.weighted_mean(scores, {"a": 1.0})


def test_structural_subscore_excludes_composition_criteria():
    refs = _refs()
    frame = pd.DataFrame([{"mfe_norm": -0.5, "ensemble_defect": 0.05,
                           "positional_entropy_mean": 0.1, "stem_fraction": 0.5,
                           "gc_fraction": 0.99}])
    structural = ["mfe_norm", "ensemble_defect", "positional_entropy_mean", "stem_fraction"]
    scores = scoring.criterion_scores(frame, refs, list(refs.values))
    sub = scoring.structural_subscore(scores, structural, {c: 1.0 for c in structural})
    expected = scores[structural].mean(axis=1)
    assert sub.iloc[0] == pytest.approx(expected.iloc[0])
