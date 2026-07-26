"""Corpus-calibrated criterion scores and the Tier 1 composite (spec §4.4, §4.9).

Two different questions, deliberately kept apart:

* **"Is this candidate normal?"** — answered against the *corpus* of validated
  aptamers, as a percentile. This is what every per-criterion bar in the UI
  shows, and it is an absolute quantity: it does not depend on what else was
  submitted alongside the candidate.
* **"Is this the best of what you gave me?"** — answered against the
  *submitted batch*, as the final ordinal rank.

Deviation from the spec's §4.9 code snippet, deliberately: the snippet computes
the composite from `rankdata` *within the batch*, which makes `tier1_score`
depend on unrelated co-submissions (a batch of one scores 1.0 on everything)
and makes evaluation E1 meaningless unless all three comparison groups are
scored in one pooled batch. The default here composes the corpus-calibrated
criterion scores directly, and the batch determines only `rank`. The literal
spec behaviour is still available as
`tier1.composite.method: batch_rank_aggregation`.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from ..config import CRITERIA
from ..errors import CorpusError


def empirical_cdf(reference: np.ndarray, values: np.ndarray) -> np.ndarray:
    """Tie-aware (mid-rank) empirical CDF of `values` within `reference`.

        F(x) = (#{r < x} + 0.5 * #{r == x}) / m

    The mid-rank convention matters: without it, a criterion whose corpus
    values are heavily tied (small integer counts) would push every tied
    candidate to one end of the distribution.
    """
    ref = np.asarray(reference, dtype=float)
    x = np.asarray(values, dtype=float)
    m = ref.size
    if m == 0:
        raise CorpusError("cannot compute a percentile against an empty reference")
    # NaN sorts to the end of the array and would silently inflate every
    # percentile; inf breaks the comparison outright. Neither may pass.
    if not np.isfinite(ref).all():
        raise CorpusError(
            f"reference distribution contains {int((~np.isfinite(ref)).sum())} "
            f"non-finite value(s); percentiles would be meaningless"
        )
    if not np.isfinite(x).all():
        raise CorpusError(
            f"cannot score {int((~np.isfinite(x)).sum())} non-finite feature "
            f"value(s); the candidate should have been failed before scoring"
        )
    ref = np.sort(ref)
    left = np.searchsorted(ref, x, side="left")
    right = np.searchsorted(ref, x, side="right")
    return (left + 0.5 * (right - left)) / m


def one_sided_score(p: np.ndarray, higher_is_better: bool) -> np.ndarray:
    """Directional criterion: percentile, flipped when lower values are better."""
    return p if higher_is_better else 1.0 - p


def two_sided_score(p: np.ndarray) -> np.ndarray:
    """Typicality: 1.0 at the corpus median, falling to 0 at either tail.

    Used where "higher is better" is simply wrong. Both extremes of GC content
    and stem fraction are bad — an all-stem candidate has nothing available to
    contact a target, an all-loop one has no stable shape.
    """
    return np.clip(1.0 - 2.0 * np.abs(p - 0.5), 0.0, 1.0)


@dataclass
class ReferenceDistributions:
    """Per-criterion reference values drawn from the validated-aptamer corpus."""

    values: dict[str, np.ndarray]
    n_sequences: int
    corpus_id: str
    is_placeholder: bool = False

    def percentile(self, criterion: str, x: Sequence[float] | np.ndarray) -> np.ndarray:
        if criterion not in self.values:
            raise CorpusError(
                f"criterion {criterion!r} is missing from the reference corpus; "
                f"available: {sorted(self.values)}"
            )
        return empirical_cdf(self.values[criterion], np.asarray(x, dtype=float))

    def score(self, criterion: str, x: Sequence[float] | np.ndarray) -> np.ndarray:
        """Corpus-calibrated score in [0,1], higher always better."""
        mode = CRITERIA[criterion]["mode"]
        p = self.percentile(criterion, x)
        if mode == "two_sided":
            return two_sided_score(p)
        return one_sided_score(p, higher_is_better=(mode == "higher"))


def criterion_scores(
    features: pd.DataFrame,
    refs: ReferenceDistributions,
    criteria: Sequence[str],
) -> pd.DataFrame:
    """A [0,1] corpus-calibrated score per criterion, one row per candidate."""
    missing = [c for c in criteria if c not in features.columns]
    if missing:
        raise CorpusError(f"feature table is missing criterion columns: {missing}")
    return pd.DataFrame(
        {c: refs.score(c, features[c].to_numpy()) for c in criteria},
        index=features.index,
    )


def weighted_mean(scores: pd.DataFrame, weights: Mapping[str, float]) -> pd.Series:
    """Σ w·s / Σ w over the given criteria.

    A *missing* weight is an error, not an implicit zero: silently dropping a
    criterion because of a typo would change what the composite means while
    still producing a plausible number.
    """
    undeclared = [c for c in scores.columns if c not in weights]
    if undeclared:
        raise CorpusError(f"no weight declared for criteria: {sorted(undeclared)}")
    used = [c for c in scores.columns if float(weights[c]) > 0.0]
    if not used:
        raise CorpusError("no criterion has a positive weight")
    w = np.array([float(weights[c]) for c in used])
    return pd.Series(
        (scores[used].to_numpy() * w).sum(axis=1) / w.sum(),
        index=scores.index,
        name="score",
    )


def rank_aggregate(scores: pd.DataFrame, weights: Mapping[str, float]) -> pd.Series:
    """Literal spec §4.9: mean of within-batch normalised ranks per criterion."""
    used = [c for c in scores.columns if float(weights.get(c, 0.0)) > 0.0]
    n = len(scores)
    ranked = pd.DataFrame(
        {c: rankdata(scores[c].to_numpy()) / n for c in used}, index=scores.index
    )
    return weighted_mean(ranked, weights)


def composite(
    scores: pd.DataFrame, weights: Mapping[str, float], method: str
) -> pd.Series:
    if method == "corpus_weighted_mean":
        return weighted_mean(scores, weights)
    if method == "batch_rank_aggregation":
        return rank_aggregate(scores, weights)
    raise CorpusError(f"unknown composite method: {method!r}")


def structural_subscore(
    scores: pd.DataFrame,
    structural_criteria: Sequence[str],
    weights: Mapping[str, float],
) -> pd.Series:
    """Absolute sub-score used for the shuffled-control comparison (§4.8).

    Must be absolute rather than batch-ranked: it compares one candidate
    against its own 20 shuffles, not against the other candidates. GC content
    and length are excluded upstream because dinucleotide shuffling preserves
    them exactly, so they would contribute an identical constant to the
    candidate and every one of its controls.
    """
    missing = [c for c in structural_criteria if c not in scores.columns]
    if missing:
        raise CorpusError(
            f"structural criteria {missing} were requested for the shuffled-control "
            f"comparison but were not scored; available: {list(scores.columns)}"
        )
    subset = list(structural_criteria)
    sub_weights = {c: float(weights[c]) for c in subset if c in weights}
    undeclared = [c for c in subset if c not in sub_weights]
    if undeclared:
        raise CorpusError(f"no weight declared for structural criteria: {undeclared}")
    return weighted_mean(scores[subset], sub_weights)


def batch_rank_fraction(values: pd.Series) -> pd.Series:
    """Where each candidate sits within the submitted batch (display only).

    Named a fraction rather than a percentile to keep it distinct from the
    corpus percentiles that carry the scientific meaning.
    """
    n = len(values)
    return pd.Series(rankdata(values.to_numpy()) / n, index=values.index)
