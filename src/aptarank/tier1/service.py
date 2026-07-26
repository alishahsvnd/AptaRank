"""Tier 1 orchestration: features -> corpus scores -> controls -> ranking.

This module owns the answer to "what order do the candidates come in?".
Tier 2 may annotate that order but never changes it (spec §6).
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Any, Callable

import numpy as np
import pandas as pd

from ..config import CRITERIA, Config
from ..ingest import IngestResult
from . import features as feature_mod
from . import scoring, shuffles
from .corpus import CorpusInfo
from .scoring import ReferenceDistributions


@dataclass
class Tier1Result:
    """Ranked candidates plus everything needed to explain the ranking."""

    table: pd.DataFrame               # one row per candidate, ranked
    criterion_scores: pd.DataFrame    # candidate_id-indexed, one column per criterion
    failures: list[dict[str, Any]] = field(default_factory=list)
    runtime_seconds: float = 0.0
    n_shuffles_used: int = 0
    composite_method: str = "corpus_weighted_mean"

    @property
    def n_scored(self) -> int:
        return len(self.table)


def run(
    cfg: Config,
    ingested: IngestResult,
    refs: ReferenceDistributions,
    corpus_info: CorpusInfo | None = None,
    progress: Callable[[int, int], None] | None = None,
) -> Tier1Result:
    started = time.perf_counter()

    criteria = cfg.active_criteria()
    weights = cfg.get("tier1.weights")
    seed = int(cfg.get("run.seed"))
    n_samples = 0 if cfg.is_fast else int(cfg.get("tier1.n_ensemble_samples"))
    n_shuffles = int(cfg.get("tier1.shuffle.n_shuffles")) if cfg.shuffles_enabled else 0

    ids = ingested.candidates["candidate_id"]
    if ids.duplicated().any():
        duplicated = sorted(ids[ids.duplicated()].unique())[:5]
        raise ValueError(
            f"candidate ids must be unique before scoring; duplicates: {duplicated}"
        )

    jobs = [
        feature_mod.FeatureJob(
            candidate_id=row.candidate_id,
            sequence=row.sequence,
            n_ensemble_samples=n_samples,
            n_shuffles=n_shuffles,
            shuffle_k=int(cfg.get("tier1.shuffle.k")),
            seed=seed,
        )
        for row in ingested.candidates.itertuples()
    ]

    results = feature_mod.compute_batch(
        jobs,
        workers=cfg.get("tier1.parallel.workers", None),
        chunk_size=int(cfg.get("tier1.parallel.chunk_size", 16)),
        progress=progress,
    )
    table, failures = feature_mod.results_to_frame(results)
    table = table.merge(
        ingested.candidates[["candidate_id", "duplicate_count"]],
        on="candidate_id",
        how="left",
    )

    # -- corpus-calibrated criterion scores (absolute, batch-independent)
    scores = scoring.criterion_scores(table, refs, criteria)
    scores.index = table["candidate_id"]

    method = cfg.get("tier1.composite.method")
    table["tier1_score"] = scoring.composite(scores, weights, method).to_numpy()
    table["batch_rank_fraction"] = scoring.batch_rank_fraction(table["tier1_score"]).to_numpy()

    for name in criteria:
        table[f"score__{name}"] = scores[name].to_numpy()

    # -- shuffled controls (§4.8)
    shuffle_rows = _shuffle_outcomes(cfg, results, refs, weights)
    if shuffle_rows is not None:
        table = table.merge(shuffle_rows, on="candidate_id", how="left")
    else:
        for column, default in (
            ("shuffle_pass", pd.NA), ("shuffle_percentile", np.nan),
            ("shuffle_p_value", np.nan), ("shuffle_margin", np.nan),
            ("shuffle_n", 0), ("shuffle_median_score", np.nan),
            ("structural_subscore", np.nan),
        ):
            table[column] = default

    # -- the ranking. `rank` is a dense ordinal broken deterministically by
    #    candidate_id so reruns are identical; `rank_min` is the statistical
    #    (competition) rank, which is the honest one to quote when scores tie.
    table = table.sort_values(
        ["tier1_score", "candidate_id"], ascending=[False, True]
    ).reset_index(drop=True)
    table.insert(0, "rank", np.arange(1, len(table) + 1))
    table["rank_min"] = (
        table["tier1_score"].rank(method="min", ascending=False).astype(int)
    )
    table["rank_is_tied"] = table["tier1_score"].duplicated(keep=False)

    scores = scores.loc[table["candidate_id"]]

    return Tier1Result(
        table=table,
        criterion_scores=scores,
        failures=failures,
        runtime_seconds=time.perf_counter() - started,
        n_shuffles_used=n_shuffles,
        composite_method=method,
    )


def _shuffle_outcomes(
    cfg: Config,
    results: list[dict[str, Any]],
    refs: ReferenceDistributions,
    weights: dict[str, float],
) -> pd.DataFrame | None:
    """Score every shuffled control and compare it to its parent candidate.

    All controls across all candidates are scored in one vectorised pass, then
    grouped back per candidate — folding happened in the workers, scoring
    happens here where the corpus lives.
    """
    if not cfg.shuffles_enabled:
        return None

    structural = [
        c for c in cfg.get("tier1.shuffle.structural_criteria") if c in CRITERIA
    ]
    alpha = float(cfg.get("tier1.shuffle.alpha"))

    control_rows: list[dict[str, Any]] = []
    real_rows: list[dict[str, Any]] = []
    for res in results:
        if res["error"] is not None or res["features"] is None:
            continue
        real_rows.append({"candidate_id": res["candidate_id"], **res["features"]})
        for control in res["shuffles"]:
            control_rows.append({"candidate_id": res["candidate_id"], **control})

    if not control_rows:
        return None

    real = pd.DataFrame(real_rows)
    controls = pd.DataFrame(control_rows)

    real_sub = scoring.structural_subscore(
        scoring.criterion_scores(real, refs, structural), structural, weights
    )
    control_sub = scoring.structural_subscore(
        scoring.criterion_scores(controls, refs, structural), structural, weights
    )
    controls = controls.assign(structural_subscore=control_sub.to_numpy())
    grouped = controls.groupby("candidate_id")
    scores_by_candidate = grouped["structural_subscore"].apply(list)
    sequences_by_candidate = grouped["sequence"].apply(list)

    expected = int(cfg.get("tier1.shuffle.n_shuffles"))
    rows = []
    for candidate_id, real_score, real_seq in zip(
        real["candidate_id"], real_sub, real["sequence"]
    ):
        control_scores = scores_by_candidate.get(candidate_id, [])
        control_seqs = sequences_by_candidate.get(candidate_id, [])
        outcome = shuffles.evaluate(float(real_score), control_scores, alpha)
        rows.append(
            {
                "candidate_id": candidate_id,
                "structural_subscore": float(real_score),
                "shuffle_pass": outcome.passed,
                "shuffle_percentile": outcome.percentile,
                "shuffle_p_value": outcome.p_value,
                "shuffle_margin": outcome.margin,
                "shuffle_median_score": outcome.median_score,
                "shuffle_n": outcome.n_shuffles,
                # A control that failed to fold would shrink the denominator
                # and quietly weaken the test; say so instead.
                "shuffle_complete": outcome.n_shuffles == expected,
                "shuffle_n_unique": len(set(control_seqs)),
                "shuffle_n_identical_to_real": sum(1 for s in control_seqs if s == real_seq),
            }
        )
    return pd.DataFrame(rows)
