"""Tier 2 orchestration: annotate the Tier 1 ranking, never reorder it.

Reads a target bundle (building it server-side first if the target is configured
by identifier rather than by file), projects the fixed calibration bank onto
that target *in the configured binding mode*, and assigns each survivor a
control-relative band. The only thing that scales with candidate count is
arithmetic on two numbers.
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any, Callable

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from ..config import Config
from ..errors import TargetError
from ..tier1.corpus import CorpusInfo
from ..tier1.service import Tier1Result
from . import bundle as bundle_mod
from . import calibration, modes


def run(
    cfg: Config,
    tier1: Tier1Result,
    corpus_table: pd.DataFrame,
    corpus_info: CorpusInfo,
    progress: Callable[[int, int], None] | None = None,
    bundle: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Score the Tier 1 survivors against one target."""
    started = time.perf_counter()

    bundle = bundle or _load_bundle(cfg)
    mode = modes.check_mode(_binding_mode(cfg, bundle))
    params = modes.parameters(cfg, mode)
    target = modes.target_measurement(bundle, mode)

    bank = calibration.build_or_load(cfg, corpus_table, corpus_info, progress=progress)

    bands = cfg.get("tier2.band_percentiles")
    moderate, strong = float(bands["moderate"]), float(bands["strong"])

    distribution = calibration.target_distribution(bank, mode, target, params)
    secondary = calibration.secondary_distributions(bank, mode, target, params)

    n_survivors = int(cfg.get("tier2.n_candidates"))
    survivors = tier1.table.head(n_survivors)

    column = modes.descriptor_column(mode, params)
    if column not in survivors.columns:
        raise TargetError(
            f"{mode} mode compares {column!r}, which is not in the Tier 1 table"
        )

    records = calibration.score_candidates(
        survivors[column].to_numpy(dtype=float),
        mode=mode,
        target=target,
        params=params,
        distribution=distribution,
        moderate=moderate,
        strong=strong,
        secondary=secondary,
    )
    for record in records:
        record["target_bundle_id"] = bundle["bundle_id"]
        record["calibration_bank_id"] = bank.bank_id
        record["parameters"] = params

    per_candidate = dict(zip(survivors["candidate_id"], records))

    return {
        "binding_mode": mode,
        "target": {
            **bundle_mod.summary(bundle),
            "bundle_path": str(_bundle_path(cfg, required=False) or ""),
            "measurement": target,
        },
        "thresholds": {
            **calibration.thresholds(distribution, moderate, strong),
            "calibration_bank": bank.meta,
            "parameters": params,
        },
        "candidates": per_candidate,
        "spearman": tier_independence(survivors, records),
        "n_evaluated": sum(1 for r in records if r["status"] == "evaluated"),
        "n_survivors": len(survivors),
        "runtime_seconds": round(time.perf_counter() - started, 3),
    }


def _binding_mode(cfg: Config, bundle: dict[str, Any]) -> str:
    """Which mode to score in, and refuse to guess when the two disagree.

    A bundle built for a surface patch carries no cavity measurement, and a
    pocket bundle carries no patch. Scoring one as the other would either fail
    obscurely later or, worse, compare against whatever number happened to be
    present.
    """
    requested = cfg.get("tier2.binding_mode")
    built_for = bundle.get("binding_mode")
    if built_for and built_for != requested:
        raise TargetError(
            f"this target was prepared for {built_for!r} mode but the run asks "
            f"for {requested!r}. Rebuild the target in {requested!r} mode, or "
            f"score it in {built_for!r}."
        )
    return requested


def tier_independence(survivors: pd.DataFrame, records: list[dict[str, Any]]) -> dict[str, Any]:
    """Spearman between Tier 1 score and the Tier 2 control percentile (§6.4).

    Correlated on `tier1_score`, not on `rank`: rank is an ordinal broken by
    candidate id, and correlating a rank (small = good) against a score
    (large = good) also flips the sign of agreement, which invites misreading.

    A high absolute correlation is worth investigating but does not prove Tier 2
    merely restates Tier 1 — both consume the same predicted structures. A low
    one does not prove independence either: the top-N cut restricts the range.
    """
    scores = survivors["tier1_score"].to_numpy(dtype=float)
    percentiles = np.array(
        [r.get("control_percentile", np.nan) for r in records], dtype=float
    )
    mask = np.isfinite(scores) & np.isfinite(percentiles)
    n = int(mask.sum())

    if n < 3 or np.unique(percentiles[mask]).size < 2 or np.unique(scores[mask]).size < 2:
        return {
            "rho": None, "p_value": None, "n": n,
            "reason": "insufficient variation or fewer than 3 evaluated survivors",
        }
    result = spearmanr(scores[mask], percentiles[mask])
    return {
        "rho": float(result.statistic),
        "p_value": float(result.pvalue),
        "n": n,
        "reason": None,
        "interpretation": (
            "Tier 2 largely restates Tier 1 — investigate before building on it"
            if abs(float(result.statistic)) > 0.7
            else "tiers carry largely independent information"
        ),
    }


def _bundle_path(cfg: Config, required: bool = True) -> Path | None:
    explicit = cfg.get("tier2.bundle_path", None)
    if explicit:
        return Path(explicit)
    identifier = cfg.get("tier2.target.id", None)
    if not identifier:
        if not required:
            return None
        raise TargetError(
            "Tier 2 is enabled but no target is configured. Set tier2.target.id "
            "(with tier2.target.source), or point tier2.bundle_path at a "
            "prepared target file."
        )
    try:
        return bundle_mod.find(
            cfg.get("tier2.bundle_dir"),
            identifier,
            cfg.get("tier2.target.chain", None),
            mode=cfg.get("tier2.binding_mode", None),
        )
    except TargetError:
        if required:
            raise
        return None


def _load_bundle(cfg: Config) -> dict[str, Any]:
    return bundle_mod.load(_bundle_path(cfg))
