"""Tier 2 orchestration: annotate the Tier 1 ranking, never reorder it.

Reads an immutable target bundle, projects the fixed calibration bank onto that
target, and assigns each survivor a control-relative band. The only thing that
scales with candidate count is arithmetic on two numbers.
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
from . import calibration


def run(
    cfg: Config,
    tier1: Tier1Result,
    corpus_table: pd.DataFrame,
    corpus_info: CorpusInfo,
    progress: Callable[[int, int], None] | None = None,
) -> dict[str, Any]:
    """Score the Tier 1 survivors against one target."""
    started = time.perf_counter()

    bundle = _load_bundle(cfg)
    pocket = bundle_mod.selected_pocket(bundle)
    d_pocket = float(pocket["geometry"]["d_pocket_A"])

    bank = calibration.build_or_load(cfg, corpus_table, corpus_info, progress=progress)

    a_per_nt = float(cfg.get("tier2.a_per_nt"))
    flex_c = float(cfg.get("tier2.flex_c"))
    sigma = float(cfg.get("tier2.sigma_A"))
    bands = cfg.get("tier2.band_percentiles")
    moderate, strong = float(bands["moderate"]), float(bands["strong"])

    dist_flex = calibration.target_distribution(
        bank, d_pocket, a_per_nt, flex_c, sigma, descriptor="flexible"
    )
    dist_ext = calibration.target_distribution(
        bank, d_pocket, a_per_nt, flex_c, sigma, descriptor="extended"
    )

    n_survivors = int(cfg.get("tier2.n_candidates"))
    survivors = tier1.table.head(n_survivors)

    records = calibration.score_candidates(
        survivors["loop_nt_median"].to_numpy(dtype=float),
        d_pocket_A=d_pocket,
        dist_flexible=dist_flex,
        dist_extended=dist_ext,
        a_per_nt=a_per_nt,
        flex_c=flex_c,
        sigma_A=sigma,
        moderate=moderate,
        strong=strong,
    )
    for record in records:
        record["target_bundle_id"] = bundle["bundle_id"]
        record["calibration_bank_id"] = bank.bank_id
        record["parameters"] = {"a_per_nt": a_per_nt, "flex_c": flex_c, "sigma_A": sigma}

    per_candidate = dict(zip(survivors["candidate_id"], records))

    return {
        "target": {
            **bundle_mod.summary(bundle),
            "bundle_path": str(_bundle_path(cfg)),
        },
        "thresholds": {
            **calibration.thresholds(dist_flex, moderate, strong),
            "calibration_bank": bank.meta,
        },
        "candidates": per_candidate,
        "spearman": tier_independence(survivors, records),
        "n_evaluated": sum(1 for r in records if r["status"] == "evaluated"),
        "n_survivors": len(survivors),
        "runtime_seconds": round(time.perf_counter() - started, 3),
    }


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
        [r.get("control_percentile_flexible", np.nan) for r in records], dtype=float
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


def _bundle_path(cfg: Config) -> Path:
    explicit = cfg.get("tier2.bundle_path", None)
    if explicit:
        return Path(explicit)
    pdb_id = cfg.get("tier2.target.pdb_id", None)
    if not pdb_id:
        raise TargetError(
            "Tier 2 is enabled but no target is configured. Set "
            "tier2.target.pdb_id, or point tier2.bundle_path at a bundle file."
        )
    return bundle_mod.find(cfg.get("tier2.bundle_dir"), pdb_id, cfg.get("tier2.target.chain", None))


def _load_bundle(cfg: Config) -> dict[str, Any]:
    return bundle_mod.load(_bundle_path(cfg))
