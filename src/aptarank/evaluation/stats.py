"""Statistics shared by the experiments.

Effect sizes with confidence intervals, not bare p-values: with a few thousand
sequences almost any difference is "significant", and a reviewer will rightly
ask how large it is.
"""

from __future__ import annotations

from typing import Any, Sequence

import numpy as np
from scipy.stats import mannwhitneyu


def probability_of_superiority(a: Sequence[float], b: Sequence[float]) -> float:
    """P(random draw from `a` > random draw from `b`), ties counted as half.

    Identical to the AUROC of `a` vs `b`, and the natural effect size to pair
    with Mann-Whitney U — it is on a scale a reader can interpret (0.5 = no
    discrimination) rather than a U statistic that scales with sample size.
    """
    x, y = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    if x.size == 0 or y.size == 0:
        return float("nan")
    u, _p = mannwhitneyu(x, y, alternative="two-sided")
    return float(u / (x.size * y.size))


def compare_groups(
    a: Sequence[float], b: Sequence[float], n_boot: int = 2000, seed: int = 0
) -> dict[str, Any]:
    """Mann-Whitney U with a bootstrap CI on the effect size."""
    x, y = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    x, y = x[np.isfinite(x)], y[np.isfinite(y)]
    if x.size < 2 or y.size < 2:
        return {"n_a": int(x.size), "n_b": int(y.size), "auc": None,
                "p_value": None, "reason": "fewer than two observations in a group"}

    u, p = mannwhitneyu(x, y, alternative="two-sided")
    auc = float(u / (x.size * y.size))

    rng = np.random.default_rng(seed)
    boot = np.empty(n_boot)
    for i in range(n_boot):
        boot[i] = probability_of_superiority(
            rng.choice(x, x.size, replace=True), rng.choice(y, y.size, replace=True)
        )
    lo, hi = np.percentile(boot[np.isfinite(boot)], [2.5, 97.5])

    return {
        "n_a": int(x.size),
        "n_b": int(y.size),
        "median_a": float(np.median(x)),
        "median_b": float(np.median(y)),
        "u_statistic": float(u),
        "p_value": float(p),
        "auc": auc,
        "auc_ci95": [float(lo), float(hi)],
        "reason": None,
    }


def compare_groups_clustered(
    a: Sequence[float],
    a_clusters: Sequence[Any],
    b: Sequence[float],
    n_boot: int = 2000,
    seed: int = 0,
) -> dict[str, Any]:
    """Mann-Whitney U with the positive group's CI bootstrapped by cluster.

    Aptamers selected against the same protein are often near-relatives, so
    resampling them individually treats correlated observations as independent
    and yields a confidence interval that is too narrow. Clusters (targets) are
    resampled instead. The control group genuinely is independent — random
    sequences and shuffles of distinct sequences — so it is resampled i.i.d.
    """
    x = np.asarray(a, dtype=float)
    y = np.asarray(b, dtype=float)
    clusters = np.asarray(a_clusters, dtype=object)
    mask = np.isfinite(x)
    x, clusters = x[mask], clusters[mask]
    y = y[np.isfinite(y)]

    base = compare_groups(x, y, n_boot=1, seed=seed)
    unique = list(dict.fromkeys(clusters.tolist()))
    if base["auc"] is None or len(unique) < 3:
        return {**base, "n_clusters_a": len(unique),
                "cluster_bootstrap": False,
                "reason": "too few clusters for a cluster bootstrap"}

    by_cluster = {c: x[clusters == c] for c in unique}
    rng = np.random.default_rng(seed)
    boot = np.empty(n_boot)
    for i in range(n_boot):
        picked = rng.choice(len(unique), len(unique), replace=True)
        resampled = np.concatenate([by_cluster[unique[j]] for j in picked])
        boot[i] = probability_of_superiority(
            resampled, rng.choice(y, y.size, replace=True)
        )
    boot = boot[np.isfinite(boot)]
    return {
        **base,
        "n_clusters_a": len(unique),
        "cluster_bootstrap": True,
        "auc_ci95": [float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))],
        "reason": None,
    }


def paired_bootstrap(
    differences: Sequence[float], n_boot: int = 2000, seed: int = 0
) -> dict[str, Any]:
    """Bootstrap CI on a mean paired difference, resampling the pairs."""
    d = np.asarray(differences, dtype=float)
    d = d[np.isfinite(d)]
    if d.size < 2:
        return {"n": int(d.size), "mean": None, "ci95": None}
    rng = np.random.default_rng(seed)
    means = np.array([rng.choice(d, d.size, replace=True).mean() for _ in range(n_boot)])
    return {
        "n": int(d.size),
        "mean": float(d.mean()),
        "median": float(np.median(d)),
        "ci95": [float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5))],
        "fraction_positive": float((d > 0).mean()),
    }


def bootstrap_auc_by_unit(
    positives_by_unit: Sequence[Sequence[float]],
    negatives_by_unit: Sequence[Sequence[float]],
    n_boot: int = 2000,
    seed: int = 0,
) -> dict[str, Any]:
    """AUROC with the bootstrap resampling *units*, not individual pairs.

    Each aptamer contributes one true-target score and several decoy scores.
    Resampling pairs would treat those decoys as independent observations and
    produce a confidence interval far too narrow.
    """
    units = [
        (np.asarray(p, dtype=float), np.asarray(n, dtype=float))
        for p, n in zip(positives_by_unit, negatives_by_unit)
        if len(p) and len(n)
    ]
    if len(units) < 2:
        return {"n_units": len(units), "auc": None, "ci95": None}

    def auc_of(selection) -> float:
        pos = np.concatenate([units[i][0] for i in selection])
        neg = np.concatenate([units[i][1] for i in selection])
        return probability_of_superiority(pos, neg)

    point = auc_of(range(len(units)))
    rng = np.random.default_rng(seed)
    boot = np.array([
        auc_of(rng.integers(0, len(units), len(units))) for _ in range(n_boot)
    ])
    boot = boot[np.isfinite(boot)]
    return {
        "n_units": len(units),
        "auc": float(point),
        "ci95": [float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))],
    }


def mean_reciprocal_rank(true_scores: Sequence[float], decoy_scores: Sequence[Sequence[float]]) -> dict[str, Any]:
    """How often the true target ranks first among its decoys.

    Reported alongside AUROC because a retrieval framing is what a biologist
    actually cares about: given this aptamer, does the tool point at the right
    protein?
    """
    reciprocal, top1 = [], []
    for true, decoys in zip(true_scores, decoy_scores):
        pool = np.asarray([true, *decoys], dtype=float)
        if not np.isfinite(pool).all() or pool.size < 2:
            continue
        # rank of the true score among all candidates, 1 = best
        rank = 1 + int((pool[1:] > pool[0]).sum()) + 0.5 * int((pool[1:] == pool[0]).sum())
        reciprocal.append(1.0 / rank)
        top1.append(rank == 1)
    if not reciprocal:
        return {"n": 0, "mrr": None, "top1_accuracy": None}
    return {
        "n": len(reciprocal),
        "mrr": float(np.mean(reciprocal)),
        "top1_accuracy": float(np.mean(top1)),
    }
