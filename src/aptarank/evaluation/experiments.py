"""Experiments E1–E5 (spec §8.2).

Each returns a JSON-serialisable result dict. `figures.py` draws from those
dicts, so every figure in the paper is regenerable from stored results without
re-running any folding.
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence

import numpy as np
import pandas as pd

from ..config import Config
from ..errors import CorpusError
from ..tier1 import features as feature_mod
from ..tier1 import scoring, shuffles
from ..tier1.corpus import CorpusInfo
from ..tier1.scoring import ReferenceDistributions
from ..tier2 import bundle as bundle_mod
from ..tier2 import calibration
from . import groups as groups_mod
from . import stats
from .groups import ComparisonGroups, target_folds


# -- shared scoring helper ----------------------------------------------


def score_sequences(
    cfg: Config,
    frame: pd.DataFrame,
    refs: ReferenceDistributions,
    n_shuffles: int = 0,
    progress: Callable[[int, int], None] | None = None,
) -> pd.DataFrame:
    """Fold and score a set of sequences against given reference distributions."""
    criteria = cfg.active_criteria()
    weights = cfg.get("tier1.weights")
    seed = int(cfg.get("run.seed"))

    jobs = [
        feature_mod.FeatureJob(
            candidate_id=row.candidate_id,
            sequence=row.sequence,
            n_ensemble_samples=int(cfg.get("tier1.n_ensemble_samples")),
            n_shuffles=n_shuffles,
            shuffle_k=int(cfg.get("tier1.shuffle.k")),
            seed=seed,
        )
        for row in frame.itertuples()
    ]
    results = feature_mod.compute_batch(
        jobs,
        workers=cfg.get("tier1.parallel.workers", None),
        chunk_size=int(cfg.get("tier1.parallel.chunk_size", 16)),
        progress=progress,
    )
    table, _failures = feature_mod.results_to_frame(results)
    if table.empty:
        raise CorpusError("every sequence in this group failed to fold")

    scores = scoring.criterion_scores(table, refs, criteria)
    table["tier1_score"] = scoring.composite(scores, weights, "corpus_weighted_mean").to_numpy()

    if n_shuffles:
        structural = [c for c in cfg.get("tier1.shuffle.structural_criteria")]
        alpha = float(cfg.get("tier1.shuffle.alpha"))
        real_sub = scoring.structural_subscore(
            scoring.criterion_scores(table, refs, structural), structural, weights
        )
        passes, p_values = [], []
        for i, res in enumerate(r for r in results if r["error"] is None):
            controls = pd.DataFrame(res["shuffles"])
            control_sub = scoring.structural_subscore(
                scoring.criterion_scores(controls, refs, structural), structural, weights
            )
            outcome = shuffles.evaluate(float(real_sub.iloc[i]), control_sub.tolist(), alpha)
            passes.append(outcome.passed)
            p_values.append(outcome.p_value)
        table["shuffle_pass"] = passes
        table["shuffle_p_value"] = p_values

    return table.merge(frame[["candidate_id"]], on="candidate_id", how="right")


def _refs_from(table: pd.DataFrame, criteria: Sequence[str], corpus_id: str,
               is_placeholder: bool) -> ReferenceDistributions:
    return ReferenceDistributions(
        values={c: table[c].to_numpy(dtype=float) for c in criteria},
        n_sequences=len(table),
        corpus_id=corpus_id,
        is_placeholder=is_placeholder,
    )


# -- E1 ------------------------------------------------------------------


def e1_discrimination(
    cfg: Config,
    corpus_features: pd.DataFrame,
    corpus_info: CorpusInfo,
    groups: ComparisonGroups,
    n_folds: int = 5,
    progress: Callable[[str, int, int], None] | None = None,
) -> dict[str, Any]:
    """Does Tier 1 separate real aptamers from controls? (spec §8.2 E1)

    Deviation from the spec, deliberately: E1 as written is circular. The
    validated aptamers *define* the corpus percentiles and are then used as the
    positive group, so they are being scored against a distribution they
    themselves produced. Here each validated sequence is scored **out of fold**
    against a reference that excludes its own target's family, and the controls
    for that fold are scored against the same reduced reference.

    The claim this supports is "Tier 1 separates aptamers it has never seen from
    composition-matched controls", which is the claim worth making.
    """
    started = time.perf_counter()
    criteria = cfg.active_criteria()
    folds = target_folds(corpus_features, n_folds, int(cfg.get("run.seed")))

    per_group: dict[str, list[float]] = {"validated": [], "random": [], "shuffled": []}
    # Cluster labels for the positive group, so its CI can be bootstrapped by
    # target rather than by sequence.
    validated_clusters: list[str] = []
    target_labels, grouping_key = groups_mod.canonical_target(corpus_features)
    fold_records = []

    for fold in folds:
        test_index = fold["test_index"]
        train = corpus_features.drop(corpus_features.index[test_index])
        if len(train) < 50:
            continue
        refs = _refs_from(train, criteria, corpus_info.corpus_id, corpus_info.is_placeholder)

        held_out = corpus_features.iloc[test_index]
        held_ids = set(held_out["candidate_id"])

        subsets = {
            "validated": groups["validated"][groups["validated"]["candidate_id"].isin(held_ids)],
            "shuffled": groups["shuffled"][groups["shuffled"]["source_id"].isin(held_ids)]
            if "source_id" in groups["shuffled"].columns else groups["shuffled"].head(0),
            "random": groups["random"].sample(
                n=min(len(held_ids), len(groups["random"])),
                random_state=int(cfg.get("run.seed")) + fold["fold"],
            ),
        }

        fold_scores = {}
        for name, subset in subsets.items():
            if subset.empty:
                continue
            scored = score_sequences(
                cfg, subset, refs,
                progress=(lambda i, n, s=f"E1 fold {fold['fold']} {name}":
                          progress(s, i, n) if progress else None),
            )
            usable = scored.dropna(subset=["tier1_score"])
            values = usable["tier1_score"].tolist()
            per_group[name].extend(values)
            if name == "validated":
                cluster_of = dict(zip(corpus_features["candidate_id"], target_labels))
                validated_clusters.extend(
                    str(cluster_of.get(cid, f"fold{fold['fold']}"))
                    for cid in usable["candidate_id"]
                )
            fold_scores[name] = {"n": len(values), "median": float(np.median(values))}

        fold_records.append(
            {"fold": fold["fold"], "grouping": fold["grouping"],
             "n_train": len(train), "n_test": len(test_index), "scores": fold_scores}
        )

    # The generated group has no fold structure: it is scored against the whole
    # corpus, exactly as a user would score it.
    generated_scores: list[float] = []
    if "generated" in groups.frames:
        refs_full = _refs_from(corpus_features, criteria, corpus_info.corpus_id,
                               corpus_info.is_placeholder)
        scored = score_sequences(
            cfg, groups["generated"], refs_full,
            progress=(lambda i, n: progress("E1 generated", i, n) if progress else None),
        )
        generated_scores = scored["tier1_score"].dropna().tolist()

    distributions = {**per_group, "generated": generated_scores}
    seed = int(cfg.get("run.seed"))
    comparisons = {
        "validated_vs_random": stats.compare_groups_clustered(
            per_group["validated"], validated_clusters, per_group["random"], seed=seed),
        "validated_vs_shuffled": stats.compare_groups_clustered(
            per_group["validated"], validated_clusters, per_group["shuffled"], seed=seed),
    }
    if generated_scores:
        comparisons["generated_vs_random"] = stats.compare_groups(
            generated_scores, per_group["random"], seed=seed)
        comparisons["validated_vs_generated"] = stats.compare_groups(
            per_group["validated"], generated_scores, seed=seed)

    return {
        "experiment": "E1",
        "title": "Tier 1 discriminates validated aptamers from matched controls",
        "design": "out-of-fold, folds grouped by target so no held-out sequence "
                  "is scored against its own target's family",
        "grouping_key": grouping_key,
        "n_unique_targets": int(pd.Series(validated_clusters).nunique())
        if validated_clusters else 0,
        "publication_eligible": not corpus_info.is_placeholder,
        "n_folds": len(fold_records),
        "folds": fold_records,
        "distributions": {k: v for k, v in distributions.items() if v},
        "comparisons": comparisons,
        "corpus": corpus_info.to_dict(),
        "runtime_seconds": round(time.perf_counter() - started, 1),
    }


# -- E2 ------------------------------------------------------------------


def e2_shuffled_controls(
    cfg: Config,
    corpus_features: pd.DataFrame,
    corpus_info: CorpusInfo,
    groups: ComparisonGroups,
    n_per_group: int = 300,
    progress: Callable[[str, int, int], None] | None = None,
) -> dict[str, Any]:
    """Is the structural score about arrangement, not composition? (E2)"""
    started = time.perf_counter()
    criteria = cfg.active_criteria()
    refs = _refs_from(corpus_features, criteria, corpus_info.corpus_id,
                      corpus_info.is_placeholder)
    n_shuffles = int(cfg.get("tier1.shuffle.n_shuffles"))
    alpha = float(cfg.get("tier1.shuffle.alpha"))

    results = {}
    for name, frame in groups.items():
        subset = frame if len(frame) <= n_per_group else frame.sample(
            n=n_per_group, random_state=int(cfg.get("run.seed"))
        )
        scored = score_sequences(
            cfg, subset, refs, n_shuffles=n_shuffles,
            progress=(lambda i, n, s=f"E2 {name}": progress(s, i, n) if progress else None),
        )
        passes = scored["shuffle_pass"].dropna()
        rate = float(passes.mean()) if len(passes) else float("nan")
        # Wilson interval: a pass rate near 0 or 1 has an asymmetric CI, and a
        # normal approximation would put the bound outside [0, 1].
        results[name] = {
            "n": int(len(passes)),
            "pass_rate": rate,
            "ci95": _wilson(int(passes.sum()), int(len(passes))),
            "median_p_value": float(scored["shuffle_p_value"].median()),
        }

    return {
        "experiment": "E2",
        "title": "Structural scores reflect arrangement, not composition",
        "alpha": alpha,
        "n_shuffles": n_shuffles,
        "note": f"pass = Monte-Carlo p <= {alpha} against {n_shuffles} "
                f"dinucleotide-preserving shuffles of the sequence itself",
        "groups": results,
        "runtime_seconds": round(time.perf_counter() - started, 1),
    }


def _wilson(successes: int, n: int, z: float = 1.96) -> list[float] | None:
    if n == 0:
        return None
    p = successes / n
    denom = 1 + z**2 / n
    centre = (p + z**2 / (2 * n)) / denom
    margin = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / denom
    return [float(max(0.0, centre - margin)), float(min(1.0, centre + margin))]


# -- E3 ------------------------------------------------------------------


def e3_matched_vs_decoy(
    cfg: Config,
    scored_validated: pd.DataFrame,
    bundles: Mapping[str, Mapping[str, Any]],
    bank: calibration.CalibrationBank,
    n_decoys: int = 5,
) -> dict[str, Any]:
    """Do aptamers score better against their true target than decoys? (E3)

    The one place with genuine labels, and the most likely to return a null
    result — which is publishable: it becomes an honest limitation and
    reinforces the plausibility-not-prediction framing.

    Scores are compared as **control percentiles**, never as raw geometric
    scores: a raw score is not comparable across targets, because each target
    has a different cavity and therefore a different control distribution.
    Pooling raw scores into one AUROC would be invalid.
    """
    started = time.perf_counter()
    a_per_nt = float(cfg.get("tier2.a_per_nt"))
    flex_c = float(cfg.get("tier2.flex_c"))
    sigma = float(cfg.get("tier2.sigma_A"))
    rng = np.random.default_rng(int(cfg.get("run.seed")))

    distributions = {}
    for pdb_id, bundle in bundles.items():
        d_pocket = float(bundle_mod.selected_pocket(bundle)["geometry"]["d_pocket_A"])
        distributions[pdb_id] = (
            d_pocket,
            calibration.target_distribution(bank, d_pocket, a_per_nt, flex_c, sigma),
        )

    available = sorted(distributions)
    if len(available) < 2:
        return {
            "experiment": "E3", "status": "skipped",
            "reason": f"need at least two target bundles, have {len(available)}",
        }

    pairs, positives, negatives, true_scores, decoy_scores = [], [], [], [], []
    for row in scored_validated.itertuples():
        true_target = str(getattr(row, "target_pdb_id", "") or "").upper()
        if true_target not in distributions:
            continue
        loop = float(row.loop_nt_median)
        if not np.isfinite(loop) or loop <= 0:
            continue

        decoy_pool = [t for t in available if t != true_target]
        chosen = list(rng.choice(decoy_pool, size=min(n_decoys, len(decoy_pool)), replace=False))

        true_p = _percentile_for(loop, distributions[true_target], a_per_nt, flex_c)
        decoys = [_percentile_for(loop, distributions[t], a_per_nt, flex_c) for t in chosen]

        pairs.append(
            {"candidate_id": row.candidate_id, "true_target": true_target,
             "true_percentile": true_p, "decoy_targets": chosen,
             "decoy_percentiles": decoys}
        )
        positives.append([true_p])
        negatives.append(decoys)
        true_scores.append(true_p)
        decoy_scores.append(decoys)

    if len(pairs) < 5:
        return {
            "experiment": "E3", "status": "skipped",
            "reason": f"only {len(pairs)} labelled aptamer/target pairs are usable",
        }

    differences = [p["true_percentile"] - float(np.mean(p["decoy_percentiles"])) for p in pairs]
    seed = int(cfg.get("run.seed"))

    return {
        "experiment": "E3",
        "status": "complete",
        "title": "True target vs decoy targets",
        "metric": "tier2_control_percentile (comparable across targets)",
        "n_pairs": len(pairs),
        "n_targets": len(available),
        "decoys_per_aptamer": n_decoys,
        "paired_difference": stats.paired_bootstrap(differences, seed=seed),
        "auroc": stats.bootstrap_auc_by_unit(positives, negatives, seed=seed),
        "retrieval": stats.mean_reciprocal_rank(true_scores, decoy_scores),
        "caveat": "Decoy targets are presumed non-matched, not confirmed "
                  "non-binders. A null result here is a limitation of the "
                  "geometric signal, not evidence that the aptamers do not bind.",
        "runtime_seconds": round(time.perf_counter() - started, 1),
    }


def _percentile_for(loop_nt: float, entry, a_per_nt: float, flex_c: float) -> float:
    from ..tier2.geometry import aptamer_dimensions

    d_pocket, dist = entry
    d_apt = aptamer_dimensions(loop_nt, a_per_nt, flex_c)["flexible"]
    return float(dist.percentile(abs(d_apt - d_pocket))[0])


# -- E4 ------------------------------------------------------------------


def e4_target_swappability(artifacts: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    """Same candidates, different targets — do the annotations actually change?

    Identical results across targets would mean the target input is not
    influencing anything, which would falsify the target-swappability claim.
    """
    from scipy.stats import spearmanr

    by_target = {}
    for artifact in artifacts:
        target = artifact.get("target")
        if not target:
            continue
        scores = {
            c["candidate_id"]: (c["tier2"] or {}).get("control_percentile_flexible")
            for c in artifact["candidates"]
            if (c["tier2"] or {}).get("status") == "evaluated"
        }
        bands = {
            c["candidate_id"]: (c["tier2"] or {}).get("band")
            for c in artifact["candidates"]
            if (c["tier2"] or {}).get("status") == "evaluated"
        }
        by_target[target["pdb_id"]] = {
            "run_id": artifact["run_id"],
            "d_pocket_A": target["selected_pocket"]["d_pocket_A"],
            "synthetic": target.get("synthetic", False),
            "scores": scores,
            "bands": bands,
        }

    if len(by_target) < 2:
        return {"experiment": "E4", "status": "skipped",
                "reason": f"need runs against at least two targets, have {len(by_target)}"}

    names = sorted(by_target)
    comparisons = []
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            shared = sorted(set(by_target[a]["scores"]) & set(by_target[b]["scores"]))
            if len(shared) < 3:
                continue
            x = [by_target[a]["scores"][c] for c in shared]
            y = [by_target[b]["scores"][c] for c in shared]
            changed = sum(
                1 for c in shared if by_target[a]["bands"][c] != by_target[b]["bands"][c]
            )
            rho = spearmanr(x, y)
            comparisons.append(
                {
                    "targets": [a, b],
                    "d_pocket_A": [by_target[a]["d_pocket_A"], by_target[b]["d_pocket_A"]],
                    "n_shared_candidates": len(shared),
                    "band_changed_fraction": changed / len(shared),
                    "spearman_rho": float(rho.statistic),
                    "spearman_p": float(rho.pvalue),
                }
            )

    return {
        "experiment": "E4",
        "status": "complete",
        "title": "Target swappability",
        "targets": {
            k: {"d_pocket_A": v["d_pocket_A"], "run_id": v["run_id"],
                "synthetic": v["synthetic"], "n_evaluated": len(v["scores"])}
            for k, v in by_target.items()
        },
        "comparisons": comparisons,
        "note": "Similar results are not necessarily a software failure: two "
                "targets with similar cavity dimensions should produce similar "
                "geometric annotations. Compare d_pocket_A before concluding.",
    }


# -- E5 ------------------------------------------------------------------


def e5_tier_independence(artifacts: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    """Spearman between the tiers, per target (spec §6.4, §8.2 E5)."""
    rows = []
    for artifact in artifacts:
        spearman = artifact["diagnostics"].get("spearman_tier1_tier2") or {}
        target = artifact.get("target") or {}
        rows.append(
            {
                "run_id": artifact["run_id"],
                "target": target.get("pdb_id"),
                "rho": spearman.get("rho"),
                "p_value": spearman.get("p_value"),
                "n": spearman.get("n"),
                "reason": spearman.get("reason"),
            }
        )
    usable = [r["rho"] for r in rows if r["rho"] is not None]
    return {
        "experiment": "E5",
        "title": "Tier independence",
        "per_run": rows,
        "max_abs_rho": float(max(abs(r) for r in usable)) if usable else None,
        "interpretation": (
            "These are weak observed monotonic associations within the "
            "surviving candidates — not evidence of statistical independence. "
            "The top-N cut restricts the range of Tier 1 scores among "
            "survivors, and both tiers consume the same predicted secondary "
            "structures, so some association is expected either way. A large "
            "|rho| would mean Tier 2 is largely restating Tier 1 and would "
            "change what the two-tier design can claim."
        ),
    }


def write_results(results: Mapping[str, Any], path: str | Path) -> Path:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(results, indent=2, default=str), encoding="utf-8")
    return out
