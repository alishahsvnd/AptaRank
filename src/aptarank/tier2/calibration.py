"""Control calibration for Tier 2 bands (spec §5.8, revised).

The raw compatibility score lives on an arbitrary scale whose meaning depends
on the target's pocket size, so a fixed cutoff would behave inconsistently the
moment the target changes — silently breaking the target-swappability claim.
Bands are therefore defined relative to shuffled controls.

Deviation from §5.8, deliberately: the spec derives the thresholds from the
shuffles of the *submitted batch*, which makes a candidate's band depend on
what else happened to be submitted alongside it, and leaves the thresholds
undefined for small batches. Instead we bank a fixed set of dinucleotide
shuffles drawn from the reference corpus, once, and reuse it for every target
and every batch.

The bank's loop descriptors are target-independent, so the expensive part —
folding and ensemble sampling — is done once and cached alongside the corpus.
Only the arithmetic that involves the pocket is redone per target.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Sequence

import numpy as np
import pandas as pd

from ..config import Config
from ..errors import CorpusError
from ..provenance import derive_seed, sha256_text
from ..tier1 import features as feature_mod
from ..tier1 import folding, shuffles
from ..tier1.corpus import CorpusInfo
from .geometry import aptamer_dimensions, compatibility

BANK_SCHEMA_VERSION = "calibration-bank-v1"


@dataclass
class CalibrationBank:
    """Fixed shuffled controls with their ensemble loop descriptors."""

    bank_id: str
    loop_nt_median: np.ndarray
    meta: dict[str, Any]

    @property
    def size(self) -> int:
        return int(self.loop_nt_median.size)


@dataclass
class ControlDistribution:
    """The bank's absolute mismatches against one specific target pocket."""

    sorted_absolute_mismatch_A: np.ndarray
    bank_id: str
    d_pocket_A: float
    descriptor: str
    n: int

    def percentile(self, absolute_mismatch: float | np.ndarray) -> np.ndarray:
        """One-sided mid-rank percentile, higher = better geometric agreement.

            P = (#{Δ_control > Δ} + 0.5 #{Δ_control = Δ}) / B

        Defined on the absolute mismatch rather than the Gaussian score: the
        Gaussian underflows for large mismatches and depends on `sigma`, and a
        band should not move because a display parameter changed.
        """
        x = np.atleast_1d(np.asarray(absolute_mismatch, dtype=float))
        ref = self.sorted_absolute_mismatch_A
        left = np.searchsorted(ref, x, side="left")
        right = np.searchsorted(ref, x, side="right")
        greater = ref.size - right
        equal = right - left
        return (greater + 0.5 * equal) / ref.size

    def quantile_mismatch(self, percentile: float) -> float:
        """The mismatch at a given control percentile — a UI reference line."""
        return float(np.quantile(self.sorted_absolute_mismatch_A, 1.0 - percentile))


def assign_band(percentile: float | None, moderate: float, strong: float) -> str:
    """Three graded bands, never a binary in/out (spec §5.8).

    A hard cutoff would imply more precision than this signal has.
    """
    if percentile is None or not np.isfinite(percentile):
        return "not_evaluated"
    if percentile >= strong:
        return "strong"
    if percentile >= moderate:
        return "moderate"
    return "weak"


# -- building the bank ---------------------------------------------------


def bank_signature(cfg: Config, corpus_info: CorpusInfo) -> str:
    relevant = {
        "schema": BANK_SCHEMA_VERSION,
        "corpus_sha256": corpus_info.corpus_sha256,
        "tool_signature": corpus_info.tool_signature,
        "model": folding.model_settings(),
        "bank_size": cfg.get("tier2.calibration.bank_size"),
        "k": cfg.get("tier1.shuffle.k"),
        "n_ensemble_samples": cfg.get("tier1.n_ensemble_samples"),
        "seed": cfg.get("run.seed"),
        "loop_definition": "max_loop_nt over h/i/m elements, forgi",
    }
    return sha256_text(json.dumps(relevant, sort_keys=True, default=str))[:16]


def build_or_load(
    cfg: Config,
    corpus_table: pd.DataFrame,
    corpus_info: CorpusInfo,
    progress: Callable[[int, int], None] | None = None,
) -> CalibrationBank:
    """Fold a fixed set of corpus-derived shuffles, once, and cache them."""
    signature = bank_signature(cfg, corpus_info)
    bank_id = f"{corpus_info.corpus_id}_{signature}"
    cache_dir = Path(cfg.get("tier2.calibration.cache_dir"))
    csv_path = cache_dir / f"{bank_id}.csv"
    meta_path = cache_dir / f"{bank_id}.meta.json"

    if csv_path.exists() and meta_path.exists():
        table = pd.read_csv(csv_path)
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        return CalibrationBank(
            bank_id=bank_id,
            loop_nt_median=table["loop_nt_median"].to_numpy(dtype=float),
            meta=meta,
        )

    if "sequence" not in corpus_table.columns:
        raise CorpusError("corpus feature table has no sequence column")

    size = int(cfg.get("tier2.calibration.bank_size"))
    k = int(cfg.get("tier1.shuffle.k"))
    n_samples = int(cfg.get("tier1.n_ensemble_samples"))
    seed = int(cfg.get("run.seed"))

    # Equal allocation across the corpus, not a truncated round robin: taking
    # exactly `size` controls in sorted order would give the lexicographically
    # early sequences an extra control each, and sequence order correlates with
    # sequence content. Every source therefore contributes the same number,
    # which overshoots `size` slightly and weights the bank evenly.
    sources = (
        corpus_table[["candidate_id", "sequence"]]
        .drop_duplicates("sequence")
        .sort_values("sequence")
        .reset_index(drop=True)
    )
    if sources.empty:
        raise CorpusError("cannot build a calibration bank from an empty corpus")

    per_source = max(1, -(-size // len(sources)))   # ceil
    total = per_source * len(sources)

    jobs, provenance = [], []
    for i in range(total):
        row = sources.iloc[i % len(sources)]
        shuffle_index = i // len(sources)
        control_seed = derive_seed(seed, row.candidate_id, "bank", shuffle_index)
        shuffled = shuffles.generate_shuffles(row.sequence, 1, k, control_seed)[0]
        control_id = f"bank{i:06d}"
        jobs.append(
            feature_mod.FeatureJob(
                candidate_id=control_id,
                sequence=shuffled,
                n_ensemble_samples=n_samples,
                n_shuffles=0,
                seed=seed,
            )
        )
        provenance.append(
            {
                "control_id": control_id,
                "source_id": row.candidate_id,
                "shuffle_index": shuffle_index,
                "shuffle_seed": control_seed,
                "identical_to_source": shuffled == row.sequence,
            }
        )

    results = feature_mod.compute_batch(
        jobs,
        workers=cfg.get("tier1.parallel.workers", None),
        chunk_size=int(cfg.get("tier1.parallel.chunk_size", 16)),
        progress=progress,
    )
    table, failures = feature_mod.results_to_frame(results)
    if table.empty:
        raise CorpusError("every calibration-bank control failed to fold")

    table = table.merge(pd.DataFrame(provenance), left_on="candidate_id",
                        right_on="control_id", how="left")
    keep = ["control_id", "source_id", "shuffle_index", "identical_to_source",
            "length", "loop_nt_median", "loop_nt_p90", "loop_nt_iqr"]
    table = table[keep]

    cache_dir.mkdir(parents=True, exist_ok=True)
    table.to_csv(csv_path, index=False)
    meta = {
        "bank_id": bank_id,
        "schema_version": BANK_SCHEMA_VERSION,
        "corpus_id": corpus_info.corpus_id,
        "corpus_sha256": corpus_info.corpus_sha256,
        "corpus_is_placeholder": corpus_info.is_placeholder,
        "requested_size": size,
        "controls_per_source": per_source,
        "allocation": "equal per source (ceil(size / n_sources) each)",
        "loop_summary_convention": "median of max_loop_nt over h/i/m forgi "
                                   "elements across sampled structures; "
                                   "numpy.median, linear interpolation",
        "n_controls": len(table),
        "n_failed": len(failures),
        "n_identical_to_source": int(table["identical_to_source"].sum()),
        "n_source_sequences": len(sources),
        "k": k,
        "n_ensemble_samples": n_samples,
        "seed": seed,
        "model": folding.model_settings(),
    }
    meta_path.write_text(json.dumps(meta, indent=2, default=str), encoding="utf-8")

    return CalibrationBank(
        bank_id=bank_id,
        loop_nt_median=table["loop_nt_median"].to_numpy(dtype=float),
        meta=meta,
    )


def target_distribution(
    bank: CalibrationBank,
    d_pocket_A: float,
    a_per_nt: float,
    flex_c: float,
    sigma_A: float,
    descriptor: str = "flexible",
) -> ControlDistribution:
    """Project the bank onto one target. Pure arithmetic — cheap per target."""
    dims = np.array(
        [
            aptamer_dimensions(float(loop), a_per_nt, flex_c)[descriptor]
            for loop in bank.loop_nt_median
        ]
    )
    mismatches = np.abs(dims - float(d_pocket_A))
    return ControlDistribution(
        sorted_absolute_mismatch_A=np.sort(mismatches),
        bank_id=bank.bank_id,
        d_pocket_A=float(d_pocket_A),
        descriptor=descriptor,
        n=int(mismatches.size),
    )


def thresholds(dist: ControlDistribution, moderate: float, strong: float) -> dict[str, Any]:
    """Band boundaries expressed in Å, for the dashboard's reference lines."""
    return {
        "bank_id": dist.bank_id,
        "n_controls": dist.n,
        "descriptor": dist.descriptor,
        "d_pocket_A": dist.d_pocket_A,
        "band_percentiles": {"moderate": moderate, "strong": strong},
        "mismatch_at_moderate_A": dist.quantile_mismatch(moderate),
        "mismatch_at_strong_A": dist.quantile_mismatch(strong),
        "control_mismatch_median_A": float(np.median(dist.sorted_absolute_mismatch_A)),
    }


def score_candidates(
    loop_nt_medians: Sequence[float],
    d_pocket_A: float,
    dist_flexible: ControlDistribution,
    dist_extended: ControlDistribution,
    a_per_nt: float,
    flex_c: float,
    sigma_A: float,
    moderate: float,
    strong: float,
) -> list[dict[str, Any]]:
    """Per-candidate Tier 2 record. The band comes from the control percentile."""
    records = []
    for loop in loop_nt_medians:
        loop = float(loop)
        if not np.isfinite(loop) or loop <= 0:
            # No accessible loop is not a zero-sized aptamer — it is an absence
            # of the measurement Tier 2 depends on.
            records.append(
                {"status": "not_evaluable_no_contact_loop", "band": "not_evaluated",
                 "loop_nt_median": loop if np.isfinite(loop) else None}
            )
            continue

        dims = aptamer_dimensions(loop, a_per_nt, flex_c)
        flex = compatibility(dims["flexible"], d_pocket_A, sigma_A)
        ext = compatibility(dims["extended"], d_pocket_A, sigma_A)
        p_flex = float(dist_flexible.percentile(flex["absolute_mismatch_A"])[0])
        p_ext = float(dist_extended.percentile(ext["absolute_mismatch_A"])[0])

        records.append(
            {
                "status": "evaluated",
                "loop_nt_median": loop,
                "d_pocket_A": float(d_pocket_A),
                "d_apt_flexible_A": dims["flexible"],
                "signed_mismatch_flexible_A": flex["signed_mismatch_A"],
                "absolute_mismatch_flexible_A": flex["absolute_mismatch_A"],
                "geometric_score_flexible": flex["geometric_score"],
                "control_percentile_flexible": p_flex,
                "d_apt_extended_A": dims["extended"],
                "signed_mismatch_extended_A": ext["signed_mismatch_A"],
                "absolute_mismatch_extended_A": ext["absolute_mismatch_A"],
                "geometric_score_extended": ext["geometric_score"],
                "control_percentile_extended": p_ext,
                # Primary, for display and for the explanation templates.
                "d_apt_A": dims["flexible"],
                "difference_A": flex["signed_mismatch_A"],
                "score": p_flex,
                "band": assign_band(p_flex, moderate, strong),
            }
        )
    return records
