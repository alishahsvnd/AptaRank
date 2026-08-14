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

The bank's descriptors are target- *and* mode-independent, so the expensive part
— folding and ensemble sampling — is done once and cached alongside the corpus.
Only the arithmetic that involves the target is redone per target, and switching
binding mode just reads a different column of the same bank.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence

import numpy as np
import pandas as pd

from ..config import Config
from ..errors import CorpusError
from ..provenance import derive_seed, sha256_text
from ..tier1 import features as feature_mod
from ..tier1 import folding, shuffles
from ..tier1 import corpus as corpus_cache
from ..tier1.corpus import CorpusInfo
from . import modes

BANK_SCHEMA_VERSION = "calibration-bank-v1"


@dataclass
class CalibrationBank:
    """Fixed shuffled controls with the descriptors every mode compares.

    Both descriptors are folded once and cached together, so switching binding
    mode costs no recomputation — only the arithmetic that involves the target
    is redone.
    """

    bank_id: str
    loop_nt_median: np.ndarray
    meta: dict[str, Any]
    length: np.ndarray | None = None
    rg_median_A: np.ndarray | None = None

    @property
    def size(self) -> int:
        return int(self.loop_nt_median.size)

    def descriptor(self, mode: str, params: Mapping[str, Any] | None = None) -> np.ndarray:
        """The control values this mode compares (refinements §4.2)."""
        column = modes.descriptor_column(modes.check_mode(mode), params or {})
        values = {
            "loop_nt_median": self.loop_nt_median,
            "length": self.length,
            "rg_median_A": self.rg_median_A,
        }.get(column)
        if values is None:
            raise CorpusError(
                f"the cached calibration bank has no {column!r} column, which "
                f"{mode} mode needs. Delete the bank cache to rebuild it."
            )
        return np.asarray(values, dtype=float)


@dataclass
class ControlDistribution:
    """The bank's disagreements against one specific target, in one mode.

    "Disagreement" is whatever the mode compares — an Å mismatch between loop
    reach and cavity width, an Å² mismatch between footprint and patch area —
    always with lower meaning better agreement, so one percentile definition
    serves every mode.
    """

    sorted_disagreement: np.ndarray
    bank_id: str
    mode: str
    units: str
    n: int
    target: dict[str, Any] = field(default_factory=dict)

    def percentile(self, disagreement: float | np.ndarray) -> np.ndarray:
        """One-sided mid-rank percentile, higher = better geometric agreement.

            P = (#{Δ_control > Δ} + 0.5 #{Δ_control = Δ}) / B

        Defined on the raw disagreement rather than the Gaussian display score:
        the Gaussian underflows for large mismatches and depends on `sigma`, and
        a band must not move because a display parameter changed.
        """
        x = np.atleast_1d(np.asarray(disagreement, dtype=float))
        ref = self.sorted_disagreement
        left = np.searchsorted(ref, x, side="left")
        right = np.searchsorted(ref, x, side="right")
        greater = ref.size - right
        equal = right - left
        return (greater + 0.5 * equal) / ref.size

    def quantile_disagreement(self, percentile: float) -> float:
        """The disagreement at a given control percentile — a UI reference line."""
        return float(np.quantile(self.sorted_disagreement, 1.0 - percentile))


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
        table = corpus_cache.read_feature_cache(csv_path)
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        return CalibrationBank(
            bank_id=bank_id,
            loop_nt_median=table["loop_nt_median"].to_numpy(dtype=float),
            length=(
                table["length"].to_numpy(dtype=float)
                if "length" in table.columns else None
            ),
            rg_median_A=(
                table["rg_median_A"].to_numpy(dtype=float)
                if "rg_median_A" in table.columns else None
            ),
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
                a_per_bp_helix=float(cfg.get("tier2.geometry.a_per_bp_helix")),
                a_per_nt_ss=float(cfg.get("tier2.geometry.a_per_nt_ss")),
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
            "length", "loop_nt_median", "loop_nt_p90", "loop_nt_iqr",
            "rg_median_A", "rg_iqr_A"]
    table = table[keep]

    cache_dir.mkdir(parents=True, exist_ok=True)
    corpus_cache.write_feature_cache(table, csv_path)
    # Same reason as the corpus cache: project the bank from the bytes that were
    # stored, so the run that builds it bands candidates exactly like the runs
    # that reuse it.
    table = corpus_cache.read_feature_cache(csv_path)
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
        "size_summary_convention": "median radius of gyration across sampled "
                                   "structures, mass-weighted over the forgi "
                                   "element graph (stems as rigid A-form rods, "
                                   "unpaired regions as Gaussian segments)",
        "geometry": {
            "a_per_bp_helix": cfg.get("tier2.geometry.a_per_bp_helix"),
            "a_per_nt_ss": cfg.get("tier2.geometry.a_per_nt_ss"),
        },
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
        length=table["length"].to_numpy(dtype=float),
        rg_median_A=table["rg_median_A"].to_numpy(dtype=float),
        meta=meta,
    )


def target_distribution(
    bank: CalibrationBank,
    mode: str,
    target: Mapping[str, Any],
    params: Mapping[str, Any],
) -> ControlDistribution:
    """Project the bank onto one target in one mode. Cheap: pure arithmetic."""
    modes.check_mode(mode)
    values = bank.descriptor(mode, params)
    disagreements = np.array(
        [
            modes.compare(mode, float(v), target, params)["disagreement"]
            for v in values
            if modes.is_evaluable(v)
        ],
        dtype=float,
    )
    if disagreements.size == 0:
        raise CorpusError(
            f"no calibration control produced a usable {mode}-mode measurement; "
            f"the bank cannot calibrate this target"
        )
    return ControlDistribution(
        sorted_disagreement=np.sort(disagreements),
        bank_id=bank.bank_id,
        mode=mode,
        units=modes.MISMATCH_UNITS[mode],
        n=int(disagreements.size),
        target=dict(target),
    )


def secondary_distributions(
    bank: CalibrationBank,
    mode: str,
    target: Mapping[str, Any],
    params: Mapping[str, Any],
) -> dict[str, ControlDistribution]:
    """Sensitivity distributions reported alongside the primary one.

    Pocket mode keeps contour length as the documented upper bound on loop
    reach; it is reported, never banded on.
    """
    if mode != modes.POCKET or params.get("primary_descriptor") != "flexible":
        return {}
    extended = dict(params)
    extended["primary_descriptor"] = "extended"
    return {"extended": target_distribution(bank, mode, target, extended)}


def thresholds(dist: ControlDistribution, moderate: float, strong: float) -> dict[str, Any]:
    """Band boundaries in the mode's own units, for the dashboard's reference lines."""
    out = {
        "bank_id": dist.bank_id,
        "n_controls": dist.n,
        "binding_mode": dist.mode,
        "units": dist.units,
        "target": dist.target,
        "band_percentiles": {"moderate": moderate, "strong": strong},
        "disagreement_at_moderate": dist.quantile_disagreement(moderate),
        "disagreement_at_strong": dist.quantile_disagreement(strong),
        "control_disagreement_median": float(np.median(dist.sorted_disagreement)),
    }
    if dist.mode == modes.POCKET:
        # Names the pocket-mode dashboard and earlier artifacts already use.
        out.update(
            d_pocket_A=dist.target.get("d_pocket_A"),
            mismatch_at_moderate_A=out["disagreement_at_moderate"],
            mismatch_at_strong_A=out["disagreement_at_strong"],
            control_mismatch_median_A=out["control_disagreement_median"],
        )
    return out


def score_candidates(
    values: Sequence[float],
    mode: str,
    target: Mapping[str, Any],
    params: Mapping[str, Any],
    distribution: ControlDistribution,
    moderate: float,
    strong: float,
    secondary: Mapping[str, ControlDistribution] | None = None,
) -> list[dict[str, Any]]:
    """Per-candidate Tier 2 record. The band comes from the control percentile."""
    modes.check_mode(mode)
    secondary = dict(secondary or {})
    records = []
    for raw in values:
        if not modes.is_evaluable(raw):
            # A missing descriptor is not a zero-sized aptamer — it is an
            # absence of the measurement this mode depends on.
            value = float(raw) if np.isfinite(_as_float(raw)) else None
            records.append(
                {
                    "status": modes.not_evaluable_status(mode),
                    "band": "not_evaluated",
                    "binding_mode": mode,
                    modes.descriptor_column(mode, params): value,
                }
            )
            continue

        result = modes.compare(mode, float(raw), target, params)
        percentile = float(distribution.percentile(result["disagreement"])[0])
        record = {
            "status": "evaluated",
            "binding_mode": mode,
            **result["fields"],
            "disagreement": result["disagreement"],
            "disagreement_units": distribution.units,
            # Paper-side vocabulary (§1.1). `score` stays as the display alias.
            "geometric_agreement_score": result["agreement"],
            "geometric_agreement_percentile": percentile,
            "control_percentile": percentile,
            "score": percentile,
            "band": assign_band(percentile, moderate, strong),
        }
        for name, dist in secondary.items():
            variant = dict(params)
            variant["primary_descriptor"] = name
            alternative = modes.compare(mode, float(raw), target, variant)
            record[f"control_percentile_{name}"] = float(
                dist.percentile(alternative["disagreement"])[0]
            )
        if mode == modes.POCKET:
            # The banded quantity, under the name the dashboard and E3/E4 use.
            record["control_percentile_flexible"] = percentile
        records.append(record)
    return records


def _as_float(value: Any) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")
