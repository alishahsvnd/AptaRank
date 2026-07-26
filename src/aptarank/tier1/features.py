"""Per-sequence feature extraction and the parallel batch driver (spec §4.4–4.8).

Workers do folding only. All scoring is done centrally and vectorised, so a
worker never needs the corpus and the corpus is never shipped to a subprocess.

Feature groups per sequence:
  composition   length, gc_fraction                                    (§4.4)
  folding       dot_bracket, mfe, mfe_norm, ensemble_defect,
                positional_entropy_mean                                (§4.5)
  elements      stem_fraction, loop counts, max_loop_nt, ...           (§4.6)
  ensemble      loop_nt_median / p90 / iqr over sampled structures     (§4.7)
"""

from __future__ import annotations

import math
import os
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from typing import Any, Callable, Iterable, Sequence

import numpy as np
import pandas as pd

from ..errors import FoldingError
from ..provenance import derive_seed
from . import elements as element_mod
from . import folding, shuffles

# Columns every feature table is guaranteed to carry.
FEATURE_COLUMNS = (
    "length",
    "gc_fraction",
    "mfe",
    "mfe_norm",
    "ensemble_free_energy",
    "ensemble_defect",
    "positional_entropy_mean",
    "n_hairpins",
    "n_interior",
    "n_multiloop",
    "n_stems",
    "stem_fraction",
    "longest_stem_bp",
    "max_loop_nt",
    "total_unpaired",
)

ENSEMBLE_COLUMNS = ("loop_nt_median", "loop_nt_p90", "loop_nt_iqr", "n_ensemble_samples")


@dataclass(frozen=True)
class FeatureJob:
    """One unit of folding work handed to a worker process."""

    candidate_id: str
    sequence: str
    n_ensemble_samples: int = 0
    n_shuffles: int = 0
    shuffle_k: int = 2
    seed: int = 0


def gc_fraction(sequence: str) -> float:
    """Pure character counting — the whole of spec §4.4's composition feature."""
    if not sequence:
        return 0.0
    return (sequence.count("G") + sequence.count("C")) / len(sequence)


def sequence_features(
    sequence: str, n_ensemble_samples: int = 0, seed: int | None = None
) -> dict[str, Any]:
    """Everything computable from one sequence without the corpus or a target."""
    fold_result, samples = folding.fold_and_sample(sequence, n_ensemble_samples, seed)
    elems = element_mod.parse_elements(fold_result.dot_bracket, sequence)

    record: dict[str, Any] = {
        "sequence": sequence,
        "length": len(sequence),
        "gc_fraction": gc_fraction(sequence),
        **fold_result.feature_dict(),
        **elems.feature_dict(),
        "element_string": elems.element_string,
    }
    record.update(ensemble_loop_stats(samples, fallback=elems.max_loop_nt))

    # A non-finite feature must fail its candidate here rather than reach the
    # percentile scorer, where NaN would silently become a plausible-looking
    # extreme score.
    bad = [
        name for name in FEATURE_COLUMNS
        if not math.isfinite(float(record[name]))
    ]
    if bad:
        raise FoldingError(f"non-finite feature value(s) {bad} for {sequence!r}")
    return record


def ensemble_loop_stats(samples: Sequence[str], fallback: int) -> dict[str, Any]:
    """Loop-size distribution over sampled structures (spec §4.7).

    Tier 2 consumes `loop_nt_median` rather than the single MFE structure's
    loop size: real molecules move between shapes, and a distribution statistic
    is a more honest input to an already-coarse geometric comparison.
    `loop_nt_iqr` says how uncertain that number itself is.
    """
    if not samples:
        return {
            "loop_nt_median": float(fallback),
            "loop_nt_p90": float(fallback),
            "loop_nt_iqr": 0.0,
            "n_ensemble_samples": 0,
        }
    sizes = np.array([element_mod.max_loop_nt(db) for db in samples], dtype=float)
    q25, q75 = np.percentile(sizes, [25, 75])
    return {
        "loop_nt_median": float(np.median(sizes)),
        "loop_nt_p90": float(np.percentile(sizes, 90)),
        "loop_nt_iqr": float(q75 - q25),
        "n_ensemble_samples": int(sizes.size),
    }


def shuffle_features(sequence: str, n: int, k: int, seed: int) -> list[dict[str, Any]]:
    """Fold each shuffled control. No ensemble sampling — too expensive (§4.8)."""
    out = []
    for i, shuffled in enumerate(shuffles.generate_shuffles(sequence, n, k, seed)):
        record = sequence_features(shuffled, n_ensemble_samples=0)
        record["shuffle_index"] = i
        out.append(record)
    return out


def run_job(job: FeatureJob) -> dict[str, Any]:
    """Worker entry point. Must stay top-level and picklable (Windows spawn)."""
    try:
        # The sequence is part of the seed material as well as the id: a reused
        # user-supplied id would otherwise give two different sequences the
        # same random stream.
        record = sequence_features(
            job.sequence,
            n_ensemble_samples=job.n_ensemble_samples,
            seed=derive_seed(job.seed, job.candidate_id, job.sequence, "ensemble"),
        )
        controls = shuffle_features(
            job.sequence,
            job.n_shuffles,
            job.shuffle_k,
            derive_seed(job.seed, job.candidate_id, job.sequence, "shuffle"),
        )
        return {
            "candidate_id": job.candidate_id,
            "features": record,
            "shuffles": controls,
            "error": None,
        }
    except Exception as exc:  # a single bad candidate must not kill the batch
        return {
            "candidate_id": job.candidate_id,
            "features": None,
            "shuffles": [],
            "error": f"{type(exc).__name__}: {exc}",
        }


def compute_batch(
    jobs: Iterable[FeatureJob],
    workers: int | None = None,
    chunk_size: int = 16,
    progress: Callable[[int, int], None] | None = None,
) -> list[dict[str, Any]]:
    """Fold a batch of sequences, in parallel across processes.

    Processes, not threads: the ViennaRNA bindings are C extensions that
    release the GIL inconsistently, so threads give almost no speed-up.
    """
    job_list = list(jobs)
    total = len(job_list)
    if total == 0:
        return []

    n_workers = workers or max(1, (os.cpu_count() or 2) - 1)
    if n_workers == 1 or total < 4:
        results = []
        for i, job in enumerate(job_list, start=1):
            results.append(run_job(job))
            if progress:
                progress(i, total)
        return results

    results = []
    with ProcessPoolExecutor(max_workers=n_workers) as pool:
        for i, result in enumerate(pool.map(run_job, job_list, chunksize=chunk_size), start=1):
            results.append(result)
            if progress:
                progress(i, total)
    return results


def results_to_frame(results: Sequence[dict[str, Any]]) -> tuple[pd.DataFrame, list[dict[str, Any]]]:
    """Split worker output into a feature table and a list of failures."""
    rows, failures = [], []
    for res in results:
        if res["error"] is not None or res["features"] is None:
            failures.append({"candidate_id": res["candidate_id"], "reason": res["error"]})
            continue
        rows.append({"candidate_id": res["candidate_id"], **res["features"]})
    frame = pd.DataFrame(rows)
    return frame, failures
