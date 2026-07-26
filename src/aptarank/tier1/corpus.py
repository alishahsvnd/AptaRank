"""Reference corpus: load, precompute features once, cache, serve distributions.

Spec §3.5.1. The corpus of experimentally validated aptamers defines what
"normal" means. No threshold in this project is hard-coded; every criterion
score is a position within these distributions.

Precomputing the corpus feature matrix and caching it to disk is the single
most important performance decision in Tier 1 — otherwise every candidate run
re-folds several thousand corpus sequences.

The cache key covers the corpus contents, the feature schema, and the versions
of the tools that produced the numbers. A stale cache would be a silent
methodological error, not a performance bug.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Sequence

import numpy as np
import pandas as pd

from ..config import CRITERIA, Config
from ..errors import CorpusError
from ..ingest import normalise_sequence, validate_sequence
from ..provenance import sha256_file, sha256_text, tool_versions
from . import features as feature_mod
from . import folding
from .scoring import ReferenceDistributions

FEATURE_SCHEMA_VERSION = "1"
REQUIRED_COLUMNS = ("sequence",)
OPTIONAL_COLUMNS = ("id", "target_name", "target_pdb_id", "source_reference")


@dataclass
class CorpusInfo:
    """Everything needed to decide whether two runs are comparable."""

    corpus_id: str
    path: str
    corpus_sha256: str
    cache_sha256: str
    n_sequences: int
    n_dropped: int
    feature_schema_version: str
    tool_signature: str
    is_placeholder: bool
    dropped_reasons: dict[str, int] = field(default_factory=dict)

    @property
    def publication_eligible(self) -> bool:
        return not self.is_placeholder

    def to_dict(self) -> dict[str, Any]:
        return {
            "corpus_id": self.corpus_id,
            "path": self.path,
            "corpus_sha256": self.corpus_sha256,
            "cache_sha256": self.cache_sha256,
            "n_sequences": self.n_sequences,
            "n_dropped": self.n_dropped,
            "dropped_reasons": self.dropped_reasons,
            "feature_schema_version": self.feature_schema_version,
            "tool_signature": self.tool_signature,
            "is_placeholder": self.is_placeholder,
            "publication_eligible": self.publication_eligible,
        }


def tool_signature() -> str:
    """Hash of everything that determines the corpus feature values.

    Tool versions alone are not enough: a change to the folding model
    (temperature, dangles, noLP, ...) changes every number in the cache while
    leaving the version strings identical. A cache that survived that would be
    a silent methodological error, so the model settings are in the key.
    """
    versions = tool_versions()
    relevant = {
        "tools": {k: versions.get(k) for k in ("viennarna", "forgi")},
        "model": folding.model_settings(),
        "features": list(feature_mod.FEATURE_COLUMNS),
    }
    return sha256_text(json.dumps(relevant, sort_keys=True, default=str))[:12]


def load_corpus_table(path: str | Path, min_length: int, max_length: int) -> tuple[pd.DataFrame, dict[str, int]]:
    """Read and clean the corpus CSV. Returns (table, dropped-reason counts)."""
    p = Path(path)
    if not p.exists():
        raise CorpusError(f"reference corpus not found: {p}")

    frame = pd.read_csv(p, dtype=str, keep_default_na=False)
    frame.columns = [c.strip().lower() for c in frame.columns]
    missing = [c for c in REQUIRED_COLUMNS if c not in frame.columns]
    if missing:
        raise CorpusError(
            f"{p}: corpus is missing required column(s) {missing}. "
            f"Expected: id, sequence, target_name, target_pdb_id, source_reference"
        )
    for col in OPTIONAL_COLUMNS:
        if col not in frame.columns:
            frame[col] = ""

    frame["sequence"] = frame["sequence"].map(normalise_sequence)
    reasons = frame["sequence"].map(lambda s: validate_sequence(s, min_length, max_length))
    dropped = reasons.dropna()
    counts: dict[str, int] = {}
    for reason in dropped:
        key = reason.split(" below ")[0].split(" above ")[0]
        key = "length out of range" if key.startswith("length") else reason
        counts[key] = counts.get(key, 0) + 1

    clean = frame[reasons.isna()].copy()
    clean = clean.drop_duplicates(subset="sequence").reset_index(drop=True)
    if clean.empty:
        raise CorpusError(f"{p}: no usable sequences after validation ({counts})")

    clean["id"] = [
        cid if cid else f"ref{i:05d}" for i, cid in enumerate(clean["id"], start=1)
    ]
    return clean, counts


def build_or_load(
    cfg: Config,
    force: bool = False,
    progress: Callable[[int, int], None] | None = None,
) -> tuple[pd.DataFrame, CorpusInfo]:
    """Return the corpus feature matrix, computing and caching it if needed."""
    path = cfg.get("corpus.path", None)
    if not path:
        raise CorpusError(
            "no reference corpus configured. Tier 1 percentile scoring is blocked "
            "without it (spec §3.5.1). Set corpus.path in your config, or pass "
            "--development-corpus PATH to run against synthetic data (development "
            "only; such runs are stamped publication_eligible: false)."
        )

    is_placeholder = bool(cfg.get("corpus.is_placeholder", False))
    if is_placeholder and not cfg.get("corpus.allow_placeholder", False):
        raise CorpusError(
            "corpus.is_placeholder is set but corpus.allow_placeholder is false. "
            "Refusing to score real candidates against synthetic reference data."
        )

    min_len, max_len = cfg.get("input.min_length"), cfg.get("input.max_length")
    corpus_sha = sha256_file(path)
    signature = tool_signature()
    key = f"{Path(path).stem}_{corpus_sha[:12]}_{FEATURE_SCHEMA_VERSION}_{signature}_{min_len}-{max_len}"

    cache_dir = Path(cfg.get("corpus.cache_dir"))
    cache_csv = cache_dir / f"{key}.csv"
    cache_meta = cache_dir / f"{key}.meta.json"

    if cache_csv.exists() and cache_meta.exists() and not force:
        table = pd.read_csv(cache_csv)
        meta = json.loads(cache_meta.read_text(encoding="utf-8"))
        meta["is_placeholder"] = is_placeholder
        meta["path"] = str(path)
        return table, CorpusInfo(**meta)

    clean, dropped = load_corpus_table(path, min_len, max_len)
    jobs = [
        feature_mod.FeatureJob(
            candidate_id=str(row.id),
            sequence=str(row.sequence),
            n_ensemble_samples=0,   # the corpus defines reference distributions
            n_shuffles=0,           # only; it needs neither sampling nor controls
            seed=int(cfg.get("run.seed")),
        )
        for row in clean.itertuples()
    ]
    results = feature_mod.compute_batch(
        jobs,
        workers=cfg.get("tier1.parallel.workers", None),
        chunk_size=int(cfg.get("tier1.parallel.chunk_size", 16)),
        progress=progress,
    )
    table, failures = feature_mod.results_to_frame(results)
    if failures:
        dropped["folding failed"] = len(failures)
    if table.empty:
        raise CorpusError(f"every corpus sequence failed to fold ({len(failures)} failures)")

    meta_columns = ["id", "target_name", "target_pdb_id", "source_reference"]
    table = table.merge(
        clean[meta_columns].rename(columns={"id": "candidate_id"}),
        on="candidate_id",
        how="left",
    )

    cache_dir.mkdir(parents=True, exist_ok=True)
    table.to_csv(cache_csv, index=False)
    info = CorpusInfo(
        corpus_id=key,
        path=str(path),
        corpus_sha256=corpus_sha,
        cache_sha256=sha256_file(cache_csv),
        n_sequences=len(table),
        n_dropped=int(sum(dropped.values())),
        dropped_reasons=dropped,
        feature_schema_version=FEATURE_SCHEMA_VERSION,
        tool_signature=signature,
        is_placeholder=is_placeholder,
    )
    payload = info.to_dict()
    payload.pop("publication_eligible")
    cache_meta.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return table, info


def reference_distributions(
    table: pd.DataFrame,
    info: CorpusInfo,
    criteria: Sequence[str] | None = None,
) -> ReferenceDistributions:
    """Turn the cached corpus feature matrix into per-criterion distributions."""
    names = list(criteria) if criteria else list(CRITERIA)
    values: dict[str, np.ndarray] = {}
    for name in names:
        if name not in table.columns:
            raise CorpusError(f"corpus feature matrix has no column {name!r}")
        column = pd.to_numeric(table[name], errors="coerce").dropna().to_numpy()
        if column.size == 0:
            raise CorpusError(f"corpus column {name!r} has no usable values")
        values[name] = column
    return ReferenceDistributions(
        values=values,
        n_sequences=len(table),
        corpus_id=info.corpus_id,
        is_placeholder=info.is_placeholder,
    )
