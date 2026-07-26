"""Run artifact assembly, writing and reading (spec §7.2).

One JSON file per run, self-describing and self-contained: structure diagrams
are embedded as SVG text rather than referenced by path, so an artifact can be
moved or attached to a paper submission without losing anything.

Any value that is unavailable is `null`. NaN never reaches the file — it is not
valid JSON and silently becomes a parse error or a string in downstream tools.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import pandas as pd

from .. import ARTIFACT_SCHEMA_VERSION, TIER2_CAVEAT
from ..config import CRITERIA, Config, REPO_ROOT
from ..errors import ArtifactError
from ..ingest import IngestResult
from ..provenance import git_state, sha256_file, tool_versions
from ..tier1.corpus import CorpusInfo
from ..tier1.service import Tier1Result

ELEMENT_FIELDS = (
    "n_hairpins", "n_interior", "n_multiloop", "n_stems", "stem_fraction",
    "longest_stem_bp", "max_loop_nt", "total_unpaired",
    "loop_nt_median", "loop_nt_p90", "loop_nt_iqr", "n_ensemble_samples",
)


def build_artifact(
    cfg: Config,
    ingested: IngestResult,
    tier1: Tier1Result,
    corpus_info: CorpusInfo,
    *,
    run_id: str | None = None,
    tier2: Mapping[str, Any] | None = None,
    diagrams: Mapping[str, str] | None = None,
    explanations: Mapping[str, str] | None = None,
    extra_diagnostics: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Assemble the complete run record. Pure function of its inputs."""
    created = datetime.now(timezone.utc)
    input_sha = sha256_file(ingested.filename) if ingested.filename else None
    run_id = run_id or _make_run_id(created, input_sha)

    tier2 = dict(tier2 or {})
    diagrams = diagrams or {}
    explanations = explanations or {}

    # A run is only publishable if BOTH its reference distributions and its
    # target evidence are real. A synthetic cavity invalidates every Tier 2
    # number just as surely as a placeholder corpus invalidates Tier 1.
    synthetic_target = bool((tier2.get("target") or {}).get("synthetic"))
    development = corpus_info.is_placeholder or synthetic_target

    artifact: dict[str, Any] = {
        "artifact_schema_version": ARTIFACT_SCHEMA_VERSION,
        "run_id": run_id,
        "created_utc": created.isoformat().replace("+00:00", "Z"),
        "run_mode": "development" if development else "standard",
        "publication_eligible": not development,
        "development_reasons": [
            reason for reason, active in (
                ("placeholder_corpus", corpus_info.is_placeholder),
                ("synthetic_target_bundle", synthetic_target),
            ) if active
        ],
        "caveat": TIER2_CAVEAT,
        "config": cfg.as_dict(),
        "config_sources": cfg.sources,
        "scoring_signature": cfg.scoring_signature(),
        "versions": tool_versions(),
        "git": git_state(REPO_ROOT),
        "corpus": corpus_info.to_dict(),
        "input": {**ingested.summary(), "sha256": input_sha},
        "target": tier2.get("target"),
        "tier2_thresholds": tier2.get("thresholds"),
        "diagnostics": {
            # §6.4 required diagnostic — computed every run, never deferred to
            # evaluation week, because a high correlation would mean Tier 2 is
            # largely restating Tier 1 and would change what the paper claims.
            "spearman_tier1_tier2": tier2.get("spearman"),
            "runtime_seconds": {
                "tier1": round(tier1.runtime_seconds, 3),
                "tier2": tier2.get("runtime_seconds"),
            },
            "n_scored": tier1.n_scored,
            "n_folding_failures": len(tier1.failures),
            "folding_failures": tier1.failures,
            "composite_method": tier1.composite_method,
            **(dict(extra_diagnostics) if extra_diagnostics else {}),
        },
        "candidates": [],
    }

    per_candidate_tier2 = tier2.get("candidates", {})
    # "not_run" (Tier 2 was disabled for this batch) and "not_evaluated" (this
    # candidate fell below the Tier 2 cut) are different facts, and neither is
    # "weak" — displaying a missing measurement as a poor one fabricates
    # evidence.
    default_status = "not_evaluated" if tier2 else "not_run"
    for row in tier1.table.to_dict(orient="records"):
        cid = row["candidate_id"]
        artifact["candidates"].append(
            _candidate_record(
                row,
                tier2=per_candidate_tier2.get(cid),
                default_status=default_status,
                svg=diagrams.get(cid),
                explanation=explanations.get(cid),
            )
        )

    return _jsonable(artifact)


def _candidate_record(
    row: Mapping[str, Any],
    tier2: Mapping[str, Any] | None,
    svg: str | None,
    explanation: str | None,
    default_status: str = "not_run",
) -> dict[str, Any]:
    criteria = {
        name: {"value": row.get(name), "score": row.get(f"score__{name}")}
        for name in CRITERIA
        if f"score__{name}" in row
    }
    shuffle = (
        {
            "pass": row.get("shuffle_pass"),
            "percentile": row.get("shuffle_percentile"),
            "p_value": row.get("shuffle_p_value"),
            "margin": row.get("shuffle_margin"),
            "n_shuffles": row.get("shuffle_n"),
            "complete": row.get("shuffle_complete"),
            "n_unique": row.get("shuffle_n_unique"),
            "n_identical_to_real": row.get("shuffle_n_identical_to_real"),
            "median_score": row.get("shuffle_median_score"),
            "structural_subscore": row.get("structural_subscore"),
        }
        if row.get("shuffle_n")
        else None
    )

    return {
        "candidate_id": row["candidate_id"],
        "sequence": row["sequence"],
        "length": row.get("length"),
        "duplicate_count": row.get("duplicate_count"),
        "rank": row.get("rank"),
        "rank_min": row.get("rank_min"),
        "rank_is_tied": row.get("rank_is_tied"),
        "tier1_score": row.get("tier1_score"),
        "batch_rank_fraction": row.get("batch_rank_fraction"),
        "criteria": criteria,
        "structure": {
            "dot_bracket": row.get("dot_bracket"),
            "element_string": row.get("element_string"),
            "svg": svg,
        },
        "elements": {k: row.get(k) for k in ELEMENT_FIELDS if k in row},
        "shuffle": shuffle,
        "tier2": dict(tier2) if tier2 else {"status": default_status, "band": "not_evaluated"},
        "explanation": explanation,
    }


def write_artifact(artifact: Mapping[str, Any], directory: str | Path) -> Path:
    """Write `<directory>/<run_id>.json` atomically."""
    out_dir = Path(directory)
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"{artifact['run_id']}.json"
    tmp = path.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(artifact, indent=2, allow_nan=False), encoding="utf-8")
    tmp.replace(path)
    return path


def read_artifact(path: str | Path) -> dict[str, Any]:
    p = Path(path)
    if not p.exists():
        raise ArtifactError(f"run artifact not found: {p}")
    with p.open("r", encoding="utf-8") as fh:
        data = json.load(fh)
    if "artifact_schema_version" not in data:
        raise ArtifactError(f"{p} is not an AptaRank run artifact")
    return data


def _make_run_id(created: datetime, input_sha: str | None) -> str:
    stamp = created.strftime("%Y%m%dT%H%M%SZ")
    suffix = (input_sha or "nohash")[:8]
    return f"run_{stamp}_{suffix}"


def _jsonable(value: Any) -> Any:
    """Recursively convert to JSON-safe types; NaN/inf/NaT become null."""
    if isinstance(value, Mapping):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_jsonable(v) for v in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        f = float(value)
        return None if (math.isnan(f) or math.isinf(f)) else f
    if isinstance(value, (np.bool_,)):
        return bool(value)
    if isinstance(value, np.ndarray):
        return [_jsonable(v) for v in value.tolist()]
    if value is pd.NA or value is pd.NaT:
        return None
    if isinstance(value, pd.Timestamp):
        return value.isoformat()
    return value
