"""Headless end-to-end orchestration.

The pipeline is fully usable without the dashboard: it consumes a candidate
file and a config, and produces one run artifact. Nothing downstream — the
dashboard, the evaluation scripts, the paper figures — computes anything of
its own; the dashboard launches *this* code as a subprocess and reads what it
writes.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .artifacts import build_artifact, write_artifact
from .artifacts.explanation import explain
from .artifacts.rendering import render_diagrams
from .config import Config
from .errors import TargetError
from .ingest import ingest
from .progress import ProgressReporter
from .tier1 import corpus as corpus_mod
from .tier1 import service as tier1_service


@dataclass
class RunOutput:
    artifact: dict[str, Any]
    path: Path | None
    tier1: tier1_service.Tier1Result


def prepare_target(cfg: Config, report: ProgressReporter) -> dict[str, Any]:
    """Get the target evidence this run needs, building it if it does not exist.

    Preparation is server-side now (refinements §3): the user supplies an
    identifier, a chain and a binding mode, and the heavy, biology-literate step
    happens here rather than on their machine. It is cached on everything that
    changes what was measured, so a repeat run is free while a changed chain or
    residue list rebuilds rather than silently reusing the old measurement.
    """
    from .tier2 import build as build_mod
    from .tier2 import bundle as bundle_mod

    explicit = cfg.get("tier2.bundle_path", None)
    if explicit:
        return bundle_mod.load(explicit)

    identifier = cfg.get("tier2.target.id", None)
    if not identifier:
        raise TargetError(
            "Tier 2 is enabled but no target is configured. Set tier2.target.id "
            "and tier2.target.source, or point tier2.bundle_path at a prepared "
            "target file."
        )

    mode = cfg.get("tier2.binding_mode")
    chain = cfg.get("tier2.target.chain", None)
    signature = build_mod.target_signature(cfg)
    directory = Path(cfg.get("tier2.bundle_dir"))

    for existing in sorted(directory.glob(f"{identifier.upper()}_*.bundle.json")):
        try:
            candidate = bundle_mod.load(existing)
        except TargetError:
            continue   # a bundle that fails its own integrity check is not a hit
        if (
            candidate.get("target_signature") == signature
            and bundle_mod.binding_mode(candidate) == mode
        ):
            return candidate

    report.warning(
        f"Preparing {identifier} chain {chain or '(first protein chain)'} for "
        f"{mode} mode — fetching the structure and measuring it. This happens "
        f"once per target.",
        code="target_building",
    )
    bundle, _path = build_mod.build_target_bundle(
        cfg,
        progress=lambda message: report.warning(message, code="target_progress"),
    )
    return bundle


def run_pipeline(
    cfg: Config,
    candidates_path: str | Path,
    write: bool = True,
    reporter: ProgressReporter | None = None,
    artifact_path: str | Path | None = None,
) -> RunOutput:
    """Score a candidate file end to end.

    `artifact_path` lets a caller fix the output location before launch, so a
    UI never has to scrape a path out of stdout to find the result.
    """
    report = reporter or ProgressReporter()
    report.run_started(candidates=str(candidates_path))

    report.stage_started("ingest")
    ingested = ingest(
        candidates_path,
        min_length=int(cfg.get("input.min_length")),
        max_length=int(cfg.get("input.max_length")),
    )
    report.stage_completed("ingest", n_valid=ingested.n_valid,
                           n_rejected=ingested.n_rejected)
    if ingested.n_rejected:
        report.warning(
            f"{ingested.n_rejected} of {ingested.n_submitted} submitted rows were "
            f"rejected and are listed in the results.",
            code="rows_rejected",
        )

    report.stage_started("corpus")
    corpus_table, corpus_info = corpus_mod.build_or_load(
        cfg, progress=report.callback("corpus", "sequences")
    )
    report.stage_completed("corpus", n_sequences=corpus_info.n_sequences)
    if corpus_info.is_placeholder:
        report.warning(
            "Placeholder reference library: scores are calibrated against "
            "synthetic data. This run cannot support any published claim.",
            code="placeholder_corpus",
        )
        print(
            "\n  !! PLACEHOLDER CORPUS — scores are calibrated against synthetic "
            "reference data.\n     This run is marked publication_eligible: false "
            "and must not be used for paper claims.\n",
            file=sys.stderr,
        )
    refs = corpus_mod.reference_distributions(
        corpus_table, corpus_info, cfg.active_criteria()
    )

    report.stage_started("tier1", total=ingested.n_valid)
    tier1 = tier1_service.run(
        cfg, ingested, refs, corpus_info, progress=report.callback("tier1", "candidates")
    )
    report.stage_completed("tier1", n_scored=tier1.n_scored,
                           n_failed=len(tier1.failures))

    tier2_payload: dict[str, Any] | None = None
    if cfg.get("tier2.enabled", False):
        from .tier2 import service as tier2_service  # imported lazily: heavy deps

        report.stage_started("target")
        bundle = prepare_target(cfg, report)
        report.stage_completed(
            "target",
            binding_mode=bundle.get("binding_mode"),
            bundle_id=bundle.get("bundle_id"),
        )

        report.stage_started("bank")
        tier2_payload = tier2_service.run(
            cfg, tier1, corpus_table, corpus_info,
            progress=report.callback("bank", "controls"),
            bundle=bundle,
        )
        report.stage_completed("bank", n_evaluated=tier2_payload["n_evaluated"])
        target_summary = tier2_payload["target"]
        for warning in (target_summary.get("selection_warnings") or []):
            report.warning(warning, code="pocket_selection")
        for warning in (target_summary.get("preparation_warnings") or []):
            report.warning(warning, code="target_preparation")
        if target_summary.get("structure_kind") == "predicted":
            report.warning(
                "Predicted structure: this target's geometry comes from an "
                "AlphaFold model, not an experiment."
                + (
                    " A predicted model may not show an interface that only forms "
                    "when a binding partner is present, so surface-mode results "
                    "here need extra caution."
                    if tier2_payload.get("binding_mode") == "surface" else ""
                ),
                code="predicted_structure",
            )
        if target_summary.get("synthetic"):
            report.warning(
                "Synthetic target evidence: the cavity was fabricated, not "
                "detected from a real structure.",
                code="synthetic_target",
            )

    report.stage_started("diagrams")
    diagrams = render_diagrams(
        tier1.table[["candidate_id", "sequence", "dot_bracket"]].to_dict("records"),
        n=int(cfg.get("output.n_diagrams")),
    ) if cfg.get("output.embed_svg", True) else {}
    report.stage_completed("diagrams", n_diagrams=len(diagrams))

    report.stage_started("artifact")
    artifact = build_artifact(
        cfg, ingested, tier1, corpus_info, tier2=tier2_payload, diagrams=diagrams
    )

    # Explanations are generated from the assembled record, so every number
    # they quote is by construction the number that was stored.
    for record in artifact["candidates"]:
        rendered = explain(record)
        record["explanation"] = rendered["text"]
        record["evidence_chips"] = rendered["chips"]
        record["rules_fired"] = rendered["rules_fired"]

    path: Path | None = None
    if write:
        if artifact_path:
            path = Path(artifact_path)
            path.parent.mkdir(parents=True, exist_ok=True)
            import json

            path.write_text(json.dumps(artifact, indent=2, allow_nan=False), encoding="utf-8")
        else:
            path = write_artifact(artifact, cfg.get("output.dir"))
    report.stage_completed("artifact", path=str(path) if path else None)
    report.run_completed(
        artifact_path=path,
        publication_eligible=artifact["publication_eligible"],
        n_candidates=len(artifact["candidates"]),
    )

    return RunOutput(artifact=artifact, path=path, tier1=tier1)
