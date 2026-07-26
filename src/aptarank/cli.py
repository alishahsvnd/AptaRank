"""Command-line interface.

    aptarank run CANDIDATES            score a candidate file, write an artifact
    aptarank corpus build              precompute and cache corpus features
    aptarank artifact summary PATH     human-readable summary of a run artifact
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .config import load_config
from .errors import AptaRankError


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="aptarank",
        description="Two-tier ranking of generated RNA aptamer candidates. "
                    "Tier 2 output is geometric plausibility, not binding prediction.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    def common(p: argparse.ArgumentParser) -> None:
        p.add_argument("-c", "--config", help="YAML config layered over configs/default.yaml")
        p.add_argument("--corpus", help="path to the validated-aptamer reference corpus CSV")
        p.add_argument(
            "--development-corpus",
            help="use a synthetic placeholder corpus; the run is stamped "
                 "publication_eligible: false and must not back paper claims",
        )
        p.add_argument(
            "--set", dest="sets", action="append", default=[], metavar="KEY=VALUE",
            help="override any config key, e.g. --set tier1.shuffle.n_shuffles=99",
        )

    run_p = sub.add_parser("run", help="score a candidate file")
    run_p.add_argument("candidates", help=".txt / .csv / .fasta file of candidate sequences")
    run_p.add_argument("-o", "--output-dir", help="directory for the run artifact")
    run_p.add_argument("--artifact-path", help="exact path to write the artifact to")
    run_p.add_argument("--fast", action="store_true", help="skip shuffled controls and sampling")
    run_p.add_argument("--target", help="target PDB ID; enables Tier 2")
    run_p.add_argument("--target-bundle", help="path to a prepared target bundle JSON")
    run_p.add_argument(
        "--progress-format", choices=("human", "jsonl"), default="human",
        help="human: carriage-return counters. jsonl: machine-readable events.",
    )
    run_p.add_argument(
        "--progress-file",
        help="append JSONL progress events here (readable while the run is live)",
    )
    run_p.add_argument("--job-id", help="identifier stamped into progress events")
    common(run_p)

    val_p = sub.add_parser(
        "validate-input", help="check a candidate file without scoring anything"
    )
    val_p.add_argument("candidates")
    val_p.add_argument("--json-out", help="write the validation summary here as JSON")
    common(val_p)

    corpus_p = sub.add_parser("corpus", help="reference corpus utilities")
    corpus_sub = corpus_p.add_subparsers(dest="corpus_command", required=True)
    build_p = corpus_sub.add_parser("build", help="precompute and cache corpus features")
    build_p.add_argument("--force", action="store_true", help="rebuild even if cached")
    common(build_p)

    target_p = sub.add_parser("target", help="target bundle utilities (Linux only)")
    target_sub = target_p.add_subparsers(dest="target_command", required=True)
    tbuild_p = target_sub.add_parser(
        "build", help="fetch a structure, run fpocket, write a target bundle"
    )
    tbuild_p.add_argument("--pdb-id", help="4-character PDB ID, e.g. 3SPU")
    tbuild_p.add_argument("--chain", help="chain identifier; default is the first protein chain")
    tbuild_p.add_argument("-o", "--output-dir", dest="bundle_out", help="where to write the bundle")
    common(tbuild_p)

    tshow_p = target_sub.add_parser("show", help="summarise a target bundle")
    tshow_p.add_argument("path")

    eval_p = sub.add_parser("evaluate", help="evaluation experiments E1-E5 (§8)")
    eval_sub = eval_p.add_subparsers(dest="evaluate_command", required=True)

    erun_p = eval_sub.add_parser("run", help="run experiments and write results JSON")
    erun_p.add_argument(
        "--experiments", default="E1,E2",
        help="comma-separated subset of E1,E2,E3,E4,E5 (default E1,E2)",
    )
    erun_p.add_argument("--generated", help="CSV/TXT of foundation-model sequences")
    erun_p.add_argument("--artifacts", nargs="*", default=[],
                        help="run artifacts for E4/E5, one per target")
    erun_p.add_argument("--folds", type=int, default=5, help="E1 out-of-fold splits")
    erun_p.add_argument("--n-per-group", type=int, default=200, help="E2 sample per group")
    erun_p.add_argument("-o", "--output", default="evaluation/results.json")
    common(erun_p)

    efig_p = eval_sub.add_parser("figures", help="regenerate paper figures from results")
    efig_p.add_argument("results", help="results JSON from `evaluate run`")
    efig_p.add_argument("-o", "--output-dir", default="evaluation/figures")
    efig_p.add_argument(
        "--allow-development", action="store_true",
        help="draw figures from placeholder/synthetic results; they are watermarked",
    )

    art_p = sub.add_parser("artifact", help="run artifact utilities")
    art_sub = art_p.add_subparsers(dest="artifact_command", required=True)
    sum_p = art_sub.add_parser("summary", help="print a summary of a run artifact")
    sum_p.add_argument("path")
    sum_p.add_argument("-n", "--top", type=int, default=10, help="rows to show")

    return parser


def _config_from_args(args: argparse.Namespace):
    overrides: dict = {}
    if getattr(args, "corpus", None):
        overrides.setdefault("corpus", {})["path"] = args.corpus
    if getattr(args, "development_corpus", None):
        overrides.setdefault("corpus", {}).update(
            {"path": args.development_corpus, "is_placeholder": True, "allow_placeholder": True}
        )
    if getattr(args, "fast", False):
        overrides.setdefault("run", {})["mode"] = "fast"
    if getattr(args, "output_dir", None):
        overrides.setdefault("output", {})["dir"] = args.output_dir
    if getattr(args, "target", None):
        overrides.setdefault("tier2", {}).update(
            {"enabled": True, "target": {"pdb_id": args.target}}
        )
    if getattr(args, "target_bundle", None):
        overrides.setdefault("tier2", {}).update(
            {"enabled": True, "bundle_path": args.target_bundle}
        )
    return load_config(getattr(args, "config", None), overrides, getattr(args, "sets", []))


def cmd_run(args: argparse.Namespace) -> int:
    from .pipeline import run_pipeline
    from .progress import ProgressReporter

    cfg = _config_from_args(args)
    reporter = ProgressReporter(
        fmt=args.progress_format, path=args.progress_file, job_id=args.job_id
    )
    try:
        out = run_pipeline(
            cfg, args.candidates, reporter=reporter, artifact_path=args.artifact_path
        )
    except AptaRankError as exc:
        # A typed error carries a message a non-expert can act on; anything
        # else is a bug and keeps its class name so we can find it in the log.
        reporter.run_failed(type(exc).__name__, str(exc))
        reporter.close()
        raise
    except Exception as exc:  # noqa: BLE001 - must still reach the progress file
        reporter.run_failed("UnexpectedError", f"{type(exc).__name__}: {exc}")
        reporter.close()
        raise
    reporter.close()

    print(f"\nRanked {len(out.tier1.table)} candidates  ->  {out.path}")
    print(_format_top(out.artifact, 10))
    return 0


def cmd_validate_input(args: argparse.Namespace) -> int:
    """Check a candidate file with exactly the ingest rules a run would apply."""
    import json as _json

    from .ingest import ingest

    cfg = _config_from_args(args)
    try:
        result = ingest(
            args.candidates,
            min_length=int(cfg.get("input.min_length")),
            max_length=int(cfg.get("input.max_length")),
        )
        summary = {"ok": True, **result.summary()}
        lengths = result.candidates["length"]
        summary["length_range"] = [int(lengths.min()), int(lengths.max())]
    except AptaRankError as exc:
        summary = {"ok": False, "error": str(exc), "error_type": type(exc).__name__}

    if args.json_out:
        Path(args.json_out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.json_out).write_text(_json.dumps(summary, indent=2), encoding="utf-8")
    print(_json.dumps(summary, indent=2))
    return 0 if summary["ok"] else 1


def cmd_corpus_build(args: argparse.Namespace) -> int:
    from .tier1 import corpus as corpus_mod

    cfg = _config_from_args(args)
    table, info = corpus_mod.build_or_load(cfg, force=args.force)
    print(
        f"\nCorpus features cached: {info.n_sequences} sequences "
        f"({info.n_dropped} dropped: {info.dropped_reasons or 'none'})\n"
        f"  corpus_id           {info.corpus_id}\n"
        f"  placeholder         {info.is_placeholder}\n"
        f"  publication_eligible {info.publication_eligible}"
    )
    return 0


def cmd_target_build(args: argparse.Namespace) -> int:
    from .tier2.build import build_target_bundle

    cfg = _config_from_args(args)
    bundle, path = build_target_bundle(
        cfg, pdb_id=args.pdb_id, chain_id=args.chain, out_dir=args.bundle_out
    )
    print(f"\nTarget bundle written: {path}")
    print(_format_target(bundle))
    return 0


def cmd_target_show(args: argparse.Namespace) -> int:
    from .tier2 import bundle as bundle_mod

    print(_format_target(bundle_mod.load(args.path)))
    return 0


def _format_target(bundle: dict) -> str:
    from .tier2 import bundle as bundle_mod

    s = bundle_mod.summary(bundle)
    pocket = s["selected_pocket"]
    lines = [
        f"\n  {s['pdb_id']} chain {s['chain']}   bundle {s['bundle_id'][:12]}",
        f"  name                 {s['name']}",
        f"  pockets detected     {s['n_pockets']}",
        f"  selection method     {s['pocket_selection']}",
        f"  selected pocket      #{pocket['index']}  "
        f"volume {pocket['volume_A3']:.1f} A^3  fpocket score {pocket['fpocket_score']:.3f}",
        f"  d_pocket (envelope)  {pocket['d_pocket_A']:.2f} A"
        f"   (centres {pocket['d_pocket_centres_A']:.2f} A, equiv sphere {pocket['d_equiv_A']:.2f} A)",
        f"  alpha spheres        {pocket['n_alpha_spheres']}"
        f"{'   [SHAPE WARNING: oddly shaped cavity]' if pocket['shape_warning'] else ''}",
        f"  retained hetero      {sorted(set(s['retained_hetero'])) or 'none'}",
        f"  electrostatics       {s['electrostatics_status']}",
    ]
    for warning in s["selection_warnings"]:
        lines.append(f"  !! {warning}")
    return "\n".join(lines)


def cmd_evaluate_run(args: argparse.Namespace) -> int:
    from .artifacts import read_artifact
    from .evaluation import build_groups, experiments
    from .tier1 import corpus as corpus_mod

    cfg = _config_from_args(args)
    wanted = {e.strip().upper() for e in args.experiments.split(",") if e.strip()}
    corpus_table, corpus_info = corpus_mod.build_or_load(cfg)
    if corpus_info.is_placeholder:
        print(
            "\n  !! PLACEHOLDER CORPUS — these results are a smoke test of the "
            "evaluation machinery,\n     not scientific findings.\n",
            file=sys.stderr,
        )

    groups = build_groups(
        corpus_table,
        seed=int(cfg.get("run.seed")),
        generated_path=args.generated,
        min_length=int(cfg.get("input.min_length")),
        max_length=int(cfg.get("input.max_length")),
    )
    print(f"\ngroups: {groups.summary()['sizes']}")

    def progress(stage: str, i: int, n: int) -> None:
        if n and (i == n or i % max(1, n // 10) == 0):
            print(f"\r  {stage:<26} {i}/{n}", end="" if i < n else "\n", flush=True)

    results: dict = {
        "config": cfg.as_dict(),
        "corpus": corpus_info.to_dict(),
        "groups": groups.summary(),
    }

    if "E1" in wanted:
        results["E1"] = experiments.e1_discrimination(
            cfg, corpus_table, corpus_info, groups, n_folds=args.folds, progress=progress
        )
    if "E2" in wanted:
        results["E2"] = experiments.e2_shuffled_controls(
            cfg, corpus_table, corpus_info, groups,
            n_per_group=args.n_per_group, progress=progress,
        )

    artifacts = [read_artifact(p) for p in args.artifacts]
    if "E4" in wanted:
        results["E4"] = experiments.e4_target_swappability(artifacts)
    if "E5" in wanted:
        results["E5"] = experiments.e5_tier_independence(artifacts)
    if "E3" in wanted:
        print("  E3 needs per-target bundles and labelled aptamer/target pairs; "
              "see evaluation/README.md", file=sys.stderr)

    path = experiments.write_results(results, args.output)
    print(f"\nwrote {path}")
    for key in sorted(k for k in results if k.startswith("E")):
        print(f"  {key}: {results[key].get('title', results[key].get('status'))}")
    return 0


def cmd_evaluate_figures(args: argparse.Namespace) -> int:
    from .evaluation import figures

    try:
        written = figures.build_all(
            args.results, args.output_dir, allow_development=args.allow_development
        )
    except ValueError as exc:
        print(f"\nerror: {exc}", file=sys.stderr)
        return 1
    for path in written:
        print(f"  wrote {path} (+ .pdf)")
    return 0 if written else 1


def cmd_artifact_summary(args: argparse.Namespace) -> int:
    from .artifacts import read_artifact

    artifact = read_artifact(args.path)
    diag = artifact["diagnostics"]
    print(
        f"\n{artifact['run_id']}   {artifact['created_utc']}\n"
        f"  mode                 {artifact['run_mode']} "
        f"(publication_eligible={artifact['publication_eligible']})\n"
        f"  input                {artifact['input']['filename']} "
        f"({artifact['input']['n_valid']} valid / {artifact['input']['n_submitted']} submitted, "
        f"{artifact['input']['n_rejected']} rejected)\n"
        f"  corpus               {artifact['corpus']['corpus_id']} "
        f"(n={artifact['corpus']['n_sequences']}, placeholder={artifact['corpus']['is_placeholder']})\n"
        f"  composite method     {diag['composite_method']}\n"
        f"  tier1 runtime        {diag['runtime_seconds']['tier1']} s\n"
        f"  spearman t1/t2       {diag['spearman_tier1_tier2']}"
    )
    print(_format_top(artifact, args.top))
    return 0


def _format_top(artifact: dict, n: int) -> str:
    header = f"\n  {'rank':>4}  {'candidate':<10} {'tier1':>6}  {'band':<14} {'shuffle':<8}  sequence"
    lines = [header, "  " + "-" * (len(header) - 3)]
    for rec in artifact["candidates"][:n]:
        shuffle = rec.get("shuffle") or {}
        passed = shuffle.get("pass")
        mark = "pass" if passed else ("fail" if passed is False else "-")
        seq = rec["sequence"]
        lines.append(
            f"  {rec['rank']:>4}  {rec['candidate_id']:<10} {rec['tier1_score']:.3f}  "
            f"{(rec['tier2'] or {}).get('band', 'not_evaluated'):<14} {mark:<8}  "
            f"{seq[:34]}{'…' if len(seq) > 34 else ''}"
        )
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    handlers = {
        ("run", None): cmd_run,
        ("validate-input", None): cmd_validate_input,
        ("corpus", "build"): cmd_corpus_build,
        ("target", "build"): cmd_target_build,
        ("target", "show"): cmd_target_show,
        ("evaluate", "run"): cmd_evaluate_run,
        ("evaluate", "figures"): cmd_evaluate_figures,
        ("artifact", "summary"): cmd_artifact_summary,
    }
    key = (
        args.command,
        getattr(args, "corpus_command", None)
        or getattr(args, "target_command", None)
        or getattr(args, "evaluate_command", None)
        or getattr(args, "artifact_command", None),
    )
    handler = handlers.get(key)
    if handler is None:
        print(f"unknown command: {key}", file=sys.stderr)
        return 2
    try:
        return handler(args)
    except AptaRankError as exc:
        print(f"\nerror: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
