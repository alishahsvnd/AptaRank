"""Sanity-check the Band 2/3 descriptors against the real demo targets.

The refinements spec (§6) asks for exactly this before the new descriptors are
relied on in evaluation figures: look at `footprint`, `planarity` and
`elongation` on the actual targets, and confirm they say something rather than
saturating, collapsing, or quietly agreeing with each other.

It answers three questions:

  1. Do the protein-side shape descriptors distinguish the two demo targets?
     A cavity and an interface should not look alike.
  2. Does the radius-of-gyration footprint separate candidates the length proxy
     cannot? If it merely reproduces length, the refinement bought nothing.
  3. Do the two footprint models put candidates in different bands? If they
     agree everywhere the choice is cosmetic; where they differ, that is the
     scientific decision being made.

    python scripts/validate_descriptors.py \
        --targets ~/aptarank-data/cache/targets/*.bundle.json \
        --candidates data/demo_candidates.csv \
        --development-corpus data/corpus/dev_placeholder_corpus.csv
"""

from __future__ import annotations

import argparse
import glob
import json
import sys
from pathlib import Path

import numpy as np

from aptarank.config import load_config
from aptarank.ingest import ingest
from aptarank.tier1 import corpus as corpus_mod
from aptarank.tier1 import features as feature_mod
from aptarank.tier2 import bundle as bundle_mod
from aptarank.tier2 import calibration, modes

FOOTPRINT_MODELS = ("radius_of_gyration", "length")


def describe_targets(paths: list[Path]) -> list[dict]:
    """Protein-side shape descriptors, per prepared target."""
    rows = []
    for path in paths:
        bundle = bundle_mod.load(path)
        summary = bundle_mod.summary(bundle)
        pocket = summary.get("selected_pocket") or {}
        patch = summary.get("patch") or {}
        rows.append(
            {
                "target": f"{summary['identifier']} {summary['chain']}",
                "mode": summary["binding_mode"],
                "n_pockets": summary["n_pockets"],
                "d_pocket_A": pocket.get("d_pocket_A"),
                "d_equiv_A": pocket.get("d_equiv_A"),
                "pocket_elongation": pocket.get("elongation"),
                "pocket_planarity_A": pocket.get("planarity_A"),
                "patch_area_A2": patch.get("patch_area_A2"),
                "patch_elongation": patch.get("elongation"),
                "patch_planarity_A": patch.get("planarity_A"),
                "n_site_residues": patch.get("n_residues"),
                "path": str(path),
            }
        )
    return rows


def candidate_descriptors(cfg, candidates_path: Path) -> "np.ndarray":
    """Fold the candidate set once and return its size descriptors."""
    ingested = ingest(
        candidates_path,
        min_length=int(cfg.get("input.min_length")),
        max_length=int(cfg.get("input.max_length")),
    )
    geometry = corpus_mod.geometry_settings(cfg)
    jobs = [
        feature_mod.FeatureJob(
            candidate_id=row.candidate_id,
            sequence=row.sequence,
            n_ensemble_samples=int(cfg.get("tier1.n_ensemble_samples")),
            n_shuffles=0,
            seed=int(cfg.get("run.seed")),
            a_per_bp_helix=geometry["a_per_bp_helix"],
            a_per_nt_ss=geometry["a_per_nt_ss"],
        )
        for row in ingested.candidates.itertuples()
    ]
    results = feature_mod.compute_batch(
        jobs, workers=cfg.get("tier1.parallel.workers", None)
    )
    table, _failures = feature_mod.results_to_frame(results)
    return table


def summarise(values, label: str, fmt: str = "{:.1f}") -> str:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return f"  {label:<26} (none)"
    parts = [fmt.format(v) for v in (values.min(), np.median(values), values.max())]
    return f"  {label:<26} min {parts[0]:>9}  median {parts[1]:>9}  max {parts[2]:>9}"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--targets", nargs="+", required=True)
    parser.add_argument("--candidates", required=True)
    parser.add_argument("-c", "--config")
    parser.add_argument("--corpus")
    parser.add_argument("--development-corpus")
    parser.add_argument("--set", dest="sets", action="append", default=[])
    parser.add_argument("--json-out")
    args = parser.parse_args(argv)

    overrides: dict = {}
    if args.corpus:
        overrides["corpus"] = {"path": args.corpus}
    if args.development_corpus:
        overrides["corpus"] = {
            "path": args.development_corpus,
            "is_placeholder": True,
            "allow_placeholder": True,
        }
    cfg = load_config(args.config, overrides, args.sets)

    paths = [Path(p) for pattern in args.targets for p in sorted(glob.glob(pattern))]
    if not paths:
        print("no prepared targets matched", file=sys.stderr)
        return 2

    report: dict = {}

    # -- 1. protein side --------------------------------------------------
    print("\n== Target shape descriptors ==\n")
    targets = describe_targets(paths)
    report["targets"] = targets
    for row in targets:
        print(f"  {row['target']}  ({row['mode']}, {row['n_pockets']} cavities detected)")
        if row["d_pocket_A"]:
            print(f"     cavity   d_pocket {row['d_pocket_A']:.1f} A   "
                  f"d_equiv {row['d_equiv_A']:.1f} A   "
                  f"elongation {row['pocket_elongation']:.2f}   "
                  f"planarity {row['pocket_planarity_A']:.1f} A")
        if row["patch_area_A2"]:
            print(f"     patch    area {row['patch_area_A2']:.0f} A^2 over "
                  f"{row['n_site_residues']} residues   "
                  f"elongation {row['patch_elongation']:.2f}   "
                  f"planarity {row['patch_planarity_A']:.1f} A")
    print()

    # -- 2. RNA side ------------------------------------------------------
    print("== Candidate size descriptors ==\n")
    table = candidate_descriptors(cfg, Path(args.candidates))
    scale = float(cfg.get("tier2.geometry.footprint_scale"))
    a_per_nt = float(cfg.get("tier2.geometry.a_per_nt_ss"))
    from aptarank.tier2.geometry import aptamer_footprint_area, footprint_area_from_radius

    rg = table["rg_median_A"].to_numpy(dtype=float)
    lengths = table["length"].to_numpy(dtype=float)
    area_rg = np.array([footprint_area_from_radius(v, scale) for v in rg])
    area_len = np.array([aptamer_footprint_area(v, a_per_nt, scale) for v in lengths])

    print(f"  {len(table)} candidates from {Path(args.candidates).name}")
    print(summarise(lengths, "length (nt)", "{:.0f}"))
    print(summarise(rg, "radius of gyration (A)"))
    print(summarise(area_rg, "footprint, Rg model (A^2)", "{:.0f}"))
    print(summarise(area_len, "footprint, length model (A^2)", "{:.0f}"))

    # The question the refinement exists to answer: does shape add anything
    # beyond length? A correlation of 1 would mean it does not.
    finite = np.isfinite(rg) & np.isfinite(lengths)
    correlation = float(np.corrcoef(rg[finite], lengths[finite])[0, 1])
    spread = [
        float(np.std(rg[finite & (lengths == n)]))
        for n in np.unique(lengths[finite])
        if (finite & (lengths == n)).sum() > 2
    ]
    print(f"\n  Rg vs length correlation      r = {correlation:+.3f}")
    if spread:
        print(f"  Rg spread among equal-length candidates   "
              f"median sd {np.median(spread):.1f} A  max sd {max(spread):.1f} A")
    report["candidates"] = {
        "n": int(len(table)),
        "rg_vs_length_correlation": correlation,
        "rg_sd_within_length_median_A": float(np.median(spread)) if spread else None,
        "length_nt": {"min": float(lengths.min()), "median": float(np.median(lengths)),
                      "max": float(lengths.max())},
        "rg_A": {"min": float(np.nanmin(rg)), "median": float(np.nanmedian(rg)),
                 "max": float(np.nanmax(rg))},
    }

    # -- 3. band comparison on each surface target ------------------------
    surface_targets = [row for row in targets if row["mode"] == "surface"]
    if surface_targets:
        print("\n== Bands, Rg model vs length model ==")
        corpus_table, corpus_info = corpus_mod.build_or_load(cfg)
        bank = calibration.build_or_load(cfg, corpus_table, corpus_info)
        bands = cfg.get("tier2.band_percentiles")
        report["band_comparison"] = []

        for row in surface_targets:
            bundle = bundle_mod.load(Path(row["path"]))
            target = modes.target_measurement(bundle, modes.SURFACE)
            per_model = {}
            for model in FOOTPRINT_MODELS:
                params = dict(modes.parameters(cfg, modes.SURFACE))
                params["footprint_model"] = model
                distribution = calibration.target_distribution(
                    bank, modes.SURFACE, target, params
                )
                column = modes.descriptor_column(modes.SURFACE, params)
                records = calibration.score_candidates(
                    table[column].to_numpy(dtype=float),
                    mode=modes.SURFACE, target=target, params=params,
                    distribution=distribution,
                    moderate=float(bands["moderate"]), strong=float(bands["strong"]),
                )
                per_model[model] = [r["band"] for r in records]

            counts = {
                model: {b: values.count(b) for b in ("strong", "moderate", "weak")}
                for model, values in per_model.items()
            }
            agree = sum(
                1 for a, b in zip(per_model["radius_of_gyration"], per_model["length"])
                if a == b
            )
            print(f"\n  {row['target']}  patch {row['patch_area_A2']:.0f} A^2")
            for model, tally in counts.items():
                print(f"     {model:<20} strong {tally['strong']:>4}  "
                      f"moderate {tally['moderate']:>4}  weak {tally['weak']:>4}")
            print(f"     same band under both  {agree}/{len(table)} "
                  f"({agree / max(len(table), 1):.0%})")
            report["band_comparison"].append(
                {"target": row["target"], "patch_area_A2": row["patch_area_A2"],
                 "counts": counts, "agreement": agree / max(len(table), 1)}
            )

    if args.json_out:
        Path(args.json_out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.json_out).write_text(json.dumps(report, indent=2), encoding="utf-8")
        print(f"\nwrote {args.json_out}")
    print()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
