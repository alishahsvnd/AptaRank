"""Paper figures, drawn from stored experiment results (spec §8.3).

Every figure is regenerable from a results JSON by a single call — no folding,
no re-scoring. That is what makes "every figure is reproducible from a stored
artifact" true rather than aspirational.

Colours are the validated categorical slots (blue, orange, aqua, yellow) for
group comparisons, and a single-hue ordinal ramp for the Tier 2 bands. Bands
deliberately avoid a green/amber/red status palette: "strong" means better
geometric agreement than most shuffled controls, not "good candidate", and a
traffic light would say otherwise.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

# Categorical slots 1-4 (validated adjacent-pair CVD ΔE 9.1, normal 22.9).
GROUP_COLORS = {
    "validated": "#2a78d6",
    "generated": "#eb6834",
    "random": "#1baf7a",
    "shuffled": "#eda100",
}
GROUP_ORDER = ["validated", "generated", "random", "shuffled"]
GROUP_LABELS = {
    "validated": "Validated\naptamers",
    "generated": "FM-generated",
    "random": "Random RNA\n(matched)",
    "shuffled": "Dinucleotide\nshuffles",
}

INK = "#0b0b0b"
SECONDARY = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"


def _style(ax) -> None:
    """Recessive chrome: the data carries the ink."""
    ax.set_facecolor("none")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color("#c3c2b7")
    ax.tick_params(colors=MUTED, labelsize=9, length=3)
    ax.yaxis.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.xaxis.label.set_color(SECONDARY)
    ax.yaxis.label.set_color(SECONDARY)


def e1_violin(results: Mapping[str, Any], path: str | Path) -> Path:
    """Tier 1 score distributions per group, with medians labelled."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    distributions = results["distributions"]
    groups = [g for g in GROUP_ORDER if distributions.get(g)]
    data = [np.asarray(distributions[g], dtype=float) for g in groups]

    fig, ax = plt.subplots(figsize=(6.0, 3.4), dpi=200)
    parts = ax.violinplot(data, showextrema=False, widths=0.75)
    for body, group in zip(parts["bodies"], groups):
        body.set_facecolor(GROUP_COLORS[group])
        body.set_alpha(0.75)
        body.set_edgecolor("none")

    for i, values in enumerate(data, start=1):
        median = float(np.median(values))
        q1, q3 = np.percentile(values, [25, 75])
        ax.vlines(i, q1, q3, color=INK, linewidth=2.0, zorder=3)
        ax.plot(i, median, "o", color="#fcfcfb", markersize=4.5,
                markeredgecolor=INK, markeredgewidth=1.2, zorder=4)
        # Direct labels are the relief for the sub-3:1 contrast of two slots.
        ax.annotate(f"{median:.2f}", (i, median), xytext=(9, -3),
                    textcoords="offset points", fontsize=8.5, color=SECONDARY)

    ax.set_xticks(range(1, len(groups) + 1))
    ax.set_xticklabels([GROUP_LABELS[g] for g in groups], fontsize=8.5, color=SECONDARY)
    ax.set_ylabel("Tier 1 score")
    ax.set_ylim(0, 1)
    _style(ax)

    comparison = results["comparisons"].get("validated_vs_random", {})
    if comparison.get("auc") is not None:
        lo, hi = comparison["auc_ci95"]
        ax.set_title(
            f"Validated vs random: AUC {comparison['auc']:.3f} "
            f"[{lo:.3f}, {hi:.3f}]   n = {comparison['n_a']}/{comparison['n_b']}",
            fontsize=9, color=SECONDARY, loc="left", pad=10,
        )
    fig.text(0.005, 0.005,
             "Out-of-fold: validated sequences scored against a corpus excluding "
             "their own target's family.",
             fontsize=7, color=MUTED)
    return _save(fig, path)


def e2_pass_rates(results: Mapping[str, Any], path: str | Path) -> Path:
    """Shuffled-control pass rate per group, with Wilson intervals."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    groups = [g for g in GROUP_ORDER if g in results["groups"]]
    rates = [results["groups"][g]["pass_rate"] for g in groups]
    errors = np.array([
        [
            results["groups"][g]["pass_rate"] - (results["groups"][g]["ci95"] or [0, 0])[0],
            (results["groups"][g]["ci95"] or [0, 0])[1] - results["groups"][g]["pass_rate"],
        ]
        for g in groups
    ]).T

    fig, ax = plt.subplots(figsize=(5.6, 3.2), dpi=200)
    positions = np.arange(len(groups))
    ax.bar(positions, rates, width=0.6,
           color=[GROUP_COLORS[g] for g in groups], zorder=2)
    ax.errorbar(positions, rates, yerr=errors, fmt="none",
                ecolor=INK, elinewidth=1.2, capsize=3, zorder=3)
    for x, rate, group in zip(positions, rates, groups):
        # Anchor the label above the error bar, not the bar: at a 59% rate the
        # CI cap sits exactly where a bar-anchored label would land.
        top = (results["groups"][group]["ci95"] or [rate, rate])[1]
        ax.annotate(f"{rate:.0%}   n={results['groups'][group]['n']}",
                    (x, top), xytext=(0, 7), textcoords="offset points",
                    ha="center", fontsize=8.5, color=SECONDARY)

    alpha = float(results["alpha"])
    ax.axhline(alpha, color=MUTED, linestyle="--", linewidth=1)
    ax.annotate(f"chance ≈ α = {alpha}", (-0.45, alpha), xytext=(0, 4),
                textcoords="offset points", ha="left", fontsize=8, color=MUTED)

    ax.set_xticks(positions)
    ax.set_xticklabels([GROUP_LABELS[g] for g in groups], fontsize=8.5, color=SECONDARY)
    ax.set_ylabel("Shuffled-control pass rate")
    ax.set_ylim(0, 1.12)
    ax.set_xlim(-0.75, len(groups) - 0.25)   # room for the chance-level label
    _style(ax)
    return _save(fig, path)


def e3_paired(results: Mapping[str, Any], path: str | Path) -> Path:
    """True-target vs decoy control percentile: paired difference + AUROC."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if results.get("status") != "complete":
        raise ValueError(f"E3 did not complete: {results.get('reason')}")

    difference = results["paired_difference"]
    auroc = results["auroc"]

    fig, ax = plt.subplots(figsize=(5.4, 3.0), dpi=200)
    for i, (label, value, ci) in enumerate(
        [
            ("Paired difference\n(true − mean decoy)", difference["mean"], difference["ci95"]),
            ("AUROC − 0.5\n(true vs decoy)",
             (auroc["auc"] - 0.5) if auroc["auc"] is not None else None,
             [c - 0.5 for c in auroc["ci95"]] if auroc.get("ci95") else None),
        ]
    ):
        if value is None:
            continue
        ax.errorbar(
            value, i,
            xerr=[[value - ci[0]], [ci[1] - value]] if ci else None,
            fmt="o", color="#2a78d6", ecolor=INK, elinewidth=1.4,
            capsize=4, markersize=8,
        )
        ax.annotate(f"{value:+.3f}", (value, i), xytext=(0, 12),
                    textcoords="offset points", ha="center",
                    fontsize=9, color=SECONDARY)

    ax.axvline(0.0, color=MUTED, linewidth=1)
    ax.set_yticks([0, 1])
    ax.set_yticklabels(
        ["Paired difference\n(true − mean decoy)", "AUROC − 0.5\n(true vs decoy)"],
        fontsize=8.5, color=SECONDARY,
    )
    ax.set_xlabel("effect size (0 = no target specificity)")
    ax.set_ylim(-0.6, 1.6)
    _style(ax)
    ax.yaxis.grid(False)
    ax.xaxis.grid(True, color=GRID, linewidth=0.8)

    retrieval = results.get("retrieval") or {}
    ax.set_title(
        f"n = {results['n_pairs']} aptamers across {results['n_targets']} targets"
        + (f" · top-1 {retrieval['top1_accuracy']:.0%}, MRR {retrieval['mrr']:.2f}"
           if retrieval.get("mrr") is not None else ""),
        fontsize=9, color=SECONDARY, loc="left", pad=10,
    )
    fig.text(0.005, 0.005,
             "CIs bootstrap the aptamers, not the individual pairs. Decoys are "
             "presumed non-matched, not confirmed non-binders.",
             fontsize=7, color=MUTED)
    return _save(fig, path)


def e4_swappability(results: Mapping[str, Any], path: str | Path) -> Path:
    """How much the Tier 2 annotation moves when the target is swapped."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    comparisons = results.get("comparisons") or []
    if not comparisons:
        raise ValueError("E4 has no target pairs to compare")

    labels = [f"{c['targets'][0]} vs {c['targets'][1]}" for c in comparisons]
    changed = [c["band_changed_fraction"] for c in comparisons]
    rhos = [c["spearman_rho"] for c in comparisons]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 3.0), dpi=200)
    positions = np.arange(len(labels))

    ax1.barh(positions, changed, height=0.55, color="#2a78d6", zorder=2)
    for y, value in zip(positions, changed):
        ax1.annotate(f"{value:.0%}", (value, y), xytext=(5, -3),
                     textcoords="offset points", fontsize=8.5, color=SECONDARY)
    ax1.set_xlim(0, 1.1)
    ax1.set_xlabel("candidates whose band changed")

    ax2.barh(positions, rhos, height=0.55, color="#eb6834", zorder=2)
    for y, value in zip(positions, rhos):
        ax2.annotate(f"{value:+.2f}", (value, y), xytext=(5 if value >= 0 else -28, -3),
                     textcoords="offset points", fontsize=8.5, color=SECONDARY)
    ax2.axvline(0, color=MUTED, linewidth=1)
    ax2.set_xlim(-1.1, 1.1)
    ax2.set_xlabel("Spearman ρ of control percentiles")

    for ax in (ax1, ax2):
        ax.set_yticks(positions)
        ax.set_yticklabels(labels, fontsize=8.5, color=SECONDARY)
        _style(ax)
        ax.yaxis.grid(False)
        ax.xaxis.grid(True, color=GRID, linewidth=0.8)

    fig.text(0.005, 0.005,
             "Similar annotations are expected when two cavities have similar "
             "dimensions; compare d_pocket before reading this as a null.",
             fontsize=7, color=MUTED)
    return _save(fig, path)


def _save(fig, path: str | Path) -> Path:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout(rect=(0, 0.02, 1, 1))
    fig.savefig(out, bbox_inches="tight", facecolor="none", transparent=True)
    # PDF alongside PNG: the paper needs vector, the dashboard needs raster.
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight", transparent=True)
    import matplotlib.pyplot as plt

    plt.close(fig)
    return out


def build_all(
    results_path: str | Path, out_dir: str | Path, allow_development: bool = False
) -> list[Path]:
    """Regenerate every figure whose experiment is present in a results file.

    Refuses development-grade results by default. A figure drawn from a
    placeholder corpus or a synthetic target bundle looks exactly like a real
    one, and the only thing standing between it and a paper is this check.
    """
    results = json.loads(Path(results_path).read_text(encoding="utf-8"))
    reasons = development_reasons(results)
    if reasons and not allow_development:
        raise ValueError(
            "refusing to draw paper figures from development-grade results: "
            + "; ".join(reasons)
            + ". Pass allow_development=True (CLI: --allow-development) to draw "
              "them anyway; they will be watermarked."
        )
    out = Path(out_dir)
    written: list[Path] = []
    builders = {
        "E1": (e1_violin, "e1_tier1_distributions.png"),
        "E2": (e2_pass_rates, "e2_shuffle_pass_rates.png"),
        "E3": (e3_paired, "e3_matched_vs_decoy.png"),
        "E4": (e4_swappability, "e4_target_swappability.png"),
    }
    for key, (builder, filename) in builders.items():
        payload = results.get(key)
        if not payload:
            continue
        try:
            path = builder(payload, out / filename)
        except (ValueError, KeyError) as exc:
            print(f"  skipped {key}: {exc}")
            continue
        if reasons:
            _watermark(path, reasons)
        written.append(path)
    return written


def development_reasons(results: Mapping[str, Any]) -> list[str]:
    """Every reason these results must not appear in the paper."""
    reasons = []
    corpus = results.get("corpus") or {}
    if corpus.get("is_placeholder"):
        reasons.append("scored against a placeholder corpus")
    for key, payload in results.items():
        if not isinstance(payload, Mapping):
            continue
        for target in (payload.get("targets") or {}).values():
            if isinstance(target, Mapping) and target.get("synthetic"):
                reasons.append(f"{key} used a synthetic target bundle")
                break
    return reasons


def _watermark(path: Path, reasons: Sequence[str]) -> None:
    """Stamp a development figure so it cannot be mistaken for a result."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.image as mpimg
    import matplotlib.pyplot as plt

    image = mpimg.imread(path)
    height, width = image.shape[0], image.shape[1]
    fig, ax = plt.subplots(figsize=(width / 200, height / 200), dpi=200)
    ax.imshow(image)
    ax.axis("off")
    ax.text(
        0.5, 0.5, "DEVELOPMENT DATA\nNOT A RESULT",
        transform=ax.transAxes, ha="center", va="center",
        fontsize=26, color="#d03b3b", alpha=0.32, rotation=24, fontweight="bold",
    )
    ax.text(
        0.5, 0.995, "; ".join(reasons), transform=ax.transAxes,
        ha="center", va="top", fontsize=6, color="#d03b3b",
    )
    fig.savefig(path, bbox_inches="tight", pad_inches=0, dpi=200)
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0)
    plt.close(fig)
