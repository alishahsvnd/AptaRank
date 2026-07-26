"""Panel renderers. Each one reads the run artifact and computes nothing.

Keeping every calculation in the pipeline is what lets the paper's figures be
generated from the same artifacts the live demo displays.
"""

from __future__ import annotations

import base64
from typing import Any, Mapping, Sequence

import altair as alt
import pandas as pd
import streamlit as st

from .theme import BAND_LABEL, band_scale, chart_config

CRITERION_LABELS = {
    "mfe_norm": "Structure stability",
    "ensemble_defect": "Fold definition",
    "positional_entropy_mean": "Fold certainty",
    "stem_fraction": "Structural composition",
    "gc_fraction": "Sequence composition",
}


# -- (a) run configuration ----------------------------------------------


def run_configuration(artifact: Mapping[str, Any]) -> None:
    cfg = artifact["config"]
    inp = artifact["input"]
    target = artifact.get("target")

    cols = st.columns(6)
    cols[0].metric("Candidates", f"{inp['n_valid']:,}", f"{inp['n_rejected']} rejected"
                   if inp["n_rejected"] else None)
    cols[1].metric("Target", (target or {}).get("pdb_id", "—"))
    cols[2].metric("Length bounds", f"{cfg['input']['min_length']}–{cfg['input']['max_length']} nt")
    cols[3].metric("Shuffled controls", cfg["tier1"]["shuffle"]["n_shuffles"])
    cols[4].metric("Seed", cfg["run"]["seed"])
    cols[5].metric("Tier 1 runtime", f"{artifact['diagnostics']['runtime_seconds']['tier1']:.0f} s")

    with st.expander("Run provenance", expanded=False):
        left, right = st.columns(2)
        left.markdown(
            f"**Input** `{inp['filename']}`  \n"
            f"sha256 `{(inp.get('sha256') or '')[:16]}…`  \n"
            f"**Corpus** `{artifact['corpus']['corpus_id']}`  \n"
            f"{artifact['corpus']['n_sequences']:,} reference sequences, "
            f"{artifact['corpus']['n_dropped']} dropped  \n"
            f"**Composite method** `{artifact['diagnostics']['composite_method']}`  \n"
            f"**Scoring signature** `{artifact['scoring_signature']}`"
        )
        versions = artifact["versions"]
        git = artifact.get("git") or {}
        right.markdown(
            f"**AptaRank** {versions.get('aptarank')} "
            f"(schema {artifact['artifact_schema_version']})  \n"
            f"**ViennaRNA** {versions.get('viennarna')} · **forgi** {versions.get('forgi')} · "
            f"**ushuffle** {versions.get('ushuffle')}  \n"
            f"**fpocket** {versions.get('fpocket') or '—'}  \n"
            f"**git** `{(git.get('commit') or 'not a checkout')[:12]}`"
            f"{' (dirty)' if git.get('dirty') else ''}"
        )


# -- (b) ranked candidate list ------------------------------------------


def ranked_table(candidates: Sequence[Mapping[str, Any]]) -> pd.DataFrame:
    rows = []
    for c in candidates:
        shuffle = c.get("shuffle") or {}
        tier2 = c.get("tier2") or {}
        rows.append(
            {
                "rank": c["rank"],
                "candidate": c["candidate_id"],
                "Tier 1": c["tier1_score"],
                "Tier 2 band": BAND_LABEL.get(tier2.get("band", "not_evaluated"), "—"),
                "shuffle": "pass" if shuffle.get("pass") else ("fail" if shuffle.get("pass") is False else "—"),
                "len": c.get("length"),
                "sequence": c["sequence"],
            }
        )
    return pd.DataFrame(rows)


#: The left-hand list stays narrow: the full sequence and every other field are
#: one click away in the detail panel, and a squeezed table hides its own
#: columns rather than scrolling them.
LIST_COLUMNS = ["rank", "candidate", "Tier 1", "Tier 2 band", "shuffle"]


def candidate_list(table: pd.DataFrame, key: str = "ranked") -> int | None:
    event = st.dataframe(
        table[LIST_COLUMNS],
        hide_index=True,
        use_container_width=True,
        height=520,
        on_select="rerun",
        selection_mode="single-row",
        key=key,
        column_config={
            "rank": st.column_config.NumberColumn("#", width="small"),
            "Tier 1": st.column_config.ProgressColumn(
                "Tier 1", min_value=0.0, max_value=1.0, format="%.3f",
                help="Corpus-calibrated composite. Absolute: it does not depend "
                     "on what else was submitted.",
            ),
            "Tier 2 band": st.column_config.TextColumn("Tier 2", width="small"),
            "shuffle": st.column_config.TextColumn("shuffle", width="small"),
        },
    )
    rows = event.selection.rows if event and event.selection else []
    return int(rows[0]) if rows else None


# -- (c) candidate detail ------------------------------------------------


def candidate_detail(candidate: Mapping[str, Any], colors: Mapping[str, Any]) -> None:
    st.markdown(f"#### Rank {candidate['rank']} · `{candidate['candidate_id']}`")

    tier2 = candidate.get("tier2") or {}
    shuffle = candidate.get("shuffle") or {}
    cols = st.columns(4)
    cols[0].metric("Tier 1 score", f"{candidate['tier1_score']:.3f}",
                   help="Corpus-calibrated composite of five named criteria. "
                        "Absolute — independent of the rest of the batch.")
    cols[1].metric(
        "Tier 2 band", BAND_LABEL.get(tier2.get("band", "not_evaluated"), "—"),
        help="Control-relative geometric size agreement with this target's "
             "cavity. 'Strong' means better agreement than ~95% of shuffled "
             "controls — not a strong candidate, and not evidence of binding.",
    )
    cols[2].metric(
        "Shuffled control",
        "pass" if shuffle.get("pass") else ("fail" if shuffle.get("pass") is False else "—"),
        f"p = {shuffle['p_value']:.3f}" if shuffle.get("p_value") is not None else None,
        # "off": a p-value is not a change, and a green arrow beside it would
        # read as a verdict the number does not carry.
        delta_color="off",
    )
    cols[3].metric("Length", f"{candidate.get('length')} nt")

    structure = candidate.get("structure") or {}
    st.markdown(
        f"<div class='apt-seq'>{structure.get('dot_bracket','')}<br>{candidate['sequence']}</div>",
        unsafe_allow_html=True,
    )

    left, right = st.columns([1, 1])
    with left:
        st.altair_chart(criterion_bars(candidate, colors), use_container_width=True)
    with right:
        svg = structure.get("svg")
        if svg:
            st.markdown(
                f"<img alt='predicted secondary structure' style='width:100%;max-width:340px' "
                f"src='data:image/svg+xml;base64,{base64.b64encode(svg.encode()).decode()}'/>",
                unsafe_allow_html=True,
            )
        else:
            st.caption("No structure diagram stored for this candidate "
                       "(diagrams are rendered for the top N only).")


def criterion_bars(candidate: Mapping[str, Any], colors: Mapping[str, Any]) -> alt.Chart:
    """One bar per named criterion — the composite decomposed (P2)."""
    rows = [
        {
            "criterion": CRITERION_LABELS.get(name, name),
            "score": entry["score"],
            "value": entry["value"],
            "raw": name,
        }
        for name, entry in (candidate.get("criteria") or {}).items()
    ]
    frame = pd.DataFrame(rows)
    base = alt.Chart(frame).encode(
        # labelOverlap=False: every criterion must be named. A composite whose
        # parts are unlabelled is the opaque score the design rules out (P2).
        y=alt.Y("criterion:N", sort=None, title=None,
                axis=alt.Axis(labelLimit=200, labelOverlap=False, ticks=False, domain=False)),
    )
    bars = base.mark_bar(
        height=10, cornerRadiusTopRight=4, cornerRadiusBottomRight=4,
        color=colors["series"],
    ).encode(
        x=alt.X("score:Q", scale=alt.Scale(domain=[0, 1]), title="corpus percentile score"),
        tooltip=[
            alt.Tooltip("criterion:N", title="Criterion"),
            alt.Tooltip("value:Q", title="Measured value", format=".4f"),
            alt.Tooltip("score:Q", title="Corpus score", format=".3f"),
        ],
    )
    labels = base.mark_text(
        align="left", dx=5, fontSize=11, color=colors["text_secondary"]
    ).encode(x=alt.X("score:Q"), text=alt.Text("score:Q", format=".2f"))
    return (
        (bars + labels)
        .properties(height=alt.Step(26))
        .configure(**chart_config(colors))
    )


# -- (d) explanation -----------------------------------------------------


def explanation_panel(candidate: Mapping[str, Any]) -> None:
    st.markdown("##### Why this candidate ranks here")
    chips = "".join(
        f"<span class='apt-chip apt-chip-{c.get('kind','neutral')}'>{c['label']}</span>"
        for c in candidate.get("evidence_chips", [])
    )
    if chips:
        st.markdown(chips, unsafe_allow_html=True)

    text = candidate.get("explanation") or ""
    caveat = "Tier 2 reflects geometric plausibility only — it is not evidence of binding."
    body = text.replace(caveat, "").strip()
    st.write(body)
    st.markdown(f"<div class='apt-caveat'>{caveat}</div>", unsafe_allow_html=True)


# -- (e) target panel ----------------------------------------------------


def target_panel(artifact: Mapping[str, Any]) -> None:
    target = artifact.get("target")
    st.markdown("##### Target")
    if not target:
        st.info(
            "Tier 1 only — no target was supplied for this run, so no "
            "target-aware evidence exists for any candidate."
        )
        return

    pocket = target["selected_pocket"]
    st.markdown(f"**{target['pdb_id']}** chain {target['chain']}")
    if target.get("name"):
        st.caption(target["name"])

    cols = st.columns(2)
    cols[0].metric("Cavities detected", target["n_pockets"])
    cols[1].metric("Selected pocket", f"#{pocket['index']}")
    cols[0].metric("Characteristic width", f"{pocket['d_pocket_A']:.1f} Å")
    cols[1].metric("Volume", f"{pocket['volume_A3']:.0f} Å³")

    method = target["pocket_selection"]
    if method == "active_site_overlap":
        st.caption("Pocket selected by overlap with literature-confirmed active-site residues.")
    else:
        st.warning(
            f"Pocket selected automatically (`{method}`) — fpocket's own score is "
            f"not guaranteed to identify the functional cavity.",
            icon="⚠️",
        )
    for warning in target.get("selection_warnings", []):
        st.warning(warning, icon="⚠️")
    if pocket.get("shape_warning"):
        st.warning(
            "This cavity is oddly shaped: its robust extent and equivalent-sphere "
            "diameter disagree by more than 2×. The geometric comparison assumes a "
            "roughly convex pocket.",
            icon="⚠️",
        )

    status = target.get("electrostatics_status")
    if status == "success":
        potential = target["electrostatic_mean_potential"]
        compatible = target["electrostatic_compatible"]
        st.metric(
            "Cavity electrostatics",
            f"{potential:+.2f} kT/e",
            "hospitable to RNA" if compatible else "repulsive to RNA",
            delta_color="off",
        )
        st.caption("Target-level only: identical for every candidate, so it is "
                   "deliberately not part of any per-candidate score.")
    else:
        st.caption(f"Electrostatics: {status or 'not computed'}.")

    if target.get("retained_hetero"):
        st.caption(
            f"Retained non-water hetero groups: {', '.join(sorted(set(target['retained_hetero'])))} "
            f"— these shape the detected cavity."
        )
    st.caption(f"Target bundle `{target['bundle_id'][:12]}`")


# -- (f) Tier 1 vs Tier 2 scatter ---------------------------------------


def tier_scatter(
    candidates: Sequence[Mapping[str, Any]],
    thresholds: Mapping[str, Any] | None,
    diagnostics: Mapping[str, Any],
    colors: Mapping[str, Any],
    selected_id: str | None = None,
) -> alt.LayerChart | None:
    rows = []
    for c in candidates:
        tier2 = c.get("tier2") or {}
        if tier2.get("status") != "evaluated":
            continue
        rows.append(
            {
                "candidate": c["candidate_id"],
                "rank": c["rank"],
                "tier1": c["tier1_score"],
                "control_percentile": tier2.get("control_percentile_flexible"),
                "mismatch": tier2.get("absolute_mismatch_flexible_A"),
                "band": BAND_LABEL.get(tier2.get("band", "not_evaluated")),
                "selected": c["candidate_id"] == selected_id,
            }
        )
    if not rows:
        return None

    frame = pd.DataFrame(rows)
    # Only the bands actually present: a legend entry with no points on the
    # plot invites the reader to hunt for a category that isn't there.
    full_domain, full_range = band_scale(colors)
    present = {b for b in frame["band"]}
    pairs = [
        (BAND_LABEL[b], colour)
        for b, colour in zip(full_domain, full_range)
        if BAND_LABEL[b] in present
    ]
    domain = [label for label, _ in pairs]
    range_ = [colour for _, colour in pairs]

    points = (
        alt.Chart(frame)
        .mark_circle(size=90, opacity=0.85, stroke=colors["surface"], strokeWidth=2)
        .encode(
            x=alt.X("tier1:Q", title="Tier 1 score (intrinsic quality)",
                    scale=alt.Scale(zero=False, nice=True),
                    axis=alt.Axis(tickCount=8)),
            y=alt.Y("control_percentile:Q", title="Tier 2 control percentile",
                    scale=alt.Scale(domain=[0, 1])),
            color=alt.Color(
                "band:N",
                title="Control-relative geometric agreement",
                scale=alt.Scale(domain=domain, range=range_),
                legend=alt.Legend(orient="top", direction="horizontal"),
            ),
            tooltip=[
                alt.Tooltip("candidate:N", title="Candidate"),
                alt.Tooltip("rank:Q", title="Rank"),
                alt.Tooltip("tier1:Q", title="Tier 1", format=".3f"),
                alt.Tooltip("control_percentile:Q", title="Control percentile", format=".3f"),
                alt.Tooltip("mismatch:Q", title="Loop/cavity mismatch (Å)", format=".1f"),
                alt.Tooltip("band:N", title="Band"),
            ],
        )
    )

    layers = [points]
    bands = (thresholds or {}).get("band_percentiles") or {}
    values = [v for v in (bands.get("moderate"), bands.get("strong")) if v is not None]
    if values:
        # Rules only, no in-plot text: the band boundaries are named in the
        # caption and on the y-axis, and a floating label in a dense scatter
        # lands on the data no matter which edge it is anchored to.
        layers.append(
            alt.Chart(pd.DataFrame({"y": values}))
            .mark_rule(strokeDash=[4, 3], color=colors["muted"], strokeWidth=1)
            .encode(y="y:Q")
        )

    if selected_id and frame["selected"].any():
        layers.append(
            alt.Chart(frame[frame["selected"]])
            .mark_point(size=260, shape="circle", strokeWidth=2,
                        stroke=colors["text_primary"], filled=False)
            .encode(x="tier1:Q", y="control_percentile:Q")
        )

    return (
        alt.layer(*layers)
        .properties(height=320)
        .configure(**chart_config(colors))
    )


def threshold_caption(thresholds: Mapping[str, Any] | None) -> str:
    """Say what the dashed lines are, and what the bands are relative to."""
    if not thresholds:
        return ""
    bands = thresholds.get("band_percentiles") or {}
    return (
        f"Dashed lines: moderate ≥ {bands.get('moderate', 0.75):.2f}, "
        f"strong ≥ {bands.get('strong', 0.95):.2f} of the control distribution. "
        f"Bands are relative to {thresholds.get('n_controls', 0):,} fixed "
        f"dinucleotide-shuffled controls scored against this same cavity "
        f"({thresholds.get('d_pocket_A', float('nan')):.1f} Å) — not to the other "
        f"candidates in this batch, so they do not shift when the batch changes."
    )


def independence_caption(diagnostics: Mapping[str, Any]) -> str:
    spearman = diagnostics.get("spearman_tier1_tier2") or {}
    if spearman.get("rho") is None:
        return f"Tier independence: not computable ({spearman.get('reason', 'no data')})."
    interpretation = spearman.get("interpretation") or ""
    if interpretation:
        interpretation = interpretation[0].upper() + interpretation[1:] + "."
    return (
        f"Spearman ρ = {spearman['rho']:+.3f} (p = {spearman['p_value']:.3g}, "
        f"n = {spearman['n']}) between Tier 1 score and the Tier 2 control percentile. "
        f"{interpretation}"
    )
