"""Panel renderers. Each one reads the run artifact and computes nothing.

Keeping every calculation in the pipeline is what lets the paper's figures be
generated from the same artifacts the live demo displays.
"""

from __future__ import annotations

import base64
from pathlib import Path
from typing import Any, Mapping, Sequence

import altair as alt
import pandas as pd
import streamlit as st

from .inputs import BINDING_MODE_DESCRIPTION, BINDING_MODE_LABEL
from .theme import BAND_LABEL, band_scale, chart_config

CRITERION_LABELS = {
    "mfe_norm": "Structure stability",
    "ensemble_defect": "Fold definition",
    "positional_entropy_mean": "Fold certainty",
    "stem_fraction": "Structural composition",
    "gc_fraction": "Sequence composition",
}

#: User-facing names for the two tiers (refinements §1.1). The artifact keeps
#: the paper's vocabulary; every screen uses these.
TIER1_NAME = "aptamer-likeness"
TIER2_NAME = "aptamer-target compatibility"

#: The same statement as `aptarank.TIER2_CAVEAT`, in the user-facing wording.
DISPLAY_CAVEAT = (
    "Aptamer-target compatibility is a geometric agreement check — it is not "
    "evidence of binding."
)

#: What the shuffled-control column is actually testing (§4.8). Shown wherever
#: the pass/fail appears, because "pass" on its own reads as a verdict on the
#: candidate rather than on what its score is made of.
SHUFFLE_HELP = (
    "Whether this candidate's structure is doing more than its letter "
    "composition. Each candidate is compared against shuffled versions of "
    "itself that keep the same nucleotide composition but scramble the order. "
    "'Pass' means the candidate scores better than ~95% of its own shuffles — "
    "its structural quality reflects how the sequence is arranged, not just "
    "which letters it contains."
)


# -- (a) run configuration ----------------------------------------------


def run_configuration(artifact: Mapping[str, Any]) -> None:
    cfg = artifact["config"]
    inp = artifact["input"]
    target = artifact.get("target")

    cols = st.columns(6)
    cols[0].metric("Candidates", f"{inp['n_valid']:,}", f"{inp['n_rejected']} rejected"
                   if inp["n_rejected"] else None)
    cols[1].metric("Target", (target or {}).get("pdb_id", "—"))
    cols[2].metric(
        "Binding mode",
        BINDING_MODE_LABEL.get(artifact.get("binding_mode"), "—"),
        help="How the expert asserted an aptamer would engage this target. "
             "AptaRank does not infer it; the comparison adapts to it.",
    )
    cols[3].metric("Shuffled controls", cfg["tier1"]["shuffle"]["n_shuffles"])
    cols[4].metric("Seed", cfg["run"]["seed"])
    cols[5].metric("Ranking runtime", f"{artifact['diagnostics']['runtime_seconds']['tier1']:.0f} s")

    with st.expander("Run provenance", expanded=False):
        left, right = st.columns(2)
        corpus = artifact["corpus"]
        # The name the user uploaded, with the hash that actually identifies the
        # content. Showing only the internal staged filename was a provenance
        # record of a name nobody recognised.
        left.markdown(
            f"**Sequences** uploaded as `{_display_name(inp)}`  \n"
            f"sha256 `{(inp.get('sha256') or '')[:16]}…`  \n"
            f"stored as `{Path(inp['filename']).name if inp.get('filename') else '—'}`  \n"
            f"**Reference library** uploaded as "
            f"`{corpus.get('original_filename') or Path(corpus['path']).name}`  \n"
            f"sha256 `{(corpus.get('corpus_sha256') or '')[:16]}…`  \n"
            f"{corpus['n_sequences']:,} reference aptamers, "
            f"{corpus['n_dropped']} dropped  \n"
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
                "aptamer-likeness": c["tier1_score"],
                "compatibility": BAND_LABEL.get(tier2.get("band", "not_evaluated"), "—"),
                "shuffle": "pass" if shuffle.get("pass") else ("fail" if shuffle.get("pass") is False else "—"),
                "len": c.get("length"),
                "sequence": c["sequence"],
            }
        )
    return pd.DataFrame(rows)


#: The left-hand list stays narrow: the full sequence and every other field are
#: one click away in the detail panel, and a squeezed table hides its own
#: columns rather than scrolling them.
LIST_COLUMNS = ["rank", "candidate", "aptamer-likeness", "compatibility", "shuffle"]


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
            "aptamer-likeness": st.column_config.ProgressColumn(
                "aptamer-likeness", min_value=0.0, max_value=1.0, format="%.3f",
                help="How much this candidate looks like a validated aptamer, "
                     "calibrated against the reference library. Absolute: it does "
                     "not depend on what else was submitted.",
            ),
            "compatibility": st.column_config.TextColumn(
                "compatibility", width="small",
                help="Aptamer-target compatibility: how well this candidate's "
                     "shape agrees with the target's binding site, relative to "
                     "shuffled controls. Not a binding prediction.",
            ),
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
    mode = tier2.get("binding_mode", "pocket")
    cols = st.columns(4)
    cols[0].metric("Aptamer-likeness", f"{candidate['tier1_score']:.3f}",
                   help="Corpus-calibrated composite of five named criteria. "
                        "Absolute — independent of the rest of the batch.")
    cols[1].metric(
        "Compatibility", BAND_LABEL.get(tier2.get("band", "not_evaluated"), "—"),
        help=(
            "How well this candidate's shape agrees with the target's binding "
            "site, relative to shuffled controls. 'Strong' means better agreement "
            "than ~95% of controls — not a strong candidate, and not evidence of "
            "binding. Mode: " + BINDING_MODE_LABEL.get(mode, mode)
        ),
    )
    cols[2].metric(
        "Shuffled control",
        "pass" if shuffle.get("pass") else ("fail" if shuffle.get("pass") is False else "—"),
        f"p = {shuffle['p_value']:.3f}" if shuffle.get("p_value") is not None else None,
        # "off": a p-value is not a change, and a green arrow beside it would
        # read as a verdict the number does not carry.
        delta_color="off",
        help=SHUFFLE_HELP,
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

    # The stored caveat is stripped by importing the exact constant rather than
    # by repeating the sentence here: a copy would stop matching the moment the
    # canonical wording changed, and the caveat would appear twice.
    from aptarank import TIER2_CAVEAT

    text = candidate.get("explanation") or ""
    body = text.replace(TIER2_CAVEAT, "").strip()
    st.write(body)
    st.markdown(f"<div class='apt-caveat'>{DISPLAY_CAVEAT}</div>", unsafe_allow_html=True)


# -- (e) target panel ----------------------------------------------------


def target_panel(artifact: Mapping[str, Any]) -> None:
    target = artifact.get("target")
    st.markdown("##### Target")
    if not target:
        st.info(
            "Aptamer-likeness only — no target was supplied for this run, so no "
            "aptamer-target compatibility evidence exists for any candidate."
        )
        return

    mode = artifact.get("binding_mode") or target.get("binding_mode") or "pocket"
    st.markdown(f"**{target.get('identifier', target['pdb_id'])}** chain {target['chain']}")
    if target.get("name"):
        st.caption(target["name"])
    st.caption(
        f"{BINDING_MODE_LABEL.get(mode, mode)} — {BINDING_MODE_DESCRIPTION.get(mode, '')}"
    )

    if target.get("structure_kind") == "predicted":
        st.warning(
            "Predicted structure (AlphaFold), not an experiment."
            + (" A predicted model may not show an interface that only forms when "
               "a binding partner is present, so treat this with extra caution."
               if mode == "surface" else ""),
            icon="⚠️",
        )

    patch = target.get("patch")
    pocket = target.get("selected_pocket")
    if mode == "surface" and patch:
        cols = st.columns(2)
        cols[0].metric("Binding-site area", f"{patch['patch_area_A2']:.0f} Å²")
        cols[1].metric("Site residues", patch["n_residues"])
        cols[0].metric("Flatness", f"{patch['planarity_A']:.1f} Å",
                       help="Thickness across the thinnest direction of the patch. "
                            "Smaller is flatter.")
        cols[1].metric("Elongation", f"{patch['elongation']:.2f}",
                       help="Longest / shortest spread. Near 1 is round; large is "
                            "groove-like.")
        st.caption(
            f"Area measured with freeSASA over the {patch['n_residues']} "
            f"binding-site residues you specified, on the isolated chain."
        )
        if patch.get("buried_residue_numbers"):
            st.warning(
                f"Binding-site residues {patch['buried_residue_numbers']} have no "
                f"exposed surface — check the residue numbering and chain.",
                icon="⚠️",
            )
        if patch.get("shape_warning"):
            st.warning(
                "This patch is not flat. Surface-mode agreement assumes a roughly "
                "planar face, so read the band with that in mind.",
                icon="⚠️",
            )
    elif pocket:
        cols = st.columns(2)
        cols[0].metric("Cavities detected", target["n_pockets"])
        cols[1].metric("Selected cavity", f"#{pocket['index']}")
        cols[0].metric("Characteristic width", f"{pocket['d_pocket_A']:.1f} Å")
        cols[1].metric("Volume", f"{pocket['volume_A3']:.0f} Å³")

        method = target.get("pocket_selection")
        if method in ("target_site_overlap", "active_site_overlap"):
            st.caption(
                "Cavity selected by overlap with literature-confirmed "
                "binding-site residues."
            )
        else:
            st.warning(
                f"Cavity selected automatically (`{method}`) — fpocket's own score "
                f"is not guaranteed to identify the functional cavity.",
                icon="⚠️",
            )
        if pocket.get("shape_warning"):
            st.warning(
                "This cavity is oddly shaped: its robust extent and equivalent-sphere "
                "diameter disagree by more than 2×. The geometric comparison assumes a "
                "roughly convex pocket.",
                icon="⚠️",
            )

    for warning in [*(target.get("selection_warnings") or []),
                    *(target.get("preparation_warnings") or [])]:
        st.warning(warning, icon="⚠️")

    status = target.get("electrostatics_status")
    if status == "success" and target.get("electrostatic_mean_potential") is not None:
        potential = target["electrostatic_mean_potential"]
        compatible = target["electrostatic_compatible"]
        st.metric(
            "Binding-site electrostatics",
            f"{potential:+.2f} kT/e",
            "hospitable to RNA" if compatible else "repulsive to RNA",
            delta_color="off",
        )
        st.caption(
            "A property of the target: identical for every candidate. In surface "
            "mode it is part of the reported agreement score, but because it is "
            "the same for everyone it can never reorder candidates or change a band."
            if mode == "surface" else
            "Target-level only: identical for every candidate, so it is "
            "deliberately not part of any per-candidate score."
        )
    else:
        st.caption(f"Electrostatics: {status or 'not computed'}.")

    if target.get("retained_hetero"):
        st.caption(
            f"Retained non-water hetero groups: {', '.join(sorted(set(target['retained_hetero'])))} "
            f"— these shape the measured site."
        )
    st.caption(f"Prepared target `{target['bundle_id'][:12]}`")


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
                "control_percentile": tier2.get(
                    "control_percentile", tier2.get("control_percentile_flexible")
                ),
                "mismatch": tier2.get("disagreement"),
                "band": BAND_LABEL.get(tier2.get("band", "not_evaluated")),
                "selected": c["candidate_id"] == selected_id,
            }
        )
    if not rows:
        return None

    units = (thresholds or {}).get("units", "Å")
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
            x=alt.X("tier1:Q", title="Aptamer-likeness",
                    scale=alt.Scale(zero=False, nice=True),
                    axis=alt.Axis(tickCount=8)),
            y=alt.Y("control_percentile:Q",
                    title="Aptamer-target compatibility (control percentile)",
                    scale=alt.Scale(domain=[0, 1])),
            color=alt.Color(
                "band:N",
                title="Compatibility, relative to shuffled controls",
                scale=alt.Scale(domain=domain, range=range_),
                legend=alt.Legend(orient="top", direction="horizontal"),
            ),
            tooltip=[
                alt.Tooltip("candidate:N", title="Candidate"),
                alt.Tooltip("rank:Q", title="Rank"),
                alt.Tooltip("tier1:Q", title="Aptamer-likeness", format=".3f"),
                alt.Tooltip("control_percentile:Q", title="Control percentile", format=".3f"),
                alt.Tooltip("mismatch:Q", title=f"Size mismatch ({units})", format=".1f"),
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
    target = thresholds.get("target") or {}
    if thresholds.get("binding_mode") == "surface":
        model = (thresholds.get("parameters") or {}).get("footprint_model", "")
        described = (
            "the folded molecule's size" if model == "radius_of_gyration"
            else "nucleotide count"
        )
        measured = (
            f"binding-site patch ({target.get('patch_area_A2', float('nan')):.0f} Å²), "
            f"with each candidate's footprint taken from {described}"
        )
    else:
        measured = f"cavity ({target.get('d_pocket_A', float('nan')):.1f} Å)"
    return (
        f"Dashed lines: moderate ≥ {bands.get('moderate', 0.75):.2f}, "
        f"strong ≥ {bands.get('strong', 0.95):.2f} of the control distribution. "
        f"Bands are relative to {thresholds.get('n_controls', 0):,} fixed "
        f"dinucleotide-shuffled controls scored against this same {measured} — "
        f"not to the other candidates in this batch, so they do not shift when "
        f"the batch changes."
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
        f"n = {spearman['n']}) between aptamer-likeness and the compatibility "
        f"control percentile. {interpretation}"
    )


def _display_name(inp: Mapping[str, Any]) -> str:
    """What the user called their file, falling back to what is on disk."""
    return inp.get("original_filename") or (
        Path(inp["filename"]).name if inp.get("filename") else "—"
    )
