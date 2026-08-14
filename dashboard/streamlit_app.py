"""AptaRank dashboard.

Three views: start an analysis, watch it run, read the results.

The dashboard is a thin orchestration and visualisation client. It stages
uploaded files, invokes the same headless `aptarank` command line used for
reproducible runs, and renders only validated artifacts. No scoring code lives
here — which is what lets the paper's figures and the live demo come from the
same runs.

    streamlit run dashboard/streamlit_app.py
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import streamlit as st

REPO_ROOT = Path(__file__).resolve().parent.parent
#: On a server the mutable state lives outside the checkout, which every deploy
#: replaces. Locally the two are the same directory.
DATA_DIR = Path(os.environ.get("APTARANK_DATA_DIR", REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT))

from dashboard import jobs, newrun, progress_view, views  # noqa: E402
from dashboard.theme import CSS, palette  # noqa: E402

VIEWS = ("New analysis", "Progress", "Results", "Recent analyses")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--artifact", help="run artifact JSON to open on startup")
    parser.add_argument("--runs-dir", default=str(REPO_ROOT / "runs"))
    known, _ = parser.parse_known_args()
    return known


@st.cache_data(show_spinner=False)
def load_artifact(path: str, _mtime: float) -> dict:
    with open(path, "r", encoding="utf-8") as fh:
        return json.load(fh)


def theme_type() -> str:
    try:
        return st.context.theme.type or "light"
    except Exception:  # pragma: no cover - older Streamlit
        return st.get_option("theme.base") or "light"


def available_artifacts(runs_dir: Path) -> list[tuple[str, str]]:
    """(label, path) for everything readable, newest first."""
    entries: list[tuple[float, str, str]] = []
    for job in jobs.discover(runs_dir, limit=50):
        if job.artifact_path.exists():
            entries.append(
                (job.artifact_path.stat().st_mtime, job.label, str(job.artifact_path))
            )
    for path in runs_dir.glob("*.json"):
        entries.append((path.stat().st_mtime, path.stem, str(path)))
    entries.sort(reverse=True)
    return [(label, path) for _mtime, label, path in entries]


def main() -> None:
    st.set_page_config(
        page_title="AptaRank", page_icon="🧬", layout="wide",
        initial_sidebar_state="expanded",
    )
    st.markdown(CSS, unsafe_allow_html=True)
    args = parse_args()
    runs_dir = Path(args.runs_dir)
    runs_dir.mkdir(parents=True, exist_ok=True)

    if args.artifact and "artifact_path" not in st.session_state:
        st.session_state["artifact_path"] = args.artifact

    # Start any run that has been waiting for a slot. Cheap, and it means the
    # queue advances whenever anyone has the page open.
    jobs.pump(runs_dir)

    with st.sidebar:
        st.markdown("## 🧬 AptaRank")
        st.caption("Interpretable ranking of generated RNA aptamers")
        running = jobs.active(runs_dir)
        default = VIEWS.index(st.session_state.get("view", "New analysis"))
        view = st.radio("Go to", VIEWS, index=default, label_visibility="collapsed")
        st.session_state["view"] = view
        if running:
            load = jobs.slots(runs_dir)
            st.info(
                f"{running.label}\n\n{load['running']}/{load['capacity']} slots busy"
                + (f", {load['queued']} waiting" if load["queued"] else ""),
                icon="⏳",
            )
        st.divider()
        st.caption(
            "Tier 1 ranks candidates on aptamer-likeness. Tier 2 annotates them "
            "with aptamer-target compatibility by performing a geometric "
            "agreement check against a target. Tier 2 is not a binding "
            "prediction, and it does not change the Tier 1 ranking."
        )

    if view == "New analysis":
        newrun.render(REPO_ROOT, runs_dir, data_dir=DATA_DIR)
        watching = st.session_state.get("watch_job")
        if watching:
            st.session_state["view"] = "Progress"
            st.rerun()

    elif view == "Progress":
        _progress_view(runs_dir)

    elif view == "Recent analyses":
        st.markdown("### Recent analyses")
        progress_view.recent(runs_dir)

    else:
        _results_view(runs_dir)


def _progress_view(runs_dir: Path) -> None:
    job_id = st.session_state.get("watch_job")
    found = {job.job_id: job for job in jobs.discover(runs_dir)}
    job = found.get(job_id) or jobs.active(runs_dir)

    if job is None:
        st.info(
            "Nothing is running. Start an analysis, or open a finished one from "
            "Recent analyses.",
            icon="ℹ️",
        )
        return

    st.session_state["watch_job"] = job.job_id
    artifact_path = progress_view.render(job, runs_dir)
    if artifact_path:
        if st.button("Open the results", type="primary"):
            st.session_state["artifact_path"] = artifact_path
            st.session_state["view"] = "Results"
            st.rerun()


def _results_view(runs_dir: Path) -> None:
    colors = palette(theme_type())
    options = available_artifacts(runs_dir)

    if not options:
        st.info(
            "No results yet. Go to **New analysis** to rank a set of sequences.",
            icon="ℹ️",
        )
        return

    paths = [path for _label, path in options]
    current = st.session_state.get("artifact_path")
    index = paths.index(current) if current in paths else 0

    header, picker = st.columns([3, 1])
    with header:
        st.markdown("## Results")
    with picker:
        labels = {path: label for label, path in options}
        chosen = st.selectbox(
            "Analysis", paths, index=index, format_func=lambda p: labels.get(p, Path(p).name),
        )
    st.session_state["artifact_path"] = chosen

    artifact = load_artifact(chosen, Path(chosen).stat().st_mtime)

    # The state badge is top-level and always present: whether a result may be
    # quoted is not a detail to be discovered in a config panel.
    if artifact.get("publication_eligible", True):
        st.success("Publication-eligible run", icon="✅")
    else:
        reasons = artifact.get("development_reasons") or ["development settings"]
        # The same wording, from the same table, as the New Analysis page shows
        # before the run — the two verdicts must never read differently.
        from dashboard.inputs import DEVELOPMENT_REASON_TEXT as pretty

        st.markdown(
            "<div class='apt-dev-banner'><b>Development run — not a result.</b> "
            + "; ".join(pretty.get(r, r) for r in reasons)
            + ". This artifact is marked <code>publication_eligible: false</code> "
              "and must not back any published claim.</div>",
            unsafe_allow_html=True,
        )

    views.run_configuration(artifact)
    st.divider()

    candidates = artifact["candidates"]
    table = views.ranked_table(candidates)
    left, centre, right = st.columns([1.15, 1.5, 0.85], gap="medium")

    with left:
        st.markdown("##### Ranked candidates")
        st.caption(f"{len(candidates):,} scored · ranked by aptamer-likeness")
        selected_row = views.candidate_list(table)

    candidate = candidates[selected_row if selected_row is not None else 0]

    with centre:
        views.candidate_detail(candidate, colors)
        views.explanation_panel(candidate)

    with right:
        views.target_panel(artifact)

    st.divider()
    st.markdown("##### Aptamer-likeness vs aptamer-target compatibility")
    chart = views.tier_scatter(
        candidates, artifact.get("tier2_thresholds"), artifact["diagnostics"],
        colors, selected_id=candidate["candidate_id"],
    )
    if chart is None:
        st.info(
            "No target was used in this analysis. Add a protein target to annotate "
            "each candidate with aptamer-target compatibility.",
            icon="ℹ️",
        )
    else:
        st.altair_chart(chart, use_container_width=True)
        st.caption(views.threshold_caption(artifact.get("tier2_thresholds")))
        st.caption(views.independence_caption(artifact["diagnostics"]))
        with st.expander("Table view of the plotted points"):
            st.dataframe(
                table[table["compatibility"] != "not evaluated"],
                hide_index=True, use_container_width=True,
            )

    st.divider()
    st.download_button(
        "Export shortlist (CSV)",
        table.to_csv(index=False).encode("utf-8"),
        file_name=f"{artifact['run_id']}_shortlist.csv",
        mime="text/csv",
    )


if __name__ == "__main__":
    main()
