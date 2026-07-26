"""Live view of a running analysis.

Reads the job's progress file from disk on a timer. Nothing here depends on the
browser session that started the run, so closing the tab, refreshing, or even
restarting Streamlit leaves the analysis untouched.

The bar shows the *current stage's* own progress, never a single blended
percentage: the stages differ in cost by more than an order of magnitude, so
one overall number would be invented rather than measured.
"""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import streamlit as st

from . import jobs
from .jobs import Job

STAGE_ORDER = ["ingest", "corpus", "tier1", "bank", "tier2", "diagrams", "artifact"]


def render(job: Job, runs_dir: Path) -> str | None:
    """Draw the progress panel. Returns the artifact path once finished."""
    state = job.state()

    st.markdown(f"### {job.label}")
    st.caption(f"Started {_when(job.created_utc)} · job `{job.job_id}`")

    if state["status"] in ("pending", "running"):
        _running(job, state)
        return None
    if state["status"] == "stalled":
        st.warning(
            "This analysis has not reported progress for a long time. It may "
            "still be working on a slow step, or it may be stuck.",
            icon="⏳",
        )
        _running(job, state)
        return None
    if state["status"] == "failed":
        _failed(job, state)
        return None

    st.success(
        f"Finished — {state.get('n_candidates') or ''} candidates ranked."
        if state.get("n_candidates") else "Finished.",
        icon="✅",
    )
    for message in state["warnings"]:
        st.warning(message, icon="⚠️")
    return state.get("artifact_path")


@st.fragment(run_every=2)
def _running(job: Job, _initial: dict[str, Any]) -> None:
    """Re-reads the job on a timer without re-running the whole page."""
    state = job.state()

    label = state.get("stage_label") or "Starting…"
    completed, total = state.get("completed"), state.get("total")
    done = [s for s in STAGE_ORDER if s in state["stages_completed"]]

    st.progress(
        (completed / total) if (completed and total) else 0.0,
        text=f"{label}" + (f" — {completed:,} of {total:,}" if completed and total else ""),
    )
    st.caption(
        f"Completed steps: {', '.join(done) if done else 'none yet'}"
        f" · last update {_ago(state.get('updated_utc'))}"
    )

    for message in state["warnings"]:
        st.warning(message, icon="⚠️")

    if st.button("Stop this analysis", key=f"cancel_{job.job_id}"):
        job.cancel()
        st.rerun()

    if state["status"] in ("completed", "failed"):
        st.rerun()   # leave the fragment and redraw the page in its final state


def _failed(job: Job, state: dict[str, Any]) -> None:
    st.error(
        f"**The analysis could not finish.**\n\n{state.get('error_message') or 'No reason recorded.'}",
        icon="🚫",
    )
    st.markdown(_next_step(state.get("error_code")))

    logs = job.logs()
    with st.expander("Technical details (for whoever maintains this tool)"):
        st.code(f"error code: {state.get('error_code')}\n\n"
                f"{logs['stderr.log'][-4000:] or logs['stdout.log'][-4000:]}",
                language="text")
    st.download_button(
        "Download the full log",
        (logs["stdout.log"] + "\n\n=== errors ===\n\n" + logs["stderr.log"]).encode("utf-8"),
        file_name=f"{job.job_id}_log.txt", mime="text/plain",
    )


def _next_step(code: str | None) -> str:
    return {
        "CorpusError": "**What to try:** choose a different reference library, or "
                       "check that the one you uploaded has a `sequence` column and "
                       "enough entries.",
        "InputError": "**What to try:** check the sequence file — the most common "
                      "causes are the wrong column name or letters other than A/C/G/U.",
        "TargetError": "**What to try:** choose a different target file. If you "
                       "uploaded one, it may be incomplete or edited.",
        "ExternalToolError": "**What to try:** this needs a tool that is not "
                             "installed on this computer. Run without a target, or "
                             "ask for a prepared target file.",
        "Cancelled": "You stopped this analysis. Nothing was saved.",
        "ProcessDied": "**What to try:** run it again. If it keeps happening, send "
                       "the log below to whoever maintains this tool.",
    }.get(code or "", "**What to try:** run it again, and send the log below to "
                      "whoever maintains this tool if it happens twice.")


def recent(runs_dir: Path, limit: int = 20) -> None:
    """A list of past analyses — the tool's memory, kept on disk."""
    found = jobs.discover(runs_dir, limit=limit)
    if not found:
        st.info("No analyses have been run on this computer yet.", icon="ℹ️")
        return

    for job in found:
        state = job.state()
        icon = {"completed": "✅", "failed": "🚫", "running": "⏳",
                "pending": "⏳", "stalled": "⚠️"}.get(state["status"], "•")
        left, right = st.columns([5, 1])
        with left:
            st.markdown(
                f"{icon} **{job.label}** · {_when(job.created_utc)} · "
                f"{state['status']}"
                + (f" · {state['stage_label']}" if state["status"] == "running" else "")
            )
        with right:
            if state["status"] == "completed" and job.artifact_path.exists():
                if st.button("Open", key=f"open_{job.job_id}"):
                    st.session_state["artifact_path"] = str(job.artifact_path)
                    st.session_state["view"] = "Results"
                    st.rerun()
            elif state["status"] in ("running", "pending", "stalled"):
                if st.button("Watch", key=f"watch_{job.job_id}"):
                    st.session_state["watch_job"] = job.job_id
                    st.session_state["view"] = "Progress"
                    st.rerun()


def _when(stamp: str | None) -> str:
    if not stamp:
        return "unknown time"
    try:
        when = datetime.fromisoformat(stamp.replace("Z", "+00:00")).astimezone()
    except ValueError:
        return stamp
    return when.strftime("%d %b %Y, %H:%M")


def _ago(stamp: str | None) -> str:
    if not stamp:
        return "just now"
    try:
        when = datetime.fromisoformat(stamp.replace("Z", "+00:00"))
    except ValueError:
        return "just now"
    seconds = (datetime.now(timezone.utc) - when).total_seconds()
    if seconds < 60:
        return f"{int(seconds)}s ago"
    return f"{int(seconds / 60)} min ago"
