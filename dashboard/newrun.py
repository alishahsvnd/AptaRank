"""The guided "new analysis" flow.

Written for someone who has never opened a terminal. Five steps, in the order
the questions actually arise, with the consequences of each choice stated in
biology rather than in configuration keys.

Two rules this module keeps:

* it never scores anything — validation and execution both go through the same
  `aptarank` command line the pipeline exposes;
* it refuses rather than warns whenever the system cannot produce the kind of
  result being asked for (see `inputs.review`).
"""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd
import streamlit as st
from aptarank.provenance import write_origin

from . import jobs
from .inputs import (
    BINDING_MODE_LABEL,
    BINDING_MODE_PREMISE,
    DEVELOPMENT,
    DEVELOPMENT_REASON_TEXT,
    UNVERIFIED,
    VERIFIED,
    ReferenceLibrary,
    TargetEvidence,
    TargetRequest,
    build_target_request,
    discover_libraries,
    discover_targets,
    inspect_library,
    review,
)

EXAMPLE_CSV = (
    "id,sequence\n"
    "candidate_1,GUUCCAUGGGCCUUGACUUGCUGUGUCAUCACCAUGGGAC\n"
    "candidate_2,CCGCUGUGAGUGUCUCACAGCGGAAUAGGUAC\n"
    "candidate_3,GGUCACUCUCACGGACUGAGAGGUGAGAGUGACCUUAGC\n"
)

STATE_BADGES = {
    VERIFIED: ("✅", "Provenance recorded"),
    UNVERIFIED: ("⚠️", "Provenance not recorded"),
    DEVELOPMENT: ("🚫", "Synthetic — testing only"),
}


def render(repo_root: Path, runs_dir: Path, data_dir: Path | None = None) -> None:
    """Draw the whole flow and, if the user confirms, launch a job.

    `data_dir` is where a deployment keeps its libraries and prepared targets —
    outside the code checkout, which a redeploy replaces.
    """
    data_dir = data_dir or repo_root
    pending = runs_dir / "pending"
    pending.mkdir(parents=True, exist_ok=True)

    st.markdown("### Start a new analysis")
    st.caption(
        "Upload your candidate sequences, choose what to compare them against, "
        "and press Run. Everything stays on this computer."
    )

    candidates_path, validation = _step_sequences(pending, repo_root, data_dir)
    library = _step_library(repo_root, data_dir, pending)
    target = _step_target(repo_root, data_dir, pending)
    preset, name = _step_settings()
    _step_review(runs_dir, candidates_path, validation, library, target, preset, name)


# -- step 1 --------------------------------------------------------------


def _step_sequences(
    pending: Path, repo_root: Path, data_dir: Path
) -> tuple[Path | None, dict[str, Any] | None]:
    st.markdown("#### 1 · Your candidate sequences")
    left, right = st.columns([3, 1])
    with left:
        uploaded = st.file_uploader(
            "RNA sequences to rank",
            type=["csv", "txt", "fasta", "fa", "tsv"],
            help="A .csv with a 'sequence' column (and optionally an 'id' column), "
                 "a plain .txt with one sequence per line, or a FASTA file.",
        )
    with right:
        st.write("")
        st.write("")
        # Two examples: a three-line file that shows the format, and the larger
        # demo batch, which is the one worth actually running - a ranking of
        # three sequences tells you nothing about whether the tool works.
        demo = _demo_candidates(repo_root, data_dir)
        if demo is not None:
            st.download_button(
                "Example file", demo.read_bytes(), file_name="demo_candidates.csv",
                mime="text/csv",
                help="199 mixed-quality sequences to try the tool on. Download, "
                     "then upload it on the left.",
            )
        else:
            st.download_button(
                "Example file", EXAMPLE_CSV, file_name="example_candidates.csv",
                mime="text/csv", help="A small file in the expected format.",
            )

    if uploaded is None:
        st.info(
            "Sequences must use the letters A, C, G and U (T is accepted and read "
            "as U) and be 20–100 letters long."
            + ("  \n\nNo sequences of your own yet? Download the example on the "
               "right and upload it back here." if demo is not None else ""),
            icon="ℹ️",
        )
        return None, None

    path = _stage(pending, uploaded, "candidates")
    validation = _validate(path)

    if not validation.get("ok"):
        st.error(_friendly_error(validation), icon="🚫")
        return path, validation

    cols = st.columns(4)
    cols[0].metric("Sequences read", f"{validation['n_submitted']:,}")
    cols[1].metric("Usable", f"{validation['n_valid']:,}")
    cols[2].metric("Excluded", f"{validation['n_rejected']:,}")
    low, high = validation.get("length_range", [0, 0])
    cols[3].metric("Length range", f"{low}–{high} nt")

    if validation["n_rejected"]:
        with st.expander(f"Why {validation['n_rejected']} row(s) were excluded", expanded=True):
            table = pd.DataFrame(validation["rejections"])
            st.dataframe(table, hide_index=True, use_container_width=True, height=200)
            st.download_button(
                "Download the exclusion report",
                table.to_csv(index=False).encode("utf-8"),
                file_name="excluded_sequences.csv", mime="text/csv",
            )
    return path, validation


# -- step 2 --------------------------------------------------------------


def _step_library(repo_root: Path, data_dir: Path, pending: Path) -> ReferenceLibrary | None:
    st.markdown("#### 2 · Reference library")
    st.caption(
        "AptaRank compares each of your sequences with experimentally validated "
        "RNA aptamers, or your own aptamer reference dataset, so that raw folding "
        "measurements become calibrated scores. This is reference data — it is "
        "not a model of binding."
    )

    libraries = [
        lib for lib in discover_libraries(
            data_dir / "data" / "corpus", data_dir / "data" / "libraries",
            repo_root / "data" / "corpus", repo_root / "data" / "libraries",
        )
        if lib.usable
    ]

    # Radio options are plain strings, not the dataclass instances themselves.
    # Streamlit stores the selected *value* in session state and matches it back
    # against the options on the next rerun; the instances are rebuilt from disk
    # every rerun, so matching them is fragile and the click appears to do
    # nothing. Strings match reliably, and the path keeps two libraries with the
    # same name distinguishable.
    by_key: dict[str, ReferenceLibrary] = {}
    for lib in libraries:
        by_key[f"{_library_label(lib)}  ·  {lib.path.parent}"] = lib
    upload_key = "Upload my own library (.csv)"

    # Never preselect synthetic data: choosing it has to be a decision.
    real = next((i for i, key in enumerate(by_key) if not by_key[key].is_placeholder), None)

    # Whichever library is actually preselected says so. Marking it here rather
    # than in its manifest keeps the label honest if a second library is added
    # and sorts ahead of it.
    if real is not None:
        key = list(by_key)[real]
        library = by_key.pop(key)
        by_key = {
            **{k: v for k, v in list(by_key.items())[:real]},
            _library_label(library, default=True) + f"  ·  {library.path.parent}": library,
            **{k: v for k, v in list(by_key.items())[real:]},
        }

    options = [*by_key, upload_key]
    default_index = real if real is not None else len(options) - 1

    choice = st.radio(
        "Which library should your sequences be compared with?",
        options,
        index=default_index,
        label_visibility="collapsed",
        key="library_choice",
    )

    if choice == upload_key:
        uploaded = st.file_uploader(
            "Reference library (.csv with id, sequence, target_name, target_pdb_id)",
            type=["csv"], key="library_upload",
        )
        if uploaded is None:
            return None
        library = inspect_library(_stage(pending, uploaded, "reference_library"))
    else:
        library = by_key[choice]

    if library.problem:
        st.error(library.problem, icon="🚫")
    elif library.state == DEVELOPMENT:
        st.error(
            "This is synthetic example data. Scores calibrated against it are a "
            "software demonstration and cannot support any scientific claim.",
            icon="🚫",
        )
    elif library.state == UNVERIFIED:
        st.warning(
            "No provenance record accompanies this library. It will still be used, "
            "but note where it came from before quoting results derived from it. "
            "(Add a `<name>.manifest.json` beside it with source, curator and "
            "curated_date to record this.)",
            icon="⚠️",
        )
    return library


def _library_label(library: ReferenceLibrary, default: bool = False) -> str:
    icon, _ = STATE_BADGES[library.state]
    marker = "  (default)" if default else ""
    return f"{icon}  {library.name}{marker}. {library.describe()}"


# -- step 3 --------------------------------------------------------------


def _step_target(repo_root: Path, data_dir: Path, pending: Path) -> TargetRequest:
    st.markdown("#### 3 · Protein target *(optional)*")
    st.caption(
        "AptaRank annotates each candidate with how well its shape agrees with "
        "measurements of the protein target's binding site. Specify a target and "
        "binding mode, and optionally target binding-site residues to choose the "
        "exact site. AptaRank does not predict binding."
    )

    targets = discover_targets(
        data_dir / "cache" / "targets", repo_root / "cache" / "targets"
    )

    # Strings as options, for the same reason as the library picker above.
    none_key = "No target — rank on aptamer-likeness only"
    new_key = "Specify a protein target"
    by_key: dict[str, TargetEvidence] = {}
    for target in targets:
        icon = "🚫" if target.synthetic else ("✅" if target.usable else "⚠️")
        by_key[
            f"{icon}  {target.pdb_id} chain {target.chain} — {target.describe()}"
        ] = target
    options = [none_key, new_key, *by_key]

    choice = st.radio(
        "Compare against a target?", options,
        label_visibility="collapsed", key="target_choice",
    )

    if choice == none_key:
        return build_target_request("none")
    if choice != new_key:
        request = build_target_request("prepared", prepared=by_key[choice])
        _target_notes(request)
        return request

    request = _target_form(pending)
    _target_notes(request)
    return request


def _target_form(pending: Path) -> TargetRequest:
    """Ask for an identifier and a chain; the server does the preparation."""
    left, right = st.columns([2, 1])
    with left:
        source = st.radio(
            "Where does the structure come from?",
            ("pdb", "alphafold"),
            format_func=lambda s: (
                "Experimental structure (PDB)" if s == "pdb"
                else "Predicted model (AlphaFold DB)"
            ),
            horizontal=True,
            key="target_source",
        )
        columns = st.columns([1, 1])
        identifier = columns[0].text_input(
            "PDB ID" if source == "pdb" else "UniProt accession",
            placeholder="7WRQ" if source == "pdb" else "P17936",
            key="target_id",
        ).strip()
        chain = columns[1].text_input(
            "Chain", value="A" if source == "alphafold" else "",
            placeholder="B",
            help="Which chain in the file is your target protein. AlphaFold "
                 "models always have a single chain, A.",
            key="target_chain",
            disabled=source == "alphafold",
        ).strip()
    with right:
        uploaded = st.file_uploader(
            "…or upload a target description", type=["txt", "yaml", "yml"],
            key="target_spec_upload",
            help="A small text file with target_name / source / id / chain / "
                 "binding_mode / target_site_residues.",
        )

    mode = st.radio(
        "How do you believe an aptamer would engage this target?",
        ("pocket", "surface"),
        format_func=lambda m: BINDING_MODE_LABEL[m],
        horizontal=True,
        key="target_binding_mode",
    )
    st.caption(BINDING_MODE_PREMISE[mode])
    st.caption(
        "You assert this; AptaRank does not infer it. The comparison it runs "
        "depends on your answer."
    )

    residues = st.text_input(
        "Binding-site residues" + (" (required for surface mode)" if mode == "surface"
                                   else " (optional)"),
        placeholder="7, 8, 9, 12, 38, 55",
        key="target_residues",
        help="Residue numbers as they are labelled in the structure, for the "
             "chain above. Write 42, not K42 — the letter is redundant, and "
             "these are the depositor's numbers, not positions in the sequence.",
    )
    partner = st.text_input(
        "Binding-partner chains to remove (optional)",
        placeholder="C",
        key="target_partner",
        help="For a complex: the chain(s) sitting on the site you want measured. "
             "They are used to confirm the site, then removed so the surface is "
             "exposed.",
    )

    if uploaded is not None:
        text = uploaded.getvalue().decode("utf-8", errors="replace")
        _stage(pending, uploaded, "target_description")
        return build_target_request("spec", spec_text=text, label=uploaded.name)

    if not identifier:
        return TargetRequest(
            kind="spec",
            problem="Enter a PDB ID or UniProt accession for the protein target.",
        )

    lines = [
        f"source: {source}",
        f"id: {identifier}",
        f"binding_mode: {mode}",
    ]
    if chain:
        lines.append(f"chain: {chain}")
    if partner.strip():
        lines.append(
            "partner_chains: ["
            + ", ".join(c.strip() for c in partner.replace(",", " ").split() if c.strip())
            + "]"
        )
    if residues.strip():
        lines.append(f"target_site_residues: [{_residue_text(residues)}]")
    return build_target_request("spec", spec_text="\n".join(lines), label="this target")


def _residue_text(raw: str) -> str:
    return ", ".join(part for part in raw.replace(",", " ").split() if part)


def _target_notes(request: TargetRequest) -> None:
    if request.problem:
        st.error(request.problem, icon="🚫")
        return
    if request.synthetic:
        st.error(
            "This target's binding site was fabricated for testing. Any "
            "compatibility annotation it produces describes nothing real.",
            icon="🚫",
        )
        return
    if request.kind == "spec":
        target = request.spec["tier2"]["target"]
        st.success(
            f"Will prepare **{target['id']}** "
            f"(chain {target.get('chain') or 'first protein chain'}) for "
            f"**{BINDING_MODE_LABEL[request.binding_mode].lower()}**"
            + (f", using {len(target.get('target_site_residues') or [])} "
               f"binding-site residues" if target.get("target_site_residues") else "")
            + ".",
            icon="✅",
        )


# -- step 4 --------------------------------------------------------------


def _step_settings() -> tuple[str, str]:
    st.markdown("#### 4 · How thorough should the analysis be?")
    preset = st.radio(
        "Analysis depth",
        list(jobs.PRESET_DESCRIPTIONS),
        format_func=lambda key: jobs.PRESET_DESCRIPTIONS[key][0],
        index=1,
        horizontal=True,
        label_visibility="collapsed",
        key="preset_choice",
    )
    st.caption(jobs.PRESET_DESCRIPTIONS[preset][1])
    name = st.text_input(
        "Name for this analysis (optional)",
        placeholder="e.g. NDM-1 first batch",
        key="analysis_name",
    )
    return preset, name


# -- step 5 --------------------------------------------------------------


def _step_review(
    runs_dir: Path,
    candidates_path: Path | None,
    validation: dict[str, Any] | None,
    library: ReferenceLibrary | None,
    target: TargetRequest | None,
    preset: str,
    name: str,
) -> None:
    st.markdown("#### 5 · Review and run")
    verdict = review(validation, library, target, preset)
    has_target = target is not None and target.kind != "none"

    summary = {
        "Sequences to rank": f"{validation['n_valid']:,}" if validation and validation.get("ok") else "—",
        "Reference library": library.name if library else "—",
        "Protein target": target.label if has_target else "None (aptamer-likeness only)",
        "Binding mode": BINDING_MODE_LABEL.get(target.binding_mode, "—") if has_target else "—",
        "Analysis depth": jobs.PRESET_DESCRIPTIONS[preset][0],
    }
    if validation and validation.get("ok"):
        low, high = jobs.estimate_runtime_s(
            validation["n_valid"], preset, with_target=has_target,
            prepare_target=has_target and target.kind == "spec",
        )
        summary["Estimated time"] = f"{_duration(low)} – {_duration(high)}"
    summary["Expected result status"] = verdict["expected_status"]

    st.table(pd.DataFrame({"": summary.keys(), " ": summary.values()}).set_index(""))

    # The eligibility verdict is decided here, from the inputs, and is the same
    # verdict the Results page will show. Promising "publication-eligible" and
    # retracting it afterwards would make the badge worthless.
    if verdict["development"]:
        reasons = verdict["development_reasons"]
        st.markdown(
            "<div class='apt-dev-banner'><b>This will be a development run — "
            "not a result.</b> "
            + ("; ".join(DEVELOPMENT_REASON_TEXT.get(r, r) for r in reasons)
               if reasons else "one or more inputs is missing or synthetic")
            + ". The results will be stamped <code>publication_eligible: false</code> "
              "and must not back any published claim.</div>",
            unsafe_allow_html=True,
        )
    elif verdict["can_run"]:
        st.success(
            "Every input carries provenance: this run will be publication-eligible.",
            icon="✅",
        )

    for message in verdict["refusals"]:
        st.error(message, icon="🚫")
    for message in verdict["warnings"]:
        st.warning(message, icon="⚠️")

    load = jobs.slots(runs_dir)
    if load["running"] >= load["capacity"]:
        st.info(
            f"{load['running']} analysis(es) already running"
            + (f", {load['queued']} waiting" if load["queued"] else "")
            + ". Yours will start automatically as soon as a slot frees — this "
              "machine is shared, so AptaRank only takes a fixed share of it.",
            icon="⏳",
        )

    if st.button(
        "Run analysis",
        type="primary",
        disabled=not verdict["can_run"],
        use_container_width=True,
    ):
        job = jobs.submit(
            runs_dir,
            candidates_path=candidates_path,
            corpus_path=library.path,
            name=name or None,
            corpus_is_placeholder=library.is_placeholder,
            target=target,
            preset=preset,
        )
        st.session_state["watch_job"] = job.job_id
        st.rerun()


# -- helpers -------------------------------------------------------------


def _demo_candidates(repo_root: Path, data_dir: Path) -> Path | None:
    """The demo batch, if this installation has one."""
    for candidate in (data_dir / "demo_candidates.csv",
                      repo_root / "data" / "demo_candidates.csv"):
        if candidate.is_file():
            return candidate
    return None


def _stage(pending: Path, uploaded, stem: str) -> Path:
    """Persist an upload under a name we choose, keyed by its content.

    Streamlit re-runs the script on every interaction; writing by content hash
    means the same upload is staged once and survives those reruns.

    The name the user gave the file is recorded in a sidecar rather than used as
    a path — a client-supplied filename is never safe to build a path from — so
    the provenance panel can show them a name they recognise alongside the hash
    that actually identifies the content.
    """
    data = uploaded.getvalue()
    digest = hashlib.sha256(data).hexdigest()[:12]
    suffix = Path(uploaded.name).suffix.lower() or ".csv"
    if suffix not in (".csv", ".txt", ".tsv", ".fasta", ".fa", ".json", ".yaml", ".yml"):
        suffix = ".csv"
    path = pending / f"{stem}_{digest}{suffix}"
    if not path.exists():
        path.write_bytes(data)
    write_origin(
        path,
        Path(uploaded.name).name,
        staged_utc=datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        sha256=hashlib.sha256(data).hexdigest(),
    )
    return path


@st.cache_data(show_spinner="Checking your sequences…")
def _validate(path: Path) -> dict[str, Any]:
    """Validate through the CLI, so the preview and the run share ingest rules."""
    result = subprocess.run(
        [sys.executable, "-m", "aptarank", "validate-input", str(path)],
        capture_output=True, text=True, shell=False,
    )
    try:
        return json.loads(result.stdout)
    except json.JSONDecodeError:
        return {
            "ok": False,
            "error": "That file could not be read as sequences.",
            "detail": (result.stderr or result.stdout)[-2000:],
        }


def _friendly_error(validation: dict[str, Any]) -> str:
    error = validation.get("error", "")
    if "no valid candidates" in error:
        return (
            "None of the rows in that file could be used. The most common causes "
            "are letters other than A/C/G/U (for example N or a protein sequence), "
            "or sequences outside the 20–100 letter range."
        )
    if "sequence' column" in error:
        return (
            "That CSV has no column called 'sequence'. Rename the column holding "
            "your sequences to 'sequence', or upload a plain .txt with one "
            "sequence per line."
        )
    if "empty" in error:
        return "That file is empty."
    return error or "That file could not be read."


def _duration(seconds: int) -> str:
    if seconds < 90:
        return f"{seconds}s"
    minutes = seconds / 60
    return f"{minutes:.0f} min" if minutes < 90 else f"{minutes / 60:.1f} h"
