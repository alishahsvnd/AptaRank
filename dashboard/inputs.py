"""Discovering and judging the inputs a run needs.

Everything here answers one of two questions: *what can the user choose?* and
*is this combination safe to run at all?*

The second question is the important one. A biologist cannot be expected to
judge when a result is untrustworthy, so the UI refuses outright rather than
warning whenever the system cannot produce the kind of result being asked for.
Warnings are reserved for choices that are legitimate but consequential.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pandas as pd

#: A reference library is only "verified" if its provenance is recorded. Correct
#: columns are not provenance: a CSV does not become validated data because it
#: parses.
MANIFEST_FIELDS = ("source", "curator", "curated_date")

VERIFIED = "verified"
UNVERIFIED = "unverified"
DEVELOPMENT = "development"


@dataclass
class ReferenceLibrary:
    """A corpus of experimentally validated aptamers, as the UI sees it."""

    path: Path
    name: str
    state: str                       # verified | unverified | development
    n_sequences: int | None = None
    manifest: dict[str, Any] = field(default_factory=dict)
    problem: str | None = None

    @property
    def is_placeholder(self) -> bool:
        return self.state == DEVELOPMENT

    @property
    def usable(self) -> bool:
        return self.problem is None

    def describe(self) -> str:
        if self.state == DEVELOPMENT:
            return "Synthetic example data — for testing the software only"
        counted = f"{self.n_sequences:,} sequences" if self.n_sequences else "unknown size"
        if self.state == VERIFIED:
            source = self.manifest.get("source", "recorded provenance")
            return f"{counted} · {source}"
        return f"{counted} · provenance not recorded"


@dataclass
class TargetEvidence:
    """A prepared target bundle, as the UI sees it."""

    path: Path
    pdb_id: str
    chain: str
    synthetic: bool
    d_pocket_A: float | None = None
    n_pockets: int | None = None
    selection_method: str | None = None
    shape_warning: bool = False
    problem: str | None = None
    name: str | None = None

    @property
    def usable(self) -> bool:
        return self.problem is None

    def describe(self) -> str:
        if self.synthetic:
            return "Fabricated cavity — for testing the software only"
        parts = [f"cavity {self.d_pocket_A:.0f} Å across"] if self.d_pocket_A else []
        if self.n_pockets:
            parts.append(f"{self.n_pockets} cavities found")
        if self.selection_method == "active_site_overlap":
            parts.append("cavity confirmed against known active-site residues")
        else:
            parts.append("cavity chosen automatically")
        return " · ".join(parts)


def discover_libraries(*directories: str | Path) -> list[ReferenceLibrary]:
    """Every reference library on disk, real ones first."""
    found: list[ReferenceLibrary] = []
    seen: set[Path] = set()
    for directory in directories:
        base = Path(directory)
        if not base.is_dir():
            continue
        for path in sorted(base.glob("*.csv")):
            resolved = path.resolve()
            if resolved in seen:
                continue
            seen.add(resolved)
            found.append(inspect_library(path))
    # Development data last, and never first in a picker.
    return sorted(found, key=lambda lib: (lib.state == DEVELOPMENT, lib.name))


def inspect_library(path: str | Path) -> ReferenceLibrary:
    """Read just enough of a CSV to describe it, without folding anything."""
    target = Path(path)
    manifest_path = target.with_suffix(".manifest.json")
    manifest: dict[str, Any] = {}
    if manifest_path.is_file():
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            manifest = {}

    problem = None
    n_sequences = None
    try:
        frame = pd.read_csv(target, dtype=str, keep_default_na=False)
        columns = {c.strip().lower() for c in frame.columns}
        if "sequence" not in columns:
            problem = "This file has no 'sequence' column, so it cannot be used " \
                      "as a reference library."
        elif not (columns & {"target_name", "target_pdb_id", "source_reference"}):
            # A bare id/sequence CSV is a candidate file, not a reference
            # library. Offering one as a library would silently calibrate a
            # user's scores against their own unvalidated sequences.
            problem = (
                "This looks like a list of candidate sequences rather than a "
                "reference library — a library also records what each aptamer "
                "was validated against (target_name / target_pdb_id)."
            )
        else:
            n_sequences = len(frame)
            if n_sequences < 50:
                problem = (
                    f"Only {n_sequences} sequences. A reference library needs "
                    f"enough entries to define a distribution — at least a few "
                    f"hundred."
                )
    except Exception as exc:  # noqa: BLE001 - a bad upload must not crash the UI
        problem = f"This file could not be read as a table ({type(exc).__name__})."

    if manifest.get("synthetic") or "placeholder" in target.stem.lower():
        state = DEVELOPMENT
    elif all(manifest.get(f) for f in MANIFEST_FIELDS):
        state = VERIFIED
    else:
        state = UNVERIFIED

    return ReferenceLibrary(
        path=target,
        name=manifest.get("name") or target.stem.replace("_", " "),
        state=state,
        n_sequences=n_sequences,
        manifest=manifest,
        problem=problem,
    )


def discover_targets(*directories: str | Path) -> list[TargetEvidence]:
    """Every prepared target bundle on disk, real ones first."""
    found: list[TargetEvidence] = []
    seen: set[Path] = set()
    for directory in directories:
        base = Path(directory)
        if not base.is_dir():
            continue
        for path in sorted(base.glob("*.bundle.json")):
            resolved = path.resolve()
            if resolved in seen:
                continue
            seen.add(resolved)
            found.append(inspect_target(path))
    return sorted(found, key=lambda t: (t.synthetic, t.pdb_id))


def inspect_target(path: str | Path) -> TargetEvidence:
    """Validate a bundle the same way the pipeline will, and describe it."""
    target = Path(path)
    try:
        from aptarank.tier2 import bundle as bundle_mod

        bundle = bundle_mod.load(target)          # validates schema + checksum
        summary = bundle_mod.summary(bundle)
        pocket = summary["selected_pocket"]
        return TargetEvidence(
            path=target,
            pdb_id=summary["pdb_id"],
            chain=summary["chain"],
            synthetic=bool(summary.get("synthetic")),
            d_pocket_A=pocket["d_pocket_A"],
            n_pockets=summary["n_pockets"],
            selection_method=summary["pocket_selection"],
            shape_warning=bool(pocket.get("shape_warning")),
            name=summary.get("name"),
        )
    except Exception as exc:  # noqa: BLE001
        return TargetEvidence(
            path=target, pdb_id=target.stem.split("_")[0], chain="?", synthetic=False,
            problem=f"This target file failed its integrity check and cannot be "
                    f"used ({type(exc).__name__}: {exc}).",
        )


# -- the gate ------------------------------------------------------------


def review(
    validation: dict[str, Any] | None,
    library: ReferenceLibrary | None,
    target: TargetEvidence | None,
    preset: str,
) -> dict[str, Any]:
    """Decide whether this combination may run, and what must be said about it.

    Returns `refusals` (hard stops) and `warnings` (proceed, but unmissably).
    """
    refusals: list[str] = []
    warnings: list[str] = []

    if not validation or not validation.get("ok"):
        refusals.append(
            validation.get("error", "No sequence file has been checked yet.")
            if validation else "Upload a file of candidate sequences first."
        )
    elif not validation.get("n_valid"):
        refusals.append(
            "None of the sequences in that file could be used. Check the "
            "rejection reasons above."
        )

    if library is None:
        refusals.append("Choose a reference library.")
    elif not library.usable:
        refusals.append(library.problem or "That reference library cannot be used.")

    if target is not None and not target.usable:
        refusals.append(target.problem or "That target file cannot be used.")

    # A rigorous preset promises statistics. Synthetic inputs cannot deliver
    # them, and quietly downgrading the preset would hide that.
    if preset == "evaluation":
        if library is not None and library.is_placeholder:
            refusals.append(
                "The rigorous setting is for numbers you intend to publish, and "
                "this is synthetic example data. Choose a real reference library, "
                "or switch to Standard analysis for a development run."
            )
        if target is not None and target.synthetic:
            refusals.append(
                "The rigorous setting cannot be used with fabricated target "
                "evidence. Use a prepared target from a real structure, or "
                "switch to Standard analysis."
            )

    # Warnings: legitimate choices with consequences the user must see.
    if validation and validation.get("n_rejected"):
        warnings.append(
            f"{validation['n_rejected']} of {validation['n_submitted']} rows will be "
            f"excluded. They are listed with reasons above and in the results."
        )
    if library is not None and library.is_placeholder:
        warnings.append(
            "Synthetic reference library: every score is calibrated against "
            "made-up data. The results are a software demonstration, not findings."
        )
    elif library is not None and library.state == UNVERIFIED:
        warnings.append(
            "This reference library has no provenance record. Note where it came "
            "from before quoting anything derived from it."
        )
    if target is None:
        warnings.append(
            "No target selected: candidates will be ranked on intrinsic structure "
            "only, with no target-aware annotation."
        )
    else:
        if target.synthetic:
            warnings.append(
                "Fabricated target cavity: the geometric annotations are a "
                "software demonstration, not measurements of a real protein."
            )
        if target.selection_method and target.selection_method != "active_site_overlap":
            warnings.append(
                "The cavity on this target was chosen automatically rather than "
                "confirmed against known active-site residues."
            )
        if target.shape_warning:
            warnings.append(
                "This cavity is an awkward shape, and the size comparison assumes "
                "a roughly rounded pocket."
            )
    if preset == "quick":
        warnings.append(
            "Quick preview skips the shuffled controls, so there is no check that "
            "a score reflects structure rather than nucleotide composition."
        )

    development = bool(
        (library is not None and library.is_placeholder)
        or (target is not None and target.synthetic)
    )
    return {
        "refusals": refusals,
        "warnings": warnings,
        "can_run": not refusals,
        "expected_status": "Development only" if development else "Publication-eligible",
        "development": development,
    }
