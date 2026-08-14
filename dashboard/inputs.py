"""Discovering and judging the inputs a run needs.

Everything here answers one of two questions: *what can the user choose?* and
*is this combination safe to run at all?*

The second question is the important one. A biologist cannot be expected to
judge when a result is untrustworthy, so the UI refuses outright rather than
warning whenever the system cannot produce the kind of result being asked for.
Warnings are reserved for choices that are legitimate but consequential.
"""

from __future__ import annotations

import hashlib
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

#: User-facing names for the binding modes (refinements §1.1, §5). One name per
#: mode, everywhere — a mode called one thing while choosing it and another
#: while reading the result is the drift §1.1 exists to stop.
BINDING_MODE_LABEL = {
    "pocket": "Pocket/groove recognition",
    "surface": "Surface-patch recognition",
}

#: Shown while *choosing* a mode: short, and about which kind of site this is,
#: because that is the judgement being made.
BINDING_MODE_PREMISE = {
    "pocket": "The aptamer inserts exposed nucleotides into a defined pocket, "
              "cleft, or groove on the protein. Common for recessed binding "
              "sites such as enzyme active sites or ligand-binding pockets.",
    "surface": "The aptamer lies across an extended, relatively shallow patch "
               "of the protein surface, making contacts over a broad interface "
               "rather than entering a discrete pocket. Common for exposed "
               "interaction surfaces such as protein–protein or "
               "protein–nucleic-acid binding sites.",
}

#: Shown while *reading* a result: the full mechanism the comparison assumed,
#: for someone deciding how much to trust the band.
BINDING_MODE_DESCRIPTION = {
    "pocket": "The aptamer presents one or more exposed nucleotides, often from "
              "a loop or bulge, into a defined pocket, cleft, or groove on the "
              "protein surface. Binding is stabilised by close shape "
              "complementarity and specific hydrogen-bonding, stacking, "
              "electrostatic, and van der Waals contacts within the confined "
              "site.",
    "surface": "The aptamer lies across an extended, relatively shallow region "
               "of the protein surface, making multiple contacts over a broad "
               "interface rather than inserting into a discrete pocket. This is "
               "common when the target epitope is a macromolecular interaction "
               "surface, such as a protein–protein, protein–nucleic-acid, or "
               "ligand-binding exosite.",
}


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

    @property
    def original_name(self) -> str:
        """What the user called it, for anything shown back to them."""
        from aptarank.provenance import original_filename

        return original_filename(self.path) or self.path.name

    def describe(self) -> str:
        if self.state == DEVELOPMENT:
            return "Synthetic example data — for testing the software only"
        # A curated library can describe itself in its manifest. The generated
        # line below is a fallback for libraries that do not, not a house style
        # to be overridden by hard-coding one library's text in the UI.
        described = str(self.manifest.get("description") or "").strip()
        if described:
            return described
        counted = f"{self.n_sequences:,} aptamers" if self.n_sequences else "unknown size"
        if self.state == VERIFIED:
            source = self.manifest.get("source", "recorded provenance")
            return f"{counted} · {source}"
        return f"{counted} · provenance not recorded"


@dataclass
class TargetEvidence:
    """A target that has already been prepared and measured, as the UI sees it."""

    path: Path
    pdb_id: str
    chain: str
    synthetic: bool
    binding_mode: str = "pocket"
    structure_kind: str = "experimental"
    d_pocket_A: float | None = None
    patch_area_A2: float | None = None
    n_site_residues: int | None = None
    n_pockets: int | None = None
    selection_method: str | None = None
    shape_warning: bool = False
    was_multi_chain: bool = False
    problem: str | None = None
    name: str | None = None

    @property
    def usable(self) -> bool:
        return self.problem is None

    @property
    def is_predicted(self) -> bool:
        return self.structure_kind == "predicted"

    def describe(self) -> str:
        if self.synthetic:
            return "Fabricated cavity — for testing the software only"
        parts = [BINDING_MODE_LABEL.get(self.binding_mode, self.binding_mode)]
        if self.binding_mode == "surface" and self.patch_area_A2:
            parts.append(
                f"binding site {self.patch_area_A2:.0f} Å² across "
                f"{self.n_site_residues or 0} residues"
            )
        elif self.d_pocket_A:
            parts.append(f"cavity {self.d_pocket_A:.0f} Å across")
        if self.is_predicted:
            parts.append("predicted structure")
        if self.binding_mode == "pocket":
            parts.append(
                "cavity confirmed against known binding-site residues"
                if self.selection_method in ("target_site_overlap", "active_site_overlap")
                else "cavity chosen automatically"
            )
        return " · ".join(parts)


@dataclass
class TargetRequest:
    """What the user asked for in step 3, before anything has been measured.

    Either a target that is already prepared on disk, or a description the
    pipeline will prepare server-side. The user never builds a file by hand.
    """

    kind: str                                   # none | prepared | spec
    label: str = ""
    binding_mode: str = "pocket"
    prepared: TargetEvidence | None = None
    spec: dict[str, Any] = field(default_factory=dict)
    spec_text: str = ""
    problem: str | None = None

    @property
    def usable(self) -> bool:
        return self.problem is None and self.kind != "none"

    @property
    def synthetic(self) -> bool:
        return bool(self.prepared and self.prepared.synthetic)

    @property
    def is_predicted(self) -> bool:
        source = (self.spec.get("tier2", {}).get("target", {}) or {}).get("source")
        return source == "alphafold" or bool(self.prepared and self.prepared.is_predicted)


def _content_key(path: Path) -> str:
    """Identify a file by what is in it, not where it is.

    The same corpus commonly exists in both the code checkout and the data
    directory; listing it twice under one name gives the user two identical
    options and no way to tell them apart.
    """
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(1 << 20):
            digest.update(block)
    return digest.hexdigest()


def discover_libraries(*directories: str | Path) -> list[ReferenceLibrary]:
    """Every distinct reference library on disk, real ones first."""
    found: list[ReferenceLibrary] = []
    seen: set[str] = set()
    for directory in directories:
        base = Path(directory)
        if not base.is_dir():
            continue
        for path in sorted(base.glob("*.csv")):
            try:
                key = _content_key(path)
            except OSError:
                continue
            if key in seen:
                continue          # same file, found in a second location
            seen.add(key)
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
    """Every distinct prepared target bundle on disk, real ones first."""
    found: list[TargetEvidence] = []
    seen: set[str] = set()
    for directory in directories:
        base = Path(directory)
        if not base.is_dir():
            continue
        for path in sorted(base.glob("*.bundle.json")):
            try:
                key = _content_key(path)
            except OSError:
                continue
            if key in seen:
                continue
            seen.add(key)
            found.append(inspect_target(path))
    return sorted(found, key=lambda t: (t.synthetic, t.pdb_id))


def inspect_target(path: str | Path) -> TargetEvidence:
    """Validate a prepared target the same way the pipeline will, and describe it."""
    target = Path(path)
    try:
        from aptarank.tier2 import bundle as bundle_mod

        bundle = bundle_mod.load(target)          # validates schema + checksum
        summary = bundle_mod.summary(bundle)
        pocket = summary.get("selected_pocket") or {}
        patch = summary.get("patch") or {}
        return TargetEvidence(
            path=target,
            pdb_id=summary["pdb_id"],
            chain=summary["chain"],
            synthetic=bool(summary.get("synthetic")),
            binding_mode=summary.get("binding_mode", "pocket"),
            structure_kind=summary.get("structure_kind", "experimental"),
            d_pocket_A=pocket.get("d_pocket_A"),
            patch_area_A2=patch.get("patch_area_A2"),
            n_site_residues=patch.get("n_residues") or len(summary.get("site_residues", [])),
            n_pockets=summary["n_pockets"],
            selection_method=summary.get("pocket_selection"),
            shape_warning=bool(pocket.get("shape_warning") or patch.get("shape_warning")),
            was_multi_chain=bool(summary.get("was_multi_chain")),
            name=summary.get("name"),
        )
    except Exception as exc:  # noqa: BLE001
        return TargetEvidence(
            path=target, pdb_id=target.stem.split("_")[0], chain="?", synthetic=False,
            problem=f"This prepared target failed its integrity check and cannot "
                    f"be used ({type(exc).__name__}: {exc}).",
        )


def build_target_request(
    kind: str,
    prepared: TargetEvidence | None = None,
    spec_text: str = "",
    label: str = "",
) -> TargetRequest:
    """Turn the step-3 choice into something the run can be launched with."""
    if kind == "none":
        return TargetRequest(kind="none", label="No target")
    if kind == "prepared":
        if prepared is None:
            return TargetRequest(kind="prepared", problem="No prepared target chosen.")
        return TargetRequest(
            kind="prepared",
            label=f"{prepared.pdb_id} chain {prepared.chain}",
            binding_mode=prepared.binding_mode,
            prepared=prepared,
            problem=prepared.problem,
        )

    from aptarank.config import parse_target_spec
    from aptarank.errors import ConfigError

    try:
        spec = parse_target_spec(spec_text, label or "this target description")
    except ConfigError as exc:
        return TargetRequest(kind="spec", spec_text=spec_text, problem=str(exc))

    target = spec["tier2"]["target"]
    mode = spec["tier2"].get("binding_mode", "pocket")
    problem = None
    if mode == "surface" and not target.get("target_site_residues"):
        # Surface mode measures the patch those residues define; without them
        # there is nothing to measure, and guessing would invent the site.
        problem = (
            "Surface mode needs binding-site residues: they define the patch "
            "being measured. Add target_site_residues, or choose pocket mode, "
            "where the cavity can be found automatically."
        )
    return TargetRequest(
        kind="spec",
        label=f"{target['id']} chain {target.get('chain') or '(first protein chain)'}",
        binding_mode=mode,
        spec=spec,
        spec_text=spec_text,
        problem=problem,
    )


# -- the gate ------------------------------------------------------------


def eligibility(
    library: ReferenceLibrary | None, target: TargetRequest | None
) -> dict[str, Any]:
    """Will this run be publication-eligible? Answered *before* it starts.

    The rule is the pipeline's, unchanged: a run is publication-eligible only if
    the reference library is real *and* carries a provenance record, and the
    target evidence is real. What changed is when it is evaluated. Promising
    eligibility on the New Analysis page and retracting it on the Results page
    is worse than either verdict alone, because it teaches the reader that the
    badge means nothing.

    Mirrors `artifacts.io.build_artifact`; the reason codes are the same
    strings, so the two surfaces cannot drift apart silently.
    """
    reasons: list[str] = []
    if library is not None and library.is_placeholder:
        reasons.append("placeholder_corpus")
    if library is not None and library.state == UNVERIFIED:
        reasons.append("unverified_corpus_provenance")
    if target is not None and target.synthetic:
        reasons.append("synthetic_target_bundle")
    return {
        "publication_eligible": not reasons and library is not None,
        "development_reasons": reasons,
        "status": "Development only" if reasons else "Publication-eligible",
    }


DEVELOPMENT_REASON_TEXT = {
    "placeholder_corpus": "the reference library is synthetic example data",
    "synthetic_target_bundle": "the target's binding site was fabricated for testing",
    "unverified_corpus_provenance":
        "the reference library has no provenance record, so where it came from "
        "cannot be cited",
}


def review(
    validation: dict[str, Any] | None,
    library: ReferenceLibrary | None,
    target: TargetRequest | None,
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

    if target is not None and target.kind != "none" and target.problem:
        refusals.append(target.problem)

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
                "evidence. Use a real protein target, or switch to Standard "
                "analysis."
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
    if target is None or target.kind == "none":
        warnings.append(
            "No target selected: candidates will be ranked on aptamer-likeness "
            "only, with no aptamer-target compatibility annotation."
        )
    else:
        if target.synthetic:
            warnings.append(
                "Fabricated target binding site: the compatibility annotations are "
                "a software demonstration, not measurements of a real protein."
            )
        if target.is_predicted:
            warnings.append(
                "Predicted structure (AlphaFold): a model, not an experiment."
                + (
                    " A predicted model may not show an interface that only forms "
                    "when a binding partner is present, so surface-mode results "
                    "here need extra caution."
                    if target.binding_mode == "surface" else ""
                )
            )
        prepared = target.prepared
        if prepared and prepared.binding_mode == "pocket" and prepared.selection_method not in (
            "target_site_overlap", "active_site_overlap"
        ):
            warnings.append(
                "The cavity on this target was chosen automatically rather than "
                "confirmed against known binding-site residues."
            )
        if prepared and prepared.shape_warning:
            warnings.append(
                "This binding site is an awkward shape. The comparison assumes a "
                "roughly rounded cavity in pocket mode, and a roughly flat face "
                "in surface mode."
            )
        if prepared and prepared.was_multi_chain:
            warnings.append(
                "This chain came from a multi-chain structure and was measured "
                "after computational removal of the partner chain(s). It retains "
                "the bound-state conformation from the complex, so the "
                "measurement reflects the exposed bound-state interface rather "
                "than a relaxed unbound state."
            )
        if target.kind == "spec":
            warnings.append(
                f"This target has not been prepared yet: the analysis will fetch "
                f"the structure and measure it first, which takes a few minutes "
                f"the first time. You have asserted the binding mode as "
                f"“{BINDING_MODE_LABEL.get(target.binding_mode, target.binding_mode)}” — "
                f"AptaRank does not infer it."
            )
    if preset == "quick":
        warnings.append(
            "Quick preview skips the shuffled controls, so there is no check that "
            "a score reflects structure rather than nucleotide composition."
        )

    verdict = eligibility(library, target)
    return {
        "refusals": refusals,
        "warnings": warnings,
        "can_run": not refusals,
        "expected_status": verdict["status"],
        "development": not verdict["publication_eligible"],
        "development_reasons": verdict["development_reasons"],
    }
