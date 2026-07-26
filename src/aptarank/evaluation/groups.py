"""Comparison groups for evaluation (spec §8.1).

Four groups, not the spec's three. The extra one is a harder negative control,
added because an IID random control is too easy: matching only single-letter
frequencies lets the tool "win" on composition rather than on structure, which
is precisely the criticism the shuffled-control check exists to answer.

    validated       real aptamers from the reference corpus (positive control)
    random          IID sampling at corpus letter frequencies, corpus lengths
    shuffled        dinucleotide shuffles of the validated aptamers (hard negative)
    generated       foundation-model output — the actual use case, unlabelled
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd

from ..errors import InputError
from ..ingest import normalise_sequence, validate_sequence
from ..provenance import derive_seed
from ..tier1 import shuffles

GROUP_NAMES = ("validated", "random", "shuffled", "generated")


@dataclass
class ComparisonGroups:
    """Sequences per group, plus how each group was constructed."""

    frames: dict[str, pd.DataFrame] = field(default_factory=dict)
    provenance: dict[str, Any] = field(default_factory=dict)

    def __getitem__(self, name: str) -> pd.DataFrame:
        return self.frames[name]

    def items(self):
        return self.frames.items()

    def summary(self) -> dict[str, Any]:
        return {
            "sizes": {name: len(frame) for name, frame in self.frames.items()},
            "provenance": self.provenance,
        }


def build_groups(
    corpus: pd.DataFrame,
    seed: int,
    generated_path: str | Path | None = None,
    n_random: int | None = None,
    min_length: int = 20,
    max_length: int = 100,
) -> ComparisonGroups:
    """Assemble every comparison group from one corpus and one optional file."""
    validated = corpus[["candidate_id", "sequence"]].copy()
    validated["group"] = "validated"
    if "target_pdb_id" in corpus.columns:
        validated["target_pdb_id"] = corpus["target_pdb_id"].to_numpy()
    if "target_name" in corpus.columns:
        validated["target_name"] = corpus["target_name"].to_numpy()

    n = n_random or len(validated)
    random_frame = random_rna_matched(validated["sequence"], n, seed)
    shuffled_frame = shuffled_negatives(validated, seed)

    frames = {
        "validated": validated,
        "random": random_frame,
        "shuffled": shuffled_frame,
    }

    provenance: dict[str, Any] = {
        "seed": seed,
        "random": "IID sampling at corpus single-letter frequencies, "
                  "lengths drawn from the corpus length distribution",
        "shuffled": "dinucleotide-preserving shuffles of the validated group",
    }

    if generated_path:
        frames["generated"] = load_generated(generated_path, min_length, max_length)
        provenance["generated"] = str(generated_path)
        overlap = set(frames["generated"]["sequence"]) & set(validated["sequence"])
        if overlap:
            # Memorised training sequences would otherwise flatter the
            # generated group with the positive group's own scores.
            frames["generated"] = frames["generated"][
                ~frames["generated"]["sequence"].isin(overlap)
            ].reset_index(drop=True)
        provenance["generated_dropped_as_corpus_duplicates"] = len(overlap)

    return ComparisonGroups(frames=frames, provenance=provenance)


def random_rna_matched(
    reference_sequences: Sequence[str], n: int, seed: int
) -> pd.DataFrame:
    """IID RNA matched to the corpus's letter frequencies and length distribution.

    Matching matters: an unmatched control would let the tool win trivially on
    composition rather than on structure, which is not the claim we want to make.
    """
    rng = np.random.default_rng(derive_seed(seed, "random_rna"))
    lengths = np.array([len(s) for s in reference_sequences])
    joined = "".join(reference_sequences)
    letters = np.array(list("ACGU"))
    freqs = np.array([joined.count(c) for c in letters], dtype=float)
    freqs = freqs / freqs.sum()

    rows = []
    for i in range(n):
        length = int(rng.choice(lengths))
        rows.append(
            {
                "candidate_id": f"rand{i:05d}",
                "sequence": "".join(rng.choice(letters, size=length, p=freqs)),
                "group": "random",
            }
        )
    return pd.DataFrame(rows)


def shuffled_negatives(validated: pd.DataFrame, seed: int) -> pd.DataFrame:
    """One dinucleotide shuffle per validated aptamer — the hard negative."""
    rows = []
    for row in validated.itertuples():
        control_seed = derive_seed(seed, row.candidate_id, "eval_shuffle")
        shuffled = shuffles.generate_shuffles(row.sequence, 1, 2, control_seed)[0]
        rows.append(
            {
                "candidate_id": f"shuf_{row.candidate_id}",
                "sequence": shuffled,
                "group": "shuffled",
                "source_id": row.candidate_id,
            }
        )
    return pd.DataFrame(rows)


def load_generated(path: str | Path, min_length: int, max_length: int) -> pd.DataFrame:
    """Read foundation-model output, normalised and validated like any input."""
    p = Path(path)
    if not p.exists():
        raise InputError(f"generated-sequence file not found: {p}")

    if p.suffix.lower() == ".csv":
        frame = pd.read_csv(p, dtype=str, keep_default_na=False)
        frame.columns = [c.strip().lower() for c in frame.columns]
        if "sequence" not in frame.columns:
            raise InputError(f"{p}: expected a 'sequence' column")
        raw = frame["sequence"].tolist()
        ids = frame["id"].tolist() if "id" in frame.columns else None
    else:
        raw = [line.strip() for line in p.read_text(encoding="utf-8").splitlines()
               if line.strip() and not line.startswith(("#", ">"))]
        ids = None

    rows = []
    for i, value in enumerate(raw):
        seq = normalise_sequence(value)
        if validate_sequence(seq, min_length, max_length) is not None:
            continue
        rows.append(
            {
                "candidate_id": (ids[i] if ids else f"gen{i:05d}"),
                "sequence": seq,
                "group": "generated",
            }
        )
    if not rows:
        raise InputError(f"{p}: no valid generated sequences")
    return pd.DataFrame(rows).drop_duplicates("sequence").reset_index(drop=True)


def canonical_target(corpus: pd.DataFrame) -> tuple[pd.Series, str]:
    """One key per *biological* target, not per structure.

    Two PDB entries can be the same protein, and grouping by PDB ID alone would
    let aptamers against the same protein land in both the reference corpus and
    the held-out fold — exactly the leak the fold split exists to prevent. The
    target name is therefore preferred as the grouping key when present, with
    the PDB ID as a fallback.
    """
    def clean(series: pd.Series) -> pd.Series:
        return (
            series.fillna("").astype(str).str.strip().str.upper()
            .str.replace(r"[^A-Z0-9]+", "-", regex=True).str.strip("-")
        )

    name = clean(corpus["target_name"]) if "target_name" in corpus.columns else None
    pdb = clean(corpus["target_pdb_id"]) if "target_pdb_id" in corpus.columns else None

    if name is not None and (name != "").any():
        key = name.where(name != "", pdb if pdb is not None else "")
        return key.replace("", "__unlabelled__"), "target_name"
    if pdb is not None and (pdb != "").any():
        return pdb.replace("", "__unlabelled__"), "target_pdb_id"
    return pd.Series(["__unlabelled__"] * len(corpus), index=corpus.index), "none"


def target_folds(
    corpus: pd.DataFrame, n_folds: int, seed: int
) -> list[dict[str, Any]]:
    """Split the corpus into folds that never split a target across folds.

    Grouping by target rather than by sequence is what makes E1 honest: two
    aptamers selected against the same protein are often near-relatives, so a
    random split would leave a held-out sequence's own family in the reference
    distribution it is being scored against.
    """
    labels, key = canonical_target(corpus)

    if key == "none":
        # No labels at all: fall back to a sequence-level split and say so
        # loudly. This is a sensitivity analysis, not a headline result.
        rng = np.random.default_rng(derive_seed(seed, "folds"))
        assignment = rng.integers(0, n_folds, len(corpus))
        return [
            {"fold": k, "grouping": "sequence (no target labels available)",
             "targets": [], "test_index": np.flatnonzero(assignment == k)}
            for k in range(n_folds)
        ]

    unique = sorted(labels.unique())
    rng = np.random.default_rng(derive_seed(seed, "folds"))
    order = rng.permutation(len(unique))
    assignment = {unique[i]: int(order[i] % n_folds) for i in range(len(unique))}

    folds = []
    for k in range(n_folds):
        members = labels.map(assignment) == k
        folds.append(
            {
                "fold": k,
                "grouping": key,
                "targets": sorted({t for t, a in assignment.items() if a == k}),
                "test_index": np.flatnonzero(members.to_numpy()),
            }
        )
    return [f for f in folds if len(f["test_index"]) > 0]
