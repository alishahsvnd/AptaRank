"""Candidate ingest and validation (spec §4.3).

Reads a `.txt` (one sequence per line), `.csv` (column `sequence`, optional
`id`) or FASTA file, normalises to RNA, validates, and deduplicates.

Rejected rows are never silently dropped: every rejection is returned with a
row number and a human-readable reason, and ends up in the run artifact.
"""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterator

import pandas as pd

from .errors import InputError

RNA_ALPHABET = frozenset("ACGU")
_WHITESPACE = re.compile(r"\s+")


@dataclass
class RawRecord:
    """One row as it appeared in the input file, before normalisation."""

    row: int
    sequence: str
    supplied_id: str | None = None


@dataclass
class IngestResult:
    """Clean candidate table plus a full account of what was thrown away."""

    candidates: pd.DataFrame  # candidate_id, sequence, length, duplicate_count, source_ids
    rejections: list[dict[str, Any]] = field(default_factory=list)
    n_submitted: int = 0
    filename: str | None = None

    @property
    def n_valid(self) -> int:
        return len(self.candidates)

    @property
    def n_rejected(self) -> int:
        return len(self.rejections)

    def summary(self) -> dict[str, Any]:
        return {
            "filename": self.filename,
            "n_submitted": self.n_submitted,
            "n_valid": self.n_valid,
            "n_rejected": self.n_rejected,
            "n_unique": int(self.candidates["duplicate_count"].eq(1).sum())
            if self.n_valid
            else 0,
            "rejections": self.rejections,
        }


def normalise_sequence(raw: str) -> str:
    """Strip whitespace, uppercase, and write DNA-style T as RNA-style U.

    The same molecule is written with T in DNA-style files and U in RNA-style
    files; this is a notation difference, not a data difference.
    """
    return _WHITESPACE.sub("", raw).upper().replace("T", "U")


def validate_sequence(seq: str, min_length: int, max_length: int) -> str | None:
    """Return a rejection reason, or None if the sequence is usable."""
    if not seq:
        return "empty sequence"
    bad = sorted(set(seq) - RNA_ALPHABET)
    if bad:
        chars = ", ".join(f"'{c}'" for c in bad[:4])
        return f"invalid character{'s' if len(bad) > 1 else ''} {chars}"
    if len(seq) < min_length:
        return f"length {len(seq)} below minimum {min_length}"
    if len(seq) > max_length:
        return f"length {len(seq)} above maximum {max_length}"
    return None


def read_records(path: str | Path) -> tuple[list[RawRecord], str]:
    """Read raw records from txt / csv / fasta. Returns (records, format)."""
    p = Path(path)
    if not p.exists():
        raise InputError(f"candidate file not found: {p}")

    suffix = p.suffix.lower()
    text = p.read_text(encoding="utf-8-sig")
    if not text.strip():
        raise InputError(f"candidate file is empty: {p}")

    if suffix == ".csv" or (suffix in {".tsv", ".txt"} and _looks_like_csv(text)):
        return list(_read_csv(text, p)), "csv"
    if text.lstrip().startswith(">"):
        return list(_read_fasta(text)), "fasta"
    return list(_read_lines(text)), "txt"


def ingest(
    path: str | Path,
    min_length: int,
    max_length: int,
    id_prefix: str = "c",
) -> IngestResult:
    """Full §4.3 pipeline: read, normalise, validate, deduplicate, assign ids."""
    records, _fmt = read_records(path)
    rejections: list[dict[str, Any]] = []

    # sequence -> first occurrence bookkeeping, preserving input order
    seen: dict[str, dict[str, Any]] = {}
    for rec in records:
        seq = normalise_sequence(rec.sequence)
        reason = validate_sequence(seq, min_length, max_length)
        if reason is not None:
            rejections.append(
                {
                    "row": rec.row,
                    "id": rec.supplied_id,
                    "sequence": _truncate(rec.sequence),
                    "reason": reason,
                }
            )
            continue
        if seq in seen:
            seen[seq]["duplicate_count"] += 1
            seen[seq]["source_ids"].append(rec.supplied_id or f"row{rec.row}")
        else:
            seen[seq] = {
                "sequence": seq,
                "length": len(seq),
                "duplicate_count": 1,
                "source_ids": [rec.supplied_id or f"row{rec.row}"],
                "first_row": rec.row,
                "supplied_id": rec.supplied_id,
            }

    rows = list(seen.values())
    width = max(5, len(str(len(rows))))
    supplied = [r["supplied_id"] for r in rows]
    use_supplied = all(s for s in supplied) and len(set(supplied)) == len(supplied)
    for i, row in enumerate(rows, start=1):
        row["candidate_id"] = (
            str(row["supplied_id"]) if use_supplied else f"{id_prefix}{i:0{width}d}"
        )

    columns = ["candidate_id", "sequence", "length", "duplicate_count", "source_ids"]
    frame = pd.DataFrame(rows, columns=columns) if rows else pd.DataFrame(columns=columns)

    if frame.empty:
        raise InputError(
            f"no valid candidates in {path}: "
            f"{len(records)} rows read, all rejected. "
            f"First reason: {rejections[0]['reason'] if rejections else 'unknown'}"
        )

    return IngestResult(
        candidates=frame.reset_index(drop=True),
        rejections=rejections,
        n_submitted=len(records),
        filename=str(path),
    )


# -- format readers ------------------------------------------------------


def _looks_like_csv(text: str) -> bool:
    head = text.lstrip().splitlines()[0].lower() if text.strip() else ""
    return "," in head and "sequence" in head


def _read_csv(text: str, path: Path) -> Iterator[RawRecord]:
    reader = csv.DictReader(text.splitlines())
    fieldnames = [f.strip().lower() for f in (reader.fieldnames or [])]
    if "sequence" not in fieldnames:
        raise InputError(
            f"{path}: CSV input requires a 'sequence' column, found {reader.fieldnames}"
        )
    lookup = {f.strip().lower(): f for f in (reader.fieldnames or [])}
    seq_col = lookup["sequence"]
    id_col = lookup.get("id")
    for row_no, row in enumerate(reader, start=2):  # row 1 is the header
        value = (row.get(seq_col) or "").strip()
        if not value:
            continue
        supplied = (row.get(id_col) or "").strip() if id_col else ""
        yield RawRecord(row=row_no, sequence=value, supplied_id=supplied or None)


def _read_lines(text: str) -> Iterator[RawRecord]:
    for row_no, line in enumerate(text.splitlines(), start=1):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        yield RawRecord(row=row_no, sequence=stripped)


def _read_fasta(text: str) -> Iterator[RawRecord]:
    header: str | None = None
    chunks: list[str] = []
    start_row = 0
    for row_no, line in enumerate(text.splitlines(), start=1):
        stripped = line.strip()
        if stripped.startswith(">"):
            if header is not None or chunks:
                yield RawRecord(row=start_row, sequence="".join(chunks), supplied_id=header)
            header = stripped[1:].split()[0] if len(stripped) > 1 else None
            chunks = []
            start_row = row_no
        elif stripped:
            chunks.append(stripped)
    if chunks:
        yield RawRecord(row=start_row, sequence="".join(chunks), supplied_id=header)


def _truncate(seq: str, limit: int = 60) -> str:
    seq = seq.strip()
    return seq if len(seq) <= limit else seq[: limit - 1] + "…"
