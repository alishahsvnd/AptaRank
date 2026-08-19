"""Progress reporting: human-readable for a terminal, JSONL for a UI.

The dashboard does not stream the pipeline's stdout through a pipe. A Streamlit
rerun can stop consuming that pipe, and once it fills the child deadlocks — the
analysis would appear to hang with no way to tell why. Instead the pipeline
writes progress events to a file, and the UI reads that file. A browser refresh,
a closed tab, or a restarted Streamlit server then has no effect on a run in
flight.

Events are versioned and append-only; `status.json` beside them is a compact
snapshot for rendering. The JSONL file is the audit trail.
"""

from __future__ import annotations

import json
import os
import sys
import time
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, TextIO

SCHEMA = "aptarank-progress-v1"

RUN_STARTED = "run_started"
STAGE_STARTED = "stage_started"
ADVANCE = "advance"
STAGE_COMPLETED = "stage_completed"
WARNING = "warning"
RUN_COMPLETED = "run_completed"
RUN_FAILED = "run_failed"

#: Plain-language stage names. A biologist reads these, not "tier1".
STAGE_LABELS = {
    "ingest": "Reading and checking sequences",
    "corpus": "Preparing the reference library",
    "tier1": "Scoring aptamer-likeness",
    "target": "Preparing the protein target",
    "bank": "Calibrating against shuffled controls",
    "tier2": "Checking aptamer-target compatibility",
    "diagrams": "Drawing structure diagrams",
    "artifact": "Writing results",
}


@dataclass
class ProgressReporter:
    """Emits progress as human text, as JSONL to a file, or both."""

    fmt: str = "human"                      # "human" | "jsonl"
    path: str | os.PathLike[str] | None = None
    job_id: str | None = None
    stream: TextIO = field(default_factory=lambda: sys.stderr)
    min_interval_s: float = 0.25            # throttle `advance` events

    _seq: int = field(default=0, init=False)
    _last_emit: dict[str, float] = field(default_factory=dict, init=False)
    _handle: TextIO | None = field(default=None, init=False)

    def __post_init__(self) -> None:
        if self.path:
            target = Path(self.path)
            target.parent.mkdir(parents=True, exist_ok=True)
            self._handle = target.open("a", encoding="utf-8")

    # -- lifecycle -------------------------------------------------------

    def run_started(self, **fields: Any) -> None:
        self.emit(RUN_STARTED, **fields)

    def stage_started(self, stage: str, total: int | None = None, **fields: Any) -> None:
        self.emit(STAGE_STARTED, stage=stage, total=total, **fields)

    def advance(self, stage: str, completed: int, total: int, unit: str = "items") -> None:
        """Throttled: a per-shuffle event would drown the log and the reader."""
        now = time.monotonic()
        final = completed >= total
        if not final and now - self._last_emit.get(stage, 0.0) < self.min_interval_s:
            return
        self._last_emit[stage] = now
        self.emit(ADVANCE, stage=stage, completed=completed, total=total, unit=unit)

    def stage_completed(self, stage: str, **fields: Any) -> None:
        self.emit(STAGE_COMPLETED, stage=stage, **fields)

    def warning(self, message: str, code: str = "warning", **fields: Any) -> None:
        self.emit(WARNING, message=message, code=code, **fields)

    def run_completed(self, artifact_path: str | os.PathLike[str] | None = None,
                      **fields: Any) -> None:
        payload: dict[str, Any] = dict(fields)
        if artifact_path:
            path = Path(artifact_path)
            payload["artifact_path"] = str(path)
            if path.exists():
                from .provenance import sha256_file

                payload["artifact_sha256"] = sha256_file(path)
        self.emit(RUN_COMPLETED, **payload)

    def run_failed(self, code: str, message: str, **fields: Any) -> None:
        self.emit(RUN_FAILED, code=code, message=message, **fields)

    def close(self) -> None:
        if self._handle:
            self._handle.close()
            self._handle = None

    def __enter__(self) -> "ProgressReporter":
        return self

    def __exit__(self, *_exc: Any) -> None:
        self.close()

    # -- emission --------------------------------------------------------

    def emit(self, event: str, **fields: Any) -> None:
        self._seq += 1
        record = {
            "schema": SCHEMA,
            "seq": self._seq,
            "created_utc": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
            "job_id": self.job_id,
            "event": event,
            **{k: v for k, v in fields.items() if v is not None},
        }
        if "stage" in record:
            record.setdefault("message", STAGE_LABELS.get(record["stage"], record["stage"]))

        if self._handle:
            # Flushed every event: a reader tailing the file must never see a
            # run as stalled purely because the writer buffered.
            self._handle.write(json.dumps(record) + "\n")
            self._handle.flush()

        if self.fmt == "human":
            self._print_human(record)

    def _print_human(self, record: dict[str, Any]) -> None:
        event = record["event"]
        if event == ADVANCE:
            stage, done, total = record["stage"], record["completed"], record["total"]
            end = "" if done < total else "\n"
            print(f"\r  {stage:<8} {done}/{total}", end=end, file=self.stream, flush=True)
        elif event == WARNING:
            print(f"\n  !! {record['message']}", file=self.stream, flush=True)
        elif event == RUN_FAILED:
            print(f"\n  error [{record['code']}]: {record['message']}",
                  file=self.stream, flush=True)

    # -- adapter ---------------------------------------------------------

    def callback(self, stage: str, unit: str = "items"):
        """A plain `(i, n)` callable, for code that predates this class."""
        def _advance(i: int, n: int) -> None:
            self.advance(stage, i, n, unit)
        return _advance


def snapshot(events: list[dict[str, Any]]) -> dict[str, Any]:
    """Collapse an event stream into the compact state a UI renders.

    Deliberately reports the *current stage's* own progress rather than a
    single overall percentage: the stages differ in cost by more than an order
    of magnitude, so one blended bar would be a fabricated number.
    """
    state: dict[str, Any] = {
        "status": "pending",
        "stage": None,
        "stage_label": None,
        "completed": None,
        "total": None,
        "unit": None,
        "warnings": [],
        "artifact_path": None,
        "artifact_sha256": None,
        "error_code": None,
        "error_message": None,
        "started_utc": None,
        "updated_utc": None,
        "stages_completed": [],
    }
    for record in events:
        event = record.get("event")
        state["updated_utc"] = record.get("created_utc", state["updated_utc"])
        if event == RUN_STARTED:
            state.update(status="running", started_utc=record.get("created_utc"))
        elif event == STAGE_STARTED:
            state.update(stage=record.get("stage"), stage_label=record.get("message"),
                         completed=0, total=record.get("total"))
        elif event == ADVANCE:
            state.update(stage=record.get("stage"), stage_label=record.get("message"),
                         completed=record.get("completed"), total=record.get("total"),
                         unit=record.get("unit"))
        elif event == STAGE_COMPLETED:
            state["stages_completed"].append(record.get("stage"))
        elif event == WARNING:
            state["warnings"].append(record.get("message"))
        elif event == RUN_COMPLETED:
            state.update(status="completed", artifact_path=record.get("artifact_path"),
                         artifact_sha256=record.get("artifact_sha256"),
                         stage=None, stage_label="Finished")
        elif event == RUN_FAILED:
            state.update(status="failed", error_code=record.get("code"),
                         error_message=record.get("message"))
    return state


def read_events(path: str | os.PathLike[str]) -> list[dict[str, Any]]:
    """Parse a progress file, tolerating a partially-written final line."""
    target = Path(path)
    if not target.exists():
        return []
    events = []
    for line in target.read_text(encoding="utf-8", errors="replace").splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            events.append(json.loads(line))
        except json.JSONDecodeError:
            continue   # the writer is mid-line; the next poll will pick it up
    return events
