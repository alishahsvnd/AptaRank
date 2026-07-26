"""Progress events are what the dashboard sees. If they lie, the UI lies."""

from __future__ import annotations

import json

from aptarank.progress import (
    ProgressReporter,
    read_events,
    snapshot,
)


def test_events_are_written_as_one_json_object_per_line(tmp_path):
    path = tmp_path / "progress.jsonl"
    with ProgressReporter(fmt="jsonl", path=path, job_id="j1") as reporter:
        reporter.run_started(candidates="x.csv")
        reporter.stage_started("tier1", total=10)
        reporter.advance("tier1", 10, 10)
        reporter.run_completed()

    lines = path.read_text(encoding="utf-8").strip().splitlines()
    assert len(lines) == 4
    records = [json.loads(line) for line in lines]
    assert [r["seq"] for r in records] == [1, 2, 3, 4]
    assert all(r["job_id"] == "j1" for r in records)


def test_stage_events_carry_a_plain_language_label(tmp_path):
    path = tmp_path / "p.jsonl"
    with ProgressReporter(fmt="jsonl", path=path) as reporter:
        reporter.stage_started("tier1")
    record = json.loads(path.read_text(encoding="utf-8").strip())
    assert record["message"] == "Scoring intrinsic structure"


def test_advance_is_throttled_but_never_drops_the_final_update(tmp_path):
    path = tmp_path / "p.jsonl"
    with ProgressReporter(fmt="jsonl", path=path, min_interval_s=60) as reporter:
        for i in range(1, 51):
            reporter.advance("tier1", i, 50)

    records = [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines()]
    assert len(records) < 50                      # throttled
    assert records[-1]["completed"] == 50         # but the last one always lands


def test_snapshot_reports_the_current_stage_not_a_blended_percentage(tmp_path):
    path = tmp_path / "p.jsonl"
    with ProgressReporter(fmt="jsonl", path=path) as reporter:
        reporter.run_started()
        reporter.stage_started("ingest")
        reporter.stage_completed("ingest")
        reporter.stage_started("tier1", total=200)
        reporter.advance("tier1", 120, 200)

    state = snapshot(read_events(path))
    assert state["status"] == "running"
    assert state["stage"] == "tier1"
    assert (state["completed"], state["total"]) == (120, 200)
    assert state["stages_completed"] == ["ingest"]


def test_snapshot_surfaces_warnings_and_failure_reasons(tmp_path):
    path = tmp_path / "p.jsonl"
    with ProgressReporter(fmt="jsonl", path=path) as reporter:
        reporter.run_started()
        reporter.warning("Placeholder reference library", code="placeholder_corpus")
        reporter.run_failed("CorpusError", "no reference corpus configured")

    state = snapshot(read_events(path))
    assert state["status"] == "failed"
    assert state["error_code"] == "CorpusError"
    assert "Placeholder reference library" in state["warnings"]


def test_completion_records_the_artifact_and_its_hash(tmp_path):
    artifact = tmp_path / "artifact.json"
    artifact.write_text('{"run_id": "x"}', encoding="utf-8")
    path = tmp_path / "p.jsonl"
    with ProgressReporter(fmt="jsonl", path=path) as reporter:
        reporter.run_completed(artifact_path=artifact)

    state = snapshot(read_events(path))
    assert state["artifact_path"] == str(artifact)
    assert len(state["artifact_sha256"]) == 64


def test_a_half_written_final_line_is_ignored_not_fatal(tmp_path):
    """The reader polls a file the writer is still appending to."""
    path = tmp_path / "p.jsonl"
    with ProgressReporter(fmt="jsonl", path=path) as reporter:
        reporter.run_started()
    with path.open("a", encoding="utf-8") as handle:
        handle.write('{"event": "advance", "compl')   # torn write

    events = read_events(path)
    assert len(events) == 1
    assert snapshot(events)["status"] == "running"
