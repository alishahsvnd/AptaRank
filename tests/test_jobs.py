"""The job runner is what stands between a biologist and a hung browser tab."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

from dashboard import jobs


def make_job(tmp_path: Path, job_id: str = "job_test", **request) -> jobs.Job:
    directory = jobs.jobs_root(tmp_path) / job_id
    (directory / "inputs").mkdir(parents=True, exist_ok=True)
    payload = {"job_id": job_id, "name": None, "created_utc": "2026-07-26T00:00:00Z",
               "pid": None, "command": [], **request}
    (directory / "request.json").write_text(json.dumps(payload), encoding="utf-8")
    return jobs.Job(job_id=job_id, directory=directory, request=payload)


def write_events(job: jobs.Job, *events: dict) -> None:
    with job.progress_file.open("a", encoding="utf-8") as handle:
        for event in events:
            handle.write(json.dumps({"schema": "aptarank-progress-v1", **event}) + "\n")


def test_a_job_with_no_pid_is_queued_not_failed(tmp_path):
    """A queued job has no process yet; that must not read as a crash."""
    job = make_job(tmp_path)
    state = job.state()
    assert state["status"] == jobs.QUEUED
    assert "Waiting" in state["stage_label"]


def test_a_dead_process_that_never_completed_is_a_failure(tmp_path):
    job = make_job(tmp_path, pid=999999)          # a PID that is not running
    write_events(job, {"event": "run_started"}, {"event": "stage_started", "stage": "tier1"})
    state = job.state()
    assert state["status"] == "failed"
    assert state["error_code"] == "ProcessDied"
    assert "stopped unexpectedly" in state["error_message"]


def test_success_without_an_artifact_is_a_failure(tmp_path):
    """Exit zero is not proof of a result."""
    job = make_job(tmp_path, pid=999999)
    write_events(job, {"event": "run_started"}, {"event": "run_completed"})
    state = job.state()
    assert state["status"] == "failed"
    assert state["error_code"] == "MissingArtifact"


def test_success_with_an_artifact_is_reported_as_complete(tmp_path):
    job = make_job(tmp_path, pid=999999)
    job.artifact_path.write_text('{"run_id": "x"}', encoding="utf-8")
    write_events(
        job,
        {"event": "run_started"},
        {"event": "run_completed", "artifact_path": str(job.artifact_path)},
    )
    assert job.state()["status"] == "completed"


def test_discovery_reads_jobs_from_disk_not_from_a_session(tmp_path):
    """The UI's memory lives on disk, so a browser refresh cannot lose a run."""
    make_job(tmp_path, "job_a")
    make_job(tmp_path, "job_b")
    found = {job.job_id for job in jobs.discover(tmp_path)}
    assert found == {"job_a", "job_b"}


def test_worker_cap_overrides_whatever_the_preset_asks_for(tmp_path, monkeypatch):
    """Server resource policy must win over user or preset configuration."""
    monkeypatch.setattr(jobs, "WORKERS_PER_JOB", "4")
    monkeypatch.setattr(jobs, "MAX_CONCURRENT_JOBS", 0)   # keep it queued, don't spawn

    job = jobs.submit(
        tmp_path,
        candidates_path=tmp_path / "c.csv",
        corpus_path=tmp_path / "corpus.csv",
        preset="standard",
        extra_sets=["tier1.parallel.workers=64"],
    )
    command = job.request["command"]
    workers = [command[i + 1] for i, a in enumerate(command)
               if a == "--set" and command[i + 1].startswith("tier1.parallel.workers")]
    assert workers == ["tier1.parallel.workers=4"]


def test_submission_beyond_capacity_queues_rather_than_oversubscribing(tmp_path, monkeypatch):
    monkeypatch.setattr(jobs, "MAX_CONCURRENT_JOBS", 0)
    job = jobs.submit(
        tmp_path, candidates_path=tmp_path / "c.csv", corpus_path=tmp_path / "corpus.csv"
    )
    assert job.request["pid"] is None
    assert job.state()["status"] == jobs.QUEUED
    assert jobs.slots(tmp_path)["queued"] == 1


def test_the_command_uses_this_interpreter_never_python_from_path(tmp_path, monkeypatch):
    monkeypatch.setattr(jobs, "MAX_CONCURRENT_JOBS", 0)
    job = jobs.submit(
        tmp_path, candidates_path=tmp_path / "c.csv", corpus_path=tmp_path / "corpus.csv"
    )
    assert job.request["command"][0] == sys.executable
    assert job.request["command"][1:4] == ["-m", "aptarank", "run"]


def test_placeholder_corpora_are_flagged_on_the_command_line(tmp_path, monkeypatch):
    """The pipeline, not the UI, decides publication eligibility."""
    monkeypatch.setattr(jobs, "MAX_CONCURRENT_JOBS", 0)
    job = jobs.submit(
        tmp_path, candidates_path=tmp_path / "c.csv", corpus_path=tmp_path / "corpus.csv",
        corpus_is_placeholder=True,
    )
    assert "--development-corpus" in job.request["command"]
    assert "--corpus" not in job.request["command"]


def test_uploads_are_never_written_under_a_client_supplied_path(tmp_path):
    class Upload:
        name = "../../evil.csv"

        def getvalue(self):
            return b"id,sequence\n"

    directory = tmp_path / "job"
    path = jobs.stage_upload(directory, Upload(), "candidates.csv")
    assert path.parent == directory / "inputs"
    assert ".." not in str(path.relative_to(directory))


def test_runtime_estimate_is_a_range_not_a_promise():
    low, high = jobs.estimate_runtime_s(1000, "standard", with_target=True)
    assert low < high
    assert jobs.estimate_runtime_s(1000, "quick", False)[1] < high


@pytest.mark.parametrize("preset", ["quick", "standard", "evaluation"])
def test_every_preset_is_described_to_the_user(preset):
    title, explanation = jobs.PRESET_DESCRIPTIONS[preset]
    assert title and explanation
    assert jobs.preset_settings(preset)
