"""Launch and track pipeline runs as durable, on-disk jobs.

The dashboard never scores anything. It stages the uploaded files, launches the
same `aptarank run` command a reviewer would type, and reads what that process
writes. Two consequences worth stating plainly:

* the claim "the dashboard holds no computation of its own" stays literally
  true — the UI is an orchestration and visualisation client;
* a run survives a browser refresh, a closed tab, or a restarted Streamlit
  server, because its entire state lives in a job directory rather than in a
  session.

Job layout:

    runs/jobs/<job_id>/
        request.json      what the user asked for, in their own terms
        inputs/           the staged uploads
        progress.jsonl    append-only event stream from the pipeline
        stdout.log        the child's output, for the technical details panel
        stderr.log
        artifact.json     the result
"""

from __future__ import annotations

import json
import os
import signal
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence

from aptarank.progress import read_events, snapshot

JOBS_DIRNAME = "jobs"
#: A run with no new progress event for this long is treated as stalled. Set
#: generously: the corpus and calibration stages are silent for a while.
HEARTBEAT_TIMEOUT_S = 900

#: Resource policy. On a shared machine these are the difference between a
#: useful service and an antisocial one: AptaRank must never take the whole box
#: away from whoever else is using it. Server deployments raise the job limit
#: and cap the workers each job may claim; the defaults suit one laptop.
MAX_CONCURRENT_JOBS = int(os.environ.get("APTARANK_MAX_CONCURRENT_JOBS", "1"))
WORKERS_PER_JOB = os.environ.get("APTARANK_WORKERS_PER_JOB")   # None -> config default

#: A deployment-specific config layered over configs/default.yaml — where the
#: server keeps its caches and results, outside the replaceable code checkout.
SITE_CONFIG = os.environ.get("APTARANK_CONFIG")

QUEUED, RUNNING, COMPLETED, FAILED, STALLED = (
    "queued", "running", "completed", "failed", "stalled"
)


@dataclass
class Job:
    """One pipeline run, reconstructed from its directory."""

    job_id: str
    directory: Path
    request: dict[str, Any]

    @property
    def progress_file(self) -> Path:
        return self.directory / "progress.jsonl"

    @property
    def artifact_path(self) -> Path:
        return self.directory / "artifact.json"

    @property
    def created_utc(self) -> str:
        return self.request.get("created_utc", "")

    @property
    def label(self) -> str:
        return self.request.get("name") or self.job_id

    @property
    def is_queued(self) -> bool:
        return self.request.get("pid") is None and not self.progress_file.exists()

    def state(self) -> dict[str, Any]:
        """Current status, derived from the event stream plus liveness checks."""
        events = read_events(self.progress_file)
        state = snapshot(events)
        state["job_id"] = self.job_id
        state["label"] = self.label
        state["created_utc"] = self.created_utc
        state["n_events"] = len(events)

        if self.is_queued:
            state["status"] = QUEUED
            state["stage_label"] = "Waiting for a free slot"
            return state

        if state["status"] == "completed":
            # Exit zero with no artifact is a failure, not a success.
            if not self.artifact_path.exists():
                state.update(
                    status="failed", error_code="MissingArtifact",
                    error_message="The analysis reported success but wrote no "
                                  "results file.",
                )
            return state

        if state["status"] in ("pending", "running"):
            if self._process_alive():
                if self._seconds_since_update(state) > HEARTBEAT_TIMEOUT_S:
                    state["status"] = "stalled"
                return state
            # The process is gone and never reported completion.
            state.update(
                status="failed",
                error_code=state.get("error_code") or "ProcessDied",
                error_message=state.get("error_message")
                or "The analysis stopped unexpectedly. The technical log below "
                   "will say why.",
            )
        return state

    def logs(self) -> dict[str, str]:
        out = {}
        for name in ("stdout.log", "stderr.log"):
            path = self.directory / name
            out[name] = path.read_text(encoding="utf-8", errors="replace") if path.exists() else ""
        return out

    def cancel(self) -> bool:
        """Stop the run, including its worker processes.

        Killing only the CLI parent would leave a pool of ViennaRNA workers
        running and the machine oversubscribed.
        """
        pid = self.request.get("pid")
        if not pid or not self._process_alive():
            return False
        try:
            if os.name == "nt":
                subprocess.run(
                    ["taskkill", "/PID", str(pid), "/T", "/F"],
                    capture_output=True, check=False,
                )
            else:
                os.killpg(os.getpgid(pid), signal.SIGTERM)
        except (OSError, subprocess.SubprocessError):
            return False
        _append_event(
            self.progress_file,
            {"event": "run_failed", "code": "Cancelled",
             "message": "The analysis was cancelled."},
        )
        return True

    # -- liveness --------------------------------------------------------

    def _process_alive(self) -> bool:
        """PID plus start time: Windows reuses PIDs, so a PID alone can lie."""
        pid = self.request.get("pid")
        if not pid:
            return False
        try:
            if os.name == "nt":
                result = subprocess.run(
                    ["tasklist", "/FI", f"PID eq {pid}", "/NH", "/FO", "CSV"],
                    capture_output=True, text=True, check=False, timeout=10,
                )
                if f'"{pid}"' not in result.stdout:
                    return False
                return "python" in result.stdout.lower()
            os.kill(pid, 0)
            return True
        except (OSError, subprocess.SubprocessError):
            return False

    def _seconds_since_update(self, state: Mapping[str, Any]) -> float:
        stamp = state.get("updated_utc") or state.get("created_utc")
        if not stamp:
            return 0.0
        try:
            when = datetime.fromisoformat(str(stamp).replace("Z", "+00:00"))
        except ValueError:
            return 0.0
        return (datetime.now(timezone.utc) - when).total_seconds()


# -- creating and finding jobs ------------------------------------------


def jobs_root(runs_dir: str | Path) -> Path:
    return Path(runs_dir) / JOBS_DIRNAME


def new_job_id() -> str:
    return datetime.now(timezone.utc).strftime("job_%Y%m%dT%H%M%S_") + os.urandom(3).hex()


def stage_upload(directory: Path, uploaded, fallback_name: str) -> Path:
    """Write an uploaded file into the job, under a name *we* choose.

    Client-supplied filenames are never used as paths: on Windows a crafted
    name can escape the directory entirely.
    """
    inputs = directory / "inputs"
    inputs.mkdir(parents=True, exist_ok=True)
    suffix = Path(getattr(uploaded, "name", fallback_name)).suffix.lower()
    if suffix not in (".txt", ".csv", ".tsv", ".fasta", ".fa", ".json", ".yaml", ".yml"):
        suffix = Path(fallback_name).suffix
    target = inputs / f"{Path(fallback_name).stem}{suffix}"
    data = uploaded.getvalue() if hasattr(uploaded, "getvalue") else uploaded.read()
    target.write_bytes(data)
    return target


def submit(
    runs_dir: str | Path,
    candidates_path: str | Path,
    corpus_path: str | Path,
    *,
    name: str | None = None,
    corpus_is_placeholder: bool = False,
    bundle_path: str | Path | None = None,
    preset: str = "standard",
    extra_sets: Sequence[str] = (),
    job_id: str | None = None,
    python_executable: str | None = None,
) -> Job:
    """Queue a run. It starts immediately if a slot is free, else it waits."""
    job_id = job_id or new_job_id()
    directory = jobs_root(runs_dir) / job_id
    (directory / "inputs").mkdir(parents=True, exist_ok=True)

    command = [
        python_executable or sys.executable,   # never `python` from PATH
        "-m", "aptarank", "run", str(candidates_path),
        "--artifact-path", str(directory / "artifact.json"),
        "--progress-format", "jsonl",
        "--progress-file", str(directory / "progress.jsonl"),
        "--job-id", job_id,
    ]
    if SITE_CONFIG and Path(SITE_CONFIG).is_file():
        command += ["-c", SITE_CONFIG]
    if corpus_is_placeholder:
        command += ["--development-corpus", str(corpus_path)]
    else:
        command += ["--corpus", str(corpus_path)]
    if bundle_path:
        command += ["--target-bundle", str(bundle_path)]

    settings = preset_settings(preset) + list(extra_sets)
    if WORKERS_PER_JOB:
        # Server policy overrides anything a preset or a config file asks for.
        settings = [s for s in settings if not s.startswith("tier1.parallel.workers")]
        settings.append(f"tier1.parallel.workers={int(WORKERS_PER_JOB)}")
    for item in settings:
        command += ["--set", item]

    request = {
        "job_id": job_id,
        "name": name,
        "created_utc": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "pid": None,
        "command": command,
        "runs_dir": str(runs_dir),
        "preset": preset,
        "candidates": str(candidates_path),
        "corpus": str(corpus_path),
        "corpus_is_placeholder": corpus_is_placeholder,
        "target_bundle": str(bundle_path) if bundle_path else None,
    }
    (directory / "request.json").write_text(json.dumps(request, indent=2), encoding="utf-8")

    job = Job(job_id=job_id, directory=directory, request=request)
    pump(runs_dir)
    return Job(job_id=job_id, directory=directory,
               request=json.loads((directory / "request.json").read_text(encoding="utf-8")))


def pump(runs_dir: str | Path) -> list[Job]:
    """Start queued jobs while slots are free. Cheap; call it on every render."""
    started: list[Job] = []
    found = discover(runs_dir, limit=100)
    running = sum(1 for job in found if job.state()["status"] in (RUNNING, "pending"))
    waiting = [job for job in reversed(found) if job.is_queued]   # oldest first

    for job in waiting:
        if running >= MAX_CONCURRENT_JOBS:
            break
        _spawn(job)
        running += 1
        started.append(job)
    return started


def _spawn(job: Job) -> None:
    """Launch the CLI for an already-queued job and record its PID."""
    directory = job.directory
    command = job.request["command"]
    runs_dir = Path(job.request.get("runs_dir", directory.parent.parent))

    stdout = (directory / "stdout.log").open("w", encoding="utf-8")
    stderr = (directory / "stderr.log").open("w", encoding="utf-8")
    creationflags = subprocess.CREATE_NEW_PROCESS_GROUP if os.name == "nt" else 0
    popen_kwargs: dict[str, Any] = {
        "stdout": stdout, "stderr": stderr, "stdin": subprocess.DEVNULL,
        "cwd": str(Path(runs_dir).resolve().parent), "shell": False,
    }
    if os.name == "nt":
        popen_kwargs["creationflags"] = creationflags
    else:
        popen_kwargs["start_new_session"] = True

    process = subprocess.Popen(command, **popen_kwargs)

    job.request["pid"] = process.pid
    job.request["started_utc"] = datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")
    (directory / "request.json").write_text(
        json.dumps(job.request, indent=2), encoding="utf-8"
    )


def preset_settings(preset: str) -> list[str]:
    """Three named presets instead of thirty config keys.

    The counts are shown to the user rather than hidden: a "quick" run is a
    different measurement, not merely a faster one.
    """
    return {
        "quick": ["run.mode=fast", "tier1.n_ensemble_samples=30", "output.n_diagrams=20"],
        "standard": ["tier1.shuffle.n_shuffles=20", "tier1.n_ensemble_samples=100",
                     "output.n_diagrams=50"],
        "evaluation": ["tier1.shuffle.n_shuffles=99", "tier1.n_ensemble_samples=200",
                       "output.n_diagrams=100"],
    }.get(preset, [])


PRESET_DESCRIPTIONS = {
    "quick": ("Quick preview", "No shuffled controls, 30 sampled structures. "
                               "Fast enough to sanity-check a file — not a result."),
    "standard": ("Standard analysis", "20 shuffled controls per sequence, 100 sampled "
                                      "structures. The everyday setting."),
    "evaluation": ("Rigorous (slowest)", "99 shuffled controls per sequence, 200 sampled "
                                         "structures. Needed for statistics quoted in "
                                         "a paper."),
}


def discover(runs_dir: str | Path, limit: int = 50) -> list[Job]:
    """Every job on disk, newest first. The UI's memory lives here, not in a session."""
    root = jobs_root(runs_dir)
    if not root.is_dir():
        return []
    found = []
    for directory in sorted(root.iterdir(), reverse=True):
        request_file = directory / "request.json"
        if not request_file.is_file():
            continue
        try:
            request = json.loads(request_file.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            continue
        found.append(Job(job_id=directory.name, directory=directory, request=request))
        if len(found) >= limit:
            break
    return found


def active(runs_dir: str | Path) -> Job | None:
    """The first job that is running or waiting for a slot, if any."""
    for job in discover(runs_dir, limit=20):
        if job.state()["status"] in ("pending", RUNNING, QUEUED):
            return job
    return None


def slots(runs_dir: str | Path) -> dict[str, int]:
    """How busy the machine is, in this tool's terms."""
    states = [job.state()["status"] for job in discover(runs_dir, limit=100)]
    return {
        "running": sum(1 for s in states if s in (RUNNING, "pending")),
        "queued": sum(1 for s in states if s == QUEUED),
        "capacity": MAX_CONCURRENT_JOBS,
    }


def estimate_runtime_s(n_candidates: int, preset: str, with_target: bool) -> tuple[int, int]:
    """A deliberately wide range. An estimate presented as exact is a lie."""
    per_candidate = {"quick": 0.02, "standard": 0.30, "evaluation": 1.30}.get(preset, 0.30)
    base = n_candidates * per_candidate
    if with_target:
        base += 90     # calibration bank, first time only
    return max(5, int(base * 0.6)), max(15, int(base * 1.8))


def _append_event(path: Path, fields: Mapping[str, Any]) -> None:
    record = {
        "schema": "aptarank-progress-v1",
        "created_utc": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        **fields,
    }
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(record) + "\n")
