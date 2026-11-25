"""SLURM and finishdir completion lookups."""

import shutil
import subprocess
from pathlib import Path

from .state import RunRecord

TERMINAL = {"COMPLETED", "FAILED", "CANCELLED", "TIMEOUT", "OUT_OF_MEMORY", "NODE_FAIL"}


def job_states(job_ids: list[str]) -> dict[str, str]:
    if not job_ids or shutil.which("sacct") is None:
        return {}
    proc = subprocess.run(
        ["sacct", "-X", "-n", "-P", "-o", "JobID,State", "-j", ",".join(job_ids)],
        capture_output=True,
        text=True,
    )
    states = {}
    for line in proc.stdout.splitlines():
        job, _, state = line.partition("|")
        states[job.strip()] = state.strip().split()[0] if state.strip() else "UNKNOWN"
    return states


def run_complete(record: RunRecord) -> bool:
    return (
        bool(record.finishdir)
        and (Path(record.finishdir) / "run_complete.out").exists()
    )


def run_status(record: RunRecord) -> str:
    """One of: submission-failed, complete, failed, running, pending, unknown."""
    if record.exit_code != 0:
        return "submission-failed"
    if run_complete(record):
        return "complete"
    if not record.slurm_ids:
        return "unknown"
    states = job_states(record.slurm_ids)
    if not states:
        return "unknown"
    values = set(states.values())
    if values <= TERMINAL:
        return "failed"
    if any(value in {"RUNNING", "COMPLETING"} for value in values):
        return "running"
    return "pending"
