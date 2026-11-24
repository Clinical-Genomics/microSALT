"""Persistent record of every command run and every check made during the validation."""

import json
from dataclasses import asdict, dataclass, field
from datetime import datetime
from pathlib import Path


def now() -> str:
    return datetime.now().isoformat(timespec="seconds")


@dataclass
class RunRecord:
    id: str
    phase: str
    test: str
    case: str
    command: str
    started: str
    duration_s: float
    exit_code: int
    log: str
    finishdir: str | None = None
    slurm_ids: list[str] = field(default_factory=list)
    note: str = ""

    @property
    def ok(self) -> bool:
        return self.exit_code == 0


class State:
    def __init__(self, path: Path):
        self.path = path
        self.data = {
            "created": now(),
            "runs": [],
            "checks": {},
            "snapshots": {},
            "manual": {},
        }
        if path.exists():
            with open(path) as stream:
                self.data = json.load(stream)

    def save(self) -> None:
        with open(self.path, "w") as stream:
            json.dump(self.data, stream, indent=2, sort_keys=True)

    def add_run(self, record: RunRecord) -> None:
        self.data["runs"] = [run for run in self.data["runs"] if run["id"] != record.id]
        self.data["runs"].append(asdict(record))
        self.save()

    def run(self, run_id: str) -> RunRecord | None:
        for run in self.data["runs"]:
            if run["id"] == run_id:
                return RunRecord(**run)
        return None

    def runs(
        self, phase: str | None = None, test: str | None = None
    ) -> list[RunRecord]:
        return [
            RunRecord(**run)
            for run in self.data["runs"]
            if (phase is None or run["phase"] == phase)
            and (test is None or run["test"] == test)
        ]

    def check(self, check_id: str, passed: bool | None, detail: str = "") -> None:
        self.data["checks"][check_id] = {
            "passed": passed,
            "detail": detail,
            "when": now(),
        }
        self.save()
        status = {True: "PASS", False: "FAIL", None: "N/A"}[passed]
        print(f"  [{status}] {check_id}: {detail}")

    def get_check(self, check_id: str) -> dict | None:
        return self.data["checks"].get(check_id)

    def checks(self, prefix: str) -> dict[str, dict]:
        return {k: v for k, v in self.data["checks"].items() if k.startswith(prefix)}

    def snapshot(self, name: str, data: dict) -> None:
        self.data["snapshots"][name] = data
        self.save()

    def get_snapshot(self, name: str) -> dict | None:
        return self.data["snapshots"].get(name)

    def manual(self, step_id: str, instruction: str) -> None:
        self.data["manual"][step_id] = {"instruction": instruction, "when": now()}
        self.save()
        print(f"  [MANUAL] {step_id}: {instruction}")
