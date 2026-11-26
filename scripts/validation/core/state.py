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

    @classmethod
    def from_dict(cls, data: dict[str, object]) -> "RunRecord":
        slurm_ids = data.get("slurm_ids", [])
        duration = data["duration_s"]
        exit_code = data["exit_code"]
        if not isinstance(duration, (int, float)):
            duration = float(str(duration))
        if not isinstance(exit_code, int):
            exit_code = int(str(exit_code))
        return cls(
            id=str(data["id"]),
            phase=str(data["phase"]),
            test=str(data["test"]),
            case=str(data["case"]),
            command=str(data["command"]),
            started=str(data["started"]),
            duration_s=float(duration),
            exit_code=exit_code,
            log=str(data["log"]),
            finishdir=(
                str(data["finishdir"]) if data.get("finishdir") is not None else None
            ),
            slurm_ids=(
                [str(item) for item in slurm_ids] if isinstance(slurm_ids, list) else []
            ),
            note=str(data.get("note", "")),
        )


class State:
    def __init__(self, path: Path):
        self.path = path
        self.data: dict[str, object] = {
            "created": now(),
            "runs": [],
            "checks": {},
            "snapshots": {},
            "manual": {},
        }
        if path.exists():
            with open(path) as stream:
                loaded = json.load(stream)
            if isinstance(loaded, dict):
                self.data = loaded

    def _runs(self) -> list[dict[str, object]]:
        runs = self.data.get("runs", [])
        if not isinstance(runs, list):
            return []
        typed: list[dict[str, object]] = []
        for run in runs:
            if not isinstance(run, dict):
                continue
            typed.append({str(key): value for key, value in run.items()})
        return typed

    def _checks(self) -> dict[str, dict[str, object]]:
        checks = self.data.get("checks", {})
        if not isinstance(checks, dict):
            return {}
        typed: dict[str, dict[str, object]] = {}
        for key, value in checks.items():
            if isinstance(key, str) and isinstance(value, dict):
                typed[key] = {
                    str(inner_key): inner_value
                    for inner_key, inner_value in value.items()
                }
        return typed

    def _snapshots(self) -> dict[str, dict]:
        snapshots = self.data.get("snapshots", {})
        if not isinstance(snapshots, dict):
            return {}
        return {
            key: value
            for key, value in snapshots.items()
            if isinstance(key, str) and isinstance(value, dict)
        }

    def _manual(self) -> dict[str, dict]:
        manual = self.data.get("manual", {})
        if not isinstance(manual, dict):
            return {}
        return {
            key: value
            for key, value in manual.items()
            if isinstance(key, str) and isinstance(value, dict)
        }

    def save(self) -> None:
        with open(self.path, "w") as stream:
            json.dump(self.data, stream, indent=2, sort_keys=True)

    def add_run(self, record: RunRecord) -> None:
        previous = self.run(record.id)
        if previous is not None:
            history = self.data.setdefault("run_history", [])
            if not isinstance(history, list):
                raise ValueError("Invalid run_history in validation state")
            history.append(asdict(previous))
        runs = [run for run in self._runs() if run.get("id") != record.id]
        runs.append(asdict(record))
        self.data["runs"] = runs
        self.save()

    def run(self, run_id: str) -> RunRecord | None:
        for run in self._runs():
            if run.get("id") == run_id:
                return RunRecord.from_dict(run)
        return None

    def runs(
        self, phase: str | None = None, test: str | None = None
    ) -> list[RunRecord]:
        return [
            RunRecord.from_dict(run)
            for run in self._runs()
            if (phase is None or run.get("phase") == phase)
            and (test is None or run.get("test") == test)
        ]

    def check(self, check_id: str, passed: bool | None, detail: str = "") -> None:
        checks = self._checks()
        if check_id in checks:
            history = self.data.setdefault("check_history", [])
            if not isinstance(history, list):
                raise ValueError("Invalid check_history in validation state")
            history.append({"id": check_id, **checks[check_id]})
        checks[check_id] = {
            "passed": passed,
            "detail": detail,
            "when": now(),
        }
        self.data["checks"] = checks
        self.save()
        status = {True: "PASS", False: "FAIL", None: "N/A"}[passed]
        print(f"  [{status}] {check_id}: {detail}")

    def get_check(self, check_id: str) -> dict[str, object] | None:
        return self._checks().get(check_id)

    def checks(self, prefix: str) -> dict[str, dict[str, object]]:
        return {k: v for k, v in self._checks().items() if k.startswith(prefix)}

    def invalidate_checks(self, prefix: str, reason: str) -> None:
        """Keep old checks in history, but mark them pending after a retry."""
        for key in list(self.checks(prefix)):
            self.check(key, None, reason)

    def snapshot(self, name: str, data: dict) -> None:
        snapshots = self._snapshots()
        snapshots[name] = data
        self.data["snapshots"] = snapshots
        self.save()

    def get_snapshot(self, name: str) -> dict | None:
        return self._snapshots().get(name)

    def manual(self, step_id: str, instruction: str) -> None:
        manual = self._manual()
        manual[step_id] = {"instruction": instruction, "when": now()}
        self.data["manual"] = manual
        self.save()
        print(f"  [MANUAL] {step_id}: {instruction}")
