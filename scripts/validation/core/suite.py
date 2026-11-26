"""Base class for validation suites and the shared context they run in."""

from dataclasses import dataclass, field
from pathlib import Path
from typing import ClassVar, TypeVar

from pydantic import BaseModel, ConfigDict

from .config import ValidationConfig
from .db import Database
from .fixtures import FixtureHelper
from .runner import MicrosaltRunner
from .state import RunRecord, State

PHASES = ("submit", "resubmit", "finish")
"""Phases in order. `submit` and `resubmit` start SLURM jobs and are followed by a wait; `finish`
runs on the login node. `collect` runs after everything and gathers results."""


@dataclass
class Context:
    cfg: ValidationConfig
    state: State
    runner: MicrosaltRunner
    fixtures: FixtureHelper
    _db: Database | None = field(default=None, repr=False)

    @property
    def db(self) -> Database:
        if self._db is None:
            self._db = Database(self.cfg.microsalt.database.SQLALCHEMY_DATABASE_URI)
        return self._db


class SuiteOptions(BaseModel):
    model_config = ConfigDict(extra="forbid")


@dataclass(frozen=True)
class RetryAction:
    """One independently retryable action exposed by an enabled suite."""

    id: str
    description: str
    phase: str
    run_id: str | None = None


OptionsT = TypeVar("OptionsT", bound=SuiteOptions)


class Suite:
    """One group of criteria. Subclasses override the phase methods they take part in.

    Check ids and run ids are namespaced by `name`, so a suite's results can be found by prefix.
    """

    name: ClassVar[str]
    title: ClassVar[str]
    requires: ClassVar[tuple[str, ...]] = ()
    Options: ClassVar[type[SuiteOptions]] = SuiteOptions

    def __init__(self, ctx: Context, options: dict):
        self.ctx = ctx
        self.cfg = ctx.cfg
        self.state = ctx.state
        self.runner = ctx.runner
        self.fixtures = ctx.fixtures
        parsed = self.Options(**options)
        self.options: SuiteOptions = parsed

    # -- phases --------------------------------------------------------------------------------

    def submit(self) -> None:
        pass

    def resubmit(self) -> None:
        pass

    def finish(self) -> None:
        pass

    def collect(self) -> None:
        pass

    def render(self) -> list[str]:
        """Markdown lines for this suite's Results subsection."""
        return []

    def retry_actions(self) -> list[RetryAction]:
        """Actions available for the configured cases in this suite."""
        return []

    def retry(self, action_id: str) -> None:
        raise ValueError(f"Suite {self.name} does not support retrying {action_id}")

    # -- helpers -------------------------------------------------------------------------------

    def check(self, key: str, passed: bool | None, detail: str = "") -> None:
        self.state.check(f"{self.name}/{key}", passed, detail)

    def get_check(self, key: str) -> dict | None:
        return self.state.get_check(f"{self.name}/{key}")

    def snapshot(self, key: str, data: dict) -> None:
        self.state.snapshot(f"{self.name}/{key}", data)

    def get_snapshot(self, key: str) -> dict | None:
        return self.state.get_snapshot(f"{self.name}/{key}")

    def run_id(self, phase: str, label: str) -> str:
        return f"{phase}/{self.name}/{label}"

    def run(self, phase: str, label: str) -> RunRecord | None:
        return self.state.run(self.run_id(phase, label))

    def analyse(
        self,
        phase: str,
        label: str,
        case: str,
        sampleinfo: Path,
        fastq: Path,
        job_name: str,
        **kwargs,
    ):
        return self.runner.analyse(
            self.run_id(phase, label),
            phase,
            self.name,
            case,
            sampleinfo,
            fastq,
            job_name,
            **kwargs,
        )

    def finish_run(
        self, label: str, case: str, sampleinfo: Path, finishdir: Path, **kwargs
    ):
        return self.runner.finish(
            self.run_id("finish", label),
            "finish",
            self.name,
            case,
            sampleinfo,
            finishdir,
            **kwargs,
        )

    def result(self) -> str:
        """Do not report PASS while an existing check is pending or unavailable."""
        checks = list(self.state.checks(f"{self.name}/").values())
        if not checks:
            return "-"
        if any(c["passed"] is False for c in checks):
            return "FAIL"
        if any(c["passed"] is None for c in checks):
            return "INCOMPLETE"
        return "PASS"
