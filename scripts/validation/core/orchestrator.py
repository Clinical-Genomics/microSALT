"""Instantiate the configured suites and drive them through the phases."""

import shutil
import sys
from pathlib import Path

from .config import ValidationConfig
from .fixtures import FixtureHelper
from .report import render_report
from .runner import MicrosaltRunner
from .slurm import TERMINAL, job_states, run_complete, run_status
from .state import State
from .suite import PHASES, Context, Suite


class Orchestrator:
    def __init__(
        self,
        cfg: ValidationConfig,
        state: State,
        registry: dict[str, type[Suite]],
        broken: dict[str, str] | None = None,
    ):
        self.cfg = cfg
        self.state = state
        self.ctx = Context(cfg, state, MicrosaltRunner(cfg, state), FixtureHelper(cfg))
        self.suites: list[Suite] = []
        for name, options in cfg.suites.items():
            if name not in registry:
                hint = ""
                if broken:
                    hint = "\nSuite modules that failed to import:\n" + "\n".join(
                        f"  {module}: {error}"
                        for module, error in sorted(broken.items())
                    )
                raise SystemExit(
                    f"unknown suite {name!r}; available: {', '.join(sorted(registry))}{hint}"
                )
            missing = [dep for dep in registry[name].requires if dep not in cfg.suites]
            if missing:
                raise SystemExit(f"suite {name!r} requires {', '.join(missing)}")
            self.suites.append(registry[name](self.ctx, options or {}))

    def _require_phases_complete(self, phases: list[str], force: bool) -> None:
        pending = [
            f"{record.id}: {run_status(record)}"
            for phase in phases
            for record in self.state.runs(phase=phase)
            if record.slurm_ids and not run_complete(record)
        ]
        if pending and not force:
            print("Cannot continue, submitted runs are not complete:")
            for line in pending:
                print(f"  {line}")
            print("Re-run with --force to continue anyway.")
            sys.exit(1)

    def start(self, phase: str, force: bool = False) -> None:
        earlier = list(PHASES[: PHASES.index(phase)])
        self._require_phases_complete(earlier, force)
        for suite in self.suites:
            getattr(suite, phase)()
        if phase == "submit":
            for step, instruction in self.cfg.manual_steps.items():
                self.state.manual(step, instruction)

    def status(self) -> None:
        rows = [
            (record.id, run_status(record), record.finishdir or "-")
            for record in self.state.runs()
            if record.slurm_ids
            or (record.phase != "finish" and record.ok and record.finishdir)
        ]
        if not rows:
            print("No submitted runs recorded.")
            return
        width = max(len(r[0]) for r in rows)
        for run_id, status, finishdir in rows:
            print(f"{run_id:<{width}}  {status:<18}  {finishdir}")

    def retry_options(self):
        return [
            (suite, action) for suite in self.suites for action in suite.retry_actions()
        ]

    def retry(self, action_id: str) -> None:
        options = {action.id: (suite, action) for suite, action in self.retry_options()}
        if action_id not in options:
            raise ValueError(
                f"{action_id!r} is not a retryable test in the configured suites"
            )
        suite, action = options[action_id]

        # Never submit a duplicate while an earlier job may still be running. Include IDs from
        # failed/partially submitted commands; do not trust only run_complete.out or exit codes.
        records = [record for record in self.state.runs() if record.slurm_ids]
        all_ids = list(
            dict.fromkeys(job for record in records for job in record.slurm_ids)
        )
        states = job_states(all_ids)
        unresolved = [
            job
            for job in all_ids
            if states.get(job, "UNKNOWN").rstrip("+") not in TERMINAL
        ]
        if unresolved:
            raise ValueError(
                "Cannot retry while recorded SLURM jobs are running or their status is unknown: "
                + ", ".join(f"{job}={states.get(job, 'UNKNOWN')}" for job in unresolved)
                + ". Wait for them to finish and ensure sacct is available."
            )
        suite.retry(action.id)

    def collect(self, copy_results: bool, output: Path) -> None:
        for suite in self.suites:
            suite.collect()
        if copy_results:
            self._copy_results()
        output.write_text(render_report(self.cfg, self.state, self.suites))
        print(f"wrote {output}")

    def _copy_results(self) -> None:
        target_root = self.cfg.validation_root / "results"
        target_root.mkdir(exist_ok=True)
        for record in self.state.runs():
            if record.phase == "finish" or not record.finishdir:
                continue
            source = Path(record.finishdir)
            target = target_root / source.name
            if not source.exists() or target.exists():
                continue
            print(f"copying {source} -> {target}")
            shutil.copytree(source, target, symlinks=True)
