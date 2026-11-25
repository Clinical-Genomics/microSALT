"""Run the microSALT CLI and record what happened."""

import shlex
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path

import yaml

from .config import ValidationConfig
from .state import RunRecord, State, now


@dataclass
class RunResult:
    command: str
    exit_code: int
    stdout: str
    stderr: str
    duration_s: float
    log: Path

    @property
    def output(self) -> str:
        return self.stdout + self.stderr


class MicrosaltRunner:
    def __init__(self, cfg: ValidationConfig, state: State):
        self.cfg = cfg
        self.state = state

    def _log_path(self, run_id: str) -> Path:
        return self.cfg.log_dir / (run_id.replace("/", "__") + ".log")

    def run(
        self, run_id: str, args: list[str], config_path: Path | None = None
    ) -> RunResult:
        command = [
            str(self.cfg.microsalt_bin),
            "--config",
            str(config_path or self.cfg.microsalt_config),
            *args,
        ]
        printable = shlex.join(command)
        print(f"$ {printable}")
        start = time.monotonic()
        proc = subprocess.run(command, capture_output=True, text=True)
        duration = time.monotonic() - start
        log = self._log_path(run_id)
        with open(log, "w") as stream:
            stream.write(
                f"# {printable}\n# exit {proc.returncode} after {duration:.1f}s\n\n"
            )
            stream.write("## stdout\n")
            stream.write(proc.stdout)
            stream.write("\n## stderr\n")
            stream.write(proc.stderr)
        print(f"  exit {proc.returncode} after {duration:.1f}s (log: {log})")
        return RunResult(
            printable, proc.returncode, proc.stdout, proc.stderr, duration, log
        )

    def _new_finishdir(self, before: set[str], job_name: str) -> Path | None:
        after = {p.name for p in self.cfg.results_root.iterdir()}
        candidates = sorted(
            (
                self.cfg.results_root / name
                for name in after - before
                if name.startswith(f"{job_name}_")
            ),
            key=lambda p: p.stat().st_mtime,
        )
        return candidates[-1] if candidates else None

    @staticmethod
    def _slurm_ids(finishdir: Path | None, job_name: str) -> list[str]:
        if finishdir is None:
            return []
        path = finishdir / f"{job_name}_slurm_ids.yaml"
        if not path.exists():
            return []
        with open(path) as stream:
            return [str(job) for job in yaml.safe_load(stream).get("jobs", [])]

    def analyse(
        self,
        run_id: str,
        phase: str,
        test: str,
        case: str,
        sampleinfo: Path,
        fastq: Path,
        job_name: str,
        *,
        qc_only: bool = False,
        dry: bool = False,
        skip_update: bool = False,
        config_path: Path | None = None,
        note: str = "",
    ) -> tuple[RunRecord, RunResult]:
        args = ["analyse", "--input", str(fastq)]
        if qc_only:
            args.append("--qc-only")
        if dry:
            args.append("--dry")
        if skip_update:
            args.append("--skip_update")
        args.append(str(sampleinfo))

        self.cfg.results_root.mkdir(parents=True, exist_ok=True)
        before = {p.name for p in self.cfg.results_root.iterdir()}
        started = now()
        result = self.run(run_id, args, config_path)
        finishdir = self._new_finishdir(before, job_name)
        record = RunRecord(
            id=run_id,
            phase=phase,
            test=test,
            case=case,
            command=result.command,
            started=started,
            duration_s=round(result.duration_s, 2),
            exit_code=result.exit_code,
            log=str(result.log),
            finishdir=str(finishdir) if finishdir else None,
            slurm_ids=self._slurm_ids(finishdir, job_name),
            note=note,
        )
        self.state.add_run(record)
        return record, result

    def finish(
        self,
        run_id: str,
        phase: str,
        test: str,
        case: str,
        sampleinfo: Path,
        finishdir: Path,
        *,
        qc_only: bool = False,
        report: str = "default",
        note: str = "",
    ) -> tuple[RunRecord, RunResult]:
        args = ["utils", "finish", "--input", str(finishdir), "--report", report]
        if qc_only:
            args.append("--qc-only")
        args.append(str(sampleinfo))
        started = now()
        result = self.run(run_id, args)
        record = RunRecord(
            id=run_id,
            phase=phase,
            test=test,
            case=case,
            command=result.command,
            started=started,
            duration_s=round(result.duration_s, 2),
            exit_code=result.exit_code,
            log=str(result.log),
            finishdir=str(finishdir),
            note=note,
        )
        self.state.add_run(record)
        return record, result


def setup_errors(output: str) -> list[str]:
    """Return the bullet list printed by AnalysisSetupError, if present."""
    if "Analysis setup failed" not in output:
        return []
    return [line[2:].strip() for line in output.splitlines() if line.startswith("- ")]
