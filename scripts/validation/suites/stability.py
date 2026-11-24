"""Core suite, part of every release validation: the standard cases run to completion and the
typing results are concordant with the previous version."""

import subprocess
import sys
from pathlib import Path

from core.report import detail, table, yn
from core.slurm import run_complete, run_status
from core.suite import Suite, SuiteOptions

COMPARE_SCRIPT = Path(__file__).resolve().parents[2] / "compare_typing_reports.py"


class StabilityOptions(SuiteOptions):
    cases: list[str] | None = None  # default: every configured case


class StabilitySuite(Suite):
    name = "stability"
    title = "Analysis stability"
    Options = StabilityOptions

    @property
    def cases(self):
        names = self.options.cases or self.cfg.case_names
        return [self.cfg.case(name) for name in names]

    def submit(self) -> None:
        print("== Full analyses ==")
        for case in self.cases:
            self.analyse(
                "submit",
                case.name,
                case.name,
                case.sampleinfo,
                case.fastq,
                case.project,
            )

    def collect(self) -> None:
        compare_dir = self.cfg.validation_root / "compare"
        compare_dir.mkdir(exist_ok=True)
        for case in self.cases:
            record = self.run("submit", case.name)
            if record is None:
                continue
            self.check(
                f"{case.name}/analysis_complete",
                run_complete(record),
                f"{run_status(record)} ({record.finishdir})",
            )
            self._compare_typing(
                case.name, case.old_typing_report, record.finishdir, compare_dir
            )

    def _compare_typing(
        self, case: str, old: Path | None, finishdir: str | None, compare_dir: Path
    ) -> None:
        key = f"{case}/typing_concordant"
        if old is None:
            self.check(key, None, "no old typing report configured")
            return
        reports = (
            sorted(
                Path(finishdir).glob("*_Typing_*.html"), key=lambda p: p.stat().st_mtime
            )
            if finishdir
            else []
        )
        if not reports:
            self.check(key, False, "no new typing report found in finishdir")
            return
        new = reports[-1]
        proc = subprocess.run(
            [sys.executable, str(COMPARE_SCRIPT), "--no-color", str(old), str(new)],
            capture_output=True,
            text=True,
        )
        (compare_dir / f"{case}.txt").write_text(proc.stdout + proc.stderr)
        self.check(
            key,
            proc.returncode == 0,
            f"compare_typing_reports exit {proc.returncode} ({old.name} vs {new.name})",
        )

    def render(self) -> list[str]:
        lines = ["#### Analysis Command", ""]
        lines += table(
            ["CaseID", "Analysis finishes without errors (Y/N)"],
            [
                [case.name, yn(self.get_check(f"{case.name}/analysis_complete"))]
                for case in self.cases
            ],
        )
        lines += ["", "#### Verify generated typing reports", ""]
        compare_dir = self.cfg.validation_root / "compare"
        for case in self.cases:
            check = self.get_check(f"{case.name}/typing_concordant")
            output = compare_dir / f"{case.name}.txt"
            lines += [
                f"##### Typing report comparison - {case.name}",
                "",
                "```plaintext",
            ]
            lines.append(
                output.read_text().rstrip() if output.exists() else detail(check)
            )
            lines += ["```", ""]
            mark = "x" if check and check["passed"] else " "
            lines += [
                f"- [{mark}] Case passed the criteria of identical typing reports.",
                "",
            ]
        return lines
