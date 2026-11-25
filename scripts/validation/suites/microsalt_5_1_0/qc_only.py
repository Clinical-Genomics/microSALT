"""QC-only workflow: `analyse --qc-only` writes QC-only job scripts, `utils finish --qc-only`
produces QC-only outputs, and the QC metrics match the full analysis."""

import re
from pathlib import Path

import yaml

from microSALT.utils.reporter import TYPING_ASSEMBLY_DELIVERABLE_TAGS

from core.db import diff_qc
from core.report import detail, table, yn
from core.state import RunRecord
from core.suite import Suite, SuiteOptions

ASSEMBLY_OR_TYPING = re.compile(r"\b(skesa|quast\.py|blastn)\b")


def runfiles(finishdir: str | None) -> list[Path]:
    if not finishdir:
        return []
    root = Path(finishdir)
    return sorted([*root.glob("runfile.sbatch"), *root.glob("*/runfile.sbatch")])


class QcOnlyOptions(SuiteOptions):
    cases: list[str] | None = None


class QcOnlySuite(Suite):
    name = "qc_only"
    title = "QC-only analysis"
    requires = ("stability",)
    Options = QcOnlyOptions
    options: QcOnlyOptions

    @property
    def cases(self):
        return [
            self.cfg.case(name) for name in (self.options.cases or self.cfg.case_names)
        ]

    def resubmit(self) -> None:
        print("== Snapshots after full analysis ==")
        for case in self.cases:
            self.snapshot(
                f"{case.name}/after_full", self.ctx.db.snapshot_project(case.project)
            )
        print("== QC-only analyses ==")
        for case in self.cases:
            record, _ = self.analyse(
                "resubmit",
                case.name,
                case.name,
                case.sampleinfo,
                case.fastq,
                case.project,
                qc_only=True,
            )
            self._check_runfiles(case.name, record)

    def _check_runfiles(self, case: str, record: RunRecord) -> None:
        files = runfiles(record.finishdir)
        bad = [p for p in files if ASSEMBLY_OR_TYPING.search(p.read_text())]
        self.check(
            f"{case}/runfile_qc_only",
            record.ok and bool(files) and not bad,
            f"{len(files)} runfile(s), {len(bad)} with skesa/quast/blastn sections",
        )
        mailjob = Path(record.finishdir) / "mailjob.sh" if record.finishdir else None
        content = mailjob.read_text() if mailjob and mailjob.exists() else ""
        self.check(
            f"{case}/mailjob_qc_track",
            "--report qc" in content and "--track qc" in content,
            "mailjob.sh finishes with --report qc --track qc",
        )

    def finish(self) -> None:
        for case in self.cases:
            print(f"== QC-only finish: {case.name} ==")
            qc_run = self.run("resubmit", case.name)
            if qc_run is None or not qc_run.finishdir:
                self.check(f"{case.name}/finish_ok", None, "Missing QC-only resubmit run or finish directory")
                continue
            finishdir = Path(qc_run.finishdir)
            record, _ = self.finish_run(
                case.name, case.name, case.sampleinfo, finishdir, qc_only=True
            )
            self.check(f"{case.name}/finish_ok", record.ok, f"exit {record.exit_code}")

            typing_html = list(finishdir.glob("*_Typing_*.html"))
            qc_html = list(finishdir.glob("*_QC_*.html"))
            leftovers = [
                p
                for pattern in ("*/assembly", "*/blast_search", "*/quast")
                for p in finishdir.glob(pattern)
            ]
            typing_tags = set()
            for deliverable in finishdir.glob("*_deliverables.yaml"):
                with open(deliverable) as stream:
                    tags = {
                        f.get("tag") for f in yaml.safe_load(stream).get("files", [])
                    }
                typing_tags |= tags & TYPING_ASSEMBLY_DELIVERABLE_TAGS
            self.check(
                f"{case.name}/outputs_qc_only",
                bool(qc_html) and not typing_html and not leftovers and not typing_tags,
                f"QC html: {len(qc_html)}, typing html: {len(typing_html)}, "
                f"assembly/typing dirs: {len(leftovers)}, typing deliverable tags: {sorted(typing_tags)}",
            )

            after = self.ctx.db.snapshot_project(case.project)
            self.snapshot(f"{case.name}/after_qc", after)
            before = self.get_snapshot(f"{case.name}/after_full")
            problems = (
                diff_qc(before, after)
                if before
                else ["no snapshot after full analysis"]
            )
            self.check(
                f"{case.name}/qc_metrics_concordant",
                not problems,
                "identical within 1%" if not problems else " | ".join(problems[:5]),
            )

    def render(self) -> list[str]:
        lines = table(
            [
                "CaseID",
                "Runfile has only QC sections (Y/N)",
                "Mailjob uses QC track (Y/N)",
                "Finish produces only QC outputs (Y/N)",
                "QC metrics concordant with full analysis (Y/N)",
            ],
            [
                [
                    case.name,
                    yn(self.get_check(f"{case.name}/runfile_qc_only")),
                    yn(self.get_check(f"{case.name}/mailjob_qc_track")),
                    yn(self.get_check(f"{case.name}/outputs_qc_only")),
                    yn(self.get_check(f"{case.name}/qc_metrics_concordant")),
                ]
                for case in self.cases
            ],
        )
        for case in self.cases:
            check = self.get_check(f"{case.name}/qc_metrics_concordant")
            if check and not check["passed"]:
                lines += ["", f"QC metric differences for {case.name}: {detail(check)}"]
        return lines
