"""QC-only runs must not alter stored typing data, must refresh metadata, and a later full finish
must still carry the typing data."""

from pathlib import Path

from core.db import diff_typing
from core.report import detail, table, yn
from core.slurm import run_complete
from core.suite import RetryAction, Suite, SuiteOptions


class QcIntegrityOptions(SuiteOptions):
    case: str
    refresh_field: str = "priority"
    refresh_value: str = "research"


class QcIntegritySuite(Suite):
    name = "qc_integrity"
    title = "QC-only data integrity"
    requires = ("stability", "qc_only")
    Options = QcIntegrityOptions
    options: QcIntegrityOptions

    LABELS = (
        (
            "typing_unchanged",
            "Sequence types, alleles and resistances unchanged after QC-only run",
        ),
        ("metadata_refreshed", "Changed metadata refreshed by QC-only finish"),
        (
            "full_report_retains_typing",
            "Full finish after QC-only run still contains typing data",
        ),
        (
            "metadata_restored",
            "Metadata restored by full finish with original sample info",
        ),
    )

    def retry_actions(self) -> list[RetryAction]:
        return [
            RetryAction(
                "finish/qc_integrity",
                f"Repeat QC integrity sequence: {self.options.case}",
                "finish",
            )
        ]

    def retry(self, action_id: str) -> None:
        case = self.options.case
        for run_id in (f"submit/stability/{case}", f"resubmit/qc_only/{case}"):
            record = self.state.run(run_id)
            if not record or not run_complete(record):
                raise ValueError(f"{run_id} must complete before retrying QC integrity")
        qc_finish = self.state.run(f"finish/qc_only/{case}")
        if (
            not qc_finish
            or not qc_finish.ok
            or not self.state.get_snapshot(f"qc_only/{case}/after_qc")
        ):
            raise ValueError(
                f"QC finish for {case} must succeed before retrying QC integrity"
            )
        self.state.invalidate_checks("qc_integrity/", "Retry in progress")
        with self.fixtures.retry_fixtures():
            self.finish()

    def finish(self) -> None:
        case = self.cfg.case(self.options.case)
        field, value = self.options.refresh_field, self.options.refresh_value
        print(f"== QC-only data integrity: {case.name} ==")
        before = self.state.get_snapshot(f"qc_only/{case.name}/after_full")
        after = self.state.get_snapshot(
            f"qc_only/{case.name}/after_qc"
        ) or self.ctx.db.snapshot_project(case.project)
        if before is None:
            self.check(
                "typing_unchanged",
                None,
                "Missing snapshot after the first full analysis",
            )
            return
        problems = diff_typing(before, after)
        self.check(
            "typing_unchanged",
            not problems,
            (
                "ST, alleles and resistances unchanged"
                if not problems
                else " | ".join(problems[:5])
            ),
        )

        qc_run = self.state.run(f"resubmit/qc_only/{case.name}")
        full_run = self.state.run(f"submit/stability/{case.name}")
        if qc_run is None or not qc_run.finishdir:
            self.check("metadata_refreshed", None, "Missing QC-only finish directory")
            return
        if full_run is None or not full_run.finishdir:
            self.check(
                "full_report_retains_typing",
                None,
                "Missing full-analysis finish directory",
            )
            return
        modified, first = self.fixtures.modified_sampleinfo(
            case,
            "refresh",
            index=0,
            **{field: value},
        )
        sample_id = first["CG_ID_sample"]
        original = before["samples"][sample_id]["metadata"].get(field)
        self.finish_run(
            "refresh",
            case.name,
            modified,
            Path(qc_run.finishdir),
            qc_only=True,
            note=f"{field} of {sample_id} set to {value!r}",
        )
        sample = self.ctx.db.sample(sample_id)
        self.check(
            "metadata_refreshed",
            sample is not None and getattr(sample, field) == value,
            f"{field}: {original!r} -> {getattr(sample, field, None)!r}",
        )

        self.finish_run(
            "full",
            case.name,
            case.sampleinfo,
            Path(full_run.finishdir),
            note="full finish on the full-analysis folder with the original sample info",
        )
        final = self.ctx.db.snapshot_project(case.project)
        self.snapshot("after_full_refinish", final)
        problems = diff_typing(before, final)
        has_typing = any(s["seq_types"] for s in final["samples"].values())
        self.check(
            "full_report_retains_typing",
            not problems and has_typing,
            (
                "typing data identical to the first full analysis"
                if not problems and has_typing
                else " | ".join(problems[:5]) or "no seq_types rows"
            ),
        )
        sample = self.ctx.db.sample(sample_id)
        self.check(
            "metadata_restored",
            sample is not None and getattr(sample, field) == original,
            f"{field} back to {original!r}",
        )

    def render(self) -> list[str]:
        return table(
            ["CaseID", "Check", "Result (Y/N)", "Detail"],
            [
                [
                    self.options.case,
                    label,
                    yn(self.get_check(key)),
                    detail(self.get_check(key)),
                ]
                for key, label in self.LABELS
            ],
        )
