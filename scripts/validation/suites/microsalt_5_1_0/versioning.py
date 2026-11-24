"""Report versioning follows the rendered roster and fingerprinted metadata; gen_reportdata
filters by sample id."""

from pathlib import Path

from core.report import detail, table, yn
from core.suite import Suite, SuiteOptions


class VersioningOptions(SuiteOptions):
    case: str
    version_field: str = "method_libprep"  # part of the report fingerprint
    version_value: str = "validation"


class VersioningSuite(Suite):
    name = "versioning"
    title = "Report versioning"
    requires = ("stability",)
    Options = VersioningOptions

    LABELS = (
        (
            "repeat_no_new_version",
            "Repeated finish without changes gives no new version",
        ),
        ("metadata_new_version", "Changed sample metadata gives new version"),
        (
            "qc_same_roster_no_new_version",
            "QC-only report over same samples/metadata gives no new version",
        ),
        ("subset_new_version", "Changed set of included samples gives new version"),
        (
            "gen_reportdata_filter",
            "gen_reportdata with sample_ids returns only those samples",
        ),
        (
            "gen_reportdata_all",
            "gen_reportdata without filter returns all project samples",
        ),
    )

    def finish(self) -> None:
        case = self.cfg.case(self.options.case)
        field, value = self.options.version_field, self.options.version_value
        db = self.ctx.db
        print(f"== Report versioning: {case.name} ==")
        full_dir = Path(self.state.run(f"submit/stability/{case.name}").finishdir)
        qc_run = self.state.run(f"resubmit/qc_only/{case.name}")
        history = []

        def versions(label: str) -> int:
            count = len(db.report_versions(case.project))
            history.append((label, count))
            return count

        def finish(
            label: str, sampleinfo: Path, finishdir: Path, qc_only=False, note=""
        ) -> int:
            self.finish_run(
                f"{case.name}/{label}",
                case.name,
                sampleinfo,
                finishdir,
                qc_only=qc_only,
                note=note,
            )
            return versions(label)

        versions("start")
        v1 = finish("full_1", case.sampleinfo, full_dir)
        v2 = finish("full_2", case.sampleinfo, full_dir)
        self.check(
            "repeat_no_new_version",
            v2 == v1,
            f"versions after first finish {v1}, after repeat {v2}",
        )

        modified, first = self.fixtures.modified_sampleinfo(
            case, "version", **{field: value}
        )
        v3 = finish(
            "full_metadata",
            modified,
            full_dir,
            note=f"{field} of {first['CG_ID_sample']} set to {value!r}",
        )
        self.check(
            "metadata_new_version",
            v3 == v2 + 1,
            f"versions {v2} -> {v3} after changing {field}",
        )
        v4 = finish(
            "full_restore",
            case.sampleinfo,
            full_dir,
            note="original sample info restored",
        )

        if qc_run:
            qc_dir = Path(qc_run.finishdir)
            v5 = finish("qc_same_roster", case.sampleinfo, qc_dir, qc_only=True)
            self.check(
                "qc_same_roster_no_new_version",
                v5 == v4,
                f"versions {v4} -> {v5} for a QC-only report over the same samples and metadata",
            )
            subset, ids = self.fixtures.subset_sampleinfo(case)
            v6 = finish("qc_subset", subset, qc_dir, qc_only=True, note=f"only {ids}")
            self.check(
                "subset_new_version",
                v6 == v5 + 1,
                f"versions {v5} -> {v6} for a QC-only report over {ids}",
            )
            finish("qc_restore", case.sampleinfo, qc_dir, qc_only=True)
        else:
            for key in ("qc_same_roster_no_new_version", "subset_new_version"):
                self.check(key, None, "qc_only suite not run for this case")
        self.snapshot("history", {"steps": history})

        all_ids = db.report_sample_ids(case.project, None)
        subset_ids = db.report_sample_ids(case.project, all_ids[:1])
        self.check(
            "gen_reportdata_filter",
            subset_ids == all_ids[:1],
            f"sample_ids={all_ids[:1]} -> {subset_ids}",
        )
        expected = db.project_sample_ids(case.project)
        self.check(
            "gen_reportdata_all",
            all_ids == expected,
            f"{len(all_ids)} of {len(expected)} project samples returned without filter",
        )

    def render(self) -> list[str]:
        lines = table(
            ["CaseID", "Check", "Result (Y/N)", "Comment"],
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
        history = self.get_snapshot("history")
        if history:
            lines += [
                "",
                f"Report version count for {self.options.case} after each step: "
                + ", ".join(f"{label}={count}" for label, count in history["steps"]),
            ]
        return lines
