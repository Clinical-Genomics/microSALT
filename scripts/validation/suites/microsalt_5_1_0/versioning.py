"""Validate the shared project report sequence introduced in microSALT 5.1.0."""

import hashlib
import json
import re
from pathlib import Path

from bs4 import BeautifulSoup

from core.report import detail, table, yn
from core.slurm import run_complete
from core.suite import RetryAction, Suite, SuiteOptions


class VersioningOptions(SuiteOptions):
    case: str


def html_version(path: Path) -> tuple[int, int | None] | None:
    """Extract the current and replaced versions from a rendered report."""
    if not path.is_file():
        return None
    text = BeautifulSoup(path.read_text(encoding="utf-8"), "html.parser").get_text(" ")
    match = re.search(
        r"\bVersion\s+(\d+)\s*,.*?(?:Ersätter\s+version\s+(\d+))?", text, re.I
    )
    if not match:
        return None
    previous = re.search(r"Ersätter\s+version\s+(\d+)", text, re.I)
    return int(match.group(1)), int(previous.group(1)) if previous else None


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class VersioningSuite(Suite):
    name = "versioning"
    title = "Report versioning"
    requires = ("stability", "qc_only")
    Options = VersioningOptions
    options: VersioningOptions

    LABELS = (
        ("default_first", "First default report allocates one shared version"),
        (
            "default_repeat",
            "Unchanged default report allocates the next shared version",
        ),
        ("qc_only", "QC-only report uses the next shared project version"),
        ("typing", "Typing-only report uses the next shared project version"),
        ("qc", "QC-only report type uses the next shared project version"),
        ("json_dump", "JSON-only report does not allocate a version"),
        (
            "previous_reports_preserved",
            "Previous report rows and generated HTML remain unchanged",
        ),
        ("gen_reportdata_filter", "Sample-ID filter returns only the selected samples"),
        ("gen_reportdata_all", "Unfiltered report data contains all project samples"),
    )

    def retry_actions(self) -> list[RetryAction]:
        return [
            RetryAction(
                "finish/versioning",
                f"Repeat report-versioning sequence: {self.options.case}",
                "finish",
            )
        ]

    def retry(self, action_id: str) -> None:
        case = self.options.case
        for run_id in (f"submit/stability/{case}", f"resubmit/qc_only/{case}"):
            record = self.state.run(run_id)
            if not record or not run_complete(record):
                raise ValueError(
                    f"{run_id} must complete before retrying report versioning"
                )
        self.state.invalidate_checks("versioning/", "Retry in progress")
        self.finish()

    def finish(self) -> None:
        case = self.cfg.case(self.options.case)
        db = self.ctx.db
        full_run = self.state.run(f"submit/stability/{case.name}")
        qc_run = self.state.run(f"resubmit/qc_only/{case.name}")
        if not full_run or not full_run.finishdir or not qc_run or not qc_run.finishdir:
            for key, _ in self.LABELS:
                self.check(
                    key,
                    None,
                    "Missing full or QC-only analysis; cannot run version checks",
                )
            return
        full_dir, qc_dir = Path(full_run.finishdir), Path(qc_run.finishdir)
        with case.sampleinfo.open() as stream:
            info = json.load(stream)
        sample = info[0] if isinstance(info, list) else info
        customer = sample["Customer_ID_project"]
        initial = db.report_versions(case.project)
        history = []
        preserved_files = {}
        reports_dir = Path(self.cfg.microsalt.folders.reports) / "analysis"
        prior_html = {}
        for directory in (reports_dir, full_dir, qc_dir):
            if directory.is_dir():
                for path in directory.glob(f"{customer}_*.html"):
                    if "_Typing_" in path.name or "_QC_" in path.name:
                        prior_html[str(path)] = file_hash(path)
        print(f"== Report versioning: {case.name} ==")

        def run_step(key, folder, report, expected_types, qc_only=False):
            before = db.report_versions(case.project)
            prior = before[-1]["version"] if before else 0
            expected = prior + 1 if expected_types else prior
            json_path = folder / f"{case.project}.json"
            json_mtime_before = (
                json_path.stat().st_mtime_ns
                if not expected_types and json_path.exists()
                else None
            )
            # A pre-existing path must not be mistaken for output from this invocation.
            paths = [
                folder / f"{customer}_{kind}_{expected}.html" for kind in expected_types
            ]
            existed = [str(path) for path in paths if path.exists()]
            record, _ = self.finish_run(
                f"{case.name}/{key}",
                case.name,
                case.sampleinfo,
                folder,
                report=report,
                qc_only=qc_only,
            )
            after = db.report_versions(case.project)
            versions = [item["version"] for item in after]
            expected_versions = [item["version"] for item in before]
            if expected_types:
                expected_versions.append(expected)
            problems = []
            if not record.ok:
                problems.append(f"finish exited {record.exit_code}; see {record.log}")
            if versions != expected_versions:
                problems.append(
                    f"versions {expected_versions} expected, got {versions}"
                )
            if after[: len(before)] != before:
                problems.append("previous database report records changed")
            if existed:
                problems.append(f"output already existed before finish: {existed}")
            for path in paths:
                parsed = html_version(path)
                if parsed != (expected, prior or None):
                    problems.append(
                        f"{path.name}: expected version/previous {(expected, prior or None)}, got {parsed}"
                    )
                if path.is_file():
                    preserved_files[str(path)] = file_hash(path)
            if not expected_types:
                if not json_path.is_file():
                    problems.append(f"JSON output missing: {json_path}")
                elif (
                    json_mtime_before is not None
                    and json_path.stat().st_mtime_ns == json_mtime_before
                ):
                    problems.append(
                        f"JSON output was not refreshed by this invocation: {json_path}"
                    )
            self.check(
                key,
                not problems,
                (
                    "; ".join(problems)
                    if problems
                    else f"versions {prior} -> {expected}; report={report}; files: {[p.name for p in paths]}"
                ),
            )
            history.append(
                {
                    "step": key,
                    "before": prior,
                    "after": versions[-1] if versions else 0,
                    "exit_code": record.exit_code,
                    "files": [str(path) for path in paths],
                }
            )

        run_step("default_first", full_dir, "default", ("Typing", "QC"))
        run_step("default_repeat", full_dir, "default", ("Typing", "QC"))
        run_step("qc_only", qc_dir, "default", ("QC",), qc_only=True)
        run_step("typing", full_dir, "typing", ("Typing",))
        run_step("qc", full_dir, "qc", ("QC",))
        run_step("json_dump", full_dir, "json_dump", ())
        final = db.report_versions(case.project)
        tracked_html = {**prior_html, **preserved_files}
        changed = [
            path
            for path, digest in tracked_html.items()
            if not Path(path).is_file() or file_hash(Path(path)) != digest
        ]
        preserved = final[: len(initial)] == initial and not changed
        self.check(
            "previous_reports_preserved",
            preserved,
            f"initial DB rows: {len(initial)}; HTML files checked: {len(tracked_html)}; "
            f"changed or missing HTML: {changed}",
        )
        self.snapshot(
            "history",
            {
                "steps": history,
                "initial_reports": initial,
                "final_reports": final,
                "prior_html_sha256": prior_html,
                "html_sha256": preserved_files,
            },
        )

        all_ids = db.report_sample_ids(case.project, None)
        expected_ids = db.project_sample_ids(case.project)
        selected = expected_ids[:1]
        subset_ids = db.report_sample_ids(case.project, selected)
        self.check(
            "gen_reportdata_filter",
            bool(selected) and subset_ids == selected,
            f"sample_ids={selected} -> {subset_ids}",
        )
        self.check(
            "gen_reportdata_all",
            bool(expected_ids) and all_ids == expected_ids,
            f"expected {expected_ids}; returned {all_ids}",
        )

    def result(self) -> str:
        checks = [self.get_check(key) for key, _ in self.LABELS]
        if any(check and check["passed"] is False for check in checks):
            return "FAIL"
        if any(not check or check["passed"] is None for check in checks):
            return "INCOMPLETE"
        return "PASS"

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
                "Versions by step: "
                + ", ".join(
                    f"{step['step']}: {step['before']} → {step['after']}"
                    for step in history["steps"]
                ),
            ]
        return lines
