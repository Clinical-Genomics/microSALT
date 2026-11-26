"""Regression tests for the release validation's project report version checks."""

import copy
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts" / "validation"))

from core.state import RunRecord, State  # noqa: E402
from suites.microsalt_5_1_0.versioning import (  # noqa: E402
    VersioningSuite,
    html_version,
)


def test_html_version_parses_current_and_replaced(tmp_path):
    report = tmp_path / "report.html"
    report.write_text(
        "<td>Version 8, 2026-01-01 (Ers&auml;tter version 7, 2025-01-01)</td>"
    )
    assert html_version(report) == (8, 7)
    report.write_text("<td>Version 1, 2026-01-01</td>")
    assert html_version(report) == (1, None)
    report.write_text("<td>no version</td>")
    assert html_version(report) is None


@pytest.mark.parametrize(
    "fault",
    [
        None,
        "no_increment",
        "failed_finish",
        "missing_html",
        "overwrite_previous",
        "json_increment",
    ],
)
def test_versioning_checks(tmp_path, fault):
    full = tmp_path / "full"
    qc = tmp_path / "qc"
    full.mkdir()
    qc.mkdir()
    info = tmp_path / "samples.json"
    info.write_text(
        json.dumps([{"Customer_ID_project": "CUST", "CG_ID_sample": "PROJA1"}])
    )
    case = SimpleNamespace(name="case", project="PROJ", sampleinfo=info)
    report_dir = tmp_path / "reports" / "analysis"
    report_dir.mkdir(parents=True)
    (report_dir / "CUST_Typing_5.html").write_text("old typing report")
    (report_dir / "CUST_QC_5.html").write_text("old QC report")
    state = State(tmp_path / "state.json")
    for key, folder in [("submit/stability/case", full), ("resubmit/qc_only/case", qc)]:
        state.add_run(
            RunRecord(
                id=key,
                phase="submit",
                test="test",
                case="case",
                command="",
                started="",
                duration_s=0,
                exit_code=0,
                log="",
                finishdir=str(folder),
            )
        )

    class FakeDB:
        def __init__(self):
            self.records = [
                {"version": 5, "fingerprint": "legacy", "date": "yesterday"}
            ]

        def report_versions(self, project):
            return copy.deepcopy(self.records)

        def project_sample_ids(self, project):
            return ["PROJA1", "PROJA2"]

        def report_sample_ids(self, project, ids):
            return ids if ids is not None else self.project_sample_ids(project)

    db = FakeDB()
    ctx = SimpleNamespace(
        cfg=SimpleNamespace(
            case=lambda name: case,
            microsalt=SimpleNamespace(
                folders=SimpleNamespace(reports=str(tmp_path / "reports"))
            ),
        ),
        state=state,
        runner=None,
        fixtures=None,
        db=db,
    )
    suite = VersioningSuite(ctx, {"case": "case"})

    def finish(label, case_name, sampleinfo, folder, report, qc_only=False):
        step = label.split("/")[-1]
        previous = db.records[-1]["version"]
        increments = report != "json_dump" and not (
            fault == "no_increment" and step == "default_repeat"
        )
        if increments or (fault == "json_increment" and report == "json_dump"):
            db.records.append(
                {"version": previous + 1, "fingerprint": None, "date": "today"}
            )
        current = db.records[-1]["version"]
        kinds = {
            "default": ["QC"] if qc_only else ["QC", "Typing"],
            "qc": ["QC"],
            "typing": ["Typing"],
            "json_dump": [],
        }[report]
        for kind in kinds:
            if fault == "missing_html" and step == "qc_only":
                continue
            (folder / f"CUST_{kind}_{current}.html").write_text(
                f"Version {current}, today (Ersätter version {previous}, yesterday)"
            )
        if report == "json_dump":
            (folder / "PROJ.json").write_text("{}")
        if fault == "overwrite_previous" and step == "qc":
            (full / "CUST_Typing_6.html").write_text("overwritten")
        ok = not (fault == "failed_finish" and step == "typing")
        return SimpleNamespace(ok=ok, exit_code=0 if ok else 1, log="failure.log"), None

    suite.finish_run = finish
    suite.finish()
    checks = state.data["checks"]
    if fault is None:
        assert suite.result() == "PASS"
        assert all(checks[f"versioning/{key}"]["passed"] for key, _ in suite.LABELS)
        assert [x["version"] for x in db.records] == [5, 6, 7, 8, 9, 10]
        assert db.records[0]["fingerprint"] == "legacy"
    else:
        assert suite.result() == "FAIL"
        key = {
            "no_increment": "default_repeat",
            "failed_finish": "typing",
            "missing_html": "qc_only",
            "overwrite_previous": "previous_reports_preserved",
            "json_increment": "json_dump",
        }[fault]
        assert checks[f"versioning/{key}"]["passed"] is False
