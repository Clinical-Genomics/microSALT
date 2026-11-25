#!/usr/bin/env python
import glob
import json
from types import SimpleNamespace
from unittest.mock import patch

import pytest
import yaml

from microSALT.store.db_manipulator import DB_Manipulator
from microSALT.utils.reporter import Reporter


def test_motif(dbm: DB_Manipulator, reporter: Reporter):
    reporter.create_subfolders()
    reporter.gen_motif(motif="resistance")
    assert len(glob.glob(f"{reporter.output}/AAA1234_resistance*")) > 0

    reporter.gen_motif(motif="expec")
    assert len(glob.glob(f"{reporter.output}/AAA1234_expec*")) > 0


def test_deliveryreport(config, dbm, reporter):
    reporter.create_subfolders()
    reporter.gen_delivery()
    assert (
        len(
            glob.glob(
                f"{config.folders.reports}/deliverables/999999_deliverables.yaml"
            )
        )
        > 0
    )


def test_jsonreport(config, dbm, reporter):
    reporter.create_subfolders()
    reporter.gen_json()
    assert len(glob.glob(f"{config.folders.reports}/json/AAA1234.json")) > 0


@pytest.mark.parametrize("qc_only", [False, True])
def test_deliveryreport_qc_only_filters_unproduced_files(config, dbm, reporter, qc_only):
    """QC-only deliverables must not list assembly, QUAST or typing outputs."""
    reporter.qc_only = qc_only
    reporter.create_subfolders()
    reporter.gen_delivery()

    with open(f"{config.folders.reports}/deliverables/999999_deliverables.yaml") as fh:
        tags = {f["tag"] for f in yaml.safe_load(fh)["files"]}

    excluded = {"microsalt-type", "assembly", "quast-results"}
    assert tags.isdisjoint(excluded) is qc_only
    assert {"microsalt-qc", "microsalt-json", "reference-alignment-deduplicated"} <= tags


@pytest.mark.parametrize("qc_only", [False, True])
def test_jsonreport_qc_only_omits_typing_and_assembly(config, reporter, qc_only):
    """QC-only JSON keeps alignment metrics and drops typing/assembly/resistance sections."""
    metrics = (
        "genome_length reference_length gc_percentage n50 contigs insert_size duplication_rate "
        "total_reads raw_reads mapped_rate average_coverage coverage_10x coverage_30x "
        "coverage_50x coverage_100x"
    ).split()
    sample = SimpleNamespace(
        CG_ID_sample="AAA1234A1",
        ST=-1,
        ST_status="Data saknas",
        threshold="-",
        resistances=[],
        **dict.fromkeys(metrics),
    )
    reporter.qc_only = qc_only
    reporter.create_subfolders()
    with patch(
        "microSALT.utils.reporter.gen_reportdata", return_value={"samples": [sample]}
    ) as report_data:
        reporter.gen_json()

    with open(f"{config.folders.reports}/json/{reporter.name}.json") as fh:
        analyses = json.load(fh)["AAA1234A1"]

    assert {"picard_markduplicate", "microsalt_samtools_stats"} <= analyses.keys()
    excluded = {"blast_pubmlst", "quast_assembly", "blast_resfinder_resistence"}
    assert analyses.keys().isdisjoint(excluded) is qc_only
    assert report_data.call_args.kwargs["sample_ids"] == (
        ["AAA1234A1"] if qc_only else None
    )


def test_reporter_sample_ids_follow_current_sample_info(reporter):
    reporter.sampleinfo = [
        {"CG_ID_sample": "AAA1234A1"},
        {"CG_ID_sample": "AAA1234A2"},
    ]
    assert reporter._sample_ids() == ["AAA1234A1", "AAA1234A2"]


def test_gen_version_allocates_project_version(reporter):
    with patch.object(reporter.db_pusher, "set_report") as set_report:
        reporter.gen_version(reporter.name)

    set_report.assert_called_once_with(reporter.name)


@pytest.mark.parametrize("qc_only", [False, True])
def test_report_qc_generates_json_and_deliverables_only_when_qc_only(reporter, qc_only):
    """A QC-only `qc` report also writes JSON + deliverables; a plain `qc` report does not."""
    reporter.qc_only = qc_only
    with (
        patch.object(reporter, "db_pusher"),
        patch.object(reporter, "gen_qc") as gen_qc,
        patch.object(reporter, "gen_json") as gen_json,
        patch.object(reporter, "gen_delivery") as gen_delivery,
        patch.object(reporter, "gen_typing") as gen_typing,
        patch.object(reporter, "mail"),
    ):
        reporter.report(type="qc")

    gen_qc.assert_called_once_with()
    gen_typing.assert_not_called()
    assert gen_json.call_count == gen_delivery.call_count == int(qc_only)


@pytest.mark.parametrize("qc_only", [False, True])
def test_default_report_skips_typing_for_qc_only(reporter, qc_only):
    reporter.qc_only = qc_only
    with (
        patch.object(reporter, "db_pusher"),
        patch.object(reporter, "gen_qc") as gen_qc,
        patch.object(reporter, "gen_json") as gen_json,
        patch.object(reporter, "gen_delivery") as gen_delivery,
        patch.object(reporter, "gen_typing") as gen_typing,
        patch.object(reporter, "mail") as mail,
    ):
        reporter.report()

    assert gen_typing.call_count == int(not qc_only)
    gen_qc.assert_called_once_with()
    gen_json.assert_called_once_with(silent=True)
    gen_delivery.assert_called_once_with()
    mail.assert_called_once_with()


def test_gen_qc_name_does_not_exist(dbm, reporter):
    reporter.name = "name_that_do_not_exist"
    with pytest.raises(SystemExit):
        reporter.gen_qc()


def test_gen_typing_name_does_not_exist(dbm, reporter):
    reporter.name = "name_that_do_not_exist"
    with pytest.raises(SystemExit):
        reporter.gen_typing()


def test_gen_motif(caplog, reporter):
    caplog.clear()
    reporter.gen_motif(motif="unrecognized")
    assert "Invalid motif type" in caplog.text
    caplog.clear()
    reporter.output = "/path/that/do/not/exists/"
    reporter.gen_motif()
    assert "Gen_motif unable to produce" in caplog.text


def test_gen_json(caplog, reporter):
    caplog.clear()
    reporter.output = "/path/that/do/not/exists/"
    reporter.folders.reports = "/path/that/do/not/exists/"
    reporter.gen_json()
    assert "Gen_json unable to produce" in caplog.text


def test_report(caplog, reporter):
    caplog.clear()
    reporter.type = "type_not_mentioned_in_list"
    with pytest.raises(Exception):
        reporter.report()
        assert "Report function recieved invalid format" in caplog.text


def test_constructor(config, logger, unpack_db_json):
    sample_info = unpack_db_json("sampleinfo_samples.json")
    reporter_obj = Reporter(
        log=logger,
        folders=config.folders,
        threshold=config.threshold,
        regex=config.regex,
        sampleinfo=sample_info,
        name="MIC1234A1",
        output="/tmp/MLST",
    )


def test_project_report_versions_increment_for_qc_full_and_qc_again(reporter):
    dbm = reporter.db_pusher
    project = "VERSIONTEST"
    dbm.add_to_session(dbm.add_project(CG_ID_project=project))
    dbm.commit_session()
    reporter.name = project
    with (
        patch.object(reporter, "create_subfolders"),
        patch.object(reporter, "gen_qc"),
        patch.object(reporter, "gen_typing"),
        patch.object(reporter, "gen_json"),
        patch.object(reporter, "gen_delivery"),
        patch.object(reporter, "mail"),
    ):
        for expected, qc_only in enumerate((True, False, True), 1):
            reporter.qc_only = qc_only
            reporter.report("default")
            assert dbm.read_report(project).version == expected


def test_report_failure_still_consumes_version(reporter):
    dbm = reporter.db_pusher
    reporter.name = "VERSIONFAIL"
    dbm.add_to_session(dbm.add_project(CG_ID_project=reporter.name))
    dbm.commit_session()
    with (
        patch.object(reporter, "create_subfolders"),
        patch.object(reporter, "gen_qc", side_effect=RuntimeError("render failed")),
        pytest.raises(RuntimeError, match="render failed"),
    ):
        reporter.report("qc")
    assert dbm.read_report(reporter.name).version == 1


@pytest.mark.parametrize("report_type", ["default", "qc", "typing"])
def test_project_html_reports_allocate_one_version(reporter, report_type):
    with (
        patch.object(reporter.db_pusher, "check_ref_lock"),
        patch.object(reporter, "create_subfolders"),
        patch.object(reporter, "gen_version") as gen_version,
        patch.object(reporter, "gen_qc"),
        patch.object(reporter, "gen_typing"),
        patch.object(reporter, "gen_json"),
        patch.object(reporter, "gen_delivery"),
        patch.object(reporter, "mail"),
    ):
        reporter.report(report_type)
    gen_version.assert_called_once_with(reporter.name)


@pytest.mark.parametrize("report_type", ["json_dump", "motif_overview", "st_update"])
def test_non_html_report_types_do_not_allocate_project_version(reporter, report_type):
    with (
        patch.object(reporter.db_pusher, "check_ref_lock"),
        patch.object(reporter, "create_subfolders"),
        patch.object(reporter, "gen_version") as gen_version,
        patch.object(reporter, "gen_json"),
        patch.object(reporter, "gen_delivery"),
        patch.object(reporter, "gen_motif"),
        patch.object(reporter, "gen_STtracker"),
        patch.object(reporter, "mail"),
    ):
        reporter.report(report_type)
    gen_version.assert_not_called()
