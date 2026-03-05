#!/usr/bin/env python

import glob
import pytest

from microSALT.utils.reporter import Reporter


def test_motif(dbm, reporter):
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
    reporter.config.folders.reports = "/path/that/do/not/exists/"
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
        config=config,
        log=logger,
        sampleinfo=sample_info,
        name="MIC1234A1",
        output="/tmp/MLST",
    )
