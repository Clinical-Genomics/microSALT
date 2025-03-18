#!/usr/bin/env python

import glob
import logging
import pytest

from microSALT.utils.reporter import Reporter
from microSALT.store.orm_models import Reports

logger = logging.getLogger("main_logger")


def test_motif(reporter):
    reporter.create_subfolders()
    reporter.gen_motif(motif="resistance")
    assert len(glob.glob("{}/AAA1234_resistance*".format(reporter.output))) > 0

    reporter.gen_motif(motif="expec")
    assert len(glob.glob("{}/AAA1234_expec*".format(reporter.output))) > 0


def test_deliveryreport(config, reporter):
    reporter.create_subfolders()
    reports = reporter.db_pusher.session.query(Reports).all()
    reporter.gen_delivery()
    assert (
        len(
            glob.glob(
                "{}/deliverables/999999_deliverables.yaml".format(config["folders"]["reports"])
            )
        )
        > 0
    )


def test_jsonreport(config, reporter):
    reporter.create_subfolders()
    reporter.gen_json()
    assert len(glob.glob("{}/json/AAA1234.json".format(config["folders"]["reports"]))) > 0


def test_gen_qc(reporter):
    reporter.name = "name_that_do_not_exist"
    with pytest.raises(AttributeError):
        reporter.gen_qc()


def test_gen_typing(reporter):
    reporter.name = "name_that_do_not_exist"
    with pytest.raises(Exception):
        reporter.gen_typing()


def test_gen_motif(caplog, reporter):
    caplog.clear()
    reporter.gen_motif(motif="unrecognized")
    assert "Invalid motif type" in caplog.text
    caplog.clear()
    reporter.output = "/path/that/do/not/exists/"
    reporter.gen_motif()
    assert "Gen_motif unable to produce" in caplog.text


def test_gen_json(caplog, config, reporter):
    caplog.clear()
    reporter.output = "/path/that/do/not/exists/"
    config["folders"]["reports"] = "/path/that/do/not/exists/"
    reporter.config = config
    reporter.gen_json()
    assert "Gen_json unable to produce" in caplog.text


def test_report(caplog, reporter):
    caplog.clear()
    reporter.type = "type_not_mentioned_in_list"
    with pytest.raises(Exception):
        reporter.report()
        assert "Report function recieved invalid format" in caplog.text


def test_restart_web(reporter):
    reporter.restart_web()


def test_constructor(app, config, dbm, logger, sampleinfo_sample):
    sample_info = [
        {
            "CG_ID_project": "AAA1234",
            "CG_ID_sample": "AAA1234A1",
            "Customer_ID_project": "999999",
            "Customer_ID_sample": "XXX0000Y1",
            "Customer_ID": "cust000",
            "application_tag": "NONE",
            "date_arrival": "0001-01-01 00:00:00",
            "date_libprep": "0001-01-01 00:00:00",
            "date_sequencing": "0001-01-01 00:00:00",
            "method_libprep": "Not in LIMS",
            "method_sequencing": "Not in LIMS",
            "organism": "Staphylococcus aureus",
            "priority": "standard",
            "reference": "AP017922.1",
        },
        {
            "CG_ID_project": "Something_else_than_AAA1234",
            "CG_ID_sample": "AAA1234A2",
            "Customer_ID_project": "999999",
            "Customer_ID_sample": "XXX0000Y1",
            "Customer_ID": "cust000",
            "application_tag": "NONE",
            "date_arrival": "0001-01-01 00:00:00",
            "date_libprep": "0001-01-01 00:00:00",
            "date_sequencing": "0001-01-01 00:00:00",
            "method_libprep": "Not in LIMS",
            "method_sequencing": "Not in LIMS",
            "organism": "Escherichia coli",
            "priority": "standard",
            "reference": "NC_011751.1",
        },
    ]
    reporter_obj = Reporter(
        app=app,
        config=config,
        dbm=dbm,
        log=logger,
        sampleinfo=sampleinfo_sample,
        name="MIC1234A1",
        output="/tmp/MLST",
    )
