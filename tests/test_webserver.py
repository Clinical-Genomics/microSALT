#!/usr/bin/env python

import pytest

from microSALT.utils.reporter import Reporter
from microSALT import config, logger
from microSALT.cli import root

@pytest.fixture
def report_obj():
  report = Reporter(config=config, log=logger)
  return report

def test_webserver(report_obj):
  report_obj.start_web()
  report_obj.kill_flask()
