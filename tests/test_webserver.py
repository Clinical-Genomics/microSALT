#!/usr/bin/env python

import pytest
import requests
import time
from unittest.mock import patch

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

def test_pages(report_obj):
  report_obj.start_web()
  #Valid pages with available data
  time.sleep(3)
  a = requests.get("http://127.0.0.1:5000/", allow_redirects=True)
  assert a.status_code == 200
  
  time.sleep(0.15)
  b = requests.get("http://127.0.0.1:5000/microSALT/", allow_redirects=True)
  assert b.status_code == 200

  time.sleep(0.15)
  c = requests.get("http://127.0.0.1:5000/microSALT/AAA1234", allow_redirects=True)
  assert c.status_code == 200

  time.sleep(0.15)
  e = requests.get("http://127.0.0.1:5000/microSALT/AAA1234/typing/all", allow_redirects=True)
  assert e.status_code == 200

  #Valid pages with unavailable data
  time.sleep(0.15)
  d = requests.get("http://127.0.0.1:5000/microSALT/AAA1234/qc", allow_redirects=True)
  assert d.status_code == 500

  time.sleep(0.15)
  f = requests.get("http://127.0.0.1:5000/microSALT/AAA1234/typing/escherichia_coli", allow_redirects=True)
  assert f.status_code == 500

  time.sleep(0.15)
  g = requests.get("http://127.0.0.1:5000/microSALT/STtracker/all", allow_redirects=True)
  assert g.status_code == 500

  time.sleep(0.15)
  h = requests.get("http://127.0.0.1:5000/microSALT/STtracker/cust000", allow_redirects=True)
  assert h.status_code == 500

  report_obj.kill_flask()

