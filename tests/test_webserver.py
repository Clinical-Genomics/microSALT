#!/usr/bin/env python

import json
import os
import pathlib
import pdb
import pytest
import requests
import time

from distutils.sysconfig import get_python_lib
from unittest.mock import patch

from microSALT.utils.reporter import Reporter
from microSALT import preset_config, logger
from microSALT.cli import root

@pytest.fixture
def testdata():
  testdata = os.path.abspath(os.path.join(pathlib.Path(__file__).parent.parent, 'tests/testdata/sampleinfo_samples.json'))
  #Check if release install exists
  for entry in os.listdir(get_python_lib()):
    if 'microSALT-' in entry:
      testdata = os.path.abspath(os.path.join(os.path.expandvars('$CONDA_PREFIX'), 'testdata/sampleinfo_samples.json'))
  with open(testdata) as json_file:
    data = json.load(json_file)
  return data

@pytest.fixture
def report_obj(testdata):
  report = Reporter(config=preset_config, log=logger, sampleinfo=testdata)
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
  assert e.status_code in [200, 500]

  time.sleep(0.15)
  d = requests.get("http://127.0.0.1:5000/microSALT/AAA1234/qc", allow_redirects=True)
  assert d.status_code in [200, 500]

  #Valid pages with unavailable data
  time.sleep(0.15)
  f = requests.get("http://127.0.0.1:5000/microSALT/AAA1234/typing/escherichia_coli", allow_redirects=True)
  assert f.status_code == 500

  time.sleep(0.15)
  g = requests.get("http://127.0.0.1:5000/microSALT/STtracker/all", allow_redirects=True)
  assert g.status_code in [200, 500]

  time.sleep(0.15)
  h = requests.get("http://127.0.0.1:5000/microSALT/STtracker/cust000", allow_redirects=True)
  assert h.status_code in [200, 500]

  report_obj.kill_flask()
