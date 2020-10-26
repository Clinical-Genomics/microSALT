#!/usr/bin/env python

import glob
import json
import logging
import os
import pathlib
import pdb
import pytest

from distutils.sysconfig import get_python_lib

from microSALT import preset_config, logger
from microSALT.utils.reporter import Reporter
from microSALT.utils.referencer import Referencer


@pytest.fixture
def testdata_prefix():
  test_path = os.path.abspath(os.path.join(pathlib.Path(__file__).parent.parent, 'tests/testdata/'))
  #Check if release install exists
  for entry in os.listdir(get_python_lib()):
    if 'microSALT-' in entry:
      test_path = os.path.abspath(os.path.join(os.path.expandvars('$CONDA_PREFIX'), 'testdata/'))
  return test_path

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
def init_references(testdata):
  ref_obj = Referencer(config=preset_config, log=logger, sampleinfo=testdata)
  ref_obj.identify_new(testdata[0].get('CG_ID_project'),project=True)
  ref_obj.update_refs()

@pytest.fixture
def reporter(testdata, init_references):
  reporter_obj = Reporter(config=preset_config, log=logger, sampleinfo=testdata[0], name="MIC1234A1", output="") 
  return reporter_obj

def test_motif(reporter, init_references):
  reporter.gen_motif(motif="resistance")
  reporter.gen_motif(motif="expec")
