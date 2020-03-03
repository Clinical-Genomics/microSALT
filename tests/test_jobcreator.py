#!/usr/bin/env python

import json
import mock
import os
import pathlib
import pdb
import pytest
import re

from distutils.sysconfig import get_python_lib
from unittest.mock import patch

from microSALT.store.db_manipulator import DB_Manipulator
from microSALT.utils.job_creator import Job_Creator
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

def fake_search(int):
  return "fake"

@patch('os.listdir')
@patch('os.stat')
@patch('gzip.open')
def test_verify_fastq(gopen, stat, listdir, testdata):
  listdir.return_value = ["ACC6438A3_HVMHWDSXX_L1_1.fastq.gz", "ACC6438A3_HVMHWDSXX_L1_2.fastq.gz", "ACC6438A3_HVMHWDSXX_L2_2.fastq.gz", "ACC6438A3_HVMHWDSXX_L2_2.fastq.gz"]
  stata = mock.MagicMock()
  stata.st_size = 2000
  stat.return_value = stata

  jc = Job_Creator(run_settings={'input':'/tmp/'}, config=preset_config, log=logger,sampleinfo=testdata)
  t = jc.verify_fastq()
  assert len(t) > 0

@patch('re.search')
@patch('microSALT.utils.job_creator.glob.glob')
def test_blast_subset(glob_search, research, testdata):
  jc = Job_Creator(run_settings={'input':'/tmp/'}, config=preset_config, log=logger,sampleinfo=testdata)
  researcha = mock.MagicMock()
  researcha.group = fake_search
  research.return_value = researcha
  glob_search.return_value = ["/a/a/a", "/a/a/b","/a/a/c"]

  jc.blast_subset('mlst', '/tmp/*')
  jc.blast_subset('other', '/tmp/*')
  outfile = open(jc.get_sbatch(), 'r')
  count = 0
  for x in outfile.readlines():
    if "blastn -db" in x:
      count = count + 1
  assert count > 0

@patch('subprocess.Popen')
def test_create_snpsection(subproc,testdata):
  #Sets up subprocess mocking
  process_mock = mock.Mock()
  attrs = {'communicate.return_value': ('output 123456789', 'error')}
  process_mock.configure_mock(**attrs)
  subproc.return_value = process_mock
  
  testdata = [testdata[0]]
  jc = Job_Creator(run_settings={'input':['AAA1234A1','AAA1234A2']}, config=preset_config, log=logger,sampleinfo=testdata)
  jc.snp_job()
  outfile = open(jc.get_sbatch(), 'r')
  count = 0
  for x in outfile.readlines():
    if "# SNP pair-wise distance" in x:
      count = count + 1
  assert count > 0

@patch('subprocess.Popen')
def test_project_job(subproc,testdata):
  #Sets up subprocess mocking
  process_mock = mock.Mock()
  attrs = {'communicate.return_value': ('output 123456789', 'error')}
  process_mock.configure_mock(**attrs)
  subproc.return_value = process_mock

  jc = Job_Creator( config=preset_config, log=logger, sampleinfo=testdata, run_settings={'pool':["AAA1234A1","AAA1234A2"], 'input':'/tmp/AAA1234'})
  jc.project_job()

def test_create_collection():
  pass

