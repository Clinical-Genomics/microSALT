#!/usr/bin/env python

import mock
import pdb
import pytest
import re

from unittest.mock import patch

from microSALT.store.db_manipulator import DB_Manipulator
from microSALT.utils.job_creator import Job_Creator
from microSALT import config, logger
from microSALT.cli import root

def fake_search(int):
  return "fake"

@pytest.fixture
def db_mani():
  db = DB_Manipulator(config, logger)
  return db

@patch('os.listdir')
@patch('os.stat')
@patch('gzip.open')
def test_verify_fastq(gopen, stat, listdir):
  listdir.return_value = ["ACC6438A3_HVMHWDSXX_L1_1.fastq.gz", "ACC6438A3_HVMHWDSXX_L1_2.fastq.gz", "ACC6438A3_HVMHWDSXX_L2_2.fastq.gz", "ACC6438A3_HVMHWDSXX_L2_2.fastq.gz"]
  stata = mock.MagicMock()
  stata.st_size = 2000
  stat.return_value = stata

  jc = Job_Creator('/tmp/', config, logger)
  t = jc.verify_fastq()
  assert len(t) > 0

@patch('re.search')
@patch('microSALT.utils.job_creator.glob.glob')
def test_blast_subset(glob_search, research):
  jc = Job_Creator('/tmp/', config, logger)
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
@patch('microSALT.store.lims_fetcher.LIMS_Fetcher')
@patch('microSALT.store.lims_fetcher.Lims.get_samples')
def test_create_snpsection(get_samples, LF, subproc):
  LF.data.return_value = {'CG_ID_project':"AAA1234",'CG_ID_sample':'AAA1234A1'}
  sample_mock = mock.MagicMock()
  sample_mock.project.id = "AAA1234"
  sample_mock.id = "AAA1234A3"
  sample_mock.name = "Trams"
  sample_mock.udf = {'Reference Genome Microbial':"NAN", 'customer':"NAN", 'Sequencing Analysis':"NAN"}
  get_samples.return_value = [sample_mock]

  #Sets up subprocess mocking
  process_mock = mock.Mock()
  attrs = {'communicate.return_value': ('output', 'error')}
  process_mock.configure_mock(**attrs)
  subproc.return_value = process_mock


  jc = Job_Creator(['a','b','c'], config, logger)
  jc.snp_job()
  outfile = open(jc.get_sbatch(), 'r')
  count = 0
  for x in outfile.readlines():
    if "# SNP pair-wise distance" in x:
      count = count + 1
  assert count > 0

def test_create_project():
  pass

def test_snp_job():
  pass

