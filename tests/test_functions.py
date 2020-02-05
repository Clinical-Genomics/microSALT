#!/usr/bin/env python

import mock
import pdb
import pytest

from unittest.mock import patch

from microSALT.store.db_manipulator import DB_Manipulator
from microSALT.utils.job_creator import Job_Creator
from microSALT import config, logger
from microSALT.cli import root


@pytest.fixture
def db_mani():
  db = DB_Manipulator(config, logger)
  return db

#@pytest.fixture
#def job_creator():
#  jc = Job_Creator(config, logger)
#  return jc


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
def test_blast_subset(research):
  jc = Job_Creator('/tmp/', config, logger)
  researcha = mock.MagicMock()
  #researcha.group(1) = 'test'
  #reserch.return_value = researcha
  #pdb.set_trace()
#
  #jc.blast_subset('mlst', '/tmp/*')
  #jc.blast_subset('other', '/tmp/*')

def test_create_snpsection():
  pass

def test_create_project():
  pass

def test_snp_job():
  pass
