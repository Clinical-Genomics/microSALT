#!/usr/bin/env python

import json
import os
import pathlib
import pdb
import pytest

from distutils.sysconfig import get_python_lib

from microSALT import preset_config, logger
from microSALT.utils.scraper import Scraper


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
def scraper(testdata):
  scrape_obj = Scraper(config=preset_config, log=logger,sampleinfo=testdata[0]) 
  return scrape_obj

def test_quast_scraping(scraper, testdata_prefix, caplog):
  scraper.scrape_quast(filename="{}/quast_results.tsv".format(testdata_prefix))


def test_blast_scraping(scraper, testdata_prefix):
  scraper.scrape_blast(type='seq_type',file_list=["{}/blast_single_loci.txt".format(testdata_prefix)])
  scraper.scrape_blast(type='resistance',file_list=["{}/blast_single_resistance.txt".format(testdata_prefix)])

def test_alignment_scraping(scraper, testdata_prefix):
  scraper.scrape_alignment(file_list="{}/*.stats.*".format(testdata_prefix))
