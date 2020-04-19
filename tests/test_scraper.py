#!/usr/bin/env python

import json
import logging
import os
import pathlib
import pdb
import pytest

from distutils.sysconfig import get_python_lib

from microSALT import preset_config, logger
from microSALT.utils.scraper import Scraper
from microSALT.utils.referencer import Referencer

@pytest.fixture
def references():
  ref_obj = Referencer(config=preset_config, log=logger)
  fixer.identify_new(project_id,project=True)
  fixer.update_refs()

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

@pytest.fixture
def init_references(testdata):
  ref_obj = Referencer(config=preset_config, log=logger, sampleinfo=testdata)
  ref_obj.identify_new(testdata[0].get('CG_ID_project'),project=True)
  ref_obj.update_refs()

def test_quast_scraping(scraper, testdata_prefix, caplog):
  scraper.scrape_quast(filename="{}/quast_results.tsv".format(testdata_prefix))

def test_blast_scraping(scraper, testdata_prefix, caplog):
  caplog.set_level(logging.DEBUG)
  scraper.scrape_blast(type='seq_type',file_list=["{}/blast_single_loci.txt".format(testdata_prefix)])
  assert "candidate" in caplog.text
  caplog.clear()
  scraper.scrape_blast(type='resistance',file_list=["{}/blast_single_resistance.txt".format(testdata_prefix)])
  assert "candidate" in caplog.text

def test_alignment_scraping(scraper, init_references, testdata_prefix):
  scraper.scrape_alignment(file_list="{}/*.stats.*".format(testdata_prefix))
