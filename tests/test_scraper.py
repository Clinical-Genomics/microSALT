#!/usr/bin/env python

import glob
import json
import logging
import os
import pathlib
import pytest

from microSALT import preset_config, logger
from microSALT.utils.scraper import Scraper
from microSALT.utils.referencer import Referencer


@pytest.fixture
def testdata_prefix():
    return os.path.abspath(os.path.join(pathlib.Path(__file__).parent.parent, 'tests/testdata/'))


@pytest.fixture
def testdata():
    testdata = os.path.abspath(
        os.path.join(pathlib.Path(__file__).parent.parent, 'tests/testdata/sampleinfo_samples.json'))
    with open(testdata) as json_file:
        data = json.load(json_file)
    return data


@pytest.fixture
def scraper(testdata):
    scrape_obj = Scraper(config=preset_config, log=logger, sampleinfo=testdata[0])
    return scrape_obj


@pytest.fixture
def init_references(testdata):
    ref_obj = Referencer(config=preset_config, log=logger, sampleinfo=testdata)
    ref_obj.identify_new(testdata[0].get('CG_ID_project'), project=True)
    ref_obj.update_refs()


@pytest.mark.xfail(reason="Can no longer fetch from databases without authenticating")
def test_quast_scraping(scraper, testdata_prefix, caplog):
    scraper.scrape_quast(filename=f"{testdata_prefix}/quast_results.tsv")


@pytest.mark.xfail(reason="Can no longer fetch from databases without authenticating")
def test_blast_scraping(scraper, testdata_prefix, caplog):
    caplog.set_level(logging.DEBUG)
    scraper.scrape_blast(type='seq_type', file_list=[f"{testdata_prefix}/blast_single_loci.txt"])
    assert "candidate" in caplog.text

    caplog.clear()
    hits = scraper.scrape_blast(type='resistance', file_list=[f"{testdata_prefix}/blast_single_resistance.txt"])
    genes = [h["gene"] for h in hits]

    assert "blaOXA-48" in genes
    assert "blaVIM-4" in genes


@pytest.mark.xfail(reason="Can no longer fetch from databases without authenticating")
def test_alignment_scraping(scraper, init_references, testdata_prefix):
    scraper.scrape_alignment(file_list=glob.glob(f"{testdata_prefix}/*.stats.*"))
