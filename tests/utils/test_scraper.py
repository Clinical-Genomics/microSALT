#!/usr/bin/env python

import glob
import logging
import pytest


@pytest.mark.xfail(reason="Can no longer fetch from databases without authenticating")
def test_quast_scraping(scraper, testdata_prefix, caplog):
    scraper.scrape_quast(filename=f"{testdata_prefix}/quast_results.tsv")


@pytest.mark.xfail(reason="Can no longer fetch from databases without authenticating")
def test_blast_scraping(scraper, testdata_prefix, caplog):
    caplog.set_level(logging.DEBUG)
    scraper.scrape_blast(type="seq_type", file_list=[f"{testdata_prefix}/blast_single_loci.txt"])
    assert "candidate" in caplog.text

    caplog.clear()
    hits = scraper.scrape_blast(
        type="resistance", file_list=[f"{testdata_prefix}/blast_single_resistance.txt"]
    )
    genes = [h["gene"] for h in hits]

    assert "blaOXA-48" in genes
    assert "blaVIM-4" in genes


@pytest.mark.xfail(reason="Can no longer fetch from databases without authenticating")
def test_alignment_scraping(scraper, init_references, testdata_prefix):
    scraper.scrape_alignment(file_list=glob.glob(f"{testdata_prefix}/*.stats.*"))
