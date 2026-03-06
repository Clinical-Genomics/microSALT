#!/usr/bin/env python
"""Tests for microSALT.utils.scraper (Scraper)."""

import glob
import logging

import pytest

from microSALT.utils.scraper import Scraper


def test_quast_scraping(
    scraper: Scraper, testdata_prefix: str, caplog: pytest.LogCaptureFixture
) -> None:
    """Scraping a QUAST results file should not raise."""
    scraper.scrape_quast(filename=f"{testdata_prefix}/quast_results.tsv")


def test_blast_scraping(
    blast_scraper: Scraper, testdata_prefix: str, caplog: pytest.LogCaptureFixture
) -> None:
    """BLAST scraping should find sequence-type candidates and resistance genes."""
    caplog.set_level(logging.DEBUG)

    blast_scraper.scrape_blast(
        type="seq_type", file_list=[f"{testdata_prefix}/blast_single_loci.txt"]
    )
    assert "candidate" in caplog.text

    caplog.clear()
    hits = blast_scraper.scrape_blast(
        type="resistance",
        file_list=[f"{testdata_prefix}/blast_single_resistance.txt"],
    )
    genes = [h["gene"] for h in hits]
    assert "blaOXA-48" in genes
    assert "blaVIM-4" in genes


def test_alignment_scraping(scraper: Scraper, testdata_prefix: str) -> None:
    """Scraping alignment stats files should not raise.

    init_references is NOT needed: scrape_alignment only reads local .stats.* files.
    """
    scraper.scrape_alignment(file_list=glob.glob(f"{testdata_prefix}/*.stats.*"))
