#!/usr/bin/env python

import glob
import logging
from unittest.mock import patch

from microSALT.utils.scraper import Scraper


class _MockLocilengths:
    """Fake locilengths mapping used in tests.

    Any ``x.startswith(prefix)`` check on keys() always returns True, and any
    key lookup returns 1000 (a realistic locus length for span computation).
    """

    class _AnyKey(str):
        def startswith(self, *args, **kwargs):
            return True

    _SENTINEL = _AnyKey(">mock_locus_1")

    def keys(self):
        return [self._SENTINEL]

    def __getitem__(self, key):
        return 1000


def test_quast_scraping(scraper, testdata_prefix, caplog):
    scraper.scrape_quast(filename=f"{testdata_prefix}/quast_results.tsv")


def test_blast_scraping(scraper, testdata_prefix, caplog):
    caplog.set_level(logging.DEBUG)
    with patch.object(Scraper, "get_locilengths", return_value=_MockLocilengths()):
        scraper.scrape_blast(
            type="seq_type", file_list=[f"{testdata_prefix}/blast_single_loci.txt"]
        )
    assert "candidate" in caplog.text

    caplog.clear()
    with patch.object(Scraper, "get_locilengths", return_value=_MockLocilengths()):
        hits = scraper.scrape_blast(
            type="resistance",
            file_list=[f"{testdata_prefix}/blast_single_resistance.txt"],
        )
    genes = [h["gene"] for h in hits]
    assert "blaOXA-48" in genes
    assert "blaVIM-4" in genes


def test_alignment_scraping(scraper, testdata_prefix):
    # init_references is NOT needed: scrape_alignment only reads local .stats.* files.
    scraper.scrape_alignment(file_list=glob.glob(f"{testdata_prefix}/*.stats.*"))
