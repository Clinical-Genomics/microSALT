#!/usr/bin/env python
"""Tests for microSALT.utils.scraper (Scraper)."""

import glob
import logging
from unittest.mock import patch

import pytest

from microSALT.utils.scraper import Scraper


class _MockLocilengths:
    """Fake locilengths mapping used in tests.

    Any ``x.startswith(prefix)`` check on keys() always returns True, and any
    key lookup returns 1000 (a realistic locus length for span computation).
    """

    class _AnyKey(str):
        def startswith(  # type: ignore[override]
            self,
            prefix: str | tuple[str, ...],
            start: int | None = None,
            end: int | None = None,
        ) -> bool:
            return True

    _SENTINEL: _AnyKey
    _SENTINEL = _AnyKey(">mock_locus_1")

    def keys(self) -> list[_AnyKey]:
        """Return the single sentinel key."""
        return [self._SENTINEL]

    def __getitem__(self, key: object) -> int:
        """Return a fixed realistic locus length for any key."""
        return 1000


def test_quast_scraping(
    scraper: Scraper, testdata_prefix: str, caplog: pytest.LogCaptureFixture
) -> None:
    """Scraping a QUAST results file should not raise."""
    scraper.scrape_quast(filename=f"{testdata_prefix}/quast_results.tsv")


def test_blast_scraping(
    scraper: Scraper, testdata_prefix: str, caplog: pytest.LogCaptureFixture
) -> None:
    """BLAST scraping should find sequence-type candidates and resistance genes."""
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


def test_alignment_scraping(scraper: Scraper, testdata_prefix: str) -> None:
    """Scraping alignment stats files should not raise.

    init_references is NOT needed: scrape_alignment only reads local .stats.* files.
    """
    scraper.scrape_alignment(file_list=glob.glob(f"{testdata_prefix}/*.stats.*"))
