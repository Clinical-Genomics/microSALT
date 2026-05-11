#!/usr/bin/env python
"""Tests for microSALT.utils.scraper (Scraper)."""

import glob
import logging
import shutil
from pathlib import Path
from typing import Generator

import pytest

from microSALT.config import MicroSALTConfig
from microSALT.store.db_manipulator import DB_Manipulator
from microSALT.utils.scraper import Scraper


class BlastScraperContext:
    """Filesystem and DB context for BLAST scraping tests.

    Attributes:
        refs_dir: Directory containing the organism-specific '.tfa' FASTA files.
        profile_path: Path to the 'staphylococcus_aureus' MLST profile file.
        combined_fsa: Path to the resistance 'combined.fsa' file.
        scraper: The configured 'Scraper' instance.
    """

    refs_dir: Path
    profile_path: Path
    combined_fsa: Path
    scraper: Scraper

    def setup_loci_fasta(self, testdata_dir: Path, references_dir: Path) -> None:
        """Write ``{references}/staphylococcus_aureus/arcC.tfa`` from blast_single_loci.txt."""
        alleles: dict[str, int] = {}
        with open(testdata_dir / "blast_single_loci.txt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 12 and parts[1] != "N/A":
                    alleles[parts[3]] = int(parts[11])

        self.refs_dir = references_dir / "staphylococcus_aureus"
        self.refs_dir.mkdir(parents=True, exist_ok=True)
        (self.refs_dir / "arcC.tfa").write_text(
            "\n".join(
                line for name, length in alleles.items() for line in (f">{name}", "A" * length)
            )
            + "\n"
        )

    def setup_profile(self, profiles_dir: Path) -> None:
        """Write a minimal ``staphylococcus_aureus`` MLST profile (ST 130 → arcC allele 3)."""
        profiles_dir.mkdir(parents=True, exist_ok=True)
        self.profile_path = profiles_dir / "staphylococcus_aureus"
        self.profile_path.write_text("ST\tarcC\n130\t3\n")

    def setup_resistance_fasta(self, testdata_dir: Path, resistances_dir: Path) -> None:
        """Write ``{resistances}/combined.fsa`` from blast_single_resistance.txt."""
        alleles: dict[str, int] = {}
        with open(testdata_dir / "blast_single_resistance.txt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 12 and parts[1] != "N/A":
                    alleles[parts[3]] = int(parts[11])

        resistances_dir.mkdir(parents=True, exist_ok=True)
        self.combined_fsa = resistances_dir / "combined.fsa"
        self.combined_fsa.write_text(
            "\n".join(
                line for name, length in alleles.items() for line in (f">{name}", "A" * length)
            )
            + "\n"
        )

    def teardown(self) -> None:
        """Remove generated filesystem artifacts created by the setup methods."""
        shutil.rmtree(self.refs_dir, ignore_errors=True)
        self.profile_path.unlink(missing_ok=True)
        self.combined_fsa.unlink(missing_ok=True)

    def _get_blast_scaper_context(
        self,
        config: MicroSALTConfig,
        logger: logging.Logger,
        testdata: list[dict],
        dbm: DB_Manipulator,
    ) -> Generator["BlastScraperContext", None, None]:
        """Return a setup BlastScraperContext instance configured for BLAST scraping tests."""
        self.setup_loci_fasta(
            Path(__file__).parent.parent / "testdata", Path(config.folders.references)
        )
        self.setup_profile(Path(config.folders.profiles))
        self.setup_resistance_fasta(
            Path(__file__).parent.parent / "testdata", Path(config.folders.resistances)
        )

        dbm.add_to_session(dbm.add_sample(CG_ID_sample="AAA1234A1", CG_ID_project="AAA1234"))
        dbm.commit_session()

        self.scraper = Scraper(
            log=logger,
            folders=config.folders,
            threshold=config.threshold,
            slurm_header=config.slurm_header,
            regex=config.regex,
            dry=False,
            config_path=config.config_path,
            pubmlst=config.pubmlst,
            pasteur=config.pasteur,
            singularity=config.singularity,
            containers=config.containers,
            sampleinfo=testdata[0],
        )

        yield self

        self.teardown()


@pytest.fixture
def blast_scraper_context(
    config: MicroSALTConfig,
    logger: logging.Logger,
    testdata: list[dict],
    dbm: DB_Manipulator,
) -> Generator[BlastScraperContext, None, None]:
    """Yield a fully set-up :class:`BlastScraperContext`, then tear it down."""
    testdata_dir = Path(__file__).parent.parent / "testdata"

    ctx = BlastScraperContext()
    ctx.setup_loci_fasta(testdata_dir, Path(config.folders.references))
    ctx.setup_profile(Path(config.folders.profiles))
    ctx.setup_resistance_fasta(testdata_dir, Path(config.folders.resistances))

    dbm.add_to_session(dbm.add_sample(CG_ID_sample="AAA1234A1", CG_ID_project="AAA1234"))
    dbm.commit_session()

    ctx.scraper = Scraper(
        log=logger,
        folders=config.folders,
        threshold=config.threshold,
        slurm_header=config.slurm_header,
        regex=config.regex,
        dry=False,
        config_path=config.config_path,
        pubmlst=config.pubmlst,
        pasteur=config.pasteur,
        singularity=config.singularity,
        containers=config.containers,
        sampleinfo=testdata[0],
    )

    yield ctx

    ctx.teardown()


def test_quast_scraping(
    scraper: Scraper, testdata_prefix: str, caplog: pytest.LogCaptureFixture
) -> None:
    """Scraping a QUAST results file should not raise."""
    scraper.scrape_quast(filename=f"{testdata_prefix}/quast_results.tsv")


def test_blast_scraping(
    blast_scraper_context: BlastScraperContext,
    testdata_prefix: str,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """BLAST scraping should find sequence-type candidates and resistance genes."""
    caplog.set_level(logging.DEBUG)

    blast_scraper_context.scraper.scrape_blast(
        type="seq_type", file_list=[f"{testdata_prefix}/blast_single_loci.txt"]
    )
    assert "candidate" in caplog.text

    caplog.clear()
    hits = blast_scraper_context.scraper.scrape_blast(
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
