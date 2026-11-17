#!/usr/bin/env python
"""Tests for microSALT.utils.scraper (Scraper)."""

import glob
import logging
import shutil
import uuid
from pathlib import Path
from typing import Generator
from unittest.mock import patch

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


def test_trimmomatic_scraping(scraper: Scraper, testdata_prefix: str) -> None:
    """Scraping a trimmomatic summary file should compute raw_reads as Input Read Pairs * 2."""
    scraper.scrape_trimmomatic(filename=f"{testdata_prefix}/trim_summary.txt")

    sample = scraper.db_pusher.get_sample_by_cg_id_sample(scraper.name)
    assert sample.raw_reads == 6508868 * 2


def _make_scraper(config, logger, sampleinfo, input_dir, qc_only=False) -> Scraper:
    return Scraper(
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
        sampleinfo=sampleinfo,
        input=str(input_dir),
        qc_only=qc_only,
    )


def test_scrape_sample_qc_only_skips_typing_for_unknown_organism(
    config: MicroSALTConfig, logger: logging.Logger, testdata: list[dict], tmp_path: Path
) -> None:
    """With qc_only, an organism without an MLST reference must not crash the scraping."""
    sampleinfo = {**testdata[0], "organism": "other"}
    input_dir = tmp_path / sampleinfo["CG_ID_sample"]
    input_dir.mkdir()
    qc_scraper = _make_scraper(config, logger, sampleinfo, input_dir, qc_only=True)

    with patch.object(Scraper, "scrape_blast", wraps=qc_scraper.scrape_blast) as blast:
        qc_scraper.scrape_sample()

    blast.assert_not_called()


def test_scrape_sample_qc_only_preserves_existing_typing_results(
    config: MicroSALTConfig,
    logger: logging.Logger,
    testdata: list[dict],
    dbm: DB_Manipulator,
    tmp_path: Path,
) -> None:
    """Standalone QC-only finish refreshes metadata while preserving typing results."""
    sampleinfo = dict(testdata[0])
    suffix = uuid.uuid4().hex[:8]
    sampleinfo["CG_ID_project"] = f"QC{suffix}"
    sampleinfo["CG_ID_sample"] = f"QC{suffix}A1"
    input_dir = tmp_path / sampleinfo["CG_ID_sample"]
    input_dir.mkdir()

    setup_scraper = _make_scraper(config, logger, sampleinfo, input_dir)
    setup_scraper.job_fallback.create_project(sampleinfo["CG_ID_project"])
    setup_scraper.job_fallback.create_sample()
    dbm.add_to_session(
        dbm.add_seq_type(
            CG_ID_sample=sampleinfo["CG_ID_sample"],
            loci="arcC",
            allele=3,
            contig_name="NODE_1",
            st_predictor=True,
        )
    )
    dbm.add_to_session(
        dbm.add_resistance(
            CG_ID_sample=sampleinfo["CG_ID_sample"],
            gene="blaOXA",
            instance="OXA",
            contig_name="NODE_1",
        )
    )
    dbm.add_to_session(
        dbm.add_expac(
            CG_ID_sample=sampleinfo["CG_ID_sample"],
            gene="virulence_gene",
            instance="ExPEC",
            contig_name="NODE_1",
        )
    )
    dbm.commit_session()
    dbm.update_sample(
        {"CG_ID_sample": sampleinfo["CG_ID_sample"]}, {"ST": 42, "pubmlst_ST": 17}
    )

    sampleinfo.update({
        "organism": "Updated organism",
        "Customer_ID_sample": "updated-customer-name",
        "reference": "UPDATED_REFERENCE",
        "date_libprep": "2026-11-01 12:00:00",
        "date_sequencing": "2026-11-02 12:00:00",
        "method_libprep": "updated-library-method",
        "method_sequencing": "updated-sequencing-method",
    })
    qc_scraper = _make_scraper(config, logger, sampleinfo, input_dir, qc_only=True)
    with (
        patch.object(qc_scraper.db_pusher, "delete_sample_results") as delete_results,
        patch.object(qc_scraper, "scrape_alignment"),
        patch.object(qc_scraper, "scrape_trimmomatic"),
        patch.object(qc_scraper, "scrape_quast"),
        patch.object(qc_scraper, "scrape_reference"),
    ):
        qc_scraper.scrape_sample()

    delete_results.assert_not_called()
    saved_sample = dbm.get_sample_by_cg_id_sample(sampleinfo["CG_ID_sample"])
    assert saved_sample.organism == sampleinfo["organism"]
    assert saved_sample.Customer_ID_sample == sampleinfo["Customer_ID_sample"]
    assert saved_sample.reference_genome == sampleinfo["reference"]
    assert saved_sample.date_libprep.isoformat(sep=" ") == sampleinfo["date_libprep"]
    assert saved_sample.date_sequencing.isoformat(sep=" ") == sampleinfo["date_sequencing"]
    assert saved_sample.method_libprep == sampleinfo["method_libprep"]
    assert saved_sample.method_sequencing == sampleinfo["method_sequencing"]
    assert (saved_sample.ST, saved_sample.pubmlst_ST) == (42, 17)
    assert [row.loci for row in saved_sample.seq_types] == ["arcC"]
    assert [row.gene for row in saved_sample.resistances] == ["blaOXA"]
    assert [row.gene for row in saved_sample.expacs] == ["virulence_gene"]


def test_scrape_project_qc_only_does_not_delete_project_samples(
    config: MicroSALTConfig,
    logger: logging.Logger,
    testdata: list[dict],
    tmp_path: Path,
) -> None:
    """QC-only project runs process current samples without deleting other project members."""
    sampleinfo = [dict(entry) for entry in testdata[:2]]
    input_dir = tmp_path / "project"
    input_dir.mkdir()
    for entry in sampleinfo:
        (input_dir / entry["CG_ID_sample"]).mkdir()

    project_scraper = _make_scraper(config, logger, sampleinfo, input_dir, qc_only=True)
    with (
        patch.object(project_scraper.db_pusher, "delete_project") as delete_project,
        patch.object(project_scraper.job_fallback, "create_project") as create_project,
        patch.object(Scraper, "scrape_sample") as scrape_sample,
    ):
        project_scraper.scrape_project()

    delete_project.assert_not_called()
    create_project.assert_called_once_with(project_scraper.name)
    assert scrape_sample.call_count == len(sampleinfo)


def test_scrape_sample_default_runs_typing(
    config: MicroSALTConfig, logger: logging.Logger, testdata: list[dict], tmp_path: Path
) -> None:
    """Without qc_only, the BLAST typing/resistance scraping still runs."""
    input_dir = tmp_path / testdata[0]["CG_ID_sample"]
    input_dir.mkdir()
    default_scraper = _make_scraper(config, logger, testdata[0], input_dir)

    with (
        patch.object(
            default_scraper.referencer, "organism2reference", return_value="staphylococcus_aureus"
        ),
        patch.object(Scraper, "scrape_blast") as blast,
    ):
        default_scraper.scrape_sample()

    called_types = [c.kwargs.get("type") for c in blast.call_args_list]
    assert "seq_type" in called_types
    assert "resistance" in called_types


def test_scrape_sample_without_mlst_scheme_keeps_other_analysis(
    config: MicroSALTConfig, logger: logging.Logger, testdata: list[dict], tmp_path: Path
) -> None:
    sampleinfo = {**testdata[0], "organism": "Unknown organism"}
    input_dir = tmp_path / sampleinfo["CG_ID_sample"]
    input_dir.mkdir()
    res_scraper = _make_scraper(config, logger, sampleinfo, input_dir)
    res_scraper.job_fallback.create_project(sampleinfo["CG_ID_project"])
    res_scraper.job_fallback.create_sample()
    res_scraper.db_pusher.update_sample(
        {"CG_ID_sample": res_scraper.name}, {"ST": 42, "pubmlst_ST": 42}
    )
    with patch.object(res_scraper, "scrape_blast", wraps=res_scraper.scrape_blast) as blast:
        res_scraper.scrape_sample()
    assert [call.kwargs["type"] for call in blast.call_args_list] == ["resistance"]
    sample = res_scraper.db_pusher.get_sample_by_cg_id_sample(res_scraper.name)
    assert sample.ST == sample.pubmlst_ST == -1
