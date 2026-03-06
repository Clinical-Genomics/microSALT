import logging
import pathlib
import shutil
from typing import Generator

import pytest

from microSALT.config import MicroSALTConfig
from microSALT.store.db_manipulator import DB_Manipulator
from microSALT.utils.referencer import Referencer
from microSALT.utils.reporter import Reporter
from microSALT.utils.scraper import Scraper


class BlastScraperContext:
    """Filesystem and DB context for BLAST scraping tests.

    Attributes:
        refs_dir: Directory containing the organism-specific '.tfa' FASTA files.
        profile_path: Path to the 'staphylococcus_aureus' MLST profile file.
        combined_fsa: Path to the resistance 'combined.fsa' file.
        scraper: The configured 'Scraper' instance.
    """

    refs_dir: pathlib.Path
    profile_path: pathlib.Path
    combined_fsa: pathlib.Path
    scraper: Scraper

    def setup_loci_fasta(self, testdata_dir: pathlib.Path, references_dir: pathlib.Path) -> None:
        """Write ``{references}/staphylococcus_aureus/arcC.tfa`` from blast_single_loci.txt.

        Allele lengths are read from the BLAST output so
        ``get_locilengths`` resolves correctly during scraping.

        Args:
            testdata_dir: Directory containing ``blast_single_loci.txt``.
            references_dir: Value of ``config.folders.references`` as a Path.
        """
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
        (self.refs_dir / "arcC.tfa").write_text("\n".join(
            line for name, length in alleles.items()
            for line in (f">{name}", "A" * length)
        ) + "\n")

    def setup_profile(self, profiles_dir: pathlib.Path) -> None:
        """Write a minimal ``staphylococcus_aureus`` MLST profile (ST 130 → arcC allele 3).

        A single-locus scheme keeps ``allele_overabundance`` at zero so
        ``alleles2st`` resolves ST 130 from the arcC=3 hit in the BLAST file.

        Args:
            profiles_dir: Value of ``config.folders.profiles`` as a Path.
        """
        profiles_dir.mkdir(parents=True, exist_ok=True)
        self.profile_path = profiles_dir / "staphylococcus_aureus"
        self.profile_path.write_text("ST\tarcC\n130\t3\n")

    def setup_resistance_fasta(self, testdata_dir: pathlib.Path, resistances_dir: pathlib.Path) -> None:
        """Write ``{resistances}/combined.fsa`` from blast_single_resistance.txt.

        Allele lengths are read from the BLAST output so resistance allele
        lengths are resolved correctly during scraping.

        Args:
            testdata_dir: Directory containing ``blast_single_resistance.txt``.
            resistances_dir: Value of ``config.folders.resistances`` as a Path.
        """
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
        self.combined_fsa.write_text("\n".join(
            line for name, length in alleles.items()
            for line in (f">{name}", "A" * length)
        ) + "\n")

    def teardown(self) -> None:
        """Remove generated filesystem artifacts created by the setup methods."""
        shutil.rmtree(self.refs_dir, ignore_errors=True)
        self.profile_path.unlink(missing_ok=True)
        self.combined_fsa.unlink(missing_ok=True)


@pytest.fixture
def testdata_prefix():
    return str(pathlib.Path(__file__).parent.parent / "testdata")


@pytest.fixture
def testdata(unpack_db_json) -> list[dict]:
    return unpack_db_json("sampleinfo_samples.json")


@pytest.fixture
def scraper(config: MicroSALTConfig, logger: logging.Logger, testdata: list[dict]) -> Scraper:
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
        sampleinfo=testdata[0],
    )


@pytest.fixture
def init_references(config: MicroSALTConfig, logger: logging.Logger, testdata: list[dict]) -> None:
    ref_obj = Referencer(
        log=logger,
        folders=config.folders,
        threshold=config.threshold,
        pubmlst=config.pubmlst,
        pasteur=config.pasteur,
        sampleinfo=testdata,
    )
    ref_obj.identify_new(testdata[0].get("CG_ID_project"), project=True)
    ref_obj.update_refs()


@pytest.fixture
def blast_scraper(
    config: MicroSALTConfig,
    logger: logging.Logger,
    testdata: list[dict],
    dbm: DB_Manipulator,
) -> Generator[BlastScraperContext, None, None]:
    """Yield a :class:`BlastScraperContext` with filesystem and DB fully set up.

    The :class:`Scraper` is created *after* the profile file is written to disk
    so its internal ``DB_Manipulator`` picks up the profile table on first use.
    """
    testdata_dir = pathlib.Path(__file__).parent.parent / "testdata"

    ctx = BlastScraperContext()
    ctx.setup_loci_fasta(testdata_dir, pathlib.Path(config.folders.references))
    ctx.setup_profile(pathlib.Path(config.folders.profiles))
    ctx.setup_resistance_fasta(testdata_dir, pathlib.Path(config.folders.resistances))

    dbm.add_rec({"CG_ID_sample": "AAA1234A1", "CG_ID_project": "AAA1234"}, "Samples")

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
        sampleinfo=testdata[0],
    )

    yield ctx

    ctx.teardown()


@pytest.fixture
def reporter(config: MicroSALTConfig, logger: logging.Logger, unpack_db_json):
    return Reporter(
        log=logger,
        folders=config.folders,
        threshold=config.threshold,
        regex=config.regex,
        sampleinfo=unpack_db_json("sampleinfo_samples.json")[0],
        name="MIC1234A1",
        output=config.folders.reports,
    )
