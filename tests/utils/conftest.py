from typing import Generator
from microSALT.store.db_manipulator import DB_Manipulator
import logging
from microSALT.config import MicroSALTConfig
import pathlib
import re
import pytest

from microSALT.utils.scraper import Scraper
from microSALT.utils.referencer import Referencer
from microSALT.utils.reporter import Reporter


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
def blast_scraper(config: MicroSALTConfig, logger: logging.Logger, testdata: list[dict], dbm: DB_Manipulator) -> Generator[Scraper, None, None]:
    """Scraper with full filesystem and DB context for BLAST scraping tests.

    Sets up:
    - ``{references}/staphylococcus_aureus/arcC.tfa`` so ``organism2reference``
      finds the organism and ``get_locilengths`` resolves allele lengths.
    - ``{profiles}/staphylococcus_aureus`` minimal TSV so ``alleles2st`` can
      resolve ST 130 from an arcC=3 hit.
    - ``{resistances}/combined.fsa`` so resistance allele lengths are resolved.
    - Project ``AAA1234`` and sample ``AAA1234A1`` in the database.

    The ``Scraper`` (and its internal ``DB_Manipulator``) is created *after*
    the profile file is written to disk so the profile table is picked up
    automatically during initialisation.
    """
    testdata_dir = pathlib.Path(__file__).parent.parent / "testdata"

    # --- 1. Build arcC.tfa from blast_single_loci.txt ---
    loci_blast = testdata_dir / "blast_single_loci.txt"
    loci_alleles: dict[str, int] = {}
    with open(loci_blast) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 12 and parts[1] != "N/A":
                loci_alleles[parts[3]] = int(parts[11])

    refs_dir = pathlib.Path(config.folders.references) / "staphylococcus_aureus"
    refs_dir.mkdir(parents=True, exist_ok=True)
    fasta_lines = []
    for allele_name, length in loci_alleles.items():
        fasta_lines.append(f">{allele_name}")
        fasta_lines.append("A" * length)
    (refs_dir / "arcC.tfa").write_text("\n".join(fasta_lines) + "\n")

    # --- 2. Build minimal staphylococcus_aureus MLST profile ---
    # Single-locus (arcC) scheme so allele_overabundance == 0 with one arcC hit.
    # ST 130 is mapped to arcC allele 3 (100 % identity hit in the blast file).
    profiles_dir = pathlib.Path(config.folders.profiles)
    profiles_dir.mkdir(parents=True, exist_ok=True)
    profile_path = profiles_dir / "staphylococcus_aureus"
    profile_path.write_text("ST\tarcC\n130\t3\n")

    # --- 3. Build combined.fsa from blast_single_resistance.txt ---
    res_blast = testdata_dir / "blast_single_resistance.txt"
    res_alleles: dict[str, int] = {}
    with open(res_blast) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 12 and parts[1] != "N/A":
                res_alleles[parts[3]] = int(parts[11])

    res_dir = pathlib.Path(config.folders.resistances)
    res_dir.mkdir(parents=True, exist_ok=True)
    res_fasta_lines = []
    for allele_name, length in res_alleles.items():
        res_fasta_lines.append(f">{allele_name}")
        res_fasta_lines.append("A" * length)
    (res_dir / "combined.fsa").write_text("\n".join(res_fasta_lines) + "\n")

    # --- 4. Seed DB: project AAA1234 is added by dbm; add sample AAA1234A1 ---
    dbm.add_rec(
        {"CG_ID_sample": "AAA1234A1", "CG_ID_project": "AAA1234"},
        "Samples",
    )

    # --- 5. Create Scraper AFTER filesystem is ready so its DB_Manipulator
    #        picks up the profile file and creates profile_staphylococcus_aureus ---
    scraper = Scraper(
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

    yield scraper

    # Cleanup filesystem (SQLite tables persist across runs, which is harmless)
    import shutil
    shutil.rmtree(str(refs_dir), ignore_errors=True)
    profile_path.unlink(missing_ok=True)
    (res_dir / "combined.fsa").unlink(missing_ok=True)


@pytest.fixture
def reporter(config, logger, unpack_db_json):
    return Reporter(
        log=logger,
        folders=config.folders,
        threshold=config.threshold,
        regex=config.regex,
        sampleinfo=unpack_db_json("sampleinfo_samples.json")[0],
        name="MIC1234A1",
        output="/tmp/MLST",
    )
