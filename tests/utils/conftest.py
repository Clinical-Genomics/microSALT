import logging
from pathlib import Path

import pytest

from microSALT.config import MicroSALTConfig
from microSALT.utils.referencer import Referencer
from microSALT.utils.reporter import Reporter
from microSALT.utils.scraper import Scraper


@pytest.fixture
def testdata_prefix():
    return str(Path(__file__).parent.parent / "testdata")


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
        singularity=config.singularity,
        containers=config.containers,
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
        singularity=config.singularity,
        containers=config.containers,
        sampleinfo=testdata,
    )
    ref_obj.identify_new()
    ref_obj.update_refs()


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
