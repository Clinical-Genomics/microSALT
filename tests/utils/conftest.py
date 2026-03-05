import pathlib
import pytest

from microSALT.utils.scraper import Scraper
from microSALT.utils.referencer import Referencer
from microSALT.utils.reporter import Reporter


@pytest.fixture
def testdata_prefix():
    return str(pathlib.Path(__file__).parent.parent / "testdata")


@pytest.fixture
def testdata(unpack_db_json):
    return unpack_db_json("sampleinfo_samples.json")


@pytest.fixture
def scraper(config, logger, testdata):
    return Scraper(config=config, log=logger, sampleinfo=testdata[0])


@pytest.fixture
def init_references(config, logger, testdata):
    ref_obj = Referencer(config=config, log=logger, sampleinfo=testdata)
    ref_obj.identify_new(testdata[0].get("CG_ID_project"), project=True)
    ref_obj.update_refs()


@pytest.fixture
def reporter(config, logger, unpack_db_json):
    return Reporter(
        config=config,
        log=logger,
        sampleinfo=unpack_db_json("sampleinfo_samples.json")[0],
        name="MIC1234A1",
        output="/tmp/MLST",
    )
