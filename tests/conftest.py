import json
import logging
import os
import pathlib
import pytest

from microSALT.config import load_config
from microSALT import setup_logger
from microSALT.store.database import initialize_database
from microSALT.store.db_manipulator import DB_Manipulator


def _config_path() -> str:
    env = os.environ.get("MICROSALT_CONFIG")
    if env:
        return env
    default = pathlib.Path(__file__).parent.parent / "configExample.json"
    return str(default)


@pytest.fixture(scope="session")
def config():
    cfg = load_config(_config_path())
    setup_logger(logging_level="INFO", log_file=cfg.folders.log_file)
    initialize_database(cfg.database.SQLALCHEMY_DATABASE_URI)
    return cfg


@pytest.fixture(scope="session")
def logger():
    return logging.getLogger("main_logger")


@pytest.fixture
def unpack_db_json():
    """Factory fixture: returns a callable that loads JSON files from tests/testdata/."""

    def _load(filename: str) -> list:
        path = pathlib.Path(__file__).parent / "testdata" / filename
        return json.loads(path.read_text())

    return _load


@pytest.fixture
def dbm(config, logger, unpack_db_json):
    """DB_Manipulator populated with the standard set of test data."""
    dbm = DB_Manipulator(log=logger, folders=config.folders, threshold=config.threshold)
    dbm.create_tables()

    for entry in unpack_db_json("sampleinfo_projects.json"):
        dbm.add_rec(entry, "Projects")
    for entry in unpack_db_json("sampleinfo_mlst.json"):
        dbm.add_rec(entry, "Seq_types")
    for entry in unpack_db_json("sampleinfo_resistance.json"):
        dbm.add_rec(entry, "Resistances")
    for entry in unpack_db_json("sampleinfo_expec.json"):
        dbm.add_rec(entry, "Expacs")
    for entry in unpack_db_json("sampleinfo_reports.json"):
        dbm.add_rec(entry, "Reports")
    return dbm


@pytest.fixture
def exp_config():
    """Expected configuration structure for config validation tests."""
    return {
        "slurm_header": {"time", "threads", "qos", "job_prefix", "project", "type"},
        "regex": {"file_pattern", "mail_recipient", "verified_organisms"},
        "folders": {
            "results",
            "reports",
            "log_file",
            "seqdata",
            "profiles",
            "references",
            "resistances",
            "genomes",
            "expec",
            "adapters",
            "credentials",
        },
        "threshold": {
            "mlst_id",
            "mlst_novel_id",
            "mlst_span",
            "motif_id",
            "motif_span",
            "total_reads_warn",
            "total_reads_fail",
            "NTC_total_reads_warn",
            "NTC_total_reads_fail",
            "mapped_rate_warn",
            "mapped_rate_fail",
            "duplication_rate_warn",
            "duplication_rate_fail",
            "insert_size_warn",
            "insert_size_fail",
            "average_coverage_warn",
            "average_coverage_fail",
            "bp_10x_warn",
            "bp_10x_fail",
            "bp_30x_warn",
            "bp_50x_warn",
            "bp_100x_warn",
        },
        "database": {"SQLALCHEMY_DATABASE_URI", "SQLALCHEMY_TRACK_MODIFICATIONS", "DEBUG"},
        "pubmlst": {"client_id", "client_secret"},
        "pasteur": {"client_id", "client_secret"},
    }
