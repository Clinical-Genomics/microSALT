import json
import pytest
import os
import shutil
import tempfile

from flask import Flask
from pathlib import Path
from typing import Generator

from microSALT.cli import initialize_logger
from microSALT.server.app import initialize_app, get_app
from microSALT.store.database import initialize_database, drop_all_tables
from microSALT.store.db_manipulator import DB_Manipulator
from microSALT.utils.reporter import Reporter


def unpack_db_json(fixture_dir, filename):
    testdata = Path(fixture_dir, filename)
    with open(testdata) as json_file:
        data = json.load(json_file)
    return data


@pytest.fixture(scope="session")
def app(config) -> Generator[Flask, None, None]:
    initialize_app(config)
    yield get_app()


@pytest.fixture(scope="session")
def fixture_dir() -> Path:
    return Path("tests", "testdata")


@pytest.fixture(scope="session")
def config_file(fixture_dir: Path) -> Path:
    return Path(fixture_dir, "config.json")


@pytest.fixture(scope="session")
def temp_dirs() -> Generator[dict, None, None]:
    temp_dirs = {
        "results": tempfile.mkdtemp(),
        "reports": tempfile.mkdtemp(),
        "log_file": tempfile.mkstemp()[1],
        "seqdata": tempfile.mkdtemp(),
        "profiles": tempfile.mkdtemp(),
        "references": tempfile.mkdtemp(),
        "resistances": tempfile.mkdtemp(),
        "genomes": tempfile.mkdtemp(),
        "pubmlst_credentials": tempfile.mkdtemp(),
        "adapters": tempfile.mkdtemp(),
        "expec": tempfile.mkdtemp(),
    }
    yield temp_dirs
    for path in temp_dirs.values():
        if os.path.isdir(path):
            shutil.rmtree(path)
        elif os.path.isfile(path):
            os.remove(path)


@pytest.fixture(scope="session")
def config(temp_dirs: dict) -> dict:
    config = {
        "_comment": "SBATCH configuration",
        "slurm_header": {
            "time": "12:00:00",
            "threads": "8",
            "qos": "normal",
            "job_prefix": "MLST",
            "project": "production",
            "type": "core",
        },
        "regex": {
            "mail_recipient": "username@suffix.com",
            "_comment": "File finding patterns. Only single capture group accepted (for reverse/forward identifier)",
            "file_pattern": "\\w{8,12}_\\w{8,10}(?:-\\d+)*_L\\d_(?:R)*(\\d{1}).fastq.gz",
            "_comment": "Organisms recognized enough to be considered stable",
            "verified_organisms": [],
        },
        "_comment": "Folders",
        "folders": {
            "results": temp_dirs["results"],
            "reports": temp_dirs["reports"],
            "log_file": temp_dirs["log_file"],
            "seqdata": temp_dirs["seqdata"],
            "profiles": temp_dirs["profiles"],
            "references": temp_dirs["references"],
            "resistances": temp_dirs["resistances"],
            "genomes": temp_dirs["genomes"],
            "pubmlst_credentials": temp_dirs["pubmlst_credentials"],
            "adapters": temp_dirs["adapters"],
            "expec": temp_dirs["expec"],
        },
        "_comment": "Database/Flask configuration",
        "database": {
            "SQLALCHEMY_DATABASE_URI": "sqlite:////tmp/microSALT.db",
            "SQLALCHEMY_TRACK_MODIFICATIONS": "False",
            "DEBUG": "True",
        },
        "_comment": "Thresholds for Displayed results",
        "threshold": {
            "_comment": "Typing thresholds",
            "mlst_id": 100,
            "mlst_novel_id": 99.5,
            "mlst_span": 90,
            "motif_id": 97,
            "motif_span": 90,
            "_comment": "Quality Control thresholds",
            "total_reads_warn": 75,
            "total_reads_fail": 70,
            "NTC_total_reads_warn": 10,
            "NTC_total_reads_fail": 20,
            "mapped_rate_warn": 50,
            "mapped_rate_fail": 30,
            "duplication_rate_warn": 20,
            "duplication_rate_fail": 80,
            "insert_size_warn": 140,
            "insert_size_fail": 100,
            "average_coverage_warn": 100,
            "average_coverage_fail": 10,
            "bp_10x_warn": 85,
            "bp_10x_fail": 75,
            "bp_30x_warn": 70,
            "bp_50x_warn": 50,
            "bp_100x_warn": 20,
        },
        "_comment": "Genologics temporary configuration file",
        "genologics": {
            "baseuri": "https://lims.facility.se/",
            "username": "limsuser",
            "password": "mypassword",
        },
        "_comment": "PubMLST credentials",
        "pubmlst": {"client_id": "", "client_secret": ""},
    }
    return config


@pytest.fixture(scope="session")
def logger(config: dict):
    """Initializes the logger"""
    yield initialize_logger(config=config)


@pytest.fixture(scope="session")
def sampleinfo_samples(fixture_dir) -> dict:
    return unpack_db_json(fixture_dir, "sampleinfo_samples.json")


@pytest.fixture(scope="session")
def sampleinfo_sample(sampleinfo_samples) -> dict:
    return sampleinfo_samples[0]


@pytest.fixture(scope="session")
def sampleinfo_projects(fixture_dir) -> dict:
    return unpack_db_json(fixture_dir, "sampleinfo_projects.json")


@pytest.fixture(scope="session")
def sampleinfo_mlst(fixture_dir) -> dict:
    return unpack_db_json(fixture_dir, "sampleinfo_mlst.json")


@pytest.fixture(scope="session")
def sampleinfo_resistance(fixture_dir) -> dict:
    return unpack_db_json(fixture_dir, "sampleinfo_resistance.json")


@pytest.fixture(scope="session")
def sampleinfo_expec(fixture_dir) -> dict:
    return unpack_db_json(fixture_dir, "sampleinfo_expec.json")


@pytest.fixture(scope="session")
def sampleinfo_reports(fixture_dir) -> dict:
    return unpack_db_json(fixture_dir, "sampleinfo_reports.json")


@pytest.fixture(scope="session")
def sample_info_files(
    sampleinfo_projects,
    sampleinfo_mlst,
    sampleinfo_resistance,
    sampleinfo_expec,
    sampleinfo_reports,
) -> dict:
    return {
        "Projects": sampleinfo_projects,
        "Seq_types": sampleinfo_mlst,
        "Resistances": sampleinfo_resistance,
        "Expacs": sampleinfo_expec,
        "Reports": sampleinfo_reports,
    }


@pytest.fixture(scope="session")
def dbm(config, sample_info_files, logger) -> Generator[DB_Manipulator, None, None]:
    initialize_database(config)
    _dbm = DB_Manipulator(config=config, log=logger)

    for table, entry in sample_info_files.items():
        for rec in entry:
            _dbm.add_rec(rec, table)

    yield _dbm
    drop_all_tables()


@pytest.fixture
def exp_config() -> dict:
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
            "pubmlst_credentials",
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
        "genologics": {"baseuri", "username", "password"},
        "pubmlst": {"client_id", "client_secret"},
        "dry": True,
    }


@pytest.fixture(scope="function")
def reporter(
    app: Flask, config, dbm, logger, sampleinfo_sample, temp_dirs
) -> Generator[Reporter, None, None]:
    yield Reporter(
        app=app,
        config=config,
        dbm=dbm,
        log=logger,
        output=temp_dirs["reports"],
        sampleinfo=sampleinfo_sample,
    )
