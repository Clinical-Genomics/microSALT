import json
import logging
import pathlib
from importlib.resources import files as resource_files

import pytest

from microSALT import setup_logger
from microSALT.config import (
    Containers,
    Database,
    Folders,
    MicroSALTConfig,
    PasteurCredentials,
    PubMLSTCredentials,
    Regex,
    Singularity,
    SlurmHeader,
    Threshold,
)
from microSALT.store.database import initialize_database
from microSALT.store.db_manipulator import DB_Manipulator


@pytest.fixture(scope="session")
def config(tmp_path_factory: pytest.TempPathFactory) -> MicroSALTConfig:
    """Session-scoped config built from tmp_path_factory so all paths are isolated."""
    base = tmp_path_factory.mktemp("microsalt")

    results = base / "results"
    reports = base / "reports"
    seqdata = base / "projects"
    profiles = base / "references" / "ST_profiles"
    references = base / "references" / "ST_loci"
    resistances = base / "references" / "resistances"
    genomes = base / "references" / "genomes"
    credentials = base / "credentials"

    for d in (results, reports, seqdata, profiles, references, resistances, genomes, credentials):
        d.mkdir(parents=True, exist_ok=True)

    db_path = base / "microsalt.db"

    cfg = MicroSALTConfig(
        slurm_header=SlurmHeader(
            time="12:00:00",
            threads="8",
            qos="normal",
            job_prefix="MLST",
            project="production",
            type="core",
        ),
        regex=Regex(
            mail_recipient="username@suffix.com",
            file_pattern=r"\w{8,12}_\w{8,10}(?:-\d+)*_L\d_(?:R)*(\d{1}).fastq.gz",
            verified_organisms=[],
        ),
        folders=Folders(
            results=str(results),
            reports=str(reports),
            seqdata=str(seqdata),
            profiles=str(profiles),
            references=str(references),
            resistances=str(resistances),
            genomes=str(genomes),
            credentials=str(credentials),
        ),
        database=Database(
            SQLALCHEMY_DATABASE_URI=f"sqlite:///{db_path}",
            SQLALCHEMY_TRACK_MODIFICATIONS="False",
            DEBUG="True",
        ),
        threshold=Threshold(),
        pubmlst=PubMLSTCredentials(),
        pasteur=PasteurCredentials(),
        singularity=Singularity(),
        containers=Containers(),
    )
    cfg.folders.expec = str(resource_files("microSALT").joinpath("unique_references", "ExPEC.fsa"))
    cfg.config_path = str(base / "config.json")

    setup_logger(logging_level="INFO")
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
def dbm(config: MicroSALTConfig, logger: logging.Logger, unpack_db_json):
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
            "seqdata",
            "profiles",
            "references",
            "resistances",
            "genomes",
            "expec",
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
        "singularity": {"binary", "bind_paths", "trimmomatic_adapters"},
        "containers": {"skesa", "blast", "bwa", "samtools", "picard", "trimmomatic", "quast"},
    }
