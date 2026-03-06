import json
import pathlib
import pytest

from microSALT.store.db_manipulator import DB_Manipulator
from microSALT import preset_config, logger


@pytest.fixture
def unpack_db_json():
    """Factory fixture: returns a callable that loads JSON files from tests/testdata/."""

    def _load(filename: str) -> list:
        path = pathlib.Path(__file__).parent / "testdata" / filename
        return json.loads(path.read_text())

    return _load


@pytest.fixture
def dbm(unpack_db_json):
    """DB_Manipulator populated with the standard set of test data."""
    dbm = DB_Manipulator(config=preset_config, log=logger)
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
