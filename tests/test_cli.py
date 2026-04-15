"""Tests for the CLI entry point — verifies config loading and propagation."""

import json
from pathlib import Path
from unittest.mock import patch

import pytest
from click.testing import CliRunner
from sqlalchemy import create_engine, inspect

from microSALT.cli import root
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
from microSALT.store.orm_models import Base


@pytest.fixture
def config_file(config: MicroSALTConfig, tmp_path: Path) -> Path:
    """Write the session config to a temporary JSON file for CLI consumption."""
    path = tmp_path / "config.json"
    path.write_text(config.model_dump_json())
    return path


@pytest.fixture
def setup_config(tmp_path: Path) -> MicroSALTConfig:
    """Fresh MicroSALTConfig whose directories do not yet exist — used to test `setup`."""
    from importlib.resources import files as resource_files

    base = tmp_path / "microsalt"
    base.mkdir()
    cfg = MicroSALTConfig(
        slurm_header=SlurmHeader(
            time="1:00:00",
            threads="4",
            qos="normal",
            job_prefix="MLST",
            project="test",
            type="core",
        ),
        regex=Regex(mail_recipient="test@test.com", file_pattern=".*", verified_organisms=[]),
        folders=Folders(
            results=str(base / "results"),
            reports=str(base / "reports"),
            seqdata=str(base / "seqdata"),
            profiles=str(base / "profiles"),
            references=str(base / "references"),
            resistances=str(base / "resistances"),
            genomes=str(base / "genomes"),
            credentials=str(base / "credentials"),
        ),
        database=Database(SQLALCHEMY_DATABASE_URI=f"sqlite:///{base / 'microsalt.db'}"),
        threshold=Threshold(),
        pubmlst=PubMLSTCredentials(),
        pasteur=PasteurCredentials(),
        singularity=Singularity(),
        containers=Containers(),
    )
    cfg.folders.expec = str(resource_files("microSALT").joinpath("unique_references", "ExPEC.fsa"))
    return cfg


@pytest.fixture
def setup_config_file(setup_config: MicroSALTConfig, tmp_path: Path) -> Path:
    """Write the setup_config to a JSON file and return its path."""
    path = tmp_path / "config.json"
    path.write_text(setup_config.model_dump_json())
    return path


def invoke_root(config_file: Path, *args: str):
    """Invoke the root CLI group with the test config and additional arguments."""
    runner = CliRunner()
    return runner.invoke(root, ["--config", str(config_file), *args])


def test_root_loads_config(config_file):
    """root command loads without error when given a valid config file."""
    result = invoke_root(config_file, "--help")
    assert result.exit_code == 0


def test_root_version(config_file):
    """--version flag returns the package version."""
    from microSALT import __version__

    result = invoke_root(config_file, "--version")
    assert result.exit_code == 0
    assert __version__ in result.output


def test_config_fields_propagated(config: MicroSALTConfig, config_file: Path):
    """Config loaded from the JSON file matches the original MicroSALTConfig fields."""
    loaded = MicroSALTConfig.model_validate(json.loads(config_file.read_text()))
    assert loaded.folders.results == config.folders.results
    assert loaded.folders.reports == config.folders.reports
    assert loaded.database.SQLALCHEMY_DATABASE_URI == config.database.SQLALCHEMY_DATABASE_URI
    assert loaded.regex.mail_recipient == config.regex.mail_recipient
    assert loaded.slurm_header.project == config.slurm_header.project


def test_setup_command(setup_config: MicroSALTConfig, setup_config_file: Path):
    """setup command creates the configured directories and reports success."""
    result = invoke_root(setup_config_file, "setup")
    assert result.exit_code == 0, result.output
    assert "Directory setup complete" in result.output
    assert "Database tables created" in result.output
    assert Path(setup_config.folders.results).exists()
    assert Path(setup_config.folders.reports).exists()


def test_setup_command_creates_tables(setup_config: MicroSALTConfig, setup_config_file: Path):
    """setup command creates every ORM-defined table in the database."""
    result = invoke_root(setup_config_file, "setup")
    assert result.exit_code == 0, result.output

    engine = create_engine(setup_config.database.SQLALCHEMY_DATABASE_URI)
    created = inspect(engine).get_table_names()
    expected = set(Base.metadata.tables.keys())
    assert expected == set(created), f"Missing tables: {expected - set(created)}"


def test_finish_exits_nonzero_when_scraper_raises(
    config: MicroSALTConfig, config_file: Path, tmp_path: Path
):
    """finish must exit with a non-zero code when the scraper raises an exception.

    A non-zero exit code causes the mailjob's ``if finish_cmd`` condition to be
    false, preventing ``run_complete.out`` from being created.
    """
    sampleinfo = tmp_path / "sampleinfo.json"
    sampleinfo.write_text(
        json.dumps(
            [
                {
                    "CG_ID_project": "TST0001",
                    "CG_ID_sample": "TST0001A1",
                    "Customer_ID_project": "100100",
                    "Customer_ID_sample": "10XY123456",
                    "Customer_ID": "cust000",
                    "application_tag": "SOMTIN100",
                    "date_arrival": "0001-01-01 00:00:00",
                    "date_libprep": "0001-01-01 00:00:00",
                    "date_sequencing": "0001-01-01 00:00:00",
                    "method_libprep": "Not in LIMS",
                    "method_sequencing": "Not in LIMS",
                    "organism": "Staphylococcus aureus",
                    "priority": "standard",
                    "reference": "None",
                }
            ]
        )
    )
    input_dir = tmp_path / "input"
    input_dir.mkdir()

    with patch("microSALT.cli.Scraper") as mock_scraper_cls:
        mock_scraper_cls.return_value.scrape_sample.side_effect = RuntimeError("scraper failed")
        result = invoke_root(
            config_file,
            "utils",
            "finish",
            str(sampleinfo),
            "--input",
            str(input_dir),
        )

    assert result.exit_code != 0, "finish should exit non-zero when scraper raises"


def test_finish_exits_nonzero_when_input_missing(
    config: MicroSALTConfig, config_file: Path, tmp_path: Path
):
    """finish must exit non-zero when the input folder does not exist.

    Mirrors the mailjob scenario where a bad path would otherwise silently
    succeed and write ``run_complete.out``.
    """
    sampleinfo = tmp_path / "sampleinfo.json"
    sampleinfo.write_text(json.dumps({"CG_ID_sample": "TST0001A1", "CG_ID_project": "TST0001"}))

    result = invoke_root(
        config_file,
        "utils",
        "finish",
        str(sampleinfo),
        "--input",
        str(tmp_path / "does_not_exist"),
    )

    assert result.exit_code != 0, "finish should exit non-zero for a missing input folder"
