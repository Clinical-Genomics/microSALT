"""Tests for the CLI entry point — verifies config loading and propagation."""

import json
import pathlib

import pytest
from click.testing import CliRunner

from microSALT.cli import root
from microSALT.config import MicroSALTConfig


@pytest.fixture
def config_file(config: MicroSALTConfig, tmp_path: pathlib.Path) -> pathlib.Path:
    """Write the session config to a temporary JSON file for CLI consumption."""
    path = tmp_path / "config.json"
    path.write_text(config.model_dump_json())
    return path


def invoke_root(config_file: pathlib.Path, *args: str):
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


def test_config_fields_propagated(config: MicroSALTConfig, config_file: pathlib.Path):
    """Config loaded from the JSON file matches the original MicroSALTConfig fields."""
    loaded = MicroSALTConfig.model_validate(json.loads(config_file.read_text()))
    assert loaded.folders.results == config.folders.results
    assert loaded.folders.reports == config.folders.reports
    assert loaded.database.SQLALCHEMY_DATABASE_URI == config.database.SQLALCHEMY_DATABASE_URI
    assert loaded.regex.mail_recipient == config.regex.mail_recipient
    assert loaded.slurm_header.project == config.slurm_header.project


def test_setup_command(tmp_path):
    """setup command creates the configured directories without error."""
    from importlib.resources import files as resource_files

    from microSALT.config import (
        Containers,
        Database,
        Folders,
        PasteurCredentials,
        PubMLSTCredentials,
        Regex,
        Singularity,
        SlurmHeader,
        Threshold,
    )

    base = tmp_path / "microsalt"
    base.mkdir()

    cfg = MicroSALTConfig(
        slurm_header=SlurmHeader(time="1:00:00", threads="4", qos="normal",
                                 job_prefix="MLST", project="test", type="core"),
        regex=Regex(mail_recipient="test@test.com", file_pattern=".*", verified_organisms=[]),
        folders=Folders(
            results=str(base / "results"),
            reports=str(base / "reports"),
            log_file=str(base / "microsalt.log"),
            seqdata=str(base / "seqdata"),
            profiles=str(base / "profiles"),
            references=str(base / "references"),
            resistances=str(base / "resistances"),
            genomes=str(base / "genomes"),
            credentials=str(base / "credentials"),
            adapters=str(base / "adapters"),
        ),
        database=Database(SQLALCHEMY_DATABASE_URI=f"sqlite:///{base}/microsalt.db"),
        threshold=Threshold(),
        pubmlst=PubMLSTCredentials(),
        pasteur=PasteurCredentials(),
        singularity=Singularity(),
        containers=Containers(),
    )
    cfg.folders.expec = str(resource_files("microSALT").joinpath("unique_references", "ExPEC.fsa"))

    config_path = tmp_path / "config.json"
    config_path.write_text(cfg.model_dump_json())

    result = invoke_root(config_path, "setup")
    assert result.exit_code == 0, result.output
    assert "Directory setup complete" in result.output
    assert (base / "results").exists()
    assert (base / "reports").exists()
