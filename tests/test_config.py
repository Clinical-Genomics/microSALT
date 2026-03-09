import pathlib

import pytest
from pydantic import ValidationError

from microSALT.config import (
    MicroSALTConfig,
    load_config,
    Folders,
    Regex,
    SlurmHeader,
    Threshold,
    Database,
)

CONFIGEXAMPLE = str(pathlib.Path(__file__).parent.parent / "configExample.json")


def test_load_config_parses_example():
    """load_config() successfully parses configExample.json into a MicroSALTConfig."""
    cfg = load_config(CONFIGEXAMPLE)
    assert isinstance(cfg, MicroSALTConfig)


def test_config_sections_present():
    """All top-level config sections are populated after parsing."""
    cfg = load_config(CONFIGEXAMPLE)
    assert isinstance(cfg.slurm_header, SlurmHeader)
    assert isinstance(cfg.regex, Regex)
    assert isinstance(cfg.folders, Folders)
    assert isinstance(cfg.database, Database)
    assert isinstance(cfg.threshold, Threshold)


def test_slurm_header_fields():
    cfg = load_config(CONFIGEXAMPLE)
    assert cfg.slurm_header.time
    assert cfg.slurm_header.threads
    assert cfg.slurm_header.qos
    assert cfg.slurm_header.job_prefix
    assert cfg.slurm_header.project
    assert cfg.slurm_header.type


def test_folders_fields():
    cfg = load_config(CONFIGEXAMPLE)
    assert cfg.folders.results
    assert cfg.folders.reports
    assert cfg.folders.log_file
    assert cfg.folders.seqdata
    assert cfg.folders.profiles
    assert cfg.folders.references
    assert cfg.folders.resistances
    assert cfg.folders.genomes
    assert cfg.folders.credentials
    assert cfg.folders.adapters


def test_expec_path_injected():
    """The expec path is derived from package data and injected by load_config."""
    cfg = load_config(CONFIGEXAMPLE)
    assert cfg.folders.expec
    assert "ExPEC.fsa" in cfg.folders.expec


def test_config_path_injected():
    """config_path is set to the resolved path of the config file."""
    cfg = load_config(CONFIGEXAMPLE)
    assert cfg.config_path
    assert cfg.config_path.endswith("configExample.json")


def test_threshold_defaults():
    cfg = load_config(CONFIGEXAMPLE)
    assert cfg.threshold.mlst_id == 100
    assert cfg.threshold.mlst_novel_id == 99.5
    assert cfg.threshold.mlst_span == 90


def test_runtime_defaults():
    """Runtime fields default to safe values before CLI sets them."""
    cfg = load_config(CONFIGEXAMPLE)
    assert cfg.dry is False
    assert cfg.config_path != ""


def test_missing_required_field_raises(tmp_path):
    """A config missing a required section raises a Pydantic ValidationError."""
    import json

    bad_config = tmp_path / "bad.json"
    bad_config.write_text(json.dumps({"slurm_header": {"time": "1:00:00"}}))
    with pytest.raises(ValidationError):
        load_config(str(bad_config))
