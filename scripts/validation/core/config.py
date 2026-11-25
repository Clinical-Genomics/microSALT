"""Validation run configuration."""

from functools import cached_property
from pathlib import Path

import yaml
from pydantic import BaseModel, model_validator

from microSALT.config import MicroSALTConfig, load_config


class Case(BaseModel):
    name: str
    project: str
    sampleinfo: Path
    fastq: Path
    old_typing_report: Path | None = None


class ValidationConfig(BaseModel):
    release: str
    microsalt_bin: Path
    microsalt_config: Path
    validation_root: Path
    cases: list[Case]
    suites: dict[str, dict] = {"stability": {}}
    manual_steps: dict[str, str] = {}

    @model_validator(mode="after")
    def _unique_cases(self) -> "ValidationConfig":
        names = [case.name for case in self.cases]
        if len(set(names)) != len(names):
            raise ValueError("case names must be unique")
        if "stability" not in self.suites:
            self.suites = {"stability": {}, **self.suites}
        return self

    def case(self, name: str) -> Case:
        for case in self.cases:
            if case.name == name:
                return case
        raise KeyError(f"unknown case {name!r}")

    @property
    def case_names(self) -> list[str]:
        return [case.name for case in self.cases]

    @cached_property
    def microsalt(self) -> MicroSALTConfig:
        return load_config(str(self.microsalt_config))

    @property
    def results_root(self) -> Path:
        return Path(self.microsalt.folders.results)

    @property
    def state_path(self) -> Path:
        return self.validation_root / "state.json"

    @property
    def log_dir(self) -> Path:
        return self.validation_root / "logs"

    @property
    def fixture_dir(self) -> Path:
        return self.validation_root / "fixtures"


def load_validation_config(path: str | Path) -> ValidationConfig:
    with open(path) as stream:
        data = yaml.safe_load(stream)
    config = ValidationConfig(**data)
    config.validation_root.mkdir(parents=True, exist_ok=True)
    config.log_dir.mkdir(exist_ok=True)
    config.fixture_dir.mkdir(exist_ok=True)
    return config
