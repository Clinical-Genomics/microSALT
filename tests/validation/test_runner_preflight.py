"""Tests for the validation runner's environment preflight and suite registry."""

import importlib.util
import sys
import types
from importlib.machinery import ModuleSpec
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts" / "validation"))

RUNNER_PATH = (
    Path(__file__).resolve().parents[2] / "scripts" / "validation" / "run_validation.py"
)
SPEC = importlib.util.spec_from_file_location("validation_runner", RUNNER_PATH)
RUNNER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(RUNNER)


def test_preflight_accepts_regular_package():
    RUNNER.preflight()


def test_preflight_rejects_namespace_package(monkeypatch):
    fake = types.ModuleType("microSALT")
    fake.__spec__ = ModuleSpec("microSALT", None, is_package=True)
    fake.__path__ = ["/site-packages/microSALT"]
    monkeypatch.setitem(sys.modules, "microSALT", fake)
    with pytest.raises(SystemExit, match="namespace package.*site-packages/microSALT"):
        RUNNER.preflight()


def test_preflight_recovers_partially_submitted_job_ids():
    from core.runner import submitted_job_ids

    output = (
        "Error: SLURM rejected submission. Already submitted job IDs: 12345, 67890. "
        "The reporting job was not submitted."
    )
    assert submitted_job_ids(output) == ["12345", "67890"]


def test_registry_reports_broken_suite_modules(tmp_path, monkeypatch):
    package = tmp_path / "suites_probe" / "microsalt_x"
    package.mkdir(parents=True)
    (tmp_path / "suites_probe" / "__init__.py").write_text(
        (
            Path(__file__).resolve().parents[2]
            / "scripts/validation/suites/__init__.py"
        ).read_text()
    )
    (tmp_path / "suites_probe" / "stability.py").write_text(
        "from suites.stability import StabilitySuite\n"
    )
    (package / "__init__.py").write_text("")
    (package / "good.py").write_text(
        "from core.suite import Suite\n\nclass GoodSuite(Suite):\n    name = 'good'\n    title = 'Good'\n"
    )
    (package / "bad.py").write_text("import does_not_exist\n")
    monkeypatch.syspath_prepend(str(tmp_path))
    probe = __import__("suites_probe")
    assert "good" in probe.REGISTRY
    assert "stability" in probe.REGISTRY
    assert list(probe.BROKEN) == ["suites_probe.microsalt_x.bad"]
    assert "does_not_exist" in probe.BROKEN["suites_probe.microsalt_x.bad"]
