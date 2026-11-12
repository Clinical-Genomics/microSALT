"""Tests for the CLI entry point — verifies config loading and propagation."""

import json
import gzip
from pathlib import Path
from unittest.mock import patch

import pytest
from click.testing import CliRunner
from sqlalchemy import create_engine, inspect

from microSALT.cli import default_sampleinfo, root
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


@pytest.mark.parametrize("qc_option", [None, "--qc-only", "--qc_only"])
@pytest.mark.parametrize("sample_count", [1, 2])
def test_analyse_passes_qc_only_to_job_creator(config_file, tmp_path, qc_option, sample_count):
    """QC-only is opt-in and reaches the job creator for samples and projects."""
    sampleinfo = tmp_path / "sampleinfo.json"
    sampleinfo.write_text(
        json.dumps(
            [
                {
                    **default_sampleinfo,
                    "CG_ID_sample": f"TST0001A{i}",
                    "CG_ID_project": "TST0001",
                    "organism": "other",
                    "reference": "NC_000913.3",
                }
                for i in range(1, sample_count + 1)
            ]
        )
    )
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    for i in range(1, sample_count + 1):
        directory = input_dir / f"TST0001A{i}" if sample_count > 1 else input_dir
        directory.mkdir(exist_ok=True)
        for mate in (1, 2):
            with gzip.open(directory / f"TST0001A{i}_HVMHWDSXX_L1_{mate}.fastq.gz", "wt") as stream:
                stream.write("@read\nACGT\n+\nIIII\n")
    args = ["analyse", str(sampleinfo), "--input", str(input_dir), "--dry"]
    if qc_option:
        args.append(qc_option)

    with patch("microSALT.cli.Job_Creator") as creator, patch("microSALT.cli.Referencer"):
        result = invoke_root(config_file, *args)

    assert result.exit_code == 0, result.output
    assert creator.call_args.kwargs["run_settings"]["qc_only"] is (qc_option is not None)
    if sample_count == 1:
        creator.return_value.project_job.assert_called_once_with(single_sample=True)
    else:
        creator.return_value.project_job.assert_called_once_with()


def test_analyse_collects_input_and_reference_errors_without_submitting(config_file, tmp_path):
    samples = [
        {**default_sampleinfo, "CG_ID_sample": "TST0001A1", "reference": "NC_1"},
        {
            **default_sampleinfo,
            "CG_ID_sample": "TST0001A2",
            "reference": "NC_2",
            "date_arrival": "bad date",
        },
    ]
    sampleinfo = tmp_path / "sampleinfo.json"
    sampleinfo.write_text(json.dumps(samples))
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    with (
        patch("microSALT.cli.Referencer") as referencer,
        patch("microSALT.cli.Job_Creator") as creator,
        patch("microSALT.utils.job_creator.subprocess.Popen") as popen,
    ):
        referencer.return_value.prepare_references.return_value = [
            "Sample TST0001A1: reference NC_1 has invalid FASTA",
            "Sample TST0001A2: reference NC_2 could not be fetched",
        ]
        result = invoke_root(config_file, "analyse", str(sampleinfo), "--input", str(input_dir))
    assert result.exit_code != 0
    assert "No SLURM jobs were submitted" in result.output
    assert "TST0001A1" in result.output and "TST0001A2" in result.output
    assert "date_arrival" in result.output
    assert "invalid FASTA" in result.output and "could not be fetched" in result.output
    creator.assert_not_called()
    popen.assert_not_called()


def test_analyse_rejects_http_200_error_and_missing_mate(config_file, tmp_path):
    import io

    samples = [
        {**default_sampleinfo, "CG_ID_sample": f"TST0001A{i}", "reference": "NC_000913.3"}
        for i in (1, 2)
    ]
    sampleinfo = tmp_path / "sampleinfo.json"
    sampleinfo.write_text(json.dumps(samples))
    input_dir = tmp_path / "input"
    for sample in samples:
        directory = input_dir / sample["CG_ID_sample"]
        directory.mkdir(parents=True)
        with gzip.open(
            directory / f"{sample['CG_ID_sample']}_HVMHWDSXX_L1_1.fastq.gz", "wt"
        ) as stream:
            stream.write("@read\nACGT\n+\nIIII\n")
    with (
        patch(
            "microSALT.utils.referencer.Entrez.efetch",
            return_value=io.StringIO("Error: Failed to understand id"),
        ) as fetch,
        patch("microSALT.utils.job_creator.subprocess.Popen") as popen,
    ):
        result = invoke_root(
            config_file, "analyse", str(sampleinfo), "--input", str(input_dir), "--qc-only"
        )
    assert result.exit_code != 0
    assert "No SLURM jobs were submitted" in result.output
    assert "missing FASTQ mate" in result.output
    assert "invalid nucleotide FASTA" in result.output
    assert "TST0001A1, TST0001A2" in result.output
    fetch.assert_called_once()
    popen.assert_not_called()


def test_analyse_collects_assembly_resolution_failure_without_submitting(config_file, tmp_path):
    import io

    samples = [
        {**default_sampleinfo, "CG_ID_sample": f"TST0001A{i}", "reference": "GCF_000262545.1"}
        for i in (1, 2)
    ]
    sampleinfo = tmp_path / "sampleinfo.json"
    sampleinfo.write_text(json.dumps(samples))
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    with (
        patch(
            "microSALT.utils.referencer.Entrez.esearch",
            return_value=io.StringIO(json.dumps({"esearchresult": {"count": "0", "idlist": []}})),
        ) as search,
        patch("microSALT.utils.job_creator.subprocess.Popen") as popen,
    ):
        result = invoke_root(
            config_file,
            "analyse",
            str(sampleinfo),
            "--input",
            str(input_dir),
            "--qc-only",
        )
    assert result.exit_code != 0
    assert "No SLURM jobs were submitted" in result.output
    assert "TST0001A1, TST0001A2" in result.output
    assert "did not resolve to exactly one NCBI Assembly record" in result.output
    search.assert_called_once()
    popen.assert_not_called()


def test_analyse_saves_normalized_reference_for_finish(config_file, tmp_path):
    from microSALT.utils.job_creator import Job_Creator

    source = {**default_sampleinfo, "reference": " NC_ 000913.3\n"}
    sampleinfo = tmp_path / "sampleinfo.json"
    sampleinfo.write_text(json.dumps(source))
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    for mate in (1, 2):
        with gzip.open(input_dir / f"XXX0000A1_HVMHWDSXX_L1_{mate}.fastq.gz", "wt") as stream:
            stream.write("@read\nACGT\n+\nIIII\n")
    output = []

    def prepare_sample(instance):
        path = Path(instance.finishdir)
        path.mkdir(parents=True, exist_ok=True)
        instance.batchfile = str(path / "runfile.sbatch")
        Path(instance.batchfile).write_text("#!/bin/sh\n")
        output.append(path)

    with (
        patch("microSALT.cli.Referencer") as referencer,
        patch.object(Job_Creator, "sample_job", autospec=True, side_effect=prepare_sample),
        patch.object(Job_Creator, "create_project"),
        patch("microSALT.utils.job_creator.subprocess.Popen") as popen,
    ):
        referencer.return_value.prepare_references.return_value = []
        result = invoke_root(
            config_file, "analyse", str(sampleinfo), "--input", str(input_dir), "--qc-only", "--dry"
        )
    assert result.exit_code == 0, result.output
    normalized = referencer.return_value.prepare_references.call_args.args[0]
    assert normalized[0]["reference"] == "NC_000913.3"
    assert json.loads((output[0] / "sampleinfo.json").read_text())["reference"] == "NC_000913.3"
    assert json.loads(sampleinfo.read_text())["reference"] == source["reference"]
    popen.assert_not_called()


def test_analyse_valid_qc_run_prepares_scripts_and_submits(config, config_file, tmp_path):
    reference = Path(config.folders.genomes) / "PREFLIGHT_TEST.fasta"
    reference.write_text(">chromosome\nACGT\n")
    for suffix in (".amb", ".ann", ".bwt", ".pac", ".sa", ".fai"):
        Path(f"{reference}{suffix}").write_text("index")
    source = {
        **default_sampleinfo,
        "CG_ID_sample": "QCINT001A1",
        "CG_ID_project": "QCINT001",
        "reference": " PREFLIGHT_TEST\n",
        "organism": "other",
    }
    sampleinfo = tmp_path / "sampleinfo.json"
    sampleinfo.write_text(json.dumps(source))
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    for mate in (1, 2):
        with gzip.open(input_dir / f"QCINT001A1_HVMHWDSXX_L1_{mate}.fastq.gz", "wt") as stream:
            stream.write("@read\nACGT\n+\nIIII\n")
    with patch("microSALT.utils.job_creator.subprocess.Popen") as popen:
        popen.return_value.communicate.return_value = (b"Submitted batch job 12345", b"")
        popen.return_value.returncode = 0
        result = invoke_root(
            config_file,
            "analyse",
            str(sampleinfo),
            "--input",
            str(input_dir),
            "--qc-only",
            "--skip_update",
        )
    assert result.exit_code == 0, result.output
    assert popen.call_count == 2  # Sample, then reporting job; no reference subprocesses.
    run_dir = next(Path(config.folders.results).glob("QCINT001A1_*"))
    script = (run_dir / "runfile.sbatch").read_text()
    assert "PREFLIGHT_TEST.fasta" in script
    assert "blastn" not in script and "skesa" not in script
    assert "--track qc" in (run_dir / "mailjob.sh").read_text()
    assert json.loads((run_dir / "sampleinfo.json").read_text())["reference"] == "PREFLIGHT_TEST"


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


@pytest.mark.parametrize("track, expected_qc_only", [("qc", True), ("default", False)])
def test_finish_track_sets_scraper_qc_only(
    config: MicroSALTConfig, config_file: Path, tmp_path: Path, track: str, expected_qc_only: bool
):
    """finish --track qc must build the Scraper in QC-only mode."""
    sampleinfo = tmp_path / "sampleinfo.json"
    sampleinfo.write_text(json.dumps({"CG_ID_sample": "TST0001A1", "CG_ID_project": "TST0001"}))
    input_dir = tmp_path / "input"
    input_dir.mkdir()

    with patch("microSALT.cli.Scraper") as mock_scraper_cls, patch("microSALT.cli.Reporter"):
        result = invoke_root(
            config_file,
            "utils",
            "finish",
            str(sampleinfo),
            "--input",
            str(input_dir),
            "--track",
            track,
        )

    assert result.exit_code == 0, result.output
    assert mock_scraper_cls.call_args.kwargs["qc_only"] is expected_qc_only
