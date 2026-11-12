#!/usr/bin/env python

from unittest import mock
from unittest.mock import patch
import gzip
from pathlib import Path

import pytest

from microSALT.utils.job_creator import Job_Creator
from microSALT.exc.exceptions import AnalysisSetupError, JobCreationError


def fake_search(int):
    return "fake"


def test_verify_fastq(config, logger, testdata, tmp_path):
    names = [
        "ACC6438A3_HVMHWDSXX_L1_1.fastq.gz",
        "ACC6438A3_HVMHWDSXX_L1_2.fastq.gz",
        "ACC6438A3_HVMHWDSXX_L2_1.fastq.gz",
        "ACC6438A3_HVMHWDSXX_L2_2.fastq.gz",
    ]
    for name in names:
        with gzip.open(tmp_path / name, "wt") as stream:
            stream.write("@read\nACGT\n+\nIIII\n")

    jc = Job_Creator(
        log=logger,
        folders=config.folders,
        slurm_header=config.slurm_header,
        regex=config.regex,
        dry=False,
        config_path=config.config_path,
        threshold=config.threshold,
        pubmlst=config.pubmlst,
        pasteur=config.pasteur,
        singularity=config.singularity,
        containers=config.containers,
        sampleinfo=testdata,
        run_settings={"input": str(tmp_path)},
    )
    t = jc.verify_fastq()
    assert len(t) > 0


@patch("re.search")
@patch("microSALT.utils.job_creator.glob.glob")
def test_blast_subset(glob_search, research, config, logger, testdata):
    jc = Job_Creator(
        log=logger,
        folders=config.folders,
        slurm_header=config.slurm_header,
        regex=config.regex,
        dry=False,
        config_path=config.config_path,
        threshold=config.threshold,
        pubmlst=config.pubmlst,
        pasteur=config.pasteur,
        singularity=config.singularity,
        containers=config.containers,
        sampleinfo=testdata,
        run_settings={"input": "/tmp/"},
    )
    researcha = mock.MagicMock()
    researcha.group = fake_search
    research.return_value = researcha
    glob_search.return_value = ["/a/a/a", "/a/a/b", "/a/a/c"]

    jc.blast_subset("mlst", "/tmp/*")
    jc.blast_subset("other", "/tmp/*")
    outfile = open(jc.get_sbatch(), "r")
    count = 0
    for x in outfile.readlines():
        if "blastn -db" in x:
            count = count + 1
    assert count > 0


@patch("subprocess.Popen")
def test_create_snpsection(subproc, config, logger, testdata):
    # Sets up subprocess mocking
    process_mock = mock.Mock()
    attrs = {"communicate.return_value": ("output 123456789", "error")}
    process_mock.configure_mock(**attrs)
    subproc.return_value = process_mock

    testdata = [testdata[0]]
    jc = Job_Creator(
        log=logger,
        folders=config.folders,
        slurm_header=config.slurm_header,
        regex=config.regex,
        dry=False,
        config_path=config.config_path,
        threshold=config.threshold,
        pubmlst=config.pubmlst,
        pasteur=config.pasteur,
        singularity=config.singularity,
        containers=config.containers,
        sampleinfo=testdata,
        run_settings={"input": ["AAA1234A1", "AAA1234A2"]},
    )
    jc.snp_job()
    outfile = open(jc.get_sbatch(), "r")
    count = 0
    for x in outfile.readlines():
        if "# SNP pair-wise distance" in x:
            count = count + 1
    assert count > 0


@patch("subprocess.Popen")
def test_project_job(subproc, config, logger, testdata, tmp_path):
    # Sets up subprocess mocking
    process_mock = mock.Mock()
    attrs = {"communicate.return_value": (b"Submitted batch job 123456789", b""), "returncode": 0}
    process_mock.configure_mock(**attrs)
    subproc.return_value = process_mock

    def prepare_sample(instance):
        Path(instance.finishdir).mkdir(parents=True, exist_ok=True)
        instance.batchfile = str(Path(instance.finishdir) / "runfile.sbatch")
        Path(instance.batchfile).write_text("#!/bin/sh\n")

    with patch.object(Job_Creator, "sample_job", autospec=True, side_effect=prepare_sample):
        jc = Job_Creator(
            log=logger,
            folders=config.folders,
            slurm_header=config.slurm_header,
            regex=config.regex,
            dry=False,
            config_path=config.config_path,
            threshold=config.threshold,
            pubmlst=config.pubmlst,
            pasteur=config.pasteur,
            singularity=config.singularity,
            containers=config.containers,
            sampleinfo=testdata,
            run_settings={"input": str(tmp_path / "input"), "finishdir": str(tmp_path / "output")},
        )
        jc.project_job()


def test_singularity_exec_binds_finishdir(config, logger, testdata):
    """finishdir is automatically added to the --bind list of every singularity exec call."""
    jc = Job_Creator(
        log=logger,
        folders=config.folders,
        slurm_header=config.slurm_header,
        regex=config.regex,
        dry=False,
        config_path=config.config_path,
        threshold=config.threshold,
        pubmlst=config.pubmlst,
        pasteur=config.pasteur,
        singularity=config.singularity,
        containers=config.containers,
        sampleinfo=testdata,
        run_settings={"input": "/tmp/", "finishdir": "/tmp/test_runfolder"},
    )
    cmd = jc._singularity_exec("blast", "blastn -help")
    assert "/tmp/test_runfolder" in cmd


def test_singularity_exec_does_not_duplicate_finishdir(config, logger, testdata):
    """finishdir is not listed twice when it already appears in singularity.bind_paths."""
    from microSALT.config import Singularity

    singularity = Singularity(bind_paths=["/tmp/test_runfolder", "/data"])
    jc = Job_Creator(
        log=logger,
        folders=config.folders,
        slurm_header=config.slurm_header,
        regex=config.regex,
        dry=False,
        config_path=config.config_path,
        threshold=config.threshold,
        pubmlst=config.pubmlst,
        pasteur=config.pasteur,
        singularity=singularity,
        containers=config.containers,
        sampleinfo=testdata,
        run_settings={"input": "/tmp/", "finishdir": "/tmp/test_runfolder"},
    )
    cmd = jc._singularity_exec("blast", "blastn -help")
    assert cmd.count("/tmp/test_runfolder") == 1


def _make_jc(config, logger, testdata, tmp_path) -> Job_Creator:
    return Job_Creator(
        log=logger,
        folders=config.folders,
        slurm_header=config.slurm_header,
        regex=config.regex,
        dry=False,
        config_path=config.config_path,
        threshold=config.threshold,
        pubmlst=config.pubmlst,
        pasteur=config.pasteur,
        singularity=config.singularity,
        containers=config.containers,
        sampleinfo=testdata,
        run_settings={"input": "/tmp/", "finishdir": str(tmp_path)},
    )


def test_build_finish_cmd_uses_existing_executable(config, logger, testdata, tmp_path):
    """The binary resolved by _build_finish_cmd must exist on the filesystem."""
    import pathlib

    jc = _make_jc(config, logger, testdata, tmp_path)

    cmd = jc._build_finish_cmd("default")

    bin_path = cmd.split()[0]
    assert pathlib.Path(bin_path).exists(), f"Binary not found on disk: {bin_path}"


def test_build_finish_cmd_contains_expected_arguments(config, logger, testdata, tmp_path):
    """_build_finish_cmd must embed the correct report flag and paths."""
    jc = _make_jc(config, logger, testdata, tmp_path)

    cmd = jc._build_finish_cmd("qc")

    assert "utils finish" in cmd
    assert "--report qc" in cmd
    assert str(tmp_path) in cmd


@pytest.mark.parametrize("qc_only", [False, True])
def test_finish_job_propagates_qc_track_and_report(config, logger, testdata, tmp_path, qc_only):
    """The submitted finish script must select QC scraping as well as the QC report."""
    jc = _make_jc(config, logger, testdata, tmp_path)
    jc.qc_only = qc_only
    with patch("microSALT.utils.job_creator.subprocess.Popen") as popen:
        popen.return_value.communicate.return_value = (b"Submitted batch job 12345", b"")
        popen.return_value.returncode = 0
        jc.finish_job(["123"])

    script = (tmp_path / "mailjob.sh").read_text()
    assert f"--report {'qc' if qc_only else 'default'}" in script
    assert ("--track qc" in script) is qc_only


@pytest.mark.parametrize("qc_only", [False, True])
def test_sample_job_qc_only_skips_assembly_and_typing(config, logger, testdata, tmp_path, qc_only):
    """Read/alignment QC is retained; assembly and BLAST are full-analysis only."""
    jc = _make_jc(config, logger, testdata, tmp_path)
    jc.qc_only = qc_only
    with (
        patch.object(jc, "create_preprocsection") as preproc,
        patch.object(jc, "create_variantsection") as alignment,
        patch.object(jc, "create_assemblysection") as assembly,
        patch.object(jc, "create_assemblystats_section") as quast,
        patch.object(jc, "create_blast_search") as blast,
        patch.object(jc, "create_sample"),
    ):
        jc.sample_job()

    preproc.assert_called_once_with()
    alignment.assert_called_once_with()
    assert assembly.call_count == quast.call_count == blast.call_count == int(not qc_only)


@pytest.mark.parametrize("qc_only", [False, True])
def test_project_job_passes_qc_only_to_samples(config, logger, testdata, tmp_path, qc_only):
    """Every child job inherits the project's QC-only setting."""
    jc = _make_jc(config, logger, testdata, tmp_path)
    jc.run_settings["qc_only"] = qc_only
    jc.dry = True
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    for sample in testdata:
        (input_dir / sample["CG_ID_sample"]).mkdir()
    jc.indir = str(input_dir)

    with (
        patch.object(Job_Creator, "sample_job", autospec=True) as sample_job,
        patch.object(jc, "create_project"),
        patch("microSALT.utils.job_creator.os.path.isfile", return_value=True),
    ):
        jc.project_job()

    assert sample_job.call_count == len(testdata)
    assert all(call.args[0].qc_only is qc_only for call in sample_job.call_args_list)


def test_project_preparation_collects_errors_before_any_submission(
    config, logger, testdata, tmp_path
):
    jc = _make_jc(config, logger, testdata, tmp_path)
    visited = []

    def prepare(instance):
        visited.append(instance.name)
        if instance.name != testdata[0]["CG_ID_sample"]:
            raise RuntimeError(f"cannot prepare {instance.name}")
        Path(instance.finishdir).mkdir(parents=True, exist_ok=True)
        instance.batchfile = str(Path(instance.finishdir) / "runfile.sbatch")
        Path(instance.batchfile).write_text("#!/bin/sh\n")

    with (
        patch.object(Job_Creator, "sample_job", autospec=True, side_effect=prepare),
        patch.object(jc, "create_project", side_effect=RuntimeError("database unavailable")),
        patch.object(jc, "_prepare_finish_job", side_effect=RuntimeError("cannot write mailjob")),
        patch("microSALT.utils.job_creator.subprocess.Popen") as popen,
        pytest.raises(AnalysisSetupError) as failure,
    ):
        jc.project_job()
    assert visited == [sample["CG_ID_sample"] for sample in testdata]
    assert "database unavailable" in str(failure.value)
    assert "cannot write mailjob" in str(failure.value)
    assert all(sample["CG_ID_sample"] in str(failure.value) for sample in testdata[1:])
    popen.assert_not_called()


def test_single_sample_setup_failure_never_submits(config, logger, testdata, tmp_path):
    jc = _make_jc(config, logger, [testdata[0]], tmp_path)
    with (
        patch.object(jc, "sample_job", side_effect=RuntimeError("missing input")),
        patch.object(jc, "create_project"),
        patch("microSALT.utils.job_creator.subprocess.Popen") as popen,
        pytest.raises(AnalysisSetupError, match="missing input"),
    ):
        jc.project_job(single_sample=True)
    popen.assert_not_called()


def test_all_jobs_are_prepared_before_submission(config, logger, testdata, tmp_path):
    jc = _make_jc(config, logger, testdata, tmp_path)
    events = []

    def prepare(instance):
        events.append(f"prepare {instance.name}")
        Path(instance.finishdir).mkdir(parents=True, exist_ok=True)
        instance.batchfile = str(Path(instance.finishdir) / "runfile.sbatch")
        Path(instance.batchfile).write_text("#!/bin/sh\n")

    def submit(command):
        events.append("submit")
        return "12345"

    with (
        patch.object(Job_Creator, "sample_job", autospec=True, side_effect=prepare),
        patch.object(jc, "create_project"),
        patch.object(
            jc, "_prepare_finish_job", side_effect=lambda: events.append("prepare report")
        ),
        patch.object(jc, "_submit_sbatch", side_effect=submit),
        patch.object(jc, "finish_job") as finish,
    ):
        jc.project_job()
    assert events[: len(testdata)] == [f"prepare {sample['CG_ID_sample']}" for sample in testdata]
    assert events[len(testdata)] == "prepare report"
    assert events[len(testdata) + 1 :] == ["submit"] * len(testdata)
    finish.assert_called_once_with(["12345"] * len(testdata), False, prepared=True)


def test_failed_slurm_submission_reports_accepted_jobs(config, logger, testdata, tmp_path):
    jc = _make_jc(config, logger, testdata, tmp_path)
    with (
        patch.object(Job_Creator, "sample_job"),
        patch("microSALT.utils.job_creator.os.path.isfile", return_value=True),
        patch.object(jc, "create_project"),
        patch.object(jc, "_prepare_finish_job"),
        patch.object(
            jc, "_submit_sbatch", side_effect=["12345", JobCreationError("queue rejected")]
        ),
        patch.object(jc, "finish_job") as finish,
        pytest.raises(JobCreationError, match="Already submitted job IDs: 12345"),
    ):
        jc.project_job()
    finish.assert_not_called()


def test_sample_setup_failure_preserves_existing_results(config, logger, testdata, tmp_path):
    jc = _make_jc(config, logger, [testdata[0]], tmp_path)
    existing = tmp_path / "user-results.txt"
    existing.write_text("keep")
    with patch.object(jc, "create_preprocsection", side_effect=RuntimeError("bad input")):
        with pytest.raises(RuntimeError, match="bad input"):
            jc.sample_job()
    assert existing.read_text() == "keep"


def test_database_failure_is_not_swallowed_by_sample_job(config, logger, testdata, tmp_path):
    jc = _make_jc(config, logger, [testdata[0]], tmp_path)
    jc.qc_only = True
    with (
        patch.object(jc, "create_preprocsection"),
        patch.object(jc, "create_variantsection"),
        patch.object(jc, "create_sample", side_effect=RuntimeError("database write failed")),
        pytest.raises(RuntimeError, match="database write failed"),
    ):
        jc.sample_job()


def test_reporting_submission_failure_is_explicit(config, logger, testdata, tmp_path):
    jc = _make_jc(config, logger, testdata, tmp_path)
    with (
        patch.object(jc, "_submit_sbatch", side_effect=JobCreationError("queue rejected")),
        pytest.raises(JobCreationError, match="Already submitted analysis job IDs: 12345"),
    ):
        jc.finish_job(["12345"])


def test_empty_analysis_job_list_never_submits_report(config, logger, testdata, tmp_path):
    jc = _make_jc(config, logger, testdata, tmp_path)
    with patch.object(jc, "_submit_sbatch") as submit:
        with pytest.raises(JobCreationError, match="No analysis job IDs"):
            jc.finish_job([])
    submit.assert_not_called()


def test_setup_logger_creates_one_log_for_project(config, logger, testdata, tmp_path):
    """project_job creates exactly one job_creator.log in finishdir via _attach_log_file;
    construction alone must not touch the filesystem.
    """
    jc = Job_Creator(
        log=logger,
        folders=config.folders,
        slurm_header=config.slurm_header,
        regex=config.regex,
        dry=False,
        config_path=config.config_path,
        threshold=config.threshold,
        pubmlst=config.pubmlst,
        pasteur=config.pasteur,
        singularity=config.singularity,
        containers=config.containers,
        sampleinfo=testdata,
        run_settings={"input": "/tmp/", "finishdir": str(tmp_path)},
    )
    assert jc.logger.name.startswith("job_creator.")
    assert not (
        tmp_path / "job_creator.log"
    ).exists(), "log file must not be created during construction"
    jc._attach_log_file()
    assert (tmp_path / "job_creator.log").exists()


def test_setup_logger_sample_reuses_project_logger(config, logger, testdata, tmp_path):
    """Sample-level Job_Creator receiving a job_creator.* logger must reuse it without
    creating an extra directory or log file.
    """
    project_dir = tmp_path / "project"
    project_dir.mkdir()
    sample_dir = tmp_path / "sample1"

    project_jc = Job_Creator(
        log=logger,
        folders=config.folders,
        slurm_header=config.slurm_header,
        regex=config.regex,
        dry=False,
        config_path=config.config_path,
        threshold=config.threshold,
        pubmlst=config.pubmlst,
        pasteur=config.pasteur,
        singularity=config.singularity,
        containers=config.containers,
        sampleinfo=testdata,
        run_settings={"input": "/tmp/", "finishdir": str(project_dir)},
    )

    # Sample instance receives the project's job_creator logger
    sample_jc = Job_Creator(
        log=project_jc.logger,
        folders=config.folders,
        slurm_header=config.slurm_header,
        regex=config.regex,
        dry=False,
        config_path=config.config_path,
        threshold=config.threshold,
        pubmlst=config.pubmlst,
        pasteur=config.pasteur,
        singularity=config.singularity,
        containers=config.containers,
        sampleinfo=testdata,
        run_settings={"input": "/tmp/", "finishdir": str(sample_dir)},
    )

    assert sample_jc.logger is project_jc.logger
    assert not sample_dir.exists(), "sample_dir must not be created by _setup_logger"
    assert not (sample_dir / "job_creator.log").exists()


def test_setup_logger_each_project_gets_fresh_logger(config, logger, testdata, tmp_path):
    """Two distinct project-level instances must each get their own fresh logger that writes
    to their own finishdir, even when run in the same process (no stale global-cache handler).
    """
    dir_a = tmp_path / "run_a"
    dir_b = tmp_path / "run_b"
    dir_a.mkdir()
    dir_b.mkdir()

    jc_a = Job_Creator(
        log=logger,
        folders=config.folders,
        slurm_header=config.slurm_header,
        regex=config.regex,
        dry=False,
        config_path=config.config_path,
        threshold=config.threshold,
        pubmlst=config.pubmlst,
        pasteur=config.pasteur,
        singularity=config.singularity,
        containers=config.containers,
        sampleinfo=testdata,
        run_settings={"input": "/tmp/", "finishdir": str(dir_a)},
    )
    jc_b = Job_Creator(
        log=logger,
        folders=config.folders,
        slurm_header=config.slurm_header,
        regex=config.regex,
        dry=False,
        config_path=config.config_path,
        threshold=config.threshold,
        pubmlst=config.pubmlst,
        pasteur=config.pasteur,
        singularity=config.singularity,
        containers=config.containers,
        sampleinfo=testdata,
        run_settings={"input": "/tmp/", "finishdir": str(dir_b)},
    )

    assert jc_a.logger is not jc_b.logger
    jc_a._attach_log_file()
    jc_b._attach_log_file()
    assert (dir_a / "job_creator.log").exists()
    assert (dir_b / "job_creator.log").exists()
