#!/usr/bin/env python

from unittest import mock
from unittest.mock import patch

from microSALT.utils.job_creator import Job_Creator


def fake_search(int):
    return "fake"


@patch("os.listdir")
@patch("os.stat")
@patch("gzip.open")
def test_verify_fastq(gopen, stat, listdir, config, logger, testdata):
    listdir.return_value = [
        "ACC6438A3_HVMHWDSXX_L1_1.fastq.gz",
        "ACC6438A3_HVMHWDSXX_L1_2.fastq.gz",
        "ACC6438A3_HVMHWDSXX_L2_2.fastq.gz",
        "ACC6438A3_HVMHWDSXX_L2_2.fastq.gz",
    ]
    stata = mock.MagicMock()
    stata.st_size = 2000
    stat.return_value = stata

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
def test_project_job(subproc, config, logger, testdata):
    # Sets up subprocess mocking
    process_mock = mock.Mock()
    attrs = {"communicate.return_value": ("output 123456789", "error")}
    process_mock.configure_mock(**attrs)
    subproc.return_value = process_mock

    with patch.dict("os.environ", {"CONDA_PREFIX": "/tmp/mock_conda"}):
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
            run_settings={"pool": ["AAA1234A1", "AAA1234A2"], "input": "/tmp/AAA1234"},
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
    assert not (tmp_path / "job_creator.log").exists(), (
        "log file must not be created during construction"
    )
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
