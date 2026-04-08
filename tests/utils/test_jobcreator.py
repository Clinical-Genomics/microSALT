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
