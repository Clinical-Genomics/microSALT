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
            sampleinfo=testdata,
            run_settings={"pool": ["AAA1234A1", "AAA1234A2"], "input": "/tmp/AAA1234"},
        )
        jc.project_job()
