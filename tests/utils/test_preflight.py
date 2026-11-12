import gzip

import pytest

from microSALT.cli import default_sampleinfo
from microSALT.utils.preflight import validate_sampleinfo, verify_fastq


@pytest.fixture
def fastq_pair(tmp_path):
    for mate in (1, 2):
        with gzip.open(tmp_path / f"sample_R{mate}.fastq.gz", "wt") as stream:
            stream.write("@read\nACGT\n+\nIIII\n")
    return tmp_path


@pytest.mark.parametrize("reference", [" NC_000913.3 ", "NC_ 000913.3", "\tNC_000913.3\n"])
def test_normalize_reference_without_modifying_source(fastq_pair, reference):
    source = {**default_sampleinfo, "reference": reference}
    samples, errors = validate_sampleinfo(
        source, default_sampleinfo, str(fastq_pair), r"sample_R(\d).fastq.gz"
    )
    assert not errors
    assert samples[0]["reference"] == "NC_000913.3"
    assert source["reference"] == reference


@pytest.mark.parametrize("reference", [" \t\n", None, "None", "../reference", "ref;touch-file"])
def test_reject_invalid_reference(fastq_pair, reference):
    _, errors = validate_sampleinfo(
        {**default_sampleinfo, "reference": reference},
        default_sampleinfo,
        str(fastq_pair),
        r"sample_R(\d).fastq.gz",
    )
    assert any("reference" in error for error in errors)


@pytest.mark.parametrize("data", [[], None, "sample", [None]])
def test_reject_invalid_sampleinfo(tmp_path, data):
    _, errors = validate_sampleinfo(
        data, default_sampleinfo, str(tmp_path), r"sample_R(\d).fastq.gz"
    )
    assert errors


def test_collect_errors_across_samples(tmp_path):
    samples = [
        {
            **default_sampleinfo,
            "CG_ID_sample": "sample1",
            "reference": " ",
            "date_arrival": "yesterday",
        },
        {**default_sampleinfo, "CG_ID_sample": "sample2", "Customer_ID": None},
    ]
    _, errors = validate_sampleinfo(
        samples, default_sampleinfo, str(tmp_path), r"sample_R(\d).fastq.gz"
    )
    summary = "\n".join(errors)
    assert "reference" in summary
    assert "date_arrival" in summary
    assert "Customer_ID" in summary
    assert "sample1" in summary and "sample2" in summary


def test_missing_fastq_mate_is_fatal(fastq_pair):
    (fastq_pair / "sample_R2.fastq.gz").unlink()
    with pytest.raises(ValueError, match="missing FASTQ mate sample_R2"):
        verify_fastq(str(fastq_pair), r"sample_R(\d).fastq.gz")


@pytest.mark.parametrize("content", ["", "@read\nACGT\n+\nIII\n", "@read\nACGT\n+"])
def test_invalid_fastq_is_fatal(fastq_pair, content):
    with gzip.open(fastq_pair / "sample_R1.fastq.gz", "wt") as stream:
        stream.write(content)
    with pytest.raises(ValueError, match="sample_R1"):
        verify_fastq(str(fastq_pair), r"sample_R(\d).fastq.gz")


def test_broken_symlink_is_fatal(fastq_pair):
    path = fastq_pair / "sample_R1.fastq.gz"
    path.unlink()
    path.symlink_to("missing.fastq.gz")
    with pytest.raises(ValueError, match="sample_R1"):
        verify_fastq(str(fastq_pair), r"sample_R(\d).fastq.gz")


def test_relative_symlink_is_supported(fastq_pair):
    path = fastq_pair / "sample_R1.fastq.gz"
    path.rename(fastq_pair / "reads.gz")
    path.symlink_to("reads.gz")
    assert len(verify_fastq(str(fastq_pair), r"sample_R(\d).fastq.gz")) == 2


def test_gzip_corruption_is_fatal(fastq_pair):
    path = fastq_pair / "sample_R1.fastq.gz"
    path.write_bytes(path.read_bytes()[:-5])
    with pytest.raises(ValueError, match="sample_R1"):
        verify_fastq(str(fastq_pair), r"sample_R(\d).fastq.gz")


def test_missing_input_does_not_inspect_working_directory():
    samples = [
        {**default_sampleinfo, "CG_ID_sample": f"sample{i}", "reference": "NC_123"} for i in (1, 2)
    ]
    _, errors = validate_sampleinfo(samples, default_sampleinfo, "", r"sample_R(\d).fastq.gz")
    assert len(errors) == 1
    assert "Sequence data folder" in errors[0]


def test_unknown_input_folder_is_reported(tmp_path):
    (tmp_path / "unexpected").mkdir()
    samples = [{**default_sampleinfo, "CG_ID_sample": f"sample{i}"} for i in (1, 2)]
    _, errors = validate_sampleinfo(
        samples, default_sampleinfo, str(tmp_path), r"sample_R(\d).fastq.gz"
    )
    assert any("unexpected" in error and "no counterpart" in error for error in errors)
