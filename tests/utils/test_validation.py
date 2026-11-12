import gzip
from pathlib import Path
from unittest.mock import patch

import pytest

from microSALT.cli import default_sampleinfo
from microSALT.utils.validation import (
    find_fastq_files,
    is_safe_identifier,
    normalize_sample,
    validate_fastq_mate,
    validate_input_file,
    validate_input_folders,
    validate_required_fields,
    validate_sample_dates,
    validate_sample_groups,
    validate_sample_identifiers,
    validate_sample_input,
    validate_sample_metadata,
    validate_sampleinfo,
    verify_fastq,
)


@pytest.fixture
def fastq_pair(tmp_path):
    for mate in (1, 2):
        with gzip.open(tmp_path / f"sample_R{mate}.fastq.gz", "wt") as stream:
            stream.write("@read\nACGT\n+\nIIII\n")
    return tmp_path


def test_find_fastq_files_filters_and_sorts(fastq_pair):
    (fastq_pair / "notes.txt").write_text("not an input")
    assert find_fastq_files(fastq_pair, r"sample_R(\d).fastq.gz") == [
        "sample_R1.fastq.gz",
        "sample_R2.fastq.gz",
    ]


def test_find_fastq_files_requires_matching_inputs(tmp_path):
    (tmp_path / "notes.txt").write_text("not an input")
    with pytest.raises(ValueError, match="No FASTQ files"):
        find_fastq_files(tmp_path, r"sample_R(\d).fastq.gz")


@pytest.mark.parametrize(
    "mate, partner", [("1", "2"), ("2", "1"), ("forward", "reverse"), ("reverse", "forward")]
)
def test_validate_fastq_mate_accepts_both_naming_conventions(mate, partner):
    names = [f"sample_{mate}.fastq.gz", f"sample_{partner}.fastq.gz"]
    validate_fastq_mate(names[0], names, r"sample_(1|2|forward|reverse).fastq.gz")


def test_validate_fastq_mate_rejects_missing_partner():
    with pytest.raises(ValueError, match="missing FASTQ mate sample_R2.fastq.gz"):
        validate_fastq_mate("sample_R1.fastq.gz", ["sample_R1.fastq.gz"], r"sample_R(\d).fastq.gz")


def test_validate_fastq_mate_rejects_unsupported_identifier():
    with pytest.raises(ValueError, match="unsupported mate identifier '3'"):
        validate_fastq_mate("sample_R3.fastq.gz", ["sample_R3.fastq.gz"], r"sample_R(\d).fastq.gz")


@pytest.mark.parametrize("pattern", [r"different_(\d).fastq.gz", r"sample_R1.fastq.gz"])
def test_validate_fastq_mate_requires_a_captured_identifier(pattern):
    with pytest.raises(ValueError, match="capture a mate identifier"):
        validate_fastq_mate("sample_R1.fastq.gz", ["sample_R1.fastq.gz"], pattern)


def test_validate_input_file_checks_access_without_reading(fastq_pair):
    with patch.object(Path, "open") as open_file:
        validate_input_file(fastq_pair / "sample_R1.fastq.gz")
    open_file.assert_called_once_with("rb")
    handle = open_file.return_value.__enter__.return_value
    handle.read.assert_not_called()
    handle.readline.assert_not_called()
    handle.readlines.assert_not_called()
    handle.__iter__.assert_not_called()


@pytest.mark.parametrize("entry_type", ["missing", "directory"])
def test_validate_input_file_rejects_non_files(tmp_path, entry_type):
    path = tmp_path / "sample_R1.fastq.gz"
    if entry_type == "directory":
        path.mkdir()
    with pytest.raises(ValueError, match="not a regular file"):
        validate_input_file(path)


def test_unreadable_inputs_are_reported_for_both_mates(fastq_pair):
    with patch.object(Path, "open", side_effect=PermissionError("permission denied")):
        with pytest.raises(ValueError) as failure:
            verify_fastq(str(fastq_pair), r"sample_R(\d).fastq.gz")
    assert "sample_R1.fastq.gz" in str(failure.value)
    assert "sample_R2.fastq.gz" in str(failure.value)
    assert "permission denied" in str(failure.value)


def test_sample_validation_never_decompresses_fastq_inputs(fastq_pair):
    sampleinfo = {**default_sampleinfo, "reference": "NC_123"}
    with patch(
        "gzip.open", side_effect=AssertionError("FASTQ decompression on login node")
    ) as decompress:
        samples, errors = validate_sampleinfo(
            sampleinfo, default_sampleinfo, str(fastq_pair), r"sample_R(\d).fastq.gz"
        )
    assert not errors
    assert len(samples) == 1
    decompress.assert_not_called()


def test_normalize_sample_preserves_other_fields():
    sample = {"reference": " NC_ 123\n", "Customer_ID_sample": "sample with spaces"}
    assert normalize_sample(sample) == {
        "reference": "NC_123",
        "Customer_ID_sample": "sample with spaces",
    }
    assert sample["reference"] == " NC_ 123\n"


@pytest.mark.parametrize("value", ["sample-1", "sample.1", "NC_123.4"])
def test_is_safe_identifier_accepts_identifiers(value):
    assert is_safe_identifier(value)


@pytest.mark.parametrize(
    "value", [None, 42, "", ".", "..", "../reference", "sample name", "ref;command"]
)
def test_is_safe_identifier_rejects_unsafe_values(value):
    assert not is_safe_identifier(value)


def test_validate_required_fields_collects_all_invalid_fields():
    assert validate_required_fields(
        {"sample": "valid", "reference": " ", "customer": None},
        ["sample", "reference", "customer", "missing"],
    ) == [
        "missing or invalid reference",
        "missing or invalid customer",
        "missing or invalid missing",
    ]


def test_validate_sample_identifiers_checks_all_path_fields():
    errors = validate_sample_identifiers(
        {"CG_ID_sample": "../sample", "CG_ID_project": "project name", "reference": "ref;command"}
    )
    assert len(errors) == 3
    assert all("unsafe" in error for error in errors)


@pytest.mark.parametrize("reference", ["None", "none", "NULL"])
def test_validate_sample_identifiers_requires_alignment_reference(reference):
    assert validate_sample_identifiers({"reference": reference}) == [
        "a reference genome is required for alignment"
    ]


def test_validate_sample_dates_collects_all_invalid_dates():
    errors = validate_sample_dates(
        {
            "date_arrival": "yesterday",
            "date_libprep": "2025-02-29 12:00:00",
            "date_sequencing": "2024-02-29 12:00:00",
        }
    )
    assert len(errors) == 2
    assert "date_arrival" in errors[0] and "date_libprep" in errors[1]


def test_validate_sample_metadata_composes_independent_checks():
    errors = validate_sample_metadata(
        {"reference": "../reference", "date_arrival": "yesterday"}, ["Customer_ID"]
    )
    assert len(errors) == 3
    assert "Customer_ID" in errors[0]
    assert "unsafe reference" in errors[1]
    assert "date_arrival" in errors[2]


def test_validate_sample_input_accepts_valid_single_sample(fastq_pair):
    assert not validate_sample_input(
        {"CG_ID_sample": "sample"},
        str(fastq_pair),
        r"sample_R(\d).fastq.gz",
        multiple_samples=False,
    )


def test_validate_sample_input_requires_project_sample_directory(fastq_pair):
    errors = validate_sample_input(
        {"CG_ID_sample": "missing"},
        str(fastq_pair),
        r"sample_R(\d).fastq.gz",
        multiple_samples=True,
    )
    assert len(errors) == 1 and "missing" in errors[0]


def test_validate_sample_input_skips_invalid_identifier(fastq_pair):
    assert not validate_sample_input(
        {"CG_ID_sample": "../sample"},
        str(fastq_pair),
        r"sample_R(\d).fastq.gz",
        multiple_samples=True,
    )


def test_validate_sample_groups_collects_duplicates_and_mixed_projects():
    errors = validate_sample_groups(
        [
            {"CG_ID_sample": "sample", "CG_ID_project": "project1"},
            {"CG_ID_sample": "sample", "CG_ID_project": "project2"},
        ]
    )
    assert errors == ["Duplicate sample ID sample", "Sample info contains mixed projects"]


def test_validate_input_folders_ignores_known_samples_and_files(tmp_path):
    (tmp_path / "known").mkdir()
    (tmp_path / "unexpected").mkdir()
    (tmp_path / "notes.txt").write_text("not a sample")
    assert validate_input_folders(str(tmp_path), [{"CG_ID_sample": "known"}]) == [
        f"Input folder {tmp_path / 'unexpected'}: no counterpart in sample info"
    ]


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


@pytest.mark.parametrize("content", [b"", b"not gzip", gzip.compress(b"@read\nACGT\n+\nIII\n")])
def test_fastq_content_validation_is_left_to_analysis(fastq_pair, content):
    (fastq_pair / "sample_R1.fastq.gz").write_bytes(content)
    assert verify_fastq(str(fastq_pair), r"sample_R(\d).fastq.gz") == [
        "sample_R1.fastq.gz",
        "sample_R2.fastq.gz",
    ]


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


def test_gzip_corruption_is_left_to_analysis(fastq_pair):
    path = fastq_pair / "sample_R1.fastq.gz"
    path.write_bytes(path.read_bytes()[:-5])
    assert len(verify_fastq(str(fastq_pair), r"sample_R(\d).fastq.gz")) == 2


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
