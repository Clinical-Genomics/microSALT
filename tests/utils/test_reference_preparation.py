"""Direct tests for resource preparation and its error-collection boundaries."""

from pathlib import Path
from unittest.mock import call, patch

import pytest

from microSALT.exc.exceptions import MLSTUnavailableError, ReferenceSetupError
from microSALT.utils.referencer import Referencer


@pytest.fixture
def genome_reference(referencer):
    path = Path(referencer.folders.genomes) / "NC_123.fasta"
    path.write_text(">chromosome\nACGT\n")
    for suffix in (".amb", ".ann", ".bwt", ".pac", ".sa", ".fai"):
        Path(f"{path}{suffix}").write_text("index")
    return path


@pytest.fixture
def mlst_reference_files(referencer, tmp_path):
    reference = "staphylococcus_aureus"
    referencer.folders.references = str(tmp_path / "loci")
    referencer.folders.profiles = str(tmp_path / "profiles")
    loci = Path(referencer.folders.references) / reference
    loci.mkdir(parents=True)
    Path(referencer.folders.profiles).mkdir()
    (loci / "locus.tfa").write_text(">locus_1\nACGT\n")
    (Path(referencer.folders.profiles) / reference).write_text("ST\tlocus\n1\t1\n")
    return reference


def test_group_samples_by_resource_deduplicates_and_preserves_sample_labels():
    samples = [
        {"CG_ID_sample": "A1", "reference": "NC_1"},
        {"CG_ID_sample": "A2", "reference": "NC_2"},
        {"CG_ID_sample": "A3", "reference": "NC_1"},
        {"reference": "NC_2"},
    ]
    assert Referencer._group_samples_by_resource(samples, "reference") == {
        "NC_1": ["A1", "A3"],
        "NC_2": ["A2", "unknown sample"],
    }


@pytest.mark.parametrize("value", [None, "", 42, "None", "NULL"])
def test_group_samples_by_resource_ignores_invalid_and_excluded_values(value):
    assert not Referencer._group_samples_by_resource(
        [{"reference": value}],
        "reference",
        excluded=("none", "null"),
    )


def test_group_samples_by_resource_can_group_organisms():
    assert Referencer._group_samples_by_resource(
        [{"CG_ID_sample": 42, "organism": "Other"}],
        "organism",
    ) == {"Other": ["42"]}


@pytest.mark.parametrize("skip_update", [False, True])
def test_prepare_genome_reference_reuses_valid_cached_genome(
    referencer, genome_reference, skip_update
):
    with (
        patch.object(referencer, "download_ncbi") as download,
        patch.object(referencer, "_index_genome") as index,
    ):
        referencer._prepare_genome_reference("NC_123", skip_update)
    download.assert_not_called()
    index.assert_not_called()


@pytest.mark.parametrize("reference", ["../reference", "bad name", ".", "..", "ref;command"])
def test_prepare_genome_reference_rejects_unsafe_identifier(referencer, reference):
    with patch.object(referencer, "download_ncbi") as download:
        with pytest.raises(ReferenceSetupError, match="invalid reference identifier"):
            referencer._prepare_genome_reference(reference, False)
    download.assert_not_called()


def test_prepare_genome_reference_downloads_missing_genome_before_validation(referencer):
    path = Path(referencer.folders.genomes) / "NC_123.fasta"
    with (
        patch.object(
            referencer,
            "download_ncbi",
            side_effect=lambda reference: path.write_text(">one\nACGT\n"),
        ) as download,
        patch.object(referencer, "_prepare_genome_indexes") as indexes,
    ):
        referencer._prepare_genome_reference("NC_123", False)
    download.assert_called_once_with("NC_123")
    indexes.assert_called_once_with(path, False)


def test_prepare_genome_reference_missing_offline_never_downloads(referencer):
    with patch.object(referencer, "download_ncbi") as download:
        with pytest.raises(ReferenceSetupError, match="missing reference.*downloads are disabled"):
            referencer._prepare_genome_reference("NC_missing", True)
    download.assert_not_called()


def test_prepare_genome_reference_corrupt_fasta_never_reindexes(referencer, genome_reference):
    genome_reference.write_text("Error: invalid id")
    with patch.object(referencer, "_prepare_genome_indexes") as indexes:
        with pytest.raises(ReferenceSetupError, match="invalid nucleotide FASTA"):
            referencer._prepare_genome_reference("NC_123", False)
    indexes.assert_not_called()


@pytest.mark.parametrize("skip_update", [False, True])
def test_prepare_genome_indexes_leaves_valid_indexes_untouched(
    referencer, genome_reference, skip_update
):
    with patch.object(referencer, "_index_genome") as index:
        referencer._prepare_genome_indexes(genome_reference, skip_update)
    index.assert_not_called()


def test_prepare_genome_indexes_rebuilds_missing_indexes(referencer, genome_reference):
    Path(f"{genome_reference}.fai").unlink()
    with patch.object(referencer, "_index_genome") as index:
        referencer._prepare_genome_indexes(genome_reference, False)
    index.assert_called_once_with(genome_reference)


def test_prepare_genome_indexes_refuses_rebuild_with_skip_update(referencer, genome_reference):
    Path(f"{genome_reference}.fai").unlink()
    with patch.object(referencer, "_index_genome") as index:
        with pytest.raises(ReferenceSetupError, match="missing or stale genome indexes.*.fai"):
            referencer._prepare_genome_indexes(genome_reference, True)
    index.assert_not_called()


def test_prepare_genome_references_collects_failures_and_attempts_later_genomes(referencer):
    with patch.object(
        referencer,
        "_prepare_genome_reference",
        side_effect=[ReferenceSetupError("fetch failed"), None],
    ) as prepare:
        errors = referencer._prepare_genome_references(
            {"NC_1": ["A1", "A2"], "NC_2": ["A3"]}, False
        )
    assert prepare.call_args_list == [call("NC_1", False), call("NC_2", False)]
    assert errors == ["Samples A1, A2, reference NC_1: fetch failed"]


@pytest.mark.parametrize("skip_update", [False, True])
def test_prepare_mlst_profile_reuses_cached_profile(referencer, skip_update):
    with (
        patch.object(referencer.db_access, "profiles", {"known": object()}),
        patch.object(referencer, "add_pubmlst") as download,
        patch.object(referencer, "create_new_profile_tables") as register,
    ):
        assert referencer._prepare_mlst_profile("Known organism", "known", skip_update) == "known"
    download.assert_not_called()
    register.assert_not_called()


def test_prepare_mlst_profile_registers_downloaded_profile(referencer):
    with (
        patch.object(referencer.db_access, "profiles", {}),
        patch.object(referencer, "add_pubmlst", return_value="new_organism") as download,
        patch.object(referencer, "create_new_profile_tables") as register,
    ):
        assert referencer._prepare_mlst_profile("New organism", None, False) == "new_organism"
    download.assert_called_once_with("New organism")
    register.assert_called_once_with(["new_organism"])


@pytest.mark.parametrize(
    "cached_reference, expected_error",
    [(None, MLSTUnavailableError), ("incomplete", ReferenceSetupError)],
)
def test_prepare_mlst_profile_distinguishes_unavailable_from_incomplete_cache(
    referencer, cached_reference, expected_error
):
    with (
        patch.object(referencer.db_access, "profiles", {}),
        patch.object(referencer, "add_pubmlst", side_effect=MLSTUnavailableError("No scheme")),
        patch.object(referencer, "create_new_profile_tables") as register,
        pytest.raises(expected_error, match="No scheme"),
    ):
        referencer._prepare_mlst_profile("Other organism", cached_reference, False)
    register.assert_not_called()


@pytest.mark.parametrize(
    "cached_reference, expected_error",
    [(None, MLSTUnavailableError), ("incomplete", ReferenceSetupError)],
)
def test_prepare_mlst_profile_offline_never_downloads(referencer, cached_reference, expected_error):
    with (
        patch.object(referencer.db_access, "profiles", {}),
        patch.object(referencer, "add_pubmlst") as download,
        pytest.raises(expected_error, match="downloads are disabled"),
    ):
        referencer._prepare_mlst_profile("Other organism", cached_reference, True)
    download.assert_not_called()


def test_prepare_mlst_profile_rejects_unexpected_empty_download_result(referencer):
    with (
        patch.object(referencer, "add_pubmlst", return_value=None),
        patch.object(referencer, "create_new_profile_tables") as register,
        pytest.raises(ReferenceSetupError, match="could not retrieve MLST resources"),
    ):
        referencer._prepare_mlst_profile("Other organism", None, False)
    register.assert_not_called()


def test_prepare_mlst_profile_propagates_download_failures(referencer):
    with (
        patch.object(referencer, "add_pubmlst", side_effect=ReferenceSetupError("network failure")),
        pytest.raises(ReferenceSetupError, match="network failure"),
    ):
        referencer._prepare_mlst_profile("Other organism", None, False)


@pytest.mark.parametrize("skip_update", [False, True])
def test_validate_mlst_reference_checks_sequences_and_indexes(
    referencer, mlst_reference_files, skip_update
):
    with patch.object(referencer, "_prepare_blast_indexes") as indexes:
        referencer._validate_mlst_reference(mlst_reference_files, skip_update)
    indexes.assert_called_once_with(
        Path(referencer.folders.references) / mlst_reference_files, ".tfa", skip_update
    )


@pytest.mark.parametrize("missing", ["locus", "profile"])
def test_validate_mlst_reference_requires_loci_and_profile(
    referencer, mlst_reference_files, missing
):
    if missing == "locus":
        path = Path(referencer.folders.references) / mlst_reference_files / "locus.tfa"
    else:
        path = Path(referencer.folders.profiles) / mlst_reference_files
    path.unlink()
    with patch.object(referencer, "_prepare_blast_indexes") as indexes:
        with pytest.raises(ReferenceSetupError, match="no MLST loci|missing MLST profile file"):
            referencer._validate_mlst_reference(mlst_reference_files, False)
    indexes.assert_not_called()


def test_validate_mlst_reference_rejects_corrupt_loci_before_indexing(
    referencer, mlst_reference_files
):
    path = Path(referencer.folders.references) / mlst_reference_files / "locus.tfa"
    path.write_text("Error: unavailable")
    with patch.object(referencer, "_prepare_blast_indexes") as indexes:
        with pytest.raises(ReferenceSetupError, match="invalid nucleotide FASTA"):
            referencer._validate_mlst_reference(mlst_reference_files, False)
    indexes.assert_not_called()


def test_prepare_mlst_references_separates_warnings_errors_and_successes(referencer):
    with (
        patch.object(referencer, "organism2reference", side_effect=[None, "broken", "valid"]),
        patch.object(
            referencer,
            "_prepare_mlst_profile",
            side_effect=[
                MLSTUnavailableError("No scheme"),
                ReferenceSetupError("download failed"),
                "valid",
            ],
        ),
        patch.object(referencer, "_validate_mlst_reference") as validate,
        patch.object(referencer.db_access.session, "rollback") as rollback,
        patch.object(referencer.logger, "warning") as warning,
    ):
        needs_expec, errors = referencer._prepare_mlst_references(
            {"Unknown": ["A1", "A2"], "Broken": ["A3"], "Valid": ["A4"]}, False
        )
    assert not needs_expec
    assert errors == ["Samples A3, organism Broken: download failed"]
    assert "A1, A2" in warning.call_args.args[0]
    assert "skipping MLST" in warning.call_args.args[0]
    rollback.assert_called_once_with()
    validate.assert_called_once_with("valid", False)


@pytest.mark.parametrize("cached_reference", [None, "escherichia_coli"])
def test_prepare_mlst_references_keeps_expec_requirement_when_mlst_validation_fails(
    referencer, cached_reference
):
    with (
        patch.object(referencer, "organism2reference", return_value=cached_reference),
        patch.object(referencer, "_prepare_mlst_profile", return_value="escherichia_coli"),
        patch.object(
            referencer, "_validate_mlst_reference", side_effect=ReferenceSetupError("bad loci")
        ),
    ):
        needs_expec, errors = referencer._prepare_mlst_references(
            {"Escherichia coli": ["A1"]}, False
        )
    assert needs_expec
    assert errors == ["Samples A1, organism Escherichia coli: bad loci"]


def test_prepare_mlst_references_keeps_cached_expec_requirement_when_profile_fails(referencer):
    with (
        patch.object(referencer, "organism2reference", return_value="escherichia_coli"),
        patch.object(
            referencer, "_prepare_mlst_profile", side_effect=ReferenceSetupError("bad profile")
        ),
    ):
        needs_expec, errors = referencer._prepare_mlst_references(
            {"Escherichia coli": ["A1"]}, False
        )
    assert needs_expec
    assert errors == ["Samples A1, organism Escherichia coli: bad profile"]


@pytest.mark.parametrize("suffix", [".fsa", ".tfa"])
@pytest.mark.parametrize("skip_update", [False, True])
def test_prepare_blast_reference_validates_matching_files_and_prepares_indexes(
    referencer, tmp_path, suffix, skip_update
):
    directory = tmp_path / "sequences"
    directory.mkdir()
    (directory / f"gene{suffix}").write_text(">gene\nACGT\n")
    (directory / "notes.txt").write_text("not a reference")
    with patch.object(referencer, "_prepare_blast_indexes") as indexes:
        referencer._prepare_blast_reference(directory, suffix, skip_update)
    indexes.assert_called_once_with(directory, suffix, skip_update)


def test_prepare_blast_reference_requires_sequences(referencer, tmp_path):
    with pytest.raises(ReferenceSetupError, match="no reference sequences"):
        referencer._prepare_blast_reference(tmp_path / "missing", ".fsa", True)


def test_prepare_blast_reference_rejects_invalid_fasta_before_indexing(referencer, tmp_path):
    (tmp_path / "gene.fsa").write_text("Error: unavailable")
    with patch.object(referencer, "_prepare_blast_indexes") as indexes:
        with pytest.raises(ReferenceSetupError, match="invalid nucleotide FASTA"):
            referencer._prepare_blast_reference(tmp_path, ".fsa", False)
    indexes.assert_not_called()


def test_prepare_blast_reference_propagates_index_failures(referencer, tmp_path):
    (tmp_path / "gene.fsa").write_text(">gene\nACGT\n")
    with patch.object(
        referencer, "_prepare_blast_indexes", side_effect=ReferenceSetupError("index failed")
    ):
        with pytest.raises(ReferenceSetupError, match="index failed"):
            referencer._prepare_blast_reference(tmp_path, ".fsa", False)


@pytest.mark.parametrize("needs_expec", [False, True])
@pytest.mark.parametrize("skip_update", [False, True])
def test_prepare_typing_references_selects_targets_and_collects_all_failures(
    referencer, tmp_path, needs_expec, skip_update
):
    resistance = tmp_path / "resistance"
    expec = tmp_path / "expec" / "ExPEC.fsa"
    referencer.folders.resistances = str(resistance)
    referencer.folders.expec = str(expec)
    expected = [call(resistance, ".fsa", skip_update)]
    if needs_expec:
        expected.append(call(expec.parent, ".fsa", skip_update))
    with patch.object(
        referencer, "_prepare_blast_reference", side_effect=ReferenceSetupError("invalid resource")
    ) as prepare:
        errors = referencer._prepare_typing_references(needs_expec, skip_update)
    assert prepare.call_args_list == expected
    assert len(errors) == len(expected)
    assert str(resistance) in errors[0]
    if needs_expec:
        assert str(expec.parent) in errors[1]


@pytest.mark.parametrize("qc_only", [False, True])
def test_prepare_references_orchestrates_resource_groups(referencer, qc_only):
    samples = [{"CG_ID_sample": "A1", "reference": "NC_1", "organism": "Other"}]
    with (
        patch.object(
            referencer, "_prepare_genome_references", return_value=["genome error"]
        ) as genomes,
        patch.object(
            referencer, "_prepare_mlst_references", return_value=(True, ["MLST error"])
        ) as mlst,
        patch.object(
            referencer, "_prepare_typing_references", return_value=["typing error"]
        ) as typing,
    ):
        errors = referencer.prepare_references(samples, qc_only=qc_only, skip_update=True)
    genomes.assert_called_once_with({"NC_1": ["A1"]}, True)
    if qc_only:
        assert errors == ["genome error"]
        mlst.assert_not_called()
        typing.assert_not_called()
    else:
        assert errors == ["genome error", "MLST error", "typing error"]
        mlst.assert_called_once_with({"Other": ["A1"]}, True)
        typing.assert_called_once_with(True, True)
