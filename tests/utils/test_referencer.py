import io
import gzip
import json
import os
import subprocess
from urllib.error import HTTPError
from pathlib import Path
from unittest.mock import Mock, patch

import pytest

from microSALT.exc.exceptions import MLSTUnavailableError, ReferenceSetupError
from microSALT.utils.pubmlst.exceptions import PubMLSTError
from microSALT.utils.referencer import Referencer


@pytest.fixture
def referencer(config, logger, tmp_path):
    folders = config.folders.model_copy(update={"genomes": str(tmp_path)})
    return Referencer(
        log=logger,
        folders=folders,
        threshold=config.threshold,
        pubmlst=config.pubmlst,
        pasteur=config.pasteur,
        singularity=config.singularity,
        containers=config.containers,
    )


def write_indexes(path):
    for suffix in (".amb", ".ann", ".bwt", ".pac", ".sa", ".fai"):
        Path(f"{path}{suffix}").write_text("index")


@pytest.mark.parametrize(
    "body",
    [
        "Error: F a i l e d  t o  u n d e r s t a n d  i d : GCF_000262545.1",
        "",
        "<html>Service unavailable</html>",
        ">empty\n",
        ">one\nACGT\n>empty\n",
        ">\nACGT\n",
        ">one\nACGT!\n",
        "ACGT\n",
    ],
)
def test_efetch_invalid_body_is_not_published_or_indexed(referencer, body):
    with (
        patch("microSALT.utils.referencer.Entrez.efetch", return_value=io.StringIO(body)),
        patch.object(referencer, "_index_genome") as index,
        pytest.raises(ReferenceSetupError),
    ):
        referencer.download_ncbi("NC_123")
    index.assert_not_called()
    assert not list(Path(referencer.folders.genomes).iterdir())


@pytest.fixture
def assembly_api():
    payload = {
        "search": {"esearchresult": {"count": "1", "idlist": ["495868"]}},
        "summary": {
            "result": {
                "495868": {
                    "assemblyaccession": "GCF_000262545.1",
                    "synonym": {"refseq": "GCF_000262545.1", "genbank": "GCA_000262545.1"},
                    "ftppath_refseq": "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/262/545/GCF_000262545.1_ASM26254v1",
                    "ftppath_genbank": "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/262/545/GCA_000262545.1_ASM26254v1",
                }
            }
        },
    }
    with (
        patch(
            "microSALT.utils.referencer.Entrez.esearch",
            side_effect=lambda **kwargs: io.StringIO(json.dumps(payload["search"])),
        ) as search,
        patch(
            "microSALT.utils.referencer.Entrez.esummary",
            side_effect=lambda **kwargs: io.StringIO(json.dumps(payload["summary"])),
        ) as summary,
    ):
        yield payload, search, summary


@pytest.mark.parametrize("prefix", ["GCF", "GCA"])
def test_assembly_download_uses_matching_database_and_complete_fasta(
    referencer, assembly_api, prefix
):
    _, search, summary = assembly_api
    reference = f"{prefix}_000262545.1"
    fasta = b">chromosome\nACGTN\n>plasmid\nTTAA\n"
    with (
        patch(
            "microSALT.utils.referencer.urllib.request.urlopen",
            return_value=io.BytesIO(gzip.compress(fasta)),
        ) as download,
        patch("microSALT.utils.referencer.Entrez.efetch") as nucleotide,
        patch.object(referencer, "_index_genome", side_effect=write_indexes),
    ):
        referencer.download_ncbi(f" {reference}\n")
    search.assert_called_once_with(
        db="assembly", term=f"{reference}[Assembly Accession]", retmax=2, retmode="json"
    )
    summary.assert_called_once_with(db="assembly", id="495868", retmode="json")
    url = f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{prefix}/000/262/545/{reference}_ASM26254v1/{reference}_ASM26254v1_genomic.fna.gz"
    download.assert_called_once_with(url, timeout=60)
    nucleotide.assert_not_called()
    path = Path(referencer.folders.genomes) / f"{reference}.fasta"
    assert path.read_bytes() == fasta
    assert not referencer._missing_genome_indexes(path)


@pytest.mark.parametrize("count, identifiers", [("0", []), ("2", ["1", "2"])])
def test_missing_or_ambiguous_assembly_is_a_setup_error(
    referencer, assembly_api, count, identifiers
):
    payload, _, summary = assembly_api
    payload["search"] = {"esearchresult": {"count": count, "idlist": identifiers}}
    with (
        patch("microSALT.utils.referencer.urllib.request.urlopen") as download,
        patch.object(referencer, "_index_genome") as index,
        pytest.raises(ReferenceSetupError, match="exactly one NCBI Assembly record"),
    ):
        referencer.download_ncbi("GCF_000262545.1")
    summary.assert_not_called()
    download.assert_not_called()
    index.assert_not_called()
    assert not list(Path(referencer.folders.genomes).iterdir())


def test_requested_assembly_version_is_not_substituted(referencer, assembly_api):
    payload, _, _ = assembly_api
    payload["summary"]["result"]["495868"]["synonym"]["refseq"] = "GCF_000262545.2"
    with patch("microSALT.utils.referencer.urllib.request.urlopen") as download:
        with pytest.raises(
            ReferenceSetupError, match="instead of requested accession GCF_000262545.1"
        ):
            referencer.download_ncbi("GCF_000262545.1")
    download.assert_not_called()


def test_missing_assembly_download_does_not_fall_back_to_paired_assembly(referencer, assembly_api):
    payload, _, _ = assembly_api
    payload["summary"]["result"]["495868"]["ftppath_refseq"] = ""
    with patch("microSALT.utils.referencer.Entrez.efetch") as nucleotide:
        with pytest.raises(ReferenceSetupError, match="No refseq genomic FASTA download"):
            referencer.download_ncbi("GCF_000262545.1")
    nucleotide.assert_not_called()


@pytest.mark.parametrize(
    "location",
    [
        "https://example.com/genomes/all/GCF_000262545.1_ASM26254v1",
        "file:///genomes/all/GCF_000262545.1_ASM26254v1",
        "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000262545.2_ASM26254v2",
    ],
)
def test_unexpected_assembly_download_location_is_rejected(referencer, assembly_api, location):
    payload, _, _ = assembly_api
    payload["summary"]["result"]["495868"]["ftppath_refseq"] = location
    with patch("microSALT.utils.referencer.urllib.request.urlopen") as download:
        with pytest.raises(ReferenceSetupError, match="Unexpected NCBI Assembly download path"):
            referencer.download_ncbi("GCF_000262545.1")
    download.assert_not_called()


@pytest.mark.parametrize(
    "content",
    [
        b"Error: download unavailable",
        gzip.compress(b"Error: no genome"),
        gzip.compress(b">empty\n"),
        gzip.compress(b">one\nACGT\n")[:-5],
    ],
)
def test_invalid_assembly_download_is_not_published(referencer, assembly_api, content):
    with (
        patch(
            "microSALT.utils.referencer.urllib.request.urlopen", return_value=io.BytesIO(content)
        ),
        patch.object(referencer, "_index_genome") as index,
        pytest.raises(ReferenceSetupError),
    ):
        referencer.download_ncbi("GCF_000262545.1")
    index.assert_not_called()
    assert not list(Path(referencer.folders.genomes).iterdir())


def test_assembly_http_failure_is_reported_without_publishing(referencer, assembly_api):
    with (
        patch(
            "microSALT.utils.referencer.urllib.request.urlopen",
            side_effect=HTTPError("https://ftp.ncbi.nlm.nih.gov", 404, "Not Found", {}, None),
        ),
        pytest.raises(ReferenceSetupError, match="GCF_000262545.1.*HTTP Error 404"),
    ):
        referencer.download_ncbi("GCF_000262545.1")
    assert not list(Path(referencer.folders.genomes).iterdir())


def test_unversioned_assembly_can_resolve_to_one_record(referencer, assembly_api):
    assert "GCF_000262545.1_ASM26254v1_genomic.fna.gz" in referencer._assembly_fasta_url(
        "GCF_000262545"
    )


def test_cached_assembly_is_used_without_network_access(referencer):
    reference = "GCF_000262545.1"
    path = Path(referencer.folders.genomes) / f"{reference}.fasta"
    path.write_text(">chromosome\nACGT\n")
    write_indexes(path)
    with (
        patch("microSALT.utils.referencer.Entrez.esearch") as search,
        patch("microSALT.utils.referencer.Entrez.efetch") as fetch,
    ):
        errors = referencer.prepare_references(
            [
                {"CG_ID_sample": "A1", "reference": reference},
            ],
            qc_only=True,
            skip_update=True,
        )
    assert not errors
    search.assert_not_called()
    fetch.assert_not_called()


def test_valid_download_is_indexed_and_published(referencer):
    with (
        patch(
            "microSALT.utils.referencer.Entrez.efetch",
            return_value=io.StringIO(">one\nACGTN\n>two\nryswkmbdhv\n"),
        ) as fetch,
        patch.object(referencer, "_index_genome", side_effect=write_indexes),
    ):
        referencer.download_ncbi(" NC_123\n")
    assert fetch.call_args.kwargs["id"] == "NC_123"
    assert fetch.call_args.kwargs["retmode"] == "text"
    path = Path(referencer.folders.genomes) / "NC_123.fasta"
    referencer.validate_fasta(path)
    assert not referencer._missing_genome_indexes(path)


def test_failed_download_preserves_existing_reference(referencer):
    path = Path(referencer.folders.genomes) / "NC_123.fasta"
    path.write_text(">original\nACGT\n")
    write_indexes(path)
    original = {file.name: file.read_bytes() for file in path.parent.iterdir()}
    with (
        patch("microSALT.utils.referencer.Entrez.efetch", return_value=io.StringIO(">new\nAAAA\n")),
        patch.object(referencer, "_index_genome", side_effect=ReferenceSetupError("index failed")),
        pytest.raises(ReferenceSetupError, match="index failed"),
    ):
        referencer.download_ncbi("NC_123")
    assert {file.name: file.read_bytes() for file in path.parent.iterdir()} == original


@pytest.mark.parametrize("skip_update", [True, False])
def test_cached_corrupt_reference_is_rejected_without_overwriting(referencer, skip_update):
    path = Path(referencer.folders.genomes) / "NC_123.fasta"
    path.write_text("Error: invalid id")
    with patch.object(referencer, "download_ncbi") as download:
        errors = referencer.prepare_references(
            [
                {"CG_ID_sample": "A1", "reference": "NC_123"},
                {"CG_ID_sample": "A2", "reference": "NC_123"},
            ],
            qc_only=True,
            skip_update=skip_update,
        )
    assert len(errors) == 1
    assert "A1, A2" in errors[0] and "invalid nucleotide FASTA" in errors[0]
    download.assert_not_called()
    assert path.read_text() == "Error: invalid id"


def test_collect_all_download_failures_once_per_reference(referencer):
    with patch.object(
        referencer, "download_ncbi", side_effect=ReferenceSetupError("fetch failed")
    ) as download:
        errors = referencer.prepare_references(
            [
                {"CG_ID_sample": "A1", "reference": "NC_1"},
                {"CG_ID_sample": "A2", "reference": "NC_2"},
                {"CG_ID_sample": "A3", "reference": "NC_1"},
            ],
            qc_only=True,
        )
    assert len(errors) == download.call_count == 2
    assert "A1, A3" in errors[0]


def test_skip_update_validates_without_downloading_or_indexing(referencer):
    path = Path(referencer.folders.genomes) / "NC_123.fasta"
    path.write_text(">one\nACGT\n")
    with (
        patch.object(referencer, "download_ncbi") as download,
        patch.object(referencer, "_index_genome") as index,
    ):
        errors = referencer.prepare_references(
            [{"CG_ID_sample": "A1", "reference": "NC_123"}], qc_only=True, skip_update=True
        )
    assert "missing or stale genome indexes" in errors[0]
    download.assert_not_called()
    index.assert_not_called()


def test_qc_only_does_not_require_mlst_or_resistance(referencer):
    path = Path(referencer.folders.genomes) / "NC_123.fasta"
    path.write_text(">one\nACGT\n")
    write_indexes(path)
    with patch.object(referencer, "add_pubmlst") as mlst:
        errors = referencer.prepare_references(
            [
                {"CG_ID_sample": "A1", "reference": "NC_123", "organism": "other"},
            ],
            qc_only=True,
            skip_update=True,
        )
    assert not errors
    mlst.assert_not_called()


def test_stale_genome_indexes_are_rebuilt(referencer):
    path = Path(referencer.folders.genomes) / "NC_123.fasta"
    path.write_text(">one\nACGT\n")
    write_indexes(path)
    for index in path.parent.glob("*.fasta.*"):
        os.utime(index, (1, 1))
    with patch.object(referencer, "_index_genome", side_effect=write_indexes) as index:
        errors = referencer.prepare_references(
            [{"CG_ID_sample": "A1", "reference": "NC_123"}], qc_only=True
        )
    assert not errors
    index.assert_called_once_with(path)


def test_nonzero_index_command_includes_stderr(referencer):
    with patch(
        "microSALT.utils.referencer.subprocess.run",
        side_effect=subprocess.CalledProcessError(
            2,
            ["bwa"],
            stderr="cannot build index",
        ),
    ):
        with pytest.raises(ReferenceSetupError, match="exit 2.*cannot build index"):
            referencer._run_index("bwa index reference.fasta", referencer.folders.genomes)


def test_success_without_genome_indexes_is_rejected(referencer):
    path = Path(referencer.folders.genomes) / "NC_123.fasta"
    path.write_text(">one\nACGT\n")
    with patch("microSALT.utils.referencer.subprocess.run") as run:
        with pytest.raises(ReferenceSetupError, match="did not produce usable files"):
            referencer._index_genome(path)
    assert run.call_count == 2


def test_blast_index_failure_is_not_swallowed(referencer, tmp_path):
    (tmp_path / "locus.tfa").write_text(">one\nACGT\n")
    with patch.object(referencer, "_run_index", side_effect=ReferenceSetupError("blast failed")):
        with pytest.raises(ReferenceSetupError, match="blast failed"):
            referencer.index_db(str(tmp_path), ".tfa")


def test_full_analysis_collects_mlst_and_resistance_failures(referencer, tmp_path):
    referencer.folders.resistances = str(tmp_path / "missing-resistance")
    with patch.object(
        referencer, "add_pubmlst", side_effect=ReferenceSetupError("download failed")
    ) as download:
        errors = referencer.prepare_references(
            [
                {"CG_ID_sample": "A1", "organism": "Unknown one"},
                {"CG_ID_sample": "A2", "organism": "Unknown two"},
            ]
        )
    assert download.call_count == 2
    assert len(errors) == 3
    assert "A1" in errors[0] and "A2" in errors[1]
    assert "missing-resistance" in errors[2]


def test_full_analysis_without_published_mlst_scheme_warns_and_continues(referencer, tmp_path):
    resistance = tmp_path / "resistance"
    resistance.mkdir()
    (resistance / "gene.fsa").write_text(">gene\nACGT\n")
    referencer.folders.resistances = str(resistance)
    with (
        patch.object(
            referencer, "add_pubmlst", side_effect=MLSTUnavailableError("No published scheme")
        ),
        patch.object(referencer, "_prepare_blast_indexes") as indexes,
        patch.object(referencer.logger, "warning") as warning,
    ):
        errors = referencer.prepare_references(
            [{"CG_ID_sample": "A1", "organism": "Unknown organism"}]
        )
    assert not errors
    indexes.assert_called_once_with(resistance, ".fsa", False)
    assert "skipping MLST" in warning.call_args.args[0]


def test_skip_update_without_cached_mlst_keeps_other_analyses(referencer, tmp_path):
    resistance = tmp_path / "resistance"
    resistance.mkdir()
    (resistance / "gene.fsa").write_text(">gene\nACGT\n")
    referencer.folders.resistances = str(resistance)
    with (
        patch.object(referencer, "add_pubmlst") as download,
        patch.object(referencer, "_prepare_blast_indexes"),
    ):
        errors = referencer.prepare_references(
            [{"CG_ID_sample": "A1", "organism": "Unknown organism"}], skip_update=True
        )
    assert not errors
    download.assert_not_called()


def test_absent_organism_is_distinct_from_failed_lookup(referencer):
    with patch.object(referencer, "query_pubmlst", return_value=[]):
        with pytest.raises(MLSTUnavailableError, match="No pubMLST database"):
            referencer.add_pubmlst("Unknown organism")
    with patch.object(referencer, "query_pubmlst", side_effect=PubMLSTError("connection failed")):
        with pytest.raises(ReferenceSetupError, match="connection failed"):
            referencer.add_pubmlst("Unknown organism")


@pytest.fixture
def mlst_client(referencer, tmp_path):
    referencer.folders.profiles = str(tmp_path / "profiles")
    referencer.folders.references = str(tmp_path / "loci")
    Path(referencer.folders.profiles).mkdir()
    Path(referencer.folders.references).mkdir()
    client = Mock()
    client.parse_url.return_value = {"db": "example_seqdef", "scheme_id": 1}
    client.retrieve_scheme_info.return_value = {
        "description": "MLST",
        "loci": ["https://example.test/loci/locus"],
    }
    client.download_profiles_csv.return_value = "ST\tlocus\n1\t1\n"
    client.download_locus.return_value = ">locus_1\nACGT\n"
    referencer.client = client
    return client


def test_available_mlst_without_version_date_can_be_downloaded(referencer, mlst_client):
    with (
        patch.object(referencer.db_access, "read_version", return_value="0"),
        patch.object(referencer, "index_db"),
    ):
        version = referencer.download_pubmlst(
            "Example species", "https://example.test/db/example_seqdef"
        )
    assert version is None
    assert (Path(referencer.folders.profiles) / "example_species").is_file()
    assert (Path(referencer.folders.references) / "example_species/locus.tfa").is_file()


def test_add_pubmlst_accepts_successful_download_without_version_date(referencer):
    catalog = [
        {
            "databases": [
                {
                    "description": "Example species sequence definitions",
                    "href": "https://example.test/db/example_seqdef",
                }
            ]
        }
    ]
    with (
        patch.object(referencer, "query_pubmlst", return_value=catalog),
        patch.object(referencer, "download_pubmlst", return_value=None),
    ):
        assert referencer.add_pubmlst("Example species") == "example_species"


def test_no_mlst_scheme_is_not_a_download_failure(referencer, mlst_client):
    mlst_client.retrieve_scheme_info.return_value = {"description": "Other scheme"}
    mlst_client.list_schemes.return_value = {"schemes": []}
    with patch.object(referencer.db_access, "read_version", return_value="0"):
        with pytest.raises(MLSTUnavailableError, match="No MLST scheme is available"):
            referencer.download_pubmlst("Example species", "https://example.test/db/example_seqdef")
    mlst_client.download_profiles_csv.assert_not_called()


def test_mlst_version_lookup_connection_failure_remains_fatal(referencer, mlst_client):
    mlst_client.retrieve_scheme_info.side_effect = PubMLSTError("connection failed")
    with pytest.raises(ReferenceSetupError, match="connection failed"):
        referencer.download_pubmlst("Example species", "https://example.test/db/example_seqdef")


@pytest.mark.parametrize("force", [False, True])
def test_cached_mlst_without_external_version_never_compares_or_stores_none(referencer, force):
    referencer.organisms = ["example_species"]
    url = "https://example.test/db/example_seqdef"
    catalog = [
        {"databases": [{"description": "Example species sequence definitions", "href": url}]}
    ]
    with (
        patch.object(referencer, "query_pubmlst", return_value=catalog),
        patch.object(referencer, "external_version", return_value=None),
        patch.object(referencer.db_access, "read_version", return_value="2025-01-01"),
        patch.object(referencer, "download_pubmlst") as download,
        patch.object(referencer.db_access, "update_version") as update,
        patch.object(referencer.db_access, "refresh_profiletable"),
    ):
        referencer.fetch_pubmlst(force=force)
    if force:
        download.assert_called_once_with("example_species", url, True)
    else:
        download.assert_not_called()
    update.assert_not_called()


def test_full_analysis_accepts_valid_cached_typing_resources(referencer, tmp_path):
    references = tmp_path / "loci"
    profiles = tmp_path / "profiles"
    resistance = tmp_path / "resistance"
    for path in (references / "staphylococcus_aureus", profiles, resistance):
        path.mkdir(parents=True)
    (profiles / "staphylococcus_aureus").write_text("ST\tlocus\n1\t1\n")
    for path in (references / "staphylococcus_aureus" / "locus.tfa", resistance / "gene.fsa"):
        path.write_text(">one\nACGT\n")
        for suffix in (".nhr", ".nin", ".nsq"):
            Path(f"{path.with_suffix('')}{suffix}").write_text("index")
    referencer.folders.references = str(references)
    referencer.folders.profiles = str(profiles)
    referencer.folders.resistances = str(resistance)
    with (
        patch.object(referencer.db_access, "profiles", {"staphylococcus_aureus": object()}),
        patch.object(referencer, "add_pubmlst") as download,
        patch.object(referencer, "index_db") as index,
    ):
        errors = referencer.prepare_references(
            [
                {"CG_ID_sample": "A1", "organism": "Staphylococcus aureus"},
            ],
            skip_update=True,
        )
    assert not errors
    download.assert_not_called()
    index.assert_not_called()


def test_blast_success_without_indexes_is_rejected(referencer, tmp_path):
    (tmp_path / "locus.tfa").write_text(">one\nACGT\n")
    with patch.object(referencer, "_run_index"):
        with pytest.raises(ReferenceSetupError, match="did not produce usable BLAST indexes"):
            referencer.index_db(str(tmp_path), ".tfa")


def test_missing_reference_is_an_error_with_skip_update(referencer):
    with patch.object(referencer, "download_ncbi") as download:
        errors = referencer.prepare_references(
            [
                {"CG_ID_sample": "A1", "reference": "NC_missing"},
            ],
            qc_only=True,
            skip_update=True,
        )
    assert "missing reference" in errors[0]
    download.assert_not_called()
