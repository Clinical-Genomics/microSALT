"""Validate analysis inputs without submitting jobs or modifying source metadata."""

import re
from collections.abc import Collection, Iterable
from datetime import datetime
from pathlib import Path


def find_fastq_files(directory: Path, file_pattern: str) -> list[str]:
    """Find matching input names, rejecting an empty or unmatched directory."""
    files = sorted(
        path.name for path in directory.iterdir() if re.fullmatch(file_pattern, path.name)
    )
    if not files:
        raise ValueError(f"No FASTQ files in {directory} match file_pattern '{file_pattern}'")
    return files


def validate_fastq_mate(name: str, files: Collection[str], file_pattern: str) -> None:
    """Require the corresponding forward/reverse filename for a matching input."""
    match = re.fullmatch(file_pattern, name)
    if match is None or not match.lastindex:
        raise ValueError(f"{name}: file_pattern must match and capture a mate identifier")
    mate = {"1": "2", "2": "1", "forward": "reverse", "reverse": "forward"}.get(match[1])
    if mate is None:
        raise ValueError(f"{name}: unsupported mate identifier {match[1]!r}")
    partner = name[: match.start(1)] + mate + name[match.end(1) :]
    if partner not in files:
        raise ValueError(f"{name}: missing FASTQ mate {partner}")


def validate_input_file(path: Path) -> None:
    """Check file type, symlink resolution and read access without reading contents."""
    if not path.is_file():
        raise ValueError("not a regular file or symlink target is missing")
    with path.open("rb"):
        pass


def verify_fastq(indir: str, file_pattern: str) -> list[str]:
    """Collect filename, mate and accessibility errors without decompressing inputs."""
    directory = Path(indir)
    files = find_fastq_files(directory, file_pattern)
    errors = []
    for name in files:
        try:
            validate_fastq_mate(name, files, file_pattern)
        except ValueError as error:
            errors.append(str(error))
        path = directory / name
        try:
            validate_input_file(path)
        except (OSError, ValueError) as error:
            errors.append(f"{path}: {error}")
    if errors:
        raise ValueError("; ".join(errors))
    return files


def normalize_sample(sample: dict) -> dict:
    """Copy metadata and remove reference whitespace without modifying the source."""
    normalized = dict(sample)
    if isinstance(normalized.get("reference"), str):
        normalized["reference"] = "".join(normalized["reference"].split())
    return normalized


def is_safe_identifier(value: object) -> bool:
    """Allow non-empty identifiers that cannot traverse paths or inject shell syntax."""
    return (
        isinstance(value, str)
        and bool(re.fullmatch(r"[\w.-]+", value))
        and value not in {".", ".."}
    )


def validate_required_fields(sample: dict, required_fields: Iterable[str]) -> list[str]:
    """Require a non-empty string for each expected metadata field."""
    return [
        f"missing or invalid {field}"
        for field in required_fields
        if not isinstance(sample.get(field), str) or not sample[field].strip()
    ]


def validate_sample_identifiers(sample: dict) -> list[str]:
    """Check path identifiers and reject placeholder alignment references."""
    errors = []
    for field in ("CG_ID_sample", "CG_ID_project", "reference"):
        value = sample.get(field)
        if isinstance(value, str) and value and not is_safe_identifier(value):
            errors.append(f"unsafe {field} {value!r}")
    reference = sample.get("reference")
    if isinstance(reference, str) and reference.lower() in {"none", "null"}:
        errors.append("a reference genome is required for alignment")
    return errors


def validate_sample_dates(sample: dict) -> list[str]:
    """Check supplied dates independently of missing-field validation."""
    errors = []
    for field in ("date_arrival", "date_libprep", "date_sequencing"):
        if sample.get(field):
            try:
                datetime.strptime(sample[field], "%Y-%m-%d %H:%M:%S")
            except (TypeError, ValueError):
                errors.append(f"invalid {field}; expected YYYY-MM-DD HH:MM:SS")
    return errors


def validate_sample_metadata(sample: dict, required_fields: Iterable[str]) -> list[str]:
    """Compose independent metadata checks; callers add the sample label."""
    return (
        validate_required_fields(sample, required_fields)
        + validate_sample_identifiers(sample)
        + validate_sample_dates(sample)
    )


def validate_sample_input(
    sample: dict, input_dir: str, file_pattern: str, multiple_samples: bool
) -> list[str]:
    """Check FASTQ inputs only when the root and sample ID are usable."""
    sample_id = sample.get("CG_ID_sample")
    root = Path(input_dir)
    if not input_dir or not root.is_dir() or not is_safe_identifier(sample_id):
        return []
    sample_dir = root / sample_id if multiple_samples else root
    try:
        verify_fastq(str(sample_dir), file_pattern)
    except Exception as error:
        return [str(error)]
    return []


def validate_sample_groups(samples: list[dict]) -> list[str]:
    """Reject duplicate sample IDs and mixed projects across the collection."""
    errors = []
    seen = set()
    projects = set()
    for sample in samples:
        sample_id = sample.get("CG_ID_sample")
        if isinstance(sample_id, str):
            if sample_id in seen:
                errors.append(f"Duplicate sample ID {sample_id}")
            seen.add(sample_id)
        project = sample.get("CG_ID_project")
        if isinstance(project, str):
            projects.add(project)
    if len(projects) > 1:
        errors.append("Sample info contains mixed projects")
    return errors


def validate_input_folders(input_dir: str, samples: list[dict]) -> list[str]:
    """Report project input directories absent from sample metadata."""
    root = Path(input_dir)
    if not input_dir or not root.is_dir():
        return []
    sample_ids = {
        sample["CG_ID_sample"] for sample in samples if isinstance(sample.get("CG_ID_sample"), str)
    }
    try:
        return [
            f"Input folder {path}: no counterpart in sample info"
            for path in sorted(root.iterdir())
            if path.is_dir() and path.name not in sample_ids
        ]
    except OSError as error:
        return [f"Input folder {root}: {error}"]


def validate_sampleinfo(
    data, required_fields: Iterable[str], input_dir: str, file_pattern: str
) -> tuple[list[dict], list[str]]:
    """Normalize samples and collect independent metadata, input and collection errors."""
    samples = [data] if isinstance(data, dict) else data
    if not isinstance(samples, list) or not samples:
        return [], ["Sample info must contain a sample object or a non-empty list of samples"]
    normalized = []
    errors = []
    required_fields = tuple(required_fields)
    if not input_dir or not Path(input_dir).is_dir():
        errors.append(f"Sequence data folder {Path(input_dir)} does not exist")
    multiple_samples = len(samples) > 1
    for index, entry in enumerate(samples, 1):
        if not isinstance(entry, dict):
            errors.append(f"Sample info entry {index} must be an object")
            continue
        sample = normalize_sample(entry)
        normalized.append(sample)
        name = sample.get("CG_ID_sample") or f"entry {index}"
        sample_errors = validate_sample_metadata(sample, required_fields)
        sample_errors.extend(
            validate_sample_input(sample, input_dir, file_pattern, multiple_samples)
        )
        errors.extend(f"Sample {name}: {error}" for error in sample_errors)
    errors.extend(validate_sample_groups(normalized))
    if multiple_samples:
        errors.extend(validate_input_folders(input_dir, normalized))
    return normalized, errors
