"""Validate analysis inputs without submitting jobs or modifying source metadata."""

import gzip
import re
from datetime import datetime
from pathlib import Path


def verify_fastq(indir: str, file_pattern: str) -> list[str]:
    """Check every matching mate and stream FASTQ records, including gzip integrity."""
    directory = Path(indir)
    files = sorted(
        path.name for path in directory.iterdir() if re.fullmatch(file_pattern, path.name)
    )
    if not files:
        raise ValueError(f"No FASTQ files in {directory} match file_pattern '{file_pattern}'")
    errors = []
    for name in files:
        match = re.fullmatch(file_pattern, name)
        mate = {"1": "2", "2": "1", "forward": "reverse", "reverse": "forward"}.get(match[1])
        if mate is None:
            errors.append(f"{name}: unsupported mate identifier {match[1]!r}")
        else:
            partner = name[: match.start(1)] + mate + name[match.end(1) :]
            if partner not in files:
                errors.append(f"{name}: missing FASTQ mate {partner}")
        path = directory / name
        try:
            with gzip.open(path, "rt") as stream:
                count = 0
                while header := stream.readline():
                    sequence = stream.readline().rstrip("\r\n")
                    separator = stream.readline()
                    quality = stream.readline().rstrip("\r\n")
                    count += 1
                    if (
                        not header.startswith("@")
                        or not separator.startswith("+")
                        or not sequence
                        or len(sequence) != len(quality)
                    ):
                        raise ValueError(f"invalid FASTQ record {count}")
                if not count:
                    raise ValueError("empty FASTQ file")
        except (OSError, EOFError, UnicodeError, ValueError) as error:
            errors.append(f"{path}: {error}")
    if errors:
        raise ValueError("; ".join(errors))
    return files


def validate_sampleinfo(data, required_fields, input_dir: str, file_pattern: str):
    """Return normalized samples and all independently detectable metadata/input errors."""
    errors = []
    samples = [data] if isinstance(data, dict) else data
    if not isinstance(samples, list) or not samples:
        return [], ["Sample info must contain a sample object or a non-empty list of samples"]
    normalized = []
    seen = set()
    projects = set()
    root = Path(input_dir)
    if not input_dir or not root.is_dir():
        errors.append(f"Sequence data folder {root} does not exist")
    for index, entry in enumerate(samples, 1):
        if not isinstance(entry, dict):
            errors.append(f"Sample info entry {index} must be an object")
            continue
        sample = dict(entry)
        if isinstance(sample.get("reference"), str):
            sample["reference"] = "".join(sample["reference"].split())
        normalized.append(sample)
        name = sample.get("CG_ID_sample") or f"entry {index}"
        for field in required_fields:
            if not isinstance(sample.get(field), str) or not sample[field].strip():
                errors.append(f"Sample {name}: missing or invalid {field}")
        for field in ("CG_ID_sample", "CG_ID_project", "reference"):
            value = sample.get(field)
            if (
                isinstance(value, str)
                and value
                and (not re.fullmatch(r"[\w.-]+", value) or value in {".", ".."})
            ):
                errors.append(f"Sample {name}: unsafe {field} {value!r}")
        reference = sample.get("reference")
        if isinstance(reference, str) and reference.lower() in {"none", "null"}:
            errors.append(f"Sample {name}: a reference genome is required for alignment")
        for field in ("date_arrival", "date_libprep", "date_sequencing"):
            if sample.get(field):
                try:
                    datetime.strptime(sample[field], "%Y-%m-%d %H:%M:%S")
                except (TypeError, ValueError):
                    errors.append(f"Sample {name}: invalid {field}; expected YYYY-MM-DD HH:MM:SS")
        sample_id = sample.get("CG_ID_sample")
        if isinstance(sample_id, str):
            if sample_id in seen:
                errors.append(f"Duplicate sample ID {sample_id}")
            seen.add(sample_id)
        project = sample.get("CG_ID_project")
        if isinstance(project, str):
            projects.add(project)
        if (
            input_dir
            and root.is_dir()
            and isinstance(sample_id, str)
            and re.fullmatch(r"[\w.-]+", sample_id)
            and sample_id not in {".", ".."}
        ):
            sample_dir = root / sample_id if len(samples) > 1 else root
            try:
                verify_fastq(str(sample_dir), file_pattern)
            except Exception as error:
                errors.append(f"Sample {name}: {error}")
    if len(projects) > 1:
        errors.append("Sample info contains mixed projects")
    if input_dir and len(samples) > 1 and root.is_dir():
        try:
            for path in sorted(root.iterdir()):
                if path.is_dir() and path.name not in seen:
                    errors.append(f"Input folder {path}: no counterpart in sample info")
        except OSError as error:
            errors.append(f"Input folder {root}: {error}")
    return normalized, errors
