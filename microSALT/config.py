import json
from importlib.resources import files as resource_files
from pathlib import Path

from pydantic import BaseModel


class SlurmHeader(BaseModel):
    time: str
    threads: str
    qos: str
    job_prefix: str
    project: str
    type: str


class Regex(BaseModel):
    mail_recipient: str
    file_pattern: str
    verified_organisms: list[str] = []


class Folders(BaseModel):
    results: str
    reports: str
    log_file: str
    seqdata: str
    profiles: str
    references: str
    resistances: str
    genomes: str
    credentials: str
    adapters: str
    expec: str = ""  # filled in after construction


class Database(BaseModel):
    SQLALCHEMY_DATABASE_URI: str
    SQLALCHEMY_TRACK_MODIFICATIONS: str = "False"
    DEBUG: str = "True"


class Threshold(BaseModel):
    mlst_id: float = 100
    mlst_novel_id: float = 99.5
    mlst_span: float = 90
    motif_id: float = 97
    motif_span: float = 90
    total_reads_warn: float = 75
    total_reads_fail: float = 70
    NTC_total_reads_warn: float = 10
    NTC_total_reads_fail: float = 20
    mapped_rate_warn: float = 50
    mapped_rate_fail: float = 30
    duplication_rate_warn: float = 20
    duplication_rate_fail: float = 80
    insert_size_warn: float = 140
    insert_size_fail: float = 100
    average_coverage_warn: float = 100
    average_coverage_fail: float = 10
    bp_10x_warn: float = 85
    bp_10x_fail: float = 75
    bp_30x_warn: float = 70
    bp_50x_warn: float = 50
    bp_100x_warn: float = 20


class BIGSdbCredentials(BaseModel):
    client_id: str = ""
    client_secret: str = ""

class PubMLSTCredentials(BIGSdbCredentials):
    pass

class PasteurCredentials(BIGSdbCredentials):
    pass


class Singularity(BaseModel):
    binary: str = "/usr/bin/singularity"
    bind_paths: list[str] = []
    trimmomatic_adapters: str = "/opt/conda/share/trimmomatic/adapters/"


class Containers(BaseModel):
    skesa: str = ""
    blast: str = ""
    bwa: str = ""
    samtools: str = ""
    picard: str = ""
    trimmomatic: str = ""
    quast: str = ""


class MicroSALTConfig(BaseModel):
    slurm_header: SlurmHeader
    regex: Regex
    folders: Folders
    database: Database
    threshold: Threshold
    pubmlst: PubMLSTCredentials = PubMLSTCredentials()
    pasteur: PasteurCredentials = PasteurCredentials()
    singularity: Singularity = Singularity()
    containers: Containers = Containers()
    # Runtime fields set by the CLI, not from the JSON file
    dry: bool = False
    config_path: str = ""


def _strip_comments(obj):
    if isinstance(obj, dict):
        return {k: _strip_comments(v) for k, v in obj.items() if k != "_comment"}
    return obj


def load_config(path: str) -> MicroSALTConfig:
    """Parse and validate a microSALT JSON config file.

    The expec reference path is derived from the package data and injected
    after parsing, so it does not need to be present in the JSON file.
    """
    with open(path) as f:
        data = json.load(f)
    data = _strip_comments(data)
    config = MicroSALTConfig(**data)
    config.folders.expec = str(
        resource_files("microSALT").joinpath("unique_references", "ExPEC.fsa")
    )
    config.config_path = str(Path(path).resolve())
    return config
