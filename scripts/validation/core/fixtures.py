"""Helpers for building modified inputs without touching real FASTQ, reference or database files.

Everything is created under <validation_root>/fixtures/<label>/ as copies or symlink farms.
"""

import json
import os
import re
import shutil
import sqlite3
from pathlib import Path

from microSALT.utils.referencer import Referencer

from .config import Case, ValidationConfig


class FixtureHelper:
    def __init__(self, cfg: ValidationConfig):
        self.cfg = cfg
        self.ms = cfg.microsalt
        self.pattern = self.ms.regex.file_pattern

    def fresh(self, label: str) -> Path:
        path = self.cfg.fixture_dir / label
        if path.exists():
            shutil.rmtree(path)
        path.mkdir(parents=True)
        return path

    @staticmethod
    def load_sampleinfo(case: Case) -> list[dict]:
        with open(case.sampleinfo) as stream:
            data = json.load(stream)
        return data if isinstance(data, list) else [data]

    @staticmethod
    def write_json(path: Path, data) -> Path:
        with open(path, "w") as stream:
            json.dump(data, stream, indent=2)
        return path

    def modified_sampleinfo(
        self, case: Case, label: str, index: int = 0, **changes
    ) -> tuple[Path, dict]:
        """Copy of the case sample info with fields of one sample changed (or deleted when None)."""
        data = self.load_sampleinfo(case)
        for key, value in changes.items():
            if value is None:
                data[index].pop(key, None)
            else:
                data[index][key] = value
        path = self.write_json(self.cfg.fixture_dir / f"{case.name}_{label}.json", data)
        return path, data[index]

    def subset_sampleinfo(self, case: Case, count: int = 1) -> tuple[Path, list[str]]:
        data = self.load_sampleinfo(case)[:count]
        path = self.write_json(
            self.cfg.fixture_dir / f"{case.name}_subset{count}.json", data
        )
        return path, [entry["CG_ID_sample"] for entry in data]

    def fastq_farm(self, root: Path, case: Case) -> Path:
        """Per-sample folders of symlinks to the case FASTQ files."""
        farm = root / "fastq"
        for sample_dir in sorted(p for p in case.fastq.iterdir() if p.is_dir()):
            target = farm / sample_dir.name
            target.mkdir(parents=True)
            for item in sample_dir.iterdir():
                os.symlink(item.resolve(), target / item.name)
        return farm

    def genomes_farm(self, root: Path) -> Path:
        """Symlinks to every file in the configured reference genome folder."""
        farm = root / "genomes"
        farm.mkdir()
        for item in Path(self.ms.folders.genomes).iterdir():
            if item.is_file():
                os.symlink(item.resolve(), farm / item.name)
        return farm

    def alt_config(
        self, root: Path, genomes: Path | None = None, db_uri: str | None = None
    ) -> Path:
        """Copy of the microSALT config with selected folders or database replaced."""
        with open(self.cfg.microsalt_config) as stream:
            data = json.load(stream)
        if genomes is not None:
            data["folders"]["genomes"] = str(genomes)
        if db_uri is not None:
            data["database"]["SQLALCHEMY_DATABASE_URI"] = db_uri
        return self.write_json(root / "microsalt_config.json", data)

    def scratch_sqlite(self, root: Path) -> tuple[Path, str] | None:
        """Copy of the sqlite database, or None when the configured database is not sqlite."""
        uri = self.ms.database.SQLALCHEMY_DATABASE_URI
        if not uri.startswith("sqlite:///"):
            return None
        scratch = root / "microsalt_scratch.db"
        shutil.copyfile(Path(uri.removeprefix("sqlite:///")), scratch)
        return scratch, f"sqlite:///{scratch}"

    @staticmethod
    def sqlite_execute(path: Path, statement: str) -> None:
        with sqlite3.connect(path) as conn:
            conn.execute(statement)

    def mates(self, sample_dir: Path) -> tuple[Path, Path]:
        """First forward FASTQ and its mate, according to the configured file pattern."""
        for item in sorted(sample_dir.iterdir()):
            match = re.fullmatch(self.pattern, item.name)
            if match and match[1] in {"1", "forward"}:
                mate = {"1": "2", "forward": "reverse"}[match[1]]
                partner = item.name[: match.start(1)] + mate + item.name[match.end(1) :]
                return item, sample_dir / partner
        raise FileNotFoundError(
            f"no forward FASTQ matching {self.pattern!r} in {sample_dir}"
        )

    def reference_of(self, case: Case) -> str:
        return "".join(self.load_sampleinfo(case)[0]["reference"].split())

    def mlst_reference(self, organism: str) -> str | None:
        """Reference folder name microSALT resolves an organism name to."""
        referencer = Referencer(
            log=None,
            folders=self.ms.folders,
            threshold=self.ms.threshold,
            pubmlst=self.ms.pubmlst,
            pasteur=self.ms.pasteur,
            singularity=self.ms.singularity,
            containers=self.ms.containers,
        )
        return referencer.organism2reference(organism)

    def synthetic_sample(
        self, base: Case, project: str, index: int, **overrides
    ) -> tuple[dict, Path]:
        """Single-sample sample info under a synthetic project, pointing at the base case FASTQ."""
        template = self.load_sampleinfo(base)[0]
        sample = dict(template)
        sample_id = f"{project}A{index}"
        sample.update(
            {
                "CG_ID_project": project,
                "CG_ID_sample": sample_id,
                "Customer_ID_project": project,
                "Customer_ID_sample": f"val-{sample_id}",
            }
        )
        sample.update(overrides)
        return sample, base.fastq / template["CG_ID_sample"]
