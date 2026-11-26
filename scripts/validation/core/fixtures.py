"""Helpers for building modified inputs without touching real FASTQ, reference or database files.

Everything is created under <validation_root>/fixtures/<label>/ as copies or symlink farms.
"""

import json
import os
import re
import shutil
import sqlite3
import tempfile
from contextlib import contextmanager
from pathlib import Path

from microSALT.utils.referencer import Referencer

from .config import Case, ValidationConfig


class FixtureHelper:
    def __init__(self, cfg: ValidationConfig):
        self.cfg = cfg
        self.ms = cfg.microsalt
        self.pattern = self.ms.regex.file_pattern
        self._retry_dir: Path | None = None

    @contextmanager
    def retry_fixtures(self):
        """Put retry fixtures in a new directory without touching earlier attempts."""
        parent = self.cfg.fixture_dir / "retries"
        parent.mkdir(parents=True, exist_ok=True)
        self._retry_dir = Path(tempfile.mkdtemp(prefix="attempt-", dir=parent))
        try:
            yield
        finally:
            self._retry_dir = None

    def fresh(self, label: str) -> Path:
        path = (self._retry_dir or self.cfg.fixture_dir) / label
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
        path = self.write_json(
            (self._retry_dir or self.cfg.fixture_dir) / f"{case.name}_{label}.json",
            data,
        )
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

    def empty_profile_sqlite(self, root: Path, scheme: str) -> tuple[Path, str] | None:
        """Create an isolated SQLite DB with the selected MLST profile table empty.

        This works even when stage uses MySQL: R7 tests the profile import guard, not the
        production DB engine, so a scratch database avoids modifying shared reference data.
        """
        scratch = root / "microsalt_scratch.db"
        profile_path = Path(self.ms.folders.profiles) / scheme
        if not profile_path.is_file():
            return None
        with profile_path.open() as stream:
            columns = stream.readline().rstrip().split("\t")[:8]
        if not columns or "ST" not in columns:
            return None

        def quote(identifier: str) -> str:
            return '"' + identifier.replace('"', '""') + '"'

        definitions = ", ".join(
            f"{quote(column)} INTEGER" + (" PRIMARY KEY" if column == "ST" else "")
            for column in columns
        )
        with sqlite3.connect(scratch) as conn:
            conn.execute(f'CREATE TABLE "profile_{scheme}" ({definitions})')
        return scratch, f"sqlite:///{scratch}"

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
