"""Validation before submission: invalid sample info, FASTQ or reference inputs are rejected with
all errors listed, no SLURM jobs submitted and no database rows written."""

from dataclasses import dataclass
from pathlib import Path

from core.report import detail, table, yn
from core.runner import setup_errors
from core.suite import Suite, SuiteOptions


@dataclass
class Fixture:
    id: str
    description: str
    sampleinfo: Path
    fastq: Path
    config_path: Path | None = None
    skip_update: bool = False
    min_errors: int = 1


class ValidationOptions(SuiteOptions):
    base_case: str
    max_validation_seconds: float = 60


class ValidationSuite(Suite):
    name = "validation"
    title = "Validation before submission"
    Options = ValidationOptions

    def _fixtures(self) -> list[Fixture]:
        fx = self.fixtures
        base = self.cfg.case(self.options.base_case)
        samples = fx.load_sampleinfo(base)
        first = samples[0]["CG_ID_sample"]
        reference = fx.reference_of(base)
        out = []

        root = fx.fresh("V1")
        path, _ = fx.modified_sampleinfo(base, "V1", organism=None)
        out.append(
            Fixture("V1", f"Remove 'organism' from sample {first}", path, base.fastq)
        )

        root = fx.fresh("V2")
        farm = fx.fastq_farm(root, base)
        fwd, _ = fx.mates(farm / first)
        fwd.rename(fwd.parent / "badname_R1.fastq.gz")
        out.append(
            Fixture(
                "V2",
                f"Rename {fwd.name} of {first} to badname_R1.fastq.gz",
                fx.write_json(root / "sampleinfo.json", samples),
                farm,
            )
        )

        root = fx.fresh("V3")
        farm = fx.fastq_farm(root, base)
        _, rev = fx.mates(farm / first)
        rev.unlink()
        out.append(
            Fixture(
                "V3",
                f"Remove mate {rev.name} of {first}",
                fx.write_json(root / "sampleinfo.json", samples),
                farm,
            )
        )

        root = fx.fresh("V4")
        farm = fx.fastq_farm(root, base)
        fwd, _ = fx.mates(farm / first)
        fwd.unlink()
        fwd.touch()
        fwd.chmod(0)
        out.append(
            Fixture(
                "V4",
                f"Make {fwd.name} of {first} unreadable (chmod 000)",
                fx.write_json(root / "sampleinfo.json", samples),
                farm,
            )
        )

        root = fx.fresh("V5")
        genomes = fx.genomes_farm(root)
        for item in genomes.glob(f"{reference}.fasta*"):
            item.unlink()
        out.append(
            Fixture(
                "V5",
                f"Remove {reference}.fasta and its indexes, run with --skip_update",
                fx.write_json(root / "sampleinfo.json", samples),
                base.fastq,
                config_path=fx.alt_config(root, genomes=genomes),
                skip_update=True,
            )
        )

        root = fx.fresh("V6")
        genomes = fx.genomes_farm(root)
        fasta = genomes / f"{reference}.fasta"
        fasta.unlink(missing_ok=True)
        fasta.write_text("<html><body>Error 500</body></html>\n")
        out.append(
            Fixture(
                "V6",
                f"Replace {reference}.fasta with an HTML error page",
                fx.write_json(root / "sampleinfo.json", samples),
                base.fastq,
                config_path=fx.alt_config(root, genomes=genomes),
                skip_update=True,
            )
        )

        root = fx.fresh("V7")
        genomes = fx.genomes_farm(root)
        (genomes / f"{reference}.fasta.bwt").unlink(missing_ok=True)
        out.append(
            Fixture(
                "V7",
                f"Remove {reference}.fasta.bwt, run with --skip_update",
                fx.write_json(root / "sampleinfo.json", samples),
                base.fastq,
                config_path=fx.alt_config(root, genomes=genomes),
                skip_update=True,
            )
        )

        root = fx.fresh("V8")
        data = fx.load_sampleinfo(base)
        del data[0]["organism"]
        farm = fx.fastq_farm(root, base)
        fwd, _ = fx.mates(farm / first)
        fwd.rename(fwd.parent / "badname_R1.fastq.gz")
        genomes = fx.genomes_farm(root)
        for item in genomes.glob(f"{reference}.fasta*"):
            item.unlink()
        out.append(
            Fixture(
                "V8",
                "Combine V1, V2 and V5",
                fx.write_json(root / "sampleinfo.json", data),
                farm,
                config_path=fx.alt_config(root, genomes=genomes),
                skip_update=True,
                min_errors=3,
            )
        )
        return out

    def submit(self) -> None:
        print("== Validation before submission (V tests) ==")
        base = self.cfg.case(self.options.base_case)
        descriptions = {}
        for fx in self._fixtures():
            descriptions[fx.id] = fx.description
            before = self.ctx.db.row_counts()
            record, result = self.analyse(
                "submit",
                fx.id,
                base.name,
                fx.sampleinfo,
                fx.fastq,
                base.project,
                skip_update=fx.skip_update,
                config_path=fx.config_path,
                note=fx.description,
            )
            errors = setup_errors(result.output)
            self.check(
                f"{fx.id}/error_raised",
                not record.ok and bool(errors),
                f"exit {record.exit_code}; {len(errors)} error(s): "
                + " | ".join(errors[:4]),
            )
            self.check(
                f"{fx.id}/no_jobs",
                not record.slurm_ids,
                f"slurm ids: {record.slurm_ids or 'none'}; finishdir: {record.finishdir or 'none'}",
            )
            after = self.ctx.db.row_counts()
            self.check(
                f"{fx.id}/no_db_rows",
                before == after,
                f"before {before}, after {after}",
            )
            if fx.min_errors > 1:
                self.check(
                    f"{fx.id}/all_errors_listed",
                    len(errors) >= fx.min_errors,
                    f"expected >= {fx.min_errors} errors in one message, got {len(errors)}",
                )

        descriptions["V9"] = (
            "Unmodified sample info, --dry (validation only, no submission)"
        )
        record, _ = self.analyse(
            "submit",
            "V9",
            base.name,
            base.sampleinfo,
            base.fastq,
            base.project,
            dry=True,
            note=descriptions["V9"],
        )
        self.check(
            "V9/accepted",
            record.ok,
            f"exit {record.exit_code} after {record.duration_s}s",
        )
        self.check(
            "V9/validation_seconds",
            record.duration_s < self.options.max_validation_seconds,
            f"{record.duration_s}s for validation + job script creation (no FASTQ decompression)",
        )
        self.snapshot("descriptions", descriptions)

    def render(self) -> list[str]:
        descriptions = self.get_snapshot("descriptions") or {}
        rows = []
        for test_id in ("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8"):
            together = self.get_check(f"{test_id}/all_errors_listed")
            rows.append(
                [
                    test_id,
                    descriptions.get(test_id, ""),
                    yn(self.get_check(f"{test_id}/error_raised")),
                    yn(together) if together else "N/A",
                    yn(self.get_check(f"{test_id}/no_jobs")),
                    yn(self.get_check(f"{test_id}/no_db_rows")),
                ]
            )
        rows.append(
            [
                "V9",
                descriptions.get("V9", ""),
                f"accepted: {yn(self.get_check('V9/accepted'))}",
                "N/A",
                "N/A",
                "N/A",
            ]
        )
        lines = table(
            [
                "Test",
                "Modification",
                "Error raised (Y/N)",
                "All errors reported together (Y/N)",
                "No jobs submitted (Y/N)",
                "No DB rows added (Y/N)",
            ],
            rows,
        )
        lines += [
            "",
            f"Validation duration for the unmodified base case: {detail(self.get_check('V9/validation_seconds'))}",
        ]
        lines += [
            "",
            "Error messages per test are stored in the run logs under `logs/submit__validation__V*.log`.",
        ]
        return lines
