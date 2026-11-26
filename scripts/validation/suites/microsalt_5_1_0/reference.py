"""Reference handling: whitespace normalisation, NCBI assembly accessions, corrupt references,
organisms without an MLST scheme and incomplete MLST profile tables."""

from dataclasses import dataclass
from pathlib import Path

from core.report import detail, table, yn
from core.runner import setup_errors
from core.slurm import run_complete, run_status
from core.suite import RetryAction, Suite, SuiteOptions

from .qc_only import runfiles


@dataclass
class Fixture:
    id: str
    description: str
    expect: str  # "error", "submit" or "skip"
    sampleinfo: Path
    fastq: Path
    job_name: str
    config_path: Path | None = None
    skip_update: bool = False
    note: str = ""


class ReferenceOptions(SuiteOptions):
    base_case: str
    project: str = "VAL5100"
    assembly_accession: str
    bogus_accession: str = "GCF_000000000.0"
    organism_without_mlst: str


class ReferenceSuite(Suite):
    name = "reference"
    title = "Reference handling"
    Options = ReferenceOptions
    options: ReferenceOptions

    def _fixtures(self, selected: set[str] | None = None) -> dict[str, Fixture]:
        fx = self.fixtures
        opt = self.options
        base = self.cfg.case(opt.base_case)
        reference = fx.reference_of(base)
        out = {}

        def add(
            test_id,
            index,
            description,
            expect,
            config_path=None,
            skip_update=False,
            **overrides,
        ):
            if selected is not None and test_id not in selected:
                return None
            root = fx.fresh(test_id)
            sample, fastq = fx.synthetic_sample(base, opt.project, index, **overrides)
            out[test_id] = Fixture(
                test_id,
                description,
                expect,
                fx.write_json(root / "sampleinfo.json", sample),
                fastq,
                sample["CG_ID_sample"],
                config_path=config_path,
                skip_update=skip_update,
            )
            return root

        add(
            "R1",
            1,
            f"Reference '  {reference}  ' with surrounding whitespace",
            "submit",
            reference=f"  {reference}  ",
        )
        add(
            "R2",
            2,
            f"NCBI assembly accession {opt.assembly_accession} as reference",
            "submit",
            reference=opt.assembly_accession,
        )
        add(
            "R3",
            3,
            f"Non-existing assembly accession {opt.bogus_accession}",
            "error",
            reference=opt.bogus_accession,
        )

        if selected is None or "R4" in selected:
            root = fx.fresh("R4")
            genomes = fx.genomes_farm(root)
            corrupt = genomes / f"{opt.assembly_accession}.fasta"
            corrupt.unlink(missing_ok=True)
            corrupt.write_text(">truncated\nACGT\n>empty\n")
            sample, fastq = fx.synthetic_sample(
                base, opt.project, 4, reference=opt.assembly_accession
            )
            out["R4"] = Fixture(
                "R4",
                f"Corrupt cached download of {opt.assembly_accession} (empty record)",
                "error",
                fx.write_json(root / "sampleinfo.json", sample),
                fastq,
                sample["CG_ID_sample"],
                config_path=fx.alt_config(root, genomes=genomes),
            )

        add(
            "R5",
            5,
            f"Organism without MLST scheme: {opt.organism_without_mlst}",
            "submit",
            organism=opt.organism_without_mlst,
        )
        add(
            "R6a",
            6,
            "Full analysis with the real organism (gives a sequence type)",
            "submit",
        )
        add(
            "R6b",
            6,
            f"Re-run same sample as {opt.organism_without_mlst} (stale ST must clear)",
            "submit",
            organism=opt.organism_without_mlst,
        )

        if selected is None or "R7" in selected:
            root = fx.fresh("R7")
            sample, fastq = fx.synthetic_sample(base, opt.project, 7)
            sampleinfo = fx.write_json(root / "sampleinfo.json", sample)
            scheme = fx.mlst_reference(fx.load_sampleinfo(base)[0]["organism"])
            scratch = fx.empty_profile_sqlite(root, scheme) if scheme else None
            if scratch is None or scheme is None:
                out["R7"] = Fixture(
                    "R7",
                    "Empty MLST profile table",
                    "skip",
                    sampleinfo,
                    fastq,
                    sample["CG_ID_sample"],
                    note=(
                        "requires a known MLST scheme and a readable profile file with an ST column; "
                        f"scheme={scheme}, profiles folder={self.cfg.microsalt.folders.profiles}"
                    ),
                )
            else:
                _, uri = scratch
                out["R7"] = Fixture(
                    "R7",
                    f"Empty profile_{scheme} in an isolated scratch SQLite database, run with --skip_update",
                    "error",
                    sampleinfo,
                    fastq,
                    sample["CG_ID_sample"],
                    config_path=fx.alt_config(root, db_uri=uri),
                    skip_update=True,
                )
        return out

    def _run(self, phase: str, fx: Fixture):
        base = self.options.base_case
        return self.analyse(
            phase,
            fx.id,
            base,
            fx.sampleinfo,
            fx.fastq,
            fx.job_name,
            skip_update=fx.skip_update,
            config_path=fx.config_path,
            note=fx.description,
        )

    def submit(self) -> None:
        print("== Reference handling (R tests) ==")
        fixtures = self._fixtures()
        descriptions = {k: v.description for k, v in fixtures.items()}
        self.snapshot("descriptions", descriptions)
        for test_id in ("R1", "R2", "R3", "R4", "R5", "R6a", "R7"):
            self._submit_fixture(fixtures[test_id])

    def _submit_fixture(self, fx: Fixture) -> None:
        test_id = fx.id
        if fx.expect == "skip":
            self.check(f"{test_id}/result", None, f"skipped: {fx.note}")
            return
        record, result = self._run("submit", fx)
        errors = setup_errors(result.output)
        if fx.expect == "error":
            self.check(
                f"{test_id}/result",
                not record.ok and not record.slurm_ids,
                f"exit {record.exit_code}; " + " | ".join(errors[:3]),
            )
        else:
            self.check(
                f"{test_id}/result",
                record.ok and bool(record.slurm_ids),
                f"exit {record.exit_code}; slurm ids {record.slurm_ids}; log: {record.log}",
            )
        if test_id == "R5" and record.ok:
            content = "".join(p.read_text() for p in runfiles(record.finishdir))
            self.check(
                "R5/mlst_skipped",
                "skipping MLST" in result.output
                and "blast_search/mlst" not in content
                and "blast_search/resistance" in content,
                "MLST warning logged; runfile has resistance but no MLST blast section",
            )

    def resubmit(self) -> None:
        if self.run("submit", "R6a") is None:
            return
        fx = self._fixtures({"R6b"})["R6b"]
        sample = self.ctx.db.sample(fx.job_name)
        self.snapshot("R6/st_before", {"ST": sample.ST if sample else None})
        self._run("resubmit", fx)

    def retry_actions(self) -> list[RetryAction]:
        actions = [
            RetryAction(
                f"submit/reference/{test_id}",
                f"Reference handling {test_id}",
                "submit",
                self.run_id("submit", test_id),
            )
            for test_id in ("R1", "R2", "R3", "R4", "R5", "R6a", "R7")
        ]
        actions.append(
            RetryAction(
                "resubmit/reference/R6b",
                "R6b: rerun without MLST",
                "resubmit",
                self.run_id("resubmit", "R6b"),
            )
        )
        actions.append(
            RetryAction(
                "finish/reference",
                "Recheck reference completion and R6 results",
                "finish",
            )
        )
        return actions

    def retry(self, action_id: str) -> None:
        if action_id == "finish/reference":
            self.finish()
            return
        phase, _, test_id = action_id.split("/")
        if test_id == "R6b":
            first = self.run("submit", "R6a")
            if not first or not run_complete(first):
                raise ValueError("R6a must complete before R6b can be retried")
            if (
                self.run("resubmit", "R6b")
                and self.get_snapshot("R6/st_before") is None
            ):
                raise ValueError(
                    "The original R6 sequence type snapshot is missing; start a new validation"
                )
        if test_id == "R6a" and self.run("resubmit", "R6b"):
            raise ValueError(
                "R6b has already run; start a new validation to repeat the complete R6 sequence"
            )
        with self.fixtures.retry_fixtures():
            fx = self._fixtures({test_id})[test_id]
            descriptions = self.get_snapshot("descriptions") or {}
            descriptions[test_id] = fx.description
            self.snapshot("descriptions", descriptions)
            self.state.invalidate_checks(
                f"reference/{'R6' if test_id == 'R6b' else test_id}/",
                "Retry in progress; result needs rechecking",
            )
            if test_id == "R6b":
                if self.get_snapshot("R6/st_before") is None:
                    sample = self.ctx.db.sample(fx.job_name)
                    self.snapshot("R6/st_before", {"ST": sample.ST if sample else None})
                self._run(phase, fx)
            else:
                self._submit_fixture(fx)

    def finish(self) -> None:
        print("== Reference handling: completion ==")
        opt = self.options
        for test_id in ("R1", "R2", "R5"):
            record = self.run("submit", test_id)
            if record:
                self.check(
                    f"{test_id}/completed",
                    run_complete(record),
                    f"{run_status(record)} ({record.finishdir})",
                )
        r2 = self.run("submit", "R2")
        if r2 and r2.ok:
            fasta = (
                Path(self.cfg.microsalt.folders.genomes)
                / f"{opt.assembly_accession}.fasta"
            )
            sample = self.ctx.db.sample(f"{opt.project}A2")
            reference_genome = (
                str(sample.reference_genome)
                if sample is not None and sample.reference_genome is not None
                else None
            )
            self.check(
                "R2/downloaded",
                fasta.exists() and reference_genome == opt.assembly_accession,
                f"{fasta} exists: {fasta.exists()}; reference_genome in db: {reference_genome}",
            )
        if self.run("resubmit", "R6b"):
            before = (self.get_snapshot("R6/st_before") or {}).get("ST")
            sample = self.ctx.db.sample(f"{opt.project}A6")
            after_value = sample.ST if sample is not None else None
            after = after_value if isinstance(after_value, int) else None
            self.check(
                "R6/result",
                before not in (None, -1) and after == -1,
                f"ST before {before}, after re-run without MLST scheme {after}",
            )

    def result(self) -> str:
        required = (
            "R1/result",
            "R1/completed",
            "R2/result",
            "R2/completed",
            "R2/downloaded",
            "R3/result",
            "R4/result",
            "R5/result",
            "R5/mlst_skipped",
            "R5/completed",
            "R6a/result",
            "R6/result",
            "R7/result",
        )
        checks = [self.get_check(key) for key in required]
        if any(check and check["passed"] is False for check in checks):
            return "FAIL"
        if any(not check or check["passed"] is None for check in checks):
            return "INCOMPLETE"
        return "PASS"

    def render(self) -> list[str]:
        descriptions = self.get_snapshot("descriptions") or {}
        rows = []
        for test_id, key, extra in (
            ("R1", "result", "completed"),
            ("R2", "result", "downloaded"),
            ("R2", "completed", None),
            ("R3", "result", None),
            ("R4", "result", None),
            ("R5", "result", "mlst_skipped"),
            ("R5", "completed", None),
            ("R6", "result", None),
            ("R7", "result", None),
        ):
            check = self.get_check(f"{test_id}/{key}")
            comment = detail(check)
            if extra and (extra_check := self.get_check(f"{test_id}/{extra}")):
                comment += f"; {extra}: {yn(extra_check)} ({detail(extra_check)})"
            label = test_id if key == "result" else f"{test_id} ({key})"
            rows.append(
                [
                    label,
                    descriptions.get(test_id) or descriptions.get(f"{test_id}b", ""),
                    yn(check),
                    comment,
                ]
            )
        return table(
            ["Test", "Description", "Result as expected (Y/N)", "Comment"], rows
        )
