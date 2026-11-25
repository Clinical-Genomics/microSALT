"""Read-only snapshots of the microSALT database used by the checks."""

from microSALT.store.database import get_session, initialize_database
from microSALT.store.orm_models import (
    Projects,
    Reports,
    Resistances,
    Samples,
    Seq_types,
)

QC_METRICS = (
    "total_reads",
    "raw_reads",
    "insert_size",
    "duplication_rate",
    "mapped_rate",
    "coverage_10x",
    "coverage_30x",
    "coverage_50x",
    "coverage_100x",
    "average_coverage",
    "reference_length",
)
METADATA = (
    "organism",
    "reference_genome",
    "priority",
    "application_tag",
    "method_libprep",
    "method_sequencing",
)
TYPING = ("ST", "pubmlst_ST", "genome_length", "n50", "contigs")


class Database:
    def __init__(self, uri: str):
        initialize_database(uri)
        self.session = get_session()

    def snapshot_project(self, project: str) -> dict:
        """Per-sample typing, metadata and QC values plus typing row counts."""
        self.session.expire_all()
        samples = {}
        rows = (
            self.session.query(Samples).filter(Samples.CG_ID_project == project).all()
        )
        for sample in rows:
            sid = sample.CG_ID_sample
            samples[sid] = {
                "typing": {k: getattr(sample, k) for k in TYPING},
                "metadata": {k: getattr(sample, k) for k in METADATA},
                "qc": {k: getattr(sample, k) for k in QC_METRICS},
                "seq_types": sorted(
                    (r.loci, r.allele, r.contig_name)
                    for r in self.session.query(Seq_types).filter(
                        Seq_types.CG_ID_sample == sid
                    )
                ),
                "resistances": sorted(
                    (r.gene, r.instance, r.contig_name)
                    for r in self.session.query(Resistances).filter(
                        Resistances.CG_ID_sample == sid
                    )
                ),
            }
        return {"samples": samples, "reports": self.report_versions(project)}

    def row_counts(self) -> dict[str, int]:
        return {
            "projects": self.session.query(Projects).count(),
            "samples": self.session.query(Samples).count(),
            "reports": self.session.query(Reports).count(),
        }

    def report_versions(self, project: str) -> list[dict]:
        self.session.expire_all()
        reports = (
            self.session.query(Reports)
            .filter(Reports.CG_ID_project == project)
            .order_by(Reports.version)
            .all()
        )
        return [
            {
                "version": r.version,
                "fingerprint": r.steps_aggregate,
                "date": str(r.date),
            }
            for r in reports
        ]

    def sample(self, sample_id: str) -> Samples | None:
        self.session.expire_all()
        return (
            self.session.query(Samples)
            .filter(Samples.CG_ID_sample == sample_id)
            .first()
        )

    def project_sample_ids(self, project: str) -> list[str]:
        self.session.expire_all()
        rows = self.session.query(Samples.CG_ID_sample).filter(
            Samples.CG_ID_project == project
        )
        return sorted(row[0] for row in rows)

    def report_sample_ids(
        self, project: str, sample_ids: list[str] | None
    ) -> list[str]:
        from microSALT.server.views import gen_reportdata

        data = gen_reportdata(project, sample_ids=sample_ids)
        return sorted(s.CG_ID_sample for s in data["samples"])


def diff_typing(before: dict, after: dict) -> list[str]:
    """Differences in ST, alleles and resistances between two project snapshots."""
    problems = []
    for sid, old in before["samples"].items():
        new = after["samples"].get(sid)
        if new is None:
            problems.append(f"{sid}: sample missing after QC-only run")
            continue
        for key in ("ST", "pubmlst_ST"):
            if old["typing"][key] != new["typing"][key]:
                problems.append(
                    f"{sid}: {key} {old['typing'][key]} -> {new['typing'][key]}"
                )
        for table in ("seq_types", "resistances"):
            if [tuple(r) for r in old[table]] != [tuple(r) for r in new[table]]:
                problems.append(
                    f"{sid}: {table} rows {len(old[table])} -> {len(new[table])} or content changed"
                )
    return problems


def diff_qc(before: dict, after: dict, rel_tol: float = 0.01) -> list[str]:
    """QC metric differences beyond a relative tolerance."""
    problems = []
    for sid, old in before["samples"].items():
        new = after["samples"].get(sid)
        if new is None:
            continue
        for key in QC_METRICS:
            a, b = old["qc"][key], new["qc"][key]
            if a is None or b is None:
                if a != b:
                    problems.append(f"{sid}: {key} {a} -> {b}")
                continue
            if a == b:
                continue
            if abs(a - b) > rel_tol * max(abs(a), abs(b), 1):
                problems.append(f"{sid}: {key} {a} -> {b}")
    return problems
