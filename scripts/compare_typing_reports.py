#!/usr/bin/env python
"""Compare typing results between two microSALT rendered HTML reports.

Usage:
    python scripts/compare_typing_reports.py report_v1.html report_v2.html

Compares:
  - Per-sample sequence types (ST calls)
  - Per-sample MLST allele assignments
  - Per-sample resistance gene calls
  - QC threshold status

Outputs a human-readable diff to stdout. Exit code 0 when identical, 1 when
differences are found, 2 on parse errors.
"""

import argparse
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

from bs4 import BeautifulSoup, Tag

# ---------------------------------------------------------------------------
# Data containers
# ---------------------------------------------------------------------------


@dataclass
class SampleSummary:
    customer_id: str
    cg_id: str
    organism: str
    sequence_type: str
    threshold: str


@dataclass
class MLSTRow:
    loci: str
    allele: str
    identity: str
    span: str


@dataclass
class ResistanceRow:
    gene: str
    group: str
    reference: str
    identity: str
    span: str


@dataclass
class SampleDetail:
    cg_id: str
    mlst: list[MLSTRow] = field(default_factory=list)
    resistances: list[ResistanceRow] = field(default_factory=list)


@dataclass
class Report:
    path: Path
    project_id: str
    report_version: str
    summaries: dict[str, SampleSummary] = field(default_factory=dict)  # keyed by CG Prov ID
    details: dict[str, SampleDetail] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------


def _clean(text: str) -> str:
    """Strip whitespace and normalise unicode spaces."""
    return text.replace("\xa0", " ").strip()


def _find_table_by_header(soup: BeautifulSoup, *header_snippets: str) -> Tag | None:
    """Return the first <table> whose <th> text contains all given snippets."""
    for table in soup.find_all("table"):
        headers = [_clean(th.get_text()) for th in table.find_all("th")]
        if all(any(s in h for h in headers) for s in header_snippets):
            return table
    return None


def _rows(table: Tag) -> list[list[str]]:
    """Return non-header rows as lists of cell text."""
    result = []
    for tr in table.find_all("tr"):
        cells = tr.find_all("td")
        if cells:
            result.append([_clean(c.get_text()) for c in cells])
    return result


def parse_report(path: Path) -> Report:
    html = path.read_text(encoding="utf-8")
    soup = BeautifulSoup(html, "html.parser")

    # --- project summary metadata ----------------------------------------
    project_id = ""
    report_version = ""

    # Look for the "Projektsammanställning" table (key/value pairs via <td>)
    for table in soup.find_all("table"):
        td_texts = [_clean(td.get_text()) for td in table.find_all("td")]
        for i, txt in enumerate(td_texts):
            if "CG Projekt ID" in txt and i + 1 < len(td_texts):
                # e.g. "ACC1234 (CG1234)" — extract the CG id inside parens
                val = td_texts[i + 1]
                m = re.search(r"\(([^)]+)\)", val)
                project_id = m.group(1) if m else val
            if "Rapport version" in txt and i + 1 < len(td_texts):
                report_version = td_texts[i + 1]

    # --- summary table ---------------------------------------------------
    summaries: dict[str, SampleSummary] = {}
    summary_table = _find_table_by_header(soup, "Sekvenstyp", "Tröskelv")
    if summary_table:
        for row in _rows(summary_table):
            if len(row) >= 5:
                s = SampleSummary(
                    customer_id=row[0],
                    cg_id=row[1],
                    organism=row[2],
                    sequence_type=row[3],
                    threshold=row[4],
                )
                summaries[s.cg_id] = s

    # --- per-sample detail tables ----------------------------------------
    details: dict[str, SampleDetail] = {}

    # Each sample detail section has an overview table containing "CG Prov ID".
    # We walk all tables and use context to pair MLST/resistance tables with
    # the sample they belong to.
    current_cg_id: str | None = None
    current_detail: SampleDetail | None = None

    for table in soup.find_all("table"):
        td_texts = [_clean(td.get_text()) for td in table.find_all("td")]
        th_texts = [_clean(th.get_text()) for th in table.find_all("th")]

        # Detect sample overview table (contains "Prov ID (CG Prov ID)" header row)
        if any("CG Prov ID" in t for t in td_texts):
            # Reset context before extraction so a failed parse doesn't leak
            # the previous sample's detail into the tables that follow.
            current_detail = None
            for i, t in enumerate(td_texts):
                if "CG Prov ID" in t and i + 1 < len(td_texts):
                    # value cell may look like "ACC5551  (CG5551)"
                    m = re.search(r"\(([^)]+)\)", td_texts[i + 1])
                    cg_id = m.group(1) if m else td_texts[i + 1].split()[-1]
                    current_detail = SampleDetail(cg_id=cg_id)
                    details[cg_id] = current_detail
                    break
            continue

        # MLST table
        if current_detail and "Loci" in th_texts and "Allel" in th_texts:
            for row in _rows(table):
                if len(row) >= 4:
                    current_detail.mlst.append(
                        MLSTRow(
                            loci=row[1],
                            allele=row[2],
                            identity=row[3],
                            span=row[4] if len(row) > 4 else "",
                        )
                    )
            continue

        # Resistance table
        if current_detail and "Gen" in th_texts and "Grupp" in th_texts:
            for row in _rows(table):
                if len(row) >= 5:
                    current_detail.resistances.append(
                        ResistanceRow(
                            gene=row[1],
                            group=row[2],
                            reference=row[3],
                            identity=row[4],
                            span=row[5] if len(row) > 5 else "",
                        )
                    )

    return Report(
        path=path,
        project_id=project_id,
        report_version=report_version,
        summaries=summaries,
        details=details,
    )


# ---------------------------------------------------------------------------
# Comparison logic
# ---------------------------------------------------------------------------

_RESET = "\033[0m"
_RED = "\033[31m"
_GREEN = "\033[32m"
_YELLOW = "\033[33m"
_BOLD = "\033[1m"


def _diff_label(a: str, b: str) -> str:
    if a == b:
        return a
    return f"{_RED}{a}{_RESET} → {_GREEN}{b}{_RESET}"


def compare_reports(r1: Report, r2: Report) -> int:
    """Print differences between two reports.  Returns number of differences."""
    diffs = 0

    print(
        f"\n{_BOLD}Report A:{_RESET} {r1.path}  (version {r1.report_version}, project {r1.project_id})"
    )
    print(
        f"{_BOLD}Report B:{_RESET} {r2.path}  (version {r2.report_version}, project {r2.project_id})\n"
    )

    all_cg_ids = sorted(set(r1.summaries) | set(r2.summaries))
    if not all_cg_ids:
        print(
            f"{_YELLOW}WARNING: No samples found — check that the HTML files are valid microSALT typing reports.{_RESET}"
        )
        return 0

    for cg_id in all_cg_ids:
        s1 = r1.summaries.get(cg_id)
        s2 = r2.summaries.get(cg_id)

        if s1 is None:
            print(f"{_GREEN}+ {cg_id}: present only in report B{_RESET}")
            diffs += 1
            continue
        if s2 is None:
            print(f"{_RED}- {cg_id}: present only in report A{_RESET}")
            diffs += 1
            continue

        sample_diffs = []

        if s1.sequence_type != s2.sequence_type:
            sample_diffs.append(f"  ST:        {_diff_label(s1.sequence_type, s2.sequence_type)}")
        if s1.organism != s2.organism:
            sample_diffs.append(f"  Organism:  {_diff_label(s1.organism, s2.organism)}")
        if s1.threshold != s2.threshold:
            sample_diffs.append(f"  Threshold: {_diff_label(s1.threshold, s2.threshold)}")

        # MLST allele comparison
        d1 = r1.details.get(cg_id)
        d2 = r2.details.get(cg_id)
        if d1 and d2:
            loci1 = {m.loci: m.allele for m in d1.mlst}
            loci2 = {m.loci: m.allele for m in d2.mlst}
            all_loci = sorted(set(loci1) | set(loci2))
            for locus in all_loci:
                a1 = loci1.get(locus, "—")
                a2 = loci2.get(locus, "—")
                if a1 != a2:
                    sample_diffs.append(f"  MLST [{locus}]: {_diff_label(a1, a2)}")

            # Resistance comparison (gene-level set diff)
            genes1 = {r.gene for r in d1.resistances}
            genes2 = {r.gene for r in d2.resistances}
            for g in sorted(genes1 - genes2):
                sample_diffs.append(f"  {_RED}- resistance: {g}{_RESET}")
            for g in sorted(genes2 - genes1):
                sample_diffs.append(f"  {_GREEN}+ resistance: {g}{_RESET}")

        if sample_diffs:
            diffs += len(sample_diffs)
            print(f"{_BOLD}{cg_id}{_RESET} ({s1.customer_id}):")
            for line in sample_diffs:
                print(line)
            print()

    if diffs == 0:
        print(f"{_GREEN}No differences found.{_RESET}")
    else:
        print(f"{_YELLOW}{diffs} difference(s) found.{_RESET}")

    return diffs


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Compare typing results between two microSALT rendered HTML reports."
    )
    parser.add_argument("report_a", type=Path, help="Path to first (older) HTML report")
    parser.add_argument("report_b", type=Path, help="Path to second (newer) HTML report")
    parser.add_argument("--no-color", action="store_true", help="Disable ANSI color output")
    args = parser.parse_args()

    if args.no_color:
        global _RESET, _RED, _GREEN, _YELLOW, _BOLD
        _RESET = _RED = _GREEN = _YELLOW = _BOLD = ""

    for p in (args.report_a, args.report_b):
        if not p.exists():
            print(f"ERROR: File not found: {p}", file=sys.stderr)
            return 2

    try:
        r1 = parse_report(args.report_a)
        r2 = parse_report(args.report_b)
    except Exception as exc:
        print(f"ERROR: Failed to parse report — {exc}", file=sys.stderr)
        return 2

    diffs = compare_reports(r1, r2)
    return 0 if diffs == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
