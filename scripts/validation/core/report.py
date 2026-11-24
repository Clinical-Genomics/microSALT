"""Render the Results section of the validation report from the recorded state."""

import getpass
from datetime import date

from .config import ValidationConfig
from .state import State
from .suite import Suite


def yn(check: dict | None) -> str:
    if check is None:
        return "-"
    return {True: "Y", False: "N", None: "N/A"}[check["passed"]]


def detail(check: dict | None) -> str:
    return check["detail"].replace("|", "/") if check else ""


def table(header: list[str], rows: list[list[str]]) -> list[str]:
    lines = [
        "| " + " | ".join(header) + " |",
        "| " + " | ".join("-" * len(h) for h in header) + " |",
    ]
    lines.extend("| " + " | ".join(str(cell) for cell in row) + " |" for row in rows)
    return lines


def render_report(cfg: ValidationConfig, state: State, suites: list[Suite]) -> str:
    today = date.today().isoformat()
    user = getpass.getuser()
    lines = ["## Results", ""]
    lines += table(
        [
            "Step of validation",
            "Date of performance",
            "Performed by",
            "Result vs criteria",
        ],
        [[suite.title, today, user, suite.result()] for suite in suites],
    )
    lines.append("")

    if cfg.manual_steps:
        lines += ["### Manual steps", ""]
        lines += table(
            ["Step", "Instruction", "Result"],
            [
                [step, instruction, "<!-- fill in -->"]
                for step, instruction in cfg.manual_steps.items()
            ],
        )
        lines.append("")

    for suite in suites:
        lines += [f"### {suite.title} - Results", ""]
        lines += suite.render()
        lines.append("")

    lines += ["### **Raw data**", ""]
    lines.append(
        f"All run folders are copied to `{cfg.validation_root / 'results'}`; command logs are in "
        f"`{cfg.log_dir}` and the full machine-readable record is `{cfg.state_path}`."
    )
    lines.append("")
    lines += table(
        ["Run", "Finish folder", "SLURM job ids"],
        [
            [record.id, record.finishdir, ", ".join(record.slurm_ids) or "-"]
            for record in state.runs()
            if record.phase != "finish" and record.finishdir
        ],
    )
    lines.append("")
    return "\n".join(lines) + "\n"
