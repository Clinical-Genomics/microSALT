#!/usr/bin/env python
"""Drive a microSALT release validation from a YAML config.

Phases (run in order, waiting for SLURM between the submitting ones):

    run_validation.py -c validation_config.yaml start submit     # first round of analyses
    run_validation.py -c validation_config.yaml status           # poll until everything is complete
    run_validation.py -c validation_config.yaml start resubmit   # second round (suites that re-run cases)
    run_validation.py -c validation_config.yaml status
    run_validation.py -c validation_config.yaml start finish     # login-node finish steps and checks
    run_validation.py -c validation_config.yaml collect          # completion, typing comparison, results.md

Which suites take part in which phase is defined by the suites enabled in the config. Run it
with the python of the microSALT environment so the microSALT package is importable.
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from core.config import load_validation_config  # noqa: E402
from core.orchestrator import Orchestrator  # noqa: E402
from core.state import State  # noqa: E402
from core.suite import PHASES  # noqa: E402
from suites import REGISTRY  # noqa: E402


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-c", "--config", help="validation YAML config (required except for `suites`)"
    )
    sub = parser.add_subparsers(dest="command", required=True)

    start = sub.add_parser("start", help="start a phase")
    start.add_argument("phase", choices=PHASES)
    start.add_argument(
        "--force",
        action="store_true",
        help="continue even if earlier runs are not complete",
    )

    sub.add_parser("status", help="show completion of submitted runs")
    sub.add_parser("suites", help="list available suites")

    collect = sub.add_parser("collect", help="gather results and write results.md")
    collect.add_argument(
        "--no-copy",
        action="store_true",
        help="do not copy run folders to the validation root",
    )
    collect.add_argument(
        "-o",
        "--output",
        help="results markdown path (default: <validation_root>/results.md)",
    )

    args = parser.parse_args(argv)
    if args.command == "suites":
        for name, suite in sorted(REGISTRY.items()):
            print(
                f"{name:<14} {suite.title}  (requires: {', '.join(suite.requires) or '-'}; module {suite.__module__})"
            )
        return 0
    if not args.config:
        parser.error("-c/--config is required")

    cfg = load_validation_config(args.config)
    state = State(cfg.state_path)
    orchestrator = Orchestrator(cfg, state, REGISTRY)
    print(f"validation root: {cfg.validation_root}; suites: {', '.join(cfg.suites)}")

    if args.command == "start":
        orchestrator.start(args.phase, force=args.force)
    elif args.command == "status":
        orchestrator.status()
    else:
        output = (
            Path(args.output) if args.output else cfg.validation_root / "results.md"
        )
        orchestrator.collect(copy_results=not args.no_copy, output=output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
