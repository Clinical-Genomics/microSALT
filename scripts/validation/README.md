# Release validation runner

Drives the validation described in the `validations` repo
(`docs/pipelines/microsalt/microSALT<version>_validation_report.md`) and collects the results in
the shape of its Results section.

The runner is split into a release-independent **core** plus a generic **stability** suite (what
every microSALT release validation has consisted of: the standard cases run to completion and the
typing reports are concordant with the previous version), and **release-specific suites** that
check the features of a particular release.

```
run_validation.py                   CLI
validation_config.example.yaml
core/                               framework
  config.py       YAML config model
  state.py        state.json: every command, exit code, log, finishdir, SLURM ids, every check
  runner.py       wraps `microsalt analyse | utils finish`
  slurm.py        run_complete.out / sacct
  db.py           database snapshots and diffs
  fixtures.py     helpers for modified inputs (sample info copies, FASTQ/genome symlink farms, ...)
  suite.py        Suite base class, Context, PHASES
  orchestrator.py drives the suites through the phases
  report.py       renders results.md
suites/
  stability.py    core suite, always on
  microsalt_5_1_0/  qc_only, qc_integrity, validation, reference, versioning
```

## Prerequisites

- The release branch is installed on stage, as `hiseq.clinical`:
    `bash /home/proj/stage/servers/resources/hasta.scilifelab.se/update-microsalt.sh [-d] -e stage -b <branch>`
- Sample info files and FASTQ folders for every case exist on stage.
- Old typing reports to compare against are copied somewhere readable.
- A copy of `validation_config.example.yaml` adjusted to the release.

Run everything with the python of the microSALT environment, since the checks import the
microSALT package to read the database:

```shell
PY=/home/proj/stage/bin/miniconda3/envs/S_microsalt/bin/python
$PY scripts/validation/run_validation.py suites                      # list available suites
$PY scripts/validation/run_validation.py -c validation_config.yaml <command>
```

## Phases

| Command          | Does                                                                                                      | Wait for SLURM after? |
| ---------------- | --------------------------------------------------------------------------------------------------------- | --------------------- |
| `start submit`   | First round of analyses (stability runs every case; other suites add their own submissions/tests).        | yes                   |
| `status`         | Prints completion of every submitted run (`run_complete.out` and `sacct`).                                | –                     |
| `start resubmit` | Second round, for suites that re-run cases after the first round completed (e.g. QC-only).                | yes                   |
| `status`         |                                                                                                           | –                     |
| `start finish`   | Login-node steps: `utils finish` sequences and database checks.                                           | no                    |
| `collect`        | Completion checks, `compare_typing_reports.py` against the old reports, copies run folders, `results.md`. | no                    |

A routine release with only `stability` enabled needs `start submit`, `status`, `collect`; the
other phases are no-ops.

`start resubmit` and `start finish` refuse to run while submitted runs of earlier phases are not
complete; pass `--force` to override.

The latest command and check results are recorded in `<validation_root>/state.json`, with earlier
attempts kept in history; command output goes to separate files under `<validation_root>/logs/`.
Re-running a whole phase still repeats all of its analyses; use `retry` to select one instead.

## Retrying a test

Use the same YAML config and validation root as the original run:

```shell
$PY scripts/validation/run_validation.py -c validation_config.yaml retry
# Or list the enabled actions, then select one by ID:
$PY scripts/validation/run_validation.py -c validation_config.yaml retry --list
$PY scripts/validation/run_validation.py -c validation_config.yaml retry --test submit/reference/R5
```

The menu is built from enabled suites and their configured cases. Select a number or `q` to
cancel, then type `yes` to confirm. `--test` skips the menu selection, not the confirmation.
Previous attempts and checks are retained in `state.json` (`run_history` and `check_history`),
logs use new attempt filenames, and modified fixtures are generated under a new
`fixtures/retries/attempt-*/` directory. The latest result is used for the report. Before retrying,
the runner requires all previously recorded SLURM jobs to have terminal states in `sacct`,
including jobs from partially failed submissions. A missing or unknown state blocks the retry.

The R6 test and the QC-integrity and report-versioning checks are sequences; they cannot always
be repeated safely as individual steps. The menu presents those checks at their safe granularity.
For example, after a QC finish has run, restarting the upstream QC analysis requires a new
validation rather than silently reusing stale downstream evidence. After retrying a reference
analysis, wait for its jobs, then select `finish/reference` to refresh completion checks before
collecting the report. A retry can modify the stage database and consume report versions.

## Writing a suite for a new release

1. Create `suites/microsalt_<version>/` with one module per criterion. Any `Suite` subclass with a
    `name` is registered automatically.
1. Subclass `core.suite.Suite`, set `name`, `title`, optional `requires` (other suite names) and an
    `Options` pydantic model for the suite's YAML options.
1. Implement the phases you need (`submit`, `resubmit`, `finish`, `collect`) and `render()`, which
    returns the markdown lines of the suite's Results subsection. Use `self.analyse(...)`,
    `self.finish_run(...)`, `self.check(...)`, `self.snapshot(...)`, `self.fixtures` (modified
    inputs) and `self.ctx.db` (snapshots).
1. Enable it in the config under `suites:` with its options.

Check ids are `<suite name>/<key>`; run ids are `<phase>/<suite name>/<label>`. Snapshots of
other suites can be read with `self.state.get_snapshot("qc_only/<case>/after_full")` etc.
