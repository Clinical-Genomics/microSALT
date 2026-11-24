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

- The release branch is installed on stage
  (`bash /home/proj/production/servers/resources/hasta.scilifelab.se/update-microsalt.sh -d stage -b <branch>`,
  add `-a` if the release includes an Alembic migration).
- Sample info files and FASTQ folders for every case exist on stage.
- Old typing reports to compare against are copied somewhere readable.
- A copy of `validation_config.example.yaml` adjusted to the release.

Run everything with the python of the microSALT environment, since the checks import the
microSALT package to read the database:

```shell
PY=/home/proj/stage/bin/miniconda3/envs/S_microSALT/bin/python
$PY scripts/validation/run_validation.py suites                      # list available suites
$PY scripts/validation/run_validation.py -c validation_config.yaml <command>
```

## Phases

| Command          | Does                                                                                                  | Wait for SLURM after? |
| ---------------- | ----------------------------------------------------------------------------------------------------- | --------------------- |
| `start submit`   | First round of analyses (stability runs every case; other suites add their own submissions/tests).    | yes                   |
| `status`         | Prints completion of every submitted run (`run_complete.out` and `sacct`).                            | –                     |
| `start resubmit` | Second round, for suites that re-run cases after the first round completed (e.g. QC-only).            | yes                   |
| `status`         |                                                                                                       | –                     |
| `start finish`   | Login-node steps: `utils finish` sequences and database checks.                                       | no                    |
| `collect`        | Completion checks, `compare_typing_reports.py` against the old reports, copies run folders, `results.md`. | no                |

A routine release with only `stability` enabled needs `start submit`, `status`, `collect`; the
other phases are no-ops.

`start resubmit` and `start finish` refuse to run while submitted runs of earlier phases are not
complete; pass `--force` to override.

Every command and check is appended to `<validation_root>/state.json`; command output goes to
`<validation_root>/logs/`. Re-running a phase overwrites the records of that phase.

## Writing a suite for a new release

1. Create `suites/microsalt_<version>/` with one module per criterion. Any `Suite` subclass with a
   `name` is registered automatically.
2. Subclass `core.suite.Suite`, set `name`, `title`, optional `requires` (other suite names) and an
   `Options` pydantic model for the suite's YAML options.
3. Implement the phases you need (`submit`, `resubmit`, `finish`, `collect`) and `render()`, which
   returns the markdown lines of the suite's Results subsection. Use `self.analyse(...)`,
   `self.finish_run(...)`, `self.check(...)`, `self.snapshot(...)`, `self.fixtures` (modified
   inputs) and `self.ctx.db` (snapshots).
4. Enable it in the config under `suites:` with its options.

Check ids are `<suite name>/<key>`; run ids are `<phase>/<suite name>/<label>`. Snapshots of
other suites can be read with `self.state.get_snapshot("qc_only/<case>/after_full")` etc.

## Notes on the 5.1.0 suites

- V5 and V7 run with `--skip_update`; without it microSALT would simply re-download or re-index
  the reference, which is correct behaviour but not an error.
- V9 is a `--dry` run of the unmodified base case; its duration is the validation-time proxy.
- R4 seeds a corrupt cached download in a reference symlink farm; the real download path is not
  intercepted.
- R7 needs a sqlite database; it works on a scratch copy with the profile table emptied.
- R8 (unreachable resource) is not reproduced on stage and is reported as covered by unit tests.
- R2 downloads a real assembly into the configured genomes folder, and all R-tests create a
  synthetic project (`VAL5100` by default) in the stage database.
