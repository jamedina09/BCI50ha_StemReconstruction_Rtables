# 2_STEM_IDENTIFICATION/run_records

The record of every engine run of `1_main_cpp_chunk_bci.R`: how the run was
set up, how to run it again, and on which code, input and machine it ran. The
driver writes it; nothing here is edited by hand.

## Table of contents

- [Contents](#contents)
- [Files of a run](#files-of-a-run)
- [Rerun a run](#rerun-a-run)
- [Resume and stop](#resume-and-stop)
- [Settings](#settings)

## Contents

| File | What it does |
| --- | --- |
| `runs_index.csv` | One line per run (below) |
| `<run folder>/` | The record of one run, named as its run folder under `BASE_OUT_DIR` (`<date>_<time>_<label>`, label explained in `../run_chunk_bci.md`) |

`runs_index.csv` columns: `run_folder`, `status` (`running`, `resumed`,
`finished`, `finished, chunks missing or failed`, `stopped`), `started` (first
start), `finished`, `duration_h` (hours from the first start to the finish,
pauses included), `chunks_done`, `chunks_failed`, `n_tags`, `git_commit`,
`code_changed` (`yes` when `dp_global/` or `2_STEM_IDENTIFICATION/` had
uncommitted changes), `input_file`, `input_md5`, `test_run` (`yes` when a test
setting was set: the merge refuses the run), `changed_settings` (settings given
other values on the command line), the key settings `USE_MEASUREMENT_ERROR`,
`DBH_ROUND_CENSUSES`, `MAX_GROWTH_FIXED`, `MAX_SHRINK_FIXED`, `DP_MAX_STATES`,
`PROB_BIRTH_DEATH`, `PROB_EXPORT_RULE`, `BB_DBH_RULE`, `HARD_PINS`,
`PIN_STOP_AT_FAST_LINK`, `DP_POSTERIOR_LABEL_FREE`, `POSTERIOR_SAMPLES`,
`POSTERIOR_SAMPLE_SEED`, and `output_dir`.

## Files of a run

The driver writes them at the start of the run and copies them to the run
folder (`<BASE_OUT_DIR>/<run folder>/`), so the engine output carries its own
record.

| File | Content |
| --- | --- |
| `run_settings.txt` | A plain-language summary of the method (measurement-error model used or not, DBH classes, growth bounds, parameter estimation, engines and routing, matcher, pins, broken-below rule, posterior, test settings), the settings changed on the command line, and every setting with its value, its default and whether it was changed (`yes`: other value on the command line; `restated`: the default given on the command line; `derived`: set from another setting), grouped as in `../run_chunk_bci.md` |
| `rerun_command.sh` | The command of the run, the git commit to check out first, and the input file with its MD5 |
| `1_main_cpp_chunk_bci.R` | Copy of the driver as it was run |
| `uncommitted_changes.patch` | `git diff` of `dp_global/` and `2_STEM_IDENTIFICATION/` against the commit, when the code had uncommitted changes |
| `run_summary.txt` | Provenance, one entry per start, resume, finish or stop: time, command line, driver (MD5), git commit and branch, uncommitted files, R and package versions (`data.table`, `Rcpp`, `arrow`, `igraph`, `MASS`, `lpSolve`), machine and number of workers, input file (path, modification time, size, MD5), chunks finished and failed, number of tags, time since the first start, status |

The engine results (parameter sets, chunk outputs) stay in the run folder
(`run_parameters_full.txt`, `run_log.txt`).

## Rerun a run

From the project root:

```sh
git checkout <git commit of rerun_command.sh>
git apply BCI_stem_reconstruction/2_STEM_IDENTIFICATION/run_records/<run folder>/uncommitted_changes.patch   # only when the file exists
sh BCI_stem_reconstruction/2_STEM_IDENTIFICATION/run_records/<run folder>/rerun_command.sh
```

The input must be the same stage-1 file (`rerun_command.sh` gives its MD5).
The rerun writes a new run folder, with a new time stamp. The copy of the
driver runs the same way with the engine code of the same commit.

## Resume and stop

A resume (`--OUT_DIR_OVERRIDE=<run folder>`) finds the record of its run and
keeps the files written at the first start. It adds a `Resume` entry to
`run_summary.txt` and stops with `❌` when the input file (MD5) or a setting
outside "Run control and output" differs from the first start, and when the
chunk size differs (it sets which tags each chunk holds), so that every chunk
of a run shares them. Cores, output folders and resume flags may change.
It also notes whether the driver and the git commit are those of the first
start. A run that ends with an error gets a `Stop (error)` entry and the
status `stopped`.

## Settings

In `1_main_cpp_chunk_bci.R`: `WRITE_RUN_RECORD` (`TRUE`; `FALSE` writes no
record) and `RUN_RECORD_DIR` (this folder). Both can be set on the command
line, e.g. `--RUN_RECORD_DIR=/tmp/records` for a test run.
