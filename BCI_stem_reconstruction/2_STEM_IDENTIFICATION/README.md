# 2_STEM_IDENTIFICATION

Stage 2 of the BCI stem-reconstruction pipeline. Runs the dp_global DP engine
on the cleaned ViewFullTable from `1_DATA_PREPARATION/` and assembles the
per-stem reconstructed identities used downstream in `3_PREPARE_R_TABLES/`.

## Scripts (run in order)

### `1_main_cpp_chunk_bci.R` — Chunked DP driver

Loads `DATA/PROCESSED/ViewFullTable_taper_corrected_growth_forms.rds`, estimates
per-species biological parameters, applies a BCI-specific TrueStemID pre-propagation
(Steps 1–3), then runs the dp_global DP engine across all multi-stem tags in
parallel chunks.

Key behaviours:

- Tags are split into single-stem (bypass, `ReconstructedStemID = StemID`) and
  multi-stem (DP reconstruction) subsets before the engine runs.
- Output is written incrementally: one Feather file per chunk to
  `BASE_OUT_DIR/<run_timestamp>/`.
- Post-engine helpers run per chunk in order:
  `maybe_add_posterior_bins` → `apply_carried_terminal_backfill` →
  `apply_orphan_stem_backfill` → `apply_broken_below_invariants` →
  `renumber_engine_minted_ids` → `finalize_posterior_paths`.
- Resumes partial runs automatically (`DP_CHUNK_RESUME=TRUE`); a chunk is
  skipped if its `_done.txt` completion marker already exists.
- All key parameters are overridable via `--KEY=VALUE` CLI flags. See the
  `CLI_REFERENCE` table inside the script for the full list.

**How to run:** see [`run_chunk_bci.txt`](run_chunk_bci.txt) for the production
command line and resume example.

### `2_merge_chunks_to_datatable.R` — Merge and assemble

Reads all Feather chunk files from a completed DP run directory, batch-converts
them to intermediate Parquet parts to control peak memory, merges into a single
dataset, reattaches single-stem records from the raw table, validates row counts
and column agreement, and saves:

- `DATA/<run_code>/merged_output.parquet` — merged multi-stem DP output
- `DATA/<run_code>/merged_output.rds` — same as RDS
- `DATA/PROCESSED/complete_dataset_with_reconstructed_stemids.rds` — full
  dataset (single-stem + multi-stem) ready for `3_PREPARE_R_TABLES/`

Set `home_dir` and `run_code` at the top of the script to point at the desired
run directory before executing.

## Other files

- [`run_chunk_bci.txt`](run_chunk_bci.txt) — Production command lines for
  starting a new run and resuming a partial run, with annotated flag descriptions.

## Data flow

```
DATA/PROCESSED/ViewFullTable_taper_corrected_growth_forms.rds
        │
        ▼
1_main_cpp_chunk_bci.R   →   BASE_OUT_DIR/<run_ts>/*_chunk_NNN.feather
                              BASE_OUT_DIR/<run_ts>/posteriors/tag_*_paths.feather
        │
        ▼
2_merge_chunks_to_datatable.R
        │
        ├──▶  DATA/<run_code>/merged_output.{parquet,rds}
        └──▶  DATA/PROCESSED/complete_dataset_with_reconstructed_stemids.rds
```
