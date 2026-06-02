# 2_STEM_IDENTIFICATION

Stage 2 of the BCI stem-reconstruction pipeline. This stage runs the DP solver
on cleaned ViewFullTable input from `BCI_stem_reconstruction/1_DATA_PREPARATION/` and produces the
reconstructed stem dataset used by `BCI_stem_reconstruction/3_PREPARE_R_TABLES/`.

## Scripts

### `1_main_cpp_chunk_bci.R`

Chunked DP driver.

- Loads `BCI_stem_reconstruction/DATA/PROCESSED/ViewFullTable_taper_corrected_growth_forms.rds`.
- Estimates per-species parameters and applies BCI-specific preprocessing.
- Runs `dp_global` on multi-stem tags in parallel chunks.
- Writes chunk outputs to `BASE_OUT_DIR/<run_timestamp>/`.
- Supports resuming interrupted runs.

### `2_merge_chunks_to_datatable.R`

Merges completed chunk outputs into final files.

- Reads Feather chunk outputs from `home_dir/run_code`.
- Converts them to temporary Parquet parts and merges them.
- Writes `merged_output.parquet` and `merged_output.rds` to
  `BCI_stem_reconstruction/DATA/<run_code>/`.

## Notes

- Refer to `BCI_stem_reconstruction/2_STEM_IDENTIFICATION/run_chunk_bci.md`
  for run and resume commands.

## Data flow

```text
DATA/PROCESSED/ViewFullTable_taper_corrected_growth_forms.rds
        │
        ▼
1_main_cpp_chunk_bci.R
        └──▶ BASE_OUT_DIR/<run_timestamp>/  (chunk Feather outputs)

2_merge_chunks_to_datatable.R
        └──▶ BCI_stem_reconstruction/DATA/<run_code>/merged_output.{parquet,rds}
```
