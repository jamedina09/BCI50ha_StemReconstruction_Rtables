# 2_STEM_IDENTIFICATION

Stage 2 of the BCI stem-reconstruction pipeline. This stage runs the DP solver
on cleaned ViewFullTable input from `BCI_stem_reconstruction/1_DATA_PREPARATION/` and produces the
reconstructed stem dataset used by `BCI_stem_reconstruction/3_PREPARE_R_TABLES/`.

## Scripts

### `1_main_cpp_chunk_bci.R`

Chunked DP driver.

- Loads `BCI_stem_reconstruction/DATA/PROCESSED/ViewFullTable_single_vs_multiple_stem_tags.rds`.
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
- Joins stems that the engine split only because a diameter change fell
  outside its hard growth bounds: the point of measurement of the trunk moved
  (e.g. 1982 diameters taken around buttresses) or one diameter was recorded
  wrongly. An ended stem is joined to the stem of the same tree that starts in
  the next census when the new stem would be an impossible recruit (it starts
  above the recruit limit) and the ended stem is the only candidate; database
  StemIDs are not used (`apply_measurement_rejoin()`; rules in
  `dp_global/R/measurement_rejoin.R` and `dp_global/README.md`;
  `ReconstructionMethod = "measurement_rejoin"`). The joined stem keeps one of
  its two IDs, so the tree's IDs skip one number per join.
- Writes the final reconstructed stem table
  `DATA/PROCESSED/complete_dataset_final_with_reconstructed_stemids.rds` (the
  input of stage 3), the joins to `DATA/PROCESSED/measurement_rejoin_audit.csv`
  and the joined observation pairs, for the posterior samples, to
  `DATA/PROCESSED/measurement_rejoin_pairs.csv`.

## Notes

- Refer to `BCI_stem_reconstruction/2_STEM_IDENTIFICATION/run_chunk_bci.md`
  for run and resume commands.

## Data flow

```text
DATA/PROCESSED/ViewFullTable_single_vs_multiple_stem_tags.rds
        │
        ▼
1_main_cpp_chunk_bci.R
        └──▶ BASE_OUT_DIR/<run_timestamp>/  (chunk Feather outputs)

2_merge_chunks_to_datatable.R
        ├──▶ BCI_stem_reconstruction/DATA/<run_code>/merged_output.{parquet,rds}
        └──▶ BCI_stem_reconstruction/DATA/PROCESSED/
               complete_dataset_final_with_reconstructed_stemids.rds
               measurement_rejoin_audit.csv, measurement_rejoin_pairs.csv
```
