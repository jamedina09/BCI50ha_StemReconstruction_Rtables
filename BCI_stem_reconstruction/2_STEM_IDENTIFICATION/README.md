# 2_STEM_IDENTIFICATION

Stage 2 of the BCI stem-reconstruction pipeline. This stage runs the
`dp_global` engines on the cleaned ViewFullTable from
`BCI_stem_reconstruction/1_DATA_PREPARATION/` and produces the reconstructed
stem dataset used by `BCI_stem_reconstruction/3_PREPARE_R_TABLES/`. Both
scripts are run from the project root.

## Scripts

### `1_main_cpp_chunk_bci.R`

Chunked DP driver.

- Loads `BCI_stem_reconstruction/DATA/PROCESSED/ViewFullTable_single_vs_multiple_stem_tags.rds`.
- Checks at load that every growth form and species code named in the routing
  settings (`DP_FALLBACK_GROWTH_FORMS`, `NON_TAPER_CORRECTED_GROWTH_FORMS`,
  `PROB_SPECIES`) occurs in the data; a value that does not (e.g. `palms` for
  the label `palm`) stops the run with `❌ CHECK FAILED`.
- Estimates the biological parameters from the 2010–2023 stems of all tags:
  one set per species with enough data and pooled sets per growth form for
  the others (all trees, all shrubs, palms + tree ferns, strangler figs). The
  growth SD is fitted as a constant when it falls with size, and the
  recruitment rate is new stems per established tree per year
  (`RECRUIT_RATE_UNIT`; see `dp_global/README.md`, *Parameter Estimation*). A
  species gets its own parameter set when its diameters cover the size range
  either from smallest to largest or over their middle 95%
  (`COVERAGE_RULE = "union"`), so a single out-of-range diameter does not
  send a well-sampled species to a pooled set.
- Pins the rows whose identity the database fixes (`TrueStemID`, Steps 1–3 in
  the script: stem tags, 2010+ records, measured break rows and the rows of
  the same StemID); the engines resolve the other rows.
- Runs `dp_global` on multi-stem tags in parallel chunks: the exact DP, or the
  probabilistic matcher for palm clumps (*Oenocarpus mapora*, *Bactris major*),
  strangler figs and trees too complex for the DP. The matcher lets every stem
  continue, die or be recruited in every census pair (`PROB_BIRTH_DEATH`).
- Applies the post-engine helper chain to every chunk (pin-track rejoin,
  backfills, terminal-to-host, broken-below invariants, renumbering of each
  tag's stems to 1..N, posterior path files; see `dp_global/README.md`).
- Writes one Feather file per chunk and the posterior path files to
  `BASE_OUT_DIR/<run folder>/` (the folder name starts with the run's time
  stamp).
- Supports resuming interrupted runs.

### `2_merge_chunks_to_datatable.R`

Merges completed chunk outputs into final files.

- Reads Feather chunk outputs from `home_dir/run_code` (both set at the top
  of the script; the run must have finished).
- Converts them to temporary Parquet parts and merges them.
- Writes `merged_output.parquet` and `merged_output.rds` to
  `BCI_stem_reconstruction/DATA/<run_code>/`.
- Compares the merged table with the stage-1 table, restores the recorded DBH
  (mm) and adds the single-stem tags, whose one stem gets
  `ReconstructedStemID = 1` (`ReconstructionMethod =
  "single_stem_tag_no_reconstructed"`).
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

## Other files

- `run_chunk_bci.md` — run and resume commands for `1_main_cpp_chunk_bci.R`.
- `tests/test_measurement_rejoin.R` — tests of the measurement rejoin
  (`testthat::test_file("BCI_stem_reconstruction/2_STEM_IDENTIFICATION/tests/test_measurement_rejoin.R")`).
- `comparissons/` (not under version control) —
  `compare_reconstruction_vs_dryad.R` and `dryad_approach_weaknesses.qmd`,
  which compare the reconstruction with the published BCI stem tables.

## Requirements

R packages: `data.table`, `parallel`, `arrow` (Feather / Parquet), and those
loaded by `dp_global/R/dp_global_main.R` (`Rcpp`, `igraph`, `MASS`,
`lpSolve`); `testthat` for the tests.

## Data flow

```text
DATA/PROCESSED/ViewFullTable_single_vs_multiple_stem_tags.rds
        │
        ▼
1_main_cpp_chunk_bci.R
        └──▶ BASE_OUT_DIR/<run folder>/  (chunk Feather outputs, posteriors/)

2_merge_chunks_to_datatable.R
        ├──▶ BCI_stem_reconstruction/DATA/<run_code>/merged_output.{parquet,rds}
        └──▶ BCI_stem_reconstruction/DATA/PROCESSED/
               complete_dataset_final_with_reconstructed_stemids.rds
               measurement_rejoin_audit.csv, measurement_rejoin_pairs.csv
```
