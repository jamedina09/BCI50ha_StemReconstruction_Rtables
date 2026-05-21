# 1_DATA_PREPARATION

Stage 1 of the BCI stem-reconstruction pipeline. Reads raw ForestGEO inputs
and produces the cleaned, harmonised tables that feed the stem-identification
DP engine.

## Scripts (run in order)

- `0_prepare_species_tables.R` — Build the per-species lookup tables
  (mnemonics, taxonomy, wood-density references) used by downstream stages.
- `1_prepare_viewfulltable.R.R` — Load the raw ViewFullTable, apply
  HOM/DBH taper corrections, fix obvious data-entry issues, and write the
  cleaned long-format census table consumed by `2_STEM_IDENTIFICATION/`.

## Subfolders

- `HELPER_FUNCTIONS/` — Shared utilities sourced by the scripts above
  (DBH unit conversions, HOM-correction routines, validation checks).

## Notes

`compare_tables1_and_tables2.R` is a local diagnostic (gitignored) used to
spot-check this stage against a reference build; it is not part of the
production pipeline.
