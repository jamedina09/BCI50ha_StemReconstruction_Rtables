# 1_DATA_PREPARATION

Stage 1 of the BCI stem-reconstruction pipeline. Reads raw ForestGEO inputs
and produces the cleaned, harmonised tables that feed the stem-identification
DP engine.

## Scripts (run in order)

- `0_prepare_species_tables.R` — Build the per-species lookup tables
  used by downstream stages.
- `1_prepare_viewfulltable.R.R` — Load the raw ViewFullTable, apply
  HOM/DBH taper corrections, fix potential data-entry issues in DBH, and write the
  cleaned long-format census table consumed by `2_STEM_IDENTIFICATION/`.
  The user will receive the raw data, with corrected stemID identification,
  but not the 'corrected' potential data entre issues in DBH.

## Subfolders

- `HELPER_FUNCTIONS/` — Shared utilities sourced by the scripts above.
