# 1_DATA_PREPARATION

This folder contains the first stage of the BCI stem reconstruction pipeline.
Its purpose is to convert raw ForestGEO exports into cleaned, harmonized input
tables for the downstream stem identification and reconstruction stages.

## What this stage does

- Loads raw ForestGEO census exports and taxonomy inputs.
- Builds standardized species lookup tables used by later pipeline stages.
- Cleans and validates DBH / HOM data and applies taper corrections.
- Writes a long-format cleaned ViewFullTable for use by
  `BCI_stem_reconstruction/2_STEM_IDENTIFICATION/`.

## Scripts (run in order)

1. `0_prepare_species_tables.R`
   - Builds the species lookup tables used by the reconstruction pipeline.
   - Produces standardized species metadata for later tables.

2. `1_prepare_viewfulltable.R.R`
   - Loads raw `ViewFullTable` data.
   - Applies HOM/DBH taper corrections and fixes common data-entry issues.
   - Harmonizes census records and outputs the cleaned long-format table.
   - Prepares the data for stem reconstruction in `2_STEM_IDENTIFICATION/`.

## Inputs

The scripts expect raw ForestGEO data exports and taxonomy reference files
from the BCI dataset. Specific input locations are defined in each script,
but typically include:

- `ViewFullTable` raw census export
- species taxonomy / lookup files

## Outputs

This stage produces:

- cleaned species lookup tables
- a standardized `ViewFullTable` long-format table
- any intermediate validation or helper outputs used by later stages

## Directory structure

- `HELPER_FUNCTIONS/` — shared utility functions and helpers used by the
  stage 1 scripts.

## Notes

- This stage does not perform final stem-level reconstruction; it only
  prepares the cleaned, harmonized inputs for the DP stem-identification
  stage.
