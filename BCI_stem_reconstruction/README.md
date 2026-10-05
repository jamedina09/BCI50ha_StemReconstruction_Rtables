# BCI Stem Reconstruction

This directory contains a four-stage workflow to reconstruct stem identities for the BCI 50-ha plot and produce analysis-ready ForestGEO-style outputs.

## What happens in this workflow

1. **Data preparation** (`1_DATA_PREPARATION/`)

- Cleans and harmonizes raw census/taxonomy inputs.
- Produces standardized species tables and a cleaned long-format ViewFullTable.

1. **Stem identification** (`2_STEM_IDENTIFICATION/`)

- Runs the chunked DP reconstruction engine on prepared inputs.
- Merges chunk outputs into the final reconstructed stem table, after joining
  trunks that the engine split only because their point of measurement moved
  or one diameter was recorded wrongly.

1. **R table creation** (`3_PREPARE_R_TABLES/`)

- Consolidates the posterior samples of the reconstruction.
- Converts the reconstructed stems into ForestGEO-compatible stem tables.
- Exports final `stemN` tables and QC/check files (the species table used in
  stage 4 is not written here; see `3_PREPARE_R_TABLES/README.md`).

1. **Example structure assessment** (`4_EXAMPLE_STRUCTURE_ASSESSMENT/`)

- Demonstrates downstream analyses (plot summaries, biomass and basal-area
  stocks/fluxes).
- Includes uncertainty propagation from posterior reconstruction paths.

## Main data flow

Raw BCI census and taxonomy data
→ cleaned/standardized inputs
→ reconstructed stem identities
→ final RTABLE products
→ example biomass and basal-area analyses

## Key output locations

- `DATA/PROCESSED/` prepared inputs from stage 1, and the final reconstructed
  stem table and rejoin files from stage 2
- `DATA/<run_code>/` merged reconstruction outputs from stage 2
- `DATA/POSTERIORS/`, `DATA/RTABLES/`, `DATA/CHECKS/` outputs from stage 3
- `4_EXAMPLE_STRUCTURE_ASSESSMENT/outputs/` example analysis products
