# 3_PREPARE_R_TABLES

Stage 3 (final) of the BCI stem-reconstruction pipeline. Takes the
reconstructed stem identities from `2_STEM_IDENTIFICATION/` and emits the
ForestGEO-format per-census R tables (`.Rdata`) along with the posterior /
diagnostic exports.

## Scripts (run in order)

- `1_prepare_posteriors_BCI.R` — Read the per-chunk arrow/feather outputs
  from stage 2, assemble the full posterior dataset, and write the
  reconstructed-StemID table consumed by the table builder.
- `2_create_R_tables_BCI.R` — Main table builder. Resolves encounter
  histories (alive / dead / gone / prior / broken-below with DBH-aware
  rules), runs the propagation and resurrection fixes, imputes missing
  coordinates and dates per tree/quadrat, and writes one ForestGEO R table
  per census plus quality-control reports.
