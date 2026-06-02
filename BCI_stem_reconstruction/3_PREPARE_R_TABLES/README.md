# 3_PREPARE_R_TABLES

Stage 3 of the BCI stem-reconstruction pipeline. This stage consumes the
reconstructed stem identities from `BCI_stem_reconstruction/2_STEM_IDENTIFICATION/` and generates the
ForestGEO-format census tables plus QC exports used for downstream analysis.

## Scripts (run in order)

- `1_prepare_posteriors_BCI.R` — Consolidates `_paths.feather` posterior files
  from a completed stage 2 run into `BCI_stem_reconstruction/DATA/POSTERIORS/posterior_sampled_paths.rds`.
- `2_create_R_tables_BCI.R` — Builds the final census tables and species table.
  It resolves encounter histories, applies DBH-aware broken-below rules,
  propagates terminal and prior states, corrects resurrection patterns,
  performs D/G remapping, imputes missing coordinates/dates, and exports
  ForestGEO-format `.Rdata` (and supporting `.csv`) tables.

## Outputs

- `BCI_stem_reconstruction/DATA/RTABLES/<site>.stemN.Rdata`
- `BCI_stem_reconstruction/DATA/RTABLES/<site>.spptable.rdata`
- `BCI_stem_reconstruction/DATA/POSTERIORS/posterior_sampled_paths.rds`
- QC exports in `BCI_stem_reconstruction/DATA/CHECKS/`

## Notes

- `1_prepare_posteriors_BCI.R` must be run before `2_create_R_tables_BCI.R`.
