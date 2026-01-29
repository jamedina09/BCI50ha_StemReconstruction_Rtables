# Diagnostics scaffolding

This folder contains consolidated diagnostic scripts and utilities for:

- sweeping biological and DP parameters
- computing per-track cost components and p(recruit)
- generating synthetic datasets for testing
- running reproducible diagnostic pipelines

Usage examples:

- Generate fake data and run example diagnostics:
  Rscript run_synthetic_diagnostics.R

- Generate parameterized synthetic datasets (single):
  Rscript generate_synthetic_data.R dp_global/dev/diagnostics/diagnostics_dataset.csv 20 "5,6,7" 42 "6,2" 5 0.5

- Run anchor plausibility grid on a specific CSV (diagnostics uses only provided datasets; Book1 is not used):
  Rscript diagnose_anchor_plausibility.R /path/to/diagnostics_dataset.csv "0.5,1,2,5,10"

- Run a generic parameter sweep (cost-only) on a dataset:
  Rscript sweep_param_generic.R Bio_Recruit_MaxDBH_unit "0.5,1,2,5,10" cost /path/to/diagnostics_dataset.csv

- Compute per-track cost time series (probe=6) on a dataset:
  Rscript track_cost_timeseries.R 6 Bio_Recruit_MaxDBH_unit "0.5,1,2,5,10" /path/to/diagnostics_dataset.csv

Files included:
- `utils.R` : helper functions (load_dataset, run_dp_on_dt, compute_cost_components, p_recruit_from_costs, save_report)
- `generate_synthetic_data.R` : parameterized dataset generator
- `diagnose_anchor_plausibility.R` : per-anchor plausibility grid
- `sweep_param_generic.R` : general parameter sweep
- `track_cost_timeseries.R` : per-track cost time series
- `run_synthetic_diagnostics.R` : convenience runner (single run)

Note: ad-hoc or one-off diagnostic/debug scripts were removed during cleanup. Use the scripts listed above for reproducible diagnostics and regression checks.

