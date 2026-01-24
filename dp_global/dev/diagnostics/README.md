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
  Rscript generate_fake_book1.R dp_global/examples/diagnostics/Book1_fake.csv 20 "5,6,7" 42 "6,2" 5 0.5

- Run anchor plausibility grid on a specific CSV:
  Rscript diagnose_anchor_plausibility.R /path/to/Book1.csv "0.5,1,2,5,10"

- Run a generic parameter sweep (cost-only) on a dataset:
  Rscript sweep_param_generic.R Bio_Recruit_MaxDBH_unit "0.5,1,2,5,10" cost /path/to/Book1.csv

- Compute per-track cost time series (probe=6) on a dataset:
  Rscript track_cost_timeseries.R 6 Bio_Recruit_MaxDBH_unit "0.5,1,2,5,10" /path/to/Book1.csv

- Generate a grid of simulated datasets and run full diagnostics on each:
  Rscript simulate_and_run_diagnostics.R
  (Outputs are placed under `dp_global/examples/diagnostics/sims/`)

Files created:
- `utils.R` : helper functions
- `generate_fake_book1.R` : parameterized dataset generator
- `diagnose_anchor_plausibility.R` : per-anchor plausibility grid
- `sweep_param_generic.R` : general parameter sweep
- `track_cost_timeseries.R` : per-track cost time series
- `run_synthetic_diagnostics.R` : convenience runner (single run)
- `simulate_and_run_diagnostics.R` : orchestrator to run diagnostics on simulated dataset grid

