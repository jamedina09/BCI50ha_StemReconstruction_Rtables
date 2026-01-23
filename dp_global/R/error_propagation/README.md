Error propagation helpers — README

Overview
--------

This directory contains helper utilities for loading posterior path summaries
produced by the DP (`*_paths.csv`, Feather, or RDS), applying full-path
reconstructions to observed DBH data, and computing per-path and per-stem
growth statistics for uncertainty propagation.

Files
-----

- `reconstruct_and_propagate.R` — Main helper script. Exposes these functions:
  - `load_posterior_paths(paths_file)` — load/parse paths file and return a
    data.table with a `recon_parsed` list-column.
  - `apply_paths_compute_growth(paths_dt, obs_dt, obs_row_id_col = NULL, strict = TRUE, plot = FALSE, plot_draws = 1000)` —
    compute per-path growth metrics and optional plot of sampled mean-growth. Accepts an optional `obs_row_id_col` (default `NULL`); when provided (or when a column named `obs_row_id` is present) matching can be done using stable observation IDs instead of relying on within-census order.
  - `sample_posterior_paths(paths_dt, n)` — sample path_sig indices according to path probability.
  - `sample_apply_growth(paths_dt, obs_dt, n)` — Monte Carlo draws of mean-growth per draw.

Design & behavior
-----------------

- The `paths` file must include columns `path_sig`, `path_prob`, and `recon`.
  The `recon` column should be a semicolon-separated list of `CensusID:ReconstructedStemID` pairs.
- Matching between recon entries and observations is done by within-census
  ordering (an `obs_pos` index). If an `obs_row_id` column is present in your
  observations (and the DP exported recon strings using `obs_row_id`), the
  helper will prefer matching by `obs_row_id` which is robust to row reordering;
  otherwise it falls back to `obs_pos`. This prevents ambiguous cartesian joins.
- If recon entries and observed counts mismatch for a census then:
  - When `strict = TRUE` (default), `apply_paths_compute_growth()` will abort with an error and show representative examples.
  - When `strict = FALSE`, the function emits a single aggregated warning summarizing examples and continues (this is useful for exploratory workflows).

Plotting
--------

- When `plot = TRUE`, the function samples mean-growth values according to `path_prob` and generates a histogram of sampled mean-growths.
  - If `ggplot2` is installed, a ggplot object is returned in `res$plot`.
  - Otherwise, a base R histogram is drawn and `res$plot` is `NULL`.
- The plot includes a vertical dashed line at the expected mean and annotated mean/sd.

Example: quick usage
--------------------

# Interactive example (R)

```r
library(data.table)
source("dp_global/R/error_propagation/reconstruct_and_propagate.R")
paths_dt <- load_posterior_paths("dp_global/output/.../posteriors/tag_20_posterior_samples_..._paths.csv")
obs <- data.table::fread("data_simulation/data/simulated_data_1.csv")[Tag == 20]
# Non-strict mode, draw a plot with 1000 sampled mean-growth values
res <- apply_paths_compute_growth(paths_dt, obs, strict = FALSE, plot = TRUE, plot_draws = 1000)
print(res$expected_mean_growth)
if (!is.null(res$plot)) print(res$plot) # ggplot object
```

# Monte-Carlo sampling route

```r
draws <- sample_apply_growth(paths_dt, obs, n = 1000)
hist(draws$mean_growth, breaks = 50, main = "Monte Carlo mean-growth draws")
```

Notes
-----

- The helpers are intentionally conservative: mismatch detection exists to
  ensure users notice when reconstructions do not align with the observed
  dataset (e.g., when you use a `paths` file for a different Tag or subset).
- If you want non-blocking behavior, use `strict = FALSE`.

Contact
-------

If you discover a workflow that produces unexpected results, please open an issue
or add unit tests demonstrating the case so we can harden the helpers further.
