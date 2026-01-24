# Error propagation helpers

## Prerequisites

R (packages: `data.table`; `ggplot2` recommended for optional plotting)

Overview
--------
This directory provides two focused, complementary scripts to work with DP posterior reconstructions:

- `process_posteriors.R` — the canonical function library that implements parsing, loading, attaching, and expansion utilities; and
- `attach_paths_to_output_run.R` — a small, interactive runner that demonstrates typical workflows and provides the convenience function `run_attach_and_expand()` for RStudio users.

The code here is built around a single deterministic matching key: ObsRowID. All reconstructions attached to the main per-row output are matched by `obs_row_id`.

What each script does
---------------------
- `process_posteriors.R` (functions):
  - `parse_recon(s)`
    - Parse a `recon` string of the form `ObsRowID:ReconstructedStemID;...` into a `data.table(ObsRowID, ReconstructedStemID)`.
    - Returns an empty table for missing or blank inputs.
  - `load_posterior_paths(paths_file)`
    - Read a per-path summary (`*_paths.csv` / feather / rds), validate required columns (`path_sig`, `path_prob`, `recon`), parse `recon` into a `recon_parsed` list-column, normalize `path_prob`, and sort by probability.
  - `attach_paths_to_output(paths, out, which = c("map","top_n","indices","sample"), ...)`
    - Attach selected paths to the main per-row reconstruction `out` (path or `data.table`) by matching `ObsRowID`.
    - Adds `DP_PathSig_k` and `DP_ReconstructedStemID_k` columns for attached paths and optionally writes a `_with_paths.csv` when `write_out = TRUE`.
  - `expand_draws(summary_dt, paths_dt, N)` and `aggregate_draws(res_dt)`
    - Sample `N` draws from the `_summary.csv` (using `sample_prob`), map each sampled `path_sig` to its `recon` string, expand into a long (Draw, ObsRowID, ReconstructedStemID) table, and aggregate per-observation frequencies into probabilities.
  - `check_map_in_paths(paths_dt, out_dt)`
    - Diagnostic helper to compare the decoded MAP joint path from `out_dt` against sampled unique paths and report the best partial match.

- `attach_paths_to_output_run.R` (runner):
  - Define the main paths here.
  - `run_attach_and_expand(paths_file, out_file, attach_n = 2, write_out = TRUE, run_expand_draws = FALSE, expand_N = 1000, expand_out_csv)` is the interactive entry point that:
    - loads the per-path summary (via `load_posterior_paths()`),
    - attaches the top-n paths to the main output using `attach_paths_to_output()` and writes `_with_paths.csv` if requested, and
    - optionally expands `N` draws and writes aggregated per-`ObsRowID` probabilities to `expand_out_csv` (when `run_expand_draws = TRUE`).
  - The runner is intentionally minimal: it calls canonical functions in `process_posteriors.R` and performs basic diagnostics and optional aggregation.

Requirements and expectations
-----------------------------
- All attachments rely on explicit ObsRowID mapping in the `recon` strings (format `ObsRowID:ReconstructedStemID;...`). Ensure your posterior-sampling step preserves ObsRowID in the exported `recon` field for deterministic matching.
- The runner is designed for interactive use (RStudio or sourcing); it does not parse CLI arguments.
- The repository no longer writes the full long per-draw posterior table by default; this workflow uses compact `paths` + `summary` files and supports on-demand expansion to reduce I/O and runtime. Use `run_expand_draws = TRUE` to generate aggregated per-ObsRowID probabilities.

Examples
--------
- Attach top 2 paths (interactive):

```r
# Source the runner (adds run_attach_and_expand)
source("dp_global/R/error_propagation/attach_paths_to_output_run.R")
# Attach and write _with_paths.csv (no draw expansion)
run_attach_and_expand(attach_n = 2, write_out = TRUE, run_expand_draws = FALSE)
```

- Expand 1000 draws and write aggregated per-ObsRowID probabilities:

```r
# Expand draws (this may take time depending on N)
run_attach_and_expand(run_expand_draws = TRUE, expand_N = 1000, expand_out_csv = "dp_global/output/expanded_draws_sample_1000.csv")
```

Files in this directory
-----------------------
- `process_posteriors.R` — canonical implementation of parsing, loading, attaching, and expansion helpers.
- `attach_paths_to_output_run.R` — interactive runner and demonstration (`run_attach_and_expand()`).

If you'd like, I can add a small example aggregated file (sample output snippet) to the README or add a simple unit test to assert ObsRowID enforcement. Let me know which you'd prefer.

Files in this directory
-----------------------
- `process_posteriors.R` — canonical implementation with parsing, loading, attaching, and expansion helpers.
- `attach_paths_to_output_run.R` — interactive runner demonstrating attachment and optional draw expansion (`run_attach_and_expand()`).

Note: `reconstruct_and_propagate.R` was renamed to `process_posteriors.R` and the old compatibility shim has been removed — update any scripts that source the old filename.