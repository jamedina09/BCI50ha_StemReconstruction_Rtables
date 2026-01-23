Minimal error-propagation helpers (README)
==========================================

Overview
--------
This directory now provides a small, focused toolkit to:

- load DP posterior `paths` summaries (per-path aggregated files), and
- attach selected posterior path mappings to an existing *main reconstruction*
  CSV by matching `ObsRowID`.

The implementation was simplified to contain only the minimal helpers required
by downstream convenience scripts (e.g., `attach_paths_to_output_run.R`).

Available functions
-------------------
- `parse_recon(s)`
  - Parse a compact `recon` string produced by the DP. Accepts:
    - explicit pairs: `"1:8;2:3;3:4"` (left value is `ObsRowID` or `CensusID`), or
    - compact path: `"7-4-8-..."` (hyphen-separated stem IDs assigned in census order).
  - Returns a `data.table` with columns `CensusID` (integer) and `ReconstructedStemID` (integer).

- `load_posterior_paths(paths_file)`
  - Read a per-path summary file (`*_paths.csv`, feather, or rds) and return a
    `data.table` with a `recon_parsed` list-column (parsed by `parse_recon`).
  - Validates presence of `path_sig`, `path_prob`, and `recon`, and normalizes `path_prob`.

- `attach_paths_to_output(paths, out, which = c("map","top_n","indices","sample"), ...)`
  - Attach selected paths as new columns to `out` (file path or `data.table`) by matching `ObsRowID`.
  - Adds columns: `DP_PathSig_k` and `DP_ReconstructedStemID_k` for each selected path `k`.
  - If `out` lacks `obs_row_id`, one is created (sequence of row numbers) before matching.
  - Modes:
    - `which = "map"` attach the single MAP path (highest `path_prob`)
    - `which = "top_n"` attach the top-`n` paths
    - `which = "indices"` attach specific path indices (pass `indices=c(...)`)
    - `which = "sample"` sample `n` paths according to `path_prob`
  - Optionally write a `_with_paths.csv` output when `write_out = TRUE`.

Notes & rationale
-----------------
- The DP posterior outputs often encode reconstructions using `ObsRowID` as the
  left-side identifier in the `recon` string. Matching by `ObsRowID` is robust
  and preferred; this helper will create an `obs_row_id` column in the main
  reconstruction CSV if one does not already exist.
- The helper intentionally leaves rows that have no mapping as `NA` for the
  `DP_ReconstructedStemID_k` column; this reflects the DP's actual sampling
  output (only observed DBH rows are included in sampled full-reconstructions).

Example: attach top 2 paths to the main reconstruction CSV
---------------------------------------------------------
Run the convenience script (example included in this directory):

```sh
Rscript dp_global/R/error_propagation/attach_paths_to_output_run.R
```

Or call the helper interactively in R:

```r
library(data.table)
source("dp_global/R/error_propagation/reconstruct_and_propagate.R")
paths_dt <- load_posterior_paths("dp_global/output/.../posteriors/tag_<TAG>_posterior_samples_<TS>_paths.csv")
out_dt <- attach_paths_to_output(paths_dt, "dp_global/output/.../stem_reconstruction_dp_global_rcpp.csv", which = "top_n", n = 2, write_out = TRUE)
```