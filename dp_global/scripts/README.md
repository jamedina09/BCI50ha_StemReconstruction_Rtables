# dp_global/scripts

This document describes the current behavior of the `dp_global` driver scripts:
- `dp_global/scripts/main_cpp.R` — the interactive/CLI driver for single-tag or targeted runs.
- `dp_global/scripts/main_cpp_chunk.R` — the chunked driver optimized for large runs.
- `dp_global/scripts/main_cpp_bci.R` — the BCI debug driver for single-tag runs on BCI census data.
- `dp_global/scripts/basal_area_uncertainty.R` — posterior-based basal area uncertainty quantification.

Both `main_cpp.R` and `main_cpp_chunk.R` accept command-line overrides of defaults using `--KEY=VALUE` flags. Keys are case-insensitive and may use `-` or `_` as separators. `main_cpp_bci.R` inherits the same CLI interface from `main_cpp.R`.

---

## Overview ✨

`main_cpp.R` is the central driver for the `dp_global` pipeline. It:

- Loads input tree census data (CSV via `input_file`).
- Ensures species information is present and optionally forces a single species label.
- Estimates biological parameters per species (via `estimate_bio_pars`).
- Runs Dynamic Programming (DP) reconstruction logic (via `run_dp_one_group` which calls the Rcpp-based DP functions).
- Optionally runs sensitivity sweeps and realism checks.
- Writes outputs (CSV, RDS, PDF) to an automatically created output directory.

You can run either script directly via Rscript:

- Single-tag or small dataset: `Rscript dp_global/scripts/main_cpp.R --INPUT_FILE=... --RUN_ALL_TAGS=TRUE`
- Large dataset (chunked, recommended): `Rscript dp_global/scripts/main_cpp_chunk.R --INPUT_FILE=...`

(When running from the `dp_global` folder: `Rscript scripts/main_cpp.R`.)

---

## Important behavior notes ⚠️

- `main_cpp.R` defaults to single-tag mode (`RUN_ALL_TAGS = FALSE`). Use `--RUN_ALL_TAGS=TRUE` to run across all tags.
- For large datasets, prefer `main_cpp_chunk.R` which divides groups (Tag + species) into chunks (`DP_CHUNK_SIZE`) and writes chunk outputs incrementally to disk to keep peak RAM low.
- When chunking, the combined `out` is not assembled in memory. Instead:
  - Each chunk is appended to `stem_reconstruction_dp_global_rcpp.csv` (if `WRITE_DP_CSV=TRUE`),
  - Each chunk is saved as `stem_reconstruction_dp_global_rcpp_chunk_###.rds` (if `WRITE_DP_RDS=TRUE`) and serves as a resume/completion marker.
- `maybe_add_posterior_bins()` is applied per-chunk when chunking; in non-chunked runs it is applied after assembling the full `out`.
- The chunked driver disables or comments out features that are not applicable to chunked runs (for example: full-run sensitivity sweeps and realism report generation).

---

## Key CLI flags & defaults 🔧

Common flags used by both drivers (case-insensitive, but use capital letters to keep it clean):

- `INPUT_FILE` — default: `data_simulation/data/simulated_data_1.csv`.
- `FORCE_ONE_SPECIES_PARAMETERS` — default: `TRUE`.
- `FORCED_SPECIES_LABEL` — default: `"all"` (used when forcing one species label).
- `SPECIES_COL` — default: `NULL` (auto-detected if not set).
- `USE_MEASUREMENT_ERROR` — default: `TRUE`.

DP / reconstruction option
- `DP_MODE` — default: `"marginals+bins"`. Allowed: 'none' 'marginals' 'marginals+bins' 'map'
- `WHICH_TAG` — character; used for single-tag runs (relevant to `main_cpp.R`). Must match the `Tag` column exactly (e.g., `"084555"` preserves leading zeros). The chunked runner processes groups (`Tag`, `species`) and does not rely on `WHICH_TAG`.
- `ANCHOR_START_CENSUS` — default: `7L`.
- `ALLOW_PROVISIONAL_DP_ANCHOR` — default: `TRUE` — when `TRUE` the DP can assign provisional anchor IDs at the last observed DBH census if the requested anchor census lacks `TrueStemID` but has DBH; set to `FALSE` to require an explicit anchored census or to fall back to the probabilistic matcher.
- `DP_VERBOSE` — default: `TRUE`.
- `DP_POSTERIOR_TOP_K` — default: `2L` - top=k posterior reconstructions to track. 
- `DP_MAX_TRACKS` — default: `NULL` (auto-computed per data when `NULL`) — optionally force max tracks; if `NULL` auto-computed.
- `DP_MAX_STATES` — default: `40000L`. Maximum number of injective assignment states allowed per census. Also controls the inter-census transition limit (`max_edges = max_states²`). The number of states at a census with $n$ observed stems and $K$ tracks is $P(K,n) = K!/(K-n)!$. When any census exceeds `max_states`, or when the cross-product of states between two adjacent censuses exceeds `max_states²`, the solver falls back to the probabilistic matcher. **Practical limits with default 40,000:** DP handles up to 6 observed stems per census exactly (P(7,6) = 5,040 < 40,000); 7+ stems trigger fallback (P(8,7) = 40,320 > 40,000). With `DP_MAX_STATES = 1,000`: max 5 stems. With `DP_MAX_STATES = 20,000`: max 6 stems. See `dp_global/README.md` for detailed tables and how to choose a value.
- `DP_SLACK_TRACKS` — default: `1L` - slack (additional) tracks allowed for DP.
- `DP_SLACK_REQUIRE_ANCHOR_RECRUITABLE` — default: `TRUE` - require anchor to be recruitable (if DBH less than max recruitment size) before granting slack.
- `DP_SLACK_REQUIRE_ANCHOR_EPS` — default: `1e-6`.
- `PROB_N_SAMPLES` — default: `200L`. Number of Gumbel-noise stochastic samples drawn by the probabilistic greedy matching fallback. Higher values produce more accurate posterior estimates at the cost of computation time.
- `PROB_LOOKAHEAD_WEIGHT` — default: `0.5`. Weight for sequential backward conditioning in the probabilistic matcher. When > 0 and $K \geq 4$, the cost matrix for each census pair is conditioned on the already-resolved forward assignment, improving path continuity. Set to `0` to disable conditioning.
- `USE_BIO_HARD_SHRINK_IN_PROB` — default: `TRUE`. Controls two hard constraints in the probabilistic matcher: (1) the bio shrink gate in pairwise edge construction (`g < Bio_Max_Shrink` → edge log-likelihood set to `-Inf`), and (2) the Layer 2 ME cumulative-shrinkage repair check in `repair_stitched_growth_violations()`. When `FALSE`, both are disabled — edges with growth below `Bio_Max_Shrink` are allowed (penalised by the soft `k_shrink` quadratic term only) and the ME cumulative check does not sever trajectories. Use `FALSE` to allow confirmed large-shrinkage events (e.g., storm damage) without forcing the matcher to split a continuous stem. **Note:** The exact DP solver is unaffected — it always applies `Bio_Max_Shrink` as a hard pruning bound regardless of this flag.
- `USE_BIO_HARD_GROWTH_IN_PROB` — default: `TRUE`. Controls the bio growth gate in pairwise edge construction in the probabilistic matcher (`g > Bio_Max_Growth_Bio` → edge log-likelihood set to `-Inf`). When `FALSE`, edges exceeding the biological growth maximum are allowed (penalised by the soft `k_growth` quadratic term only). The exact DP solver is unaffected by this flag.
- `PIN_TRUESTEMID` — default: `TRUE`. When `TRUE`, observations with a known `TrueStemID` at non-anchor censuses are pinned to their field-observed identity track in the probabilistic matcher, reducing effective state space and preventing re-identification of labelled stems.
- `DP_FALLBACK_GROWTH_FORMS` — default: `character(0)`; comma- or
  semicolon-separated list of values in the `growth_form` column that should
  trigger an immediate probabilistic fallback and prevent the DP solver from running
  on that tag. The driver automatically splits the string into a vector.
- `NON_TAPER_CORRECTED_GROWTH_FORMS` — default: `c("palm", "strangler_fig", "tree_fern")`;
  growth forms whose DBH is NOT taper-corrected. These forms exhibit both real
  biological DBH growth (palms: 1–3 cm/yr; strangler figs: variable as they
  encircle hosts) and apparent DBH shifts when the measurement height (HOM)
  changes between censuses. Because their trunk geometry does not allow taper
  correction, any HOM change produces an uncompensated apparent DBH change.
  The DP solver replaces the general pruning bounds with wide base bounds
  (see below) and optionally applies HOM-proportional widening on top.
- `NON_TAPER_CORRECTED_PRUNE_MIN_GROWTH` — default: `-0.625` cm/year
  (`1.25 × MAX_SHRINK_FIXED`); lower annual growth bound used during pruning
  for non-taper-corrected growth forms. Replaces the general effective pruning
  minimum when the tag's `growth_form` matches `NON_TAPER_CORRECTED_GROWTH_FORMS`.
  Note: in the default driver configuration the general `prune_min_growth` is
  also `1.25 × MAX_SHRINK_FIXED`, so this override is only active when
  `prune_use_bio_bounds = TRUE` or narrower general bounds are set.
- `NON_TAPER_CORRECTED_PRUNE_MAX_GROWTH` — default: `6.25` cm/year
  (`1.25 × MAX_GROWTH_FIXED`); upper annual growth bound for non-taper-corrected
  growth forms during pruning. Same interplay with the general bounds as above.
- `HOM_TOLERANCE_SCALE` — default: `2.0` (cm/year per meter of HOM deviation);
  when a non-taper-corrected tag has a `hom` (or `HOM`) column in its data,
  the prune bounds are widened for each census pair by
  `hom_tolerance_scale × max(|HOM − 1.3|) / interval_years`.
  This is where the effective differentiation between regular trees and
  non-taper-corrected forms occurs in the default configuration. NA HOM values
  are treated as 1.3 (zero contribution). Set to `0` to disable
  HOM widening while keeping the wide base bounds.

**Anchor scoping and post-anchor preservation:** If observations exist after the requested `ANCHOR_START_CENSUS`, the DP is scoped to censuses <= `ANCHOR_START_CENSUS` and post-anchor rows are preserved and appended to the output. Post-anchor rows with non-NA `DBH` and a non-NA `TrueStemID` are set to `ReconstructedStemID = TrueStemID` and `ReconstructionMethod = "given"`. Remaining post-anchor rows (missing `TrueStemID` or `DBH`) receive `ReconstructionMethod = "none_after_anchor"`. If scoping removes all pre-anchor observations, the original rows are returned and anchor-census `TrueStemID` values are treated as `"given"` while other rows are labeled `"none_after_anchor"`. **Note:** Pre-anchor rows (censuses before the anchor) with non-NA `TrueStemID` are processed by the solver and can receive different identity assignments based on biological likelihood — however, when `PIN_TRUESTEMID=TRUE` (default), they are constrained to their field-observed identity.

**Provisional anchor behavior:** When a requested anchor census lacks `TrueStemID` but contains DBH observations and `ALLOW_PROVISIONAL_DP_ANCHOR=TRUE`, the DP will assign provisional anchor IDs at the last-observed DBH census and mark those anchor rows with `ReconstructionMethod = "provisional_dp"`.

Posterior sampling:
- `POSTERIOR_SAMPLES` — default: `200L` (set to `0` to disable sampling). When `>0`, the DP engine will draw full-path posterior samples and write them into a per-run `posteriors/` subdirectory.
- `POSTERIOR_SAMPLES_FORMAT` — default: `"csv"` (options: `rds`, `feather`, `csv`).
- `POSTERIOR_SAMPLES_PATH` — default: `NULL`. If `NULL` the run's `out_dir` is used and the DP writes to `<out_dir>/posteriors/`. If you supply a path that itself ends in `posteriors`, the script strips that suffix to avoid creating nested `posteriors/posteriors` folders (the DP will create the `posteriors/` subdirectory itself).
- `POSTERIOR_SAMPLE_SEED` — default: `NULL`. If sampling is enabled and the seed is unset, the chunked runner defaults the seed to `123L` to improve reproducibility; you can override this with `--POSTERIOR_SAMPLE_SEED=<int>`.

Example posterior output files (per tag/run): `tag_11_posterior_samples__summary.csv`, `tag_11_posterior_samples__paths.csv`.

Output controls:
- `WRITE_DP_CSV` — default: `TRUE` - write incremental/combined CSV output - memory heavy.
- `WRITE_DP_RDS` — default: `TRUE` - write per-chunk RDS or combined for non-chunk runs.
- `WRITE_DP_FEATHER` — default: `FALSE` (requires the `arrow` package) - write per-chunk feather (.feather) files or combined for non-chunk runs.
- `WRITE_DP_PDF` — default: `TRUE` - write pdf visualizations per tag - memory heavy.
- `DP_PDF_INCLUDE_REFERENCE` — `main_cpp.R` default: `TRUE`; `main_cpp_chunk.R` default: `TRUE` — include biologically-informed reference lines in PDFs (useful with simulated data).
- `WRITE_DP_PDF_PER_CHUNK` — default in `main_cpp_chunk.R`: `TRUE` (controls per-chunk PDFs).

Parallel & chunking controls (chunked runner specific):
- `DP_CHUNK_SIZE` — default in `main_cpp_chunk.R`: `7L` (set `<= 0` to disable chunking behavior when applicable).
- `DP_CHUNK_RESUME` — default: `TRUE` (skip chunks whose `_done.txt` completion marker exists) — allows stopping and resuming runs. A chunk is considered complete only when its `_done.txt` file is present; partial RDS files from interrupted runs are re-processed.
- `OUT_DIR_OVERRIDE` — default: `NULL`. When set, bypasses automatic output directory creation and writes into the specified path directly. Use this to resume into an existing output directory (e.g., `--OUT_DIR_OVERRIDE=dp_global/output/<previous_run_dir>`).
- `DP_CHUNK_OVERWRITE` — default: `FALSE` (when `TRUE`, overwrite existing chunk RDS files).
- `DP_CHUNK_START`, `DP_CHUNK_END` — default: `NULL` (limit chunk range for tests).
- `RUN_ALL_TAGS` — default: `FALSE` (when `TRUE` run across all tags; chunking is recommended for large runs).
- `MANUAL_CORES` & `MANUAL_CORES_VALUE` — default: `TRUE` and `1L` respectively.

Notes on CLI differences:
- `main_cpp_chunk.R` exposes a reduced `CLI_REFERENCE` relative to `main_cpp.R` (it omits `WHICH_TAG` and leaves sensitivity/realism flags commented out) to reflect the chunked runner's intent; however it still accepts CLI overrides for many run-level parameters when invoked via Rscript.

Helpful post-run utilities:
- Merge chunk RDS/Feather files into a single CSV (run this in R or source the script and call the helper):

```r
# from R in project root
source("dp_global/scripts/main_cpp_chunk.R")
merge_chunks_to_csv("dp_global/output/<your_run_dir>")
```

This streams each chunk file to a single CSV to avoid loading the full dataset into memory.

`main_cpp.R` runs the non-chunked workflow (single-tag or parallelized tags) and does not perform per-chunk writing.

Misc:
- `USE_MEASUREMENT_ERROR` (default: `TRUE`) — enable measurement-error-aware parameter estimation.

Output directory & naming:
- `PROJECT_ROOT` (default: project root via `here::here()`) — override to set a different project root and thus change where `dp_global/output/` is created.
- `base_out_dir` (default: `dp_global/output`) — base directory where run-specific output directories are created.
- `CONFIG_NAME` (default: `NULL`) — optional string used when assembling the run-specific output directory name.
- The final `out_dir` is automatically constructed from timestamp, config name, DP mode, and other key parameters. The script writes a `run_parameters_full.txt` file into `out_dir` documenting the run configuration.

Files produced as run markers/logs:
- `run_started.txt` and `run_finished.txt` — small timestamp files written at start and finish to allow job watchers to detect progress.
- `run_parameters_full.txt` — text file capturing all important run-level variables for reproducibility.
- `run_log.txt` — appended by `log_msg()` throughout the run; writes performed via `maybe_write()` ensure directories exist and the script records success/failure messages here (for example: `Wrote RDS chunk 2: <path>`).

PDF & plotting controls:
- `WRITE_DP_PDF` (default: `TRUE`) — control whether PDFs are generated via `plot_tag_to_pdf()`.
- `DP_PDF_INCLUDE_REFERENCE` (default: `TRUE`) — include biologically-informed reference lines in PDFs.
- `PLOT_PDF_ONE_TAG_ONLY` (main: `TRUE` when `RUN_ALL_TAGS=FALSE`; not used by `main_cpp_chunk.R`) — when `TRUE` produce PDFs only for `WHICH_TAG` (useful for single-tag runs).

Sensitivity & realism flags (available in `main_cpp.R`):
- `SENSITIVITY_MODE` (default: `"none"`) — Options: `"none"`, `"run"`, `"run+write"`, `"run+write+pdf"`. Controls whether sensitivity sweeps are executed and if results are written.
- `WRITE_OUTPUTS` (derived from `SENSITIVITY_MODE`) — internal flag to control writing sensitivity outputs when requested.
- `MAKE_ALL_SWEEPS_PDF` (derived) — whether to render all sweeps to PDF when `SENSITIVITY_MODE="run+write+pdf"`.
- `RUN_REALISM_REPORT` (default: `FALSE`) — when `TRUE` the script will generate a realism report for a representative species.
- `RUN_K_SWEEP_DEMO` (default: `FALSE`) — optional demo mode for k-sweep visualizations.

Note: the chunked runner (`main_cpp_chunk.R`) disables or comments out these options because per-chunk processing does not assemble a full `out` object for full-run sensitivity/realism processing.

Biological realism settings (defaults in script):
- `MAX_GROWTH_HARD_SOURCE = "fixed"`, `MAX_GROWTH_FIXED = 7.5` (`main_cpp.R`) / `5` (`main_cpp_chunk.R`)
- `MAX_SHRINK_HARD_SOURCE = "fixed"`, `MAX_SHRINK_FIXED = -0.5`
- `K_SHRINK_SOURCE = "fixed"`, `K_SHRINK_FIXED = 0`
- `K_GROWTH_SOURCE = "fixed"`, `K_GROWTH_FIXED = 0`
- `RECRUIT_MAX_SOURCE = "fixed"`, `RECRUIT_MAX_FIXED = (MAX_GROWTH_FIXED * 5) + 0.9999`
- `USE_MEASUREMENT_ERROR = TRUE`

Notes about chunking & downstream outputs:
- When chunking is active (`DP_CHUNK_SIZE > 0`), the script processes and writes chunk outputs incrementally and intentionally sets the in-memory `out` object to `NULL` to avoid excessive memory use.
- Because `out` is not assembled in memory for chunked runs, downstream steps that expect a combined `out` (e.g., writing a single combined RDS `stem_reconstruction_dp_global_rcpp.rds`, generating per-run PDFs from a combined `out`, or creating the realism report from `out`) will be skipped. Instead, you can work with the incremental CSV or per-chunk RDS files produced by the run.

---

## Chunking details & recommendations ✅

- When `RUN_ALL_TAGS=TRUE` and `DP_CHUNK_SIZE > 0`, groups are constructed from `unique(xrun[, .(Tag, species)])` and split into chunks of size `DP_CHUNK_SIZE`.
- For each chunk `ci`:
  - The script runs the DP in parallel over the groups in the chunk (each child sets `data.table::setDTthreads(1L)` to limit thread contention).
  - The chunk's results are combined into `out_chunk` and annotated with `DP_Chunk = ci`.
  - `maybe_add_posterior_bins(out_chunk)` is applied to add `DP_PosteriorBin` if requested.
- The chunk rows include a `run_out_dir` column set to the basename of the run `out_dir`.
- If `WRITE_DP_CSV=TRUE`, `out_chunk` is appended to `stem_reconstruction_dp_global_rcpp.csv` (the script writes the header only on the first write).
- If `WRITE_DP_RDS=TRUE`, `out_chunk` is saved as `stem_reconstruction_dp_global_rcpp_chunk_###.rds` (used as a resume marker).
- If `WRITE_DP_PDF=TRUE` and `WRITE_DP_PDF_PER_CHUNK=TRUE`, the script will attempt to generate a per-chunk PDF `stem_reconstruction_dp_global_rcpp_chunk_###.pdf` using `plot_tag_to_pdf()`; PDF generation errors are logged but will not abort the run.

Memory-saving recommendations:
- Prefer the chunking + incremental CSV approach for very large datasets (keeps peak RAM low).
- Keep `WRITE_DP_CSV=TRUE` so you get a single on-disk CSV that grows incrementally (append is memory-friendly).
- Keep per-chunk RDS files (`WRITE_DP_RDS=TRUE`) as reliable completion markers for resume and for reproducibility. These RDS files contain chunk results and can be merged later using `data.table::rbindlist(lapply(chunk_files, readRDS), use.names=TRUE, fill=TRUE)` on a machine with enough RAM or processed in streaming fashion.
- Avoid assembling a full `out` in memory if your dataset is large — the script intentionally does `out <- NULL` when chunking to prevent a memory spike.

---

## Where `maybe_add_posterior_bins` is applied

- If the script runs in **non-chunked** mode (i.e., `RUN_ALL_TAGS=TRUE` but `DP_CHUNK_SIZE <= 0` or chunking disabled), the full `out` is assembled in memory and `maybe_add_posterior_bins(out)` is applied immediately after assembly.
- If **chunking** is active, `maybe_add_posterior_bins()` is applied to each `out_chunk` before the chunk is written to disk.

This ensures posterior-bin computation is done while memory per-chunk is small and avoids double-processing.

Runner integration

- `main_cpp.R` and `main_cpp_chunk.R` are the primary entrypoints and are designed to be invoked directly with `Rscript`. External orchestrators can build CLI flags using the canonical names in `CLI_REFERENCE` (defined in each script); keys are case-insensitive and can use `-` or `_`.

## Post-reconstruction notes

At the anchor census (`CensusID == ANCHOR_START_CENSUS`), `TrueStemID` values serve as hard constraints and `ReconstructedStemID` will equal `TrueStemID` for those rows (`ReconstructionMethod = "given"`). Pre-anchor rows with `TrueStemID` are processed by the solver (DP or probabilistic) and can receive different identity assignments based on biological likelihood.

### Hard-invariant sweep and the `SweepAuditOverride` column

When `PIN_TRUESTEMID = TRUE` (default), an idempotent **hard-invariant sweep** runs at three sites — inside `finalize_out()` (DP path), inside `match_stems_probabilistic()` (probabilistic fallback), and at the script-level inside `run_dp_one_group()` — forcing every row with a non-NA `TrueStemID` to `ReconstructedStemID = TrueStemID` and `ReconstructionMethod = "given"`. This guarantees the invariant on rows the engines never visit (NA-DBH terminal rows anchored by the BCI driver's Steps 2/3, MF re-insertion edge cases, and probabilistic-fallback leaks).

In the rare case where the engine had already assigned a non-NA `ReconstructedStemID` that disagrees with `TrueStemID`, the sweep silently overrides it. Each such row is flagged `TRUE` in the **`SweepAuditOverride`** boolean column (FALSE elsewhere) and a `[audit]` log line is emitted naming the tag and override count. Downstream uncertainty consumers should treat any `SweepAuditOverride == TRUE` row as observed (P=1, entropy=0); the `DP_PosteriorTop*` columns on those rows describe the engine's overridden choice, not the final `ReconstructedStemID`.

The engine's pre-sweep `ReconstructedStemID` is preserved in the **`ReconstructedStemID_PreSweep`** column (snapshot taken once by the first sweep layer that fires). On `SweepAuditOverride == FALSE` rows it equals `ReconstructedStemID`; on `SweepAuditOverride == TRUE` rows it carries the engine's original (overridden) ID, allowing direct comparison without re-running the engine.

### Post-engine `carried_terminal` backfill (all drivers)

After the DP/probabilistic engine returns and `maybe_add_posterior_bins()` has run, **every driver** (`main_cpp.R`, `main_cpp_chunk.R`, `main_cpp_bci.R`) calls the shared helper `apply_carried_terminal_backfill()` defined in `dp_global/R/dp_global_main.R`.

The helper finds rows where **all three** conditions hold:

- `is.na(ReconstructedStemID)` (engine produced no assignment),
- `is.na(DBH)` (no measurement to match against), and
- `Status %in% c("dead", "stem dead", "broken below")` (a terminal event for the stem).

For each such row, after sorting by `(Tag, OriginalStemID, CensusID)`, it copies the most recent prior non-NA `ReconstructedStemID` from the same `(Tag, OriginalStemID)` group (LOCF) and sets `ReconstructionMethod = "carried_terminal"`. Biologically, a death/break row ends the trajectory of the most recent prior identity that shared the `OriginalStemID` — without this fill those rows would be dropped from any downstream trajectory.

Where it fires:

- `main_cpp.R` — Step 5.5b, immediately after `maybe_add_posterior_bins(out)`.
- `main_cpp_chunk.R` — inside the per-chunk parallel block, applied to each `out_chunk` with `verbose = FALSE` to keep multi-tag chunked logs quiet.
- `main_cpp_bci.R` — Step 9b. The BCI driver also performs an upstream pre-DP `TrueStemID` propagation (Steps 1–3, see the BCI section below) that is *not* applicable to non-BCI inputs; the post-engine `carried_terminal` step is the only piece that is universal across all three drivers.

**Warning messages from the probabilistic matcher** (sample-level repair counts, ME cumulative-shrinkage breaks, growth-aware resolver diagnostics) are emitted via `message()` on stderr and are also printed to stdout via `cat()` when `DP_VERBOSE=TRUE`. To capture all warnings in a log file, redirect both streams: `Rscript ... > log.txt 2>&1`.


---

## BCI Debug Driver (`main_cpp_bci.R`)

`main_cpp_bci.R` is a debug driver for single-tag runs on BCI (Barro Colorado Island) multi-stem census data. It sources `main_cpp.R` to inherit all helper functions and CLI handling, then overrides the input file and a few BCI-specific defaults, and adds a **pre-DP `TrueStemID` reconstruction** step that is specific to BCI's identity conventions.

### Pre-DP TrueStemID propagation (BCI-specific, Steps 1–3)

Before calling the DP, the BCI driver writes `TrueStemID` on every row whose biological identity is unambiguous from the BCI database conventions:

- **Step 1a** — any row with a non-NA `StemTag` (the field crew physically tagged the stem) gets `TrueStemID = OriginalStemID`.
- **Step 1b** — any row at `CensusID >= 7` (BCI's systematic re-tagging campaign from 2010 onward) gets `TrueStemID = OriginalStemID`.
- **Step 2a/b/c** — within each `(Tag, OriginalStemID)` group, after the last live DBH measurement, terminal-event rows (`Status %in% c("dead", "stem dead", "broken below")` or R-family resprout codes in `ListOfTSM`) are anchored to their own `OriginalStemID`; remaining post-last-DBH gaps are filled by bidirectional LOCF/NOCB.
- **Step 3a** — same direct anchoring of any remaining terminal-event rows, dropping the post-last-DBH guard.
- **Step 3b** — within each `(Tag, OriginalStemID)` group, if all non-NA `TrueStemID` values agree, fill remaining NAs with that value; conflicting groups are left alone and counted in a diagnostic message.

Pre-anchor rows without an unambiguous identity remain NA and are resolved by the DP. The combination of this pre-DP propagation and the DP/probabilistic engines' hard-invariant sweep (see below) guarantees `ReconstructedStemID == TrueStemID` on every row that Steps 1–3 anchored.

### Post-DP carried_terminal backfill (Step 9b)

After `run_dp_one_group()` returns and `maybe_add_posterior_bins()` has been applied, the BCI driver invokes the shared helper `apply_carried_terminal_backfill()` (Step 9b). This is the same helper called by `main_cpp.R` (Step 5.5b) and `main_cpp_chunk.R`; see the *Post-engine `carried_terminal` backfill* section above for the rule and rationale. The helper is the only post-engine piece that is identical across all three drivers — the upstream Steps 1–3 above are BCI-input-specific and have no analogue in the simulator (which already ships `TrueStemID`).

### BCI-specific defaults

- `INPUT_FILE`: `bci_data/multistem_tags.rds` (loaded via `readRDS()`; not tracked by git)
- `WHICH_TAG`: `"115203"` (a multi-stem debugging tag)
- `FORCE_ONE_SPECIES_PARAMETERS`: `FALSE` (use real BCI species)
- `MAX_GROWTH_FIXED`: `5.0`, `MAX_SHRINK_FIXED`: `-0.5`, `DP_MAX_STATES`: `1039L`
- `ANCHOR_START_CENSUS`: `7L`
- `PIN_TRUESTEMID`: `TRUE`
- Output written to `dp_global/output/<timestamp>_BCI_tag<WHICH_TAG>_*` under the BCI-specific run directory naming scheme

### Example invocations

```bash
# Default tag
Rscript dp_global/scripts/main_cpp_bci.R

# Specific tag
Rscript dp_global/scripts/main_cpp_bci.R --WHICH_TAG=123375

# Quiet mode
Rscript dp_global/scripts/main_cpp_bci.R --WHICH_TAG=187064 --DP_VERBOSE=FALSE

# Tight state-space cap (force probabilistic fallback for testing)
Rscript dp_global/scripts/main_cpp_bci.R --WHICH_TAG=000184 --DP_MAX_STATES=2
```

### Prerequisites

- R packages: `data.table`, `here` (in addition to the standard DP prerequisites)
- `bci_data/multistem_tags.rds` (or whichever `INPUT_FILE` is set to) must exist; not tracked by git

---

## Example invocations

- Run a single tag interactively:

```bash
Rscript dp_global/scripts/main_cpp.R --INPUT_FILE=data_simulation/data/simulated_data_1.csv --WHICH_TAG=20
```

- Run all tags with chunking (use the chunked driver):

Note: `dp_global/scripts/main_cpp_chunk.R` is the chunked driver intended for large runs. You can edit configuration variables at the top of the file or pass overrides via CLI flags (e.g., `--DP_CHUNK_SIZE=7`) when invoking it with `Rscript`. Run it directly via:

```
Rscript dp_global/scripts/main_cpp_chunk.R
```

- Run all tags without chunking (not recommended for very large datasets):

```
Rscript dp_global/scripts/main_cpp.R --RUN_ALL_TAGS=TRUE --DP_CHUNK_SIZE=0
```

---

## Post-run: combining chunk files or troubleshooting

- By default, the merge helpers write to `stem_reconstruction_dp_global_rcpp_merged.csv` in the run `out_dir`. They locate chunk files matching `stem_reconstruction_dp_global_rcpp_chunk_\\d{3}\\.rds` or `stem_reconstruction_dp_global_rcpp_chunk_\\d{3}\\.feather` (you can use `prefer = 'feather'` with `merge_chunks_to_csv()` to prefer Feather sources).

To merge per-chunk RDS files into a single RDS (only do this on a machine with enough RAM):

```r
library(data.table)
chunk_files <- list.files("<out_dir>", pattern = "stem_reconstruction_dp_global_rcpp_chunk_\\d{3}\\.rds$", full.names = TRUE)
all <- rbindlist(lapply(chunk_files, readRDS), use.names = TRUE, fill = TRUE)
```

- If `stem_reconstruction_dp_global_rcpp.csv` already exists when you resume with `DP_CHUNK_RESUME=TRUE`, the script will append new chunk rows and skip chunks whose `_done.txt` completion marker exists. Use `--OUT_DIR_OVERRIDE=<path>` to point at the existing output directory.

---

## Basal Area Uncertainty (`basal_area_uncertainty.R`)

Post-processing script that uses posterior path samples from a completed run to quantify how identity uncertainty propagates into individual-level basal area (BA) estimates.

### Usage

```bash
Rscript dp_global/scripts/basal_area_uncertainty.R \
  --RUN_DIR=dp_global/output/<run_dir>
```

The script reads:
- `<RUN_DIR>/stem_reconstruction_dp_global_rcpp.csv` (main reconstruction)
- `<RUN_DIR>/posteriors/tag_*_posterior_samples__paths.csv` (posterior paths)

### Outputs

Two CSV files and one PDF written to `<RUN_DIR>/`:

| File | Rows | Description |
|------|------|-------------|
| `basal_area_tag_census.csv` | Tag × Census | Total BA (m²), stem count, year |
| `basal_area_tag_change.csv` | Tag × Census interval | BA change decomposed into survivor growth, mortality loss, and recruitment gain — MAP values plus posterior mean, SD, and 95% CI for each component (all in m²) |
| `basal_area_figures.pdf` | — | Multi-page PDF: summary page (all tags), per-tag detail pages (2×2 panels: BA trajectory, stem count, decomposition bars with uncertainty whiskers, stem demographics), posterior uncertainty histograms, and **posterior density plots** (kernel densities of Growth/Loss/Gain/DeltaBA with weighted-mean vertical lines — one page pooled across all census intervals, then one page per interval) |

### Key insight

Tag-level total BA per census is **invariant** to identity assignment — the same DBH values are summed regardless of which stem identity each observation receives. The **decomposition** of BA change into growth (surviving stems), loss (mortality), and gain (recruitment) **is** identity-dependent. Different posterior path samples assign different stems as survivors vs. deaths vs. recruits, producing uncertainty in the attribution of BA change to these demographic components.

---