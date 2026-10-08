# dp_global/scripts

This document describes the current behavior of the `dp_global` driver scripts:

- `dp_global/scripts/main_cpp.R` — the interactive/CLI driver for single-tag or targeted runs.
- `dp_global/scripts/main_cpp_chunk.R` — the chunked driver optimized for large runs.
- `dp_global/scripts/main_cpp_bci.R` — the BCI debug driver for single-tag runs on BCI census data.
- `dp_global/scripts/basal_area_uncertainty.R` — posterior-based basal area uncertainty quantification.

Both `main_cpp.R` and `main_cpp_chunk.R` accept command-line overrides of defaults using `--KEY=VALUE` flags. Keys are case-insensitive and may use `-` or `_` as separators (write the keys whose values must stay text, such as `WHICH_TAG`, with `_`). Any top-level variable of the script can be overridden, not only the keys listed by `--help`; an unknown key gives a warning and is ignored. `main_cpp_bci.R` inherits the same CLI interface from `main_cpp.R`.

Run every script from the **project root**: `dp_global/R/dp_global_main.R` resolves the modules and the C++ file from the working directory.

R packages: `data.table`, `here`, `Rcpp`, `igraph`, `MASS` (loaded by `dp_global_main.R`), `lpSolve` (probabilistic matcher), `parallel`; optional `arrow` (Feather output), `ggplot2` and `cowplot` (PDF plots), `reshape2` and `scales` (bio-parameter report, k-sweep plot); `ggplot2`, `patchwork` and `scales` for `basal_area_uncertainty.R`.

---

## Overview ✨

`main_cpp.R` is the central driver for the `dp_global` pipeline. It:

- Loads input tree census data (CSV via `INPUT_FILE`).
- Ensures species information is present and optionally forces a single species label.
- Estimates biological parameters per species (via `estimate_bio_pars`).
- Runs Dynamic Programming (DP) reconstruction logic (via `run_dp_one_group` which calls the Rcpp-based DP functions).
- Optionally runs sensitivity sweeps and realism checks.
- Writes outputs (CSV, RDS, PDF) to an automatically created output directory, plus `bio_pars_report.pdf` (plots of the estimated parameters).

You can run either script directly via Rscript:

- Single tag: `Rscript dp_global/scripts/main_cpp.R --INPUT_FILE=... --WHICH_TAG=...`
- All tags of a small dataset: `Rscript dp_global/scripts/main_cpp.R --INPUT_FILE=... --RUN_ALL_TAGS=TRUE --PLOT_PDF_ONE_TAG_ONLY=FALSE`
- Large dataset (chunked, recommended): `Rscript dp_global/scripts/main_cpp_chunk.R --INPUT_FILE=...`

## Stem Identity Renumbering

All drivers (`main_cpp.R`, `main_cpp_chunk.R`, `main_cpp_bci.R`) use a universal post-engine helper chain: `maybe_add_posterior_bins()`, `apply_pin_track_rejoin()`, `apply_carried_terminal_backfill()`, `apply_orphan_stem_backfill()`, `apply_terminal_to_host()`, `apply_broken_below_invariants()`, `renumber_engine_minted_ids()`, and finally `finalize_posterior_paths()`. After these steps, all `ReconstructedStemID` values are renumbered **sequentially from 1 to N within each tag**, ordered by the earliest census in which each stem appears. If multiple stems first appear in the same census, the largest DBH at that census gets the lower ID, with ties broken by original ID. **Negative or zero IDs are never produced.**

---

## Important behavior notes ⚠️

- `main_cpp.R` defaults to single-tag mode (`RUN_ALL_TAGS = FALSE`, `WHICH_TAG = "20"`). Use `--RUN_ALL_TAGS=TRUE` to run across all tags. `PLOT_PDF_ONE_TAG_ONLY` is set from the default `RUN_ALL_TAGS` before the command line is applied, so add `--PLOT_PDF_ONE_TAG_ONLY=FALSE` to get every tag in the PDF.
- `main_cpp_chunk.R` always processes every (Tag, species) group; there `RUN_ALL_TAGS` and `WHICH_TAG` only enter the name of the output folder (`allT` or `T0`).
- For large datasets, prefer `main_cpp_chunk.R` which divides groups (Tag + species) into chunks (`DP_CHUNK_SIZE`) and writes chunk outputs incrementally to disk to keep peak RAM low.
- When chunking, the combined `out` is not assembled in memory. Instead:
  - Each chunk is appended to `stem_reconstruction_dp_global_rcpp.csv` (if `WRITE_DP_CSV=TRUE`),
  - Each chunk is saved as `stem_reconstruction_dp_global_rcpp_chunk_###.rds` (if `WRITE_DP_RDS=TRUE`) and serves as a resume/completion marker.
- `maybe_add_posterior_bins()` is applied per-chunk when chunking; in non-chunked runs it is applied after assembling the full `out`.
- The chunked driver disables or comments out features that are not applicable to chunked runs (for example: full-run sensitivity sweeps and realism report generation).

---

## Key CLI flags & defaults 🔧

Common flags used by both drivers (case-insensitive, but use capital letters to keep it clean):

- `INPUT_FILE` — default: `data_simulation/data/simulated_data_1.csv`. A CSV (read with `data.table::fread()`) with at least `Tag`, `CensusID`, `DBH` (cm), `ExactDate` and `TrueStemID`; `OriginalStemID` (or `StemID`), `Status`, `ListOfTSM`, `growth_form` and `HOM` are used when present.
- `FORCE_ONE_SPECIES_PARAMETERS` — default: `TRUE`.
- `FORCED_SPECIES_LABEL` — default: `"all"` (used when forcing one species label).
- `SPECIES_COL` — default: `NULL` (auto-detected if not set).
- `USE_MEASUREMENT_ERROR` — default: `TRUE`. Applies to the parameter estimation and to the DP growth likelihood (the probabilistic matcher has no measurement-error mixture).

DP / reconstruction option

- `DP_MODE` — default: `"marginals+bins"`. Allowed: 'none' 'marginals' 'marginals+bins' 'map'. `'map'` and `'marginals'` run the same DP (MAP assignment plus marginals) and differ only in the folder label; `'marginals+bins'` adds the `DP_PosteriorBin` column; `'none'` skips the DP in `main_cpp.R` only (the chunked runner runs it for every mode).
- `WHICH_TAG` — character; used for single-tag runs (relevant to `main_cpp.R`). Must match the `Tag` column exactly (e.g., `"084555"` preserves leading zeros). The chunked runner processes groups (`Tag`, `species`) and does not rely on `WHICH_TAG`.
- `ANCHOR_START_CENSUS` — default: `7L`.
- `ALLOW_PROVISIONAL_DP_ANCHOR` — default: `TRUE` — when `TRUE` the DP can assign provisional anchor IDs at the last observed DBH census if the requested anchor census lacks `TrueStemID` but has DBH; set to `FALSE` to require an explicit anchored census or to fall back to the probabilistic matcher.
- `DP_VERBOSE` — default: `TRUE`.
- `DP_POSTERIOR_TOP_K` — default: `2L` - number of most probable identities reported per observation (columns `DP_PosteriorTop<k>ID` / `DP_PosteriorTop<k>Prob`).
- `DP_MAX_TRACKS` — default: `NULL` (auto-computed per data when `NULL`) — optionally force max tracks; if `NULL` auto-computed.
- `DP_MAX_STATES` — default: `40000L`. Maximum number of injective assignment states allowed per census. Also controls the inter-census transition limit (`max_edges = max_states²`). The number of states at a census with $n$ observed stems and $K$ tracks is $P(K,n) = K!/(K-n)!$. When any census exceeds `max_states`, or when the cross-product of states between two adjacent censuses reaches `max_states²`, the solver falls back to the probabilistic matcher. (When pins or growth bounds restrict the tracks of a census, `enumerate_states_constrained()` keeps the first `max_states` assignments instead of signalling the overflow, so such a census is solved on a truncated state set.) **Practical limits with default 40,000** (for $K = n + 1$ tracks): DP handles up to 6 observed stems per census exactly (P(7,6) = 5,040 < 40,000); 7+ stems trigger fallback (P(8,7) = 40,320 > 40,000). With `DP_MAX_STATES = 1,000`: max 5 stems. With `DP_MAX_STATES = 20,000`: max 6 stems. With more tracks the limits are lower. See `dp_global/README.md` for detailed tables and how to choose a value.
- `DP_SLACK_TRACKS` — default: `1L` - slack (additional) tracks allowed for DP.
- `DP_SLACK_REQUIRE_ANCHOR_RECRUITABLE` — default: `TRUE` - require anchor to be recruitable (if DBH less than max recruitment size) before granting slack.
- `DP_SLACK_REQUIRE_ANCHOR_EPS` — default: `1e-6`.
- `PROB_N_SAMPLES` — default: `200L`. Number of Gumbel-noise stochastic samples drawn by the probabilistic matcher. Higher values produce more accurate posterior estimates at the cost of computation time.
- `PROB_SPECIES` — default: `character(0)`. Values of the `Species` column sent straight to the probabilistic matcher (comma- or semicolon-separated on the command line).
- `PROB_N_SIGMA_ME` — `main_cpp.R` default: `3` (not a variable of `main_cpp_chunk.R`, where the engine default of 3 applies). Threshold of the cumulative-shrinkage repair of the probabilistic matcher, in SD of the measurement error; `Inf` turns that layer off.
- `PROB_LOOKAHEAD_WEIGHT` — default: `1`. Weight for sequential backward conditioning in the probabilistic matcher. When > 0 and $K \geq 4$, the cost matrix for each census pair is conditioned on the already-resolved forward assignment, improving path continuity. Set to `0` to disable conditioning.
- `USE_BIO_HARD_SHRINK_IN_PROB` — default: `TRUE`. Controls two hard constraints in the probabilistic matcher: (1) the bio shrink gate in pairwise edge construction (`g < Bio_Max_Shrink` → edge log-likelihood set to `-Inf`), and (2) the Layer 2 ME cumulative-shrinkage repair check in `repair_stitched_growth_violations()`. When `FALSE`, both are disabled — edges with growth below `Bio_Max_Shrink` are allowed (penalised by the soft `k_shrink` quadratic term only) and the ME cumulative check does not sever trajectories. Use `FALSE` to allow confirmed large-shrinkage events (e.g., storm damage) without forcing the matcher to split a continuous stem. **Note:** The exact DP solver is unaffected — its transition cost always charges the hard penalty for growth below `Bio_Max_Shrink`, regardless of this flag.
- `USE_BIO_HARD_GROWTH_IN_PROB` — default: `TRUE`. Controls the bio growth gate in pairwise edge construction in the probabilistic matcher (`g > Bio_Max_Growth` → edge log-likelihood set to `-Inf`). When `FALSE`, edges exceeding the biological growth maximum are allowed (penalised by the soft `k_growth` quadratic term only). The exact DP solver is unaffected by this flag.
- `PIN_TRUESTEMID` — default: `TRUE`. When `TRUE`, observations with a known `TrueStemID` at non-anchor censuses are pinned to their field-observed identity track in the DP (state enumeration) and in the probabilistic matcher, reducing effective state space and preventing re-identification of labelled stems. Two stems with different pins are never joined: when the pins forbid links that a census pair's death/recruit slots relied on, the pair gets extra slots (`pin_masked_pair()`, see `dp_global/README.md`).
- Birth-death assignment in the probabilistic matcher: these drivers use the
  engine default (`prob_birth_death = TRUE` in
  `match_stems_dp_global_backward_marginals_batch()`, `birth_death = TRUE` in
  `match_stems_probabilistic()`): in every census pair each stem may continue,
  die or be recruited, scored with the same biological parameters as the DP. The BCI
  driver exposes it as `PROB_BIRTH_DEATH` (`FALSE` = count-based slots: with a
  stable stem count every stem must continue; see
  `dp_global/README.md`, *Probabilistic Matching Fallback*).
- DBH recorded in classes (`dbh_round_censuses` of the engine) is not set by
  these drivers; the BCI driver exposes it as `DBH_ROUND_CENSUSES`.
- `DP_FALLBACK_GROWTH_FORMS` — default: `character(0)`; comma- or
  semicolon-separated list of values in the `growth_form` column that should
  trigger an immediate probabilistic fallback and prevent the DP solver from running
  on that tag. The driver automatically splits the string into a vector. Values
  are matched exactly against the data's labels (the BCI driver checks them at
  load and stops on a value that is not a label).
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
- `NON_TAPER_CORRECTED_PRUNE_MAX_GROWTH` — default: `1.25 × MAX_GROWTH_FIXED`
  (`9.375` cm/year in `main_cpp.R`, `6.25` in `main_cpp_chunk.R`); upper annual growth bound for non-taper-corrected
  growth forms during pruning. Same interplay with the general bounds as above.
- `HOM_TOLERANCE_SCALE` — default: `2.0` (cm of DBH per meter of HOM deviation);
  when a non-taper-corrected tag has a `hom` (or `HOM`) column in its data,
  the prune bounds are widened for each census pair by
  `hom_tolerance_scale × max(|HOM − 1.3|) / interval_years`.
  This is where the effective differentiation between regular trees and
  non-taper-corrected forms occurs in the default configuration. NA HOM values
  are treated as 1.3 (zero contribution). Set to `0` to disable
  HOM widening while keeping the wide base bounds.

**Anchor scoping and post-anchor preservation:** If observations exist after the requested `ANCHOR_START_CENSUS`, the DP is scoped to censuses <= `ANCHOR_START_CENSUS` and post-anchor rows are preserved and appended to the output. Post-anchor rows with a non-NA `TrueStemID` (with or without a `DBH`) are set to `ReconstructedStemID = TrueStemID` and `ReconstructionMethod = "given"`. Remaining post-anchor rows receive `ReconstructionMethod = "none_after_anchor"`. If scoping removes all pre-anchor observations, the original rows are returned and anchor-census `TrueStemID` values are treated as `"given"` while other rows are labeled `"none_after_anchor"`. **Note:** Pre-anchor rows (censuses before the anchor) with non-NA `TrueStemID` are processed by the solver and can receive different identity assignments based on biological likelihood — however, when `PIN_TRUESTEMID=TRUE` (default), they are constrained to their field-observed identity.

**Provisional anchor behavior:** When a requested anchor census lacks `TrueStemID` but contains DBH observations and `ALLOW_PROVISIONAL_DP_ANCHOR=TRUE`, the DP will assign provisional anchor IDs at the last-observed DBH census and mark those anchor rows with `ReconstructionMethod = "provisional_dp"`.

Posterior sampling:

- `POSTERIOR_SAMPLES` — default: `200L` (set to `0` to disable sampling). When `>0`, each engine (DP and probabilistic) draws full-path posterior samples and writes a per-run `posteriors/` subdirectory under `out_dir`.
- `POSTERIOR_SAMPLES_FORMAT` — default: `"csv"` (options: `rds`, `feather`, `csv`). `feather` requires the `arrow` package; if `arrow` is missing the writer silently falls back to `rds`. `basal_area_uncertainty.R` reads CSV path files only.
- `POSTERIOR_SAMPLES_PATH` — default: `NULL`. If `NULL` the run's `out_dir` is used and posteriors are written to `<out_dir>/posteriors/`. If you supply a path that itself ends in `posteriors`, the script strips that suffix to avoid creating nested `posteriors/posteriors` folders (the engine creates the `posteriors/` subdirectory itself).
- `POSTERIOR_SAMPLE_SEED` — default: `NULL`. If sampling is enabled and the seed is unset, both drivers set the seed to `123L`; you can override this with `--POSTERIOR_SAMPLE_SEED=<int>`. The same seed is used by the probabilistic matcher.

**Posterior writing pipeline (staging architecture).** The engines do not write the final paths file themselves. They stage the raw per-sample reconstruction table (in engine ID space, with a `logp` weight column from the DP) to:

```
<out_dir>/posteriors/.staging/tag_<Tag>_samples_raw_<BATCH_TS>.rds
```

Each driver script then runs the standard post-engine helper chain — `maybe_add_posterior_bins()` → `apply_pin_track_rejoin()` → `apply_carried_terminal_backfill()` → `apply_orphan_stem_backfill()` → `apply_terminal_to_host()` → `apply_broken_below_invariants()` — and finally calls, in order:

1. `renumber_engine_minted_ids(out_chunk, posterior_top_k = DP_POSTERIOR_TOP_K, ...)` — renumbers `ReconstructedStemID` from 1 to N within each tag, in order of first appearance, and returns a `Tag/old_id/new_id` mapping table. It writes no file; its `posterior_samples_path` and `mapping_format` arguments are not used.
2. `finalize_posterior_paths(out_chunk, posterior_samples_path = out_dir, mapping = .renum$mapping, ...)` — for each staging file, translates `ReconstructedStemID` via the mapping, applies the `TrueStemID` pins per sample (`apply_pins_to_samples()`), runs `apply_bb_invariants_to_samples()` in the renumbered ID space, computes `path_sig` / `path_count` / `path_prob` / `recon`, writes the final `tag_<Tag>_posterior_samples_<BATCH_TS>_paths.<feather|rds|csv>` file, and deletes the staging file on success.

A healthy completed run therefore leaves `<out_dir>/posteriors/.staging/` empty and `<out_dir>/posteriors/` populated with one `*_paths.<ext>` file per tag for which posteriors were drawn. The `*_paths.<ext>` file is the only posterior artefact written per tag. `path_count` is the number of samples with a path; `path_prob` is that count divided by the number of samples when the samples carry no `logp` (probabilistic matcher), and the share of `exp(logp)` of the path's samples when they do (DP).

Example posterior output file (per tag/run, BCI defaults with `POSTERIOR_SAMPLES_FORMAT="feather"`): `posteriors/tag_11_posterior_samples_20260528_215527_paths.feather`. When `BATCH_TS` is empty (the default for the dp_global scripts), the timestamp segment collapses, e.g. `posteriors/tag_11_posterior_samples__paths.csv`.

Output controls:

- `WRITE_DP_CSV` — default: `TRUE` - write incremental/combined CSV output - memory heavy.
- `WRITE_DP_RDS` — default: `TRUE` - write per-chunk RDS or combined for non-chunk runs.
- `WRITE_DP_FEATHER` — default: `FALSE` (requires the `arrow` package) - write per-chunk feather (.feather) files or combined for non-chunk runs.
- `WRITE_DP_PDF` — default: `TRUE` - write pdf visualizations per tag - memory heavy.
- `DP_PDF_INCLUDE_REFERENCE` — `main_cpp.R` default: `TRUE`; `main_cpp_chunk.R` default: `TRUE` — add a second panel with the trajectories grouped by `OriginalStemID` to each PDF page (needs that column; useful with simulated data).
- `WRITE_DP_PDF_PER_CHUNK` — default in `main_cpp_chunk.R`: `TRUE` (controls per-chunk PDFs).

Parallel & chunking controls (chunked runner specific):

- `DP_CHUNK_SIZE` — default in `main_cpp_chunk.R`: `7L`; number of (Tag, species) groups per chunk (`0` puts all groups in one chunk).
- `DP_CHUNK_RESUME` — default: `TRUE` (skip chunks whose `_done.txt` completion marker exists) — allows stopping and resuming runs. A chunk is considered complete only when its `_done.txt` file is present; partial RDS files from interrupted runs are re-processed.
- `OUT_DIR_OVERRIDE` — default: `NULL`. When set, bypasses automatic output directory creation and writes into the specified path directly. Use this to resume into an existing output directory (e.g., `--OUT_DIR_OVERRIDE=dp_global/output/<previous_run_dir>`).
- `DP_CHUNK_OVERWRITE` — default: `FALSE` (when `TRUE`, chunks that already have a `_done.txt` marker are run again and their files overwritten).
- `DP_CHUNK_START`, `DP_CHUNK_END` — default: `NULL` (limit chunk range for tests).
- `RUN_ALL_TAGS` — default: `FALSE` (`main_cpp.R`: when `TRUE` run across all tags, in parallel over `MC_CORES`; the chunked runner processes all groups whatever its value).
- `MANUAL_CORES` & `MANUAL_CORES_VALUE` — default: `TRUE` and `1L` respectively. With `MANUAL_CORES=FALSE` the number of workers is `parallel::detectCores() - 1`. Workers are forked with `parallel::mclapply()`, which supports more than one worker only on Unix-alikes.

Notes on CLI differences:

- `main_cpp_chunk.R` exposes a reduced `CLI_REFERENCE` relative to `main_cpp.R` (it omits `WHICH_TAG`, `PROB_N_SIGMA_ME` and the sensitivity/realism flags, and adds the `DP_CHUNK_*` and `WRITE_DP_PDF_PER_CHUNK` keys); like `main_cpp.R` it accepts an override for any variable it defines.

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
- `DP_PDF_INCLUDE_REFERENCE` (default: `TRUE`) — add the `OriginalStemID` panel to each PDF page.
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

- The chunked runner processes and writes chunk outputs incrementally and removes each chunk from memory before the next one.
- Because no combined `out` is assembled in memory for chunked runs, it writes no combined RDS (`stem_reconstruction_dp_global_rcpp.rds`), no run-level PDF and no realism report. Instead, you can work with the incremental CSV or per-chunk RDS files produced by the run.

---

## Chunking details & recommendations ✅

- In `main_cpp_chunk.R`, groups are constructed from `unique(xrun[, .(Tag, species)])` and split into chunks of size `DP_CHUNK_SIZE`.
- For each chunk `ci`:
  - The script runs the DP in parallel over the groups in the chunk (each child sets `data.table::setDTthreads(1L)` to limit thread contention).
  - The chunk's results are combined into `out_chunk` and annotated with `DP_Chunk = ci`.
  - `maybe_add_posterior_bins(out_chunk)` is applied to add `DP_PosteriorBin` if requested, followed by the rest of the post-engine helper chain and `finalize_posterior_paths()`.
  - A chunk that raises an error leaves `stem_reconstruction_dp_global_rcpp_chunk_###_failed.txt` and the run continues with the next chunk.
- The chunk rows include a `run_out_dir` column set to the basename of the run `out_dir`.
- If `WRITE_DP_CSV=TRUE`, `out_chunk` is appended to `stem_reconstruction_dp_global_rcpp.csv` (the script writes the header only on the first write).
- If `WRITE_DP_RDS=TRUE`, `out_chunk` is saved as `stem_reconstruction_dp_global_rcpp_chunk_###.rds`. The resume marker is the `..._chunk_###_done.txt` file written at the end of each chunk.
- If `WRITE_DP_PDF=TRUE` and `WRITE_DP_PDF_PER_CHUNK=TRUE`, the script will attempt to generate a per-chunk PDF `stem_reconstruction_dp_global_rcpp_chunk_###.pdf` using `plot_tag_to_pdf()`; PDF generation errors are logged but will not abort the run.

Memory-saving recommendations:

- Prefer the chunking + incremental CSV approach for very large datasets (keeps peak RAM low).
- Keep `WRITE_DP_CSV=TRUE` so you get a single on-disk CSV that grows incrementally (append is memory-friendly).
- Keep per-chunk RDS files (`WRITE_DP_RDS=TRUE`) for reproducibility and for rebuilding the CSV. These RDS files contain chunk results and can be merged later using `data.table::rbindlist(lapply(chunk_files, readRDS), use.names=TRUE, fill=TRUE)` on a machine with enough RAM or processed in streaming fashion.
- Avoid assembling a full `out` in memory if your dataset is large — the chunked runner never holds more than one chunk.

---

## Where `maybe_add_posterior_bins` is applied

- `main_cpp.R` (non-chunked) assembles the full `out` in memory and applies `maybe_add_posterior_bins(out)` immediately after assembly.
- `main_cpp_chunk.R` applies `maybe_add_posterior_bins()` to each `out_chunk` before the chunk is written to disk.

This ensures posterior-bin computation is done while memory per-chunk is small and avoids double-processing.

Runner integration

- `main_cpp.R` and `main_cpp_chunk.R` are the primary entrypoints and are designed to be invoked directly with `Rscript`. External orchestrators can build CLI flags using the canonical names in `CLI_REFERENCE` (defined in each script); keys are case-insensitive and can use `-` or `_`.

## Post-reconstruction notes

At the anchor census (`CensusID == ANCHOR_START_CENSUS`), `TrueStemID` values serve as hard constraints and `ReconstructedStemID` will equal `TrueStemID` for those rows (`ReconstructionMethod = "given"`). Pre-anchor rows with `TrueStemID` are pinned to that identity when `PIN_TRUESTEMID = TRUE` (default); with `PIN_TRUESTEMID = FALSE` the solver (DP or probabilistic) is free to give them another identity.

### Hard-invariant sweep and the `SweepAuditOverride` column

When `PIN_TRUESTEMID = TRUE` (default), an idempotent **hard-invariant sweep** runs at three sites — inside `finalize_out()` (fallback paths and resprout segment splits of the DP function), inside `match_stems_probabilistic()` (probabilistic matcher), and at the script-level inside `run_dp_one_group()` (every output; the only sweep for tags solved by the plain DP path, which enforces the pins through its state enumeration) — forcing every row with a non-NA `TrueStemID` to `ReconstructedStemID = TrueStemID` and `ReconstructionMethod = "given"`. This guarantees the invariant on rows the engines never visit (NA-DBH terminal rows anchored by the BCI driver's Steps 2/3, MF re-insertion edge cases, and probabilistic-fallback leaks).

In the rare case where the engine had already assigned a non-NA `ReconstructedStemID` that disagrees with `TrueStemID`, the sweep silently overrides it. Each such row is flagged `TRUE` in the **`SweepAuditOverride`** boolean column (FALSE elsewhere) and a `[audit]` log line is emitted naming the tag and override count. Downstream uncertainty consumers should treat any `SweepAuditOverride == TRUE` row as observed (P=1, entropy=0); the `DP_PosteriorTop*` columns on those rows describe the engine's overridden choice, not the final `ReconstructedStemID`.

The engine's pre-sweep `ReconstructedStemID` is preserved in the **`ReconstructedStemID_PreSweep`** column (snapshot taken once by the first sweep layer that fires). On `SweepAuditOverride == FALSE` rows it equals `ReconstructedStemID`; on `SweepAuditOverride == TRUE` rows it carries the engine's original (overridden) ID, allowing direct comparison without re-running the engine.

#### Duplicate-aware pinning and `SweepRollbackToPreSweep`

The script-level sweep additionally enforces `ReconstructedStemID` uniqueness within each `(Tag, CensusID)`. Pins are processed in stable row order; before applying `ReconstructedStemID := TrueStemID` on a row, the sweep checks whether that integer value is already present on another row at the same `(Tag, CensusID)` in the working `ReconstructedStemID` column (engine baseline + any pins committed earlier in this pass). If a collision would result, the sweep **respects the engine's reconstruction**: the row's `ReconstructedStemID` is restored to its `ReconstructedStemID_PreSweep` value, `ReconstructionMethod` is left untouched, and the row is flagged in the **`SweepRollbackToPreSweep`** boolean column (FALSE elsewhere). A `[audit]` log line is emitted naming the tag and rollback count.

`SweepAuditOverride` continues to mark the disagreement; `SweepRollbackToPreSweep` records the chosen resolution. The canonical case is BCI tag `258411` C6, where retag-campaign reuse of `OriginalStemID = 995110` would otherwise put two distinct stems on the same `ReconstructedStemID` at C6; the engine's fresh id (`995113`) is retained instead.

### Post-engine `carried_terminal` backfill (all drivers)

After the DP/probabilistic engine returns and `maybe_add_posterior_bins()` and `apply_pin_track_rejoin()` have run, **every driver** (`main_cpp.R`, `main_cpp_chunk.R`, `main_cpp_bci.R`) calls the shared helper `apply_carried_terminal_backfill()` defined in `dp_global/R/dp_global_main.R`.

The helper finds rows where **all three** conditions hold:

- `is.na(ReconstructedStemID)` (engine produced no assignment),
- `is.na(DBH)` (no measurement to match against), and
- `Status %in% c("dead", "stem dead", "broken below", "missing")` (a terminal event for the stem, or a missing record after the stem's earlier records).

For each such row, after sorting by `(Tag, source ID, CensusID)` — the source ID is `StemID`, or `OriginalStemID` when the table has no `StemID` — it copies the most recent prior non-NA `ReconstructedStemID` from the same `(Tag, source ID)` group (LOCF) and sets `ReconstructionMethod = "carried_terminal"`. Biologically, a death/break row ends the trajectory of the most recent prior identity that shared the source ID — without this fill those rows would be dropped from any downstream trajectory.

Where it fires:

- `main_cpp.R` — Step 5.5b, after `maybe_add_posterior_bins(out)` and `apply_pin_track_rejoin(out)`.
- `main_cpp_chunk.R` — after the chunk's groups have been combined into `out_chunk`, with `verbose = FALSE` to keep multi-tag chunked logs quiet.
- `main_cpp_bci.R` — section 9b. The BCI driver also performs an upstream pre-DP `TrueStemID` propagation (Steps 1–3, see the BCI section below) that is *not* applicable to non-BCI inputs; the post-engine helper chain is shared by all three drivers.

### Post-engine `given_orphan` backfill (all drivers)

Immediately after `apply_carried_terminal_backfill()`, every driver also calls `apply_orphan_stem_backfill()` (defined alongside it in `dp_global/R/dp_global_main.R`). The helper finds rows where **all four** conditions hold:

- `is.na(ReconstructedStemID)` (engine produced no assignment, and `apply_carried_terminal_backfill()` did not fill it either),
- `is.na(TrueStemID)` (no upstream anchor),
- `is.na(DBH)` (no measurement to match against), and
- the source-id column (`StemID` if present, otherwise `OriginalStemID`) is non-NA.

For each such row it copies the source id into `ReconstructedStemID` and sets `ReconstructionMethod = "given_orphan"`. This handles "born-orphan" stems whose source identifier first appears with no DBH and no upstream anchor (e.g. a brand-new stem id first recorded as broken-below at C7+). The source id is unambiguous, the DP has no signal to disambiguate without DBH, and there is no collision risk because the DP never reached these rows. Where it fires:

- `main_cpp.R` — Step 5.5c, immediately after `apply_carried_terminal_backfill()`.
- `main_cpp_chunk.R` — applied to `out_chunk` with `verbose = FALSE` immediately after the chunk's `carried_terminal` backfill.
- `main_cpp_bci.R` — section 9b, immediately after `apply_carried_terminal_backfill()`.

`apply_terminal_to_host()` and `apply_broken_below_invariants()` follow in every driver (rules in `dp_global/README.md` and in the comments of `dp_global/R/dp_global_main.R`).

**Warning messages from the probabilistic matcher** (sample-level repair counts, ME cumulative-shrinkage breaks, growth-aware resolver diagnostics) are emitted via `message()` on stderr and are also printed to stdout via `cat()` when `DP_VERBOSE=TRUE`. To capture all warnings in a log file, redirect both streams: `Rscript ... > log.txt 2>&1`.

---

## BCI Debug Driver (`main_cpp_bci.R`)

`main_cpp_bci.R` is a debug driver for single-tag runs on BCI (Barro Colorado Island) multi-stem census data. It sources `main_cpp.R` to inherit all helper functions and CLI handling, then overrides the input file and a few BCI-specific defaults, and adds a **pre-DP `TrueStemID` reconstruction** step that is specific to BCI's identity conventions. It estimates no biological parameters: the input file already carries the `Bio_*` columns. The production BCI run uses `BCI_stem_reconstruction/2_STEM_IDENTIFICATION/1_main_cpp_chunk_bci.R`, not this script.

### Pre-DP TrueStemID propagation (BCI-specific, Steps 1–3)

Before calling the DP, the BCI driver writes `TrueStemID` on every row whose biological identity is unambiguous from the BCI database conventions:

- **Step 1a** — any row with a non-NA `StemTag` (the field crew physically tagged the stem) gets `TrueStemID = OriginalStemID`.
- **Step 1b** — any row at `CensusID >= 7` (BCI's systematic re-tagging campaign from 2010 onward) gets `TrueStemID = OriginalStemID`, except rows without a DBH whose `Status` is `"dead"`, `"stem dead"` or `"broken below"` (the end of a trajectory, matched by the engine and the post-engine backfill).
- **Step 2a/b/c** — within each `(Tag, OriginalStemID)` group, the rows after the last DBH measurement form the terminal phase. Step 2b is written to anchor measured break / resprout rows of that phase, but a row after the last DBH cannot have a DBH, so it anchors nothing. Step 2c fills the terminal-phase rows that are still NA from an anchored terminal-phase row of the same group (LOCF, then NOCB).
- **Step 3a** — any remaining row **with a DBH** whose `Status` is `"broken below"` or whose `ListOfTSM` has an R-family resprout code gets `TrueStemID = OriginalStemID` (start of a new trajectory). Death / break rows without a DBH are not anchored.
- **Step 3a.5** — same-`OriginalStemID` continuity anchor. Within each `(Tag, OriginalStemID)` group, if (a) every row is currently still unanchored (no non-NA `TrueStemID`), (b) there is at least one alive DBH-bearing row, and (c) there is at least one NA-DBH terminal row (`Status ∈ {"dead", "stem dead", "broken below"}`), then the alive DBH-bearing rows of the group are anchored to their own `OriginalStemID`. In the BCI data 99.8% of dead/stem-dead NA-DBH terminals and 55.6% of broken-below NA-DBH terminals have an earlier record with the same `OriginalStemID` (computed in `data_simulation/sample_data_BCI/general_data/dead_pattern.qmd` and `broken_below_pattern.qmd`). Without this step, an unanchored alive history can be reassigned by the DP to a competing newly-born stem and produce duplicate `ReconstructedStemID` values at later censuses (BCI tags `060145`, `233660`, `606162`, `639010`, `739002`).
- **Step 3b** — within each `(Tag, OriginalStemID)` group, if the rows that carry a DBH hold exactly one non-NA `TrueStemID` value, fill remaining NAs with that value; conflicting groups are left alone and counted in a diagnostic message.

Pre-anchor rows without an unambiguous identity remain NA and are resolved by the DP. The combination of this pre-DP propagation and the DP/probabilistic engines' hard-invariant sweep (see below) guarantees `ReconstructedStemID == TrueStemID` on every row that Steps 1–3 anchored.

### Post-DP helpers (section 9b)

After `run_dp_one_group()` returns and `maybe_add_posterior_bins()` has been applied, the BCI driver invokes the shared helpers `apply_pin_track_rejoin()`, `apply_carried_terminal_backfill()`, `apply_orphan_stem_backfill()`, `apply_terminal_to_host()`, `apply_broken_below_invariants()`, `renumber_engine_minted_ids()` and `finalize_posterior_paths()`. These are the same helpers called by `main_cpp.R` and `main_cpp_chunk.R`; see the *Post-engine `carried_terminal` backfill*, *Post-engine `given_orphan` backfill*, and *Posterior writing pipeline* sections above for the rules and rationale. Together with the script-level hard-invariant + duplicate-aware sweep (which also runs in this driver via the inherited `run_dp_one_group()`), these are the only post-engine pieces that are identical across all three drivers — the upstream Steps 1–3 above are BCI-input-specific and have no analogue in the simulator (which already ships `TrueStemID`).

### BCI-specific defaults

- `INPUT_FILE`: `data_simulation/sample_data_BCI/multistem_tags.rds` (loaded via `readRDS()`; BCI multi-stem tags with `Species`, `Tag`, `StemTag`, `OriginalStemID`, `CensusID`, `ExactDate`, `DBH` in cm, `ListOfTSM`, `HOM`, `Status` and the `Bio_*` columns; it has no `growth_form` column)
- `WHICH_TAG`: `"080297"` (a multi-stem debugging tag)
- `FORCE_ONE_SPECIES_PARAMETERS`: `FALSE` (use real BCI species)
- `MAX_GROWTH_FIXED`: `5.0`, `MAX_SHRINK_FIXED`: `-0.5`, `DP_MAX_STATES`: `1039L`, `PRUNE_BOUND_FACTOR`: `5`
- `USE_MEASUREMENT_ERROR`: `FALSE`
- `POSTERIOR_SAMPLES`: `0L` (no posterior samples unless `--POSTERIOR_SAMPLES=<n>` is given)
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
- `data_simulation/sample_data_BCI/multistem_tags.rds` (or whichever `INPUT_FILE` is set to) must exist

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
Rscript dp_global/scripts/main_cpp.R --RUN_ALL_TAGS=TRUE --PLOT_PDF_ONE_TAG_ONLY=FALSE
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

- `<RUN_DIR>/stem_reconstruction_dp_global_rcpp.csv` (main reconstruction; required)
- `<RUN_DIR>/posteriors/tag_*_posterior_samples_*_paths.csv` (posterior paths; only CSV path files are read, so the run must use `POSTERIOR_SAMPLES_FORMAT="csv"`, the default of `main_cpp.R` and `main_cpp_chunk.R`). Without them only the MAP outputs are produced. Paths are weighted by their `path_prob`.

It needs the packages `data.table`, `ggplot2`, `patchwork` and `scales`.

### Outputs

Two CSV files and one PDF written to `<RUN_DIR>/`:

| File | Rows | Description |
|------|------|-------------|
| `basal_area_tag_census.csv` | Tag × Census | Total BA (m²), stem count, year |
| `basal_area_tag_change.csv` | Tag × Census interval | BA change decomposed into survivor growth, mortality loss, and recruitment gain — MAP values plus posterior mean, SD, and 95% CI for each component (all in m²) |
| `basal_area_figures.pdf` | — | Multi-page PDF: summary page (all tags), per-tag detail pages (2×2 panels: BA trajectory, stem count, decomposition bars with uncertainty whiskers, stem demographics), posterior uncertainty histograms, and **posterior density plots** (kernel densities of Growth/Loss/Gain/DeltaBA with weighted-mean vertical lines — one page pooled across all census intervals, then one page per interval) |

### Key insight

Tag-level total BA per census is **invariant** to identity assignment — the same DBH values are summed regardless of which stem identity each observation receives. The **decomposition** of BA change into growth (surviving stems), loss (mortality), and gain (recruitment) **is** identity-dependent. Different posterior path samples assign different stems as survivors vs. deaths vs. recruits, producing uncertainty in the attribution of BA change to these demographic components.
