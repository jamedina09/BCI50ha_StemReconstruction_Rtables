# Global Dynamic Programming Stem-ID Reconstruction with Biological Costs

**Author:** José A. Medina-Vega

---

## Table of Contents

1. [Overview](#overview)
2. [Quickstart](#quickstart)
3. [Problem Description](#problem-description)
4. [Algorithm Architecture](#algorithm-architecture)
5. [Data Requirements](#data-requirements)
6. [Biological Model & Parameters](#biological-model--parameters)
7. [Core Algorithm Details](#core-algorithm-details)
8. [Uncertainty Quantification](#uncertainty-quantification)
9. [Outputs & Diagnostics](#outputs--diagnostics)
10. [Parameter Estimation](#parameter-estimation)
11. [Fallback Mechanisms](#fallback-mechanisms)
12. [ForestGEO Code Handling](#forestgeo-code-handling)
13. [Pruning & Conservative Guards](#pruning--conservative-guards)
14. [Workflows & Usage Patterns](#workflows--usage-patterns)
15. [Implementation Reference](#implementation-reference)
16. [Notes & Common Issues](#notes--common-issues)
17. [References](#references)
18. [Basal Area Uncertainty Quantification](#basal-area-uncertainty-quantification)

---

## Overview

Reconstruct stable stem identities across censuses using a biologically informed global dynamic programming solver that enforces life-cycle constraints and provides uncertainty quantification.

### Key Features

- **Global optimization** across multiple censuses
- **Biological realism** (growth, mortality, recruitment)
- **Measurement error handling** following Chave et al. (2004)
- **Exact uncertainty quantification** via posterior marginals
- **Probabilistic matcher** for the tags the DP does not solve (Gumbel-perturbed assignments per census pair, with per-observation posterior probabilities)
- **C++ implementation** of the transition costs and feasibility checks (required by the DP)

### Stem Identity Renumbering Workflow

All reconstructed stem identities (`ReconstructedStemID`) are assigned **sequentially from 1 to N within each tag**, ordered by the earliest census in which each stem appears. If multiple stems first appear in the same census, the largest DBH at that census gets the lower ID, with ties broken by original ID. This renumbering is always applied after the engine and all post-processing helpers, regardless of which driver or fallback is used. **Negative or zero IDs are never produced.** The BCI merge step's `apply_measurement_rejoin()` runs after this renumbering and does not renumber: a joined stem keeps one of its two IDs, so a tree it changes skips one number per join.

All drivers (`main_cpp.R`, `main_cpp_chunk.R`, `main_cpp_bci.R`) use the same post-engine helper chain: `maybe_add_posterior_bins()`, `apply_pin_track_rejoin()`, `apply_carried_terminal_backfill()`, `apply_orphan_stem_backfill()`, `apply_terminal_to_host()`, `apply_broken_below_invariants()`, `renumber_engine_minted_ids()`, and finally `finalize_posterior_paths()`. Posterior path files are always written after renumbering, using the staging architecture. All downstream outputs, including posterior path files, use this renumbered ID space.

---

## Quickstart

### Prerequisites

R with the packages `data.table`, `Rcpp`, `igraph` and `MASS` (checked by `dp_global/R/dp_global_main.R`), `here` (driver scripts), `lpSolve` (probabilistic matcher) and `parallel`; a C++ compiler (the DP needs the compiled functions of `dp_global/src/transition_cost_rcpp.cpp`); optional: `ggplot2`, `cowplot`, `reshape2`, `scales`, `patchwork` (plots and reports), `arrow` (Feather files), `testthat` (tests).

Run everything from the project root (the folder that contains `dp_global/`).

### Running the Code

**Interactive (single-tag or small datasets):**

```r
source("dp_global/scripts/main_cpp.R")   # defines the settings and functions; does not start a run
res <- run_main()
```

**Command line (single-tag or small datasets):**

```bash
Rscript dp_global/scripts/main_cpp.R
```

**Command line (large datasets — chunked incremental output):**

```bash
Rscript dp_global/scripts/main_cpp_chunk.R
```

Notes:

- For large datasets prefer `main_cpp_chunk.R` which processes groups (Tag + species) in configurable chunks, writes incremental CSV/RDS output, and supports resume (`DP_CHUNK_RESUME=TRUE`).
- Census timing comes from the `ExactDate` column: the interval of a census pair is the difference of the mean dates of the two censuses (per stem in the parameter estimation). There is no interval argument.
- `.rds` outputs are written by default (`--WRITE_DP_RDS=TRUE`).
- Tests: `testthat::test_file("dp_global/tests/test_estimate_bio_pars.R")` and `testthat::test_file("dp_global/tests/test_probabilistic_birth_death.R")`.

### File Structure

```
dp_global/
├── README.md
├── R/
│   ├── dp_global_main.R              # Source loader (compiles the C++ file, sources the engine modules) + post-engine helpers
│   ├── dp_global_bio.R               # Biological parameter estimation (estimate_bio_pars) + R cost breakdown
│   ├── dp_global_dp.R                # Core DP solver (match_stems_dp_global_backward_marginals_batch)
│   ├── dp_global_states.R            # State enumeration and track-DBH helpers
│   ├── dp_global_matchers.R          # Stepwise igraph matcher (not called by the DP workflow)
│   ├── dp_probabilistic_matching.R   # Probabilistic matcher (match_stems_probabilistic)
│   ├── dp_global_utils.R             # State-count helper and interval helpers
│   ├── dp_global_diag.R              # Diagnostics and PDF plotting
│   ├── measurement_rejoin.R          # Post-engine rejoin of stems split by a moved point of measurement
│   ├── naming_helpers.R              # Output directory naming helpers
│   ├── sensitivity_transition_cost_bio.R  # Sensitivity sweep helpers
│   ├── realism_calibration.R         # Realism calibration helpers
│   ├── k_tuning_viz.R                # K-tuning visualisation
│   ├── check_functions.r             # Bio-parameter inspection and plotting utilities
│   ├── complexity/
│   │   ├── estimate_dp_complexity_function.R  # DP complexity estimator
│   │   └── profile_complexity.R               # Complexity profiling
│   └── dpglobal_bundle/
│       ├── README.md                     # How to build and deploy the bundle
│       ├── dpglobal_bundle_loader.R      # Bundle builder (creates RData + manifest)
│       ├── package_bundle.sh             # Packaging helper (creates tarball in dist/)
│       ├── verify_bundle.R               # Basic smoke-test for deployed bundles
│       └── dist/                         # Generated tarballs
├── scripts/
│   ├── README.md                     # Driver options and outputs
│   ├── main_cpp.R                    # Interactive / single-tag driver
│   ├── main_cpp_chunk.R              # Chunked driver for large runs
│   ├── main_cpp_bci.R               # BCI debug driver (single-tag, RDS input, sources main_cpp.R for helpers)
│   └── basal_area_uncertainty.R      # Posterior-based basal area uncertainty quantification
├── src/
│   ├── README.md                     # C++ functions and their R wrappers
│   ├── transition_cost_rcpp.cpp      # C++ transition cost, phase feasibility and assignment functions
│   └── transition_cost_rcpp.R        # R wrapper for Rcpp-compiled functions
├── tests/
│   ├── test_estimate_bio_pars.R      # testthat: growth-SD refit and recruitment rate units
│   └── test_probabilistic_birth_death.R  # testthat: birth-death matcher and assignment solver
```

Output artifacts (not tracked by git) are written to `dp_global/output/` at runtime.

---

## Problem Description

### Context: Forest Census Data

Forest census measurements track multiple stems per tree over time:

- **Tag:** Plot/tree grouping identifier
- **CensusID:** Temporal sequence (1, 2, 3, ...)
- **DBH:** Diameter at breast height (may be missing = not observed)
- **Species:** Taxonomic identity
- **Anchor Census:** Later census with known `TrueStemID` for each observed stem

### The Challenge

**Goal:** Assign earlier observations to anchor identities such that:

1. DBH trajectories are biologically plausible
2. Identity assignments are globally optimal
3. Life-cycle constraints are respected (no resurrection)

**Why it's hard:**

- Multiple stems per tag creates combinatorial assignment space
- Missing observations require handling of recruitment/mortality
- Measurement error complicates growth assessment
- Local greedy matching can trap in suboptimal solutions

Note on measurement heights (HOM):
If DBH measurements were taken at differing heights, include a `HOM` column (or `hom`; detection is case-insensitive) and convert DBH to a common reference height (1.3 m) prior to running the workflow for taper-corrected growth forms (trees, shrubs, figs, unknown). For **non-taper-corrected** growth forms (palms, strangler figs, tree ferns) where taper correction is not applicable, the DP solver uses wide base pruning bounds (default 1.25× the standard growth/shrink limits) because these forms exhibit real biological DBH growth and, when the measurement point (HOM) changes between censuses, the recorded DBH can shift substantially even if the true diameter at a fixed height has not changed. When a `HOM` column is present, the solver additionally widens bounds per census pair in proportion to the worst-case HOM deviation from 1.3 m (controlled by `hom_tolerance_scale`, default 2.0 cm of DBH per meter of deviation, divided by the census interval). NA HOM values are treated as 1.3 m (zero deviation). If no HOM column is present, HOM widening is disabled and only the wide base bounds apply.

---

## Algorithm Architecture

### Track-Based State Space

The DP uses **K latent tracks** (identity slots):

- Each observed stem must be assigned to exactly one track (injective mapping)
- Tracks can be empty (`NA`) in any census
- Track occupancy implies identity continuity across censuses

**State at census t:** An injective assignment vector of length `n_obs(t)`, where each element indicates which of K tracks that observation occupies.

**Number of states:**
$$P(K, n_{obs}) = K \times (K-1) \times \cdots \times (K-n_{obs}+1)$$

### Understanding `max_states` (DP_MAX_STATES)

`DP_MAX_STATES` (passed as `max_states` to the solver) is the single parameter that controls when the DP solver falls back to the probabilistic matcher. It imposes two limits:

1. **Per-census enumeration (`enum_exceeded`):** If $P(K, n_{obs})$ exceeds `max_states` at any single census, the solver cannot enumerate the assignment states and falls back. This holds for a census whose observations may use every track. When pins or growth bounds restrict the tracks of a census, `enumerate_states_constrained()` is used: it stops at `max_states` and returns the first `max_states` assignments instead of signalling an overflow, so a constrained census with more feasible assignments than `max_states` is solved on a truncated state set.
2. **Inter-census transitions (`edge_count_exceeded`):** If the cross-product of state counts between two adjacent censuses reaches `max_states`² (called `max_edges`), the transition matrix is too large and the solver falls back.

In practice, the per-census enumeration limit is reached first because state counts grow factorially.

#### Fallback thresholds by `DP_MAX_STATES` value

The tables below show when the per-census enumeration forces fallback. "Stems observed" ($n$) is the number of stems with non-NA DBH in a single census. $K$ is the number of identity tracks.

**With $K = n + 1$ (minimum realistic: `slack_tracks = 1`, no births needed):**

| Stems ($n$) | Tracks ($K$) | States $P(K,n)$ | 1,000 | 20,000 | 40,000 |
|:-:|:-:|--:|:-:|:-:|:-:|
| 2 | 3 | 6 | DP | DP | DP |
| 3 | 4 | 24 | DP | DP | DP |
| 4 | 5 | 120 | DP | DP | DP |
| 5 | 6 | 720 | DP | DP | DP |
| 6 | 7 | 5,040 | fallback | DP | DP |
| 7 | 8 | 40,320 | fallback | fallback | fallback |
| 8 | 9 | 362,880 | fallback | fallback | fallback |

**With $K = n + 2$ (typical: `slack_tracks = 1` plus 1 birth track needed):**

| Stems ($n$) | Tracks ($K$) | States $P(K,n)$ | 1,000 | 20,000 | 40,000 |
|:-:|:-:|--:|:-:|:-:|:-:|
| 2 | 4 | 12 | DP | DP | DP |
| 3 | 5 | 60 | DP | DP | DP |
| 4 | 6 | 360 | DP | DP | DP |
| 5 | 7 | 2,520 | fallback | DP | DP |
| 6 | 8 | 20,160 | fallback | fallback | DP |
| 7 | 9 | 181,440 | fallback | fallback | fallback |
| 8 | 10 | 1,814,400 | fallback | fallback | fallback |

**Summary — maximum stems per census handled by exact DP:**

| `DP_MAX_STATES` | `max_edges` (= `DP_MAX_STATES`²) | Max stems ($K = n+1$) | Max stems ($K = n+2$) |
|--:|--:|:-:|:-:|
| 1,000 | 1,000,000 | 5 | 4 |
| 20,000 | 400,000,000 | 6 | 5 |
| 40,000 | 1,600,000,000 | 6 | 6 |

**Key insight:** With the default `DP_MAX_STATES = 40,000` of `main_cpp.R` and `main_cpp_chunk.R`, the DP handles tags with up to **6 observed stems per census** exactly (for these values of $K$). Tags with **7 or more stems** in any census are routed to the probabilistic matcher. (The BCI driver uses `DP_MAX_STATES = 10,000`.)

#### Inter-census transition budget

The `max_edges` limit (= `DP_MAX_STATES`²) guards the cross-product of states between adjacent censuses. For the values above, when both censuses are within the per-census state limit, the cross-product does not exceed `max_edges` — so in practice the per-census enumeration limit is the binding constraint. The edge guard becomes relevant when `DP_MAX_STATES` is set high enough to pass per-census enumeration but the product of two large censuses overflows the transition budget.

#### How to choose a value

```r
# 1. Find the most complex tags in your dataset
library(data.table)
dt <- fread("your_data.csv")
obs_per_census <- dt[!is.na(DBH), .N, by = .(Tag, CensusID)]
max_obs <- obs_per_census[, .(max_n = max(N)), by = Tag][order(-max_n)]
head(max_obs, 10)  # top 10 most complex tags

# 2. Compute states for a specific stem count
n <- 6   # max observed stems in any census
K <- 8   # n + 2 (slack + 1 birth)
states <- prod(K:(K - n + 1))  # P(8, 6) = 20,160
cat("States:", states, "\n")

# 3. Set DP_MAX_STATES above that to guarantee exact DP
# Rscript dp_global/scripts/main_cpp_chunk.R --DP_MAX_STATES=25000
```

**Trade-off:** Higher values → exact DP for more tags (slower, more memory). Lower values → more tags use the probabilistic matcher (faster, approximate; it uses the same `Bio_*` parameters and hard pruning bounds).

### Life-Cycle Phase System

Each track maintains a phase variable:

- **Phase 0 (prebirth):** Track not yet occupied, can transition to Phase 1
- **Phase 1 (alive):** Track occupied, can remain alive or transition to Phase 2
- **Phase 2 (dead / unused):** Track empty from here on (its stem died earlier, or the track is never used), cannot return to Phase 1

**Enforced constraint:** `OBS → NA → OBS` is forbidden (no resurrection)

### Choosing K (Number of Tracks)

K must accommodate:

1. Number of unique anchor stem IDs
2. Maximum observed stems in any census
3. Births needed to explain stem count increases

```text
births_needed = Σ max(0, n_obs(t+1) - n_obs(t))
K_base = max(#anchor_IDs, max_obs, n_obs(1) + births_needed)
K = min(K_base + slack_tracks, max_tracks)
K = max(K, #anchor_IDs + #TrueStemIDs pinned before the anchor that are not at the anchor)
```

**Slack tracks:** Extra tracks allowing simultaneous death+birth in constant-count intervals. With `slack_require_anchor_recruitable = TRUE` they are granted only when at least one anchor stem is no larger than `Bio_Recruit_MaxDBH_unit`. If `K` is smaller than the largest stem count of a census, the tag goes to the probabilistic matcher (`K_too_small`).

### DP Objective

Minimize total cost across all transitions:

$$\text{TotalCost} = \sum_{t=1}^{\text{anchor}-1} \text{TransitionCost}(t \to t+1)$$

Subject to:

- Injective assignments at each census
- Valid phase transitions per track
- Fixed assignments at anchor census (from `TrueStemID`)

---

## Data Requirements

### Input Columns (Required)

| Column | Type | Description |
|--------|------|-------------|
| `Tag` | integer/character | Plot/tree grouping identifier |
| `CensusID` | integer | Census sequence (1, 2, 3, ...) |
| `DBH` | numeric | Diameter in cm (`NA` = not observed) |
| `ExactDate` | Date | Measurement date; gives the census intervals |
| `species` | character | Species identifier (parameter sets; the drivers create it) |
| `TrueStemID` | integer | Known stem ID at anchor census (and at any earlier row to pin) |

Used when present: `Status` and `ListOfTSM` (resprout, broken-below, dead and missing-from-field records), `growth_form` and `HOM` (pruning bounds and routing), `Species` (routing with `prob_species`), `StemID` / `OriginalStemID` (post-engine helpers).

### Biological Parameter Columns (Required)

The parameters must be present in the dataset, with one value per tag, before running DP. These are typically added by calling `estimate_bio_pars()` and joining results. `Bio_Gamma_Growth`, `Bio_Max_Growth`, `Bio_Max_Growth_Soft` and `Bio_K_Growth` are optional (defaults 0, no bound, no cap and 0); the DP stops when any other column is missing, `NA` or not constant within the tag:

#### Growth Parameters

- `Bio_Mu_Growth` ($\mu_{\text{const}}$): Intercept of mean annual growth (cm/year)
- `Bio_Gamma_Growth` ($\mu_{\gamma}$): Log-DBH slope in mean growth model (cm/year per unit of log DBH)
- `Bio_Sigma0_Growth` ($\sigma_0$): Baseline growth process SD (cm/year)
- `Bio_Sigma1_Growth` ($\sigma_1$): Growth SD slope vs DBH ((cm/year)/cm)

#### Shrinkage Penalties

- `Bio_Max_Shrink` (`max_shrink`): Hard lower bound on annual growth (cm/year, negative)
- `Bio_K_Shrink` ($k_{\text{shrink}}$): Soft shrinkage penalty weight (1/cm²)

#### Extreme Growth Penalties

- `Bio_Max_Growth` (`max_growth`): Hard upper bound on annual growth (cm/year)
- `Bio_Max_Growth_Soft` (`max_growth_soft`): Soft growth cap (cm/year) above which the soft penalty applies
- `Bio_K_Growth` ($k_{\text{growth}}$): Soft extreme growth penalty weight (1/cm²)

#### Mortality Parameters

- `Bio_H0_Mortality` ($h_0$): Baseline hazard parameter
- `Bio_Beta_Mortality` ($\beta$): Size effect on mortality hazard

#### Recruitment Parameters

- `Bio_Recruit_Meanlog`: LogNormal meanlog for recruit size
- `Bio_Recruit_Sdlog`: LogNormal sdlog for recruit size
- `Bio_Recruit_MaxDBH_unit`: Maximum plausible recruit DBH (cm)
- `Bio_Recruitment_lambda`: Recruitment rate per empty track (1/year)

### Output Columns (Added by Solver)

| Column | Description |
|--------|-------------|
| `ReconstructedStemID` | Assigned stem identity |
| `ReconstructionMethod` | One of: `"given"`, `"dp"`, `"probabilistic"`, `"provisional_dp"`, `"dp_mf_inferred"`, `"pin_track"`, `"carried_terminal"`, `"given_orphan"`, `"terminal_to_host"`, `"measurement_rejoin"`, `"bb_split"`, `"bb_split_carry"`, `"bb_post_terminator_split"`, `"bb_post_terminator_split_carry"`, `"none_after_anchor"`, `"skipped_no_data"` — see *Primary Outputs* below |
| `ConstraintViolation` | Post-hoc diagnostic flag (DP-solved tags) |
| `DP_KUsed` | Number of tracks used |
| `DP_MaxStatesPerCensus` | Largest per-census state count $P(K, n_{obs})$ of the tag, before any track constraint |
| `DP_MaxStatesCensusID` | Census with that largest count |
| `DP_FallbackReason` | Why the tag went to the probabilistic matcher (`NA` for tags solved by the DP); see *DP fallback reason codes* |
| `obs_row_id` | Row number within the tag; the key used by the posterior path files |
| `SweepAuditOverride` | `TRUE` on rows where the hard-invariant sweep in `finalize_out()` (and the equivalent sweep in the probabilistic matcher / driver scripts) overrode an engine-assigned `ReconstructedStemID` to match a known `TrueStemID`. Use to surface engine-vs-pin disagreements that would otherwise be silently corrected. See the *Hard-invariant sweep and the `SweepAuditOverride` column* section in `dp_global/scripts/README.md`. |
| `ReconstructedStemID_PreSweep` | Snapshot of the engine's `ReconstructedStemID` *before* the hard-invariant sweep ran. It equals the engine's original choice on rows where `SweepAuditOverride == TRUE`. Post-engine helpers that relabel rows after the sweep (`apply_pin_track_rejoin()`, `apply_terminal_to_host()`, the broken-below pass) do not update it, so it can also differ from the final `ReconstructedStemID` on other rows. The column is populated once by the first sweep that fires (engine `finalize_out` → probabilistic matcher → script-level backstop) and preserved unchanged by later sweeps. |
| `SweepRollbackToPreSweep` | `TRUE` where the script-level sweep left a pin unapplied to avoid a duplicate ID in one census, or where `apply_terminal_to_host()` restored the engine's choice. |
| `DP_Chunk`, `run_out_dir` | Chunk number (chunked drivers) and name of the run folder, added by the drivers. |

Notes on post-anchor output semantics:

- When DP is scoped to pre-anchor censuses (because there are observations after the requested anchor), post-anchor rows are preserved and appended to the DP output. Post-anchor rows with a non-NA `TrueStemID` (with or without a `DBH`) will have `ReconstructedStemID = TrueStemID` and `ReconstructionMethod = "given"`.
- On tags solved by the DP, unmeasured broken-below / R-coded rows without a `TrueStemID` in the censuses that directly follow the anchor (up to the first later measurement) take the ID of the anchor stem when exactly one stem is measured at the anchor (`"given"`).
- Remaining post-anchor rows are labeled `ReconstructionMethod = "none_after_anchor"` and leave the engine with `ReconstructedStemID = NA`; the post-engine backfills may fill them later.
- Pre-anchor rows (censuses before the anchor) with non-NA `TrueStemID` are constrained to that identity when `pin_truestemid = TRUE` (default) and labeled `"given"`. With `pin_truestemid = FALSE` the solver is free to assign them another identity.
- If no anchored census with `TrueStemID` exists, and `allow_provisional_anchor = TRUE` (default), the DP can assign provisional anchor IDs at the last observed DBH census; those anchor rows are labeled `"provisional_dp"` and treated as anchors for the reconstruction.

### Posterior Uncertainty Columns (Optional)

When running `match_stems_dp_global_backward_marginals_batch()`:

- `DP_PosteriorTop1ID`, `DP_PosteriorTop1Prob`
- `DP_PosteriorTop2ID`, `DP_PosteriorTop2Prob` (one pair per `posterior_top_k`)
- `DP_PosteriorEntropy`
- `DP_PosteriorReconstructedProb`
- `DP_PosteriorUnlinkedProb` (DP-solved tags only; the probabilistic matcher leaves it `NA`)
- `DP_PosteriorBin` (if using `add_dp_posterior_bins()`)

Posterior path summaries (`*_paths.{feather|rds|csv}`) encode the `recon` column as `ObsRowID:ReconstructedStemID` pairs. The attachment helpers expect ObsRowID-based encodings.

### Posterior Paths File Format

When `POSTERIOR_SAMPLES > 0`, every tag for which posteriors are drawn ends up with one final file at:

```
<out_dir>/posteriors/tag_<Tag>_posterior_samples_<BATCH_TS>_paths.<feather|rds|csv>
```

The extension follows `POSTERIOR_SAMPLES_FORMAT` (`feather` falls back to `rds` when the `arrow` package is unavailable). When `BATCH_TS` is empty (the script default), the timestamp segment collapses to an empty string and the filename contains a double underscore — this is expected and filenames remain unique per tag.

Columns:

| Column | Type | Description |
|--------|------|-------------|
| `path_sig` | character | Dash-separated per-sample stem labels of all measured observations up to the anchor, ordered by census |
| `path_count` | integer | Number of posterior samples that produced this exact label sequence. `path_count / n_samples` is its sampled frequency (use this for sampling and Monte Carlo) |
| `path_prob` | numeric | Reference only (sums to 1 across all rows). DP: draws re-weighted by their own sampling probability (`logp`), which double counts it; probabilistic engine: equals `path_count / n_samples` |
| `recon` | character | Compact mapping of `ObsRowID:ReconstructedStemID` pairs, semicolon-separated |

Each row is a distinct **labelled reconstruction** of the tag. Posterior samples with identical label sequences are aggregated into `path_count` / `path_prob`. The labels are per-sample labels, not the exported `ReconstructedStemID` values (pins are re-encoded per tag by `apply_pins_to_samples()`), and the same grouping of observations into stems can appear in several rows under different labels (for example when a stem that ends before the anchor sits on a different spare track in two DP samples). The frequency of a grouping is therefore the sum of `path_count` over the rows that share it; for Monte Carlo use every row is simply drawn with weight `path_count`. The `ObsRowID` values in the `recon` column correspond to the `obs_row_id` column in the main reconstruction CSV, providing the join key between posterior paths and per-observation data (Tag, CensusID, OriginalStemID, DBH).

**How these files are produced (staging architecture).** The DP and probabilistic engines do *not* write the final paths files themselves. Instead, each engine stages the raw per-sample reconstruction table (in engine ID space, before any post-hoc renaming) to:

```
<out_dir>/posteriors/.staging/tag_<Tag>_samples_raw_<BATCH_TS>.rds
```

After the engine returns and the post-engine helpers run (`apply_pin_track_rejoin()`, `apply_carried_terminal_backfill()`, `apply_orphan_stem_backfill()`, `apply_terminal_to_host()`, `apply_broken_below_invariants()`, `renumber_engine_minted_ids()`), the driver calls `finalize_posterior_paths()` (defined in `dp_global/R/dp_global_main.R`). That function reads each staging file, translates `ReconstructedStemID` via the renumber mapping, applies the database pins to each sample as the export applies them after the engine (`apply_pins_to_samples()`: pin sweep with duplicate rollback, then the track rejoin of `apply_pin_track_rejoin()`), re-runs `apply_bb_invariants_to_samples()` so per-sample bb-minted IDs are derived from the renumbered track IDs, computes `path_sig` / `path_count` / `path_prob` / `recon`, writes the final paths file, and deletes the staging file on success. If `posterior_samples == 0`, no staging file is created and no paths file is written.

`path_prob` uses logp-weighted sample weights when the engine attaches a `logp` column to the staged samples (DP path); when no `logp` is present (probabilistic path), `path_prob = path_count / n_samples`.

**The two engines' draws, and how to use them together.**

- *DP:* the samples are drawn by backward sampling from the exact DP posterior, so `path_count / n_samples` already estimates the posterior probability of each path. Likely trajectories repeat (median 4 unique paths per multi-path tree; about half of all paths are drawn more than once). Because each draw was already chosen with its own probability `exp(logp)`, `path_prob` (which re-weights by `exp(logp)`) is roughly proportional to the square of the path's probability and must not be used for sampling.
- *Probabilistic engine (`dp_probabilistic_matching.R`):* the samples are approximate. Each is a noisy (Gumbel-perturbed) assignment per census pair, stitched across censuses, repaired for growth violations and filtered for pin consistency. With many stems and censuses almost every draw differs somewhere, so nearly every path has `path_count = 1` and weight `1 / n_samples` (0.005 with 200 samples). There is no single most probable path; the exported reconstruction is the most representative draw (`select_consensus_trajectory()`). The spread is a heuristic for identity uncertainty, not a calibrated posterior.
- *Pins in the samples:* engines that do not honour every pin while sampling (the DP with a provisional anchor, the pre segment of a resprout split) would otherwise give samples that group pinned rows differently from the export. `finalize_posterior_paths()` therefore applies the export's pin treatment to every sample (`apply_pins_to_samples()`), before R1/R2. Only pinned links change; unpinned links keep their per-sample variation, so the posterior spread is preserved and the export becomes one of the sampled groupings.
- *Using both together:* weight every path by `path_count / sum(path_count)` within its tree, whatever the engine. Drawing one path per tree with these weights is the same as picking one of the tree's draws uniformly, so trees from both engines can be sampled together (as `BCI_stem_reconstruction/4_EXAMPLE_STRUCTURE_ASSESSMENT/basal_area_uncertainty.R` does). `DP_PosteriorReconstructedProb` in the reconstruction table has the same meaning for both engines: the share of samples that give the observation its exported ID.

### MAP vs posterior-sampled paths 🔀

**What these two outputs represent**

- **`ReconstructedStemID` (main output)** is, for tags solved by the DP, the *MAP joint assignment* (MAP — Maximum a posteriori) decoded by the DP (a deterministic Viterbi-style backtrace of the most probable full path), and for tags solved by the probabilistic matcher the consensus sample; in both cases after the pin sweep and the post-engine helpers. This is written per-observation in the main `stem_reconstruction_*.csv`.

- **Per-path posterior summary (`*_paths.<feather|rds|csv>`)** is an *empirical* summary of full reconstructions produced by the posterior sampler (only generated when `posterior_samples > 0`). Each row is a unique path observed among draws; `path_count / n_samples` is its posterior probability (`path_prob` is a reference column, see above).

**Why they can differ**

- Posterior sampling is finite and stochastic: the MAP joint path may have non-zero posterior mass yet still not be drawn among the finite samples. Consequently, the exported grouping may *not* appear among the paths of the paths file.
- The labels differ: a path carries per-sample labels, so the exported reconstruction and a path are compared by how they group the observations, not by `path_sig`.

**Practical recommendations**

- Increase `posterior_samples` to raise the chance the MAP path is drawn and therefore present in the paths file.

**Quick check example (R)**

```r
library(data.table)
dp <- fread("dp_global/output/<run_dir>/stem_reconstruction_dp_global_rcpp.csv")
tg <- dp[Tag == 11]
# Load via the format used for the run (csv is the default of the dp_global drivers)
paths <- fread("dp_global/output/<run_dir>/posteriors/tag_11_posterior_samples__paths.csv")
grouping <- function(ids) match(ids, ids)   # same grouping <=> identical vectors
has_export <- vapply(strsplit(paths$recon, ";", fixed = TRUE), function(p) {
  kv <- do.call(rbind, strsplit(p, ":", fixed = TRUE))
  exported <- tg$ReconstructedStemID[match(as.integer(kv[, 1]), tg$obs_row_id)]
  identical(grouping(exported), grouping(kv[, 2]))
}, logical(1))
paths[has_export]  # empty => the exported grouping was not sampled
```

---

## Biological Model & Parameters

### Growth Model

**Mean growth** (size-dependent):
$\mu(D) = \alpha + \gamma \log(D)$

Where:

- $\alpha$ = `Bio_Mu_Growth` (intercept, cm/year)
- $\gamma$ = `Bio_Gamma_Growth` (slope on log-DBH, cm/year per unit of log DBH)
- If $\gamma = 0$, reduces to constant mean growth

**Process variability** (heteroskedastic):
$\sigma(D) = \sigma_0 + \sigma_1 D$

Where:

- $\sigma_0$ = `Bio_Sigma0_Growth` (baseline SD, cm/year)
- $\sigma_1$ = `Bio_Sigma1_Growth` (SD slope vs DBH, (cm/year)/cm)

### Mortality Model

**Hazard function:**
$\text{hazard}(D) = h_0 \exp(\beta D)$

**Death probability over interval:**
$P_{\text{death}}(D, \Delta t) = 1 - \exp(-\text{hazard}(D) \cdot \Delta t)$

Where:

- $h_0$ = `Bio_H0_Mortality` (baseline hazard)
- $\beta$ = `Bio_Beta_Mortality` (size effect parameter)
- Probability clamped to $[10^{-12}, 1-10^{-12}]$ for numerical stability

### Recruitment Model

**Size distribution:**
$D_{\text{recruit}} \sim \text{LogNormal}(\text{meanlog}, \text{sdlog})$

**Rate per empty track:**
$P_{\text{recruit}}(\Delta t) = 1 - \exp(-\lambda_{\text{recruit}} \cdot \Delta t)$

Where:

- `Bio_Recruit_Meanlog`, `Bio_Recruit_Sdlog`: LogNormal parameters
- $\lambda_{\text{recruit}}$ = `Bio_Recruitment_lambda` (rate per year; `estimate_bio_pars()` estimates it per established tree by default, see *Recruitment Parameter Estimation*)
- Probability clamped to $[10^{-12}, 1-10^{-12}]$ for numerical stability

### Measurement Error Model (Chave et al. 2004)

When enabled (`use_measurement_error = TRUE`), DBH observations include remeasurement noise:

$$D_{\text{obs}} = D_{\text{true}} + \varepsilon$$

**Error distribution** (2-component mixture):
$$\varepsilon \sim (1-p)\,\mathcal{N}(0, \text{SD1}(D)^2) + p\,\mathcal{N}(0, \text{SD2}^2)$$

where:
$$\text{SD1}(D) = a \cdot D + b$$

**Parameters:**

- `meas_sd1_a`: Small-error slope (a)
- `meas_sd1_b`: Small-error intercept (b)
- `meas_sd2`: Large-error SD (constant)
- `meas_p_big`: Probability of large error (p)

**Effect on likelihood:** For DBH→DBH transitions, observed growth becomes a **4-component mixture** (combining measurement errors at t₀ and t₁), evaluated via log-sum-exp.

### DBH recorded in classes (rounded down)

In BCI 1982 and 1985 (censuses 1 and 2), saplings were measured in 5 mm increments, rounded down. Sources: the CTFS R Package growth tutorial — *"The argument rnd indicates that dbhs<50 mm are rounded down to 5-mm, necessary because saplings at BCI in 1982 and 1985 were measuring in 5-mm increments"* ([CTFS tutorial: growth changes](https://ctfs.si.edu/ctfsdev/CTFSRPackageNew/index.php/web/tutorials/GrowthChange/index.html)) — and Piponiot et al. 2024 (Appendix S1 R code: 1985 stems with DBH < 5.5 cm recorded in 5 mm classes, rounded down), which `BCI_stem_reconstruction/4_EXAMPLE_STRUCTURE_ASSESSMENT/biomass_stocks_fluxes.R` follows (Section 10b). The data agree: 100% (1982) and 99.9% (1985) of DBH values below 50 mm are multiples of 5 mm, against about 20% (chance) from 1990 on; between 50 and 55 mm, 93% (1982) and 73% (1985), mostly recorded 50 mm (a floored 50–54.9 mm), against about 22% later; the mean 1985→1990 increment of small stems (5.4 mm) is about twice that of later intervals (2.0–2.5 mm), the bias of a floored 1985 value, while 1982→1985 (both floored) shows none. A recorded value below 55 mm (the classes 10, 15, …, 50 mm) is therefore treated as rounded. With the narrow Gaussian growth model (σ₀ ≈ 0.04 cm/yr for many species) a 5 mm class step read as growth is a 4σ+ event, so without the rule below the engines would split one database stem into a death and a recruit.

When a census is flagged (`dbh_round_censuses`; BCI driver `DBH_ROUND_CENSUSES = "1,2"`), a DBH $d$ below `dbh_round_max` (5.5 cm) at that census is read as a true size in $[d, d + w)$ with $w$ = `dbh_round_width` (0.5 cm). For a link with rounded flags $r_0, r_1 \in \{0, 1\}$:

$$g_{\text{mid}} = g + \frac{(r_1 - r_0)\, w/2}{\Delta t}, \qquad \sigma_{\text{eff}}^2 = \sigma(D_0 + r_0 w/2)^2 + \frac{(r_0 + r_1)\, w^2/12}{\Delta t^2}$$

and the growth likelihood (step D5, and every component of the measurement-error mixture) uses $g_{\text{mid}}$, $\mu(D_0 + r_0 w/2)$ and $\sigma_{\text{eff}}$. The hard growth limits (D1, D2) stay on the measured growth $g$: rounding changes how likely a link is, never which links are allowed. For censuses that are not flagged, and stems at or above `dbh_round_max`, the rule has no effect on the cost. The same rule is applied by the probabilistic matcher (`compute_pairwise_log_likelihood()`). Only the BCI driver sets `dbh_round_censuses`; it checks at run time that every flagged census shows the rounding in the data (`✓ DBH rounded down to 5 mm classes …`, otherwise ❌ and stop).

---

## Core Algorithm Details

### Transition Cost Function

The transition cost is the negative log-likelihood for transitioning from census t to t+1, obtained by summing per-track costs across four mutually exclusive cases. It is implemented in C++ (`transition_cost_paired_rcpp_cpp()` and `transition_cost_tracks_bio_batch_rcpp_cpp()` in `dp_global/src/transition_cost_rcpp.cpp`; the DP calls the paired version once per census pair); `transition_cost_tracks_bio_components()` in `dp_global/R/dp_global_bio.R` gives the same cost split into its terms, for diagnostics.

**Key principle:** Each track contributes independently to the total cost. The function processes all K tracks sequentially.

**Inputs:**

- `track_dbh_t`: length-K numeric vector (DBH per track at census t, NA if unoccupied)
- `track_dbh_tp1`: length-K numeric vector (DBH per track at census t+1, NA if unoccupied)
- `interval_years`: $\Delta t$ between censuses
- Biological parameters (growth, mortality, recruitment, shrinkage)
- `eps_tiebreak`: deterministic tie-break weight (default $10^{-6}$)
- `hard_penalty`: cost of a forbidden transition (default $10^6$)

**Output:** Single scalar cost per pair of track vectors (sum across all tracks + tie-break term)

### Recruitment Probability

Computed once per transition (shared across all tracks):

$P_{\text{recruit}}(\Delta t) = 1 - \exp(-\lambda_{\text{recruit}} \cdot \Delta t)$

Clamped to $[10^{-12}, 1-10^{-12}]$ for numerical stability.

### Per-Track Transition Cases

#### Case A: NA → NA (No recruitment)

Track remains empty throughout interval.

$\text{cost} += -\log(1 - P_{\text{recruit}})$

**Interpretation:** Penalizes leaving track empty based on recruitment probability. Higher recruitment rates make long empty histories more costly.

#### Case B: NA → DBH (Recruitment)

Track transitions from empty to occupied.

**Hard constraint:** If D_1 > `recruit_max_dbh` (or D_1 is not positive):
`cost += hard_penalty` ($10^6$)

**Otherwise (valid recruit):**
$\text{cost} += -\log(P_{\text{recruit}}) - \log f_{\text{LogNormal}}(D_1; \text{meanlog}, \text{sdlog})$

where:
$f_{\text{LogNormal}}(x) = \frac{1}{x \cdot \text{sdlog} \cdot \sqrt{2\pi}} \exp\left(-\frac{(\log x - \text{meanlog})^2}{2 \cdot \text{sdlog}^2}\right)$

Implementation: `dlnorm(d1, meanlog, sdlog, log=TRUE)`

#### Case C: DBH → NA (Mortality)

Track transitions from occupied to empty (disappearance/death).

**Mortality probability:**
$P_{\text{death}}(D_0, \Delta t) = 1 - \exp(-h_0 \exp(\beta D_0) \cdot \Delta t)$

Clamped to $[10^{-12}, 1-10^{-12}]$.

$\text{cost} += -\log(P_{\text{death}})$

**Interpretation:** Disappearance is "cheap" only when mortality probability is high. Low mortality at size $D_0$ makes DBH→NA costly.

#### Case D: DBH → DBH (Growth)

Track remains occupied; DBH changes from $D_0$ to $D_1$.

**Annualized growth:**
$g = \frac{D_1 - D_0}{\Delta t}$

**D1. Hard shrinkage constraint:** If g < `max_shrink`: `cost += hard_penalty`, skip to next track

**D2. Hard extreme growth constraint:** If g > `max_growth`: `cost += hard_penalty`, skip to next track

**D3. Heteroskedastic growth variance:**
$\sigma(D_0) = \max(\sigma_0 + \sigma_1 D_0, \, 10^{-6})$

**D4. Mean growth (size-dependent):**
$\mu(D_0) = \begin{cases}
\mu_{\text{const}} + \mu_{\gamma} \log(D_0) & \text{if } \mu_{\gamma} \neq 0 \text{ and } D_0 > 0\\
\mu_{\text{const}} & \text{otherwise}
\end{cases}$

**D5. Growth likelihood:**

**Without measurement error:**
$\text{cost} += \frac{(g - \mu(D_0))^2}{2\sigma(D_0)^2} + \log \sigma(D_0) + \frac{1}{2}\log(2\pi)$

(At a census flagged as rounded, $g$, $D_0$ and $\sigma$ are replaced by $g_{\text{mid}}$, $D_0 + r_0 w/2$ and $\sigma_{\text{eff}}$; see *DBH recorded in classes* above.)

**With measurement error (Chave et al. 2004):**

Four-component Normal mixture for observed growth:

Component parameters:
$\text{SD}_{\text{meas},i} = \frac{\sqrt{\text{SD}_a^2 + \text{SD}_b^2}}{\Delta t}$

where $({\text{SD}_a, \text{SD}_b})$ from:

1. Small-small: (SD1($D_0$), SD1($D_1$))
2. Small-big: (SD1($D_0$), SD2)
3. Big-small: (SD2, SD1($D_1$))
4. Big-big: (SD2, SD2)

Component weights:
$w_i \in \{(1-p)^2, (1-p)p, p(1-p), p^2\}$

Total SD per component:
$\text{SD}_{\text{tot},i} = \sqrt{\sigma(D_0)^2 + \text{SD}_{\text{meas},i}^2}$

Log-likelihood:
$\ell_i = \log w_i + \log \phi(g; \mu(D_0), \text{SD}_{\text{tot},i})$

where $\phi$ is Normal PDF.

Stable log-sum-exp:
$m = \max(\ell_1, \ell_2, \ell_3, \ell_4)$
$\text{cost} += -\left(m + \log\sum_{i=1}^4 \exp(\ell_i - m)\right)$

**D6. Soft shrinkage penalty:**
$\text{if } D_1 < D_0: \quad \text{cost} += k_{\text{shrink}} (D_0 - D_1)^2$

**Units:** $(D_0 - D_1)$ in cm, $k_{\text{shrink}}$ in 1/cm²

**D7. Soft extreme growth penalty:**
d_{1,cap} = D_0 + `max_growth_soft` * Δt
If D_1 > d_{1,cap}: `cost += k_growth * (D_1 - d_{1,cap})^2`

**Units:** Excess in cm, $k_{\text{growth}}$ in 1/cm²

**D8. Deterministic tie-break (non-biological):**

After summing all track costs:

`r_0` = rank(`track_dbh_t`, ties.method="first")
`r_1` = rank(`track_dbh_{t+1}`, ties.method="first")
`both_obs` = {tracks where both t and t+1 observed}
`cost` += `eps_tiebreak` * sum_{k in `both_obs`} |r_0[k] - r_1[k]|

(ranks are taken among the observed tracks of each census)

Default: $\varepsilon_{\text{tiebreak}} = 10^{-6}$

**Purpose:** Discourages unnecessary rank crossings when biological costs are tied.

### Backward DP Recursion

**Initialization:** Anchor census assignments fixed by `TrueStemID`

**Backward loop:** For each census from anchor-1 down to 1:

1. Enumerate all valid injective states at current census (restricted to the tracks that the pins and the growth bounds allow)
2. For each state, pair it with the reachable next states (from census t+1)
3. Derive valid phase transitions and apply the hard pruning (`derive_phase_prev_batch_rcpp()`)
4. Compute the transition cost of all remaining pairs in one call + future cost from DP table
5. Keep minimum cost and store backpointer (Viterbi), and accumulate the log-sum-exp weights used for the marginals

**Path decoding:** Follow backpointers from best start state through to anchor (among equal costs the first enumerated state wins)

**Complexity:** Managed via `max_states` parameter; falls back if exceeded

---

## Uncertainty Quantification

### Posterior Distribution

The DP defines a Gibbs-like posterior over full reconstructions:

$$P(\text{path}) \propto \exp\left(-\frac{\text{TotalCost}(\text{path})}{\tau}\right)$$

where $\tau$ is the temperature parameter (the drivers pass 1):

- $\tau < 1$: Sharper posterior (more confident, concentrates on near-MAP paths)
- $\tau > 1$: Flatter posterior (less confident, spreads mass across alternatives)

**Posterior samples (optional):** If `POSTERIOR_SAMPLES` > 0, the engines draw full-path posterior samples. Use `POSTERIOR_SAMPLES_FORMAT` to select `rds`, `feather` (arrow), or `csv` for the final path files, which are written into a `posteriors/` subdirectory under the run's output directory (or under `POSTERIOR_SAMPLES_PATH` if that option is supplied); see *Posterior Paths File Format*.

### Marginal Computation

`match_stems_dp_global_backward_marginals_batch()` computes exact marginals via:

1. **Backward pass:** Log-sum-exp recursion computing log partition function
2. **Forward pass:** Normalization to obtain marginal probabilities
3. **Per-observation summaries:** Probability of assignment to each track (anchor IDs and spare tracks)

### Posterior Summary Statistics

- **Top1/Top2 IDs and probabilities:** Most likely assignments
- **Entropy:** $H = -\sum_i p_i \log p_i$ (higher = more uncertain)
- **Reconstructed probability:** Posterior mass on chosen `ReconstructedStemID`
- **Unlinked probability:** Posterior mass on not matching any anchor ID

### Posterior Binning

`add_dp_posterior_bins()` creates categorical labels:

- **confident:** Probability of the reconstructed ID at least `confident_prob` (0.95 in the drivers)
- **ambiguous:** Everything else with a posterior
- **unlinked-likely:** Probability of being unlinked at least `unlinked_prob` (0.5 in the drivers); takes precedence over confident

---

## Outputs & Diagnostics

### Primary Outputs

1. **ReconstructedStemID:** Assigned identity for each observation
2. **ReconstructionMethod:** One of:
   - `"given"`: identity copied from `TrueStemID` (anchor row, pre-anchor row pinned to its track, post-anchor row with known `TrueStemID`, or row corrected by the hard-invariant sweep)
   - `"dp"`: Assigned by the exact DP solver
   - `"probabilistic"`: Assigned by the probabilistic matcher (`match_stems_probabilistic()`)
   - `"provisional_dp"`: Provisional anchor assigned by the DP at the last observed DBH census when the requested anchor lacked `TrueStemID` (`allow_provisional_anchor = TRUE`)
   - `"dp_mf_inferred"`: Row of a missing-from-field census, set aside before the engine and re-inserted after it (see *Missing-from-Field (MF) Detection*). It takes the ID of the track that is assigned in the nearest census before and after and has no row at that census, when exactly one such track exists; otherwise its `ReconstructedStemID` stays `NA`.
   - `"carried_terminal"`: Orphan terminal-event row (`Status` ∈ {`dead`, `stem dead`, `broken below`, `missing`} with `DBH = NA`) that the engine left without an ID. The post-engine helper `apply_carried_terminal_backfill()` (defined in `dp_global/R/dp_global_main.R` and invoked by all four driver scripts) copies the most recent prior `ReconstructedStemID` of the same `(Tag, source ID)` group via LOCF; the source ID is `StemID`, or `OriginalStemID` when the table has no `StemID` column (the simulated data). Biologically these rows close out the trajectory of the most recent prior identity carrying that source ID. Rows with no earlier record of their source ID stay `NA`.
   - `"given_orphan"`: "Born-orphan" stem row that the engine left without an ID. The source identifier (`StemID`, or `OriginalStemID` when the table has no `StemID` column) is non-NA, but `DBH`, `TrueStemID`, and the engine's `ReconstructedStemID` are all NA — typically a brand-new stem id first recorded as broken-below at C7+ with no DBH and no upstream `TrueStemID`. The post-engine helper `apply_orphan_stem_backfill()` (defined in `dp_global/R/dp_global_main.R` and invoked after `apply_carried_terminal_backfill()` by all four driver scripts) copies the source identifier into `ReconstructedStemID`.
   - `"pin_track"`: Unpinned row that the engine had linked to pinned rows of the same stem but that the TrueStemID sweep left on another ID (the sweep moves pinned rows only). Typical cases: the database renumbered the stem before 2010 and only the later StemID carries a pin (Step 3a.5), or the DP resprout split offset a pinned track. `apply_pin_track_rejoin()` (first post-engine helper) moves the whole left-behind ID onto the pin when the engine track segment holds exactly one pin value, all measured rows of that ID belong to the segment, none of its rows is pinned, and the pin has no row in those censuses. Segments end at every broken-below row with a DBH (R1) and after every NA-R barrier census, so those splits stay. Pinned rows never move; `TrueStemID` is unchanged.
   - `"terminal_to_host"`: Unmeasured `dead` / `stem dead` / `broken below` record that came *before* the first evidence of life of the identity it carried (or sat on an identity that never lives). The BCI database usually stores a stem's break record under the new StemID of the resprout measured later, so pins and the sweep put it on the resprout, where R2 would leave it as a stem that never lives. `apply_terminal_to_host()` (after `apply_orphan_stem_backfill()`, before `apply_broken_below_invariants()`) moves it to the stem that ended: (a) the engine's own pre-sweep choice when that identity lived before the record — or, if a later step relabelled that engine track (e.g. `apply_pin_track_rejoin()`), the current ID of the track's last measured row before the record; the row then also gets `SweepRollbackToPreSweep = TRUE`; (b) otherwise the only identity of the tag whose last evidence of life is in the census just before, when exactly one such record group starts then. A record is never moved onto an identity that has a row in one of its censuses, measured rows never move, and ambiguous cases stay as they are.
   - `"measurement_rejoin"`: Row of a stem that the engine had ended only because the next measurement of the same trunk fell outside the hard growth bounds after its point of measurement moved, or because one diameter was recorded wrongly (BCI: 1982 diameters taken around buttresses, POM raised in later censuses; the taper correction cannot remove buttress flare). The pruning removes such a link before costs are compared while an oversized recruit is only penalised, so the engine books a death plus an impossible recruit. `apply_measurement_rejoin()` (`dp_global/R/measurement_rejoin.R`, run by the BCI merge script; idempotent) joins an ended stem to the stem of the same tree that starts in the next census when the new stem is an impossible recruit (it starts at or above the recruit limit in a tree measured before), the ended stem is the only stem of the tree whose measurements end in the census before (the only candidate), and the link lies outside the hard growth bounds; the later measurement has no break/resprout code (R1), the earlier stem was >= 10 cm, the taper-corrected size ratio is between 0.4 and 1.5, the end/start is clean (no census gets two rows) and the two stems carry at most one pin. Database StemIDs are not used (before the anchor they are the identity the reconstruction replaces). A split the engine chose inside the bounds is left as it is. The relabelled stem takes the id of the stem carrying the pin, else of the later stem (IDs are not renumbered afterwards, so the tree's IDs skip one number per join); the export joins are written to `measurement_rejoin_audit.csv`, and their observation pairs (`measurement_rejoin_pairs.csv`) are joined by `apply_measurement_rejoin_to_paths()` in each posterior sample where they are split with a clean end/start and at most one pin.
   - `"bb_split"`: Measured broken-below row (`Status == "broken below"` with a DBH) that had an earlier row on its reconstructed trajectory (**R1**, split-on-break). The post-engine helper `apply_broken_below_invariants()` (defined in `dp_global/R/dp_global_main.R`) starts a new `ReconstructedStemID` at that row. See *Broken-Below Invariants* below.
   - `"bb_split_carry"`: Later rows of the same trajectory that carry the ID started at a `"bb_split"` row (up to the next measured broken-below row). The separate label lets downstream code tell them from the splitting row.
   - `"bb_post_terminator_split"`: First measured row that followed a stump (`broken below` without a DBH) on the same reconstructed trajectory (**R2**, terminate-on-stump). `apply_broken_below_invariants()` starts a new `ReconstructedStemID` at that row.
   - `"bb_post_terminator_split_carry"`: Later rows of the same trajectory that carry the ID started at a `"bb_post_terminator_split"` row.
   - `"none_after_anchor"`: Post-anchor row with no `TrueStemID` and no engine assignment
   - `"skipped_no_data"`: Tag/segment skipped because there were no usable observations
3. **SweepAuditOverride:** Boolean — `TRUE` flags rows where the hard-invariant sweep overrode an engine-assigned `ReconstructedStemID` to enforce a known `TrueStemID`. See *Hard-invariant sweep and the `SweepAuditOverride` column* in `dp_global/scripts/README.md` for downstream-uncertainty implications.
4. **ReconstructedStemID_PreSweep:** Snapshot of the engine's `ReconstructedStemID` before the sweep. On `SweepAuditOverride == TRUE` rows it carries the original engine choice that the sweep overrode (or, when `SweepRollbackToPreSweep == TRUE`, the value that the sweep rolled back *to*). The post-engine helpers that relabel rows afterwards do not update it.
5. **SweepRollbackToPreSweep:** Boolean — `TRUE` flags rows where the script-level sweep refused to pin `ReconstructedStemID := TrueStemID` because doing so would have produced a duplicate `ReconstructedStemID` at the same `(Tag, CensusID)`, and the engine's `ReconstructedStemID_PreSweep` value was retained instead, or where `apply_terminal_to_host()` restored the engine's own choice for an unmeasured terminal record the sweep had pinned elsewhere (`ReconstructionMethod = "terminal_to_host"`). `SweepAuditOverride` will also be `TRUE` on duplicate-rollback rows (the disagreement was detected); `SweepRollbackToPreSweep` records the chosen resolution. `ReconstructionMethod` is left untouched on a rollback so it still reflects the engine's original assignment path.

### Diagnostic Outputs

- **ConstraintViolation:** Post-hoc flag for growth violations (`add_constraint_violation()`: annual growth of a reconstructed stem outside `[min_growth, max_growth]`)
- **DP_KUsed:** Actual number of tracks used
- **DP_MaxStatesPerCensus:** Worst-case state space size, $P(K, n_{obs})$ of the census with most stems
- **DP_MaxStatesCensusID:** Census achieving maximum states
- **DP_FallbackReason:** Why the tag was solved by the probabilistic matcher (see *DP fallback reason codes*)
- **`attr(out, "DP_PruneInfo")`:** Pruning counts and effective bounds (see *Diagnostics & reproducibility*)

### Broken-Below Invariants (Post-Engine Pass)

Two life-cycle invariants are enforced **after** the DP / probabilistic engine and **after** `apply_pin_track_rejoin()`, the two backfills and `apply_terminal_to_host()`, by `apply_broken_below_invariants()` (defined in `dp_global/R/dp_global_main.R`). The pass is deterministic, idempotent, and operates per reconstructed trajectory: the rows of one `(Tag, ReconstructedStemID)` in census order. The trajectory, not the `StemTag`, is the unit because StemTags only exist from the 2010 census: a StemTag group would keep a stem's pre-2010 rows apart from its later rows, so a break measured in 2010 would never be compared with the pre-2010 trunk it continues.

- **R1 (split-on-break):** A row with `Status == "broken below"` **and a DBH** that has an earlier row on its trajectory starts a new identity. A broken stem measured again is the resprout, a new physical stem. The new ID is carried forward through the later rows of the trajectory until the next measured broken-below row, which starts another new ID.
- **R2 (terminate-on-stump):** A row with `Status == "broken below"` and **no DBH** (a stump) ends its trajectory. The first later row of the same trajectory that has a DBH starts a new identity, carried forward as in R1. Unmeasured rows between the stump and that measurement (later stump, `dead` or `stem dead` records) keep the old ID.

`dead` and `stem dead` rows trigger neither rule. New IDs are taken above the largest `ReconstructedStemID` of the table (the drivers renumber every tag to 1..N afterwards, see *Stem Identity Renumbering Workflow*). The row that starts the new ID is tagged `"bb_split"` (R1) or `"bb_post_terminator_split"` (R2); the later rows that carry it are tagged `"bb_split_carry"` / `"bb_post_terminator_split_carry"`.

**Pin-override semantics.** A pin does not protect a row: `apply_broken_below_invariants()` also splits rows where `ReconstructedStemID = TrueStemID` (e.g. the BCI drivers pin measured broken-below rows to their own database StemID, which the database sometimes reuses from the trunk). On the relabelled rows whose `TrueStemID` equals the replaced ID, `TrueStemID` is rewritten to the new ID so the two columns stay consistent; the function logs the number of such overrides. `ReconstructedStemID_PreSweep` keeps the engine's value. (The script-level `TrueStemID` backstop sweep is intentionally **not** re-run after the BB pass.)

**Posterior samples.** The same R1/R2 operator is applied to each posterior sample by `apply_bb_invariants_to_samples()` (`dp_global/R/dp_global_main.R`). `finalize_posterior_paths()` calls it after it has translated the staged labels with the renumber mapping and applied the database pins per sample (`apply_pins_to_samples()`). The path signature (`paste0(ReconstructedStemID, collapse = "-")` per sample) is computed *after* this relabel, so samples with the same labelled reconstruction are counted together and `sum(path_prob) == 1` holds.

### Visualization

**Plot reconstructed trajectories:**

```r
plot_tag_to_pdf(
  out,
  pdf_file = "output/trajectories.pdf",
  include_reference = TRUE,   # adds the TrueStemID panel (needs a TrueStemID column)
  tag = NULL                  # one tag, or NULL for one page per tag
)
```

**Sensitivity analysis:**

```r
source("dp_global/R/sensitivity_transition_cost_bio.R")
# build_all_sweeps(), plot_all_sweeps_to_pdf(): cost of fixed scenarios along each parameter
```

**Realism calibration:**

```r
source("dp_global/R/realism_calibration.R")
# realism_report_from_reconstruction(): summary, per-group table and tuning suggestions
```

### Files written by the driver scripts

When you run `dp_global/scripts/main_cpp.R` or `dp_global/scripts/main_cpp_chunk.R`, outputs are written to a run-specific directory under `dp_global/output/` (the driver prints `out_dir` on startup). Common files written by the run include:

- `run_started.txt`, `run_finished.txt` — simple markers indicating job start/finish timestamps
- `run_parameters_full.txt` — full parameter dump and command-line overrides for reproducibility
- `run_log.txt` — appended log lines from `log_msg()` throughout the run
- `stem_reconstruction_dp_global_rcpp.csv` — main reconstruction CSV (written when `--WRITE_DP_CSV=TRUE`; chunked runner appends one chunk at a time)
- `stem_reconstruction_dp_global_rcpp_chunk_NNN.rds` / `*.feather` — per-chunk binary outputs from the chunked runner; `*_done.txt` companion markers gate `DP_CHUNK_RESUME=TRUE`, and a chunk that raised an error leaves `*_failed.txt`
- `stem_reconstruction_dp_global_rcpp.rds` / `.feather` — binary copies of the full reconstruction (non-chunked runs only, written when `--WRITE_DP_RDS=TRUE` / `--WRITE_DP_FEATHER=TRUE`)
- `stem_reconstruction_dp_global_rcpp.pdf` / `stem_reconstruction_dp_global_rcpp_chunk_NNN.pdf` — per-tag reconstruction plots (if enabled)
- `posteriors/tag_<Tag>_posterior_samples_<BATCH_TS>_paths.<feather|rds|csv>` — final per-tag posterior path summaries (written when `POSTERIOR_SAMPLES > 0`)
- `posteriors/.staging/` — internal staging directory for raw per-sample reconstructions emitted by the engines. `finalize_posterior_paths()` consumes and deletes each staging file after writing the corresponding final paths file; a healthy completed run leaves this directory empty
- `bio_pars_report.pdf` — plots of the estimated biological parameters (`main_cpp.R`)
- `tag_<which_tag>_realism_summary_rcpp.csv`, `tag_<which_tag>_realism_by_tag_rcpp.csv`, `tag_<which_tag>_realism_tuning_suggestions_rcpp.csv` — realism report tables when `--RUN_REALISM_REPORT=TRUE` (`main_cpp.R`)
- `simulated_all_transition_cost_sweeps_rcpp.rds`, `simulated_all_transition_cost_sweeps.csv`, `simulated_all_transition_cost_sweep_jumps_rcpp.csv`, `simulated_all_transition_cost_jumps_rcpp.csv`, `simulated_all_transition_cost_jumps_rcpp.rds` — sensitivity sweep outputs when `--SENSITIVITY_MODE` enables write (`main_cpp.R`); `simulated_all_transition_cost_sweeps_rcpp.pdf` with `--SENSITIVITY_MODE=run+write+pdf`
- `k_sweep_join_vs_split_demo_rcpp.pdf` — k-sweep demo plot when `--RUN_K_SWEEP_DEMO=TRUE` (`main_cpp.R`)

Note: enabling `--WRITE_DP_RDS=TRUE` (or `--WRITE_DP_FEATHER=TRUE`) is recommended when you plan to post-process reconstructions in R (it preserves types and attributes without re-parsing CSVs). The directory name is built from `--BATCH_TS` and the run settings; to write into a directory of your choice (for example to resume a run) pass `--OUT_DIR_OVERRIDE=<path>`. See `dp_global/scripts/README.md` for all options.

---

## Parameter Estimation

### Overview

`estimate_bio_pars()` (`dp_global/R/dp_global_bio.R`) derives the biological parameters from stems whose identity is known (`TrueStemID`). The drivers call it once per species and attach the result to every row as `Bio_*` columns. Each guardrail and soft-penalty weight has its own source argument (`recruit_max_source`, `max_shrink_source`, `max_growth_source`, `k_shrink_source`, `k_growth_source`):

- **`"data"`:** estimated from the observations
- **`"fixed"`:** a user-specified constant (`*_fixed`)

**Input:** a table with `Tag`, `CensusID`, `DBH` (cm), `TrueStemID`, `ExactDate` and optionally `species`. Only rows with `DBH > 0` and a non-NA `TrueStemID` are used.

**Censuses used (`anchor_start_census`):** the estimation grid holds every stem at every census with `CensusID >= anchor_start_census`, so only the census pairs from that census on enter the growth, mortality and recruitment estimates. The default reads the global `ANCHOR_START_CENSUS` of the calling script.

**Intervals:** the function has no interval argument. The interval of a census pair is computed per stem from `ExactDate` (date at $t_1$ minus date at $t_0$, in years of 365.25 days). A stem without a record at a census takes its tag's mean date at that census, or else the mean date of the census. The intervals of the growth pairs are returned in `res$interval$per_pair_intervals`; `res$interval$pairs_candidate_count` is the number of growth pairs plus mortality at-risk rows.

**Errors:** the function stops when no row is usable, when fewer than two censuses or five growth pairs are available, and when a `"fixed"` source or an `enforce_*` option is given an invalid value.


### Growth Parameter Estimation

**Data filtering:** Keep rows with non-NA `DBH > 0` and non-NA `TrueStemID`, on the grid of censuses from `anchor_start_census`

**Wide format conversion:**

```r
dw <- dcast(anchor_data_complete, Tag + TrueStemID + species ~ CensusID, value.var = "DBH")       # DBH
iw <- dcast(anchor_data_complete, Tag + TrueStemID + species ~ CensusID, value.var = "ExactDate") # dates
```

**Annual growth increments:**
For each adjacent census pair $(t_0, t_1)$:
$g = \frac{D_1 - D_0}{T_i}$ where $T_i$ is the stem's own interval (years) for that pair, from its two dates (see *Intervals* above).

**Mean model fitting:**
$\mu(D) = \alpha + \gamma \log(D)$

Fits `lm(g ~ log(d0))` when:

- At least 10 valid observations exist
- Variance of `log(d0)` > 0

Otherwise falls back to constant mean: $\alpha = \bar{g}$, $\gamma = 0$

**Predicted mean for each observation:**
$\mu_{\text{pred}}(i) = \begin{cases}
\alpha + \gamma \log(D_{0,i}) & \text{if fit successful and } D_{0,i} > 0\\
\alpha & \text{otherwise}
\end{cases}$

**Process SD estimation:**

Uses robust SD proxy from absolute residuals:
$\text{residual}_i = |g_i - \mu_{\text{pred}}(i)|$

For Normal$(0, \sigma^2)$: $\mathbb{E}[|X|] = \sigma\sqrt{2/\pi}$, therefore:
$\hat{\sigma}_{\text{total},i} = \text{residual}_i \times \sqrt{\pi/2}$

**Measurement error correction (if enabled):**

Expected measurement variance for annualized increment:
$\text{Var}_{\text{meas}}(g_i) = \frac{\text{Var}(\varepsilon | D_{0,i}) + \text{Var}(\varepsilon | D_{1,i})}{\Delta t^2}$

where mixture variance at diameter $D$:
$\text{Var}(\varepsilon | D) = (1-p)\text{SD1}(D)^2 + p \cdot \text{SD2}^2$
$\text{SD1}(D) = \max(a \cdot D + b, 10^{-6})$

Process SD in quadrature:
$\hat{\sigma}_{\text{proc},i} = \sqrt{\max(\hat{\sigma}_{\text{total},i}^2 - \text{Var}_{\text{meas}}(g_i), 10^{-8})}$

**Without measurement error:**
$\hat{\sigma}_{\text{proc},i} = \max(\hat{\sigma}_{\text{total},i}, 10^{-6})$

**Linear model for heteroskedasticity:**
Fits `lm(sd_proc_hat ~ d0_all)`:

- slope ≥ 0: $\sigma_0$ = intercept (at least 0.01), $\sigma_1$ = slope
- slope < 0 (the SD falls with size, e.g. palms): a line cannot be kept, since it would reach zero for large stems, and keeping its intercept with the slope set to 0 would give every stem the SD extrapolated to DBH 0. The SD is therefore refitted as a constant, $\sigma_0$ = mean of the per-pair proxies (at least 0.01) and $\sigma_1 = 0$ (tested in `dp_global/tests/test_estimate_bio_pars.R`).

Final model:
$\sigma(D_0) = \sigma_0 + \sigma_1 D_0$

### Optional enforcement: user-specified growth bounds

You can optionally enforce user-specified annual growth bounds in `estimate_bio_pars()` to remove extreme increments prior to parameter estimation. These options are independent of the guardrails returned by the function and affect only the *data used to fit* the growth mean and variance.

- `enforce_growth_bounds` (logical, default `FALSE`): when `TRUE`, observations with annualized growth $g$ (cm/year) outside the provided bounds will be dropped before fitting the mean and variance models.
- `growth_min_fixed`, `growth_max_fixed` (numeric, cm/year): one-sided bounds are allowed (set one of them to `NA` to have only a single-sided filter).

Behavior:

- Dropped observations produce a warning stating how many points were removed.
- If enforcing the bounds causes too few growth observations to remain (the function requires at least 5 total growth observations), `estimate_bio_pars()` will raise an error.

Example (enforce growth between -0.5 and 7.5 cm/year):

```r
estimate_bio_pars(x, enforce_growth_bounds = TRUE, growth_min_fixed = -0.5, growth_max_fixed = 7.5)
```

### Penalty Parameter Estimation

#### Soft Shrinkage Penalty ($k_{\text{shrink}}$)

**Data mode with measurement error** (`use_measurement_error = TRUE` and at least 10 growth pairs):

Typical measurement SD for raw DBH difference:
$s_{\text{typ}} = \text{median}\left(\sqrt{\text{Var}(\varepsilon|D_0) + \text{Var}(\varepsilon|D_1)}\right)$

Penalty weight:
$k_{\text{shrink}} = \frac{1}{2 s_{\text{typ}}^2}$

Clamped to $[10^{-6}, 10^6]$ for numerical stability.

**Data mode without measurement error** (or with fewer than 10 growth pairs):

Uses variance of observed shrinkage increments (in cm):
$\Delta D_{\text{shrink}} = D_0 - D_1 \quad \text{for pairs where } D_1 < D_0$

If at least 5 shrinkage events exist:
$k_{\text{shrink}} = \frac{1}{2 \cdot \text{Var}(\Delta D_{\text{shrink}})}$

Otherwise defaults to 50.

**Fixed mode:**

```r
k_shrink_source = "fixed"
k_shrink_fixed = 50  # or 0 to disable
```

**Interpretation:** To make shrinkage of $s$ cm cost ≈1 unit:
$k \approx \frac{1}{2s^2}$

Example: K=50 corresponds to $s \approx \sqrt{1/(2 \times 50)} = 0.1$ cm

**Application in cost function:**
$\text{cost} += k_{\text{shrink}} (D_0 - D_1)^2 \quad \text{when } D_1 < D_0$

#### Hard Shrinkage Guardrail (`max_shrink`)

**Data mode:**

Empirical lower quantile (default 0.1%):
`max_shrink_data` = Q_{0.001}(g_all)

**Measurement-only lower quantile (if measurement error enabled):**

Four-component mixture for $(e_1 - e_0)/\Delta t$:

Component SDs:
`SD_meas_i` = sqrt(SD_a^2 + SD_b^2) / Delta_t

where (SD_a, SD_b) ∈ {(SD1(D_0), SD1(D_1)), (SD1(D_0), SD2), (SD2, SD1(D_1)), (SD2, SD2)}

Component weights: {(1-p)^2, (1-p)p, p(1-p), p^2}

Solve for quantile q = `shrink_hard_prob` (default 10^{-4}) via:
CDF(x) = sum_{i=1}^4 w_i Phi(x; 0, SD_meas_i)

using typical diameter D_typ = median(D_0) and, as Delta_t, the median interval of the last census pair.

**Final value:**

- With measurement error: `max_shrink` = min(`max_shrink_data`, `max_shrink_meas`)
- Without: `max_shrink` = `max_shrink_data`

**Fixed mode:**

```r
max_shrink_source = "fixed"
max_shrink_fixed = -0.5  # cm/year (any finite number is accepted; a shrink bound is negative)
```

**Application in cost function:**
If g < `max_shrink`: `cost += hard_penalty` ($10^6$)

#### Soft Extreme Growth Penalty ($k_{\text{growth}}$)

**Data mode with measurement error** (`use_measurement_error = TRUE` and at least 10 growth pairs): Same calibration as $k_{\text{shrink}}$:
$k_{\text{growth}} = \frac{1}{2 s_{\text{typ}}^2}$

Clamped to $[10^{-6}, 10^6]$.

**Data mode without measurement error** (or with fewer than 10 growth pairs):

Uses variance of observed positive increments (in cm):
$\Delta D_{\text{pos}} = D_1 - D_0 \quad \text{for pairs where } D_1 > D_0$

If at least 5 positive events exist:
$k_{\text{growth}} = \frac{1}{2 \cdot \text{Var}(\Delta D_{\text{pos}})}$

Otherwise defaults to 50.

**Fixed mode:**

```r
k_growth_source = "fixed"
k_growth_fixed = 50  # or 0 to disable
```

**Interpretation:** Same rule as shrinkage: $k \approx 1/(2s^2)$

**Application in cost function:**
d_{1,cap} = D_0 + `max_growth_soft` * Δt
If D_1 > d_{1,cap}: `cost += k_growth * (D_1 - d_{1,cap})^2`

#### Hard Extreme Growth Guardrail (`max_growth`)

**Data mode:**

Empirical upper quantile (default 99.9%):
`max_growth_data` = Q_{0.999}(g_all)

**Measurement-only upper quantile (if measurement error enabled):**

Uses same four-component mixture as shrinkage, solving for the upper quantile $q = 1 -$ `growth_hard_prob` (default $1 - 10^{-4}$, the 99.99th percentile).

**Final value:**

- With measurement error: `max_growth` = max(`max_growth_data`, `max_growth_meas`)
- Without: `max_growth` = `max_growth_data`

**Soft growth cap:**
`max_growth_soft` = min(`max_growth`, Q_{0.99}(g_all))

**Fixed mode:**

```r
max_growth_source = "fixed"
max_growth_fixed = 7.5  # cm/year, must be positive
```

**Application in cost function:**
If g > `max_growth`: `cost += hard_penalty` ($10^6$)

### Mortality Parameter Estimation

**Data preparation:**
For each adjacent census pair $(t_0, t_1)$:

- At-risk stems: observed at census $t_0$ (non-NA DBH)
- Died indicator: 1 if NA at $t_1$, 0 if still observed (a stem that is not measured at $t_1$ counts as dead, also when it is measured again later)
- $\Delta t$: the stem's own interval for the pair

**Model:**
$\text{hazard}(D) = h_0 \exp(\beta D)$
$P_{\text{death}}(D, \Delta t) = 1 - \exp(-\text{hazard}(D) \cdot \Delta t)$

Probability clamped to $[10^{-12}, 1-10^{-12}]$ for numerical stability.

**Estimation:** Maximum likelihood via `optim()`:

Negative log-likelihood:
$-\sum_i \left[\text{died}_i \cdot \log P_{\text{death}}(D_{0,i}) + (1-\text{died}_i) \cdot \log(1-P_{\text{death}}(D_{0,i}))\right]$

Optimization on unconstrained scale:

- Parameter vector: $[\log(h_0), \beta]$
- Starting values: `mortality_start = c(log(0.01), 0)`
- Method: BFGS

**Output:** $h_0 = \exp(\text{par}[1])$, $\beta = \text{par}[2]$

### Recruitment Parameter Estimation

**Identifying recruitment events:**
For each adjacent census pair:

- At-risk tracks: NA at $t_0$
- Recruited: NA at $t_0$ AND observed positive DBH at $t_1$

**Size distribution:**

Fit LogNormal to recruited DBH values using `MASS::fitdistr`:
$D_{\text{recruit}} \sim \text{LogNormal}(\text{meanlog}, \text{sdlog})$

Fallback (if < 2 recruits): meanlog = log(2), sdlog = 0.5

**Maximum recruit DBH:**
Upper quantile guardrail (`recruit_max_source = "data"`, default 99.9%):
`recruit_max_dbh` = Q_{0.999}(D_recruited)

Fallback (no recruit observed): 5 cm. With `recruit_max_source = "fixed"` the value is `recruit_max_fixed`.

**Optional enforcement: recruit max DBH**
`estimate_bio_pars()` supports an optional pre-fit filter for recruit sizes:

- `enforce_recruit_max` (logical, default `FALSE`): when `TRUE`, recruits with DBH strictly greater than `recruit_max_fixed` (cm) are removed from the sample prior to fitting the lognormal recruit-size distribution.
- When `enforce_recruit_max = TRUE`, you must provide a positive finite `recruit_max_fixed` value. The function will warn if any recruits were dropped.

Example (enforce a 37.5 cm recruit cap):

```r
estimate_bio_pars(x, enforce_recruit_max = TRUE, recruit_max_source = "fixed", recruit_max_fixed = 37.5)
```

**Recruitment rate** (`recruit_rate_unit`, driver flag `RECRUIT_RATE_UNIT`):

- `"tree"` (default): new stems per established tree per year,
  $\lambda_{\text{recruit}} = \frac{n_{\text{new}} + 0.5}{\sum_{\text{trees, intervals}} \Delta t}$,
  where a tree is established in an interval when it has a measured stem at its start, $n_{\text{new}}$ counts the stems of those trees without a DBH at the start and with one at the end, and $\Delta t$ is each established tree's interval (the mean over its stems). The 0.5 pseudo-recruit keeps the rate positive for sets without new stems in established trees. This is the probability per year that a tree gains a stem, which is what both engines use it for (the DP's empty track, the matcher's recruitment cell).
- `"slot"`: recruits per empty slot per year, $\lambda = n_{\text{recruits}} / \sum \Delta t$ over the grid cells without a DBH at the start. Those cells are dead stems and stems that recruit later, so the ratio barely depends on the species: about 0.08/yr for nearly every BCI set (5–95%: 0.048–0.097), while recruitment observed in 2010–2023 ranges from 0.008 to 0.74 new stems per stem-year (correlation 0.27). With `"tree"` the BCI rates range from 0.0011 to 0.076 (5–95%; median 0.010; *Oenocarpus mapora* 0.149, *Hybanthus prunifolius* 0.019, *Faramea occidentalis* 0.0065). This value is also used when the per-tree rate cannot be computed (no established tree).

Both values are returned (`recruitment$lambda_tree`, `recruitment$lambda_slot`, with `lambda_unit`, `n_new_stems_established_trees` and `established_tree_years`).

### Quantile Selection

| Parameter | Default Quantile | Purpose |
|-----------|------------------|---------|
| `max_shrink_data` | 0.1% (0.001) | Lower tail for hard shrink constraint |
| `max_growth_data` | 99.9% (0.999) | Upper tail for hard growth constraint |
| `max_growth_soft` | 99% (0.99) | Upper tail for soft growth penalty |
| `recruit_max_dbh` | 99.9% (0.999) | Maximum plausible recruit size |

All quantiles are configurable via `estimate_bio_pars()` arguments (`shrink_data_quantile`, `growth_data_quantile`, `growth_soft_quantile`, `recruit_max_quantile`; the measurement-noise tail probabilities are `shrink_hard_prob` and `growth_hard_prob`).

---

## Fallback Mechanisms

### When Fallback Occurs

The DP function can be told to **bypass the solver** for whole groups of tags:

- **Growth form** (`fallback_growth_forms`, driver flag `DP_FALLBACK_GROWTH_FORMS`): a tag whose rows up to the anchor hold at least one `growth_form` value of the list goes to the matcher with reason `growth_form_forced`.
- **Species** (`prob_species`, driver flag `PROB_SPECIES`): a tag whose `Species` column holds a listed value goes to the matcher with reason `species_forced_probabilistic`.

Both arguments take a character vector or one comma- or semicolon-separated string, which the DP function splits.

The DP solver also falls back to the probabilistic matcher when:

1. No DBH is observed up to the anchor, or no usable anchor exists (see *Anchor fallback behavior*)
2. A census has too many injective states for `max_states` → reason: `enum_exceeded` (see *Understanding `max_states`* for the constrained-enumeration caveat)
3. The cross-product of adjacent census states reaches `max_states²` → reason: `edge_count_exceeded`
4. K is insufficient ($K < \max$ observed stems) → reason: `K_too_small`
5. The backward recursion, the decoding or the forward pass finds no feasible state or transition
6. The solver raises an error (handled by the drivers, see *Error fallback*)

**Fallback routing:** The `do_fallback()` helper routes all fallback reasons to `match_stems_probabilistic()`. The matcher receives the same growth bounds (`min_growth`, `max_growth`) and prune bounds (`prune_min_growth`, `prune_max_growth`, `prune_recruit_max_dbh`) as the DP solver. It also receives the `use_bio_hard_shrink_in_prob` and `use_bio_hard_growth_in_prob` flags (propagated from `USE_BIO_HARD_SHRINK_IN_PROB` and `USE_BIO_HARD_GROWTH_IN_PROB`), which control whether the bio hard gates and ME cumulative-shrinkage check are active in the probabilistic matcher. After the matcher returns, `do_fallback()` applies two rules to its output:

- **NA-R barrier:** a census without a measured stem that holds an unmeasured R-coded / broken-below row cuts every track that crosses it (see *ForestGEO Code Handling*).
- **R-boundary splitting** (`split_live_resprout_tracks()`): a measured row with a resprout code or `Status == "broken below"` starts a new stem, so only that stem's track is split, and its rows **before** the R census get a new ID. The R row and its continuation keep the matcher's ID, which is the anchor `TrueStemID` that the final sweep restores. (If the rows after the R census were renamed instead, that sweep would re-join the anchor to the stem that broke and leave the censuses in between as one-census stems counted twice in stock.) The tree's other stems keep the matcher's links, as they do in the probabilistic posterior samples, and a track is not split when one of its rows before the break is pinned. The anchor census is left to `apply_broken_below_invariants()` (R1). Unlike the DP solver's resprout segment split, which ends every stem of the tree at its first R census, the fallback splits only the broken stem.

**Error fallback:** `run_dp_one_group()` (in the driver scripts) wraps the DP call in a `tryCatch` block. If the solver throws a runtime error (e.g., memory exhaustion), the error handler runs the probabilistic matcher with the same R-boundary splitting (`split_live_resprout_tracks()`) rather than returning `NA` or crashing the run. The error message is logged via `log_msg()` and stored in `DP_FallbackReason` as `error:<message>` (first 200 characters).

### Cross-Product Edge Guard

Before evaluating backward transitions, the DP checks whether the cross-product of states between adjacent censuses reaches `max_states²`. If `n_states(t) × n_reachable_states(t+1) >= max_states²`, the solver triggers `do_fallback("edge_count_exceeded")` which routes to the probabilistic matcher. This means `max_states` is the single parameter controlling both per-census enumeration limits and inter-census transition budgets.

### Segment Split Consistency

When a tag is split into sub-problems at an R-event boundary (resprout segment split), the solver first evaluates the **whole-tag** state space ($P(K, n_{obs})$ of every census, with the pre-split K before slack tracks). If the whole-tag space already exceeds `max_states`, both sub-calls are forced to `max_states = 0`, guaranteeing both segments use the probabilistic matcher and avoiding a situation where the smaller sub-segment K would otherwise pass the threshold and use the DP inconsistently.

However, when the whole-tag check passes, each sub-segment is solved **independently** with the original `max_states` and computes its own K from its sub-segment anchor and stem counts. It is therefore possible for the pre-resprout and post-resprout segments to use **different algorithms** — for example, the post-resprout segment (anchored at the known census with fewer stems) solving cleanly with the exact DP while the pre-resprout segment (anchored at an earlier census with more stems and a larger K) falls back to the probabilistic matcher. Both segments' outputs are combined into a single reconstruction for the tag with their respective `ReconstructionMethod` labels.

#### Anchor fallback behavior

If the tag's rows end before the requested `anchor_start`, the last census with a DBH becomes the anchor. If the requested anchor census has no row for the tag, or **all its rows have NA for both `DBH` and `TrueStemID`**, the algorithm searches backwards and uses the most recent earlier census that has at least one row with a non-NA `DBH` and a non-NA `TrueStemID` as the anchor instead of immediately falling back. If no such earlier census exists, a provisional anchor is used when allowed (see the note below); otherwise the tag goes to the probabilistic matcher (`anchor_missing_truestem`).

#### Anchor extension (forward search)

When the tag has measurements after the nominal anchor census and that census has 0 living stems (all DBH values are `NA`), the algorithm extends the anchor search **forward** to post-anchor censuses. It selects the first census after the anchor that has a living stem with a non-NA `TrueStemID` (or, when there is none, the first census with a living stem). This prevents tags where the anchor census happens to have only dead stems from unnecessarily falling back when a later census has reliable identity information.

Note: In addition, when the anchor census contains DBH observations but some lack a `TrueStemID`, `allow_provisional_anchor = TRUE` (the default of the DP function; driver flag `ALLOW_PROVISIONAL_DP_ANCHOR`, default `TRUE`) lets the DP give every measured row of that census a provisional `TrueStemID`/`ReconstructedStemID` (marked with `ReconstructionMethod = "provisional_dp"`) and proceed instead of falling back. The same is done at the last census with a DBH when no census has both a DBH and a `TrueStemID`.

### DP fallback reason codes (`DP_FallbackReason`)

To aid diagnostics the DP sets a `DP_FallbackReason` string on returned rows when the solver falls back to the probabilistic matcher (`NA` on tags it solves). Codes and suggested remedies:

| Code | Meaning | Suggested action |
|------|---------|------------------|
| `growth_form_forced` | A `growth_form` value of the tag is in `fallback_growth_forms`. | None — the routing was requested. |
| `species_forced_probabilistic` | A `Species` value of the tag is in `prob_species`. | None — the routing was requested. |
| `no_obs_up_to_anchor` | No DBH observations exist at or prior to the requested anchor. | Check `anchor_start` or ensure DBH observations exist before the anchor. |
| `anchor_missing_truestem` | Anchor census has no DBH and no `TrueStemID`, no earlier census has both, and no provisional anchor could be set. | Supply `TrueStemID` or enable `allow_provisional_anchor=TRUE`. |
| `anchor_missing_obs` | Anchor census has no DBH observations. | Ensure DBH at anchor or select a different `anchor_start`. |
| `anchor_missing_truestem_prov_disabled` | Anchor has DBH but `TrueStemID` missing and provisional anchors are disabled. | Set `allow_provisional_anchor=TRUE` or provide `TrueStemID`. |
| `anchor_ids_missing` | No anchor IDs were found after attempting to locate an anchor census. | Inspect data or supply anchor IDs. |
| `K_too_small` | Chosen `K` is smaller than the maximum per-census observed stems. | Increase `max_tracks` or `slack_tracks`. |
| `enum_exceeded` | State enumeration of a census returned no state matrix: more injective assignments than `max_states`. Routes to probabilistic matcher. | Increase `max_states` or reduce `K`. |
| `edge_count_exceeded` | Cross-product of adjacent census states reached `max_states²`. Routes to probabilistic matcher. | Increase `DP_MAX_STATES`. |
| `anchor_truestem_not_found` | Anchor `TrueStemID` values could not be mapped to track indices. | Check for unexpected `TrueStemID` values or duplicates. |
| `next_assign_row_mismatch` | Internal mismatch mapping assignment keys to enumerated rows. | Likely internal error — reproduce and report with a minimal example. |
| `no_reachable_next_states` | No reachable next full-states in the backward recursion (all pruned/filtered). | Relax pruning bounds (`prune_hard=FALSE`, widen `prune_*`) or check data. |
| `no_states_produced` | No DP states produced at a census after pruning/constraints. | Relax pruning or widen growth bounds. |
| `no_feasible_edges` | No feasible edges found between states after pruning/cost checks. | Relax pruning and growth bounds. |
| `decode_failure` | Failed to obtain a valid MAP start index. | Check pruning/parameter choices; report if persistent. |
| `viterbi_decode_failure` | Invalid pointer during Viterbi backtrace. | Likely a bug — reproduce and report. |
| `assign_mismatch` | MAP assignment length does not match observed rows. | Reproduce and report (enumeration bug). |
| `forward_edges_missing` | Missing edges during forward pass (unexpected). | Check edge construction and pruning. |
| `forward_no_alpha` | No finite forward alpha (numerical underflow or fully pruned transitions). | Relax pruning, adjust temperature, or report. |
| `error:<message>` | (drivers) The DP call raised an error; the tag was solved by the probabilistic matcher. | Read the message; reproduce with a single-tag run. |

> Note: When a tag with rows after the anchor falls back, its post-anchor rows are appended with the same `DP_FallbackReason`. If the output held several reasons, they are concatenated with `;` in the propagated column.

**How to check for fallback reasons (R example):**

```r
res <- match_stems_dp_global_backward_marginals_batch(dt, anchor_start = 5, ...)
unique(na.omit(res$DP_FallbackReason))
```

**Important:** `min_growth` and `max_growth` are passed to the probabilistic matcher, together with `prune_min_growth`, `prune_max_growth` and `prune_recruit_max_dbh`. When a prune bound is given, the matcher uses it in place of the corresponding base bound (`prune_recruit_max_dbh` in place of `Bio_Recruit_MaxDBH_unit`).

### Probabilistic Matcher: `match_stems_probabilistic()`

When the DP cannot be used for a tree (any reason in the fallback list above: a state space that is too large, a dead end with no feasible state, species or growth-form routing, …), `do_fallback()` routes it to `match_stems_probabilistic()`.

**Algorithm:**

1. **Pairwise log-likelihoods** (`compute_pairwise_log_likelihood()`): For each adjacent census pair, computes a log-likelihood matrix between all observed stems from the `Bio_*` parameters the DP uses: the Gaussian growth likelihood with size-dependent mean and variance, the log survival probability of the earlier stem, the soft shrinkage penalty, and a soft growth penalty above `Bio_Max_Growth`. The DP's measurement-error mixture is not implemented in the matcher. A link outside the growth bounds, or (with the `use_bio_hard_*_in_prob` flags) outside `Bio_Max_Shrink` / `Bio_Max_Growth`, is forbidden. DBH rounded down to classes at flagged censuses (`dbh_round_censuses`) is treated as in the DP (see *DBH recorded in classes*).

2. **Cost matrix augmentation** (`augment_cost_matrix()`): Expands the pairwise likelihood matrix to K×K with mortality slots (for stems disappearing) and recruitment slots (for new stems appearing), using the same mortality hazard and recruitment size/rate distributions as the DP. Two modes (`birth_death`, set from the DP's `prob_birth_death`, driver flag `PROB_BIRTH_DEATH`):
   - **Birth-death (default, `TRUE`).** K = n_curr + n_next. Each current stem has one death cell (`[i, n_next + i]`, log P(death)) and each next stem one recruitment cell (`[n_curr + j, j]`, log P(recruit) + log f(size)); the cells that pair a recruitment row with a death column carry no event and score 0. Every combination of survivals, deaths and recruitments is an assignment, scored by the sum of its events, so a link is taken only when it is more likely than the death of the earlier stem plus the recruitment of the later one — the choice the DP makes. The matrix carries attribute `"bd"`.
   - **Count-based slots (`FALSE`).** K = max(n_curr, n_next), raised only until an assignment without forbidden links exists (with M the largest set of allowed survival links, K ≥ n_curr + n_next − M), so every allowed assignment keeps M survivals: with stable stem counts no stem can die and none can be recruited. In a test against the 2010–2023 stem tags (identities hidden, 2023 pinned), this forced survival recovers 30% of the links of *Oenocarpus* clumps and 65% of their deaths, and before 2010 it gives 0.8–4.4 palm replacements (a death and a recruit in one tree and interval) per 100 trees per year, against 4.6–6.9 with certain identity. The birth-death mode recovers 45% of the links and 82–85% of the deaths, and gives 5.0–7.2 replacements per 100 palm trees per year in 1985–2010 (1982–85: 9.3, when palms with a single stem, whose deaths are certain, also died about three times faster than later). The mode only concerns tags solved by the matcher.

3. **Stochastic assignment** (`greedy_assignment_gumbel()`): Draws `n_samples` (default 200) stochastic assignments by adding Gumbel(0, temperature) noise to the log-likelihoods. Birth-death mode: the noise goes on every event cell (none on the empty cells) and each sample is the exact best assignment of the perturbed pair (`hungarian_min_rcpp()`, Kuhn–Munkres in `dp_global/src/transition_cost_rcpp.cpp`; one cell per event, so deaths and recruitments are not favoured by having several equivalent slots). Count-based slots: rows are assigned greedily to the best available column, and a pair whose greedy assignment uses a forbidden link is solved again exactly (`enforce_feasible_assignment()`). This approximates sampling from the Gibbs distribution over assignments.

   **TrueStemID pins while sampling.** Two masks constrain each pair's cost matrix. `apply_pin_mask()` sends an observation pinned to a stem present at the anchor to the column that carries that anchor stem. `apply_track_pin_mask()` covers the pins of every stem, including stems that end before the anchor (Step 3a.5 / Step 3b StemID pins, trees last measured before 2010): while sampling backward each next-census observation carries the pin of its track (its own pin, or one inherited from a later census of the same sample, `propagate_track_pin()`), and a pinned observation must join the observation that carries its pin and may not join one that carries a different pin. Unpinned observations are free, so a stem's earlier rows under an older, renumbered StemID can still be linked to it. The samples therefore group pinned observations as the pin sweep groups them in the exported table.

   **Pins never force a link (`pin_masked_pair()`).** The masks are applied to each sample's pair matrices by `pin_masked_pair()`. With count-based slots, step 2 sizes a pair for the growth limits only, so the masks can forbid links that its slots relied on: one stem at C6 pinned to A and one at C7 pinned to B form a 1×1 matrix (the growth between them is allowed, so there is no death slot) whose only link the mask forbids; without another slot every sample would have to join A to B. When the masked survival links (largest allowed set M) need more death/recruit slots than the pair has, `pin_masked_pair()` rebuilds the pair with K = n_curr + n_next − M (`augment_cost_matrix(K_min = …)`) and masks it again, so one pinned stem ends and the other starts. Pairs with enough slots keep their size (a birth-death matrix always has them); the verbose log counts the draws that needed extra slots ("Pins: … extra death/recruit slots").

4. **Sequential backward conditioning** (when `prob_lookahead_weight > 0`, the tag has at least three measured censuses and $K \geq 4$): When stitching per-pair assignments backward from the anchor, the cost matrix for census pair $(t, t+1)$ is conditioned on the already-resolved assignment at pair $(t+1, t+2)$. Specifically, `condition_cost_matrix()` computes a continuity bonus for each candidate column $j$ based on the log-likelihood that $j$'s forward assignment (the stem it maps to at $t+1$) continues plausibly to its assignment at $t+2$. This bonus is normalised (best column gets 0, others receive negative penalties capped at $-2$ log units) and weighted by `prob_lookahead_weight` (function default 0.5; the drivers pass `PROB_LOOKAHEAD_WEIGHT`, default 1). Conditioning is gated on $K \geq 4$ to avoid distorting results for simple tags with few stems.

5. **Backward stitching**: Starting from the anchor census (where identities are known), stitches per-pair assignments backward to build full per-sample track assignments across all censuses.

6. **Sample-level growth violation repair** (`repair_stitched_growth_violations()`): Before computing marginals, each sample's trajectories are walked and links violating growth constraints are severed. This ensures that `compute_marginals_from_samples()` only counts biologically valid paths, so posterior probabilities reflect post-constraint uncertainty. Two layers of defense:

   - **Hard-rate check:** Annualised growth outside `[min_rate, max_rate]` is severed immediately. Pinned observations are never severed: if the earlier observation of the link is pinned, the later one is severed instead (unless it is pinned too or sits at the anchor), and a link between two pinned observations is kept, as the pin sweep keeps it in the export.
   - **ME-informed cumulative-shrinkage check** (active when `use_bio_hard_shrink_in_prob = TRUE` and `n_sigma_me` is finite): Tracks cumulative shrinkage along each trajectory. Even when each consecutive pair passes the hard rate, a long run of small decreases can accumulate more shrinkage than measurement error can explain. The threshold is derived from the small-error component of the BCI measurement-error model: $\text{SD}(D) = 0.0062 \times D + 0.0904$ cm. When cumulative shrinkage exceeds $n_\sigma \times \sqrt{\text{SD}(d_\text{start})^2 + \text{SD}(d_\text{curr})^2}$ ($n_\sigma$ = `n_sigma_me`, default 3), the link at the start of the shrinkage run is severed. The BCI driver sets `PROB_N_SIGMA_ME = Inf`, which turns this check off; `USE_BIO_HARD_SHRINK_IN_PROB=FALSE` also disables it and allows confirmed large-shrinkage events without forced splitting.

   **Relationship to DP:** The exact DP solver has no explicit trajectory repair: its global cost minimisation sums each per-step growth likelihood over the entire trajectory, so a run of decreases pays a compounding cost even when each individual step is within hard bounds. The probabilistic matcher works per census pair and lacks this global cost accumulation. The ME cumulative-shrinkage check provides an analogous trajectory-level constraint for it.

7. **Pin-consistent sample filter** (`filter_pin_consistent_samples()`): samples in which an observation pinned to a stem present at the anchor is not on that stem's track are dropped, so the marginals are conditioned on those pins. If fewer than `min(10, n_samples / 4)` samples would remain, all are kept (soft pins) and a warning is logged. This happens when no sample can satisfy a pin, e.g. a stem measured at C3 and again at the anchor but not in between, while other stems are measured in between: the matcher has no alive-but-unmeasured state, so that stem's track cannot bridge the gap (the pin sweep, step 11, joins the two rows in the export).

8. **Marginal posterior computation** (`compute_marginals_from_samples()`): Two-pass approach:
   - **Pass 1**: Computes the full marginal posterior distribution for every observation by counting track assignments across samples (Top-K IDs and probabilities, entropy).
   - **Pass 2**: Growth-aware greedy conflict resolution — censuses are processed from the anchor outward (anchor first, then ±1, ±2, …). Within each census, observations are sorted by confidence (descending); each observation gets its MAP track ID unless it conflicts with an already-assigned ID **or** it would create a growth-rate violation against the nearest already-resolved adjacent census. When a candidate ID is rejected (by conflict or growth check), the next-best posterior alternative is tried. `match_stems_probabilistic()` passes `intervals`, `min_rate` and `max_rate`, which turns the growth check on; when one of them is `NULL`, the censuses are resolved in order and only ID conflicts are checked.

   The IDs resolved in Pass 2 are replaced by the consensus export (step 9); the Top-K and entropy columns stay marginal.

9. **Consensus export** (`select_consensus_trajectory()`): `ReconstructedStemID` is the most representative sample: the one whose predecessor links agree most with the other samples' links (maximum expected accuracy; ties are broken by how often its whole partition occurs). The exported reconstruction is therefore a trajectory that exists in a sample. `DP_PosteriorReconstructedProb` is the share of samples that give the observation the exported track.

10. **Diagnostic growth-violation check**: A scan walks each exported stem trajectory and counts residual hard-rate violations. The check is diagnostic only — it does not modify `tree_data` or posteriors. Any violations are logged as warnings.

11. **Pin sweep** (when `pin_truestemid = TRUE`, the default): every row with a database pin (non-NA `TrueStemID`, not a provisional anchor) gets `ReconstructedStemID = TrueStemID` and `ReconstructionMethod = "given"`. The engine's own value is kept in `ReconstructedStemID_PreSweep`, and `SweepAuditOverride = TRUE` flags the rows the sweep changed. Unpinned rows keep the engine's ID; after the engine, `apply_pin_track_rejoin()` moves unpinned rows that the sweep left behind on a track whose pinned rows it moved (one pin per track segment).

12. **Posterior export** (`export_probabilistic_posteriors()`): Stages the per-tag samples in the same format as the DP posterior sampler; `finalize_posterior_paths()` writes the final path file after the post-engine helpers.

**Parameters:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `n_samples` | `200` | Number of Gumbel-noise stochastic samples (DP argument `prob_n_samples`) |
| `temperature` | `1.0` | Temperature for Gumbel sampling (shared with DP) |
| `posterior_top_k` | `2` | Number of top posterior alternatives to report |
| `prob_lookahead_weight` | `0.5` | Weight for sequential backward conditioning (0 = disabled) |
| `n_sigma_me` | `3` | Number of ME standard deviations for cumulative shrinkage threshold (`Inf` = off; DP argument `prob_n_sigma_me`). The small-error SD model of this check, SD(D) = 0.0062 × D + 0.0904 cm, is fixed in the call to `repair_stitched_growth_violations()` |
| `birth_death` | `TRUE` | Birth-death slots (step 2); `FALSE` = count-based slots (DP argument `prob_birth_death`) |
| `pin_truestemid` | `TRUE` | Honour `TrueStemID` pins while sampling and in the pin sweep |
| `dbh_round_censuses`, `dbh_round_max`, `dbh_round_width` | `integer(0)`, `5.5`, `0.5` | DBH recorded in classes (see *DBH recorded in classes*) |
| `use_bio_hard_shrink_in_prob` | `TRUE` | Apply `Bio_Max_Shrink` hard gate in pairwise edge construction and ME cumulative-shrinkage check in trajectory repair. When `FALSE`, edges with growth below `Bio_Max_Shrink` are allowed (penalised by soft `k_shrink` only) and the Layer 2 cumulative-shrinkage check is skipped. Propagated from `USE_BIO_HARD_SHRINK_IN_PROB`. |
| `use_bio_hard_growth_in_prob` | `TRUE` | Apply `Bio_Max_Growth` hard gate in pairwise edge construction. When `FALSE`, edges exceeding `Bio_Max_Growth` are allowed (penalised by soft `k_growth` only). Propagated from `USE_BIO_HARD_GROWTH_IN_PROB`. |

**Output columns:** The probabilistic matcher produces the same output schema as the DP solver (including `DP_PosteriorTop1ID`, `DP_PosteriorTop1Prob`, `DP_PosteriorEntropy`, etc.; `DP_PosteriorUnlinkedProb` stays `NA`). `ReconstructionMethod` is `"probabilistic"` for unpinned rows, `"given"` for rows with a database pin (pin sweep, step 11) and `"provisional_dp"` for provisional anchor rows. Anchor-census rows have `DP_PosteriorReconstructedProb = 1` (every sample gives them their anchor ID).

### DP vs Probabilistic: How Shrinkage is Handled

Both engines read the same `Bio_*` parameters and receive the same growth and prune bounds, but they apply them differently:

| Aspect | DP (global) | Probabilistic (per census pair) |
|--------|-------------|----------------------|
| **Scope** | Minimises total cost over the entire trajectory | Stochastic assignment per census pair, stitched backward from the anchor |
| **Shrinkage defence** | Global cost accumulation: each step's growth likelihood compounds across the trajectory, penalising runs of decreases | Explicit trajectory repair: ME cumulative-shrinkage check severs runs exceeding the measurement-error threshold (when `use_bio_hard_shrink_in_prob = TRUE` and `n_sigma_me` is finite) |
| **Hard bounds** | Transitions outside the prune bounds are removed before costing (`prune_hard = TRUE`); a transition outside `Bio_Max_Shrink` / `Bio_Max_Growth` that is not pruned costs `hard_penalty` ($10^6$). Unaffected by `USE_BIO_HARD_SHRINK_IN_PROB` | Links outside the prune bounds (or `min_growth` / `max_growth`) are forbidden, and so are links outside `Bio_Max_Shrink` / `Bio_Max_Growth` when `use_bio_hard_shrink_in_prob` / `use_bio_hard_growth_in_prob` are `TRUE`; the hard-rate check of the sample-level repair uses the same rate bounds |
| **Soft penalty** | `k_shrink` on every decrease; `k_growth` above the soft cap `Bio_Max_Growth_Soft` (both 0 in the drivers' default configuration) | `k_shrink` on every decrease; `k_growth` above `Bio_Max_Growth`, so it acts only when the bio growth gate is off |
| **Survival term** | Mortality cost only on the track that dies | Log survival probability added to every link |
| **Measurement error** | 4-component mixture model when `use_measurement_error = TRUE` | No mixture; the ME coefficients 0.0062 and 0.0904 serve as a ruler for the cumulative threshold regardless of `USE_MEASUREMENT_ERROR` |
| **Asymmetry** | Bio bounds always act (as `hard_penalty`) | With `use_bio_hard_shrink_in_prob = FALSE`, probabilistic is **more permissive** than the DP for the same tree; use with care |


---

## ForestGEO Code Handling

The engines read two columns that encode field-recorded events: the code list `ListOfTSM` and `Status`. They influence transition feasibility, segment splitting and the handling of unmeasured censuses.

### Resprout Codes (R, RP, RF, RT, QR, OR)

**Recognized codes:** `R`, `RP`, `RF`, `RT`, `QR` and `OR`, matched as whole codes in `ListOfTSM`. A row with `Status == "broken below"` is treated the same way. A **resprout observation** is such a row *with a DBH*.

When the tag has a resprout observation at a census before the anchor, the mechanisms below apply.

#### 1. Resprout segment split

The DP splits the problem into two independent sub-problems at the first census (not the first census of the tag) that holds a resprout observation and lies before the anchor:

- **Pre-segment**: from the earliest census up to and including the R census, anchored at its last measured census (with provisional anchor IDs where it has no `TrueStemID`). In the DP the R-coded measurement is the last record of the stem that broke.
- **Post-segment**: from the census after the R census to the anchor. All stems at its first census are treated as new stems (recruits) that are independent of the stems in earlier censuses.

The split is biologically motivated: R means the tree resprouted, so the stems measured after that census have no identity continuity with the stems before the break. The exported identity of the R-coded row itself is settled after the engine: the pin sweep puts pinned rows on their database identity, and `apply_broken_below_invariants()` (R1) starts a new ID at a measured broken-below row that still has earlier rows on its trajectory.

Pre-segment IDs are offset by the maximum post-segment ID to prevent clashes. **Posterior samples:** each sub-call draws its own samples and hands them back (`posterior_return_samples = TRUE`, or `return_samples` in the probabilistic matcher) instead of writing a staging file, which the other segment's file would overwrite. Because every stem ends at the split, the two segments are independent and draw *k* of one segment paired with draw *k* of the other is an exact draw of the whole tag (`pair_segment_posterior_samples()`); the pre segment uses `posterior_sample_seed + 1` so the paired draws are independent, a segment with fewer draws is recycled, a segment that drew none (e.g. one census) contributes its exported grouping to every draw, and pre-segment labels are shifted above all post-segment labels. The paired samples are staged once per tag (`stage_posterior_samples()`). Sub-calls set `allow_segment_split = FALSE` to prevent cascading splits; a later resprout observation inside a segment is handled by the R-recruit constraint (3). Virtual census IDs (gaps in `seq.int()` not present in the sub-call's data) are filtered out.

#### 2. Post-segment recruit continuity constraint

In the post-segment, all stems at the first census (the census after the R census) are known recruits. When the number of stems at the first census ≤ the number at the next census, the constraint requires that every first-census stem be on a track that is also occupied at the next census. The rationale: a freshly resprouted stem dying immediately while a new independent recruit simultaneously appears on a different track is far less parsimonious than all recruits continuing.

This constraint prunes candidate assignments where any recruit's track is empty at the next census.

**Dead-end safety guard:** If applying the recruit continuity constraint would remove *all* candidate transitions at a census (creating a dead end in the backward recursion), the constraint is skipped for that census with a diagnostic message. This prevents the DP from becoming infeasible due to overly restrictive pruning in edge cases where no transition can simultaneously satisfy recruit continuity and the biological cost model.

#### 3. R-recruit constraint

A resprout observation at census $t+1$ must **continue a track that is occupied at census $t$**: in the DP it is the last record of the stem that broke, so it cannot sit on an empty (not yet recruited) track. It therefore needs no track of its own, and K is not raised for resprout observations.

The constraint is enforced for every candidate transition in `derive_phase_prev_batch_rcpp()` (C++): a transition that places a resprout observation on a track that is empty at the prior census is infeasible. The backward loop repeats the test as an R-level filter, which is skipped when census $t$ has no measured stem; the C++ test has no such exception, so in that case no transition is feasible and the tag (or segment) goes to the probabilistic matcher (`no_states_produced`).

#### 4. NA-R barrier

A census before the anchor with **no measured stem** and at least one unmeasured row that carries a resprout code or `Status == "broken below"` is a barrier: the whole tree died back, so no identity crosses it. After decoding, a track with rows on both sides of the barrier is cut and its unpinned rows before the barrier take a new ID. A track is left whole when its ID is a `TrueStemID` known after the barrier and no row before the barrier carries a different `TrueStemID`; pinned rows are never relabelled. The unmeasured R-coded rows of the barrier census take the IDs of the stems of the preceding census, in row order. `do_fallback()` applies the same cut to the probabilistic matcher's output (every crossing track; the pin sweep then restores the pinned rows).

### Missing-from-Field (MF) Detection

Rows of stems that were alive but not measured are set aside before the engine runs and re-inserted afterwards (`ReconstructionMethod = "dp_mf_inferred"`), so the gap is not read as a death followed by a recruitment. Two cases are detected:

- **Explicit:** a row without a DBH whose `ListOfTSM` holds the code `MF` or whose `Status` is `"missing"`. The episode extends through the following consecutive censuses in which no row of the tag has a DBH.
- **Implicit:** a census in which no row of the tag has a DBH, lying between the tag's first and last measured census.

**R-code exclusion:** A census whose rows carry a resprout code or `Status == "broken below"` is not an implicit MF census. Such a record marks a break, not a field absence.

A re-inserted row takes an ID only when exactly one track is assigned in the nearest census before and after and has no row at the MF census; otherwise its `ReconstructedStemID` stays `NA`.

### Regex Implementation

All code matching uses `grepl()` with `perl = TRUE` to ensure `\b` word-boundary anchors work correctly. The resprout regex pattern is:

```r
"\\b(R|RP|RF|RT|QR|OR)\\b"
```

This matches whole-word `R`, `RP`, `RF`, `RT`, `QR`, or `OR` codes in fields that may contain multiple semicolon-separated codes (e.g., `"R;B"` matches, `"NORMAL"` does not). The MF pattern is `"\\bMF\\b"`. The labelling of post-anchor rows uses the pattern without `QR`.


## Pruning & Conservative Guards

### Motivation — factorial growth of the state space

The DP enumerates injective assignment states per census. For a census with `n_obs` observed stems and `K` tracks the number of assignment states is the falling-permutation

```
P(K, n_obs) = K × (K-1) × ... × (K - n_obs + 1)
```

The total number of full-paths (and candidate transitions between adjacent census assignment states) grows multiplicatively across censuses, which quickly becomes intractable (factorial-like explosion). To keep the DP practical, *conservative, cheap* pre-filters remove biologically impossible or extremely implausible candidate transitions before the full (expensive) transition cost is evaluated.

These filters are intentionally conservative: they are meant to remove clearly impossible candidates (hard physical/biological limits) while preserving borderline cases for the full cost evaluation.

### Where pruning runs

- Pruning happens in the backward recursion, in `derive_phase_prev_batch_rcpp()` (C++), which tests every pair of a state at census $t$ and a reachable state at $t+1$: first the life-cycle phase rules and the R-recruit constraint, then the hard growth and recruit-size bounds. Only the pairs that pass are costed, in one call to `transition_cost_paired_rcpp()`. This avoids the costly likelihood evaluation for obviously impossible transitions.
- The growth and recruit-size bounds are only applied when `prune_hard = TRUE` (default). With `prune_hard = FALSE` the phase rules and the R-recruit constraint still apply, and every other candidate transition is passed to the cost routine (slower).
- The prune bounds also restrict the state enumeration (`enumerate_states_constrained()`): a track that holds an anchor stem is open to an observation only if an observation of the next census that may use that track is growth-compatible with it (for non-taper growth forms the test uses twice the non-taper bounds, to leave room for the HOM widening). Tracks that are empty at the anchor are open to every observation, and a pinned observation is restricted to its pinned track. A census can therefore have far fewer states than $P(K, n_{obs})$.


### Parameters & exact semantics

We separate pruning thresholds from the biological (`Bio_*`) parameters so you can control pruning behaviour without changing the biological model used in the transition-cost calculations.

- `prune_min_growth` (numeric | NULL): explicit lower bound (cm/year) used for pruning. If `NULL` the value `min_growth` passed to the DP is used.
- `prune_max_growth` (numeric | NULL): explicit upper bound (cm/year) used for pruning. If `NULL` the value `max_growth` passed to the DP is used.
- `prune_use_bio_bounds` (logical, default TRUE): whether to intersect the user-specified prune bounds with the biological hard bounds found in the data (`Bio_Max_Shrink`, `Bio_Max_Growth`). If TRUE (default):
  - eff_min_grow = max(user_min, Bio_Max_Shrink)
  - eff_max_grow = min(user_max, Bio_Max_Growth)
  where `user_min` = `prune_min_growth` if provided else `min_growth`, and similarly for `user_max`.
  If `prune_use_bio_bounds = FALSE` the `eff_*` values equal `user_*` directly.
- `prune_recruit_max_dbh` (numeric | NULL): override for the recruit-size cut used during pruning. If `NULL` the biological `Bio_Recruit_MaxDBH_unit` is used. If `prune_use_bio_recruit = TRUE` (default) and `Bio_Recruit_MaxDBH_unit` is finite, the more conservative value `min(prune_recruit_max_dbh, Bio_Recruit_MaxDBH_unit)` is used; otherwise the explicit override is used. A recruit larger than this cut is pruned; a recruit above `Bio_Recruit_MaxDBH_unit` that is not pruned costs `hard_penalty`.
- `prune_use_bio_recruit` (logical, default TRUE): controls whether the biological recruit bound is intersected with `prune_recruit_max_dbh` (the minimum of the two is used). `Bio_Recruit_MaxDBH_unit` is either the 0.999 quantile of the observed recruit DBHs or a fixed value, as chosen when the parameters were estimated (`recruit_max_source`).

#### Non-taper-corrected growth form override

Growth forms that are **not taper-corrected** (default: `"palm"`, `"strangler_fig"`, `"tree_fern"`) get their own pruning bounds, for two reasons:

1. **Real biological growth.** Palms, strangler figs, and tree ferns grow in DBH — palms can add 1–3 cm/yr, and strangler figs can change diameter substantially as they encircle or replace their host.
2. **Apparent DBH variation from HOM changes.** When the measurement height (HOM) shifts between censuses, the recorded DBH can change dramatically even if the true diameter at a fixed height remained constant. Because these growth forms lack the tapered trunk geometry that allows taper correction, any HOM shift produces an uncompensated apparent DBH change.

For these reasons, the non-taper-corrected override **replaces** the general effective pruning bounds with its own pair of bounds, which are meant to be at least as wide as the standard growth/shrink limits. The function defaults are −0.625 and 6.25 cm/year; the drivers pass `PRUNE_BOUND_FACTOR ×` the standard limits (`MAX_SHRINK_FIXED` and `MAX_GROWTH_FIXED`). An optional HOM-proportional widening layer then extends the bounds further on a per-census-pair basis when HOM data are available.

**Parameters:**

- `non_taper_corrected_growth_forms` (character vector, default `c("palm", "strangler_fig", "tree_fern")`): growth forms whose DBH measurements are not taper-corrected. When a tag's `growth_form` matches any entry in this list (exact match; a tag with several `growth_form` values uses the most common one, with a warning), the non-taper override is activated. Accepts comma- or semicolon-separated strings. The values must be the `growth_form` labels of the data: the BCI driver assigns `palm`, `strangler` and `fern` and therefore sets `c("palm", "strangler", "fern")`; it stops with a visible check when a value listed here, in `DP_FALLBACK_GROWTH_FORMS` or in `PROB_SPECIES` does not occur in the data.
- `non_taper_corrected_prune_min_growth` (numeric, default `-0.625`): lower prune bound (cm/year) that replaces the general effective minimum for matching growth forms. The drivers pass `PRUNE_BOUND_FACTOR × MAX_SHRINK_FIXED`.
- `non_taper_corrected_prune_max_growth` (numeric, default `6.25`): upper prune bound (cm/year) that replaces the general effective maximum for matching growth forms. The drivers pass `PRUNE_BOUND_FACTOR × MAX_GROWTH_FIXED`.
- `hom_tolerance_scale` (numeric, default `2.0`): additional DBH tolerance (cm) per meter of HOM deviation from 1.3 m, spread over the census interval. When a `hom` (or `HOM`) column (in m) is present in the data and the tag is non-taper-corrected, the prune bounds are widened for each census pair by:

  ```
  hom_tol = hom_tolerance_scale × max(|HOM − 1.3|, na.rm=TRUE) / interval_years
  eff_min_grow_pair = non_taper_corrected_prune_min_growth − hom_tol
  eff_max_grow_pair = non_taper_corrected_prune_max_growth + hom_tol
  ```

  where `max(|HOM − 1.3|)` is the worst-case deviation across all rows of the tag at the two censuses being compared. NA HOM values are treated as 1.3 m (zero deviation contribution). Set `hom_tolerance_scale = 0` to disable HOM widening.

**How the override layers interact:**

1. **User layer**: general prune bounds are set from `prune_min/max_growth` (or `min/max_growth` if NULL).
2. **Bio layer**: if `prune_use_bio_bounds = TRUE`, general bounds are tightened by intersecting with biological hard limits (`Bio_Max_Shrink`, `Bio_Max_Growth`).
3. **Non-taper override**: if the tag's `growth_form` matches `non_taper_corrected_growth_forms`, the effective bounds from steps 1–2 are **replaced** with `non_taper_corrected_prune_min/max_growth`. This ensures non-taper forms always get wide bounds regardless of how tight the bio-constrained general bounds may be.
4. **HOM widening**: if a `HOM` column is present and `hom_tolerance_scale > 0`, the bounds from step 3 are **widened symmetrically** for each census pair based on the worst-case HOM deviation.

**Note on the default driver configuration:** The driver scripts pass `prune_min_growth = PRUNE_BOUND_FACTOR × MAX_SHRINK_FIXED` and `prune_max_growth = PRUNE_BOUND_FACTOR × MAX_GROWTH_FIXED` with `prune_use_bio_bounds = FALSE`, and the same two values as the non-taper bounds. The general bounds therefore already equal the non-taper bounds, so in the default configuration the non-taper override at step 3 is a no-op — all trees (taper-corrected or not) receive the same base prune window: `[-0.625, 9.375]` cm/year in `main_cpp.R` (factor 1.25, `MAX_GROWTH_FIXED = 7.5`), `[-0.625, 6.25]` in `main_cpp_chunk.R` (factor 1.25, `MAX_GROWTH_FIXED = 5`) and `[-2.5, 25]` in the BCI drivers (factor 5). The effective differentiation for non-taper forms comes from step 4 (HOM widening), and the non-taper override provides an independent lever when users tighten the general bounds (e.g. via `prune_use_bio_bounds = TRUE` or narrower `prune_min/max_growth`).

Notes:

- These pruning parameters affect only the *pre-filtering* step (they do not change the transition cost function or post-hoc diagnostics beyond recording the effective prune values).
- If the interval between two censuses is not finite or not positive, the growth and recruit-size pruning are both skipped for that pair.
- You can always define extra margins for these parameters. For example:

```r
prune_min_growth = MAX_SHRINK_FIXED * 2.5 # very wide fixed bounds
prune_max_growth = MAX_GROWTH_FIXED * 1.5 # very wide fixed bounds
prune_use_bio_bounds = FALSE # use fixed prune bounds instead of biological ones
prune_recruit_max_dbh = RECRUIT_MAX_FIXED * 1.2 # very high recruit max dbh
prune_use_bio_recruit = FALSE # use prune_recruit_max_dbh as given, not min(prune_recruit_max_dbh, Bio_Recruit_MaxDBH_unit)
```

### Diagnostics & reproducibility

- The DP exposes `attr(out, "DP_PruneInfo")` with:
  - `total_examined`, `total_pruned` (counts of candidate transitions; `total_pruned` counts every transition removed before costing: phase-infeasible ones, pruned ones and those removed by the R-code constraints)
  - `per_census` (removed transitions per census pair, named by the earlier census)
  - `eff_min_growth`, `eff_max_growth`, `eff_recruit_max` (the effective thresholds used, before HOM widening)
  - `is_non_taper_corrected`, `use_hom_relax` (whether the non-taper override and the HOM widening were active)
- The counts are only filled when `prune_hard = TRUE` and the tag is solved by the DP.
- The effective prune bounds are also `vcat`-logged at the start of the backward pass for transparency.

### Practical guidance and examples

- Function defaults: with `prune_min_growth = NULL` and `prune_max_growth = NULL` the code uses `user_min = min_growth`, `user_max = max_growth` and intersects them with the `Bio_*` values. This preserves biological hard limits while letting you set study-level `min_growth`/`max_growth` to narrow behaviour across the run. (The driver scripts pass explicit prune bounds, see the note above.)
- To apply *wider* pruning (less aggressive rejection) than the biological bounds allow (for example, when the DP cost model is trusted to handle extreme cases), set `prune_use_bio_bounds = FALSE` and set `prune_min_growth`/`prune_max_growth` to the desired wide range (e.g., `-10`..`25`).
- To apply *stricter* recruit-size pruning, provide a small `prune_recruit_max_dbh` (and set `prune_use_bio_recruit = FALSE` to use it without intersecting the biological bound).
- Beware of overly tight pruning: if pruning removes all feasible transitions at a census the DP will fall back to the probabilistic matcher (or produce no DP states). If you see frequent fallbacks for reasonable data, relax the prune bounds or set `prune_hard = FALSE` for a diagnostic run.

### Why this approach?

- The combinatorial nature of injective assignments makes per-census state counts grow factorially with `n_obs`. Pruning cheaply removes impossibilities and reduces the number of pairwise assignment evaluations from `O(n_states_cc × n_states_{cc+1})` to a manageable number while retaining feasible candidates for the full probabilistic scoring.
- Because pruning is conservative and logged, it is auditable: you can use `DP_PruneInfo` to quantify how many candidate transitions were removed and where.

---

## Workflows & Usage Patterns

### Single-Tag Debug

```r
# In dp_global/scripts/main_cpp.R (or via CLI: --RUN_ALL_TAGS=FALSE --WHICH_TAG=19):
RUN_ALL_TAGS <- FALSE
WHICH_TAG <- "19"

WRITE_DP_PDF <- TRUE
```

### Single-Tag with Posterior Uncertainty

```r
DP_MODE <- "marginals+bins"  # default; also adds posterior bins
POSTERIOR_SAMPLES <- 200L    # default; 0 disables the sampled path files
# the drivers pass temperature = 1 to match_stems_dp_global_backward_marginals_batch()
```

### Full Parallel Run

```r
RUN_ALL_TAGS <- TRUE
MANUAL_CORES <- TRUE
MANUAL_CORES_VALUE <- 8L

# Adjust if needed:
DP_MAX_STATES <- 1e5
DP_MAX_TRACKS <- 50L
```

### Sensitivity Analysis Only

```r
source("dp_global/R/dp_global_main.R")
source("dp_global/R/sensitivity_transition_cost_bio.R")
# build_all_sweeps() / plot_all_sweeps_to_pdf(); or run main_cpp.R with --SENSITIVITY_MODE=run+write+pdf
```

### Example: enforce growth/recruit bounds

Units: growth bounds are **cm/year**; recruit caps are **cm**.

```r
# Enforce growth bounds (0.05 - 0.75 cm/yr) and cap recruits at 37.5 cm
bio <- estimate_bio_pars(
  x,
  anchor_start_census = 7L,
  enforce_growth_bounds = TRUE,
  growth_min_fixed = 0.05,
  growth_max_fixed = 0.75,
  enforce_recruit_max = TRUE,
  recruit_max_source = "fixed",
  recruit_max_fixed = 37.5
)
```

### Realism Report Only

```r
source("dp_global/R/realism_calibration.R")
# realism_report_from_reconstruction(out, interval_years, base_args); or run main_cpp.R with --RUN_REALISM_REPORT=TRUE
```

### Parameter Estimation Workflow

The example below runs on the simulated data (`data_simulation/data/simulated_data_1.csv`) from the project root. The driver scripts do the same and then apply the post-engine helper chain (see *Stem Identity Renumbering Workflow*).

```r
library(data.table)
library(here)
source(here("dp_global", "R", "dp_global_main.R"))

xraw <- fread(here("data_simulation", "data", "simulated_data_1.csv"))
xraw[, species := as.character(Species)]
USE_MEASUREMENT_ERROR <- TRUE

# 1. Estimate parameters per species (intervals come from ExactDate)
bio_pars_list <- list()
for (sp in unique(xraw$species)) {
  bio_pars_list[[sp]] <- estimate_bio_pars(
    xraw[species == sp],
    anchor_start_census = 7L,    # first census used (the default reads the global ANCHOR_START_CENSUS)
    use_measurement_error = USE_MEASUREMENT_ERROR,
    # Hard constraint sources
    max_shrink_source = "data",  # or "fixed"
    max_shrink_fixed = -0.5,
    max_growth_source = "data",  # or "fixed"
    max_growth_fixed = 7.5,
    # Soft penalty sources
    k_shrink_source = "data",    # or "fixed"
    k_shrink_fixed = 50,
    k_growth_source = "data",    # or "fixed"
    k_growth_fixed = 50,
    # Recruit size cap
    recruit_max_source = "data", # or "fixed"
    recruit_max_fixed = 5,
    # Quantile configuration (optional)
    shrink_data_quantile = 0.001,
    shrink_hard_prob = 1e-4,
    growth_data_quantile = 0.999,
    growth_hard_prob = 1e-4,
    growth_soft_quantile = 0.99,
    recruit_max_quantile = 0.999
  )
}

# 2. Attach as columns to tree_data
xraw[, Bio_Mu_Growth := bio_pars_list[[species]]$growth$alpha, by = species]
xraw[, Bio_Gamma_Growth := bio_pars_list[[species]]$growth$gamma, by = species]
xraw[, Bio_Sigma0_Growth := bio_pars_list[[species]]$growth$sigma0, by = species]
xraw[, Bio_Sigma1_Growth := bio_pars_list[[species]]$growth$sigma1, by = species]
xraw[, Bio_Max_Shrink := bio_pars_list[[species]]$shrinkage$max_shrink, by = species]
xraw[, Bio_K_Shrink := bio_pars_list[[species]]$shrinkage$k_shrink, by = species]
xraw[, Bio_Max_Growth := bio_pars_list[[species]]$growth$max_growth, by = species]
xraw[, Bio_Max_Growth_Soft := bio_pars_list[[species]]$growth$max_growth_soft, by = species]
xraw[, Bio_K_Growth := bio_pars_list[[species]]$growth$k_growth, by = species]
xraw[, Bio_H0_Mortality := bio_pars_list[[species]]$mortality$h0, by = species]
xraw[, Bio_Beta_Mortality := bio_pars_list[[species]]$mortality$beta, by = species]
xraw[, Bio_Recruit_Meanlog := bio_pars_list[[species]]$recruitment$meanlog, by = species]
xraw[, Bio_Recruit_Sdlog := bio_pars_list[[species]]$recruitment$sdlog, by = species]
xraw[, Bio_Recruit_MaxDBH_unit := bio_pars_list[[species]]$recruitment$recruit_max_dbh, by = species]
xraw[, Bio_Recruitment_lambda := bio_pars_list[[species]]$recruitment$lambda, by = species]

# 3. Run the DP solver, one tag at a time (marginals are always computed)
out <- rbindlist(lapply(unique(xraw$Tag), function(tg) {
  match_stems_dp_global_backward_marginals_batch(
    tree_data = copy(xraw[Tag == tg]),
    min_growth = -2,           # growth bounds: pruning (when prune_* are NULL), matcher, ConstraintViolation
    max_growth = 10,
    anchor_start = 7,
    max_tracks = 30,
    max_states = 50000,
    slack_tracks = 1,
    temperature = 1,
    posterior_top_k = 2,
    use_measurement_error = USE_MEASUREMENT_ERROR,
    verbose = TRUE
  )
}), fill = TRUE)

# 4. Add posterior bins
out <- add_dp_posterior_bins(
  out,
  confident_prob = 0.95,
  unlinked_prob = 0.50,
  use_reconstructed_prob = TRUE
)
```

In the driver scripts the recruit cap comes from `RECRUIT_MAX_SOURCE` / `RECRUIT_MAX_FIXED` (default `"fixed"`, `MAX_GROWTH_FIXED × 5 + 0.9999` cm), the guardrails from `MAX_SHRINK_FIXED` / `MAX_GROWTH_FIXED` (source `"fixed"`), and both soft penalties are fixed at 0.

**Important workflow notes:**

1. **min_growth/max_growth in DP calls:** These parameters are **not part of the transition cost**. They are used as:
   - The prune bounds of the DP when `prune_min_growth` / `prune_max_growth` are `NULL` (intersected with `Bio_Max_Shrink` / `Bio_Max_Growth` when `prune_use_bio_bounds = TRUE`)
   - The growth bounds of the probabilistic matcher (replaced by the prune bounds when those are given)
   - The bounds of the post-hoc `ConstraintViolation` diagnostic

2. **Constraints inside the DP cost:** Come from Bio_* columns:
   - Hard shrinkage: `Bio_Max_Shrink`
   - Hard growth: `Bio_Max_Growth`
   - Soft shrinkage: `Bio_K_Shrink`
   - Soft growth: `Bio_K_Growth` and `Bio_Max_Growth_Soft`

3. **Measurement error:** Must be set consistently between `estimate_bio_pars()` and DP solver calls.

4. **Parameter sources:** Each parameter's source (data/fixed) and estimated value are stored in the returned list structure for provenance tracking.

5. **Bio_* columns:** Each required `Bio_*` column (see *Biological Parameter Columns*) must hold one non-NA value per tag; the DP stops with an error otherwise.

---


## Implementation Reference

### Function Map

| Task | Primary Function | Location |
|------|------------------|----------|
| State enumeration | `enumerate_states_injective()`, `enumerate_states_constrained()` | `dp_global/R/dp_global_states.R` |
| Track DBH vector of one state (helper; the DP builds its track-DBH matrices directly) | `state_to_track_dbh()` | `dp_global/R/dp_global_states.R` |
| Transition costs (C++) | `transition_cost_paired_rcpp()` (called by the DP), `transition_cost_tracks_bio_batch_rcpp()` | `dp_global/src/transition_cost_rcpp.cpp`, R wrappers in `dp_global/src/transition_cost_rcpp.R` |
| Phase feasibility and hard pruning (C++) | `derive_phase_prev_batch_rcpp()` | `dp_global/src/transition_cost_rcpp.cpp` |
| Cost breakdown (debug) | `transition_cost_tracks_bio_components()` | `dp_global/R/dp_global_bio.R` |
| Posterior marginals (production) | `match_stems_dp_global_backward_marginals_batch()` | `dp_global/R/dp_global_dp.R` |
| Posterior binning | `add_dp_posterior_bins()` | `dp_global/R/dp_global_diag.R` |
| Probabilistic matcher | `match_stems_probabilistic()` | `dp_global/R/dp_probabilistic_matching.R` |
| R-boundary splitting (fallback paths) | `split_live_resprout_tracks()` | `dp_global/R/dp_global_dp.R` |
| Posterior samples of a resprout-split tag (paired segments) | `pair_segment_posterior_samples()`, `segment_fixed_assignment()`, `stage_posterior_samples()` | `dp_global/R/dp_global_dp.R` |
| Pins of every stem while sampling (probabilistic) | `apply_track_pin_mask()`, `propagate_track_pin()` | `dp_global/R/dp_probabilistic_matching.R` |
| Extra death/recruit slots when pins forbid a pair's links (probabilistic, count-based slots) | `pin_masked_pair()`, `augment_cost_matrix(K_min = …)` | `dp_global/R/dp_probabilistic_matching.R` |
| Birth-death assignment: one death and one recruitment cell per stem, exact perturbed assignment (probabilistic) | `augment_cost_matrix(birth_death = TRUE)`, `greedy_assignment_gumbel()`, `hungarian_min_rcpp()`; argument `prob_birth_death`, tests `dp_global/tests/test_probabilistic_birth_death.R` | `dp_global/R/dp_probabilistic_matching.R`, `dp_global/src/transition_cost_rcpp.cpp` |
| Growth SD refit when it falls with size; recruitment rate per established tree | `estimate_bio_pars(recruit_rate_unit = …)`, tests `dp_global/tests/test_estimate_bio_pars.R` | `dp_global/R/dp_global_bio.R` |
| DBH recorded in classes, rounded down (BCI 1982/1985 < 55 mm) | `transition_cost_paired_rcpp(round_t, round_tp1, …)`, `compute_pairwise_log_likelihood(round_curr, round_next, …)`; argument `dbh_round_censuses` | `dp_global/src/transition_cost_rcpp.cpp`, `dp_global/R/dp_global_dp.R`, `dp_global/R/dp_probabilistic_matching.R` |
| Carried-terminal backfill (post-engine) | `apply_carried_terminal_backfill()` | `dp_global/R/dp_global_main.R` |
| Orphan-stem backfill (post-engine) | `apply_orphan_stem_backfill()` | `dp_global/R/dp_global_main.R` |
| Rows left behind by the pin sweep rejoin their track's pin (post-engine) | `apply_pin_track_rejoin()` | `dp_global/R/dp_global_main.R` |
| Terminal records back to their stem (post-engine) | `apply_terminal_to_host()` | `dp_global/R/dp_global_main.R` |
| Stems split by a moved point of measurement rejoined (export and posterior samples, BCI merge / posterior prep) | `apply_measurement_rejoin()`, `apply_measurement_rejoin_to_paths()` | `dp_global/R/measurement_rejoin.R` |
| Broken-below invariants (post-engine) | `apply_broken_below_invariants()` | `dp_global/R/dp_global_main.R` |
| Per-sample database pins (used by `finalize_posterior_paths`) | `apply_pins_to_samples()` | `dp_global/R/dp_global_main.R` |
| Per-sample BB invariants (used by `finalize_posterior_paths`) | `apply_bb_invariants_to_samples()` | `dp_global/R/dp_global_main.R` |
| Renumbering of every tag's IDs to 1..N (post-engine; audit labels with no final ID get numbers after N) | `renumber_engine_minted_ids()` | `dp_global/R/dp_global_main.R` |
| Finalize posterior path files (post-renumber) | `finalize_posterior_paths()` | `dp_global/R/dp_global_main.R` |
| Parameter estimation | `estimate_bio_pars()` | `dp_global/R/dp_global_bio.R` |
| DP complexity estimate | `estimate_dp_complexity()` | `dp_global/R/complexity/estimate_dp_complexity_function.R` |
| Plotting | `plot_tag_to_pdf()` | `dp_global/R/dp_global_diag.R` |
| Driver (interactive / single-tag) | `run_main()`; per-tag call `run_dp_one_group()` | `dp_global/scripts/main_cpp.R` |
| Driver (chunked / large runs) | `run_main_chunked()` | `dp_global/scripts/main_cpp_chunk.R` |
| Driver (BCI debug, single-tag) | (sources `main_cpp.R`) | `dp_global/scripts/main_cpp_bci.R` |

### match_stems_dp_global_backward_marginals_batch — Function reference (implementation details)

This function is the DP solver: it returns the MAP reconstruction, the posterior marginals and, on request, posterior samples of one tag. Below are the key computations, the helper functions it uses, and the semantics of important parameters (including the pruning controls).

Key computations and helpers:

- Anchor selection
  - If the requested `anchor_start` has no DBH/TrueStemID, the function searches backward for the most recent census with at least one row having non-NA DBH and non-NA TrueStemID and uses that as the anchor. If the tag has measurements after the anchor and the anchor has 0 living stems, the function instead searches forward for the first post-anchor census with living stems (preferring one with a non-NA TrueStemID). If no valid anchor is found, a provisional anchor is used when allowed; otherwise it falls back to the probabilistic matcher.

- State enumeration
  - `enumerate_states_injective(K, n_obs, max_states)` enumerates injective assignments (permutation-based states) and returns `NULL` when there are more than `max_states`, which sends the tag to the probabilistic matcher. `enumerate_states_constrained(K, n_obs, allowed_tracks, max_states)` is used when pins or growth bounds restrict the tracks of a census; it stops at `max_states` assignments (see *Understanding `max_states`*).
  - `count_injective_states(K, n_obs)` computes theoretical counts used for diagnostics and for the segment split consistency check.

- Track DBH matrices
  - For each census the function builds an `n_states × K` matrix holding the DBH on each track in each state (`track_dbh_by_state`), used for the feasibility test and the cost evaluation.

- Phase constraints
  - `derive_phase_prev_batch_rcpp(...)` (C++, `dp_global/src/transition_cost_rcpp.cpp`) checks phase-transition feasibility, the R-recruit constraint and the hard pruning bounds for all (i, j) state pairs in batch; returns `from_i`, `to_j`, and the `phase_t` matrix of the feasible pairs.

- Interval computation
  - Uses per-census mean `ExactDate` to compute `interval_val` (years) between census pairs. If `interval_val` is NA, non-finite or not positive, the growth and recruit-size pruning are skipped for that pair.

- Conservative pruning (pre-filters)
  - Controlled by: `prune_hard` (logical); `prune_min_growth`, `prune_max_growth`, `prune_use_bio_bounds`, `prune_recruit_max_dbh`, `prune_use_bio_recruit`; `non_taper_corrected_growth_forms`, `non_taper_corrected_prune_min/max_growth`, `hom_tolerance_scale`.
  - Effective pruning thresholds computed as

    user_min = prune_min_growth (if given) else min_growth

    user_max = prune_max_growth (if given) else max_growth

    if prune_use_bio_bounds:
      eff_min_grow = max(user_min, Bio_Max_Shrink)
      eff_max_grow = min(user_max, Bio_Max_Growth)
    else:
      eff_min_grow = user_min
      eff_max_grow = user_max

    # Non-taper-corrected override: replace the effective bounds with the
    # non-taper bounds so palms/strangler figs/tree ferns are not spuriously
    # pruned. This step matters when prune_use_bio_bounds=TRUE (which
    # would tighten the general bounds); in the default driver config
    # (prune_use_bio_bounds=FALSE, both pairs of bounds = PRUNE_BOUND_FACTOR ×
    # the fixed limits) it is a no-op because the values are the same

    if growth_form in non_taper_corrected_growth_forms:
      eff_min_grow = non_taper_corrected_prune_min_growth  (override)
      eff_max_grow = non_taper_corrected_prune_max_growth  (override)

    # HOM widening (non-taper growth forms only): per census pair, widen
    # the bounds symmetrically by the worst-case HOM deviation across both
    # censuses. This is the main source of differentiation for non-taper
    # forms in the default config

    Per census pair, if non-taper-corrected, HOM column present and hom_tolerance_scale > 0:
      hom_tol = hom_tolerance_scale × max(|HOM − 1.3|) / interval_years
      eff_min_grow_pair = eff_min_grow − hom_tol
      eff_max_grow_pair = eff_max_grow + hom_tol

    eff_recruit_max chosen from `prune_recruit_max_dbh` and `Bio_Recruit_MaxDBH_unit` depending on flags.

  - For each candidate (assignment-state pair) and each track with DBH at both times compute

    g = (D_{t+1} - D_t) / Δt

    and prune if g ∉ [eff_min_grow, eff_max_grow]; for recruits (NA→DBH) prune if DBH > eff_recruit_max.

  - Pruning updates `prune_stats` and avoids costing the pruned candidates.

- Transition cost computation
  - `transition_cost_paired_rcpp(tdbh0_mat, tdbh1_mat, interval_years, ...)` computes, in one call per census pair, the cost of every feasible state pair (sum across tracks) including growth likelihood, mortality, recruitment, and hard penalties (1e6) for impossible transitions.

- Backward recursion & Viterbi
  - Uses `log_sum_exp` to accumulate marginal weights and a Viterbi update to compute MAP path (per-step `vit_cost`, `vit_ptr` arrays).

- Forward pass & posterior marginals
  - Builds per-observation posterior distributions (`DP_PosteriorTop{k}ID` / `DP_PosteriorTop{k}Prob`, entropy, the probability of the reconstructed ID and the unlinked probability) via normalized state weights.

- Posterior sampling
  - Optional sample drawing from the DP graph (`posterior_samples`), with backward sampling and per-sample `logp` weights.
  - The engine writes the raw per-sample data.table (in engine ID space) to `<out_dir>/posteriors/.staging/tag_<Tag>_samples_raw_<BATCH_TS>.rds` and attaches `attr(out, "DP_Posterior_Staging_File")` / `attr(out, "DP_Sampling_Profile")` for diagnostics. The final `*_paths.<ext>` file is written downstream by `finalize_posterior_paths()` after `renumber_engine_minted_ids()` has translated engine IDs into the renumbered space (see `dp_global/R/dp_global_main.R`). The probabilistic engine's `export_probabilistic_posteriors()` writes a structurally identical staging file.

- Diagnostics & fallbacks
  - `attr(out, "DP_PruneInfo")` contains pruning diagnostics (counts and effective thresholds), and `attr(out, "DP_Compute_Profile")` the number and time of the transition-cost calls. The function falls back to `match_stems_probabilistic()` on anchor failures, enumeration exhaustion, K insufficiency or if DP produces no feasible states after pruning (see *DP fallback reason codes*). These fallback return values always include `DP_PruneInfo`.

Notes:

- `min_growth`/`max_growth` are the prune bounds of the DP when `prune_min_growth`/`prune_max_growth` are `NULL`, and they are passed to the probabilistic matcher as base growth bounds together with `prune_min_growth`/`prune_max_growth`/`prune_recruit_max_dbh`. They also define the post-hoc `ConstraintViolation` check.

---

### Key Implementation Details

**NA interpretation:** inside the DP an empty track (`NA`) means "no stem," not "alive but missed"; censuses detected as missing-from-field are taken out before the DP (see *Missing-from-Field (MF) Detection*)

**Large penalties:** `1e6` encodes hard constraints (impossible transitions)

**Measurement error:** Controlled by `use_measurement_error` flag, affects:

- Likelihood computation (4-component mixture)
- Parameter estimation (variance subtraction, quantile selection)

**Phase variables:** Extended DP key format: `"<assignment>|<phases>"`, where `<assignment>` is the comma-separated track of each observation and `<phases>` one digit (0/1/2) per track

**Anchor initialization:** Tracks without an anchor stem start in phase 2 (dead / unused), allowing reconstruction of earlier deaths

### Debugging Tools

**Transition cost breakdown:**

```r
components <- transition_cost_tracks_bio_components(
  track_dbh_t, track_dbh_tp1, interval_years, ...
)
# components$per_track: one row per track with its case and cost terms (cost_recruit,
#   cost_no_recruit, cost_mortality, cost_growth_lik, cost_shrink_soft, cost_growth_soft,
#   cost_hard, total_track); components$tiebreak; components$total; components$p_recruit
# (no DBH-rounding arguments: equals the compiled cost for censuses not flagged as rounded)
```

**Posterior analysis:**

```r
# After running marginal solver:
high_uncertainty <- out[DP_PosteriorEntropy > 1.5]
ambiguous <- out[DP_PosteriorBin == "ambiguous"]
```

### Complexity Estimation

The DP algorithm's computational complexity depends on the number of observations per census and the resulting state space size. Use the complexity estimator to identify which tags will take the longest to process before running the full workflow.

**Load the complexity estimation functions** (from the project root; the file needs the `here` package attached):

```r
library(here)
source(here("dp_global", "R", "complexity", "estimate_dp_complexity_function.R"))
```

**Estimate complexity for all tags:**

```r
complexity <- estimate_dp_complexity(here("data_simulation", "data", "simulated_data_1.csv"))
print(complexity)
```

**Get detailed analysis for a specific tag:**

```r
details <- get_tag_complexity_details(here("data_simulation", "data", "simulated_data_1.csv"), tag = 11)
print(details$summary)     # the tag's row of estimate_dp_complexity()
print(details$per_census)  # CensusID, N_Obs, N_States = P(K, N_Obs)
```

**Output columns:**

- `Tag`: The tag identifier
- `Species`: Species name
- `census_range`, `n_censuses`: The censuses with a DBH up to the anchor
- `K`: Number of tracks assumed by the estimator
- `max_obs`: Maximum number of observations in any single census
- `max_states_per_census`: Maximum states in any single census
- `total_states`: Total number of states across all censuses
- `estimated_edges_unpruned`: Number of candidate transitions between adjacent censuses (proxy for runtime)
- `estimated_edges_pruned`: The same after the hard pruning bounds
- `estimated_fallback`, `fallback_reason`: Whether a census exceeds `max_states` (the DP would send the tag to the probabilistic matcher)
- `predicted_seconds`, `predicted_hours`: Runtime predicted from the edge count by a fitted log-log polynomial

**How it works:**

1. **Number of tracks (K)**: `max(#anchor IDs, max_obs, n_obs(1) + births)` plus one track per resprout observation plus slack tracks. The DP adds no track for resprout observations and adds tracks for `TrueStemID`s pinned before the anchor, so the estimator's K can differ from the K the DP uses.
2. **States per census**: P(K, n_obs) = K! / (K - n_obs)! where n_obs is observations in that census
3. **Transitions**: For each pair of consecutive censuses, all pairs of states are counted (and, for the pruned count, tested against the growth and recruit-size bounds)

Tags without fallback come first, sorted by the number of candidate transitions in descending order, so the slowest tags appear first.

**Example output for the simulated data** (`simulated_data_1.csv`, defaults): of 131 tags, 35 are flagged as fallback. Tag 19 has 107,956,800 candidate transitions and tag 11 has 57,153,600 (both K = 7, up to 6 stems per census); every other tag has at most 2,462,400.

---


## Notes & Common Issues

### Interpretation Notes

1. **DBH → NA is mortality** inside the engines, not a missing measurement (except in censuses detected as missing-from-field, which are set aside before the engine)
2. **K selection** is critical: too small fails, too large is inefficient
3. **Slack tracks** enable realistic death+birth dynamics in constant-count transitions
4. **Temperature** in posterior controls confidence (not a physical parameter)

### Parameter Roles Summary Table

| Parameter Type | Parameter Name | Used By | Purpose | Typical Value |
|----------------|---------------|---------|---------|---------------|
| **Hard Constraint (DP)** | `Bio_Max_Shrink` | DP cost function | Hard lower bound on $g$ | -0.5 cm/year (`MAX_SHRINK_FIXED`) |
| **Hard Constraint (DP)** | `Bio_Max_Growth` | DP cost function | Hard upper bound on $g$ | 5 or 7.5 cm/year (`MAX_GROWTH_FIXED`) |
| **Soft Penalty (DP)** | `Bio_K_Shrink` | DP cost function | Quadratic shrinkage weight | 0 (drivers) to 50 (1/cm²) |
| **Soft Penalty (DP)** | `Bio_K_Growth` | DP cost function | Quadratic growth excess weight | 0 (drivers) to 50 (1/cm²) |
| **Soft Threshold (DP)** | `Bio_Max_Growth_Soft` | DP cost function | Growth excess threshold | min(`Bio_Max_Growth`, 99% quantile of observed growth) |
| **Base bounds** | `min_growth` | DP pruning (when `prune_min_growth` is `NULL`), probabilistic matcher | Base lower growth bound | -0.5 cm/year |
| **Base bounds** | `max_growth` | DP pruning (when `prune_max_growth` is `NULL`), probabilistic matcher | Base upper growth bound | 7.5 cm/year |
| **Prune bounds** | `prune_min_growth` | DP pruning, probabilistic matcher | Hard lower pruning bound (`PRUNE_BOUND_FACTOR` × base) | -0.625 cm/year |
| **Prune bounds** | `prune_max_growth` | DP pruning, probabilistic matcher | Hard upper pruning bound (`PRUNE_BOUND_FACTOR` × base) | 9.375 cm/year |
| **Prune bounds** | `prune_recruit_max_dbh` | DP pruning, probabilistic matcher | Hard recruit size cap (`PRUNE_BOUND_FACTOR` × `RECRUIT_MAX_FIXED`) | 48.12 cm |
| **Diagnostic Only** | `min_growth` | `ConstraintViolation` | Post-hoc violation flag | Same as base bounds |
| **Diagnostic Only** | `max_growth` | `ConstraintViolation` | Post-hoc violation flag | Same as base bounds |
| **Probabilistic matcher** | `prob_n_samples` | Probabilistic matcher | Number of Gumbel-noise samples | 200 |
| **Probabilistic matcher** | `prob_lookahead_weight` | Probabilistic matcher | Sequential backward conditioning weight | 0.5 (function default), 1 (drivers) |

The typical values are those of `main_cpp.R` (`MAX_GROWTH_FIXED = 7.5`, `PRUNE_BOUND_FACTOR = 1.25`, `RECRUIT_MAX_FIXED = 7.5 × 5 + 0.9999`).


### Quantile Configuration Reference

Complete table of all quantiles used in parameter estimation:

| Quantile Parameter | Default Value | Tail | What It Computes | Used For | Configurable |
|-------------------|---------------|------|------------------|----------|--------------|
| `shrink_data_quantile` | 0.001 (0.1%) | Lower | Empirical $g$ lower bound | `max_shrink_data` | Yes, via `estimate_bio_pars()` |
| `shrink_hard_prob` | $10^{-4}$ (0.01%) | Lower | Measurement-only $g$ lower bound | `max_shrink_meas` | Yes, via `estimate_bio_pars()` |
| `growth_data_quantile` | 0.999 (99.9%) | Upper | Empirical $g$ upper bound | `max_growth_data` | Yes, via `estimate_bio_pars()` |
| `growth_hard_prob` | $10^{-4}$ (0.01% tail) | Upper | Measurement-only $g$ upper bound (99.99%) | `max_growth_meas` | Yes, via `estimate_bio_pars()` |
| `growth_soft_quantile` | 0.99 (99%) | Upper | Empirical $g$ soft threshold | `max_growth_soft_data` | Yes, via `estimate_bio_pars()` |
| `recruit_max_quantile` | 0.999 (99.9%) | Upper | Maximum recruit DBH | `recruit_max_dbh` | Yes, via `estimate_bio_pars()` |

**Key distinctions:**

1. **Empirical quantiles** (data_quantile): Computed from observed $g_{\text{all}}$ values
2. **Measurement-only quantiles** (hard_prob): Computed from 4-component mixture at typical diameter
3. **Combination logic:**
   - Shrinkage: `min()` = more conservative (less negative = more restrictive)
   - Growth: `max()` = more permissive (larger = less restrictive)

### Common Confusions

1. **"Why huge penalties?"**
   - $10^6$ encodes hard biological impossibilities
   - Makes infeasible transitions effectively infinite cost
   - DP will never choose these paths unless no alternative exists

2. **"Why track-based?"**
   - Allows global optimization over entire history
   - Prevents local greedy decisions from causing global ID swaps
   - Enables exact life-cycle constraint enforcement

3. **"Data vs fixed mode?"**
   - **Data mode:** Use when you have reliable known IDs for parameter estimation
   - **Fixed mode:** Use for literature values or when data is too sparse
   - Can mix (e.g., data hard, fixed soft)

4. **"Measurement error impact?"**
   - Makes penalties more conservative (won't over-penalize measurement noise)
   - Adds 4-component mixture to growth likelihood (more robust)
   - Separates process from measurement variance in estimation

5. **"min_growth vs Bio_Max_Shrink?"**
   - `Bio_Max_Shrink`: Used in DP cost (primary reconstruction constraint)
   - `min_growth`: Used for pruning (when no prune bound is given), by the probabilistic matcher and for the `ConstraintViolation` diagnostic
   - They can have different values

6. **"Why do soft penalties exist if hard constraints are enforced?"**
   - Hard: Binary (allowed/forbidden), creates discrete cost landscape
   - Soft: Continuous discouragement, helps DP choose among many feasible paths
   - Together: Hard prevents absurdity, soft guides toward realism

### Performance Considerations

- State space grows factorially: $P(K, n_{obs})$ states per census, at most $K^{n_{obs}}$
- Set `DP_MAX_STATES` based on memory availability (default 40,000 in `main_cpp.R` and `main_cpp_chunk.R`); also controls the cross-product edge limit (`max_states²`)
- Set `PROB_N_SAMPLES` to control probabilistic matching accuracy vs speed (default 200)
- Parallel processing via `MANUAL_CORES_VALUE` (controlled by `MANUAL_CORES=TRUE`) for multiple tags
- Consider fallback threshold adjustment for large datasets
- The DP costs all feasible transitions of a census pair in one C++ call (`transition_cost_paired_rcpp()`)

### Numerical Stability Details

| Operation | Stability Mechanism | Value |
|-----------|-------------------|-------|
| Growth SD | Floor clamp | $\sigma(D) \geq 10^{-6}$ |
| Measurement SD1 | Floor clamp | $\text{SD1}(D) \geq 10^{-6}$ |
| Mortality probability | Range clamp | $[10^{-12}, 1-10^{-12}]$ |
| Recruitment probability | Range clamp | $[10^{-12}, 1-10^{-12}]$ |
| Log-sum-exp mixture | Subtract max before exp | Prevents overflow/underflow |
| K estimation (data mode) | Range clamp | $k \in [10^{-6}, 10^6]$ |
| Rank tie-breaking | Ranks from a stable sort (equal DBHs keep track order) | `std::stable_sort` |

---

## References

**Measurement error model:**
Chave, J., Condit, R., Aguilar, S., Hernandez, A., Lao, S., & Perez, R. (2004). Error propagation and scaling for tropical forest biomass estimates. *Philosophical Transactions of the Royal Society of London. Series B: Biological Sciences*, 359(1443), 409-420.

**DBH recorded in 5 mm classes (BCI 1982 and 1985):**
CTFS R Package tutorial, growth changes: <https://ctfs.si.edu/ctfsdev/CTFSRPackageNew/index.php/web/tutorials/GrowthChange/index.html>
Piponiot, C., Condit, R., Hubbell, S. P., Pérez, R., Lao, S., Aguilar, S., & Muller-Landau, H. (2024). Woody Biomass Stocks and Fluxes in the Barro Colorado Island 50-ha Plot. (Appendix S1 R code.)

**Additional documentation:**

- Driver scripts and options: `dp_global/scripts/README.md`
- C++ functions: `dp_global/src/README.md`
- Sensitivity analysis: `dp_global/R/sensitivity_transition_cost_bio.R`
- Realism calibration: `dp_global/R/realism_calibration.R`
- Basal area uncertainty: `dp_global/scripts/basal_area_uncertainty.R`
- API documentation: See function headers in the relevant R files (`dp_global/R/dp_global_dp.R`, `dp_global/R/dp_global_bio.R`, etc.)

---

## Basal Area Uncertainty Quantification

The `basal_area_uncertainty.R` script uses posterior path samples to quantify how identity uncertainty propagates into basal area (BA) estimates.

### Key Insight

**Tag-level BA per census is invariant to identity assignment** — the total BA (sum of $\pi/4 \cdot (\text{DBH}/100)^2$ across all living stems, in m²) does not change because the same DBH observations are always present regardless of which stem identity they are assigned to. Identity uncertainty affects:

1. **Per-stem BA trajectories**: Different identity assignments allocate different DBH values to each stem across censuses
2. **BA growth rates**: Different trajectories produce different consecutive BA differences and thus different growth rates
3. **Demographic accounting**: Mortality and recruitment events change with identity assignment

### Usage

```bash
Rscript dp_global/scripts/basal_area_uncertainty.R \
  --RUN_DIR=dp_global/output/<run_dir>
```

The script reads `stem_reconstruction_dp_global_rcpp.csv` and the posterior path files `posteriors/tag_<Tag>_posterior_samples_*_paths.csv` of the run directory. Only CSV path files are read (`POSTERIOR_SAMPLES_FORMAT = "csv"`, the default of the dp_global drivers); without them only the MAP columns are produced. Paths are weighted by their `path_prob`. DBH is in cm.

### Outputs

Two CSV files and one PDF are written to the run directory:

| File | Contents |
|------|----------|
| `basal_area_tag_census.csv` | Per-tag × census: total BA (m²), stem count, year |
| `basal_area_tag_change.csv` | Per-tag BA change between consecutive censuses: MAP decomposition into growth (surviving stems), loss (mortality), and gain (recruitment), plus posterior mean, SD, and 95% CI for each component |
| `basal_area_figures.pdf` | Multi-page PDF with summary pages, per-tag detail panels (BA trajectory, stem count, decomposition bars with uncertainty whiskers, stem demographics), uncertainty histograms, and **posterior density plots** (kernel densities of Growth/Loss/Gain/DeltaBA with weighted-mean vertical lines — one page pooled across all census intervals, then one page per interval) |

### Column Reference

**Tag census** (`basal_area_tag_census.csv`):

| Column | Description |
|--------|-------------|
| `Tag` | Individual tree identifier |
| `CensusID` | Census identifier |
| `Year` | Calendar year (from mean ExactDate) |
| `TotalBA_m2` | Total basal area = Σ(π/4 × (DBH/100)²) across all stems (m²) |
| `NumStems` | Number of stems with non-NA DBH |

**Tag change** (`basal_area_tag_change.csv`):

| Column | Description |
|--------|-------------|
| `Tag` | Individual tree identifier |
| `CensusID_from`, `CensusID_to` | Census pair defining the interval |
| `Interval_yr` | Interval length in years (from mean ExactDate) |
| `Growth_BA` | MAP: BA change attributable to surviving stems (m²) |
| `Loss_BA` | MAP: BA removed by stem mortality (negative, m²) |
| `Gain_BA` | MAP: BA added by stem recruitment (m²) |
| `DeltaBA_total` | Total BA change = Growth + Loss + Gain (invariant to identity) |
| `NumSurvivors`, `NumDeaths`, `NumRecruits` | MAP stem counts per demographic category |
| `Growth_mean`, `Growth_sd`, `Growth_q025`, `Growth_q975` | Posterior summary of growth component |
| `Loss_mean`, `Loss_sd`, `Loss_q025`, `Loss_q975` | Posterior summary of loss component |
| `Gain_mean`, `Gain_sd`, `Gain_q025`, `Gain_q975` | Posterior summary of gain component |
| `DeltaBA_check` | Posterior weighted mean of total BA change (should match `DeltaBA_total`) |
| `NumSurvivors_mean`, `NumDeaths_mean`, `NumRecruits_mean` | Posterior weighted mean stem counts per demographic category |
| `NumPaths` | Number of posterior paths used for uncertainty estimation |

---

## Building This Documentation

Use `pandoc` (if available) to render `README.md` to HTML, or use your preferred tooling.
