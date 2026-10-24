# dp_global

Global dynamic-programming reconstruction of stem identities with biological
costs.

In each census a tree (tag) has a set of stem measurements, but which
measurement belongs to which stem is known only from some census on (the
anchor census: field tags, database IDs). The engine links the earlier
measurements to stems under a biological model of growth, mortality and
recruitment:

- an exact dynamic program (DP) finds the most probable reconstruction of the
  tree, its posterior marginals and posterior samples of whole
  reconstructions;
- a probabilistic matcher with the same model solves the trees the DP cannot
  (too many stems, no usable anchor, or routed there by species or growth
  form);
- post-engine helpers apply known identities (pins) and the ForestGEO field
  codes (breaks, resprouts, stumps), then number the stems 1..N per tree.

`BCI_stem_reconstruction/2_STEM_IDENTIFICATION/` runs the engine on the BCI
50-ha plot.

Key features:

- global optimisation over all censuses of a tree, not census by census;
- biological realism: size-dependent growth, mortality and recruitment;
- life-cycle constraints: a stem that dies does not come back;
- measurement error following Chave et al. (2004), and DBH recorded in
  classes;
- exact posterior marginals and posterior samples (DP), and per-observation
  posterior probabilities from perturbed samples (matcher);
- transition costs, feasibility checks and the assignment solver in C++.

Author: José A. Medina-Vega

## Contents

| File | What it does |
| --- | --- |
| `R/dp_global_main.R` | Loads the engine (compiles the C++ file, sources the modules) and defines the post-engine helpers |
| `R/dp_global_dp.R` | DP solver `match_stems_dp_global_backward_marginals_batch()`, resprout split, posterior staging |
| `R/dp_probabilistic_matching.R` | Probabilistic matcher `match_stems_probabilistic()` |
| `R/dp_global_bio.R` | Parameter estimation `estimate_bio_pars()` and the cost breakdown `transition_cost_tracks_bio_components()` |
| `R/dp_global_states.R` | State enumeration and track-DBH helpers |
| `R/dp_global_utils.R` | State counts and census intervals |
| `R/dp_global_diag.R` | Posterior bins, `ConstraintViolation` and trajectory plots |
| `R/measurement_rejoin.R` | Rejoin of stems split by a moved point of measurement (run by the BCI merge) |
| `R/naming_helpers.R` | Run folder names |
| `R/dp_global_matchers.R` | Stepwise igraph matcher (not used by the DP workflow) |
| `R/sensitivity_transition_cost_bio.R` | Sensitivity sweeps of the transition cost |
| `R/realism_calibration.R` | Realism report of a reconstruction |
| `R/k_tuning_viz.R` | Join-versus-split cost as the soft-penalty weight varies |
| `R/check_functions.r` | Inspection and plots of estimated parameters |
| `R/complexity/` | DP state counts and run-time estimate per tag |
| `R/dpglobal_bundle/` | Packaging of the engine for other machines |
| `scripts/` | Drivers (`main_cpp.R`, `main_cpp_chunk.R`, `main_cpp_bci.R`) and `basal_area_uncertainty.R` |
| `src/` | C++ transition costs, feasibility checks and assignment solver, with R wrappers |
| `tests/` | testthat tests |

```text
dp_global/
├── README.md
├── R/
│   ├── dp_global_main.R              # loader + post-engine helpers
│   ├── dp_global_dp.R                # DP solver
│   ├── dp_probabilistic_matching.R   # probabilistic matcher
│   ├── dp_global_bio.R               # parameter estimation + R cost breakdown
│   ├── dp_global_states.R            # state enumeration
│   ├── dp_global_utils.R             # state counts, intervals
│   ├── dp_global_diag.R              # bins, ConstraintViolation, plots
│   ├── measurement_rejoin.R          # moved point of measurement rejoin
│   ├── naming_helpers.R              # run folder names
│   ├── dp_global_matchers.R          # stepwise igraph matcher
│   ├── sensitivity_transition_cost_bio.R
│   ├── realism_calibration.R
│   ├── k_tuning_viz.R
│   ├── check_functions.r
│   ├── complexity/
│   │   ├── estimate_dp_complexity_function.R
│   │   └── profile_complexity.R
│   └── dpglobal_bundle/              # bundle builder; dist/ is built, not tracked
├── scripts/
│   ├── README.md
│   ├── main_cpp.R                    # single-pass driver
│   ├── main_cpp_chunk.R              # chunked, resumable driver
│   ├── main_cpp_bci.R                # BCI single-tag debug driver
│   └── basal_area_uncertainty.R
├── src/
│   ├── README.md
│   ├── transition_cost_rcpp.cpp
│   └── transition_cost_rcpp.R
└── tests/
    ├── test_enumerate_states.R
    ├── test_estimate_bio_pars.R
    └── test_probabilistic_birth_death.R
```

Run outputs go to `dp_global/output/` (not tracked).

## How to run

Run everything from the project root (`dp_global_main.R` finds its files from
the working directory).

```sh
Rscript dp_global/scripts/main_cpp.R          # one tag or a small data set
Rscript dp_global/scripts/main_cpp_chunk.R    # large data sets: chunks, incremental output, resume
```

```r
source("dp_global/scripts/main_cpp.R")   # settings and functions; starts no run
res <- run_main()
```

`scripts/README.md` lists the driver settings. The BCI run uses its own
driver, `BCI_stem_reconstruction/2_STEM_IDENTIFICATION/1_main_cpp_chunk_bci.R`.
Census timing comes from `ExactDate`; there is no interval argument.

### Calling the engine directly

The example runs on the simulated data (`data_simulation/data/simulated_data_1.csv`).
The drivers do the same, then apply the post-engine helpers (section
"Post-engine helpers").

```r
library(data.table)
library(here)
source(here("dp_global", "R", "dp_global_main.R"))

xraw <- fread(here("data_simulation", "data", "simulated_data_1.csv"))
xraw[, species := as.character(Species)]
USE_MEASUREMENT_ERROR <- TRUE

# 1. Parameters per species (intervals come from ExactDate)
bio_pars_list <- list()
for (sp in unique(xraw$species)) {
  bio_pars_list[[sp]] <- estimate_bio_pars(
    xraw[species == sp],
    anchor_start_census = 7L,          # first census used
    use_measurement_error = USE_MEASUREMENT_ERROR,
    max_shrink_source = "data", max_shrink_fixed = -0.5,   # hard bounds: "data" or "fixed"
    max_growth_source = "data", max_growth_fixed = 7.5,
    k_shrink_source = "data", k_shrink_fixed = 50,         # soft penalties: "data" or "fixed"
    k_growth_source = "data", k_growth_fixed = 50,
    recruit_max_source = "data", recruit_max_fixed = 5,    # recruit size cap
    shrink_data_quantile = 0.001, shrink_hard_prob = 1e-4, # quantiles (defaults)
    growth_data_quantile = 0.999, growth_hard_prob = 1e-4,
    growth_soft_quantile = 0.99, recruit_max_quantile = 0.999
  )
}

# 2. Bio_* columns (one value per tag)
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

# 3. The DP, one tag at a time
out <- rbindlist(lapply(unique(xraw$Tag), function(tg) {
  match_stems_dp_global_backward_marginals_batch(
    tree_data = copy(xraw[Tag == tg]),
    min_growth = -2, max_growth = 10,  # pruning (when prune_* are NULL), matcher, ConstraintViolation
    anchor_start = 7, max_tracks = 30, max_states = 50000, slack_tracks = 1,
    temperature = 1, posterior_top_k = 2,
    use_measurement_error = USE_MEASUREMENT_ERROR
  )
}), fill = TRUE)

# 4. Posterior bins
out <- add_dp_posterior_bins(out, confident_prob = 0.95, unlinked_prob = 0.50,
                             use_reconstructed_prob = TRUE)
```

In the drivers the hard bounds come from `MAX_SHRINK_FIXED` /
`MAX_GROWTH_FIXED` (source `"fixed"`), the recruit cap from
`RECRUIT_MAX_FIXED` (`MAX_GROWTH_FIXED × 5 + 0.9999` cm), and both soft
penalties are fixed at 0.

Things to keep in mind:

1. `min_growth` / `max_growth` are not part of the transition cost. They are
   the DP's pruning bounds when `prune_min_growth` / `prune_max_growth` are
   `NULL`, the matcher's growth bounds (replaced by the prune bounds when
   those are given), and the bounds of the `ConstraintViolation` diagnostic.
2. The constraints inside the cost come from the `Bio_*` columns:
   `Bio_Max_Shrink`, `Bio_Max_Growth` (hard), `Bio_K_Shrink`, `Bio_K_Growth`
   and `Bio_Max_Growth_Soft` (soft).
3. `use_measurement_error` must be the same in `estimate_bio_pars()` and in
   the engine.
4. The returned parameter list stores each parameter's source and value.

### Usage patterns

```r
# One tag, with plots (main_cpp.R or --RUN_ALL_TAGS=FALSE --WHICH_TAG=19)
RUN_ALL_TAGS <- FALSE; WHICH_TAG <- "19"; WRITE_DP_PDF <- TRUE

# Posterior uncertainty (defaults)
DP_MODE <- "marginals+bins"; POSTERIOR_SAMPLES <- 200L   # 0 = no path files

# All tags in parallel
RUN_ALL_TAGS <- TRUE; MANUAL_CORES <- TRUE; MANUAL_CORES_VALUE <- 8L
DP_MAX_STATES <- 1e5; DP_MAX_TRACKS <- 50L

# Sensitivity sweeps only (or main_cpp.R --SENSITIVITY_MODE=run+write+pdf)
source("dp_global/R/dp_global_main.R")
source("dp_global/R/sensitivity_transition_cost_bio.R")   # build_all_sweeps(), plot_all_sweeps_to_pdf()

# Realism report only (or main_cpp.R --RUN_REALISM_REPORT=TRUE)
source("dp_global/R/realism_calibration.R")   # realism_report_from_reconstruction(out, interval_years, base_args)

# Growth bounds (cm/yr) and recruit cap (cm) enforced before the fit
bio <- estimate_bio_pars(x, anchor_start_census = 7L,
  enforce_growth_bounds = TRUE, growth_min_fixed = 0.05, growth_max_fixed = 0.75,
  enforce_recruit_max = TRUE, recruit_max_source = "fixed", recruit_max_fixed = 37.5)
```

## Inputs

One row per stem record (a stem in a census).

| Column | Type | Content |
| --- | --- | --- |
| `Tag` | integer / character | Tree: the unit the engine solves |
| `CensusID` | integer | Census number 1, 2, 3, … |
| `DBH` | numeric | Diameter in cm (`NA` = not measured) |
| `ExactDate` | Date | Measurement date; the interval of a census pair is the difference of the mean census dates of the tag |
| `species` | character | Parameter set of the row (the drivers create it) |
| `TrueStemID` | integer | Pin: known stem identity, at the anchor census and at any earlier row to pin; `NA` where the engine is free |

Used when present:

| Column | Used for |
| --- | --- |
| `Status`, `ListOfTSM` | Field status and codes: breaks, resprouts, stumps, dead and missing stems |
| `growth_form`, `hom` / `HOM` | Pruning bounds of non-taper-corrected forms; routing to the matcher |
| `Species` | Routing species to the matcher (`prob_species`) |
| `PrimaryStem` | Main-stem rule (`"main"`) |
| `StemTag` | Field tags (hard pins) |
| `StemID` or `OriginalStemID` | Database stem ID: post-engine backfills and grouping of terminal records |

If DBH was measured at different heights, give the taper-corrected DBH at
1.3 m for the taper-corrected growth forms. Forms that cannot be
taper-corrected (palms, strangler figs, tree ferns) get wider pruning bounds
and HOM widening (section "Pruning").

### Biological parameter columns

One value per tag, usually from `estimate_bio_pars()`.
`Bio_Gamma_Growth`, `Bio_Max_Growth`, `Bio_Max_Growth_Soft` and `Bio_K_Growth`
are optional (defaults 0, no bound, no cap and 0); the DP stops when another
column is missing, `NA` or not constant within the tag.

| Column | Symbol | Content |
| --- | --- | --- |
| `Bio_Mu_Growth` | $\mu_{\text{const}}$ (α) | Intercept of the mean annual growth (cm/yr) |
| `Bio_Gamma_Growth` | $\mu_\gamma$ (γ) | Slope of the mean growth on log DBH (cm/yr per unit log DBH) |
| `Bio_Sigma0_Growth` | $\sigma_0$ | Baseline SD of the growth process (cm/yr) |
| `Bio_Sigma1_Growth` | $\sigma_1$ | Slope of the growth SD on DBH ((cm/yr)/cm) |
| `Bio_Max_Shrink` | `max_shrink` | Hard lower bound of annual growth (cm/yr, negative) |
| `Bio_K_Shrink` | $k_{\text{shrink}}$ | Weight of the soft shrinkage penalty (1/cm²) |
| `Bio_Max_Growth` | `max_growth` | Hard upper bound of annual growth (cm/yr) |
| `Bio_Max_Growth_Soft` | `max_growth_soft` | Soft growth cap (cm/yr) above which the soft growth penalty applies |
| `Bio_K_Growth` | $k_{\text{growth}}$ | Weight of the soft growth penalty (1/cm²) |
| `Bio_H0_Mortality` | $h_0$ | Baseline mortality hazard (1/yr) |
| `Bio_Beta_Mortality` | $\beta$ | Size effect on the mortality hazard (1/cm) |
| `Bio_Recruit_Meanlog`, `Bio_Recruit_Sdlog` | | Lognormal size of a new stem |
| `Bio_Recruit_MaxDBH_unit` | | Largest plausible DBH of a new stem (cm) |
| `Bio_Recruitment_lambda` | $\lambda$ | Recruitment rate of an empty track (1/yr) |

## Outputs

Columns added to every row:

| Column | Content |
| --- | --- |
| `ReconstructedStemID` | Reconstructed stem, numbered 1..N within each tag |
| `ReconstructionMethod` | How the row got its stem (below) |
| `DP_PosteriorTop1ID`, `DP_PosteriorTop1Prob`, `DP_PosteriorTop2ID`, `DP_PosteriorTop2Prob` | The most probable stems of the row and their probabilities (`posterior_top_k` pairs) |
| `DP_PosteriorEntropy` | Entropy of the row's posterior, $H = -\sum_i p_i \log p_i$ |
| `DP_PosteriorReconstructedProb` | Posterior probability that the row belongs to its exported stem |
| `DP_PosteriorUnlinkedProb` | Posterior probability that the row belongs to no anchor stem (DP only; `NA` from the matcher) |
| `DP_PosteriorBin` | `confident`, `ambiguous` or `unlinked-likely` (`add_dp_posterior_bins()`) |
| `DP_KUsed` | Number of tracks |
| `DP_MaxStatesPerCensus`, `DP_MaxStatesCensusID` | Largest state count P(K, n) of a census, before track restrictions, and that census |
| `DP_FallbackReason` | Why the tag went to the probabilistic matcher; `NA` when the DP solved it |
| `ConstraintViolation` | Annual growth of the reconstructed stem outside `[min_growth, max_growth]` (diagnostic, `add_constraint_violation()`) |
| `obs_row_id` | Row number within the tag; the key of the posterior path files |
| `SweepAuditOverride` | `TRUE` where the pin sweep replaced the stem the engine had assigned |
| `ReconstructedStemID_PreSweep` | The engine's stem before the pin sweep. Set once by the first sweep that runs; the later helpers do not update it, so it can differ from the final stem on other rows too |
| `SweepRollbackToPreSweep` | `TRUE` where the engine's stem was kept: the pin would have put two rows of one census on one stem (`SweepAuditOverride` is then `TRUE` too, and `ReconstructionMethod` keeps the engine's label), or `apply_terminal_to_host()` restored the engine's choice |
| `DP_Chunk`, `run_out_dir` | Chunk number and run folder (added by the drivers) |

`ReconstructionMethod`:

| Value | Meaning |
| --- | --- |
| `given` | Stem from the pin (`TrueStemID`): anchor row, pinned row before the anchor, row after the anchor with a pin, or row set by the pin sweep |
| `dp` | Assigned by the DP |
| `probabilistic` | Assigned by the probabilistic matcher |
| `provisional_dp` | Provisional anchor stem (section "Anchor census") |
| `dp_mf_inferred` | Row of a census in which the stem was alive but not measured, set aside before the engine; it takes the stem of the one track that spans it, otherwise `ReconstructedStemID` stays `NA` |
| `pin_track` | Unpinned row moved onto the pin of its engine track (`apply_pin_track_rejoin()`) |
| `carried_terminal` | Unmeasured dead, stem-dead, broken-below or missing row given the stem of the previous record with its `StemID` (`apply_carried_terminal_backfill()`) |
| `given_orphan` | Unmeasured row of a stem never measured and never pinned, given its `StemID` (`apply_orphan_stem_backfill()`) |
| `terminal_to_host` | Unmeasured dead, stem-dead or broken-below record moved to the stem that ended (`apply_terminal_to_host()`) |
| `bb_split` | Measured broken-below row that starts a new stem by rule R1 |
| `bb_split_carry` | Later rows of a stem started by R1 |
| `bb_post_terminator_split` | First measured row after a stump that starts a new stem by rule R2 |
| `bb_post_terminator_split_carry` | Later rows of a stem started by R2 |
| `measurement_rejoin` | Row of a stem joined by `apply_measurement_rejoin()` (BCI merge) |
| `none_after_anchor` | Row after the anchor census without a pin or an engine stem |
| `skipped_no_data` | Tag without usable measurements (drivers) |

Rows after the anchor: a row with a pin takes it (`given`, measured or not);
unmeasured broken-below / resprout-coded rows without a pin that directly
follow the anchor (up to the first later measurement) take the anchor stem
when exactly one stem is measured at the anchor (`given`); the other rows stay
`none_after_anchor` with `ReconstructedStemID = NA`, and the post-engine
backfills may fill them.

### Posterior path files

With `posterior_samples > 0` every tag gets one file

```text
<out_dir>/posteriors/tag_<Tag>_posterior_samples_<BATCH_TS>_paths.<feather|rds|csv>
```

(`feather` needs `arrow`, else `rds`; an empty `BATCH_TS` gives a double
underscore). One row per distinct sampled reconstruction:

| Column | Type | Content |
| --- | --- | --- |
| `path_sig` | character | The per-sample stem label of every measured observation up to the anchor, in census order, joined by `-` |
| `path_count` | integer | Number of samples with this reconstruction; `path_count / n_samples` is its posterior probability |
| `path_prob` | numeric | Reference only (sums to 1). DP: each sample re-weighted by its own sampling probability `exp(logp)`, which counts that probability twice (≈ p²); matcher: equals `path_count / n_samples` |
| `recon` | character | `obs_row_id:label` pairs, separated by `;` |

The labels are per-sample labels, not the exported `ReconstructedStemID`, so
one grouping of observations can appear in several rows (e.g. a stem that ends
before the anchor on a different spare track in two DP samples). The
frequency of a grouping is the sum of `path_count` over its rows. `obs_row_id`
joins the paths to the reconstruction table.

How to use the samples:

- **DP**: drawn by backward sampling from the exact posterior, so
  `path_count / n_samples` estimates each reconstruction's probability;
  likely reconstructions repeat.
- **Matcher**: approximate. Each sample is a Gumbel-perturbed assignment per
  census pair, stitched across censuses, repaired and filtered for pins. With
  many stems nearly every sample differs (`path_count = 1`, weight
  1/`n_samples`). The spread measures identity uncertainty heuristically, not
  as a calibrated posterior.
- **Both together**: weight every path by `path_count / sum(path_count)`
  within its tree. Drawing a path with these weights is drawing one of the
  tree's samples uniformly, so trees of both engines can be sampled together
  (as `BCI_stem_reconstruction/4_EXAMPLE_STRUCTURE_ASSESSMENT/basal_area_uncertainty.R`
  does).

**Staging.** The engines do not write the path files. Each stages its raw
samples (engine ID space) in
`<out_dir>/posteriors/.staging/tag_<Tag>_samples_raw_<BATCH_TS>.rds`. After
the post-engine helpers, `finalize_posterior_paths()` reads each staging file
and:

1. translates the labels with the renumbering of `renumber_engine_minted_ids()`;
2. applies the database pins to each sample as the export applies them
   (`apply_pins_to_samples()`: pin sweep with duplicate rollback, then the
   track rejoin of `apply_pin_track_rejoin()`). Engines that do not honour
   every pin while sampling (a provisional anchor, the pre-segment of a
   resprout split) would otherwise group pinned rows differently from the
   export. Only pinned links change, so the spread of the unpinned links is
   kept;
3. applies the broken-below rules to each sample (`apply_bb_invariants_to_samples()`);
4. computes `path_sig`, `path_count`, `path_prob` and `recon`, writes the
   path file and deletes the staging file.

A completed run leaves `.staging/` empty.

**MAP versus samples.** The exported `ReconstructedStemID` is the MAP
reconstruction (DP: the Viterbi backtrace of the most probable full path) or
the most representative sample (matcher), after the pin sweep and the
helpers. The samples are finite and random, so the exported grouping can have
posterior mass and still be absent from the path file; more samples make it
more likely to appear. Compare a path with the export by grouping:

```r
library(data.table)
dp <- fread("dp_global/output/<run_dir>/stem_reconstruction_dp_global_rcpp.csv")
tg <- dp[Tag == 11]
paths <- fread("dp_global/output/<run_dir>/posteriors/tag_11_posterior_samples__paths.csv")
grouping <- function(ids) match(ids, ids)   # same grouping <=> identical vectors
has_export <- vapply(strsplit(paths$recon, ";", fixed = TRUE), function(p) {
  kv <- do.call(rbind, strsplit(p, ":", fixed = TRUE))
  exported <- tg$ReconstructedStemID[match(as.integer(kv[, 1]), tg$obs_row_id)]
  identical(grouping(exported), grouping(kv[, 2]))
}, logical(1))
paths[has_export]  # empty: the exported grouping was not sampled
```

### Files written by the drivers

In the run folder (`dp_global/output/<run>/` for the generic drivers; the
driver prints it at start):

| File | Content |
| --- | --- |
| `run_started.txt`, `run_finished.txt` | Start and end time stamps |
| `run_log.txt` | Log lines of the run |
| `run_parameters_full.txt` | Settings, command line, overrides and parameter structure |
| `stem_reconstruction_dp_global_rcpp.csv` / `.rds` / `.feather` / `.pdf` | Reconstruction and plots (single-pass run; `WRITE_DP_*`) |
| `stem_reconstruction_dp_global_rcpp_chunk_NNN.*`, `_done.txt`, `_failed.txt` | Chunked runs: per-chunk outputs; `_done.txt` lets a resumed run skip the chunk |
| `posteriors/` | Path files (`.staging/` empty after a completed run) |
| `bio_pars_report.pdf` | Estimated parameters (`main_cpp.R`) |
| `tag_<tag>_realism_*_rcpp.csv` | Realism report (`main_cpp.R`, `RUN_REALISM_REPORT`) |
| `simulated_all_transition_cost_*` | Sensitivity sweeps (`main_cpp.R`, `SENSITIVITY_MODE`) |
| `k_sweep_join_vs_split_demo_rcpp.pdf` | k-sweep plot (`main_cpp.R`, `RUN_K_SWEEP_DEMO`) |

The folder name is built from `BATCH_TS` and the settings
(`build_out_dir_name()`); `OUT_DIR_OVERRIDE` writes into a given folder (e.g.
to resume a chunked run).

## What it does

### The problem

A census records, per tree, the stems measured that year: `Tag`, `CensusID`,
`DBH` (missing when a stem was not measured), species, and from the anchor
census on a known stem identity (`TrueStemID`). The task is to assign the
earlier measurements to stems so that:

1. DBH trajectories are biologically plausible;
2. the assignment is globally optimal over all censuses of the tree;
3. life-cycle constraints hold: a stem is born once and does not come back
   after it dies.

It is hard because several stems per tree create a combinatorial assignment
space, unmeasured stems must be explained by mortality and recruitment,
measurement error blurs growth, and matching census by census can lock in
swaps that a global view would avoid.

### Overview of a run

For each tag:

1. The driver estimates the biological parameters per species and attaches
   them as `Bio_*` columns; the BCI driver also sets the pins (`TrueStemID`)
   from field tags and database IDs.
2. Rows of censuses in which a stem was alive but not measured are set aside
   (missing from the field).
3. The anchor census is chosen and the DP builds K tracks.
4. If the state space fits within `max_states`, the DP solves the censuses up
   to the anchor: MAP reconstruction, marginals and samples. Otherwise, or
   when the tag is routed there, the probabilistic matcher solves them.
5. The pin sweep sets pinned rows to their pin; rows after the anchor keep
   their pin or stay unassigned.
6. The post-engine helpers fix rows the pin sweep left behind, attach
   unmeasured terminal records, apply the broken-below rules, number the stems
   1..N and write the posterior path files.

### Track-based state space

The DP uses K latent tracks (identity slots):

- each measured observation of a census goes to exactly one track, and a
  track holds at most one observation per census (injective assignment);
- a track can be empty (`NA`) in any census;
- occupying the same track in successive censuses is identity continuity.

A **state** of census t is the injective assignment vector of its
$n_{\text{obs}}(t)$ observations to tracks. The number of states is

$$P(K, n_{\text{obs}}) = \frac{K!}{(K-n_{\text{obs}})!} = K (K-1) \cdots (K-n_{\text{obs}}+1)$$

The tracks of the anchor census carry the pins of its stems, so the
reconstruction is expressed in anchor identities: an earlier observation on an
anchor track belongs to that stem; an observation on another track belongs to
a stem that ends before the anchor.

### Life-cycle phases

Each track has a phase:

- **0 (before birth)**: not yet occupied; can move to phase 1;
- **1 (alive)**: occupied; stays alive or moves to phase 2;
- **2 (dead or unused)**: empty from here on (its stem died, or the track is
  never used); cannot return to phase 1.

Allowed changes are 0 → 0, 0 → 1 (recruitment), 1 → 1 (growth), 1 → 2
(death) and 2 → 2, so `OBS → NA → OBS` on one track is forbidden (no
resurrection). A stem that was alive but not measured in a census is handled
by the missing-from-field rule, which removes that census from the DP. The DP
key of a state is `"<assignment>|<phases>"`: the track of each observation,
comma-separated, and one phase digit per track. At the anchor, tracks without
an anchor stem are in phase 2 (dead or unused), which lets the DP reconstruct
stems that died before the anchor.

### Choosing K

K must hold the stems at the anchor, the largest stem count of any census and
the new stems needed to explain increases in stem count:

```text
births_needed = Σ_t max(0, n_obs(t+1) − n_obs(t))
K_base = max(#anchor stems, max_t n_obs(t), n_obs(first census) + births_needed)
K      = min(K_base + slack_tracks + main-rule breaks, max_tracks)
K      = max(K, #anchor stems + #pins before the anchor that are not at the anchor)
```

- **Slack tracks** let a stem die and a new one appear in the same interval
  when the stem count does not change. With
  `slack_require_anchor_recruitable = TRUE` they are granted only when an
  anchor stem is no larger than `Bio_Recruit_MaxDBH_unit` (plus
  `slack_require_anchor_eps`).
- **Main-rule breaks**: with the main-stem rule (section "ForestGEO field
  codes") each census with a break after a census with a main-stem record
  adds one track, because the broken record may be a new stem. Resprout
  observations under rule R1 need no track of their own.
- **Pins before the anchor**: every pinned stem that ends before the anchor
  needs its own track, or its pin could not be honoured.

When K is smaller than the largest stem count of a census, the tag goes to the
matcher (`K_too_small`).

**Interchangeable tracks and the label-free posterior.** A track with no stem
at the anchor and no pinned observation carries no identity. A reconstruction
with m stems on S such tracks therefore occurs under S!/(S−m)! labelings, all
with the same cost:

- the MAP compares single labelings, so the exported reconstruction is not
  affected;
- summing over labelings weights each reconstruction by its number of
  labelings, so the marginals and samples depend on the number of spare
  tracks (slack tracks, tracks added for breaks) and favour reconstructions
  with more stems ending before the anchor;
- a row on such a track has its marginal split across the equivalent labels
  (about 0.5 with two spare tracks, 0.33 with three).

With option `dp_global.posterior_label_free = TRUE` (BCI driver
`DP_POSTERIOR_LABEL_FREE`), each step back from the anchor that puts $a$ more
interchangeable tracks into use, while $u$ are already in use, is divided by
$(S-u)!/(S-u-a)!$. Each path is then weighted by $1/P(S, m)$, so each
reconstruction counts once in the marginals and samples, and the rows on
interchangeable tracks get a label-free `DP_PosteriorReconstructedProb`
(section "Uncertainty quantification").

### Understanding `max_states`

`max_states` (driver setting `DP_MAX_STATES`) is the one parameter that
decides when the DP hands a tag to the probabilistic matcher. It sets two
limits:

1. **Per census** (`enum_exceeded`): a census with more allowed states than
   `max_states`. When pins or growth bounds restrict the tracks of a census,
   `enumerate_states_constrained()` counts only the allowed assignments and
   returns `NULL` past `max_states` (never a truncated set).
2. **Per census pair** (`edge_count_exceeded`): the states of a census × the
   reachable states of the next census reach `max_states²` (`max_edges`).

The per-census limit binds first, because state counts grow factorially.
Fallback by stem count n, with every track allowed:

With K = n + 1 (`slack_tracks = 1`, no births needed):

| Stems (n) | Tracks (K) | States P(K, n) | 1,000 | 20,000 | 40,000 |
| :-: | :-: | --: | :-: | :-: | :-: |
| 2 | 3 | 6 | DP | DP | DP |
| 3 | 4 | 24 | DP | DP | DP |
| 4 | 5 | 120 | DP | DP | DP |
| 5 | 6 | 720 | DP | DP | DP |
| 6 | 7 | 5,040 | fallback | DP | DP |
| 7 | 8 | 40,320 | fallback | fallback | fallback |
| 8 | 9 | 362,880 | fallback | fallback | fallback |

With K = n + 2 (`slack_tracks = 1` plus one birth):

| Stems (n) | Tracks (K) | States P(K, n) | 1,000 | 20,000 | 40,000 |
| :-: | :-: | --: | :-: | :-: | :-: |
| 2 | 4 | 12 | DP | DP | DP |
| 3 | 5 | 60 | DP | DP | DP |
| 4 | 6 | 360 | DP | DP | DP |
| 5 | 7 | 2,520 | fallback | DP | DP |
| 6 | 8 | 20,160 | fallback | fallback | DP |
| 7 | 9 | 181,440 | fallback | fallback | fallback |
| 8 | 10 | 1,814,400 | fallback | fallback | fallback |

Largest stem count per census solved exactly:

| `max_states` | `max_edges` | K = n + 1 | K = n + 2 |
| --: | --: | :-: | :-: |
| 1,000 | 1,000,000 | 5 | 4 |
| 10,000 | 100,000,000 | 6 | 5 |
| 20,000 | 400,000,000 | 6 | 5 |
| 40,000 | 1,600,000,000 | 6 | 6 |

Restrictions from pins and pruning often leave far fewer states than P(K, n),
so tags above these counts can still be solved. Defaults: 50,000 in the DP
function, 40,000 in `main_cpp.R` and `main_cpp_chunk.R`, 10,000 in the BCI
driver.

Choosing a value:

```r
library(data.table)
dt <- fread("your_data.csv")
max_obs <- dt[!is.na(DBH), .N, by = .(Tag, CensusID)][, .(max_n = max(N)), by = Tag][order(-max_n)]
head(max_obs, 10)              # the most complex tags
n <- 6; K <- 8                 # stems in one census; n + 2 tracks
prod(K:(K - n + 1))            # P(8, 6) = 20,160: set DP_MAX_STATES above it
```

Higher values solve more tags exactly (slower, more memory); lower values send
more tags to the matcher (faster, approximate; same `Bio_*` parameters and
prune bounds). `estimate_dp_complexity()` gives the state counts of every tag
before a run (section "Checks and diagnostics").

### DP objective

The DP minimises the total cost of the transitions up to the anchor A:

$$\text{TotalCost} = \sum_{t=1}^{A-1} \text{TransitionCost}(s_t \to s_{t+1})$$

subject to injective assignments in every census, valid phase changes on
every track, and the anchor state fixed by the pins.

### Biological model

**Growth** (size-dependent mean, heteroskedastic variance). The annual growth
g = (D₁ − D₀)/Δt is normal with

$$\mu(D) = \alpha + \gamma \log D, \qquad \sigma(D) = \sigma_0 + \sigma_1 D$$

(α = `Bio_Mu_Growth`, γ = `Bio_Gamma_Growth`; γ = 0 gives a constant mean;
σ₀ = `Bio_Sigma0_Growth`, σ₁ = `Bio_Sigma1_Growth`).

**Mortality.** Hazard and probability of death over an interval:

$$h(D) = h_0 \exp(\beta D), \qquad P_{\text{death}}(D, \Delta t) = 1 - \exp\left(-h(D)\,\Delta t\right)$$

(h₀ = `Bio_H0_Mortality`, β = `Bio_Beta_Mortality`).

**Recruitment.** A new stem on an empty track, and its size:

$$P_{\text{recruit}}(\Delta t) = 1 - \exp(-\lambda\,\Delta t), \qquad D_{\text{recruit}} \sim \text{LogNormal}(\text{meanlog}, \text{sdlog})$$

(λ = `Bio_Recruitment_lambda`; `estimate_bio_pars()` estimates it per
established tree by default, section "Parameter estimation").

Probabilities are clamped to $[10^{-12}, 1 - 10^{-12}]$.

**Measurement error** (`use_measurement_error = TRUE`; Chave et al. 2004). A
recorded DBH carries remeasurement noise, $D_{\text{obs}} = D_{\text{true}} + \varepsilon$, with
a two-component mixture

$$\varepsilon \sim (1-p)\,\mathcal{N}\left(0, \mathrm{SD1}(D)^2\right) + p\,\mathcal{N}\left(0, \mathrm{SD2}^2\right), \qquad \mathrm{SD1}(D) = aD + b$$

(`meas_sd1_a` = a = 0.0062, `meas_sd1_b` = b = 0.0904 cm, `meas_sd2` = SD2 =
4.64 cm, `meas_p_big` = p = 0.05). The observed growth of a link then follows
a four-component mixture (small or large error at each census), evaluated with
log-sum-exp (case D below).

**DBH recorded in classes (rounded down).** In BCI 1982 and 1985 (censuses 1
and 2) saplings were measured in 5 mm increments, rounded down (CTFS R package
growth tutorial: "dbhs<50 mm are rounded down to 5-mm, necessary because
saplings at BCI in 1982 and 1985 were measuring in 5-mm increments"; Piponiot
et al. 2024, Appendix S1: DBH < 5.5 cm). The data agree: 100% (1982) and 99.9%
(1985) of DBH values below 50 mm are multiples of 5 mm, against about 20%
(chance) from 1990 on; between 50 and 55 mm, 93% and 73%, mostly a recorded
50 mm (a floored 50–54.9 mm). With a narrow growth distribution (σ₀ ≈ 0.04
cm/yr for many species) a 5 mm class step read as growth is a 4σ+ event, so
without the rule below one stem would be split into a death and a recruit.

When a census is flagged (`dbh_round_censuses`; BCI driver
`DBH_ROUND_CENSUSES = "1,2"`), a DBH $d$ below `dbh_round_max` (5.5 cm) at
that census is read as a true size in $[d, d + w)$, $w$ = `dbh_round_width`
(0.5 cm). For a link with rounding flags $r_0, r_1 \in \{0, 1\}$:

$$g_{\text{mid}} = g + \frac{(r_1 - r_0)\, w/2}{\Delta t}, \qquad \sigma_{\text{eff}}^2 = \sigma(D_0 + r_0 w/2)^2 + \frac{(r_0 + r_1)\, w^2/12}{\Delta t^2}$$

and the growth likelihood (step D5, and every component of the
measurement-error mixture) uses $g_{\text{mid}}$, $\mu(D_0 + r_0 w/2)$ and
$\sigma_{\text{eff}}$. The hard growth bounds (D1, D2) stay on the recorded
growth $g$: rounding changes how likely a link is, never which links are
allowed. Censuses that are not flagged, and stems at or above
`dbh_round_max`, are not affected. The matcher applies the same rule
(`compute_pairwise_log_likelihood()`); the BCI driver checks at run time that
every flagged census shows the rounding in the data.

### Transition cost

The transition cost from census t to t+1 is the negative log-likelihood,
summed over the K tracks, plus a tie-break. It is computed in C++
(`transition_cost_paired_rcpp_cpp()` and
`transition_cost_tracks_bio_batch_rcpp_cpp()` in `src/transition_cost_rcpp.cpp`;
the DP calls the paired version once per census pair for all candidate
pairs). `transition_cost_tracks_bio_components()` (`R/dp_global_bio.R`) gives
the same cost split into its terms.

Each track contributes independently. Inputs: the DBH of each track at t and
t+1 (`NA` when empty), Δt (years), the `Bio_*` parameters, `eps_tiebreak`
(10⁻⁶) and `hard_penalty` (10⁶). The recruitment probability
$P_{\text{recruit}}(\Delta t)$ is computed once per transition and shared by
all tracks.

**Case A: NA → NA** (the track stays empty):

$$\text{cost} \mathrel{+}= -\log(1 - P_{\text{recruit}})$$

A high recruitment rate makes long empty histories costly.

**Case B: NA → DBH** (recruitment). When $D_1$ > `Bio_Recruit_MaxDBH_unit`
or $D_1 \le 0$: `cost += hard_penalty`. Otherwise

$$\text{cost} \mathrel{+}= -\log P_{\text{recruit}} - \log f_{\text{LogNormal}}(D_1;\ \text{meanlog}, \text{sdlog}), \qquad f_{\text{LogNormal}}(x) = \frac{1}{x\,\text{sdlog}\sqrt{2\pi}} \exp\left(-\frac{(\log x - \text{meanlog})^2}{2\,\text{sdlog}^2}\right)$$

**Case C: DBH → NA** (death):

$$\text{cost} \mathrel{+}= -\log P_{\text{death}}(D_0, \Delta t)$$

A disappearance is cheap only where mortality is high.

**Case D: DBH → DBH** (growth), with $g = (D_1 - D_0)/\Delta t$:

- **D1** hard shrinkage: $g$ < `Bio_Max_Shrink` → `cost += hard_penalty`; next track.
- **D2** hard growth: $g$ > `Bio_Max_Growth` → `cost += hard_penalty`; next track.
- **D3** growth SD: $\sigma(D_0) = \max(\sigma_0 + \sigma_1 D_0,\ 10^{-6})$.
- **D4** mean growth: $\mu(D_0) = \mu_{\text{const}} + \mu_\gamma \log D_0$ when
  $\mu_\gamma \ne 0$ and $D_0 > 0$, else $\mu_{\text{const}}$.
- **D5** growth likelihood. Without measurement error:

  $$\text{cost} \mathrel{+}= \frac{(g - \mu(D_0))^2}{2\sigma(D_0)^2} + \log \sigma(D_0) + \tfrac{1}{2}\log(2\pi)$$

  With measurement error, a four-component normal mixture. Component $i$ has
  the error SDs $(\text{SD}_a, \text{SD}_b)$ at the two censuses:
  (SD1($D_0$), SD1($D_1$)), (SD1($D_0$), SD2), (SD2, SD1($D_1$)), (SD2, SD2),
  with weights $w_i \in \{(1-p)^2, (1-p)p, p(1-p), p^2\}$, and

  $$\text{SD}_{\text{meas},i} = \frac{\sqrt{\text{SD}_a^2 + \text{SD}_b^2}}{\Delta t}, \qquad \text{SD}_{\text{tot},i} = \sqrt{\sigma(D_0)^2 + \text{SD}_{\text{meas},i}^2}, \qquad \ell_i = \log w_i + \log \phi\left(g;\ \mu(D_0), \text{SD}_{\text{tot},i}\right)$$

  $$\text{cost} \mathrel{+}= -\left(m + \log \sum_{i=1}^{4} \exp(\ell_i - m)\right), \qquad m = \max_i \ell_i$$

  (at a rounded census: $g_{\text{mid}}$, $D_0 + r_0 w/2$ and $\sigma_{\text{eff}}$).
- **D6** soft shrinkage: when $D_1 < D_0$,
  $\text{cost} \mathrel{+}= k_{\text{shrink}} (D_0 - D_1)^2$ (cm; $k$ in 1/cm²).
- **D7** soft growth: with $D_{1,\text{cap}} = D_0 + g_{\text{soft}}\,\Delta t$
  ($g_{\text{soft}}$ = `Bio_Max_Growth_Soft`), when $D_1 > D_{1,\text{cap}}$,
  $\text{cost} \mathrel{+}= k_{\text{growth}} (D_1 - D_{1,\text{cap}})^2$.
- **D8** tie-break (not biological), after summing the tracks: with $r_0$, $r_1$
  the DBH ranks of the measured tracks at t and t+1 (stable ranks, ties in track
  order),

  $$\text{cost} \mathrel{+}= \varepsilon_{\text{tiebreak}} \sum_{k \in \text{both measured}} |r_0[k] - r_1[k]|$$

  Among otherwise equal costs it favours the assignment that keeps the size
  order of the stems.

A recruitment is costed only on an empty track (phase 0 → 1) and a death only
on a track that leaves phase 1, so a stem that keeps growing is never charged a
death and a recruit.

### Anchor census

`anchor_start` is the census whose identities are trusted (BCI: census 7,
2010, the first census that tagged every stem). The DP solves the censuses up
to it.

1. A tag whose last measurement is before `anchor_start` is anchored at its
   last measured census.
2. A tag measured after the anchor is solved up to the anchor and its later
   rows are appended (section "Outputs"). When the anchor census has no
   measured stem, the anchor moves forward to the first later census with a
   measured stem, preferring one with a pin.
3. When the anchor census has no row with a DBH or a pin, the latest earlier
   census with a measured, pinned row becomes the anchor.
4. Provisional anchor (`allow_provisional_anchor = TRUE`): when no census has
   both, the last measured census becomes the anchor and its measured rows get
   new provisional IDs (`provisional_dp`); the same happens when the anchor
   census has measured rows without a pin.
5. Otherwise the tag goes to the matcher (`anchor_*`, `no_obs_up_to_anchor`).

### Backward DP recursion

Let $s_t$ be a state of census t, $c(s_t, s_{t+1})$ the transition cost and A
the anchor census, whose state is fixed by the pins. Going back from A, census
by census, the DP:

1. enumerates the allowed states of census t (restricted to the tracks that
   pins and growth bounds allow, below);
2. pairs each state with the reachable states of census t+1 and removes the
   pairs that break the phase rules, the code rules or the pruning bounds (C++
   `derive_phase_prev_batch_rcpp()`);
3. costs the remaining pairs in one call and keeps, for each state, the
   cheapest continuation (Viterbi value and back-pointer) and the log-sum-exp
   of all continuations:

$$V_t(s) = \min_{s'} \left[ c(s, s') + V_{t+1}(s') \right], \qquad \log Z_t(s) = \log \sum_{s'} \exp\left( -\frac{c(s, s')}{\tau} + \log Z_{t+1}(s') \right)$$

with $V_A = 0$ and $\log Z_A = 0$ at the anchor state.

The MAP reconstruction follows the back-pointers from the cheapest start
state to the anchor (among equal costs the first enumerated state wins).

### Uncertainty quantification

The DP defines a Gibbs posterior over full reconstructions:

$$P(\text{path}) \propto \exp\left(-\frac{\text{TotalCost}(\text{path})}{\tau}\right)$$

with temperature τ (`temperature`; 1 in the drivers: the costs are negative
log-likelihoods, so τ = 1 is the Bayesian posterior; τ < 1 sharpens it towards
the MAP, τ > 1 flattens it).

- **Marginals**: the backward pass gives the log partition function; a forward
  pass normalises the state weights into the marginal probability of every
  observation on every track (anchor tracks and spare tracks).
- **Samples** (`posterior_samples` > 0): full reconstructions drawn by backward
  sampling from the same weights (seed `posterior_sample_seed`), each with its
  `logp`, staged for the path file.

Summaries per observation:

- `DP_PosteriorTop{k}ID` / `Prob`: the most likely tracks.
- `DP_PosteriorEntropy`: $H = -\sum_i p_i \log p_i$ (higher = more uncertain).
- `DP_PosteriorReconstructedProb`: posterior probability that the row belongs
  to its exported stem. For a stem on an anchor or pinned track, the marginal
  of that track (the anchor or pinned observation fixes the label). For a
  stem on an interchangeable track (it ends before the anchor and holds no
  pin): with the label-free posterior, the share of samples in which the row's
  stem ends at the same observation as the exported stem (`NA` without
  samples); otherwise the track marginal, split across equivalent labels.
- `DP_PosteriorUnlinkedProb`: posterior mass on the tracks without an anchor
  stem.

Top-k, entropy and the unlinked probability stay per-track quantities.

`add_dp_posterior_bins()` labels each row:

- **unlinked-likely**: unlinked probability ≥ `unlinked_prob` (0.5); takes
  precedence over confident;
- **confident**: probability of the reconstructed ID ≥ `confident_prob`
  (0.95);
- **ambiguous**: everything else with a posterior.

### Pruning and conservative guards

The state space grows factorially with the stem count, and the candidate
transitions between two censuses grow with the product of their state counts.
Cheap, conservative pre-filters remove biologically impossible transitions
before the expensive cost is evaluated, while keeping borderline cases for the
full cost.

**Where pruning runs.** In the backward recursion,
`derive_phase_prev_batch_rcpp()` (C++) tests every pair of a state at t and a
reachable state at t+1: first the phase rules and the R-recruit constraint,
then the hard growth and recruit-size bounds. Only the pairs that pass are
costed. The growth and recruit-size bounds apply only with `prune_hard = TRUE`
(default); without it the phase rules and the code rules still apply.

The prune bounds also restrict the state enumeration
(`enumerate_states_constrained()`): a track holding an anchor stem is open to
an observation only if an observation of the next census that may use that
track is growth-compatible with it (twice the non-taper bounds for those
forms, leaving room for the HOM widening). Tracks empty at the anchor are open
to every observation, and a pinned observation is restricted to its pin's
track. A census can therefore have far fewer states than P(K, n).

**Effective bounds**, built in layers (all in cm/yr):

```text
user_min = prune_min_growth  (NULL: min_growth)
user_max = prune_max_growth  (NULL: max_growth)

if prune_use_bio_bounds:                 # 2. biological layer
    eff_min = max(user_min, Bio_Max_Shrink)
    eff_max = min(user_max, Bio_Max_Growth)
else:
    eff_min = user_min; eff_max = user_max

if growth_form in non_taper_corrected_growth_forms:    # 3. non-taper override
    eff_min = non_taper_corrected_prune_min_growth
    eff_max = non_taper_corrected_prune_max_growth

per census pair, for those forms, with a hom column and hom_tolerance_scale > 0:   # 4. HOM widening
    hom_tol      = hom_tolerance_scale × max(|HOM − 1.3|) / interval_years
    eff_min_pair = eff_min − hom_tol
    eff_max_pair = eff_max + hom_tol
```

A link with $g \notin$ [eff_min, eff_max] is pruned; a recruit (NA → DBH)
above `eff_recruit_max` is pruned, where `eff_recruit_max` is
`prune_recruit_max_dbh` (default `Bio_Recruit_MaxDBH_unit`), capped at
`Bio_Recruit_MaxDBH_unit` when `prune_use_bio_recruit = TRUE` and that value
is finite. A transition inside the prune bounds but outside the `Bio_*` hard
bounds is not pruned and costs `hard_penalty`. A census pair whose interval is
not finite or not positive is not pruned.

**Non-taper-corrected growth forms** (default `"palm"`, `"strangler_fig"`,
`"tree_fern"`; a tag uses its most common `growth_form`, exact match, with a
warning when it has several) get their own bounds, for two reasons: they grow
in diameter (palms can add 1–3 cm/yr; strangler figs change diameter as they
encircle or replace their host), and they cannot be taper-corrected, so a
change in the height of measurement changes the recorded DBH. Their bounds
replace the general ones (defaults −0.625 and 6.25 cm/yr; the drivers pass
`PRUNE_BOUND_FACTOR` × the fixed limits). HOM widening then adds
`hom_tolerance_scale` (2 cm of DBH per m of deviation from 1.3 m) per census
pair, using the worst deviation over the tag's rows in the two censuses (`NA`
HOM counts as 1.3 m; 0 turns it off). The values must be `growth_form` labels
of the data: the BCI driver assigns `palm`, `strangler` and `fern` and checks
at run time that every value it lists occurs in the data.

**Driver configuration.** The drivers pass `prune_min_growth` /
`prune_max_growth` = `PRUNE_BOUND_FACTOR` × `MAX_SHRINK_FIXED` /
`MAX_GROWTH_FIXED` with `prune_use_bio_bounds = FALSE`, the same two values as
non-taper bounds, and `prune_recruit_max_dbh` = `PRUNE_BOUND_FACTOR` ×
`RECRUIT_MAX_FIXED` with `prune_use_bio_recruit = FALSE`. Every tag therefore
gets the same base window: [−0.625, 9.375] in `main_cpp.R` (factor 1.25,
`MAX_GROWTH_FIXED` = 7.5), [−0.625, 6.25] in `main_cpp_chunk.R` (1.25, 5) and
[−2.5, 25] in the BCI drivers (5, 5). What differs for non-taper forms is the
HOM widening; the override matters when the general bounds are tightened
(e.g. `prune_use_bio_bounds = TRUE`).

Practical guidance:

- with the function defaults (`prune_*` `NULL`), `min_growth` / `max_growth`
  are intersected with the `Bio_*` bounds;
- for wider pruning than the biological bounds, set
  `prune_use_bio_bounds = FALSE` and wide `prune_min_growth` /
  `prune_max_growth` (e.g. −10 to 25);
- for stricter recruit-size pruning, give a small `prune_recruit_max_dbh` with
  `prune_use_bio_recruit = FALSE`;
- pruning that removes every feasible transition of a census sends the tag to
  the matcher; frequent fallbacks on reasonable data call for wider bounds, or
  `prune_hard = FALSE` in a diagnostic run.

`attr(out, "DP_PruneInfo")` records `total_examined`, `total_pruned` (every
transition removed before costing: phase-infeasible, pruned, or removed by the
code rules), `per_census` (removed transitions per census pair, named by the
earlier census), `eff_min_growth`, `eff_max_growth`, `eff_recruit_max` (before
HOM widening), `is_non_taper_corrected` and `use_hom_relax`. The counts are
filled only with `prune_hard = TRUE` on tags the DP solves; the effective
bounds are also logged at the start of the backward pass.

### ForestGEO field codes

The engines read two columns of field-recorded events, the code list
`ListOfTSM` and `Status`. They affect transition feasibility, segment
splitting and the handling of unmeasured censuses.

**Codes.** `R`, `RP`, `RF`, `RT`, `QR` and `OR`, matched as whole codes in
`ListOfTSM` (`grepl()` with `perl = TRUE`, pattern `\b(R|RP|RF|RT|QR|OR)\b`:
`"R;B"` matches, `"NORMAL"` does not), and `Status == "broken below"`, mark a
break. A **resprout observation** is such a row with a DBH; a **stump** is a
broken-below row without a DBH. The missing-from-field pattern is `\bMF\b`;
the labelling of rows after the anchor uses the code pattern without `QR`.
When a tag has a resprout observation before the anchor, the mechanisms below
apply.

#### 1. Resprout segment split

Argument `allow_segment_split` of the DP function (default `TRUE`; the BCI
driver sets `DP_ALLOW_SEGMENT_SPLIT = FALSE`: before 2010 a break is written on
the tree's main-stem record only, and the tree's other stems keep growing, so
ending every stem at the break would cut real links). With the split, the DP
solves two independent problems at the first census before the anchor (not the
tag's first census) that holds a resprout observation:

- **pre-segment**: from the first census up to and including the R census,
  anchored at its last measured census (provisional IDs where it has no pin).
  In the DP the R-coded measurement is the last record of the stem that broke;
- **post-segment**: from the census after the R census to the anchor. Every
  stem of its first census is a new stem, independent of the stems before.

R means the tree resprouted, so the stems after that census have no identity
continuity with those before it. The exported identity of the R-coded row is
settled after the engine (pin sweep; rule R1). Pre-segment IDs are shifted
above the post-segment IDs.

Posterior samples: each sub-call draws its own samples and hands them back
(`posterior_return_samples = TRUE`, or `return_samples` in the matcher)
instead of staging them. Because every stem ends at the split, the two
segments are independent, and draw k of one segment paired with draw k of the
other is an exact draw of the whole tag (`pair_segment_posterior_samples()`).
The pre-segment uses `posterior_sample_seed + 1`, a segment with fewer draws
is recycled, a segment that drew none (e.g. one census) contributes its
exported grouping (`segment_fixed_assignment()`), and pre-segment labels are
shifted above post-segment labels. The paired samples are staged once per tag
(`stage_posterior_samples()`). Sub-calls set `allow_segment_split = FALSE`, so
a later resprout observation inside a segment follows the R-recruit
constraint (3) or the main-stem rule (5); census IDs absent from a sub-call's
data are filtered out.

#### 2. Post-segment recruit continuity

In the post-segment, every stem of the first census is a known recruit. When
that census has no more stems than the next, every one of them must be on a
track that is also occupied at the next census: a freshly resprouted stem
dying at once while an independent recruit appears on another track is far
less parsimonious than all recruits continuing. If the constraint would remove
every candidate transition of the census, it is skipped there (with a
message), so the DP does not dead-end.

#### 3. R-recruit constraint

Only with rule R1 (option `dp_global.bb_dbh_rule = "r1"`, the engine
default). A resprout observation at census t+1 must continue a track occupied
at census t: in the DP it is the last record of the stem that broke, so it
cannot sit on a not-yet-recruited track. It needs no track of its own, and K
is not raised for it. The constraint is enforced for every candidate
transition in `derive_phase_prev_batch_rcpp()`. The backward loop repeats the
test as an R filter, skipped when census t has no measured stem; the C++ test
has no such exception, so then no transition is feasible and the tag (or
segment) goes to the matcher (`no_states_produced`).

#### 4. NA-R barrier

A census before the anchor with no measured stem and at least one unmeasured
row with a break code is a barrier: the whole tree died back, so no identity
crosses it. After decoding, a track with rows on both sides is cut and its
unpinned rows before the barrier take a new ID. A track is left whole when its
ID is a `TrueStemID` known after the barrier and no row before the barrier
carries another `TrueStemID`; pinned rows are never relabelled. The unmeasured
break rows of the barrier census take the stems of the census before, in row
order. `do_fallback()` applies the same cut to the matcher's output (every
crossing track; the pin sweep then restores the pinned rows).

#### 5. Main-stem rule (option `dp_global.bb_dbh_rule = "main_ends"`)

BCI field convention up to 2010: a break is written on the tree's main-stem
record (`PrimaryStem == "main"`; in the census before there is exactly one
main record in 99.6–100% of the cases, always the largest stem). The record
carries the DBH of the largest stem still alive when any stem of the tree is
alive (in a single-stem tree the resprout; in a multiple-stem tree usually the
former second stem), and no DBH when nothing is alive. The rule follows the
record, with no size threshold:

- the main stem of the census before ends at the break;
- the broken record with a DBH is matched like any other measurement, but it
  cannot continue that main stem. With one stem before the break it is always
  a new stem; with several it continues another stem when one fits.

**DP**: at a transition into a census with a resprout observation, the track
of the main-stem observation of the census before may hold no observation of
the next census (main-stem filter in the backward loop); the resprout flags are
not passed to `derive_phase_prev_batch_rcpp()` (no forced continuation), the
R-recruit filter is skipped, and K gets one more track per such break.
**Matcher**: in that census pair the main stem's link cells are set to −Inf
(`apply_main_end_mask()`; only its death cell stays open), in the same masking
step as the pins (`pin_masked_pair()`), and `split_live_resprout_tracks()` is
not applied. **After the engine**: R1 only checks the rule, and
`apply_pin_track_rejoin()` / `apply_pins_to_samples()` cut a track at a
broken-below row with a DBH only where it would continue the main stem
(`bb_main_violation()`). The rule acts on the same records in the export and in
every posterior sample. From 2015 BCI broken-below records carry no DBH, so
the rule has no effect after 2010.

#### 6. Missing from the field

Rows of stems that were alive but not measured are set aside before the engine
and re-inserted afterwards, so the gap is not read as a death followed by a
recruitment:

- **explicit**: a row without a DBH whose `ListOfTSM` holds `MF` or whose
  `Status` is `"missing"`; the episode extends through the following
  consecutive censuses in which no row of the tag has a DBH;
- **implicit**: a census in which no row of the tag has a DBH, between the
  tag's first and last measured census, unless a row of it carries a break
  code (a break is an event, not a field absence).

A census left without rows is dropped from the DP, so tracks link the censuses
on both sides directly. After the engine each set-aside row takes the one
track assigned in the nearest censuses before and after and free at its census
(`dp_mf_inferred`); with several candidates it stays `NA`.

### Pins, the pin sweep and hard pins

A pin (`TrueStemID`) is a known stem identity: at the anchor census, and at
any earlier row the driver pins (BCI: field tags, database IDs from 2010 and
database IDs propagated to earlier records; the pin steps are in
`BCI_stem_reconstruction/2_STEM_IDENTIFICATION/README.md`).

**While solving** (`pin_truestemid = TRUE`, default): the DP enumerates only
the states that put each pinned observation on its pin's track; the matcher
constrains every sampled pair (section "Probabilistic matcher").

**Pin sweep**: every row with a pin (not a provisional anchor) gets
`ReconstructedStemID = TrueStemID` and `ReconstructionMethod = "given"`. The
engine's value goes to `ReconstructedStemID_PreSweep` and changed rows are
flagged in `SweepAuditOverride`. The sweep runs inside the DP's fallback paths
(`finalize_out()`), in the matcher and in the drivers' `run_dp_one_group()`,
which is the only sweep on the normal DP path (the DP honours the pins through
its state enumeration). The drivers' sweep processes rows in order and keeps
the engine's stem when the pin would put two rows of one census on one stem
(`SweepRollbackToPreSweep`). Unpinned rows keep the engine's stem;
`apply_pin_track_rejoin()` handles unpinned rows the sweep left behind.

**Hard pins** (option `dp_global.hard_pins = TRUE`, BCI driver `HARD_PINS`):
pins also hold after the engine.

- `apply_broken_below_invariants()` (and its per-sample version) never cuts
  between two rows with the same pin. When the rows it would move hold the
  current ID as their pin, the other rows of the stem take the new ID instead,
  so no `TrueStemID` is rewritten.
- `apply_terminal_to_host()` leaves pinned records where they are.
- `renumber_engine_minted_ids()` renumbers `TrueStemID` only where the engine
  minted it (`provisional_dp`); database pins are kept.
- **Field tags**: the `StemTag`s from the census of option
  `dp_global.field_tag_first_census` on (`field_tag_rows()`; default `Inf`,
  none; BCI 7, the first census that tagged every stem). A broken-below record
  without a DBH that carries a field tag is that stem with no measurable DBH
  at 1.3 m, so it stays with its tag's stem, also before the tag's first
  measurement (a 2010 stump whose tag is measured in 2015) and between two
  measurements (a stem measured, broken below and measured again under its
  tag). Field-tagged stems are pinned before the engine, so each stays one
  stem in the export and in every posterior sample.
- **Exception** (`pin_prelife_rows()`), so that a stump stays the end record
  of the stem that broke: an unmeasured dead, stem-dead, broken-below or
  missing record that comes before any evidence of life of its own pin is not
  bound by the pin when it carries no field tag, or when its pin never lives
  (all its records then go to the stem that broke together). Before 2010 the
  BCI database stores a break under a new StemID, the one of the resprout
  measured later, and shows that StemID's later tag on its earlier records;
  those are not field tags. `apply_terminal_to_host()` moves such a record to
  the stem that broke, and R2 still cuts the resprout from it.

### Fallback mechanisms

**Routing.** Whole groups of tags can bypass the DP:

- growth form (`fallback_growth_forms`, driver `DP_FALLBACK_GROWTH_FORMS`): a
  tag whose rows up to the anchor hold a listed `growth_form` goes to the
  matcher (`growth_form_forced`);
- species (`prob_species`, driver `PROB_SPECIES`): a tag whose `Species`
  holds a listed value (`species_forced_probabilistic`).

Both take a character vector or one comma- or semicolon-separated string.

The DP also falls back when no DBH is observed up to the anchor or no usable
anchor exists, when a census exceeds `max_states` (`enum_exceeded`), when a
census pair reaches `max_states²` (`edge_count_exceeded`), when K is too small
(`K_too_small`), when the recursion, the decoding or the forward pass finds no
feasible state or transition, and when the DP raises an error (drivers).

**`do_fallback()`** sends every reason to `match_stems_probabilistic()`, with
the DP's growth bounds (`min_growth`, `max_growth`), prune bounds
(`prune_min_growth`, `prune_max_growth`, `prune_recruit_max_dbh`; a prune
bound replaces the corresponding base bound, and `prune_recruit_max_dbh`
replaces `Bio_Recruit_MaxDBH_unit`) and the `use_bio_hard_*` flags. After the
matcher it applies:

- the **NA-R barrier** (every crossing track);
- with rule R1, **R-boundary splitting** (`split_live_resprout_tracks()`): a
  measured broken-below or resprout row starts a new stem, so only that stem's
  track is split and its rows before the break take a new ID. The break row
  and its continuation keep the matcher's ID, the anchor `TrueStemID` that the
  sweep restores (renaming the rows after the break instead would let the sweep
  re-join the anchor to the stem that broke, leaving one-census stems counted
  twice in the stock). The tree's other stems keep their links, as in the
  matcher's samples, and a track with a pinned row before the break is not
  split. The anchor census is left to R1. Unlike the DP's segment split, which
  ends every stem at the first R census, this splits only the broken stem.
  With the main-stem rule it is not applied: the matcher already ends the main
  stem and matches the broken record freely.

**Error fallback.** The drivers' `run_dp_one_group()` wraps the DP call in
`tryCatch`. On an error it runs the matcher on the rows up to the anchor with
the same R-boundary splitting, labels the later rows (`given` /
`none_after_anchor`), logs the message and stores `error:<message>` (first
200 characters) in `DP_FallbackReason`.

**Segment split consistency.** When a tag is split at a resprout census, the
DP first evaluates the whole-tag state space (P(K, n) of every census with the
pre-split K, before slack tracks). If it exceeds `max_states`, both
sub-calls get `max_states = 0`, so both segments use the matcher; otherwise
each segment is solved independently with its own K and can use a different
engine (e.g. the post-segment by the DP, the pre-segment by the matcher), each
with its own `ReconstructionMethod`.

Reason codes (`DP_FallbackReason`; post-anchor rows of a fallback tag carry
the same value, several reasons are joined by `;`):

| Code | Meaning | What to do |
| --- | --- | --- |
| `growth_form_forced` | A `growth_form` of the tag is in `fallback_growth_forms` | Requested routing |
| `species_forced_probabilistic` | A `Species` value of the tag is in `prob_species` | Requested routing |
| `no_obs_up_to_anchor` | No DBH at or before the requested anchor | Check `anchor_start` |
| `anchor_missing_truestem` | No census has both a DBH and a pin, and no provisional anchor could be set | Give pins or allow provisional anchors |
| `anchor_missing_obs` | The anchor census has no DBH | Choose another `anchor_start` |
| `anchor_missing_truestem_prov_disabled` | The anchor has DBH without pins and provisional anchors are off | `allow_provisional_anchor = TRUE` or give pins |
| `anchor_ids_missing` | No anchor IDs after locating an anchor census | Inspect the data |
| `K_too_small` | K is smaller than the largest stem count of a census | Raise `max_tracks` or `slack_tracks` |
| `enum_exceeded` | A census has more allowed states than `max_states` | Raise `max_states` |
| `edge_count_exceeded` | A census pair reaches `max_states²` state pairs | Raise `DP_MAX_STATES` |
| `anchor_truestem_not_found` | Anchor pins could not be mapped to tracks | Check for unexpected or duplicated pins |
| `next_assign_row_mismatch` | Internal mismatch of assignment keys | Report with a minimal example |
| `no_reachable_next_states` | No reachable next state in the recursion | Relax pruning or check the data |
| `no_states_produced` | No state left at a census after pruning and constraints | Relax pruning or widen the growth bounds |
| `no_feasible_edges` | No feasible transition between states | Relax pruning and growth bounds |
| `decode_failure` | No valid MAP start state | Check pruning and parameters |
| `viterbi_decode_failure` | Invalid back-pointer in the Viterbi backtrace | Report |
| `assign_mismatch` | MAP assignment length differs from the observed rows | Report |
| `forward_edges_missing` | Missing transitions in the forward pass | Check the transitions and pruning |
| `forward_no_alpha` | No finite forward weight (underflow or fully pruned) | Relax pruning or change the temperature |
| `error:<message>` | The DP raised an error (drivers) | Read the message; rerun the tag alone |

```r
res <- match_stems_dp_global_backward_marginals_batch(dt, anchor_start = 5, ...)
unique(na.omit(res$DP_FallbackReason))
```

### Probabilistic matcher

`match_stems_probabilistic()` works census pair by census pair and stitches
the pairs into full reconstructions:

1. **Link scores** (`compute_pairwise_log_likelihood()`): for every stem of
   census t and every stem of t+1, from the `Bio_*` parameters of the DP: the
   Gaussian growth log-likelihood with size-dependent mean and variance, the
   log survival probability of the earlier stem, the soft shrinkage penalty,
   and a soft growth penalty above `Bio_Max_Growth`. The measurement-error
   mixture is not used. A link outside the growth bounds, or outside
   `Bio_Max_Shrink` / `Bio_Max_Growth` when `use_bio_hard_shrink_in_prob` /
   `use_bio_hard_growth_in_prob` are `TRUE`, is forbidden. DBH classes are read
   as in the DP.
2. **Assignment matrix** (`augment_cost_matrix()`), with the mortality hazard
   and recruitment size and rate of the DP. Two modes (`birth_death`, DP
   argument `prob_birth_death`, driver `PROB_BIRTH_DEATH`):
   - **birth-death** (default): K = n_curr + n_next. Each current stem has one
     death cell (`[i, n_next + i]`, log P_death) and each next stem one
     recruitment cell (`[n_curr + j, j]`, log P_recruit plus the log size
     density); cells pairing a recruitment row with a death column carry no
     event and score 0. Every combination of survivals, deaths and
     recruitments is an assignment scored by the sum of its events, so a link
     is taken only when it is more likely than the death of the earlier stem
     plus the recruitment of the later one, the choice the DP makes. The
     matrix carries attribute `"bd"`;
   - **count-based slots** (`FALSE`): K = max(n_curr, n_next), raised only
     until an assignment without forbidden links exists (with M the largest
     set of allowed survival links, K ≥ n_curr + n_next − M). Every allowed
     assignment keeps M survivals, so with a stable stem count no stem can die
     and none can be recruited.
3. **Stochastic assignment** (`greedy_assignment_gumbel()`): `n_samples` (200)
   assignments per pair. Gumbel(0, `temperature`) noise is added to the
   log-likelihoods. Birth-death mode: noise on every event cell (none on the
   empty cells), and each sample is the exact best assignment of the perturbed
   pair (`hungarian_min_rcpp()`, Kuhn–Munkres in C++; one cell per event, so
   deaths and recruitments are not favoured by having several equivalent
   slots). Count-based mode: rows assigned greedily to the best free column;
   a pair whose greedy assignment uses a forbidden link is solved again
   exactly (`enforce_feasible_assignment()`). Taking the best assignment of a
   Gumbel-perturbed matrix approximates sampling from the Gibbs distribution
   over assignments.
4. **Pins while sampling.** `apply_pin_mask()` sends an observation pinned to a
   stem present at the anchor to the column that carries that stem.
   `apply_track_pin_mask()` covers the pins of every stem, also of stems that
   end before the anchor: while sampling backward each next-census observation
   carries the pin of its track (its own, or one inherited from a later census
   of the same sample, `propagate_track_pin()`), and a pinned observation must
   join the observation that carries its pin and may not join one that
   carries another pin. Unpinned observations are free, so a stem's earlier
   rows under an older StemID can still be linked to it. **Pins never force a
   link** (`pin_masked_pair()`): with count-based slots a pair is sized for the
   growth limits only, so the masks can forbid every link its slots relied on
   (one stem at C6 pinned to A, one at C7 pinned to B: a 1×1 matrix whose only
   link the mask forbids). When the masked survival links need more death and
   recruitment slots than the pair has, the pair is rebuilt with
   K = n_curr + n_next − M (`augment_cost_matrix(K_min = …)`) and masked
   again, so one pinned stem ends and the other starts. With the main-stem
   rule the main stem's links are masked in the same step.
5. **Look-ahead** (`prob_lookahead_weight` > 0, at least three censuses and
   K ≥ 4): when stitching backward from the anchor, the matrix of pair
   (t, t+1) is conditioned on the already-resolved pair (t+1, t+2)
   (`condition_cost_matrix()`): each candidate column gets a continuity bonus
   from the log-likelihood that its stem continues plausibly to t+2,
   normalised (best column 0, others down to −2 log units) and weighted
   (function default 0.5; drivers 1). The K ≥ 4 gate avoids distorting simple
   tags.
6. **Stitching** the per-pair assignments backward from the anchor into full
   samples.
7. **Repair** (`repair_stitched_growth_violations()`), so the marginals count
   only biologically valid paths:
   - hard-rate check: a link with annual growth outside [min_rate, max_rate]
     is cut. Pinned observations are never cut off: if the earlier observation
     is pinned the later one is cut instead (unless it is pinned too or at the
     anchor), and a link between two pinned observations is kept;
   - cumulative-shrinkage check (`use_bio_hard_shrink_in_prob = TRUE` and
     finite `n_sigma_me`): a run of small decreases can add up to more
     shrinkage than measurement error explains. With SD(D) = 0.0062 D + 0.0904
     cm (the small-error component), a run whose cumulative shrinkage exceeds
     $n_\sigma \sqrt{\text{SD}(d_{\text{start}})^2 + \text{SD}(d_{\text{curr}})^2}$
     is cut at its start ($n_\sigma$ = `n_sigma_me`, default 3; the BCI driver
     sets `Inf`, off). The DP needs no such check: each step pays its growth
     likelihood, so a run of decreases accumulates cost.
8. **Pin filter** (`filter_pin_consistent_samples()`): samples in which an
   observation pinned to an anchor stem is not on that stem's track are
   dropped, so the marginals are conditioned on those pins. If fewer than
   min(10, `n_samples`/4) would remain, all are kept (soft pins) with a
   warning, e.g. when a stem measured at C3 and at the anchor but not in
   between cannot be bridged (the matcher has no alive-but-unmeasured state;
   the pin sweep joins the two rows in the export).
9. **Marginals** (`compute_marginals_from_samples()`): pass 1 counts the track
   of every observation across samples (Top-k, entropy). Pass 2 resolves
   conflicts greedily from the anchor outward, by confidence, taking the next
   alternative when an ID is taken or would break the growth rates against
   the nearest resolved census; its IDs are replaced by the export of step 10.
10. **Export** (`select_consensus_trajectory()`): the sample whose predecessor
    links agree most with the links of the other samples (maximum expected
    accuracy; ties broken by how often its whole partition occurs). With
    option `dp_global.prob_export_rule = "fewest_me"` the samples with the
    fewest links whose recorded DBH shrinks beyond measurement error (CTFS rule,
    $D_1 \le D_0 - 4\,(0.006214\,D_0 + 0.9036)$ in mm; `count_beyond_me_links()`)
    are considered first. The samples are unchanged and the export is always a
    sampled reconstruction.
11. **Diagnostic check**: residual hard-rate violations of the export are
    counted and logged; nothing is changed.
12. **Pin sweep** (`pin_truestemid = TRUE`) and **posterior export**
    (`export_probabilistic_posteriors()`): the samples are staged in the DP's
    format; `finalize_posterior_paths()` writes the path file later.

`DP_PosteriorReconstructedProb` is the share of samples that give the row its
exported stem (1 at the anchor); `DP_PosteriorUnlinkedProb` stays `NA`.
`ReconstructionMethod` is `probabilistic`, `given` (pinned rows) or
`provisional_dp`.

| Parameter | Default | Content |
| --- | --- | --- |
| `n_samples` | 200 | Gumbel samples (DP argument `prob_n_samples`) |
| `temperature` | 1 | Gumbel scale (shared with the DP) |
| `posterior_top_k` | 2 | Top-k alternatives reported |
| `prob_lookahead_weight` | 0.5 | Look-ahead weight (0 = off) |
| `n_sigma_me` | 3 | Cumulative-shrinkage threshold in ME SDs (`Inf` = off; DP argument `prob_n_sigma_me`) |
| `birth_death` | `TRUE` | Birth-death cells; `FALSE` = count-based slots (DP argument `prob_birth_death`) |
| `pin_truestemid` | `TRUE` | Pins while sampling and in the sweep |
| `dbh_round_censuses`, `dbh_round_max`, `dbh_round_width` | none, 5.5, 0.5 | DBH classes |
| `use_bio_hard_shrink_in_prob` | `TRUE` | `Bio_Max_Shrink` gate on links and the cumulative-shrinkage check; `FALSE` allows such links (soft `k_shrink` only) |
| `use_bio_hard_growth_in_prob` | `TRUE` | `Bio_Max_Growth` gate on links; `FALSE` allows them (soft `k_growth` only) |

#### DP and matcher compared

| Aspect | DP (global) | Matcher (per census pair) |
| --- | --- | --- |
| Scope | Minimises the total cost over the whole tag | Perturbed assignment per pair, stitched from the anchor |
| Runs of shrinkage | Each step pays its growth likelihood, so runs accumulate cost | Explicit cumulative-shrinkage cut (when on) |
| Hard bounds | Pruned before costing; inside the prune bounds but outside `Bio_*`: `hard_penalty`. Not affected by `USE_BIO_HARD_*_IN_PROB` | Links outside the prune bounds (or `min_growth` / `max_growth`) forbidden; `Bio_*` bounds with the `use_bio_hard_*` flags |
| Soft penalties | `k_shrink` on every decrease; `k_growth` above `Bio_Max_Growth_Soft` | `k_shrink` on every decrease; `k_growth` above `Bio_Max_Growth` (acts only when the growth gate is off) |
| Survival | Mortality cost on the track that dies | Log survival added to every link |
| Measurement error | Four-component mixture (optional) | No mixture; 0.0062 and 0.0904 only scale the cumulative check |
| Posterior | Exact marginals and samples | Heuristic spread of perturbed samples |

With `use_bio_hard_shrink_in_prob = FALSE` the matcher is more permissive than
the DP for the same tree.

### Post-engine helpers

The drivers apply, to each chunk and in this order (all in
`R/dp_global_main.R`):

1. **`maybe_add_posterior_bins()`**: `DP_PosteriorBin`.
2. **`apply_pin_track_rejoin()`**: the pin sweep moves pinned rows only, so an
   unpinned stretch the engine had linked to them can be left on another ID
   (e.g. a stem renumbered by the database before 2010 with only the later
   StemID pinned, or a pinned track offset by the resprout split). The whole
   left-behind ID moves onto the pin when the engine track segment holds
   exactly one pin value, all measured rows of the ID belong to the segment,
   none of its rows is pinned and the pin has no row in those censuses
   (`pin_track`). Segments end at every broken-below row with a DBH (with the
   main-stem rule only where it would continue the main stem) and after every
   NA-R barrier census, so those splits stay. Pinned rows never move;
   `TrueStemID` is unchanged.
3. **`apply_carried_terminal_backfill()`**: an unmeasured row with status
   dead, stem dead, broken below or missing that the engine left without a
   stem takes the latest earlier `ReconstructedStemID` of its
   `(Tag, StemID)` (LOCF; `OriginalStemID` without `StemID`), closing that
   stem's trajectory (`carried_terminal`). Rows with no earlier record of
   their StemID stay `NA`.
4. **`apply_orphan_stem_backfill()`**: a row with a StemID but no DBH, no pin
   and no engine stem (typically a new StemID first recorded broken below
   without a DBH) takes its StemID (`given_orphan`).
5. **`apply_terminal_to_host()`**: the BCI database usually stores a stem's
   break record under the new StemID of the resprout measured later, so pins
   and the sweep put it on the resprout, where R2 would leave it as a stem
   that never lives. An unmeasured dead, stem-dead or broken-below record that
   comes before the first evidence of life (DBH or status alive) of its stem,
   or sits on a stem that never lives, moves to the stem that ended
   (`terminal_to_host`), per group of records (Tag, stem, StemID):
   - (a) the engine's own choice before the sweep, when that stem was last
     alive before the record; if a later step relabelled that track (e.g.
     `apply_pin_track_rejoin()`), the current ID of its last measured row
     before the record. `SweepRollbackToPreSweep` is set;
   - (b) else the one stem of the tag whose last evidence of life is in the
     census before, when exactly one record group starts in that census;
   - (c) else, when the group's first record is the tree's main record
     (`PrimaryStem == "main"`), the main stem of the census before: the one
     measured main record of that census, if its last evidence of life is
     there. Up to 2010 a break is written on the main stem's record, so when
     the whole tree broke (several stems end together and (b) finds no single
     one) the stump is the main stem's end record;
   - otherwise the records stay (a stem that is never alive).

   A missing record of the same StemID that follows the group's first record
   moves with it (e.g. a 1985 stump recorded missing in 1990). A record never
   moves onto a stem with a row in one of its censuses; measured rows never
   move; each record keeps its own `Status`; `TrueStemID` is unchanged (the
   sweep is not rerun). With hard pins a pinned record moves only under the
   exception of `pin_prelife_rows()`. Idempotent.
6. **`apply_broken_below_invariants()`**: two life-cycle rules, deterministic
   and idempotent, per reconstructed stem (rows of one Tag and
   `ReconstructedStemID`) in census order. The stem, not the `StemTag`, is the
   unit, because StemTags exist only from 2010: a StemTag group would keep a
   stem's pre-2010 rows apart from its later rows.
   - **R1 (split on break)**: a broken-below row with a DBH that has earlier
     rows on its stem starts a new stem: a broken stem measured again is the
     resprout, a new physical stem. The new ID is carried to the later rows
     until the next such row (`bb_split`, `bb_split_carry`). With the
     main-stem rule R1 is a check: only such a row whose previous measured row
     is the main stem of the tag's previous measured census starts a new stem;
     the engines do not make these links.
   - **R2 (terminate on stump)**: a stump ends its stem; the first later
     measured row starts a new stem, carried forward the same way
     (`bb_post_terminator_split`, `bb_post_terminator_split_carry`).
     Unmeasured rows between (later stumps, dead or stem-dead records) keep
     the old ID.

   Dead and stem-dead rows trigger neither rule. New IDs are taken above the
   largest ID of the table. Without hard pins a pin does not protect a row
   (the BCI driver pins measured broken-below rows to their own StemID, which
   the database sometimes reuses from the trunk): a relabelled row whose
   `TrueStemID` equals the replaced ID gets the new ID as its `TrueStemID`,
   and the count is logged. With hard pins, see "Pins". The sweep is not rerun
   after this pass. `apply_bb_invariants_to_samples()` applies the same rules
   to each posterior sample.
7. **`renumber_engine_minted_ids()`**: numbers the stems of each tag 1..N in
   order of first census; in a tie the largest DBH in that census gets the
   lower number, then the original ID. No ID is zero or negative. The mapping
   is applied to `ReconstructedStemID`, `ReconstructedStemID_PreSweep` and the
   Top-k IDs (in the last two, engine labels that are no longer a final ID get
   numbers after N), and to `TrueStemID` only where the engine minted it
   (`provisional_dp`, and without hard pins the R1/R2 rows). It returns the
   mapping for the next step.
8. **`finalize_posterior_paths()`**: writes the posterior path files in the
   renumbered ID space (section "Posterior path files").

The generic drivers (`main_cpp.R`, `main_cpp_chunk.R`, `main_cpp_bci.R`) and
the BCI driver use this chain, whatever engine solved the tag.

### Measurement rejoin

`apply_measurement_rejoin()` (`R/measurement_rejoin.R`) runs later, in the
BCI merge, after the renumbering. When the point of measurement (POM) of a
trunk moves (buttresses growing up the trunk, deformities, strangler figs;
the new height recorded or not), the next DBH can fall outside the hard growth
bounds. The taper correction cannot remove buttress flare, and pruning removes
such a link before costs are compared while an oversized recruit is only
penalised, so the engines book a stem loss plus an impossible recruit for one
trunk. The rejoin joins a stem that is an impossible recruit (it starts at or
above the recruit limit in a tree measured before) to a stem of the same tree
whose measurements end in the census before (or two censuses before, when that
stem has no measurement in between) only when all of these hold:

- **Not a break or a death**: the later measurement has no break or resprout
  code and the earlier stem has no dead, stem-dead or broken-below record in
  between; the end and start are clean (no census gets two rows).
- **The size change is that of a moved POM**: between consecutive censuses
  the taper-corrected growth is outside the hard bounds (a split the engine
  chose inside the bounds stays); a smaller DBH down to a taper-corrected ratio
  of 0.4 (earlier stem ≥ 10 cm); a larger DBH only when the recorded HOM
  changed and the recorded DBH grew inside the bounds, ratio at most 1.5;
  across a census without a measurement, only a recorded HOM rise with a
  smaller DBH, with no other measurement of the tree in between that could be
  the same trunk.
- **The tree shows that its POM moves** (any census): a measurement away from
  1.3 m, a buttress code (B), a growth form without a stable POM
  (`irregular_forms`; BCI: strangler figs), or a species that usually has
  buttresses (`mr_buttress_species()`: ≥ 80% of its trees ≥ 30 cm carry B)
  when neither stem is measured in a census that recorded buttress codes (BCI:
  from 1990). A census without HOM records is not a sign on its own.
- **One identity**: the two stems carry at most one pin, and the pair is
  clearly the best (same HOM first, then the size ratio closest to 1, with
  less than half the deviation of the next candidate); consecutive pairs are
  resolved before pairs across a census, and when no pairing is clearly best
  nothing is joined. Database StemIDs are not used.

The joined stem takes the ID of the stem with the pin, else of the later stem
(IDs are not renumbered afterwards, so a changed tree skips one number per
join; `measurement_rejoin`). The joins go to `measurement_rejoin_audit.csv`
(with the signs found in `evidence`), and their observation pairs
(`measurement_rejoin_pairs.csv`) are joined by
`apply_measurement_rejoin_to_paths()` in each posterior sample where they are
split cleanly with at most one pin. The file header numbers and justifies
every rule.

## Parameter estimation

`estimate_bio_pars()` (`R/dp_global_bio.R`) estimates the parameters of one
parameter set from stems of known identity; the drivers call it per species
(and for pooled sets) and attach the result as `Bio_*` columns.

- **Input**: `Tag`, `CensusID`, `DBH` (cm), `TrueStemID`, `ExactDate` and
  optionally `species`. Only rows with `DBH > 0` and a `TrueStemID` are used.
- **Censuses** (`anchor_start_census`): the estimation grid holds every stem
  at every census from `anchor_start_census` on, so only those census pairs
  enter the estimates (default: the global `ANCHOR_START_CENSUS`).
- **Intervals**: computed per stem from `ExactDate` (date at t₁ minus date at
  t₀, in years of 365.25 days); a stem without a record at a census takes its
  tag's mean date there, else the census mean. `res$interval$per_pair_intervals`
  holds the intervals of the growth pairs.
- **Sources**: each guardrail and soft penalty has `*_source` `"data"`
  (estimated) or `"fixed"` (`*_fixed`): `recruit_max_source`,
  `max_shrink_source`, `max_growth_source`, `k_shrink_source`,
  `k_growth_source`.
- **Errors**: the function stops when no row is usable, when fewer than two
  censuses or five growth pairs remain, and when a fixed source or an
  `enforce_*` option gets an invalid value.

### Growth

Wide tables of DBH and dates per stem (`dcast` of Tag + TrueStemID + species
by CensusID); for each adjacent census pair

$$g_i = \frac{D_{1,i} - D_{0,i}}{T_i}$$

with $T_i$ the stem's own interval.

**Mean**: `lm(g ~ log(d0))` when at least 10 pairs exist and log D₀ varies,
giving α and γ; otherwise α = mean(g), γ = 0. The predicted mean is
$\mu_{\text{pred}}(i) = \alpha + \gamma \log D_{0,i}$ (α when the fit failed or
$D_{0,i} \le 0$).

**Process SD.** A robust SD proxy from the absolute residuals: for a normal
$\mathbb{E}|X| = \sigma\sqrt{2/\pi}$, so

$$\hat\sigma_{\text{total},i} = |g_i - \mu_{\text{pred}}(i)| \sqrt{\pi/2}$$

With measurement error the error variance of the annual increment is removed
in quadrature:

$$\text{Var}_{\text{meas}}(g_i) = \frac{\text{Var}(\varepsilon \mid D_{0,i}) + \text{Var}(\varepsilon \mid D_{1,i})}{\Delta t^2}, \qquad \text{Var}(\varepsilon \mid D) = (1-p)\,\text{SD1}(D)^2 + p\,\text{SD2}^2, \qquad \text{SD1}(D) = \max(aD + b,\ 10^{-6})$$

$$\hat\sigma_{\text{proc},i} = \sqrt{\max\left(\hat\sigma_{\text{total},i}^2 - \text{Var}_{\text{meas}}(g_i),\ 10^{-8}\right)}$$

and without it $\hat\sigma_{\text{proc},i} = \max(\hat\sigma_{\text{total},i}, 10^{-6})$.

**Heteroskedasticity**: `lm(sd_proc_hat ~ d0)`:

- slope ≥ 0: σ₀ = intercept (at least 0.01), σ₁ = slope;
- slope < 0 (the SD falls with size, e.g. palms): a line would reach zero for
  large stems, and keeping its intercept with slope 0 would give every stem the
  SD extrapolated to DBH 0. The SD is refitted as a constant: σ₀ = mean of the
  per-pair proxies (at least 0.01), σ₁ = 0.

The model is $\sigma(D_0) = \sigma_0 + \sigma_1 D_0$.

**Enforced growth bounds** (`enforce_growth_bounds = TRUE`, with
`growth_min_fixed` / `growth_max_fixed`, cm/yr; one side may be `NA`): the
increments outside the bounds are dropped before the mean and SD are fitted
(with a warning giving the count). They affect only the data of the fit, not
the returned guardrails. The function stops when fewer than five pairs remain.
The drivers enforce the fixed growth bounds.

### Soft penalties

**k_shrink**, data mode:

- with measurement error and at least 10 growth pairs: with
  $s_{\text{typ}} = \text{median}\sqrt{\text{Var}(\varepsilon \mid D_0) + \text{Var}(\varepsilon \mid D_1)}$,
  the typical measurement SD of a DBH difference,
  $k_{\text{shrink}} = 1/(2 s_{\text{typ}}^2)$;
- otherwise from the observed decreases $\Delta D = D_0 - D_1$ (cm, pairs with
  $D_1 < D_0$): $k_{\text{shrink}} = 1/(2\,\text{Var}(\Delta D))$ with at least
  five decreases, else 50.

Clamped to $[10^{-6}, 10^6]$. Fixed mode: `k_shrink_fixed` (0 disables it).
A shrinkage of $s$ cm costs about one unit when $k \approx 1/(2s^2)$; k = 50
corresponds to s ≈ 0.1 cm. In the cost: $k_{\text{shrink}}(D_0 - D_1)^2$ when
$D_1 < D_0$.

**k_growth**: the same calibration with the positive increments
$\Delta D = D_1 - D_0$ ($D_1 > D_0$); in the cost
$k_{\text{growth}}(D_1 - D_{1,\text{cap}})^2$ above the soft cap. Fixed mode:
`k_growth_fixed`. The drivers fix both weights at 0.

### Hard bounds

**max_shrink**, data mode: the empirical lower quantile
$Q_{0.001}(g)$ (`shrink_data_quantile`). With measurement error also the
measurement-only quantile: the $q$ = `shrink_hard_prob` (10⁻⁴) quantile of the
mixture of $(e_1 - e_0)/\Delta t$,

$$F(x) = \sum_{i=1}^{4} w_i\, \Phi\left(x;\ 0,\ \text{SD}_{\text{meas},i}\right)$$

with the four components of case D, evaluated at the typical diameter
$D_{\text{typ}}$ = median($D_0$) and, as Δt, the median interval of the last
census pair. Then `max_shrink` = min(data, measurement) (the more permissive,
more negative value), or the data value alone. Fixed mode: `max_shrink_fixed`
(cm/yr).

**max_growth**, data mode: $Q_{0.999}(g)$ (`growth_data_quantile`); with
measurement error also the $1 -$ `growth_hard_prob` (99.99%) quantile of the
same mixture; `max_growth` = max(data, measurement). Fixed mode:
`max_growth_fixed` (cm/yr, positive). The soft cap is
`max_growth_soft` = min(`max_growth`, $Q_{0.99}(g)$) (`growth_soft_quantile`).

In the cost both are hard: $g$ outside [`max_shrink`, `max_growth`] costs
`hard_penalty`. The drivers fix both and keep the data quantile for the soft
cap.

### Mortality

For each adjacent census pair, the stems measured at t₀ are at risk; a stem
not measured at t₁ died (also when it is measured again later); Δt is the
stem's own interval. The model is
$P_{\text{death}}(D, \Delta t) = 1 - \exp(-h_0 e^{\beta D} \Delta t)$, fitted by
maximum likelihood (`optim`, BFGS, on $[\log h_0, \beta]$, start
`mortality_start` = (log 0.01, 0)) of

$$-\sum_i \left[ d_i \log P_{\text{death}}(D_{0,i}) + (1 - d_i) \log\left(1 - P_{\text{death}}(D_{0,i})\right) \right]$$

giving $h_0 = \exp(\text{par}_1)$, $\beta = \text{par}_2$.

### Recruitment

- **Events**: a grid cell without a DBH at t₀ and with a positive DBH at t₁.
- **Size**: LogNormal(meanlog, sdlog) fitted to the recruit DBH
  (`MASS::fitdistr`); with fewer than two recruits meanlog = log 2,
  sdlog = 0.5. `enforce_recruit_max = TRUE` drops recruits above
  `recruit_max_fixed` before the fit (a positive value is required).
- **Limit** `recruit_max_dbh`: data mode $Q_{0.999}$ of the recruit DBH
  (`recruit_max_quantile`; 5 cm without recruits), fixed mode
  `recruit_max_fixed`.
- **Rate** (`recruit_rate_unit`, driver `RECRUIT_RATE_UNIT`):
  - `"tree"` (default): new stems per established tree per year,
    $\lambda = (n_{\text{new}} + 0.5) / \sum \Delta t$. A tree is established
    in an interval when it has a measured stem at its start; $n_{\text{new}}$
    counts its stems without a DBH at the start and with one at the end; Δt is
    each established tree's interval (mean over its stems). The 0.5
    pseudo-recruit keeps the rate positive. This is the yearly probability
    that a tree gains a stem, which is what both engines use λ for (the DP's
    empty track, the matcher's recruitment cell);
  - `"slot"`: recruits per empty grid cell per year,
    $\lambda = n_{\text{recruits}} / \sum \Delta t$ over the cells without a
    DBH at the start. Those cells are mostly dead stems and stems that recruit
    later, so this rate barely varies between species (about 0.08/yr for
    nearly every BCI parameter set) although observed recruitment varies far
    more. It is also used when no tree is established.

  Both values are returned (`recruitment$lambda_tree`, `lambda_slot`,
  `lambda_unit`, `n_new_stems_established_trees`, `established_tree_years`).

### Quantiles

| Argument | Default | Tail | Computes | Used for |
| --- | --- | --- | --- | --- |
| `shrink_data_quantile` | 0.001 | lower | Empirical lower bound of g | `max_shrink` (data) |
| `shrink_hard_prob` | 10⁻⁴ | lower | Measurement-only lower bound | `max_shrink` (measurement) |
| `growth_data_quantile` | 0.999 | upper | Empirical upper bound of g | `max_growth` (data) |
| `growth_hard_prob` | 10⁻⁴ | upper | Measurement-only upper bound (99.99%) | `max_growth` (measurement) |
| `growth_soft_quantile` | 0.99 | upper | Empirical soft threshold | `max_growth_soft` |
| `recruit_max_quantile` | 0.999 | upper | Largest recruit DBH | `recruit_max_dbh` |

Empirical quantiles come from the observed g; measurement-only quantiles from
the four-component mixture at the typical diameter. Shrinkage combines them
with `min()` and growth with `max()`.

## Settings

Options read by the engine (`options(dp_global.<name> = …)`; the BCI driver
sets all of them):

| Option | Default | What it does |
| --- | --- | --- |
| `dp_global.bb_dbh_rule` | `"r1"` | `"main_ends"`: main-stem rule for broken-below records with a DBH |
| `dp_global.hard_pins` | `FALSE` | Pins hold after the engine |
| `dp_global.field_tag_first_census` | `Inf` | First census whose `StemTag`s are field tags |
| `dp_global.posterior_label_free` | `FALSE` | Label-free DP posterior |
| `dp_global.prob_export_rule` | `"mea"` | Export of the matcher: `"mea"` or `"fewest_me"` |

Arguments of `match_stems_dp_global_backward_marginals_batch()`:

| Argument | Default | What it does |
| --- | --- | --- |
| `anchor_start` | | Anchor census |
| `max_states` | 50,000 | State-space limit |
| `max_tracks`, `slack_tracks` | 30, 1 | Cap on K; slack tracks |
| `slack_require_anchor_recruitable`, `slack_require_anchor_eps` | `FALSE`, 10⁻⁶ | Slack tracks only with a recruitable anchor stem |
| `min_growth`, `max_growth` | −Inf, Inf | Base growth bounds: pruning when no prune bound is given, matcher, `ConstraintViolation` |
| `prune_hard`, `prune_min_growth`, `prune_max_growth`, `prune_use_bio_bounds` | `TRUE`, `NULL`, `NULL`, `TRUE` | Growth pruning |
| `prune_recruit_max_dbh`, `prune_use_bio_recruit` | `NULL`, `TRUE` | Recruit-size pruning |
| `non_taper_corrected_growth_forms`, `non_taper_corrected_prune_min_growth`, `non_taper_corrected_prune_max_growth`, `hom_tolerance_scale` | palm, strangler_fig, tree_fern; −0.625; 6.25; 2 | Non-taper bounds and HOM widening |
| `pin_truestemid`, `allow_provisional_anchor` | `TRUE`, `TRUE` | Pins and provisional anchors |
| `allow_segment_split` | `TRUE` | Resprout segment split |
| `fallback_growth_forms`, `prob_species` | none | Tags sent to the matcher |
| `use_measurement_error`, `meas_sd1_a`, `meas_sd1_b`, `meas_sd2`, `meas_p_big` | `FALSE`, 0.0062, 0.0904, 4.64, 0.05 | Measurement-error mixture |
| `dbh_round_censuses`, `dbh_round_max`, `dbh_round_width` | none, 5.5, 0.5 | DBH classes |
| `temperature`, `eps_tiebreak` | 1, 10⁻⁶ | Posterior temperature; tie-break weight |
| `posterior_top_k`, `posterior_samples`, `posterior_samples_format`, `posterior_samples_path`, `posterior_sample_seed` | 2, 0, `"rds"`, `NULL`, `NULL` | Posterior summaries and samples |
| `prob_n_samples`, `prob_lookahead_weight`, `prob_n_sigma_me`, `prob_birth_death` | 200, 0.5, 3, `TRUE` | Matcher |
| `use_bio_hard_shrink_in_prob`, `use_bio_hard_growth_in_prob` | `TRUE`, `TRUE` | Matcher hard gates |
| `verbose`, `chunk_id` | `FALSE`, `NULL` | Log output |

Roles of the bounds and penalties (typical values of `main_cpp.R`:
`MAX_GROWTH_FIXED` = 7.5, `PRUNE_BOUND_FACTOR` = 1.25,
`RECRUIT_MAX_FIXED` = 7.5 × 5 + 0.9999):

| Type | Parameter | Used by | Purpose | Typical value |
| --- | --- | --- | --- | --- |
| Hard constraint | `Bio_Max_Shrink` | DP cost | Lower bound of g | −0.5 cm/yr |
| Hard constraint | `Bio_Max_Growth` | DP cost | Upper bound of g | 5 or 7.5 cm/yr |
| Soft penalty | `Bio_K_Shrink` | DP cost | Quadratic shrinkage weight | 0 (drivers) to 50 1/cm² |
| Soft penalty | `Bio_K_Growth` | DP cost | Quadratic growth-excess weight | 0 (drivers) to 50 1/cm² |
| Soft threshold | `Bio_Max_Growth_Soft` | DP cost | Growth-excess threshold | min(`Bio_Max_Growth`, 99% quantile of g) |
| Base bound | `min_growth`, `max_growth` | Pruning (no prune bound), matcher, `ConstraintViolation` | Base growth bounds | −0.5, 7.5 cm/yr |
| Prune bound | `prune_min_growth`, `prune_max_growth` | Pruning, matcher | `PRUNE_BOUND_FACTOR` × base | −0.625, 9.375 cm/yr |
| Prune bound | `prune_recruit_max_dbh` | Pruning, matcher | `PRUNE_BOUND_FACTOR` × `RECRUIT_MAX_FIXED` | 48.12 cm |
| Matcher | `prob_n_samples` | Matcher | Gumbel samples | 200 |
| Matcher | `prob_lookahead_weight` | Matcher | Look-ahead weight | 0.5 (function), 1 (drivers) |

The driver settings are in `scripts/README.md` and in the BCI stage-2 README.

## Checks and diagnostics

- The DP stops when `temperature` is not positive or a required `Bio_*` column
  is missing, `NA` or not constant within the tag; `estimate_bio_pars()` stops
  on unusable data or invalid fixed values.
- `attr(out, "DP_PruneInfo")` (section "Pruning"), `attr(out, "DP_Compute_Profile")`
  (number and time of the cost calls), `DP_Sampling_Profile` when samples are
  drawn.
- `ConstraintViolation` and `DP_FallbackReason` on every row.
- Posterior screening:

  ```r
  high_uncertainty <- out[DP_PosteriorEntropy > 1.5]
  ambiguous <- out[DP_PosteriorBin == "ambiguous"]
  ```

- Cost of one transition, by term and track:

  ```r
  components <- transition_cost_tracks_bio_components(track_dbh_t, track_dbh_tp1, interval_years, ...)
  # components$per_track: one row per track with its case and terms (cost_recruit, cost_no_recruit,
  #   cost_mortality, cost_growth_lik, cost_shrink_soft, cost_growth_soft, cost_hard, total_track);
  # components$tiebreak, components$total, components$p_recruit
  # (no DBH-rounding arguments: equals the compiled cost at censuses not flagged as rounded)
  ```

- `plot_tag_to_pdf(out, pdf_file, include_reference = TRUE, tag = NULL)`:
  reconstructed trajectories per tag (with the `TrueStemID` panel).
- `R/sensitivity_transition_cost_bio.R` (`build_all_sweeps()`,
  `plot_all_sweeps_to_pdf()`): the cost of fixed scenarios along each
  parameter; an opt-in demo runs on the simulated data (file header).
  `R/realism_calibration.R` (`realism_report_from_reconstruction()`): summary,
  per-group table and tuning suggestions. `R/k_tuning_viz.R`: join-versus-split
  cost as k varies.

### Complexity estimation

Before a run, `R/complexity/estimate_dp_complexity_function.R` (needs `here`)
estimates the work per tag:

```r
library(here)
source(here("dp_global", "R", "complexity", "estimate_dp_complexity_function.R"))
complexity <- estimate_dp_complexity(here("data_simulation", "data", "simulated_data_1.csv"))
details <- get_tag_complexity_details(here("data_simulation", "data", "simulated_data_1.csv"), tag = 11)
details$summary      # the tag's row of estimate_dp_complexity()
details$per_census   # CensusID, N_Obs, N_States = P(K, N_Obs)
```

Columns: `Tag`, `Species`, `census_range`, `n_censuses` (censuses with a DBH
up to the anchor), `K`, `max_obs`, `max_states_per_census`, `total_states`,
`estimated_edges_unpruned` (candidate transitions between adjacent censuses, a
proxy for run time), `estimated_edges_pruned` (after the hard pruning bounds),
`estimated_fallback`, `fallback_reason` (a census above `max_states`),
`predicted_seconds`, `predicted_hours` (from the edge count, fitted log-log
polynomial). The estimator's K is
`max(#anchor IDs, max_obs, n_obs(1) + births)` plus one track per resprout
observation plus slack tracks, so it can differ from the DP's K (which adds no
track for resprout observations and adds tracks for pins before the anchor).
States per census are P(K, n_obs); transitions count all pairs of states of
consecutive censuses. Tags without fallback come first, by decreasing number of
transitions, so the slowest tags appear first.

### Numerical stability

| Operation | Mechanism | Value |
| --- | --- | --- |
| Growth SD | Floor | σ(D) ≥ 10⁻⁶ |
| Measurement SD1 | Floor | SD1(D) ≥ 10⁻⁶ |
| Mortality and recruitment probabilities | Range clamp | [10⁻¹², 1 − 10⁻¹²] |
| Log-sum-exp (mixture, recursion) | Maximum subtracted before `exp` | No overflow or underflow |
| k in data mode | Range clamp | [10⁻⁶, 10⁶] |
| Rank tie-break | Stable ranks (equal DBHs keep track order) | `std::stable_sort` |

### Performance

- States grow factorially: P(K, n) per census, at most Kⁿ.
- Set `DP_MAX_STATES` for the available memory; it also sets the edge limit
  `max_states²`.
- `PROB_N_SAMPLES` trades matcher accuracy for speed (200).
- Tags run in parallel with `MANUAL_CORES = TRUE` and `MANUAL_CORES_VALUE`.
- The DP costs all feasible transitions of a census pair in one C++ call.

## Implementation reference

### Function map

| Task | Function | Location |
| --- | --- | --- |
| State enumeration | `enumerate_states_injective()`, `enumerate_states_constrained()` | `R/dp_global_states.R` |
| State counts | `count_injective_states()` | `R/dp_global_utils.R` |
| Track DBH vector of one state (helper) | `state_to_track_dbh()` | `R/dp_global_states.R` |
| Transition costs (C++) | `transition_cost_paired_rcpp()` (DP), `transition_cost_tracks_bio_batch_rcpp()` | `src/transition_cost_rcpp.cpp`, wrappers in `src/transition_cost_rcpp.R` |
| Phase feasibility and pruning (C++) | `derive_phase_prev_batch_rcpp()` | `src/transition_cost_rcpp.cpp` |
| Cost breakdown | `transition_cost_tracks_bio_components()` | `R/dp_global_bio.R` |
| DP solver | `match_stems_dp_global_backward_marginals_batch()` | `R/dp_global_dp.R` |
| Fallback routing | `do_fallback()` (inside the DP function) | `R/dp_global_dp.R` |
| R-boundary splitting | `split_live_resprout_tracks()` | `R/dp_global_dp.R` |
| Samples of a split tag | `pair_segment_posterior_samples()`, `segment_fixed_assignment()`, `stage_posterior_samples()` | `R/dp_global_dp.R` |
| Probabilistic matcher | `match_stems_probabilistic()` | `R/dp_probabilistic_matching.R` |
| Link scores; assignment matrix | `compute_pairwise_log_likelihood()`, `augment_cost_matrix()` | `R/dp_probabilistic_matching.R` |
| Perturbed assignment | `greedy_assignment_gumbel()`, `enforce_feasible_assignment()`, `hungarian_min_rcpp()` | `R/dp_probabilistic_matching.R`, `src/transition_cost_rcpp.cpp` |
| Pins while sampling | `apply_pin_mask()`, `apply_track_pin_mask()`, `propagate_track_pin()`, `pin_masked_pair()`, `apply_main_end_mask()` | `R/dp_probabilistic_matching.R` |
| Look-ahead, repair, filter | `condition_cost_matrix()`, `repair_stitched_growth_violations()`, `filter_pin_consistent_samples()` | `R/dp_probabilistic_matching.R` |
| Marginals and export | `compute_marginals_from_samples()`, `select_consensus_trajectory()`, `count_beyond_me_links()`, `export_probabilistic_posteriors()` | `R/dp_probabilistic_matching.R` |
| DBH classes | `transition_cost_paired_rcpp(round_t, round_tp1, …)`, `compute_pairwise_log_likelihood(round_curr, round_next, …)` | `src/transition_cost_rcpp.cpp`, `R/dp_probabilistic_matching.R` |
| Posterior bins; growth diagnostic | `add_dp_posterior_bins()`, `add_constraint_violation()` | `R/dp_global_diag.R` |
| Hard pins and field tags | `hard_pins_active()`, `field_tag_rows()`, `pin_prelife_rows()`, `bb_main_violation()` | `R/dp_global_main.R` |
| Post-engine helpers | `apply_pin_track_rejoin()`, `apply_carried_terminal_backfill()`, `apply_orphan_stem_backfill()`, `apply_terminal_to_host()`, `apply_broken_below_invariants()`, `renumber_engine_minted_ids()`, `finalize_posterior_paths()` | `R/dp_global_main.R` |
| Per-sample pins and rules | `apply_pins_to_samples()`, `apply_bb_invariants_to_samples()` | `R/dp_global_main.R` |
| Measurement rejoin | `apply_measurement_rejoin()`, `apply_measurement_rejoin_to_paths()`, `mr_buttress_species()` | `R/measurement_rejoin.R` |
| Parameter estimation | `estimate_bio_pars()` | `R/dp_global_bio.R` |
| Complexity estimate | `estimate_dp_complexity()`, `get_tag_complexity_details()` | `R/complexity/estimate_dp_complexity_function.R` |
| Plots | `plot_tag_to_pdf()` | `R/dp_global_diag.R` |
| Drivers | `run_main()` and `run_dp_one_group()` (`main_cpp.R`), `run_main_chunked()` (`main_cpp_chunk.R`), `main_cpp_bci.R` (sources `main_cpp.R`) | `scripts/` |

### Inside the DP function

- **Anchor selection**: section "Anchor census".
- **State enumeration**: `enumerate_states_injective(K, n_obs, max_states)`
  returns `NULL` above `max_states` (fallback);
  `enumerate_states_constrained(K, n_obs, allowed_tracks, max_states)` when
  pins or growth bounds restrict tracks. `count_injective_states()` gives the
  theoretical counts used for diagnostics and the segment split check.
- **Track DBH matrices**: per census an `n_states × K` matrix of the DBH on
  each track in each state (`track_dbh_by_state`), used for the feasibility
  test and the costs.
- **Phase constraints**: `derive_phase_prev_batch_rcpp()` tests phase
  feasibility, the R-recruit constraint and the pruning bounds for all state
  pairs in one batch and returns the feasible pairs (`from_i`, `to_j`) with
  their phase matrix.
- **Intervals**: from the per-census mean `ExactDate`; a pair with an `NA`,
  non-finite or non-positive interval is not pruned.
- **Costs**: `transition_cost_paired_rcpp(tdbh0_mat, tdbh1_mat, …)`, one call
  per census pair.
- **Recursion**: Viterbi values and back-pointers (`vit_cost`, `vit_ptr`) and
  log-sum-exp weights.
- **Forward pass**: normalised state weights give the per-observation
  posteriors.
- **Sampling**: backward sampling with per-sample `logp`; staged in
  `.staging/` with attributes `DP_Posterior_Staging_File` and
  `DP_Sampling_Profile`.
- **NA in the DP**: an empty track means "no stem", not "alive but missed";
  missing-from-field censuses are taken out before the DP.
- **Hard constraints**: `hard_penalty` = 10⁶ encodes impossible transitions;
  the DP takes such a path only when no alternative exists.

### Common questions

1. **Why a track-based global DP?** It optimises over the whole history, so
   local choices cannot cause later identity swaps, and it enforces the
   life-cycle constraints exactly.
2. **Why both hard bounds and soft penalties?** Hard bounds are binary
   (allowed or forbidden) and prevent absurd links; soft penalties
   discourage implausible links continuously and help choose among many
   feasible paths.
3. **Data or fixed mode?** Data mode when there are reliable known
   identities; fixed mode for literature values or sparse data. They can be
   mixed (e.g. data hard bounds, fixed soft penalties).
4. **What does measurement error change?** The growth likelihood becomes a
   four-component mixture (more robust to noisy measurements), the estimation
   separates process from measurement variance, and the bounds become more
   conservative.
5. **`min_growth` versus `Bio_Max_Shrink`?** `Bio_Max_Shrink` is in the DP
   cost; `min_growth` is a pruning bound (when no prune bound is given), the
   matcher's bound and the `ConstraintViolation` bound. They can differ.
6. **Is DBH → NA a missing measurement?** No: inside the engines it is a
   death, except in censuses detected as missing from the field. K matters:
   too small fails, too large is slow; temperature controls confidence, not
   biology.

## Basal area uncertainty (`scripts/basal_area_uncertainty.R`)

The script propagates identity uncertainty into basal area (BA). The BA of a
tag in a census, Σ π/4 (DBH/100)² over its measured stems (m²), does not
depend on which stem a measurement belongs to. What depends on identity is the
per-stem BA trajectories and the split of each BA change into growth of
surviving stems, loss by mortality and gain by recruitment.

```sh
Rscript dp_global/scripts/basal_area_uncertainty.R --RUN_DIR=dp_global/output/<run_dir>
```

It reads `stem_reconstruction_dp_global_rcpp.csv` and the CSV path files of
the run (`POSTERIOR_SAMPLES_FORMAT = "csv"`, the drivers' default; without
them only the MAP columns are produced). Each path is weighted by its share of
the tag's samples, `path_count / sum(path_count)`. DBH is in cm. The path
files hold the observations up to the anchor only. An observation outside the
census window of the paths (after the anchor, or before the paths' first
census) joins, in every path, the path stem that holds its exported stem-mate
at the nearest window edge (the last in-path observation of the same exported
stem for later censuses, the first for earlier ones); an exported stem without
in-path observations keeps its own stem. Every census pair is thus decomposed
from one consistent set of stems, and the weighted total change
(`DeltaBA_check`) equals `DeltaBA_total`. When a measured observation inside
the window is missing from the paths, the tag's paths are replaced by the
exported reconstruction. The script reports both counts.

| File | Content |
| --- | --- |
| `basal_area_tag_census.csv` | Per tag and census: `Tag`, `CensusID`, `Year` (mean ExactDate), `TotalBA_m2`, `NumStems` |
| `basal_area_tag_change.csv` | Per tag and census pair: `Interval_yr`; MAP `Growth_BA`, `Loss_BA` (negative), `Gain_BA`, `DeltaBA_total` (identity-free), `NumSurvivors`, `NumDeaths`, `NumRecruits`; posterior mean, SD and 2.5% / 97.5% quantiles of growth, loss and gain; `DeltaBA_check` (weighted mean of the total change); mean stem counts; `NumPaths` |
| `basal_area_figures.pdf` | Summary pages, per-tag panels (BA trajectory, stem count, decomposition with uncertainty, demographics), uncertainty histograms and posterior densities |

## Tests

From the project root:

```r
testthat::test_file("dp_global/tests/test_estimate_bio_pars.R")          # growth SD refit, recruitment rate units
testthat::test_file("dp_global/tests/test_probabilistic_birth_death.R")  # birth-death matcher, assignment solver
testthat::test_file("dp_global/tests/test_enumerate_states.R")           # state enumeration above max_states
```

## Requirements

R packages `data.table`, `Rcpp`, `igraph`, `MASS` (checked by
`dp_global_main.R`), `lpSolve` (matcher), `parallel`, `here` (drivers);
optional `arrow` (Feather), `ggplot2`, `cowplot`, `reshape2`, `scales`,
`patchwork` (plots), `testthat`. A C++ compiler: `dp_global_main.R` compiles
`src/transition_cost_rcpp.cpp` with `Rcpp::sourceCpp()`.

Further documentation: `scripts/README.md` (drivers), `src/README.md` (C++
functions), `R/dpglobal_bundle/README.md` (bundle), and the function headers
in the R files.

## References

- Chave J, Condit R, Aguilar S, Hernandez A, Lao S, Perez R (2004). Error
  propagation and scaling for tropical forest biomass estimates.
  *Philosophical Transactions of the Royal Society B* 359: 409–420.
- CTFS R package, growth tutorial:
  <https://ctfs.si.edu/ctfsdev/CTFSRPackageNew/index.php/web/tutorials/GrowthChange/index.html>
- Piponiot C, Condit R, Hubbell SP, Pérez R, Lao S, Aguilar S, Muller-Landau H
  (2024). Woody biomass stocks and fluxes in the Barro Colorado Island 50-ha
  plot (Appendix S1).
