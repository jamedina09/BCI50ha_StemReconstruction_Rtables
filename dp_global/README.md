# dp_global

Engine that reconstructs the stem identities of multiple-stem trees across
forest censuses. In each census a tree (tag) has a set of stem measurements,
but which measurement belongs to which stem is known only from some census on
(the anchor census: field tags, database IDs). The engine links the earlier
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

Author: José A. Medina-Vega

## Contents

| File | What it does |
| --- | --- |
| `R/dp_global_main.R` | Loads the engine (compiles the C++ file, sources the modules) and defines the post-engine helpers |
| `R/dp_global_dp.R` | DP solver `match_stems_dp_global_backward_marginals_batch()`, resprout split, posterior staging |
| `R/dp_probabilistic_matching.R` | Probabilistic matcher `match_stems_probabilistic()` |
| `R/dp_global_bio.R` | Parameter estimation `estimate_bio_pars()` and the cost breakdown `transition_cost_tracks_bio_components()` |
| `R/dp_global_states.R` | State enumeration (`enumerate_states_injective()`, `enumerate_states_constrained()`) |
| `R/dp_global_utils.R` | State counts and census intervals |
| `R/dp_global_diag.R` | Posterior bins, `ConstraintViolation` and trajectory plots |
| `R/measurement_rejoin.R` | Rejoin of stems split by a moved point of measurement (run by the BCI merge) |
| `R/naming_helpers.R` | Run folder names |
| `R/dp_global_matchers.R` | Stepwise igraph matcher (not used by the DP workflow) |
| `R/sensitivity_transition_cost_bio.R`, `R/realism_calibration.R`, `R/k_tuning_viz.R`, `R/check_functions.r` | Sensitivity sweeps, realism report, k-tuning plots, parameter inspection |
| `R/complexity/` | DP state counts and run-time estimate per tag |
| `R/dpglobal_bundle/` | Packaging of the engine for other machines |
| `scripts/` | Drivers (`main_cpp.R`, `main_cpp_chunk.R`, `main_cpp_bci.R`) and `basal_area_uncertainty.R` |
| `src/` | C++ transition costs, feasibility checks and assignment solver, with R wrappers |
| `tests/` | testthat tests |

## How to run

From the project root (`dp_global_main.R` finds its files from the working
directory):

```sh
Rscript dp_global/scripts/main_cpp.R          # one tag or a small data set
Rscript dp_global/scripts/main_cpp_chunk.R    # large data sets, in chunks, resumable
```

```r
source("dp_global/scripts/main_cpp.R")   # settings and functions; starts no run
res <- run_main()
```

`scripts/README.md` lists the driver settings. The BCI run uses its own
driver, `BCI_stem_reconstruction/2_STEM_IDENTIFICATION/1_main_cpp_chunk_bci.R`.

The engine can also be called directly. A minimal workflow for one tag:

```r
source("dp_global/R/dp_global_main.R")

# 1. Parameters of one species, from the stems of known identity
bio <- estimate_bio_pars(stems_sp, anchor_start_census = 7L)

# 2. Bio_* columns on the rows of that species (the drivers' attach_bio_columns())
# 3. The DP, one tag at a time
out <- match_stems_dp_global_backward_marginals_batch(
  tree_data = tag_rows, anchor_start = 7L,
  min_growth = -0.5, max_growth = 5,
  max_states = 10000L, slack_tracks = 1L,
  posterior_samples = 200L, posterior_samples_path = out_dir
)

# 4. Post-engine helpers, in the drivers' order (section "Post-engine helpers")
out <- add_dp_posterior_bins(out)
out <- apply_pin_track_rejoin(out)
out <- apply_carried_terminal_backfill(out)
out <- apply_orphan_stem_backfill(out)
out <- apply_terminal_to_host(out)
out <- apply_broken_below_invariants(out)
ren <- renumber_engine_minted_ids(out, posterior_samples_path = out_dir)
finalize_posterior_paths(ren$out, posterior_samples_path = out_dir, mapping = ren$mapping)
```

## Inputs

One row per stem record (a stem in a census). Required columns:

| Column | Content |
| --- | --- |
| `Tag` | Tree: the unit the engine solves |
| `CensusID` | Census number 1, 2, 3, … |
| `DBH` | Diameter in cm (`NA` = not measured) |
| `ExactDate` | Measurement date; the interval of a census pair is the difference of the mean census dates of the tag |
| `species` | Parameter set of the row |
| `TrueStemID` | Pin: known stem identity, at the anchor census and at any earlier row to pin; `NA` where the engine is free |

Used when present:

| Column | Used for |
| --- | --- |
| `Status`, `ListOfTSM` | Field status and codes: breaks, resprouts, stumps, missing stems |
| `growth_form`, `hom` / `HOM` | Pruning bounds of non-taper-corrected forms and routing to the matcher |
| `Species` | Routing species to the matcher (`prob_species`) |
| `PrimaryStem` | Main-stem rule (`"main"`) |
| `StemTag` | Field tags (hard pins) |
| `StemID` or `OriginalStemID` | Database stem ID: post-engine backfills and grouping of terminal records |

Biological parameter columns, one value per tag (`estimate_bio_pars()` gives
them; section "Parameter estimation"):

| Column | Symbol | Content |
| --- | --- | --- |
| `Bio_Mu_Growth`, `Bio_Gamma_Growth` | α, γ | Mean annual growth μ(D) = α + γ log D (cm/yr) |
| `Bio_Sigma0_Growth`, `Bio_Sigma1_Growth` | σ₀, σ₁ | Growth SD σ(D) = σ₀ + σ₁ D (cm/yr) |
| `Bio_Max_Shrink`, `Bio_Max_Growth` | | Hard bounds of annual growth (cm/yr) |
| `Bio_Max_Growth_Soft` | | Soft cap of annual growth (cm/yr) |
| `Bio_K_Shrink`, `Bio_K_Growth` | k_shrink, k_growth | Weights of the soft shrinkage and growth penalties (1/cm²) |
| `Bio_H0_Mortality`, `Bio_Beta_Mortality` | h₀, β | Mortality hazard h₀ exp(β D) (1/yr) |
| `Bio_Recruit_Meanlog`, `Bio_Recruit_Sdlog` | | Lognormal size of a new stem |
| `Bio_Recruit_MaxDBH_unit` | | Largest DBH of a new stem (cm) |
| `Bio_Recruitment_lambda` | λ | Recruitment rate (1/yr) |

`Bio_Gamma_Growth`, `Bio_Max_Growth`, `Bio_Max_Growth_Soft` and `Bio_K_Growth`
are optional (defaults 0, no bound, no cap, 0). The DP stops when another
column is missing, `NA` or not constant within the tag.

## Outputs

Columns added to every row:

| Column | Content |
| --- | --- |
| `ReconstructedStemID` | Reconstructed stem, numbered 1..N within each tag |
| `ReconstructionMethod` | How the row got its stem (table below) |
| `DP_PosteriorTop1ID`, `DP_PosteriorTop1Prob`, `DP_PosteriorTop2ID`, `DP_PosteriorTop2Prob` | The most probable stems of the row and their probabilities (`posterior_top_k` pairs) |
| `DP_PosteriorEntropy` | Entropy of the row's posterior, −Σ p log p |
| `DP_PosteriorReconstructedProb` | Posterior probability that the row belongs to its exported stem |
| `DP_PosteriorUnlinkedProb` | Posterior probability that the row belongs to no anchor stem (DP only; `NA` from the matcher) |
| `DP_PosteriorBin` | `confident`, `ambiguous` or `unlinked-likely` (`add_dp_posterior_bins()`) |
| `DP_KUsed`, `DP_MaxStatesPerCensus`, `DP_MaxStatesCensusID` | Tracks used, the largest state count P(K, n) of a census (before track restrictions) and that census |
| `DP_FallbackReason` | Why the tag went to the probabilistic matcher; `NA` when the DP solved it |
| `ConstraintViolation` | Annual growth of the reconstructed stem outside `[min_growth, max_growth]` (diagnostic) |
| `obs_row_id` | Row number within the tag; the key of the posterior path files |
| `SweepAuditOverride` | `TRUE` where the pin sweep replaced the engine's stem |
| `ReconstructedStemID_PreSweep` | The engine's stem before the pin sweep (not updated by the later helpers) |
| `SweepRollbackToPreSweep` | `TRUE` where the engine's stem was kept: the pin would have put two rows of one census on one stem, or `apply_terminal_to_host()` restored the engine's choice |
| `DP_Chunk`, `run_out_dir` | Chunk number and run folder (added by the drivers) |

`ReconstructionMethod`:

| Value | Meaning |
| --- | --- |
| `given` | Stem from the pin (`TrueStemID`): anchor row, pinned row before the anchor, row after the anchor with a pin, or row set by the pin sweep |
| `dp` | Assigned by the DP |
| `probabilistic` | Assigned by the probabilistic matcher |
| `provisional_dp` | Provisional anchor stem, for a tag whose anchor census has measured rows without a pin |
| `dp_mf_inferred` | Row of a census in which the stem was alive but not measured, set aside before the engine and given the stem of the one track that spans it (`ReconstructedStemID` `NA` when no single track does) |
| `pin_track` | Unpinned row moved onto the pin of its engine track (`apply_pin_track_rejoin()`) |
| `carried_terminal` | Unmeasured dead, stem-dead, broken-below or missing row given the stem of the previous record with its `StemID` (`apply_carried_terminal_backfill()`) |
| `given_orphan` | Unmeasured row of a stem never measured, given its `StemID` (`apply_orphan_stem_backfill()`) |
| `terminal_to_host` | Unmeasured dead or broken-below record moved to the stem that ended (`apply_terminal_to_host()`) |
| `bb_split`, `bb_split_carry` | Measured broken-below row that starts a new stem by rule R1, and the later rows of that stem |
| `bb_post_terminator_split`, `bb_post_terminator_split_carry` | First measured row after a stump that starts a new stem by rule R2, and the later rows of that stem |
| `measurement_rejoin` | Row of a stem joined by `apply_measurement_rejoin()` (BCI merge) |
| `none_after_anchor` | Row after the anchor census without a pin or an engine stem |
| `skipped_no_data` | Tag without usable measurements (drivers) |

### Posterior path files

With `posterior_samples > 0` every tag gets one file
`<out_dir>/posteriors/tag_<Tag>_posterior_samples_<BATCH_TS>_paths.<feather|rds|csv>`
(`BATCH_TS` empty gives a double underscore; `feather` needs `arrow`, else
`rds`). One row per distinct sampled reconstruction:

| Column | Content |
| --- | --- |
| `path_sig` | The per-sample stem label of every measured observation up to the anchor, in census order, joined by `-` |
| `path_count` | Number of samples with this reconstruction; `path_count / n_samples` is its posterior probability |
| `path_prob` | Reference only: for the DP each sample is re-weighted by its own probability, which counts that probability twice (≈ p²); for the matcher it equals `path_count / n_samples` |
| `recon` | `obs_row_id:label` pairs, separated by `;` |

The labels are per-sample labels, so one grouping of observations can appear
in several rows; compare a path with the export by grouping, not by label. To
sample a tree, draw a row with weight `path_count`. Trees of both engines can
be sampled together this way: drawing a row by `path_count` is drawing one of
the tree's samples uniformly.

- DP samples are drawn by backward sampling from the exact posterior; likely
  reconstructions repeat.
- Matcher samples are approximate (perturbed assignments per census pair,
  stitched and repaired); with many stems nearly every sample is unique
  (`path_count = 1`). Their spread measures identity uncertainty
  heuristically, not as a calibrated posterior.

The engines stage their raw samples in `<out_dir>/posteriors/.staging/`. After
the post-engine helpers, `finalize_posterior_paths()` translates the labels
with the renumbering, applies the pins (`apply_pins_to_samples()`: pin sweep
with duplicate rollback, then the track rejoin) and the broken-below rules
(`apply_bb_invariants_to_samples()`) to each sample, writes the path file and
deletes the staging file. Only pinned links change, so the spread of the
unpinned links is kept and the export's pin treatment holds in every sample.

The exported `ReconstructedStemID` is the MAP reconstruction (DP) or the most
representative sample (matcher), after the pin sweep and the helpers. With a
finite number of samples the exported grouping can be missing from the path
file.

### Files written by the drivers

In the run folder (`dp_global/output/<run>/` for the generic drivers):

- `run_started.txt`, `run_finished.txt`, `run_log.txt`, `run_parameters_full.txt`
  (settings, command line and parameter structure);
- `stem_reconstruction_dp_global_rcpp.csv` / `.rds` / `.feather` / `.pdf`
  (single-pass run) or `stem_reconstruction_dp_global_rcpp_chunk_NNN.*` with
  `_done.txt` or `_failed.txt` markers (chunked runs);
- `posteriors/` (path files; `.staging/` is empty after a completed run);
- `main_cpp.R` only: `bio_pars_report.pdf`, realism tables
  (`RUN_REALISM_REPORT`), sensitivity sweeps (`SENSITIVITY_MODE`) and the
  k-sweep plot (`RUN_K_SWEEP_DEMO`).

## What it does

### Overview

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
5. The pin sweep sets pinned rows to their pin. Rows after the anchor keep
   their pin or stay unassigned.
6. The post-engine helpers fix rows the pin sweep left behind, attach
   unmeasured terminal records, apply the broken-below rules, number the stems
   1..N and write the posterior path files.

### Tracks, states and phases

The DP assigns the measured observations of each census to K tracks (stem
slots), at most one observation per track. A **state** of a census is such an
injective assignment; a census with n observations has

$$P(K, n) = \frac{K!}{(K-n)!} = K (K-1) \cdots (K-n+1)$$

states. Each track has a **phase**: 0 before its stem appears, 1 while the
stem is alive, 2 after it died (or when the track is never used). Allowed
changes are 0 → 0, 0 → 1 (recruitment), 1 → 1 (growth), 1 → 2 (death) and
2 → 2. A track never returns from 2 to 1: a stem measured, then not measured,
then measured again is not one stem, unless the census in between is detected
as missing from the field (it is then left out of the DP, below).

The tracks of the anchor census carry the pins of its stems, so the
reconstruction is expressed in anchor identities: an earlier observation on an
anchor track belongs to that stem; an observation on another track belongs to
a stem that ends before the anchor.

### Number of tracks K

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
  anchor stem is no larger than `Bio_Recruit_MaxDBH_unit` (a tree of only large
  stems is unlikely to have had a recent recruit).
- **Main-rule breaks**: with the main-stem rule (ForestGEO codes) each census
  with a break after a census with a main-stem record adds one track, because
  the broken record may be a new stem.
- **Pins before the anchor**: every pinned stem that ends before the anchor
  needs its own track.

When K is smaller than the largest stem count of a census, the tag goes to the
matcher (`K_too_small`).

### State-space limit (`max_states`)

`max_states` (driver setting `DP_MAX_STATES`) bounds the DP in two ways:

- **per census** (`enum_exceeded`): a census with more allowed states than
  `max_states`. When pins or growth bounds restrict the tracks of a census,
  only the allowed assignments are counted (`enumerate_states_constrained()`,
  which returns `NULL` past `max_states`, never a truncated set);
- **per census pair** (`edge_count_exceeded`): states of a census × reachable
  states of the next census ≥ `max_states²`.

In practice the per-census limit binds first. Largest number of measured stems
in one census that the DP solves when every track is allowed:

| `max_states` | K = n + 1 | K = n + 2 |
| --- | --- | --- |
| 1,000 | 5 (720 states) | 4 (360) |
| 10,000 | 6 (5,040) | 5 (2,520) |
| 40,000 | 6 (5,040) | 6 (20,160) |

Restrictions from pins and pruning often leave far fewer states than P(K, n),
so tags above these counts can still be solved. Defaults: 50,000 in the DP
function, 40,000 in `main_cpp.R` and `main_cpp_chunk.R`, 10,000 in the BCI
driver. `estimate_dp_complexity()` (`R/complexity/`) gives the state counts of
every tag before a run.

### Biological model

| Process | Model |
| --- | --- |
| Growth | g = (D₁ − D₀)/Δt ~ Normal(μ(D₀), σ(D₀)²), μ(D) = α + γ log D (μ = α when γ = 0), σ(D) = max(σ₀ + σ₁ D, 10⁻⁶) |
| Mortality | P_death(D₀) = 1 − exp(−h₀ exp(β D₀) Δt) |
| Recruitment | P_recruit = 1 − exp(−λ Δt) for an empty track; the new stem's DBH ~ Lognormal(meanlog, sdlog) |

Δt is in years; probabilities are clamped to [10⁻¹², 1 − 10⁻¹²].

**Measurement error** (`use_measurement_error = TRUE`; Chave et al. 2004): a
recorded DBH is D_obs = D_true + ε with

$$\varepsilon \sim (1-p)\,\mathcal{N}(0, \mathrm{SD1}(D)^2) + p\,\mathcal{N}(0, \mathrm{SD2}^2), \qquad \mathrm{SD1}(D) = aD + b$$

(`meas_sd1_a` = 0.0062, `meas_sd1_b` = 0.0904 cm, `meas_sd2` = 4.64 cm,
`meas_p_big` = 0.05). The observed growth is then a mixture of four normals
(small or large error at each census, weights (1−p)², (1−p)p, p(1−p), p²),
each with variance σ(D₀)² + (SD_a² + SD_b²)/Δt², evaluated with log-sum-exp.

**DBH recorded in classes**: at a census listed in `dbh_round_censuses`, a DBH
d below `dbh_round_max` (default 5.5 cm) was rounded down to a class of width
w = `dbh_round_width` (default 0.5 cm), so the true size is in [d, d + w). For
a link with rounding flags r₀, r₁ ∈ {0, 1}:

$$g_{\text{mid}} = g + \frac{(r_1 - r_0)\,w/2}{\Delta t}, \qquad \sigma_{\text{eff}}^2 = \sigma(D_0 + r_0 w/2)^2 + \frac{(r_0 + r_1)\,w^2/12}{\Delta t^2}$$

and the growth likelihood uses g_mid, μ(D₀ + r₀w/2) and σ_eff (also in every
component of the measurement-error mixture). The hard growth bounds stay on
the recorded values: rounding changes how likely a link is, not which links
are allowed. Without this rule a one-class step (e.g. 10 → 15 mm) read as
growth is so unlikely under a narrow growth distribution that one stem would
be split into a death and a recruit. BCI recorded stems below 55 mm in 5 mm
classes in 1982 and 1985 (CTFS R package growth tutorial; Piponiot et al.
2024); the matcher applies the same rule.

### Transition cost

The cost of a transition between a state at census t and a state at t+1 is the
negative log-likelihood, summed over the K tracks, plus a tie-break. It is
computed in C++ (`transition_cost_paired_rcpp()`, one call per census pair
for all candidate pairs); `transition_cost_tracks_bio_components()` gives the
same cost split into its terms.

| Track at t → t+1 | Cost |
| --- | --- |
| A. empty → empty | −log(1 − P_recruit) |
| B. empty → DBH (recruitment) | −log P_recruit − log f_Lognormal(D₁); `hard_penalty` when D₁ > `Bio_Recruit_MaxDBH_unit` or D₁ ≤ 0 |
| C. DBH → empty (death) | −log P_death(D₀) |
| D. DBH → DBH (growth) | steps D1–D7 below |

Case D, with g = (D₁ − D₀)/Δt:

1. **D1** g < `Bio_Max_Shrink`: `hard_penalty` (10⁶); next track.
2. **D2** g > `Bio_Max_Growth`: `hard_penalty`; next track.
3. **D3** σ(D₀) = max(σ₀ + σ₁ D₀, 10⁻⁶).
4. **D4** μ(D₀) = α + γ log D₀ (α when γ = 0 or D₀ ≤ 0).
5. **D5** growth likelihood: without measurement error
   (g − μ)²/(2σ²) + log σ + ½ log 2π; with it, −log of the four-component
   mixture (rounded censuses: g_mid and σ_eff).
6. **D6** soft shrinkage: + k_shrink (D₀ − D₁)² when D₁ < D₀.
7. **D7** soft growth: + k_growth (D₁ − D_cap)² when D₁ > D_cap =
   D₀ + `Bio_Max_Growth_Soft` × Δt.

**Tie-break** (D8): + `eps_tiebreak` (10⁻⁶) × Σ |rank_t(k) − rank_t+1(k)| over
the tracks measured at both censuses (ranks of DBH among the measured tracks
of each census). Among otherwise equal costs it favours the assignment that
keeps the size order of the stems.

A recruitment is costed only on an empty track (phase 0 → 1), a death only on
a track that leaves phase 1, so a stem that keeps growing is never charged a
death and a recruit.

### Anchor census

`anchor_start` is the census whose identities are trusted (BCI: 7 = 2010, the
first census that tagged every stem). The DP solves the censuses up to it;
rows after it are not solved.

1. A tag whose last measurement is before `anchor_start` is anchored at its
   last measured census.
2. When the tag is measured after the anchor but the anchor census has no
   measured stem, the anchor moves forward to the first later census with a
   measured stem (one with a pin, when there is one).
3. When the anchor census has no row with a DBH or a pin, the latest earlier
   census with both is used.
4. When no census has both and `allow_provisional_anchor = TRUE`, the last
   measured census becomes a provisional anchor: its measured rows get new
   IDs (`provisional_dp`). The same happens when the anchor census has
   measured rows without a pin.
5. Otherwise the tag goes to the matcher (reasons `anchor_*`,
   `no_obs_up_to_anchor`).

Rows after the anchor: a row with a pin takes it (`given`); unmeasured
broken-below / resprout-coded rows that directly follow the anchor (up to the
first later measurement) take the anchor stem when exactly one stem is
measured at the anchor; the other rows stay `none_after_anchor`.

### DP recursion, MAP and posterior

Let s_t be a state of census t, c(s_t, s_t+1) the transition cost and A the
anchor census, whose state is fixed by the pins. Going back from A census by
census, the DP:

1. enumerates the allowed states of census t (pins and growth bounds restrict
   the tracks an observation may use);
2. pairs each state with the reachable states of census t+1 and removes the
   pairs that break the phase rules, the code rules or the pruning bounds
   (C++ `derive_phase_prev_batch_rcpp()`);
3. costs the remaining pairs and keeps, for each state, the cheapest
   continuation and the log-sum-exp of all continuations:

$$V_t(s) = \min_{s'} \left[ c(s, s') + V_{t+1}(s') \right], \qquad \log Z_t(s) = \log \sum_{s'} \exp\!\left( -\frac{c(s, s')}{\tau} + \log Z_{t+1}(s') \right)$$

with V_A = 0 and log Z_A = 0 at the anchor state.

- **MAP**: the reconstruction decoded from the cheapest start state through
  the back-pointers (among equal costs the first enumerated state wins). This
  is the exported `ReconstructedStemID` (`dp`).
- **Posterior**: P(path) ∝ exp(−TotalCost / τ), τ = `temperature` (1 in the
  drivers: the costs are negative log-likelihoods, so τ = 1 is the Bayesian
  posterior; τ < 1 sharpens it, τ > 1 flattens it). A forward pass gives the
  marginal probability of every observation on every track.
- **Samples**: `posterior_samples` full reconstructions drawn by backward
  sampling from the same weights (seed `posterior_sample_seed`).

**Interchangeable tracks.** A track with no stem at the anchor and no pinned
observation carries no identity, so a reconstruction with m stems on S such
tracks occurs under S!/(S−m)! labelings of equal cost. The MAP compares single
labelings and is not affected, but summing over labelings would weight each
reconstruction by its number of labelings and make the marginals depend on the
number of spare tracks. With option `dp_global.posterior_label_free = TRUE`
each step back that puts a more interchangeable tracks into use, while u are
already in use, is divided by (S−u)!/(S−u−a)!, so each reconstruction counts
once in the marginals and samples.

**Posterior summaries**:

- `DP_PosteriorReconstructedProb`: for a row on an anchor or pinned track, the
  marginal of that track; for a row on an interchangeable track with the
  label-free posterior, the share of samples in which the row's stem ends at
  the same observation as in the export (`NA` without samples).
- `DP_PosteriorTop{k}ID` / `Prob`, `DP_PosteriorEntropy`: per-track
  marginals.
- `DP_PosteriorUnlinkedProb`: posterior mass of the row on tracks without an
  anchor stem.
- `add_dp_posterior_bins()`: `unlinked-likely` when the unlinked probability
  is at least `unlinked_prob` (0.5); otherwise `confident` when
  `DP_PosteriorReconstructedProb` is at least `confident_prob` (0.95);
  otherwise `ambiguous`.

### Pruning

Before costing, the DP removes the transitions outside conservative bounds,
so that the expensive likelihood is evaluated only for possible links
(`prune_hard = TRUE`). The bounds are built in layers:

1. **User bounds**: `prune_min_growth` / `prune_max_growth` (cm/yr; default
   `min_growth` / `max_growth`).
2. **Biological bounds**: with `prune_use_bio_bounds = TRUE` (function
   default) they are narrowed to `Bio_Max_Shrink` / `Bio_Max_Growth`.
3. **Non-taper-corrected growth forms**: a tag whose `growth_form` is in
   `non_taper_corrected_growth_forms` (default palm, strangler fig, tree fern;
   the most common form of the tag) takes
   `non_taper_corrected_prune_min_growth` /
   `non_taper_corrected_prune_max_growth` (defaults −0.625 and 6.25) instead.
   These forms grow in diameter and cannot be taper-corrected, so a change in
   the height of measurement changes the recorded DBH.
4. **HOM widening**: for those forms, with a `hom` column (m), each census
   pair widens the bounds by `hom_tolerance_scale` (2 cm per m) ×
   max |HOM − 1.3| over the rows of the two censuses / Δt (`NA` HOM counts as
   1.3 m; 0 turns it off).
5. **Recruit size**: `prune_recruit_max_dbh` (default
   `Bio_Recruit_MaxDBH_unit`), capped at `Bio_Recruit_MaxDBH_unit` when
   `prune_use_bio_recruit = TRUE`.

A pruned transition is removed; a transition inside the prune bounds but
outside the `Bio_*` hard bounds is kept and costs `hard_penalty`. The bounds
also restrict the state enumeration: a track holding an anchor stem is open to
an observation only if an observation of the next census that may use that
track can follow it (twice the non-taper bounds for those forms, leaving room
for the HOM widening). A pair of censuses with a non-positive interval is not
pruned. The drivers pass `PRUNE_BOUND_FACTOR` × the fixed growth bounds and
recruit limit with `prune_use_bio_bounds = FALSE` and
`prune_use_bio_recruit = FALSE`, and the same values as non-taper bounds, so
all tags share one wide window and the HOM widening is what differs for
non-taper forms.

When pruning leaves no feasible transition, the tag goes to the matcher
(`no_states_produced`, `no_feasible_edges`, …). `attr(out, "DP_PruneInfo")`
records the transitions examined and removed (per census pair) and the
effective bounds.

### ForestGEO field codes

Resprout codes `R`, `RP`, `RF`, `RT`, `QR` and `OR` in `ListOfTSM` (whole
codes, regex `\b(R|RP|RF|RT|QR|OR)\b`) and `Status == "broken below"` mark a
break. A **resprout observation** is such a row with a DBH; a **stump** is a
broken-below row without a DBH.

**R-recruit constraint** (rule R1, option `dp_global.bb_dbh_rule = "r1"`, the
engine default). A resprout observation at census t+1 must continue a track
occupied at census t: in the engine it is the last record of the stem that
broke, so it needs no track of its own. Enforced in
`derive_phase_prev_batch_rcpp()`; after the engine, R1 starts a new stem at
it (post-engine helpers).

**Main-stem rule** (option `dp_global.bb_dbh_rule = "main_ends"`, BCI driver
`BB_DBH_RULE`). Up to 2010 the BCI field writes a break on the tree's
main-stem record (`PrimaryStem == "main"`), with the DBH of the largest stem
still alive, or with no DBH when nothing is alive. Under this rule:

- the main stem of the census before ends at the break;
- the broken record with a DBH is matched like any other measurement, but it
  cannot continue that main stem. With one stem before the break it is always
  a new stem; with several it continues another stem when one fits.

The DP forbids the main stem's track to hold an observation of the next
census (and adds one track per such break); the matcher masks the main stem's
links in that pair (`apply_main_end_mask()`); after the engine R1 only checks
the rule. From 2015 broken-below records in BCI carry no DBH, so the rule has
no effect after 2010.

**Resprout segment split** (`allow_segment_split = TRUE`, the default of the
DP function; the BCI driver sets `DP_ALLOW_SEGMENT_SPLIT = FALSE`). The DP
solves the tag as two independent problems at the first census before the
anchor (not the tag's first census) that holds a resprout observation:

- the pre-segment, up to and including that census, anchored at its last
  measured census (provisional IDs where it has no pin);
- the post-segment, from the next census to the anchor, in which every stem of
  the first census is a new stem. When that census has no more stems than the
  next, every one of them must continue to the next census (recruit
  continuity; skipped where it would leave no feasible transition).

Pre-segment IDs are shifted above the post-segment IDs. Each segment draws
its own samples; because every stem ends at the split, draw k of one segment
paired with draw k of the other is an exact draw of the whole tag
(`pair_segment_posterior_samples()`; the pre-segment uses
`posterior_sample_seed + 1`). When the whole tag exceeds `max_states`, both
segments use the matcher. Segments do not split again: a later resprout
observation follows the R-recruit constraint or the main-stem rule.

**NA-R barrier**. A census before the anchor with no measured stem and an
unmeasured break record means the whole tree died back: no identity crosses
it. After decoding, a track with rows on both sides is cut and its unpinned
rows before the barrier take a new ID, unless the track's pin is known after
the barrier and no row before carries another pin. The unmeasured break rows
of the barrier census take the stems of the census before, in row order.

**Missing from the field**. Rows of a stem that was alive but not measured are
set aside before the engine, so the gap is not read as a death plus a
recruitment:

- explicit: a row without a DBH with code `MF` or status `missing`, extended
  through the following consecutive censuses in which no row of the tag has a
  DBH;
- implicit: a census in which no row of the tag has a DBH, between the tag's
  first and last measured census, unless a row of it carries a break code.

A census left without rows is dropped from the DP, so tracks link the
censuses on both sides directly. After the engine each set-aside row takes the
one track assigned in the nearest censuses before and after and free at its
census (`dp_mf_inferred`); with several candidates it stays `NA`.

### Pins

A pin (`TrueStemID`) is a known stem identity: at the anchor census, and at
any earlier row the driver pins (BCI: field tags, database IDs from 2010 and
database IDs propagated to earlier records; the pin steps are in
`BCI_stem_reconstruction/2_STEM_IDENTIFICATION/README.md`).

- **While solving** (`pin_truestemid = TRUE`): the DP enumerates only the
  states that put each pinned observation on its pin's track; the matcher
  constrains every sampled pair (below).
- **Pin sweep**: every pinned row is set to its pin (`given`); the engine's
  choice is kept in `ReconstructedStemID_PreSweep` and changed rows are
  flagged in `SweepAuditOverride`. The sweep runs inside the DP's fallback
  path, in the matcher and in the drivers' `run_dp_one_group()`, which is the
  only sweep on the normal DP path. The drivers' sweep processes rows in order
  and keeps the engine's stem when the pin would put two rows of one census on
  one stem (`SweepRollbackToPreSweep`).
- **Hard pins** (option `dp_global.hard_pins = TRUE`, BCI driver
  `HARD_PINS`): pins also hold after the engine. The broken-below rules never
  cut between two rows with the same pin (when the rows they would move hold
  the current ID as their pin, the other rows take the new ID instead) and
  never rewrite a pin, and `apply_terminal_to_host()` leaves pinned records in
  place.
- **Field tags**: the `StemTag`s from census
  `dp_global.field_tag_first_census` on (default `Inf`, none; BCI 7). A stump
  that carries a field tag is that stem with no measurable DBH at 1.3 m, so it
  stays with its tag's stem.
- **Exception** (`pin_prelife_rows()`): an unmeasured dead, stem-dead,
  broken-below or missing record that comes before any life of its own pin is
  not bound by the pin when it carries no field tag, or when its pin never
  lives. Such a record is the end of the stem that broke (in BCI before 2010
  the database stores a break under the new StemID of the resprout measured
  later), so `apply_terminal_to_host()` can move it there.

### Fallback to the probabilistic matcher

| `DP_FallbackReason` | Cause |
| --- | --- |
| `growth_form_forced` | A `growth_form` of the tag (rows up to the anchor) is in `fallback_growth_forms` |
| `species_forced_probabilistic` | A `Species` value of the tag is in `prob_species` |
| `enum_exceeded` | A census has more allowed states than `max_states` |
| `edge_count_exceeded` | Two adjacent censuses reach `max_states²` state pairs |
| `K_too_small` | K is smaller than the largest stem count of a census |
| `no_obs_up_to_anchor`, `anchor_missing_obs`, `anchor_missing_truestem`, `anchor_missing_truestem_prov_disabled`, `anchor_ids_missing`, `anchor_truestem_not_found` | No usable anchor census |
| `no_states_produced`, `no_reachable_next_states`, `no_feasible_edges`, `forward_edges_missing`, `forward_no_alpha` | No feasible state or transition after pruning and the code rules |
| `decode_failure`, `viterbi_decode_failure`, `assign_mismatch`, `next_assign_row_mismatch` | Internal decoding failure |
| `error:<message>` | The DP raised an error (drivers: first 200 characters) |

`do_fallback()` passes the matcher the DP's growth bounds, prune bounds and
`use_bio_hard_*` flags. After the matcher it applies the NA-R barrier (every
crossing track) and, with rule R1, `split_live_resprout_tracks()`: a measured
broken-below or resprout row starts a new stem, so only the broken stem's
track is split and its rows before the break take a new ID (the break row
keeps the matcher's ID, which the pin sweep ties to the anchor). The tree's
other stems keep their links, and a track with a pinned row before the break
is not split. With the main-stem rule the matcher already ends the main stem,
so this split is not applied. When the DP raises an error, the drivers run the
matcher on the rows up to the anchor in the same way and label the later rows
(`given` / `none_after_anchor`).

### Probabilistic matcher

`match_stems_probabilistic()` works census pair by census pair and stitches
the pairs into full reconstructions:

1. **Link scores** (`compute_pairwise_log_likelihood()`): for every stem of
   census t and every stem of t+1, the growth log-likelihood of the DP model
   (no measurement-error mixture; DBH classes as in the DP), the log survival
   probability of the earlier stem and the soft shrinkage penalty, plus a soft
   growth penalty above `Bio_Max_Growth`. A link outside the growth bounds,
   or outside `Bio_Max_Shrink` / `Bio_Max_Growth` when
   `use_bio_hard_shrink_in_prob` / `use_bio_hard_growth_in_prob` are `TRUE`,
   is forbidden.
2. **Assignment matrix** (`augment_cost_matrix()`):
   - birth-death (`birth_death = TRUE`, default): one death cell per earlier
     stem (log P_death) and one recruitment cell per later stem (log P_recruit
     plus the log size density), so every combination of survivals, deaths and
     recruitments is an assignment. A link is taken only when it is more
     likely than the death of the earlier stem plus the recruitment of the
     later one, as in the DP;
   - count-based (`FALSE`): only as many death and recruitment slots as the
     stem counts and the forbidden links require, so with a stable stem count
     no stem can die and none can be recruited.
3. **Sampling**: `n_samples` (200) assignments per pair, each the best
   assignment of the pair after adding Gumbel noise of scale `temperature`
   to the event cells: solved exactly in birth-death mode
   (`hungarian_min_rcpp()`), greedily with a feasibility repair in
   count-based mode (`enforce_feasible_assignment()`). Adding Gumbel noise and
   taking the best assignment approximates sampling from the Gibbs
   distribution over assignments.
4. **Pins while sampling**: `apply_pin_mask()` sends an observation pinned to
   an anchor stem to that stem's column; `apply_track_pin_mask()` covers every
   pin, also of stems that end before the anchor: a pinned observation must
   join the observation that carries its pin (its own or one inherited from a
   later census of the same sample) and may not join one that carries another
   pin. When the masks forbid every link a count-based pair relied on,
   `pin_masked_pair()` rebuilds the pair with enough death and recruitment
   slots, so a pin never forces a link. With the main-stem rule the main
   stem's links into a census with a break are masked in the same step.
5. **Look-ahead** (`prob_lookahead_weight` > 0, at least three censuses and
   K ≥ 4): when stitching backward from the anchor, the matrix of pair
   (t, t+1) gets a continuity bonus for the stems whose assignment at t+1
   continues plausibly to t+2 (best column 0, others down to −2 log units,
   times the weight; function default 0.5, drivers 1).
6. **Stitching** the pairs backward from the anchor into full samples.
7. **Repair** (`repair_stitched_growth_violations()`): links outside the hard
   growth rates are cut (pinned observations are never cut off; a link
   between two pinned observations is kept). With `n_sigma_me` finite and
   `use_bio_hard_shrink_in_prob = TRUE`, a run of decreases whose cumulative
   shrinkage exceeds `n_sigma_me` × √(SD(d_start)² + SD(d_now)²), SD(D) =
   0.0062 D + 0.0904 cm, is cut at its start (default 3; the BCI driver sets
   `Inf`, off).
8. **Pin filter** (`filter_pin_consistent_samples()`): samples that break a
   pin to an anchor stem are dropped; all are kept (with a warning) when fewer
   than min(10, `n_samples`/4) would remain.
9. **Marginals** (`compute_marginals_from_samples()`): Top-k, entropy and the
   share of samples per track for every observation.
10. **Export** (`select_consensus_trajectory()`): the sample whose links agree
    most with the links of the other samples (maximum expected accuracy; ties
    broken by how often its whole partition occurs). With option
    `dp_global.prob_export_rule = "fewest_me"` the samples with the fewest
    links that shrink beyond measurement error (CTFS rule:
    D₁ ≤ D₀ − 4 (0.006214 D₀ + 0.9036), mm) are considered first. The export
    is always a sampled reconstruction.
11. A diagnostic count of the hard-rate violations left in the export
    (reported only; the posteriors are not changed), the pin sweep, and the
    staging of the samples for the path file.

`DP_PosteriorReconstructedProb` is the share of samples that give the row its
exported stem (1 at the anchor); `DP_PosteriorUnlinkedProb` stays `NA`.

#### DP and matcher compared

| | DP | Matcher |
| --- | --- | --- |
| Scope | Exact over the whole tag | Per census pair, stitched |
| Hard bounds | Pruned before costing; inside the prune bounds, outside `Bio_*`: `hard_penalty` | Forbidden links (prune bounds; `Bio_*` with the `use_bio_hard_*` flags) |
| Runs of shrinkage | Every step pays its growth likelihood, so runs add up | Optional cumulative-shrinkage cut |
| Survival | Mortality cost on the track that dies | Log survival on every link |
| Soft growth penalty | Above `Bio_Max_Growth_Soft` | Above `Bio_Max_Growth` (acts only when the growth gate is off) |
| Measurement error | Four-component mixture (optional) | Not modelled in the link scores |
| Posterior | Exact marginals and samples | Heuristic spread of perturbed samples |

### Post-engine helpers

The drivers apply, to each chunk and in this order (all in
`R/dp_global_main.R`):

1. **`maybe_add_posterior_bins()`**: `DP_PosteriorBin`.
2. **`apply_pin_track_rejoin()`**: the pin sweep moves pinned rows only, so an
   unpinned stretch of rows that the engine had linked to them can be left on
   another ID (e.g. a stem renumbered by the database before 2010, with only
   the later StemID pinned). Such an ID moves onto the pin of its engine track
   when the track segment holds one pin value, all measured rows of the ID
   belong to the segment, none of its rows is pinned and the pin has no row in
   its censuses (`pin_track`). Segments end at every broken-below row with a
   DBH (with the main-stem rule only where it would continue the main stem)
   and after every NA-R barrier, so those splits stay.
3. **`apply_carried_terminal_backfill()`**: an unmeasured dead, stem-dead,
   broken-below or missing row without a stem takes the stem of the previous
   record with the same `StemID` (`OriginalStemID` without `StemID`)
   (`carried_terminal`).
4. **`apply_orphan_stem_backfill()`**: an unmeasured row of a stem never
   measured, without a stem or pin, takes its `StemID` (`given_orphan`).
5. **`apply_terminal_to_host()`**: an unmeasured dead, stem-dead or
   broken-below record that comes before any life of its stem (or sits on a
   stem that never lives) moves to the stem that ended (`terminal_to_host`):
   - (a) the engine's own choice before the sweep, when that stem was last
     alive before the record (if a later step relabelled that track, the
     current ID of its last measured row before the record);
     `SweepRollbackToPreSweep` is then set;
   - (b) else the one stem of the tag last alive in the census before, when
     exactly one such record group starts in that census;
   - (c) else, when the record is the tree's main record, the main stem of the
     census before (up to 2010 a break is written on the main stem's record,
     so when the whole tree broke the stump is the main stem's end);
   - otherwise it stays (a stem that is never alive).

   A missing record of the same StemID that follows moves with it. A record
   never moves onto a stem with a row in one of its censuses; measured rows
   never move; each record keeps its own `Status`.
6. **`apply_broken_below_invariants()`**, per reconstructed stem (Tag ×
   `ReconstructedStemID`) in census order:
   - **R1**: a broken-below row with a DBH that has earlier rows on its stem
     starts a new stem; the new ID is carried to the later rows until the next
     such row (`bb_split`, `bb_split_carry`). With the main-stem rule, only
     where it would continue the main stem of the census before (a check: the
     engines do not make such links).
   - **R2**: after a stump, the first measured row of the stem starts a new
     stem, carried forward the same way (`bb_post_terminator_split`,
     `bb_post_terminator_split_carry`); unmeasured rows between keep the old
     ID.

   Dead and stem-dead rows trigger neither rule. New IDs are taken above the
   largest ID of the table. Without hard pins, a relabelled row whose
   `TrueStemID` equals the replaced ID gets the new ID as its `TrueStemID`.
7. **`renumber_engine_minted_ids()`**: numbers the stems of each tag 1..N by
   first census, then by decreasing DBH in that census, then by old ID. The
   mapping is applied to `ReconstructedStemID`, `ReconstructedStemID_PreSweep`
   and the Top-k IDs, and to `TrueStemID` only where the engine minted it
   (`provisional_dp`, and without hard pins the R1/R2 rows); database pins are
   not renumbered.
8. **`finalize_posterior_paths()`**: writes the posterior path files (see
   Outputs).

### Measurement rejoin

`apply_measurement_rejoin()` (`R/measurement_rejoin.R`) runs later, in the
BCI merge. When the point of measurement of a trunk moves (buttresses,
deformities, strangler figs), the next DBH can fall outside the hard growth
bounds; the engines then end the stem and start an impossible recruit (above
the recruit limit) for the same trunk. The rejoin joins the two stems when:

- nothing marks a break or a death between them, and the end and start are
  clean (no census gets two rows);
- the size change is that of a moved point of measurement (outside the hard
  bounds; a smaller DBH down to a taper-corrected ratio of 0.4, a larger one
  only with a recorded HOM change, ratio at most 1.5);
- the tree shows that its point of measurement moves: a HOM away from 1.3 m,
  a buttress code, a growth form without a stable point of measurement, or a
  buttressed species in censuses without buttress codes;
- the two stems carry at most one pin, and the pair is clearly the best one.

The joined stem keeps one of its IDs (`measurement_rejoin`);
`apply_measurement_rejoin_to_paths()` makes the same joins in the posterior
samples. The file header numbers and justifies every rule.

## Parameter estimation

`estimate_bio_pars()` (`R/dp_global_bio.R`) estimates the parameters of one
parameter set from stems of known identity: rows with `DBH > 0` and a
`TrueStemID`, on a grid of every stem at every census from
`anchor_start_census` on. The interval of a census pair is computed per stem
from its two dates (a missing date takes the mean date of its tag in that
census, else of the census). Each guardrail and soft penalty has a source,
`"data"` (estimated) or `"fixed"` (`*_fixed`). The function stops when no row
is usable, or fewer than two censuses or five growth pairs remain.

**Growth.** Annual increments g = (D₁ − D₀)/Δt of consecutive measurements.

- Mean: `lm(g ~ log(d0))` with at least 10 pairs and variation in log D₀
  (α, γ); otherwise α = mean(g), γ = 0.
- SD: per pair, |g − μ(D₀)| × √(π/2) (E|X| = σ√(2/π) for a normal). With
  measurement error, the error variance of the increment,
  [Var(ε|D₀) + Var(ε|D₁)]/Δt² with Var(ε|D) = (1−p) SD1(D)² + p SD2², is
  subtracted in quadrature. `lm(sd ~ d0)` then gives σ₀ (at least 0.01) and
  σ₁; when the slope is negative (the SD falls with size, e.g. palms) σ₀ is
  the mean SD and σ₁ = 0, so the SD cannot reach zero for large stems.
- `enforce_growth_bounds = TRUE` drops increments outside
  `growth_min_fixed` / `growth_max_fixed` before the fit (the drivers do, with
  the fixed growth bounds).

**Hard bounds.** `max_shrink` = the `shrink_data_quantile` (0.001) quantile
of g; with measurement error, the smaller of it and the `shrink_hard_prob`
(10⁻⁴) quantile of the measurement-error distribution of an increment.
`max_growth` = the `growth_data_quantile` (0.999) quantile; with measurement
error, the larger of it and the 1 − `growth_hard_prob` quantile.
`max_growth_soft` = min(`max_growth`, the `growth_soft_quantile` (0.99)
quantile). The drivers fix `max_shrink` and `max_growth` and keep the data
quantiles for the soft cap.

**Soft penalties.** k = 1/(2s²), so that a change of s cm costs about one
unit: with measurement error (at least 10 pairs) s is the median measurement
SD of a DBH difference; otherwise s² is the variance of the observed decreases
(`k_shrink`) or increases (`k_growth`), with k = 50 below five events. k is
clamped to [10⁻⁶, 10⁶]; the drivers fix both at 0.

**Mortality.** A stem measured at t₀ and not at t₁ died. h₀ and β maximise

$$\sum_i \left[ d_i \log P_{\text{death}}(D_{0,i}) + (1 - d_i) \log\left(1 - P_{\text{death}}(D_{0,i})\right) \right]$$

(`optim`, BFGS, on log h₀ and β; start `mortality_start` = (log 0.01, 0)).

**Recruitment.**

- Size: lognormal fitted to the DBH of stems measured at t₁ and not at t₀
  (`MASS::fitdistr`; meanlog = log 2, sdlog = 0.5 with fewer than two
  recruits). `enforce_recruit_max = TRUE` drops recruits above
  `recruit_max_fixed` before the fit.
- Limit `recruit_max_dbh`: the `recruit_max_quantile` (0.999) quantile of the
  recruit DBH (5 cm without recruits) or `recruit_max_fixed`.
- Rate (`recruit_rate_unit`):
  - `"tree"` (default): new stems per established tree per year,
    λ = (n_new + 0.5)/Σ Δt, over the trees with a measured stem at the start
    of an interval; n_new counts their stems without a DBH at the start and
    with one at the end. This is the probability per year that a tree gains a
    stem, which is what both engines use λ for (the DP's empty track, the
    matcher's recruitment cell);
  - `"slot"`: recruits per empty grid cell per year. The empty cells are
    mostly dead stems and stems that recruit later, so this rate barely
    varies between species. It is also used when no tree is established.

| Quantile argument | Default | Gives |
| --- | --- | --- |
| `shrink_data_quantile` | 0.001 | `max_shrink` |
| `growth_data_quantile` | 0.999 | `max_growth` |
| `growth_soft_quantile` | 0.99 | `max_growth_soft` |
| `recruit_max_quantile` | 0.999 | `recruit_max_dbh` |
| `shrink_hard_prob`, `growth_hard_prob` | 10⁻⁴ | Measurement-error tails |

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

Arguments of `match_stems_dp_global_backward_marginals_batch()` (function
defaults):

| Argument | Default | What it does |
| --- | --- | --- |
| `anchor_start` | | Anchor census |
| `max_states` | 50,000 | State-space limit |
| `max_tracks`, `slack_tracks` | 30, 1 | Cap on K and slack tracks |
| `slack_require_anchor_recruitable`, `slack_require_anchor_eps` | `FALSE`, 10⁻⁶ | Slack tracks only with a recruitable anchor stem |
| `min_growth`, `max_growth` | −Inf, Inf | Base growth bounds: pruning when no prune bound is given, matcher, `ConstraintViolation` |
| `prune_hard`, `prune_min_growth`, `prune_max_growth`, `prune_use_bio_bounds` | `TRUE`, `NULL`, `NULL`, `TRUE` | Growth pruning |
| `prune_recruit_max_dbh`, `prune_use_bio_recruit` | `NULL`, `TRUE` | Recruit-size pruning |
| `non_taper_corrected_growth_forms`, `non_taper_corrected_prune_min_growth`, `non_taper_corrected_prune_max_growth`, `hom_tolerance_scale` | palm, strangler_fig, tree_fern; −0.625; 6.25; 2 | Non-taper bounds and HOM widening |
| `pin_truestemid`, `allow_provisional_anchor` | `TRUE`, `TRUE` | Pins and provisional anchors |
| `allow_segment_split` | `TRUE` | Resprout segment split |
| `fallback_growth_forms`, `prob_species` | none | Tags sent to the matcher (vector or comma-separated string) |
| `use_measurement_error`, `meas_sd1_a`, `meas_sd1_b`, `meas_sd2`, `meas_p_big` | `FALSE`, 0.0062, 0.0904, 4.64, 0.05 | Measurement-error mixture |
| `dbh_round_censuses`, `dbh_round_max`, `dbh_round_width` | none, 5.5, 0.5 | DBH classes |
| `temperature`, `eps_tiebreak` | 1, 10⁻⁶ | Posterior temperature and tie-break |
| `posterior_top_k`, `posterior_samples`, `posterior_samples_format`, `posterior_samples_path`, `posterior_sample_seed` | 2, 0, `"rds"`, `NULL`, `NULL` | Posterior summaries and samples |
| `prob_n_samples`, `prob_lookahead_weight`, `prob_n_sigma_me`, `prob_birth_death` | 200, 0.5, 3, `TRUE` | Matcher |
| `use_bio_hard_shrink_in_prob`, `use_bio_hard_growth_in_prob` | `TRUE`, `TRUE` | Matcher hard gates |
| `verbose`, `chunk_id` | `FALSE`, `NULL` | Log output |

`Bio_*` columns enter the cost; `min_growth` / `max_growth` do not. The driver
values are in `scripts/README.md` and in the BCI stage-2 README.

## Checks and diagnostics

- The DP stops when `temperature` is not positive or a required `Bio_*` column
  is missing, `NA` or not constant within the tag;
  `estimate_bio_pars()` stops on unusable data or invalid fixed values.
- `attr(out, "DP_PruneInfo")`: transitions examined and removed (total and per
  census pair) and the effective bounds; `attr(out, "DP_Compute_Profile")`:
  number and time of the cost calls; `DP_Sampling_Profile` when samples are
  drawn.
- `ConstraintViolation`: reconstructed links outside `[min_growth,
  max_growth]`.
- `DP_FallbackReason`: why a tag went to the matcher.
- `estimate_dp_complexity()` and `get_tag_complexity_details()`
  (`R/complexity/estimate_dp_complexity_function.R`): states and transitions
  per tag, fallback flag and predicted run time, before a run.
- `transition_cost_tracks_bio_components()`: cost of a transition by term
  and track.
- `plot_tag_to_pdf()`: reconstructed trajectories per tag (with the
  `TrueStemID` panel when `include_reference = TRUE`).
- `R/sensitivity_transition_cost_bio.R` (`build_all_sweeps()`),
  `R/realism_calibration.R` (`realism_report_from_reconstruction()`),
  `R/k_tuning_viz.R` (join-versus-split cost as k varies).

Numerical safeguards: σ(D) and SD1(D) ≥ 10⁻⁶; probabilities in
[10⁻¹², 1 − 10⁻¹²]; log-sum-exp with the maximum subtracted; data-mode k in
[10⁻⁶, 10⁶]; stable ranks for the tie-break.

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

## References

- Chave J, Condit R, Aguilar S, Hernandez A, Lao S, Perez R (2004). Error
  propagation and scaling for tropical forest biomass estimates.
  *Philosophical Transactions of the Royal Society B* 359: 409–420.
- CTFS R package, growth tutorial:
  <https://ctfs.si.edu/ctfsdev/CTFSRPackageNew/index.php/web/tutorials/GrowthChange/index.html>
- Piponiot C, Condit R, Hubbell SP, Pérez R, Lao S, Aguilar S, Muller-Landau H
  (2024). Woody biomass stocks and fluxes in the Barro Colorado Island 50-ha
  plot (Appendix S1).
