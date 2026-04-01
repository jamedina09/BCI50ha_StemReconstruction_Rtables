# Global Dynamic Programming Stem-ID Reconstruction with Biological Costs

**Date:** 2026-03-30  
**Author:** José A. Medina-Vega

---

## Table of Contents

1. [Overview](#overview)
2. [Installation & Quickstart](#installation--quickstart)
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

---

## Overview

Reconstruct stable stem identities across censuses using a biologically informed global dynamic programming solver that enforces life-cycle constraints and provides uncertainty quantification.

### Key Features

- **Global optimization** across multiple censuses
- **Biological realism** (growth, mortality, recruitment)
- **Measurement error handling** following Chave et al. (2004)
- **Exact uncertainty quantification** via posterior marginals
- **C++ acceleration** for performance and fallbacks when DP is intractable

---

## Quickstart

### Prerequisites

R (packages: `data.table`, `igraph`; optional: `ggplot2`, `cowplot`) 

### Running the Code

**Interactive (single-tag or small datasets):**
```r
source("dp_global/scripts/main_cpp.R")
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
- To specify census timing, pass `interval_years` (scalar) or ensure a per-row interval column (e.g., `Bio_IntervalYears`) is present and set `interval_years = NULL` to enable automatic per-pair interval detection.
- Request `.rds` outputs by passing `--WRITE_DP_RDS=TRUE`.

### File Structure

```
dp_global/
├── README.md
├── R/
│   ├── dp_global_main.R              # Source loader: sources all R modules below
│   ├── dp_global_bio.R               # Biological parameter estimation (estimate_bio_pars)
│   ├── dp_global_dp.R                # Production DP solver (match_stems_dp_global_backward_marginals_batch)
│   ├── dp_global_states.R            # State enumeration and track-DBH helpers
│   ├── dp_global_matchers.R          # Fallback igraph matcher (match_stems_optimal_backward)
│   ├── dp_global_utils.R             # Shared utilities (interval resolution, etc.)
│   ├── dp_global_diag.R              # Diagnostics and PDF plotting
│   ├── naming_helpers.R              # Output directory naming helpers
│   ├── sensitivity_transition_cost_bio.R  # Sensitivity sweep helpers
│   ├── realism_calibration.R         # Realism calibration helpers
│   ├── k_tuning_viz.R                # K-tuning visualisation
│   ├── check_functions.r             # Mortality parameter inspection and plotting utilities
│   ├── complexity/
│   │   ├── estimate_dp_complexity_function.R  # DP complexity estimator
│   │   ├── test_complexity_estimator.R        # Complexity estimator test/demo
│   │   ├── validate_run.R                     # Run validation (compare two outputs)
│   │   ├── validate_tiebreaker.R              # Tie-breaker determinism validation
│   │   ├── validate_m_implementation.R        # M-code pinning validation
│   │   └── check_symmetry.R                   # Cost symmetry checks
│   └── dpglobal_bundle/
│       ├── dpglobal_bundle_loader.R      # Bundle builder (creates RData + manifest)
│       ├── package_bundle.sh             # Packaging helper (creates tarball in dist/)
│       ├── verify_bundle.R               # Basic smoke-test for deployed bundles
│       └── dist/                         # Generated tarballs (not tracked by git)
├── scripts/
│   ├── main_cpp.R                    # Interactive / single-tag driver
│   ├── main_cpp_chunk.R              # Chunked driver for large runs
│   └── main_cpp_bci.R               # BCI debug driver (single-tag, RDS input, withr bundle sourcing)
├── src/
│   ├── transition_cost_rcpp.cpp      # C++ transition cost + phase feasibility functions
│   └── transition_cost_rcpp.R        # R wrapper for Rcpp-compiled functions
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
If DBH measurements were taken at differing heights, include a `HOM` column (or `hom`; detection is case-insensitive) and convert DBH to a common reference height (1.3 m) prior to running the workflow for taper-corrected growth forms (trees, shrubs, figs, unknown). For **non-taper-corrected** growth forms (palms, strangler figs, tree ferns) where taper correction is not applicable, the DP solver uses wide base pruning bounds (default 1.25× the standard growth/shrink limits) because these forms exhibit real biological DBH growth and, when the measurement point (HOM) changes between censuses, the recorded DBH can shift substantially even if the true diameter at a fixed height has not changed. When a `HOM` column is present, the solver additionally widens bounds per census pair in proportion to the worst-case HOM deviation from 1.3 m (controlled by `hom_tolerance_scale`, default 2.0 cm/yr per meter of deviation). NA HOM values are treated as 1.3 m (zero deviation). If no HOM column is present, HOM widening is disabled and only the wide base bounds apply.

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

### Life-Cycle Phase System

Each track maintains a phase variable:

- **Phase 0 (prebirth):** Track not yet occupied, can transition to Phase 1
- **Phase 1 (alive):** Track occupied, can remain alive or transition to Phase 2
- **Phase 2 (dead):** Track was occupied but now empty, cannot return to Phase 1

**Enforced constraint:** `OBS → NA → OBS` is forbidden (no resurrection)

### Choosing K (Number of Tracks)

K must accommodate:
1. Number of unique anchor stem IDs
2. Maximum observed stems in any census
3. Births needed to explain stem count increases

```text
births_needed = Σ max(0, n_obs(t+1) - n_obs(t))
K_base = max(#anchor_IDs, max_obs, n_obs(1) + births_needed)
K_final = K_base + slack_tracks
```

**Slack tracks:** Extra tracks allowing simultaneous death+birth in constant-count intervals

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
| `species` | character | Species identifier |
| `TrueStemID` | integer | Known stem ID at anchor census |

### Biological Parameter Columns (Required)

All parameters must be present in the dataset before running DP. These are typically added by calling `estimate_bio_pars()` and joining results:

#### Growth Parameters
- `Bio_Mu_Growth` ($\mu_{\text{const}}$): Intercept of mean annual growth (cm/year)
- `Bio_Gamma_Growth` ($\mu_{\gamma}$): Log-DBH slope in mean growth model (cm/year)
- `Bio_Sigma0_Growth` ($\sigma_0$): Baseline growth process SD (cm/year)
- `Bio_Sigma1_Growth` ($\sigma_1$): Growth SD slope vs DBH ((cm/year)/cm)

#### Shrinkage Penalties
- `Bio_Max_Shrink` (`max_shrink`): Hard lower bound on annual growth (cm/year, negative)
- `Bio_K_Shrink` ($k_{\text{shrink}}$): Soft shrinkage penalty weight (1/cm²)

#### Extreme Growth Penalties
- `Bio_Max_Growth` (`max_growth`): Hard upper bound on annual growth (cm/year)
- `Bio_K_Growth` ($k_{\text{growth}}$): Soft extreme growth penalty weight (1/cm²)

#### Mortality Parameters
- `Bio_H0_Mortality` ($h_0$): Baseline hazard parameter
- `Bio_Beta_Mortality` ($\beta$): Size effect on mortality hazard

#### Recruitment Parameters
- `Bio_Recruit_Meanlog`: LogNormal meanlog for recruit size
- `Bio_Recruit_Sdlog`: LogNormal sdlog for recruit size
- `Bio_Recruit_MaxDBH_unit`: Maximum plausible recruit DBH (cm)
- `Bio_Recruitment_lambda`: Recruitment Poisson rate (1/year)

### Output Columns (Added by Solver)

| Column | Description |
|--------|-------------|
| `ReconstructedStemID` | Assigned stem identity |
| `ReconstructionMethod` | One of: `"given"`, `"dp"`, `"igraph"`, `"provisional_dp"`, `"provisional_igraph"`, `"none_after_anchor"` — see notes below |
| `ConstraintViolation` | Post-hoc diagnostic flag |
| `DP_KUsed` | Number of tracks used |
| `DP_MaxStatesPerCensus` | Largest state space encountered |
| `DP_MaxStatesCensusID` | Census with largest state space |

Notes on post-anchor output semantics:
- When DP is scoped to pre-anchor censuses (because there are observations after the requested anchor), post-anchor rows are preserved and appended to the DP output. Post-anchor rows with non-NA `DBH` and a `TrueStemID` that was actually used by the DP will have `ReconstructedStemID = TrueStemID` and `ReconstructionMethod = "given"`.
- Remaining post-anchor rows without a DP assignment are labeled `ReconstructionMethod = "none_after_anchor"` and keep `ReconstructedStemID = NA`.
- If no anchored census with `TrueStemID` exists, and `allow_provisional_anchor = TRUE` (default), the DP can assign provisional anchor IDs at the last observed DBH census; those anchor rows are labeled `"provisional_dp"` and treated as anchors for the reconstruction. If DP falls back to the igraph matcher and assigns provisional IDs, those are labeled `"provisional_igraph"`.

### Posterior Uncertainty Columns (Optional)

When running `match_stems_dp_global_backward_marginals()`:

- `DP_PosteriorTop1ID`, `DP_PosteriorTop1Prob`
- `DP_PosteriorTop2ID`, `DP_PosteriorTop2Prob`
- `DP_PosteriorEntropy`
- `DP_PosteriorReconstructedProb`
- `DP_PosteriorUnlinkedProb`
- `DP_PosteriorBin` (if using `add_dp_posterior_bins()`)

Posterior path summaries (`*_paths.csv`) encode the `recon` column as `ObsRowID:ReconstructedStemID` pairs. The attachment helpers expect ObsRowID-based encodings.

### MAP vs posterior-sampled paths 🔀

**What these two outputs represent**


- **`ReconstructedStemID` (main output)** is the *MAP joint assignment* (MAP — Maximum a posteriori) decoded by the DP (a deterministic Viterbi-style backtrace of the most probable full path). This is written per-observation in the main `stem_reconstruction_*.csv` as the best joint reconstruction under the model.

- **Per-path posterior summary (`*_paths.csv`)** is an *empirical* summary of full reconstructions produced by the posterior sampler (only generated when `posterior_samples > 0`). Each row is a unique path observed among draws and `path_prob` is the normalized sampling weight for that unique path (sums to 1 across sampled unique paths).

**Why they can differ**

- Posterior sampling is finite and stochastic: the MAP joint path may have non-zero posterior mass yet still not be drawn among the finite samples. Consequently, the concatenation of per-row MAPs (`ReconstructedStemID`) may *not* appear as any `path_sig` in the `*_paths.csv` file.

**Practical recommendations**

- Increase `posterior_samples` to raise the chance the MAP path is drawn and therefore present in `*_paths.csv`.
- If you require the MAP path to be represented in per-path summaries, you can explicitly insert the MAP signature into `paths_summary` after sampling (and mark it).

**Quick check example (R)**

```r
library(data.table)
dp <- fread("dp_global/output/.../stem_reconstruction_dp_global_rcpp.csv")
sig <- paste0(dp[ReconstructionMethod == "dp", ReconstructedStemID], collapse = "-")
paths <- fread("dp_global/output/.../posteriors/..._paths.csv")
paths[path_sig == sig]  # empty => MAP path wasn't sampled
```

---

## Biological Model & Parameters

### Growth Model

**Mean growth** (size-dependent):
$\mu(D) = \alpha + \gamma \log(D)$

Where:
- $\alpha$ = `Bio_Mu_Growth` (intercept, cm/year)
- $\gamma$ = `Bio_Gamma_Growth` (slope on log-DBH, cm/year)
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
- $\lambda_{\text{recruit}}$ = `Bio_Recruitment_lambda` (rate per year)
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

---

## Core Algorithm Details

### Transition Cost Function

`transition_cost_tracks_bio()` computes the negative log-likelihood for transitioning from census t to t+1 by summing per-track costs across four mutually exclusive cases.

**Key principle:** Each track contributes independently to the total cost. The function processes all K tracks sequentially.

**Inputs:**
- `track_dbh_t`: length-K numeric vector (DBH per track at census t, NA if unoccupied)
- `track_dbh_tp1`: length-K numeric vector (DBH per track at census t+1, NA if unoccupied)
- `interval_years`: $\Delta t$ between censuses
- Biological parameters (growth, mortality, recruitment, shrinkage)
- `eps_tiebreak`: deterministic tie-break weight (default $10^{-6}$)

**Output:** Single scalar cost (sum across all tracks + tie-break term)

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

**Hard constraint:** If D_1 > `recruit_max_dbh`:
`cost += 10^6`

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

**D1. Hard shrinkage constraint:** If g < `max_shrink`: `cost += 10^6`, skip to next track

**D2. Hard extreme growth constraint:** If g > `max_growth`: `cost += 10^6`, skip to next track

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
`cost` += `epsilon_tiebreak` * sum_{k in `both_obs`} |r_0[k] - r_1[k]|

Default: $\varepsilon_{\text{tiebreak}} = 10^{-6}$

**Purpose:** Discourages unnecessary rank crossings when biological costs are tied.

### Backward DP Recursion

**Initialization:** Anchor census assignments fixed by `TrueStemID`

**Backward loop:** For each census from anchor-1 down to 1:
1. Enumerate all valid injective states at current census
2. For each state, iterate over valid next states (from census t+1)
3. Derive valid phase transitions (enforce life-cycle constraints)
4. Compute transition cost + future cost from DP table
5. Keep minimum cost and store backpointer

**Path decoding:** Follow backpointers from best start state through to anchor

**Complexity:** Managed via `max_states` parameter; falls back if exceeded

---

## Uncertainty Quantification

### Posterior Distribution

The probabilistic solver defines a Gibbs-like posterior:

**Posterior samples (optional):** If `POSTERIOR_SAMPLES` > 0, the DP can draw full-path posterior samples and write them to disk. Use `POSTERIOR_SAMPLES_FORMAT` to select `rds`, `feather` (arrow), or `csv`. By default, samples are written into a `posteriors/` subdirectory under the DP output directory (or to `POSTERIOR_SAMPLES_PATH` if that option is supplied).
$$P(\text{path}) \propto \exp\left(-\frac{\text{TotalCost}(\text{path})}{\tau}\right)$$

where $\tau$ is the temperature parameter:
- $\tau < 1$: Sharper posterior (more confident, concentrates on near-MAP paths)
- $\tau > 1$: Flatter posterior (less confident, spreads mass across alternatives)

### Marginal Computation

`match_stems_dp_global_backward_marginals()` computes exact marginals via:
1. **Backward pass:** Log-sum-exp recursion computing log partition function
2. **Forward pass:** Normalization to obtain marginal probabilities
3. **Per-observation summaries:** Probability of assignment to each anchor ID

### Posterior Summary Statistics

- **Top1/Top2 IDs and probabilities:** Most likely assignments
- **Entropy:** $H = -\sum_i p_i \log p_i$ (higher = more uncertain)
- **Reconstructed probability:** Posterior mass on chosen `ReconstructedStemID`
- **Unlinked probability:** Posterior mass on not matching any anchor ID

### Posterior Binning

`add_dp_posterior_bins()` creates categorical labels:
- **confident:** High probability on single ID
- **ambiguous:** Probability split across multiple IDs
- **unlinked-likely:** High probability of being unlinked

---

## Outputs & Diagnostics

### Primary Outputs

1. **ReconstructedStemID:** Assigned identity for each observation
2. **ReconstructionMethod:**
   - `"given"`: From input `TrueStemID`
   - `"dp"`: Assigned by DP solver
   - `"igraph"`: Assigned by fallback method

### Diagnostic Outputs

- **ConstraintViolation:** Post-hoc flag for growth violations
- **DP_KUsed:** Actual number of tracks used
- **DP_MaxStatesPerCensus:** Worst-case state space size
- **DP_MaxStatesCensusID:** Census achieving maximum states

### Visualization

**Plot reconstructed trajectories:**
```r
plot_tag_to_pdf(
  out, 
  pdf_file = "output/trajectories.pdf",
  include_reference = TRUE
)
```

**Sensitivity analysis:**
```r
source("dp_global/R/sensitivity_transition_cost_bio.R")
# Run parameter sweeps
```

**Realism calibration:**
```r
source("dp_global/R/realism_calibration.R")
# Analyze reconstruction quality
```

### Files written by the driver scripts

When you run `dp_global/scripts/main_cpp.R` or `dp_global/scripts/main_cpp_chunk.R`, outputs are written to a run-specific directory under `dp_global/output/` (the driver prints `out_dir` on startup). Common files written by the run include:

- `run_started.txt`, `run_finished.txt` — simple markers indicating job start/finish timestamps
- `run_parameters_full.txt` — full parameter dump and command-line overrides for reproducibility
- `stem_reconstruction_dp_global_rcpp.csv` — main reconstruction CSV (also written as RDS when `--WRITE_DP_RDS=TRUE`)
- `stem_reconstruction_dp_global_rcpp.rds` — binary R object of the reconstruction (written when `--WRITE_DP_RDS=TRUE`)
- `stem_reconstruction_dp_global_rcpp.pdf` — per-tag reconstruction plots (if enabled)
- `tag_<which_tag>_realism_summary_rcpp.csv`, `tag_<which_tag>_realism_by_tag_rcpp.csv`, `tag_<which_tag>_realism_tuning_suggestions_rcpp.csv` — realism report tables when `--RUN_REALISM_REPORT=TRUE`
- `simulated_all_transition_cost_sweeps_rcpp.rds`, `simulated_all_transition_cost_jumps_rcpp.csv`, `simulated_all_transition_cost_jumps_rcpp.rds` — sensitivity sweep outputs when `--SENSITIVITY_MODE` enables write

Note: enabling `--WRITE_DP_RDS=TRUE` is recommended when you plan to post-process reconstructions in R (it preserves types and attributes without re-parsing CSVs). If you need a different output location, set `--PROJECT_ROOT` and `--BATCH_TS` or supply `OUT_DIR_NAME` via an override to control the output directory naming.

---

## Parameter Estimation

### Overview

`estimate_bio_pars()` derives biological parameters from known `TrueStemID` data. Supports two modes:
- **data mode:** Estimates from empirical observations
- **fixed mode:** Uses user-specified constants

**Interval years handling:** `interval_years` may be a scalar (e.g., `5`) applied to all census pairs, or set to `NULL` to enable per-row (per-pair) interval detection from the input table. When `NULL`, the function searches for interval columns (default candidates: `"Bio_IntervalYears"`, `"IntervalYears"`, `"interval_years"`, `"census_interval_years"`, `"CensusIntervalYears"`) and uses per-pair intervals when available. Missing per-pair values are coalesced (prefer `t1`, then `t0`, then scalar); if a scalar is not provided the function may infer a representative scalar (median) from available pairs. Diagnostic information is returned in `res$interval` (see below).

### Per-row interval support

When per-census intervals vary across tags or measurement pairs you can enable per-pair interval inference by calling `estimate_bio_pars()` with `interval_years = NULL` (or by passing an explicit `interval_col_candidates` vector). Behavior summary:

- Searches the original input for interval columns (default candidates listed above) and builds a wide per-pair interval table aligned with the DBH wide table.
- For each pair, uses the interval value from the later census (`t1`) when present, otherwise falls back to `t0`, then to a supplied scalar `interval_years` if provided.
- If no scalar is provided and some pairs are missing intervals, the function infers a representative scalar (median of available per-pair intervals) and records whether rows were filled or dropped.
- Diagnostics returned in `res$interval`: `inferred_interval_years`, `per_pair_intervals` (numeric vector), `pairs_candidate_count`, `pairs_filled_with_scalar_count`, `pairs_dropped_count`.

### Growth Parameter Estimation

**Data filtering:** Keep rows with non-NA `DBH > 0` and non-NA `TrueStemID`

**Wide format conversion:**
```r
dw <- dcast(Tag + TrueStemID + species ~ CensusID, value.var = "DBH")
```

**Annual growth increments:**
For each adjacent census pair $(t_0, t_1)$:
$g = \frac{D_1 - D_0}{T_i}$ where $T_i$ is the interval (years) for that specific pair. $T_i$ may be a constant scalar (`interval_years`) or vary per pair when a per-row interval column (e.g., `Bio_IntervalYears`) is present. When per-row intervals are used, `estimate_bio_pars()` coalesces per-pair `t1` → `t0` → scalar and returns diagnostics in `res$interval`.

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
Fits `lm(sd_proc_hat ~ d0_all)` with constraints:
- $\sigma_0 \geq 10^{-4}$
- $\sigma_1 \geq 0$

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

**Data mode with measurement error:**

Typical measurement SD for raw DBH difference:
$s_{\text{typ}} = \text{median}\left(\sqrt{\text{Var}(\varepsilon|D_0) + \text{Var}(\varepsilon|D_1)}\right)$

Penalty weight:
$k_{\text{shrink}} = \frac{1}{2 s_{\text{typ}}^2}$

Clamped to $[10^{-6}, 10^6]$ for numerical stability.

**Data mode without measurement error:**

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

Solve for quantile q = 10^{-4} via:
CDF(x) = sum_{i=1}^4 w_i Phi(x; 0, SD_meas_i)

using typical diameter D_typ = median(D_0)

**Final value:**
- With measurement error: `max_shrink` = min(`max_shrink_data`, `max_shrink_meas`)
- Without: `max_shrink` = `max_shrink_data`
**Fixed mode:**
```r
max_shrink_source = "fixed"
max_shrink_fixed = -0.5  # cm/year, must be negative
```

**Application in cost function:**
If g < `max_shrink`: `cost += 10^6`

#### Soft Extreme Growth Penalty ($k_{\text{growth}}$)

**Data mode:** Same calibration as $k_{\text{shrink}}$:
$k_{\text{growth}} = \frac{1}{2 s_{\text{typ}}^2}$

Clamped to $[10^{-6}, 10^6]$.

**Data mode without measurement error:**

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

Uses same four-component mixture as shrinkage, solving for upper quantile $q = 1 - 10^{-4}$ (99.99th percentile).

**Final value:**
- With measurement error: `max_growth` = max(`max_growth_data`, `max_growth_meas`)
- Without: `max_growth` = `max_growth_data`

**Soft growth cap:**
`max_growth_soft` = min(`max_growth`, Q_{0.99}(g_all))

**Fixed mode:**
```r
max_growth_source = "fixed"
max_growth_fixed = 7.5  # cm/year
```

**Application in cost function:**
If g > `max_growth`: `cost += 10^6`

### Mortality Parameter Estimation

**Data preparation:**
For each adjacent census pair $(t_0, t_1)$:
- At-risk stems: observed at census $t_0$ (non-NA DBH)
- Died indicator: 1 if NA at $t_1$, 0 if still observed

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
Upper quantile guardrail (default 99.9%):
`recruit_max_dbh` = Q_{0.999}(D_recruited)

Fallback: 5 cm

**Optional enforcement: recruit max DBH**
`estimate_bio_pars()` supports an optional pre-fit filter for recruit sizes:

- `enforce_recruit_max` (logical, default `FALSE`): when `TRUE`, recruits with DBH strictly greater than `recruit_max_fixed` (cm) are removed from the sample prior to fitting the lognormal recruit-size distribution.
- When `enforce_recruit_max = TRUE`, you must provide a positive finite `recruit_max_fixed` value. The function will warn if any recruits were dropped.

Example (enforce a 37.5 cm recruit cap):

```r
estimate_bio_pars(x, enforce_recruit_max = TRUE, recruit_max_source = "fixed", recruit_max_fixed = 37.5)
```

**Recruitment rate:**

Poisson rate per empty slot per year:
$\lambda_{\text{recruit}} = \frac{n_{\text{recruits}}}{n_{\text{at-risk}} \cdot \sum_{\text{intervals}} \Delta t}$

where:
- $n_{\text{recruits}}$: total recruitment events across all intervals
- $n_{\text{at-risk}}$: total empty-slot census-years

Fallback: 0

### Quantile Selection

| Parameter | Default Quantile | Purpose |
|-----------|------------------|---------|
| `max_shrink_data` | 0.1% (0.001) | Lower tail for hard shrink constraint |
| `max_growth_data` | 99.9% (0.999) | Upper tail for hard growth constraint |
| `max_growth_soft` | 99% (0.99) | Upper tail for soft growth penalty |
| `recruit_max_dbh` | 99.9% (0.999) | Maximum plausible recruit size |

All quantiles are configurable via `estimate_bio_pars()` arguments.

---

## Fallback Mechanisms

### When Fallback Occurs

In addition to the structural and numerical conditions below, the DP routine
can be instructed to **bypass the solver based on a data column named
`growth_form`**.  If the caller supplies a non-empty character vector via the
`fallback_growth_forms` argument (propagated from the CLI variable
`DP_FALLBACK_GROWTH_FORMS`), any tag whose pre-anchor subset contains at least
one row whose `growth_form` value matches the vector will immediately fall back
with reason `growth_form_forced`.  The CLI flag accepts comma‑ or
semicolon‑separated lists, which are split automatically by the DP loader.


The DP solver automatically falls back to `match_stems_optimal_backward()` when:

1. Anchor census missing or has no observed stems
2. Any census has too many injective states (P(K, n_obs) > `max_states`)
3. K insufficient ($K < \max$ observed stems)
4. DP recursion yields no feasible keys

#### Anchor fallback behavior
If the user requests an `anchor_start` that exists in the dataset but **all rows at that census have NA for both `DBH` and `TrueStemID`**, the algorithm will search backwards and select the most recent earlier census that has at least one row with a non-NA `DBH` and a non-NA `TrueStemID` and use that census as the anchor instead of immediately falling back to the igraph matcher. If no such earlier census exists, the algorithm falls back to `match_stems_optimal_backward()`.

#### Anchor extension (forward search)
When the nominal anchor census has 0 living stems (all DBH values are `NA`), the algorithm extends the anchor search **forward** to post-anchor censuses. It selects the first census after the anchor that has at least one living stem with a non-NA `TrueStemID`. This prevents tags where the anchor census happens to have only dead stems from unnecessarily falling back to the igraph matcher when a later census has reliable identity information.

Note: In addition, when the requested anchor census contains DBH observations but lacks `TrueStemID` values, setting the `ALLOW_PROVISIONAL_DP_ANCHOR` flag (default: `TRUE` in the chunk runner and the DP function) allows the DP to assign provisional `TrueStemID`/`ReconstructedStemID` values at the last-observed DBH census (marked with `ReconstructionMethod = "provisional_dp"`) and proceed with the DP instead of falling back to the igraph matcher.

### Fallback Method: `match_stems_optimal_backward()`

**Algorithm:**
1. Works backward one census at a time
2. Uses `igraph::max_bipartite_match()` with deterministic edge weights
3. Edges allowed only when growth ∈ [min_growth, max_growth]
4. If more stems at t+1 than t, smallest future stems treated as recruits
5. Unmatched stems receive new IDs

**Marking:** Rows assigned by fallback have `ReconstructionMethod="igraph"`

### DP fallback reason codes (`DP_FallbackReason`)

To aid diagnostics the DP sets a `DP_FallbackReason` string on returned rows when the solver falls back to the igraph matcher. Common codes and suggested remedies:

| Code | Meaning | Suggested action |
|------|---------|------------------|
| `no_obs_up_to_anchor` | No DBH observations exist at or prior to the requested anchor. | Check `anchor_start` or ensure DBH observations exist before the anchor. |
| `anchor_missing_truestem` | Anchor census missing `TrueStemID` values and no earlier anchor found. | Supply `TrueStemID` or enable `allow_provisional_anchor=TRUE`. |
| `anchor_missing_obs` | Anchor census has no DBH observations. | Ensure DBH at anchor or select a different `anchor_start`. |
| `anchor_missing_truestem_prov_disabled` | Anchor has DBH but `TrueStemID` missing and provisional anchors are disabled. | Set `allow_provisional_anchor=TRUE` or provide `TrueStemID`. |
| `provisional_igraph_anchor_assigned` | The igraph matcher assigned provisional anchor IDs at the last observed census. | Verify provisional assignments or provide true anchors. |
| `anchor_ids_missing` | No anchor IDs were found after attempting to locate an anchor census. | Inspect data or supply anchor IDs. |
| `K_too_small` | Chosen `K` is smaller than the maximum per-census observed stems. | Increase `max_tracks` or `slack_tracks`. |
| `enum_exceeded` | State enumeration exceeded `max_states` for a census (too many injective assignments). | Increase `max_states` or reduce `K`. |
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
| `no_dbh` | (matcher) No DBH observations present; nothing to assign. | Ensure DBH data are present for the Tag. |
| `igraph_assigned` | (matcher) Rows explicitly assigned by the igraph fallback. | None — indicates fallback assignment occurred. |

> Note: When DP is scoped to pre-anchor censuses and post-anchor rows are appended, any `DP_FallbackReason` present on the DP output will be propagated to matching post-anchor rows. If multiple reasons apply, they will be concatenated with `;` in the propagated column.

**How to check for fallback reasons (R example):**

```r
res <- match_stems_dp_global_backward_marginals_batch(dt, anchor_start = 5, ...)
unique(na.omit(res$DP_FallbackReason))
```

**Important:** `min_growth` and `max_growth` constraints affect **only** the fallback method and post-hoc diagnostics, not the DP objective.

---

## ForestGEO Code Handling

The DP solver recognizes ForestGEO stem status codes that encode field-recorded events. These codes influence state enumeration and transition feasibility.

### Resprout Codes (R, OR)

**Recognized codes:** `R` (resprout from main stem) and `OR` (other breakage resprout).

**R-recruit constraint:** When a living observation at census $t+1$ carries an R or OR code, the DP forces that stem to appear as a **new recruit** (empty track at census $t$). This prevents the algorithm from chaining a resprout observation to a pre-existing identity track, which would violate the biological meaning of the code: the stem broke and re-grew, so it is a new physical entity.

The constraint is implemented as a hard filter in the backward transition loop: any candidate assignment that places an R/OR-coded living stem on a track that was already occupied at the prior census is removed from the feasible set.

### Main-Stem Code (M)

**Recognized code:** `M` (main stem).

**M-pin constraint:** At branching events (censuses where $n_{obs} > 1$ at both the current and next census), any observation carrying an M code is **pinned to track 1**. This ensures the main stem retains identity priority and prevents the DP from swapping it with secondary stems during optimization.

The constraint is applied only when:
- The prior census has at least one observed stem ($n_{obs,t} > 0$)
- Both the current and next census have multiple stems

When triggered, the state enumeration is filtered to retain only assignments where the M-coded observation occupies track 1.

### Missing-from-Field (MF) Detection

The DP detects implicit "missing from field" events: censuses where all stems have `NA` DBH, sandwiched between censuses with observed stems. These are treated as temporary absences rather than mortality.

**R-code exclusion:** Censuses where all observations are `NA` but carry resprout codes (R or OR) are excluded from MF detection. A resprout code on an unobserved stem indicates a known break event, not a temporary field absence.

### Regex Implementation

All ForestGEO code matching uses `grepl()` with `perl = TRUE` to ensure `\b` word-boundary anchors work correctly. The resprout regex pattern is:

```r
"\\b(R|OR)\\b"
```

This matches whole-word `R` or `OR` codes in fields that may contain multiple semicolon-separated codes (e.g., `"R;B"` matches, `"NORMAL"` does not).

## Pruning & Conservative Guards

### Motivation — factorial growth of the state space
The DP enumerates injective assignment states per census. For a census with `n_obs` observed stems and `K` tracks the number of assignment states is the falling-permutation

```
P(K, n_obs) = K × (K-1) × ... × (K - n_obs + 1)
```

The total number of full-paths (and candidate transitions between adjacent census assignment states) grows multiplicatively across censuses, which quickly becomes intractable (factorial-like explosion). To keep the DP practical, we apply *conservative, cheap* pre-filters that remove biologically impossible or extremely implausible candidate transitions before evaluating the full (expensive) transition cost.

These filters are intentionally conservative: they are meant to remove clearly impossible candidates (hard physical/biological limits) while preserving borderline cases for the full cost evaluation.

### Where pruning runs
- Pruning happens in the backward recursion immediately *before* calling `transition_cost_tracks_bio_batch_rcpp` for a candidate transition. This avoids the costly likelihood evaluation for obviously impossible transitions.
- Pruning is only applied when `prune_hard = TRUE` (default). If `prune_hard = FALSE`, no pre-filtering is performed and all candidate transitions are passed to the cost routine (slower but conservative).

### Parameters & exact semantics
We separate pruning thresholds from the biological (`Bio_*`) parameters so you can control pruning behaviour without changing the biological model used in the transition-cost calculations.

- `prune_min_growth` (numeric | NULL): explicit lower bound (cm/year) used for pruning. If `NULL` the value `min_growth` passed to the DP is used.
- `prune_max_growth` (numeric | NULL): explicit upper bound (cm/year) used for pruning. If `NULL` the value `max_growth` passed to the DP is used.
- `prune_use_bio_bounds` (logical, default TRUE): whether to intersect the user-specified prune bounds with the biological hard bounds found in the data (`Bio_Max_Shrink`, `Bio_Max_Growth`). If TRUE (default):
  - eff_min_grow = max(user_min, Bio_Max_Shrink)
  - eff_max_grow = min(user_max, Bio_Max_Growth)
  where `user_min` = `prune_min_growth` if provided else `min_growth`, and similarly for `user_max`.
  If `prune_use_bio_bounds = FALSE` the `eff_*` values equal `user_*` directly.
- `prune_recruit_max_dbh` (numeric | NULL): override for the recruit-size cut used during pruning. If `NULL` the biological `Bio_Recruit_MaxDBH_unit` is used. If `prune_use_bio_recruit = TRUE` (default) and both are finite, we use the more conservative value `min(prune_recruit_max_dbh, Bio_Recruit_MaxDBH_unit)`; otherwise the explicit override is used.
- `prune_use_bio_recruit` (logical, default TRUE): controls whether the biological recruit bound should be intersected with `prune_recruit_max_dbh`. If TRUE, we select the minimum as min(prune_recruit_max_dbh, maximum recruitment dbh). Note: maximum recruitment dbh could be either the 0.999 quantile of recruitment DBH or a fixed value defined earlier when estimating the parameters.

#### Non-taper-corrected growth form override

Growth forms that are **not taper-corrected** (default: `"palm"`, `"strangler_fig"`, `"tree_fern"`) need wider pruning bounds than standard trees for two reasons:

1. **Real biological growth.** Palms, strangler figs, and tree ferns grow in DBH — palms can add 1–3 cm/yr, and strangler figs can change diameter substantially as they encircle or replace their host.
2. **Apparent DBH variation from HOM changes.** When the measurement height (HOM) shifts between censuses, the recorded DBH can change dramatically even if the true diameter at a fixed height remained constant. Because these growth forms lack the tapered trunk geometry that allows taper correction, any HOM shift produces an uncompensated apparent DBH change.

For these reasons, the non-taper-corrected override **replaces** the general effective pruning bounds with values that are wider than (or at least as wide as) the standard growth/shrink limits, rather than tighter. The default bounds are 1.25× the standard limits (`MAX_SHRINK_FIXED` and `MAX_GROWTH_FIXED`). An optional HOM-proportional widening layer then extends the bounds further on a per-census-pair basis when HOM data are available.

**Parameters:**

- `non_taper_corrected_growth_forms` (character vector, default `c("palm", "strangler_fig", "tree_fern")`): growth forms whose DBH measurements are not taper-corrected. When a tag's `growth_form` matches any entry in this list, the non-taper override is activated. Accepts comma- or semicolon-separated strings.
- `non_taper_corrected_prune_min_growth` (numeric, default `-0.625`): lower prune bound (cm/year) that replaces the general effective minimum for matching growth forms. Default is `1.25 × MAX_SHRINK_FIXED`.
- `non_taper_corrected_prune_max_growth` (numeric, default `6.25`): upper prune bound (cm/year) that replaces the general effective maximum for matching growth forms. Default is `1.25 × MAX_GROWTH_FIXED`.
- `hom_tolerance_scale` (numeric, default `2.0`): additional annual DBH tolerance (cm/yr) per meter of HOM deviation from 1.3 m. When a `hom` (or `HOM`) column is present in the data and the tag is non-taper-corrected, the prune bounds are widened for each census pair by:

  ```
  hom_tol = hom_tolerance_scale × max(|HOM − 1.3|, na.rm=TRUE) / interval_years
  eff_min_grow_pair = non_taper_corrected_prune_min_growth − hom_tol
  eff_max_grow_pair = non_taper_corrected_prune_max_growth + hom_tol
  ```

  where `max(|HOM − 1.3|)` is the worst-case deviation across all stems at the two censuses being compared. NA HOM values are treated as 1.3 m (zero deviation contribution). Set `hom_tolerance_scale = 0` to disable HOM widening.

**How the override layers interact:**

1. **User layer**: general prune bounds are set from `prune_min/max_growth` (or `min/max_growth` if NULL).
2. **Bio layer**: if `prune_use_bio_bounds = TRUE`, general bounds are tightened by intersecting with biological hard limits (`Bio_Max_Shrink`, `Bio_Max_Growth`).
3. **Non-taper override**: if the tag's `growth_form` matches `non_taper_corrected_growth_forms`, the effective bounds from steps 1–2 are **replaced** with `non_taper_corrected_prune_min/max_growth`. This ensures non-taper forms always get wide bounds regardless of how tight the bio-constrained general bounds may be.
4. **HOM widening**: if a `HOM` column is present and `hom_tolerance_scale > 0`, the bounds from step 3 are **widened symmetrically** for each census pair based on the worst-case HOM deviation.

**Note on the default driver configuration:** The driver scripts (`main_cpp_chunk.R`, `main_cpp.R`) pass `prune_min_growth = 1.25 × MAX_SHRINK_FIXED` and `prune_max_growth = 1.25 × MAX_GROWTH_FIXED` with `prune_use_bio_bounds = FALSE`. This means the general bounds already match the non-taper defaults, so in the default configuration the non-taper override at step 3 is a no-op — all trees (taper-corrected or not) receive the same base prune window of `[-0.625, 6.25]`. The effective differentiation for non-taper forms comes from step 4 (HOM widening), and the non-taper override provides an independent lever when users tighten the general bounds (e.g. via `prune_use_bio_bounds = TRUE` or narrower `prune_min/max_growth`).

Notes:
- These pruning parameters affect only the *pre-filtering* step (they do not change the transition cost function or post-hoc diagnostics beyond recording the effective prune values).
- If interval length between censuses is invalid or not finite for a pair, growth-based pruning is skipped for that pair and only recruit-size pruning (if applicable) is used.
- You can always define extra margins for these parameters. For example:

```r
prune_min_growth = MAX_SHRINK_FIXED * 2.5 # very wide fixed bounds
prune_max_growth = MAX_GROWTH_FIXED * 1.5 # very wide fixed bounds
prune_use_bio_bounds = FALSE # use fixed prune bounds instead of biological ones
prune_recruit_max_dbh = RECRUIT_MAX_FIXED * 1.2 # very high recruit max dbh
prune_use_bio_recruit = FALSE# use fixed prune bounds instead of biological ones
```

### Diagnostics & reproducibility
- The DP exposes `attr(out, "DP_PruneInfo")` with:
  - `total_examined`, `total_pruned` (counts)
  - `per_census` (count per census pair)
  - `eff_min_growth`, `eff_max_growth`, `eff_recruit_max` (the effective thresholds used)
- The effective prune bounds are also `vcat`-logged at the start of the backward pass for transparency.

### Practical guidance and examples
- Default behaviour (safe): leave `prune_min_growth = NULL` and `prune_max_growth = NULL`. The code will use `user_min = min_growth`, `user_max = max_growth` and intersect them with `Bio_*` values. This preserves biological hard limits while letting you set study-level `min_growth`/`max_growth` to narrow behaviour across the run.
- To apply *wider* pruning (less aggressive rejection) than the biological bounds allow (for example, when the DP cost model is trusted to handle extreme cases), set `prune_use_bio_bounds = FALSE` and set `prune_min_growth`/`prune_max_growth` to the desired wide range (e.g., `-10`..`25`).
- To apply *stricter* recruit-size pruning, provide a small `prune_recruit_max_dbh` (and set `prune_use_bio_recruit = FALSE` to use it without intersecting the biological bound).
- Beware of overly tight pruning: if pruning removes all feasible transitions at a census the DP will fallback to the `igraph` matcher (or produce no DP states). If you see frequent fallbacks for reasonable data, relax the prune bounds or set `prune_hard = FALSE` for a diagnostic run.

### Why this approach?
- The combinatorial nature of injective assignments makes per-census state counts grow factorially with `n_obs`. Pruning cheaply removes impossibilities and reduces the number of pairwise assignment evaluations from `O(n_states_cc × n_states_{cc+1})` to a manageable number while retaining feasible candidates for the full probabilistic scoring.
- Because pruning is conservative and logged, it is auditable: you can use `DP_PruneInfo` to quantify how many candidate transitions were removed and where.

---

## Workflows & Usage Patterns

### Single-Tag Debug

```r
# In dp_global/scripts/main_cpp.R (or via CLI):
RUN_ALL_TAGS <- FALSE
WHICH_TAG <- 123456L

WRITE_DP_PDF <- TRUE
```

### Single-Tag with Posterior Uncertainty

```r
DP_MODE <- "marginals+bins"  # default; also adds posterior bins
# temperature is passed directly to match_stems_dp_global_backward_marginals_batch()
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
# Run parameter sweeps on costs
```

### Example: enforce growth/recruit bounds
Units: growth bounds are **cm/year**; recruit caps are **cm**.

```r
# Enforce growth bounds (0.05 - 0.75 cm/yr) and cap recruits at 37.5 cm
bio <- estimate_bio_pars(
  x,
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
# Analyze reconstruction quality metrics
```

### Parameter Estimation Workflow

```r
# 1. Estimate parameters per species
bio_pars_list <- list()
for (sp in unique(xraw$species)) {
  sp_data <- xraw[species == sp]
  # To use per-row/per-pair interval years from the dataset, set `interval_years = NULL` and
  # ensure an interval column exists (e.g., `Bio_IntervalYears`). Example: `interval_col_candidates = "Bio_IntervalYears"`.
  bio_pars_list[[sp]] <- estimate_bio_pars(
    sp_data,
    interval_years = NULL, # prefer per-row/per-pair interval columns (e.g., Bio_IntervalYears); set scalar like 5 to force a constant interval
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
# You may also supply a per-row `Bio_Recruit_MaxDBH_unit` column in the input to override the default.
# Additionally, you can force a fixed recruit max via runtime flags:
#   RECRUIT_MAX_SOURCE = "fixed"  # or "data"
#   RECRUIT_MAX_FIXED  = 6.0        # cm when RECRUIT_MAX_SOURCE="fixed"
# By default RECRUIT_MAX_SOURCE="data" (estimate from observed recruits).
xraw[, Bio_Recruitment_lambda := bio_pars_list[[species]]$recruitment$lambda, by = species]

# 3. Run DP solver (marginals + posterior bins)
out <- xraw[, match_stems_dp_global_backward_marginals_batch(
  .SD,
  min_growth = -2,           # Used by fallback and diagnostics only
  max_growth = 10,           # Used by fallback and diagnostics only
  interval_years = NULL,     # set to NULL to prefer per-row interval column (e.g., Bio_IntervalYears) or pass a scalar like 5
  anchor_start = 7,
  max_tracks = 30,
  max_states = 50000,
  slack_tracks = 1,
  dp_mode = "marginals+bins",
  temperature = 0.5,
  posterior_top_k = 2,
  use_measurement_error = USE_MEASUREMENT_ERROR,
  verbose = TRUE
), by = .(Tag, species)]

# 4. Add posterior bins (applied per-group or to full out)
out <- add_dp_posterior_bins(
  out,
  confident_prob = 0.95,
  unlinked_prob = 0.50,
  use_reconstructed_prob = TRUE
)
```

**Important workflow notes:**

1. **min_growth/max_growth in DP calls:** These parameters are **NOT used** by the DP objective function. They only affect:
   - The fallback `match_stems_optimal_backward()` (stepwise igraph matcher)
   - Post-hoc `ConstraintViolation` diagnostics

2. **Actual constraints used by DP:** Come from Bio_* columns:
   - Hard shrinkage: `Bio_Max_Shrink`
   - Hard growth: `Bio_Max_Growth`
   - Soft shrinkage: `Bio_K_Shrink`
   - Soft growth: `Bio_K_Growth` and `Bio_Max_Growth_Soft`

3. **Measurement error:** Must be set consistently between `estimate_bio_pars()` and DP solver calls.

4. **Parameter sources:** Each parameter's source (data/fixed) and estimated value are stored in the returned list structure for provenance tracking.

5. **Per-row intervals:** You can enable per-pair intervals for parameter estimation by providing a per-row interval column (e.g., `Bio_IntervalYears`) and calling `estimate_bio_pars()` with `interval_years = NULL` (or setting `interval_col_candidates` explicitly). The function will prefer per-pair `t1` values, fall back to `t0`, then to a scalar if provided, and will report diagnostic counters in `res$interval`.

---

## Implementation Reference

### Function Map

| Task | Primary Function | Location |
|------|------------------|----------|
| State enumeration | `enumerate_states_injective()` | `dp_global/R/dp_global_states.R` |
| Track DBH vectors | `state_to_track_dbh()` | `dp_global/R/dp_global_states.R` |
| Transition costs (C++) | `transition_cost_tracks_bio_batch_rcpp_cpp()` | `dp_global/src/transition_cost_rcpp.cpp` |
| Cost breakdown (debug) | `transition_cost_tracks_bio_components()` | `dp_global/R/dp_global_bio.R` |
| Posterior marginals (production) | `match_stems_dp_global_backward_marginals_batch()` | `dp_global/R/dp_global_dp.R` |
| Posterior binning | `add_dp_posterior_bins()` | `dp_global/R/dp_global_diag.R` |
| Fallback matcher | `match_stems_optimal_backward()` | `dp_global/R/dp_global_matchers.R` |
| Parameter estimation | `estimate_bio_pars()` | `dp_global/R/dp_global_bio.R` |
| DP complexity estimate | `estimate_dp_complexity()` | `dp_global/R/complexity/estimate_dp_complexity_function.R` |
| Plotting | `plot_tag_to_pdf()` | `dp_global/R/dp_global_diag.R` |
| Driver (interactive / single-tag) | `run_dp_one_group()` | `dp_global/scripts/main_cpp.R` |
| Driver (chunked / large runs) | `run_main_chunked()` | `dp_global/scripts/main_cpp_chunk.R` |
| Driver (BCI debug, single-tag) | (sources `main_cpp.R`) | `dp_global/scripts/main_cpp_bci.R` |

### match_stems_dp_global_backward_marginals_batch — Function reference (implementation details)

This implementation is the production-grade, batch-capable marginal DP solver. Below are the key computations, helper functions it uses, and the semantics of important parameters (including the new pruning controls).

Key computations and helpers:

- Anchor selection
  - If requested `anchor_start` has no DBH/TrueStemID, the function searches backward for the most recent census with at least one row having non-NA DBH and non-NA TrueStemID and uses that as the anchor. If the anchor has 0 living stems, the function also searches forward for the first post-anchor census with living stems and a non-NA TrueStemID. If no valid anchor is found in either direction, it falls back to `match_stems_optimal_backward()`.

- State enumeration
  - `enumerate_states_injective(K, n_obs, max_states)` enumerates injective assignments (permutation-based states). If enumeration exceeds `max_states` it falls back to igraph.
  - `count_injective_states(K, n_obs)` computes theoretical counts used for diagnostics.

- Track DBH vectors
  - `state_to_track_dbh(assign_vec, obs_dbh, K)` constructs length-K DBH vectors for each state used in cost evaluation.

- Phase constraints
  - `derive_phase_prev_batch_rcpp(...)` (C++, `dp_global/src/transition_cost_rcpp.cpp`) checks phase-transition feasibility for all (i, j) assignment pairs in batch; returns `from_i`, `to_j`, and `phase_t` matrix for feasible pairs.

- Interval computation
  - Uses per-census mean `ExactDate` to compute `interval_val` (years) between census pairs. If `interval_val` is NA or non-finite, growth-based pruning is skipped for that pair.

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

    # Non-taper-corrected override: replace effective bounds with wide
    # base values so palms/strangler figs/tree ferns are not spuriously
    # pruned. This step matters when prune_use_bio_bounds=TRUE (which
    # would tighten the general bounds); in the default driver config
    # (prune_use_bio_bounds=FALSE, general bounds already 1.25×) it is
    # a no-op because the values are the same.
    if growth_form in non_taper_corrected_growth_forms:
      eff_min_grow = non_taper_corrected_prune_min_growth  (override)
      eff_max_grow = non_taper_corrected_prune_max_growth  (override)

    # HOM widening: per census pair, widen bounds symmetrically by the
    # worst-case HOM deviation across both censuses. This is the main
    # source of differentiation for non-taper forms in the default config.
    Per census pair, if HOM column present and hom_tolerance_scale > 0:
      hom_tol = hom_tolerance_scale × max(|HOM − 1.3|) / interval_years
      eff_min_grow_pair = eff_min_grow − hom_tol
      eff_max_grow_pair = eff_max_grow + hom_tol

    eff_recruit_max chosen from `prune_recruit_max_dbh` and `Bio_Recruit_MaxDBH_unit` depending on flags.

  - For each candidate (assignment-state pair) and each track with DBH at both times compute

    g = (D_{t+1} - D_t) / Δt

    and prune if g ∉ [eff_min_grow, eff_max_grow]; for recruits (NA→DBH) prune if DBH > eff_recruit_max.

  - Pruning updates `prune_stats` and avoids calling the expensive `transition_cost_tracks_bio_batch_rcpp()` for pruned candidates.

- Transition cost computation
  - `transition_cost_tracks_bio_batch_rcpp(track_dbh_t, track_dbh_tp1, interval_years, ...)` computes per-candidate cost (sum across tracks) including growth likelihood, mortality, recruitment, and hard-penalties (1e6) for impossible transitions.

- Backward recursion & Viterbi
  - Uses `log_sum_exp` to accumulate marginal weights and a Viterbi update to compute MAP path (per-step `vit_cost`, `vit_ptr` arrays).

- Forward pass & posterior marginals
  - Builds per-observation posterior distributions (`DP_PosteriorTopK*` columns) via normalized state weights.

- Posterior sampling
  - Optional sample drawing from the DP graph (`posterior_samples`), with backward sampling and per-sample `logp` weights converted to `sample_prob`.

- Diagnostics & fallbacks
  - `attr(out, "DP_PruneInfo")` contains pruning diagnostics (counts and effective thresholds). The function falls back to `match_stems_optimal_backward()` on anchor failures, enumeration exhaustion, K insufficiency or if DP produces no feasible states after pruning. These fallback return values always include `DP_PruneInfo`.

Notes:
- `min_growth`/`max_growth` still primarily control the **fallback matcher** and post-hoc `ConstraintViolation` checks. They can be defined as `Bio_Max_Shrink`/`Bio_Max_Growth`. 

---

### Key Implementation Details

**NA interpretation:** `NA` means "not observed," not "alive but missed"

**Large penalties:** `1e6` encodes hard constraints (impossible transitions)

**Measurement error:** Controlled by `use_measurement_error` flag, affects:
- Likelihood computation (4-component mixture)
- Parameter estimation (variance subtraction, quantile selection)

**Phase variables:** Extended DP key format: `"<assignment>|<phase0,phase1,...,phaseK>"`

**Anchor initialization:** Non-anchor tracks start in phase 2 (dead), allowing reconstruction of earlier deaths

### Debugging Tools

**Transition cost breakdown:**
```r
components <- transition_cost_tracks_bio_components(
  track_dbh_t, track_dbh_tp1, interval_years, ...
)
# Returns per-track: growth_ll, shrink_soft, growth_soft, mortality, recruit
```

**Posterior analysis:**
```r
# After running marginal solver:
high_uncertainty <- out[DP_PosteriorEntropy > 1.5]
ambiguous <- out[DP_PosteriorBin == "ambiguous"]
```

### Complexity Estimation

The DP algorithm's computational complexity depends on the number of observations per census and the resulting state space size. Use the complexity estimator to identify which tags will take the longest to process before running the full workflow.

**Load the complexity estimation functions:**
```r
source("R/estimate_dp_complexity_function.R")
```

**Estimate complexity for all tags:**
```r
complexity <- estimate_dp_complexity("../data_simulation/data/simulated_data_1.csv")
print(complexity)
```

**Get detailed analysis for a specific tag:**
```r
details <- get_tag_complexity_details("../data_simulation/data/simulated_data_1.csv", tag = 11)
print(details$observations_per_census)  # Observations per census
print(details$states_per_census)        # States per census
cat("Total transition computations:", format(details$transition_computations, big.mark = ","), "\n")
```

**Output columns:**
- `Tag`: The tag identifier
- `Species`: Species name  
- `MaxObs`: Maximum number of observations in any single census
- `K`: Number of tracks used by the DP algorithm
- `TotalStates`: Total number of states across all censuses
- `MaxStatesPerCensus`: Maximum states in any single census
- `TransitionComputations`: Estimated number of transition cost calculations (proxy for runtime)

**How it works:**
1. **Number of tracks (K)**: Determined by the maximum observations per census plus slack tracks
2. **States per census**: P(K, n_obs) = K! / (K - n_obs)! where n_obs is observations in that census
3. **Transitions**: For each pair of consecutive censuses, all pairs of states must be evaluated

Tags are sorted by `TransitionComputations` in descending order, so the slowest tags appear first.

**Example output for your dataset:**
Tag 11 requires ~89 million transition computations, while most other tags require <1 million.

---

## Notes & Common Issues

### Interpretation Notes

1. **DBH → NA is mortality**, not missing measurement
2. **K selection** is critical: too small fails, too large is inefficient
3. **Slack tracks** enable realistic death+birth dynamics in constant-count transitions
4. **Temperature** in posterior controls confidence (not a physical parameter)
5. **NA means "not observed"**, not "alive but missed"

### Parameter Roles Summary Table

| Parameter Type | Parameter Name | Used By | Purpose | Typical Value |
|----------------|---------------|---------|---------|---------------|
| **Hard Constraint (DP)** | `Bio_Max_Shrink` | DP cost function | Hard lower bound on $g$ | -0.5 to -1.0 cm/year |
| **Hard Constraint (DP)** | `Bio_Max_Growth` | DP cost function | Hard upper bound on $g$ | 5 to 10 cm/year |
| **Soft Penalty (DP)** | `Bio_K_Shrink` | DP cost function | Quadratic shrinkage weight | 0 to 50 (1/cm²) |
| **Soft Penalty (DP)** | `Bio_K_Growth` | DP cost function | Quadratic growth excess weight | 0 to 50 (1/cm²) |
| **Soft Threshold (DP)** | `Bio_Max_Growth_Soft` | DP cost function | Growth excess threshold | 95% of `Bio_Max_Growth` |
| **Fallback Only** | `min_growth` | Fallback matcher | Edge constraint for stepwise | -2 cm/year |
| **Fallback Only** | `max_growth` | Fallback matcher | Edge constraint for stepwise | 10 cm/year |
| **Diagnostic Only** | `min_growth` | `ConstraintViolation` | Post-hoc violation flag | Same as fallback |
| **Diagnostic Only** | `max_growth` | `ConstraintViolation` | Post-hoc violation flag | Same as fallback |

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
   - `min_growth`: Used in fallback and diagnostics only
   - They can have different values (usually min_growth is more permissive)

6. **"Why do soft penalties exist if hard constraints are enforced?"**
   - Hard: Binary (allowed/forbidden), creates discrete cost landscape
   - Soft: Continuous discouragement, helps DP choose among many feasible paths
   - Together: Hard prevents absurdity, soft guides toward realism

### Performance Considerations

- State space grows factorially: $O(K^{n_{obs}})$
- Set `DP_MAX_STATES` based on memory availability (default 40,000)
- Parallel processing via `MANUAL_CORES_VALUE` (controlled by `MANUAL_CORES=TRUE`) for multiple tags
- Consider fallback threshold adjustment for large datasets
- Batch cost computation (`match_stems_dp_global_backward_marginals_batch`) is faster for posterior inference

### Numerical Stability Details

| Operation | Stability Mechanism | Value |
|-----------|-------------------|-------|
| Growth SD | Floor clamp | $\sigma(D) \geq 10^{-6}$ |
| Measurement SD1 | Floor clamp | $\text{SD1}(D) \geq 10^{-6}$ |
| Mortality probability | Range clamp | $[10^{-12}, 1-10^{-12}]$ |
| Recruitment probability | Range clamp | $[10^{-12}, 1-10^{-12}]$ |
| Log-sum-exp mixture | Subtract max before exp | Prevents overflow/underflow |
| K estimation (data mode) | Range clamp | $k \in [10^{-6}, 10^6]$ |
| Rank tie-breaking | Deterministic first index | `ties.method="first"` |

---

## References

**Measurement error model:**
Chave, J., Condit, R., Aguilar, S., Hernandez, A., Lao, S., & Perez, R. (2004). Error propagation and scaling for tropical forest biomass estimates. *Philosophical Transactions of the Royal Society of London. Series B: Biological Sciences*, 359(1443), 409-420.

**Additional documentation:**
- Sensitivity analysis: `dp_global/R/sensitivity_transition_cost_bio.R`
- Realism calibration: `dp_global/R/realism_calibration.R`
- API documentation: See function headers in the relevant R files (`dp_global/R/dp_global_dp.R`, `dp_global/R/dp_global_bio.R`, etc.)

---

## Building This Documentation

Use `pandoc` (if available) to render `README.md` to HTML, or use your preferred tooling.