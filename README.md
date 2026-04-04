# STEM_IDENTIFICATION_TEST

## Overview

Biologically informed dynamic-programming (DP) solver that reconstructs stem identities across forest censuses backward in time from a known anchor census. Given multi-stem tree measurements and a late-census anchor with trusted `TrueStemID`, the algorithm assigns each earlier observation to a latent identity track by minimising negative log-likelihood costs that encode growth, mortality, and recruitment biology. Uncertainty is quantified via forward-backward marginals and optional posterior sampling.

When the exact DP state space is too large (combinatorial explosion from many stems per tree), the solver automatically falls back to a **probabilistic greedy matching** module that uses the same biological cost model with Gumbel-noise stochastic sampling to produce approximate reconstructions with per-observation posterior probabilities.

## Directory Layout

```
├── dp_global/                 # Core algorithm, drivers, and C++ acceleration
│   ├── R/                     # R modules (sourced by dp_global_main.R)
│   │   ├── dp_global_main.R           # Module loader (sources all R modules in order)
│   │   ├── dp_global_bio.R            # Biological parameter estimation
│   │   ├── dp_global_states.R         # State enumeration & track-DBH helpers
│   │   ├── dp_global_matchers.R       # Fallback igraph matcher
│   │   ├── dp_probabilistic_matching.R # Probabilistic greedy matching fallback
│   │   ├── dp_global_dp.R            # Core DP solver (backward/forward pass, marginals)
│   │   ├── dp_global_utils.R         # Shared utilities
│   │   ├── dp_global_diag.R          # Diagnostics & PDF plotting
│   │   ├── naming_helpers.R          # Output directory naming
│   │   ├── complexity/               # DP complexity estimator
│   │   └── dpglobal_bundle/          # Portable deployment bundle builder
│   ├── scripts/               # CLI driver scripts (main_cpp*.R)
│   └── src/                   # C++ transition cost (Rcpp)
├── data_simulation/           # Simulated forest-census data generator
│   └── data/                  # Generated test datasets (CSV)
├── bci_data/                  # BCI census data (not tracked by git)
└── Makefile                   # Convenience targets (smoke test)
```

**Not tracked by git** (see `.gitignore`): `dp_global/output/`, `dp_global/ForestGEO_codes/`, `dp_global/examples/`, `bci_data/`, `posteriors/`, `tests/`, `*.rds`, `*.pdf`, `*.log`.

## Prerequisites

R ≥ 4.0 with packages: `data.table`, `igraph`, `Rcpp`, `here`.
Optional: `ggplot2`, `cowplot` (plotting), `arrow` (feather output), `withr` (bundle sourcing).

## Quickstart

```bash
# Verify modules load
make smoke

# Single-tag run on simulated data
Rscript dp_global/scripts/main_cpp.R --WHICH_TAG=20

# Full run (all tags, chunked output with resume support)
Rscript dp_global/scripts/main_cpp_chunk.R --RUN_ALL_TAGS=TRUE --DP_CHUNK_SIZE=7

# BCI data (single tag)
Rscript dp_global/scripts/main_cpp_bci.R --WHICH_TAG=123375
```

Key CLI parameters for controlling solver behavior:

| Flag | Default | Purpose |
|------|---------|---------|
| `--DP_MAX_STATES` | `40000` | Max injective states per census before probabilistic fallback |
| `--PROB_N_SAMPLES` | `200` | Number of Gumbel-noise samples for probabilistic matching |
| `--POSTERIOR_SAMPLES` | `200` | Number of posterior path samples (0 to disable) |

### Understanding `DP_MAX_STATES`

`DP_MAX_STATES` controls the maximum number of assignment states the DP solver will enumerate at any single census. The number of states at a census is the number of ways to assign $n$ observed stems to $K$ identity tracks:

$$P(K, n) = K \times (K-1) \times \cdots \times (K - n + 1)$$

This grows factorially. If any census exceeds `DP_MAX_STATES`, or if the cross-product between two adjacent censuses exceeds `DP_MAX_STATES²`, the solver falls back to the probabilistic matcher.

**What the default value of 40,000 means in practice:**

| Stems observed | Tracks (K) | States $P(K,n)$ | Fits in 40,000? |
|:-:|:-:|--:|:-:|
| 2 | 4 | 12 | Yes |
| 3 | 5 | 60 | Yes |
| 4 | 6 | 360 | Yes |
| 5 | 7 | 2,520 | Yes |
| 6 | 8 | 20,160 | Yes |
| 7 | 9 | 181,440 | No → probabilistic |
| 8 | 10 | 1,814,400 | No → probabilistic |

**How to choose a value for your data:**
1. Check how many stems your most complex tags have (e.g., `max_obs = max(table(data$Tag, data$CensusID))`)
2. With `K = max_obs + slack_tracks + births`, compute $P(K, max\_obs)$
3. Set `DP_MAX_STATES` above that number to use exact DP, or below it to use the faster probabilistic matcher for those tags
4. Higher values use more memory and time; lower values route more tags through the probabilistic fallback (which is faster but approximate)

See `dp_global/scripts/README.md` for the full CLI flag reference.

## Key Documentation

| Document | Contents |
|----------|----------|
| `dp_global/README.md` | Algorithm details, cost model, data requirements, parameter estimation, fallback mechanisms |
| `dp_global/scripts/README.md` | CLI flags, chunking, resume, example invocations |
| `dp_global/src/README.md` | C++ acceleration API and validation |
| `data_simulation/README.md` | Simulation parameters, biological models, output format |

## Reconstruction Methods

Each observation in the output receives a `ReconstructionMethod` label indicating how its `ReconstructedStemID` was determined:

| Method | Description |
|--------|-------------|
| `given` | Identity known from input `TrueStemID` (anchor census) |
| `dp` | Assigned by the exact DP solver |
| `probabilistic` | Assigned by the probabilistic greedy matching fallback |
| `igraph` | Assigned by the igraph bipartite matching fallback |
| `provisional_dp` | Provisional anchor assigned by DP |
| `provisional_igraph` | Provisional anchor assigned by igraph fallback |
| `none_after_anchor` | Post-anchor row without assignment |

## Conventions

- Driver scripts: `dp_global/scripts/` — `main_cpp.R` (single-tag/small), `main_cpp_chunk.R` (large chunked), `main_cpp_bci.R` (BCI-specific with `withr` bundle sourcing).
- Module internals: `dp_global/R/` — sourced in order by `dp_global_main.R`.
- Output directories: auto-created under `dp_global/output/` (not tracked by git).
