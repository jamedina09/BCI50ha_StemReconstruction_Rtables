# STEM_IDENTIFICATION_TEST

## Overview

Biologically informed dynamic-programming (DP) solver that reconstructs stem identities across forest censuses backward in time from a known anchor census. Given multi-stem tree measurements and a late-census anchor with trusted `TrueStemID`, the algorithm assigns each earlier observation to a latent identity track by minimising negative log-likelihood costs that encode growth, mortality, and recruitment biology. Uncertainty is quantified via forward-backward marginals and optional posterior sampling.

When the exact DP state space is too large (combinatorial explosion from many stems per tree), the solver automatically falls back to a **probabilistic greedy matching** module that uses the same biological cost model and hard pruning bounds with Gumbel-noise stochastic sampling to produce approximate reconstructions with per-observation posterior probabilities.

## Directory Layout

```
├── dp_global/                 # Core algorithm, drivers, and C++ acceleration
│   ├── R/                     # R modules (sourced by dp_global_main.R)
│   │   ├── dp_global_main.R           # Module loader (sources all R modules in order)
│   │   ├── dp_global_bio.R            # Biological parameter estimation
│   │   ├── dp_global_states.R         # State enumeration & track-DBH helpers
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

R ≥ 4.0 with packages: `data.table`, `Rcpp`, `here`.
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

`DP_MAX_STATES` controls the maximum number of assignment states the DP solver will enumerate at any single census before falling back to the probabilistic greedy matcher. It also controls the inter-census transition budget: the cross-product of states between any two adjacent censuses must not exceed `DP_MAX_STATES²`.

#### How states are counted

At each census, the DP enumerates all injective (one-to-one) assignments of $n$ observed stems to $K$ identity tracks. The number of such assignments is the falling factorial:

$$P(K, n) = K \times (K-1) \times \cdots \times (K - n + 1) = \frac{K!}{(K-n)!}$$

where $K$ is the number of tracks (determined by the anchor stem count, births needed, and slack). Typically $K = n + 1$ (with `slack_tracks = 1` and no births) or $K = n + 2$ (with 1 birth track added).

This grows **factorially**, so even modest increases in stem count cause explosive growth in the state space.

#### Two fallback triggers

1. **Per-census enumeration (`enum_exceeded`):** If $P(K, n) > \text{DP\_MAX\_STATES}$ at any single census, the solver cannot enumerate states and falls back.
2. **Inter-census transitions (`edge_count_exceeded`):** If $P(K, n_1) \times P(K, n_2) > \text{DP\_MAX\_STATES}^2$ for any pair of adjacent censuses, the transition matrix is too large and the solver falls back.

In practice, the per-census limit is reached first because the state counts grow so rapidly.

#### Fallback thresholds by `DP_MAX_STATES` value

The tables below show when fallback occurs for different `DP_MAX_STATES` values. "Stems observed" is the number of stems with non-NA DBH in a single census. $K = n + 1$ assumes `slack_tracks = 1` with no birth tracks needed.

**With $K = n + 1$ (minimum realistic tracks):**

| Stems ($n$) | Tracks ($K$) | States $P(K,n)$ | 1,000 | 20,000 | 40,000 |
|:-:|:-:|--:|:-:|:-:|:-:|
| 2 | 3 | 6 | DP | DP | DP |
| 3 | 4 | 24 | DP | DP | DP |
| 4 | 5 | 120 | DP | DP | DP |
| 5 | 6 | 720 | DP | DP | DP |
| 6 | 7 | 5,040 | fallback | DP | DP |
| 7 | 8 | 40,320 | fallback | fallback | fallback |
| 8 | 9 | 362,880 | fallback | fallback | fallback |

**With $K = n + 2$ (when 1 birth track is needed):**

| Stems ($n$) | Tracks ($K$) | States $P(K,n)$ | 1,000 | 20,000 | 40,000 |
|:-:|:-:|--:|:-:|:-:|:-:|
| 2 | 4 | 12 | DP | DP | DP |
| 3 | 5 | 60 | DP | DP | DP |
| 4 | 6 | 360 | DP | DP | DP |
| 5 | 7 | 2,520 | DP | DP | DP |
| 6 | 8 | 20,160 | fallback | DP | DP |
| 7 | 9 | 181,440 | fallback | fallback | fallback |
| 8 | 10 | 1,814,400 | fallback | fallback | fallback |

**Summary — maximum stems per census handled by exact DP:**

| `DP_MAX_STATES` | `max_edges` ($= \text{DP\_MAX\_STATES}^2$) | Max stems ($K = n+1$) | Max stems ($K = n+2$) |
|--:|--:|:-:|:-:|
| 1,000 | 1,000,000 | 5 | 5 |
| 20,000 | 400,000,000 | 6 | 6 |
| 40,000 | 1,600,000,000 | 6 | 6 |

**Key insight:** With the default `DP_MAX_STATES = 40,000`, the DP handles tags with up to **6 observed stems per census** exactly. Tags with **7 or more stems** in any census are routed to the probabilistic greedy matcher. Increasing `DP_MAX_STATES` to 50,000 would not help — the next factorial step (40,320 for 7 stems with $K=8$) requires `DP_MAX_STATES ≥ 40,321` AND the inter-census product must fit, which it does since $40{,}320^2 = 1.6 \times 10^9 < 40{,}321^2$.

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

**Trade-off:** Higher values → exact DP for more tags (slower, more memory). Lower values → more tags use the probabilistic fallback (faster, approximate but uses the same biological model and pruning bounds).

See `dp_global/README.md` for the full algorithm description and `dp_global/scripts/README.md` for the CLI flag reference.

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
| `provisional_dp` | Provisional anchor assigned by DP |
| `dp_mf_inferred` | Missing-from-field census identity inferred from flanking DP assignments |
| `none_after_anchor` | Post-anchor row without assignment |
| `skipped_no_data` | Tag had no usable data for reconstruction |

## Conventions

- Driver scripts: `dp_global/scripts/` — `main_cpp.R` (single-tag/small), `main_cpp_chunk.R` (large chunked), `main_cpp_bci.R` (BCI-specific with `withr` bundle sourcing).
- Module internals: `dp_global/R/` — sourced in order by `dp_global_main.R`.
- Output directories: auto-created under `dp_global/output/` (not tracked by git).
