# STEM_IDENTIFICATION_TEST

## Overview

This repository has three components:

- **`dp_global/`** — a biologically informed dynamic-programming (DP) engine that reconstructs multi-stem tree identities across forest censuses. Given measurements from a long-term plot and a late "anchor" census with confirmed stem labels, the engine assigns each earlier observation to a latent identity track by minimising negative log-likelihood costs that encode growth, mortality, recruitment, and measurement error. Posterior path sampling provides uncertainty estimates on all downstream derived quantities.

- **`data_simulation/`** — a synthetic forest-census generator used to develop and validate the `dp_global` engine. It produces biologically plausible multi-species, multi-stem datasets with controlled ground truth, including hardcoded edge-case and regression-test tags derived from BCI field data.

- **`BCI_stem_reconstruction/`** — an end-to-end pipeline that applies `dp_global` to the Barro Colorado Island (BCI) 50-ha permanent plot across nine censuses (1982–2022/3). It covers ForestGEO data preparation, chunked DP stem reconstruction, posterior consolidation, ForestGEO-format R-table assembly, and estimation of aboveground biomass (AGB) stocks and fluxes and basal area (BA) uncertainty.

---

## Repository Structure

```
├── dp_global/                 # Core algorithm, drivers, and C++ acceleration
│   ├── R/                     # R modules (sourced by dp_global_main.R)
│   │   ├── dp_global_main.R           # Module loader (sources all R modules in order)
│   │   ├── dp_global_bio.R            # Biological parameter estimation
│   │   ├── dp_global_states.R         # State enumeration & track-DBH helpers
│   │   ├── dp_probabilistic_matching.R # Probabilistic greedy matching fallback
│   │   ├── dp_global_dp.R             # Core DP solver (backward/forward pass, marginals)
│   │   ├── dp_global_utils.R          # Shared utilities
│   │   ├── dp_global_diag.R           # Diagnostics & PDF plotting
│   │   ├── naming_helpers.R           # Output directory naming
│   │   ├── complexity/                # DP complexity estimator
│   │   └── dpglobal_bundle/           # Portable deployment bundle builder
│   ├── output/                 # Runtime outputs (not tracked by git)
│   ├── scripts/               # CLI driver scripts, batch runner, and BA uncertainty
│   └── src/                   # C++ transition cost (Rcpp)
├── BCI_stem_reconstruction/   # Full BCI 50-ha pipeline
│   ├── 1_DATA_PREPARATION/    # Species tables and cleaned ViewFullTable
│   ├── 2_STEM_IDENTIFICATION/ # Chunked DP runner on BCI data + chunk merger
│   ├── 3_PREPARE_R_TABLES/    # Posterior consolidation and ForestGEO-format R tables
│   ├── 4_EXAMPLE_STRUCTURE_ASSESSMENT/  # AGB stocks/fluxes and BA uncertainty
│   └── DATA/                  # Generated BCI outputs (not tracked by git)
├── data_simulation/           # Simulated forest-census data generator
│   ├── data/                  # Generated test datasets (CSV + diagnostic PDFs)
│   ├── sample_data_BCI/       # Example BCI-formatted sample data
│   └── simulate_data.R        # Simulation driver script
└── Makefile                   # Convenience targets (smoke test)
```

---

## `dp_global/` — Stem Reconstruction Engine

### Prerequisites

R ≥ 4.0 with packages: `data.table`, `Rcpp`, `here`.
Optional: `ggplot2`, `cowplot` (plotting), `arrow` (feather output), `withr` (bundle sourcing).

### Quickstart

```bash
# Verify all R modules load
make smoke

# Single-tag run on simulated data
Rscript dp_global/scripts/main_cpp.R --WHICH_TAG=20

# Full run — chunked, with resume support
Rscript dp_global/scripts/main_cpp_chunk.R --RUN_ALL_TAGS=TRUE --DP_CHUNK_SIZE=7

# BCI single-tag debug run
Rscript dp_global/scripts/main_cpp_bci.R --WHICH_TAG=123375
```

### Key CLI flags

| Flag | Default | Purpose |
|------|---------|---------|
| `--DP_MAX_STATES` | `40000` | Max injective states per census before probabilistic fallback |
| `--PROB_N_SAMPLES` | `200` | Gumbel-noise samples for probabilistic matching |
| `--PROB_LOOKAHEAD_WEIGHT` | `1` | Sequential backward conditioning weight (0 = disabled) |
| `--POSTERIOR_SAMPLES` | `200` | Posterior path samples drawn per tag (0 to disable) |
| `--USE_BIO_HARD_SHRINK_IN_PROB` | `TRUE` | Hard shrink gate in probabilistic matcher; set `FALSE` for confirmed large-shrinkage events |
| `--USE_BIO_HARD_GROWTH_IN_PROB` | `TRUE` | Hard growth gate in probabilistic matcher; set `FALSE` to allow exceptional growth |

### How the two algorithms work

#### The problem

In long-term forest census plots, individual trees can have multiple stems measured every few years, but **stem identity labels are only reliable at one late census** (the "anchor"). For all earlier censuses, we need to determine which measurement belongs to which stem — a problem compounded by stem death, new recruitment, and measurement noise.

#### Exact DP solver

The DP solver works backward from the anchor census, evaluating every possible assignment of earlier measurements to known stem identities. Each candidate is scored using biology: size-dependent growth rates, mortality hazard, recruitment probability, and measurement error. The algorithm picks the single jointly optimal assignment across all censuses.

Because it examines every possibility, the DP always finds the mathematically optimal answer. The downside is that the number of possibilities grows factorially with stem count: the solver is exact for tags with **≤ 6 stems per census** (with `DP_MAX_STATES = 40,000`), and falls back for tags with 7 or more stems.

#### Probabilistic greedy matcher (fallback)

When the state space is too large for exact enumeration, the probabilistic matcher draws hundreds of Gumbel-noise-perturbed samples, stitches them backward from the anchor, repairs biological constraint violations, and votes across surviving samples. The vote share becomes the posterior probability per observation. The same biological cost model and pruning bounds used by the DP apply here.

#### When each algorithm runs

| Scenario | Algorithm |
|----------|-----------|
| ≤ 6 observed stems per census | Exact DP |
| 7+ observed stems in any census | Probabilistic matcher |
| Species / growth forms in `FALLBACK_GROWTH_FORMS` | Probabilistic matcher |
| DP hits a runtime error | Probabilistic matcher (automatic fallback) |

For tags split by an R-event (resprout/breakage codes R, RP, RF, RT, QR, OR), each segment chooses its algorithm independently. See `dp_global/README.md` for the full algorithm reference including the DP_MAX_STATES state-space tables, biological cost model, measurement error model, and posterior path format.

---

## `data_simulation/` — Test Dataset

`data_simulation/simulate_data.R` generates a synthetic multi-species tropical forest census dataset used to develop and regression-test the `dp_global` engine. The output contains 136 tags:

- **42 simulated trees** (Tags 1–42) across 3 species with species-specific growth scaling.
- **46 hardcoded diagnostic tags** (Tags 43–88) derived from real BCI multi-stem patterns for regression testing.
- **3 M-code test tags** (Tags 901–903) that validate the M-coded main-stem constraint.
- **45 row-count invariant edge-case tags** (Tags 9901–9945) covering every combination of census span, stem count, DBH availability, and special flags.

```bash
Rscript data_simulation/simulate_data.R
```

Outputs are written to `data_simulation/data/`. See `data_simulation/README.md` for the full simulation parameters and column schema.

---

## `BCI_stem_reconstruction/` — BCI 50-ha Pipeline

A sequential four-stage pipeline that takes raw BCI ForestGEO exports through to biomass flux estimates.

### Stage 1 — Data Preparation (`1_DATA_PREPARATION/`)

Converts raw ForestGEO census exports into a cleaned, harmonized ViewFullTable. Builds species lookup tables, applies Cushman et al. 2014 taper corrections, fixes common data-entry issues, and assigns growth forms used by the DP engine.

Scripts (run in order): `0_prepare_species_tables.R` → `1_prepare_viewfulltable.R.R`

### Stage 2 — Stem Identification (`2_STEM_IDENTIFICATION/`)

Runs `dp_global` on the cleaned ViewFullTable in parallel chunks and merges outputs.

- `1_main_cpp_chunk_bci.R` — chunked DP driver for BCI data; writes feather chunk outputs with resume support.
- `2_merge_chunks_to_datatable.R` — merges chunk feathers into `merged_output.parquet` and `.rds`.

See `BCI_stem_reconstruction/2_STEM_IDENTIFICATION/run_chunk_bci.md` for run and resume commands.

### Stage 3 — Prepare R Tables (`3_PREPARE_R_TABLES/`)

Consolidates posterior path files and builds ForestGEO-format census and species R tables.

- `1_prepare_posteriors_BCI.R` — aggregates `_paths.feather` files into `posterior_sampled_paths.rds`.
- `2_create_R_tables_BCI.R` — resolves encounter histories, applies broken-below rules, imputes missing data, and exports `<site>.stemN.Rdata` and `<site>.spptable.rdata`.

### Stage 4 — Biomass Stocks and Fluxes (`4_EXAMPLE_STRUCTURE_ASSESSMENT/`)

Two independent analysis scripts; all outputs are written to `outputs/`.

**`biomass_stocks_fluxes.R`** — estimates AGB stocks, productivity, mortality, and net AGB change across nine BCI censuses using Chave et al. 2014 allometry with Martinez-Cano et al. 2019 height model (trees) and Goodman et al. 2013 (palms). Applies optional strangler-fig removal, palm DBH correction, taper correction, DBH interpolation, 1985 rounding-bias correction, size-class stratification, and Kohyama et al. 2019 productivity/mortality bias correction. Outputs: `outputs/plot_agb_dynamics.png`, `outputs/plot_agb_by_size.png`.

**`basal_area_uncertainty.R`** — propagates stem-identity uncertainty from `dp_global` posterior paths into basal area stocks and fluxes via Monte Carlo realizations. Reports MAP estimates and empirical 95 % CIs; uncertainty is non-zero only for pre-anchor census intervals. Outputs: `outputs/fig1_BA_stock.pdf`, `outputs/fig2_BA_fluxes.pdf`, `outputs/fig3_BA_trajectories.pdf`, plus MAP and MC feather tables.

---

## Engine Output Reference

After the engine and all post-processing helpers run, `ReconstructedStemID` values are renumbered sequentially from 1 to N within each tag, ordered by the earliest census in which each stem appears (ties broken by largest DBH, then original ID). IDs are always positive and contiguous.

---

## Key Documentation

| Document | Contents |
|----------|----------|
| `dp_global/README.md` | Algorithm details, cost model, data requirements, parameter estimation, fallback mechanisms |
| `dp_global/scripts/README.md` | CLI flags, chunking, resume, example invocations, basal area uncertainty, batch runner |
| `dp_global/src/README.md` | C++ acceleration API and validation |
| `dp_global/R/dpglobal_bundle/README.md` | Portable bundle builder for deploying the algorithm on other machines |
| `data_simulation/README.md` | Simulation parameters, biological models, output format |
| `BCI_stem_reconstruction/1_DATA_PREPARATION/README.md` | Build species tables and cleaned ViewFullTable from ForestGEO exports |
| `BCI_stem_reconstruction/2_STEM_IDENTIFICATION/README.md` | Chunked BCI DP runner and chunk merger |
| `BCI_stem_reconstruction/3_PREPARE_R_TABLES/README.md` | Consolidate posteriors and build ForestGEO-format R tables |
| `BCI_stem_reconstruction/4_EXAMPLE_STRUCTURE_ASSESSMENT/README.md` | AGB stocks/fluxes and basal-area uncertainty for the BCI 50-ha plot |
