# Work in progress

# Stem Reconstruction and Rtables for the Barro Colorado 50-ha Forest Dynamics Plot

## Overview

This repository has three components:

- **`dp_global/`** — a biologically informed dynamic-programming (DP) engine, with a probabilistic matcher for the tags the DP cannot handle, that reconstructs multi-stem tree identities across forest censuses. Given measurements from a long-term plot and a late "anchor" census from which stem labels are confirmed, the engine assigns each earlier observation to a latent identity track by minimising negative log-likelihood costs that encode growth, mortality, and recruitment (or ingrowth). Posterior path sampling provides uncertainty estimates on downstream derived quantities.

- **`data_simulation/`** — a synthetic forest-census generator used to develop and validate the `dp_global` engine. It produces biologically plausible multi-species, multi-stem datasets with controlled ground truth, including hardcoded edge-case and regression-test tags derived from BCI field data itself.

- **`BCI_stem_reconstruction/`** — an end-to-end pipeline that applies `dp_global` to the Barro Colorado Island (BCI) 50-ha permanent plot across nine censuses (1982–2022/3). It covers ForestGEO data preparation, chunked DP stem reconstruction, posterior consolidation, ForestGEO-format R-table assembly, and estimation of aboveground biomass (AGB) stocks and fluxes and basal area (BA) uncertainty.

---

## Repository Structure

```
├── dp_global/                 # Core algorithm, drivers, and C++ acceleration
│   ├── R/                     # R modules
│   │   ├── dp_global_main.R           # Module loader (compiles the C++ file, sources the engine modules) + post-engine helpers
│   │   ├── dp_global_bio.R            # Biological parameter estimation
│   │   ├── dp_global_states.R         # State enumeration & track-DBH helpers
│   │   ├── dp_probabilistic_matching.R # Probabilistic matcher
│   │   ├── dp_global_dp.R             # Core DP solver (backward/forward pass, marginals)
│   │   ├── dp_global_matchers.R       # Stepwise igraph matcher (not called by the DP workflow)
│   │   ├── dp_global_utils.R          # State-count and interval helpers
│   │   ├── dp_global_diag.R           # Diagnostics & PDF plotting
│   │   ├── measurement_rejoin.R       # Post-engine rejoin of stems split by a moved point of measurement
│   │   ├── naming_helpers.R           # Output directory naming
│   │   ├── check_functions.r, k_tuning_viz.R, realism_calibration.R,
│   │   │   sensitivity_transition_cost_bio.R  # Inspection, tuning and sensitivity helpers
│   │   ├── complexity/                # DP complexity estimator
│   │   └── dpglobal_bundle/           # Portable deployment bundle builder
│   ├── output/                # Runtime outputs (not tracked by git)
│   ├── scripts/               # CLI driver scripts and BA uncertainty
│   ├── src/                   # C++ transition cost (Rcpp)
│   └── tests/                 # testthat tests
├── BCI_stem_reconstruction/   # Full BCI 50-ha pipeline
│   ├── 1_DATA_PREPARATION/    # Species tables and cleaned ViewFullTable
│   ├── 2_STEM_IDENTIFICATION/ # Chunked runner on BCI data + chunk merger
│   ├── 3_PREPARE_R_TABLES/    # Posterior consolidation and ForestGEO-format R tables
│   ├── 4_EXAMPLE_STRUCTURE_ASSESSMENT/  # Plot summaries, AGB stocks/fluxes and BA uncertainty
│   └── DATA/                  # Raw inputs and generated BCI outputs (mostly not tracked by git)
├── data_simulation/           # Simulated forest-census data generator
│   ├── data/                  # Generated test dataset (CSV + diagnostic PDFs)
│   ├── sample_data_BCI/       # BCI multi-stem tags table (input of dp_global/scripts/main_cpp_bci.R)
│   └── simulate_data.R        # Simulation driver script
└── Makefile                   # Convenience target (smoke test)
```

---

## `dp_global/` — Stem Reconstruction Engine

### Prerequisites

R with the packages `data.table`, `Rcpp`, `igraph`, `MASS`, `here`, `lpSolve` and `parallel`, and a C++ compiler (the DP needs the compiled functions of `dp_global/src/transition_cost_rcpp.cpp`).
Optional: `ggplot2`, `cowplot`, `reshape2`, `scales`, `patchwork` (plots and reports), `arrow` (feather output), `testthat` (tests).

Run every command from the project root.

### Quickstart

```bash
# Verify all R modules load
make smoke

# Single-tag run on simulated data
Rscript dp_global/scripts/main_cpp.R --WHICH_TAG=20

# All tags of the simulated data — chunked, with resume support
# (the chunked driver always processes every tag)
Rscript dp_global/scripts/main_cpp_chunk.R --DP_CHUNK_SIZE=7

# BCI single-tag debug run (reads data_simulation/sample_data_BCI/multistem_tags.rds)
Rscript dp_global/scripts/main_cpp_bci.R --WHICH_TAG=123375
```

### Key CLI flags

Defaults are those of `main_cpp.R` and `main_cpp_chunk.R` unless a driver is named.

| Flag | Default | Purpose |
| ------ | --------- | --------- |
| `--DP_MAX_STATES` | `40000` (BCI driver `1_main_cpp_chunk_bci.R`: `10000`) | Max injective states per census before probabilistic fallback |
| `--PROB_N_SAMPLES` | `200` | Gumbel-noise samples for probabilistic matching |
| `--PROB_LOOKAHEAD_WEIGHT` | `1` | Sequential backward conditioning weight (0 = disabled) |
| `--POSTERIOR_SAMPLES` | `200` | Posterior path samples drawn per tag (0 to disable) |
| `--USE_BIO_HARD_SHRINK_IN_PROB` | `TRUE` | Hard shrink gate in probabilistic matcher; set `FALSE` for confirmed large-shrinkage events |
| `--USE_BIO_HARD_GROWTH_IN_PROB` | `TRUE` | Hard growth gate in probabilistic matcher; set `FALSE` to allow exceptional growth |
| `--DBH_ROUND_CENSUSES` | `1,2` (BCI driver `1_main_cpp_chunk_bci.R`) | Censuses whose small-stem DBH (< 55 mm) was recorded in 5 mm classes, rounded down (BCI 1982, 1985); both engines read such a DBH as a size within its class. `none` turns it off |

### How the two algorithms work

#### The problem

In long-term forest census plots, individual trees can have multiple stems measured every few years, but **stem identity labels are only reliable from one late census on** (the "anchor"). For all earlier censuses, we need to determine which measurement belongs to which stem — a problem compounded by stem death and new recruitment.

#### Exact DP solver

The DP solver works backward from the anchor census, evaluating the assignments of earlier measurements to identity tracks (the anchor's stems plus spare tracks for stems that died earlier). Each candidate is scored using biology: size-dependent growth rates, mortality hazard, and recruitment probability. The algorithm picks the single jointly optimal assignment across all censuses.

The number of assignments grows factorially with stem count. With `DP_MAX_STATES = 40,000` the DP handles tags with up to about **6 measured stems per census** (the exact limit depends on the number of tracks; see *Understanding `max_states`* in `dp_global/README.md`, which also describes the case in which the DP works on a truncated set of assignments); larger tags go to the probabilistic matcher.

#### Probabilistic greedy matcher (fallback)

When the DP cannot be used (the state space is too large for exact enumeration, or no assignment satisfies its constraints), the probabilistic matcher draws hundreds of Gumbel-noise-perturbed samples, stitches them backward from the anchor and cuts the links of each sample that break the hard growth limits. In every census pair each stem may continue, die or be recruited (birth-death assignment, `PROB_BIRTH_DEATH = TRUE`): a link is kept only when it is more likely than the death of the earlier stem plus the recruitment of the later one, as in the DP, and each sample is the exact best assignment of a perturbed pair. Database pins (`TrueStemID`) constrain every sample, and two stems with different pins are never joined. The exported reconstruction is the most representative sample (the one whose links agree most with the other samples), and the share of samples that agree with it is the posterior probability per observation. The matcher scores links with the same growth, mortality and recruitment parameters and the same hard growth bounds as the DP; the differences between the two engines are listed in `dp_global/README.md`.

#### When each algorithm runs

| Scenario | Algorithm |
| ---------- | ----------- |
| Assignments per census within `DP_MAX_STATES` (up to about 6 measured stems at 40,000) | Exact DP |
| More assignments than `DP_MAX_STATES` in any census | Probabilistic matcher |
| Species in `PROB_SPECIES` / growth forms in `DP_FALLBACK_GROWTH_FORMS` | Probabilistic matcher |
| DP finds no feasible assignment (e.g. pins it cannot honour) | Probabilistic matcher (automatic fallback) |
| DP hits a runtime error | Probabilistic matcher (automatic fallback) |

For tags that the DP splits at a resprout census (codes R, RP, RF, RT, QR, OR), each segment chooses its algorithm independently. See `dp_global/README.md` for the full algorithm reference including the DP_MAX_STATES state-space tables, biological cost model, measurement error model, and posterior path format.

---

## `data_simulation/` — Test Dataset

`data_simulation/simulate_data.R` generates a synthetic multi-species tropical forest census dataset used to develop and regression-test the `dp_global` engine. The output contains 142 tags:

- **42 simulated-style trees** (Tags 1–42): 38 randomly generated trees across 3 species with species-specific trait scaling, and 4 hand-written tags (39–42).
- **52 hardcoded diagnostic tags** (17 tags in 43–88 and 35 BCI tag numbers) copied from BCI multi-stem patterns for regression testing.
- **3 M-code test tags** (Tags 901–903) whose main stem carries the `M` code (the engines do not read this code).
- **45 row-count invariant edge-case tags** (Tags 9901–9945) covering combinations of census span, stem count, DBH availability, and special flags.

```bash
Rscript data_simulation/simulate_data.R
```

Outputs are written to `data_simulation/data/` (the folder must exist). See `data_simulation/README.md` for the full simulation parameters and column schema.

---

## `BCI_stem_reconstruction/` — BCI 50-ha Pipeline

A sequential four-stage pipeline that takes raw BCI ForestGEO exports through to biomass flux estimates.

### Stage 1 — Data Preparation (`1_DATA_PREPARATION/`)

Converts raw ForestGEO census exports into the species table and a cleaned ViewFullTable. Checks the species names (TNRS), keeps one measurement per stem and census, flags likely DBH entry errors and writes a corrected value in a new column (the recorded `DBH` is kept), applies the Cushman et al. 2014 taper correction (the corrected DBH is the one the engine parameters are estimated from), adds each species' growth form (`Lifeform`, which the stage-2 driver maps to the engine's growth forms), and labels tags as single- or multiple-stem.

Scripts (run in order): `0_prepare_species_tables.R` → `1_prepare_viewfulltable.R.R`

### Stage 2 — Stem Identification (`2_STEM_IDENTIFICATION/`)

Runs `dp_global` on the multiple-stem tags of the cleaned ViewFullTable, chunk by chunk (the tags of a chunk run in parallel), and merges the outputs with the single-stem tags.

- `1_main_cpp_chunk_bci.R` — chunked driver for BCI data; writes feather chunk outputs with resume support.
- `2_merge_chunks_to_datatable.R` — merges chunk feathers into `merged_output.parquet` and `.rds`, joins trunks that the engine split only because their point of measurement moved or one diameter was recorded wrongly (`apply_measurement_rejoin()`, `dp_global/R/measurement_rejoin.R`), and writes the final table `DATA/PROCESSED/complete_dataset_final_with_reconstructed_stemids.rds`.

See `BCI_stem_reconstruction/2_STEM_IDENTIFICATION/run_chunk_bci.md` for run and resume commands.

### Stage 3 — Prepare R Tables (`3_PREPARE_R_TABLES/`)

Consolidates posterior path files and builds ForestGEO-format census R tables.

- `1_prepare_posteriors_BCI.R` — aggregates `_paths.feather` files into `posterior_sampled_paths.rds`, applying the merge step's joins to every sample.
- `2_create_R_tables_BCI.R` — resolves encounter histories into the corrected status (`Rstatus`), gives every tree one location, fills the `date` column (days since 1960-01-01; `ExactDate`, `dbh` and `DFstatus` stay exactly as recorded), and exports `<site>.stemN.Rdata`. The species table is not written here (see `BCI_stem_reconstruction/3_PREPARE_R_TABLES/README.md`).

### Stage 4 — Example Structure Assessment (`4_EXAMPLE_STRUCTURE_ASSESSMENT/`)

Three independent analysis scripts; `biomass_stocks_fluxes.R` and `basal_area_uncertainty.R` write their outputs to `outputs/`, and `general_plot_information.R` prints its summaries to the console.

**`biomass_stocks_fluxes.R`** — estimates AGB stocks, productivity, mortality, and net AGB change across nine BCI censuses using Chave et al. 2014 allometry with Martinez-Cano et al. 2019 height model (trees) and Goodman et al. 2013 (palms). Applies optional strangler-fig removal, palm DBH correction, taper correction, DBH interpolation, 1985 rounding-bias correction, size-class stratification, and Kohyama et al. 2019 productivity/mortality bias correction. Outputs: `outputs/plot_agb_dynamics.png`, `outputs/plot_agb_by_size.png`.

**`basal_area_uncertainty.R`** — propagates stem-identity uncertainty from `dp_global` posterior paths into basal area stocks and fluxes via Monte Carlo realizations. Reports the exported reconstruction and the mean and 2.5–97.5 % quantile interval of the realizations; identities from the anchor census (2010) on are fixed, so the realizations differ before it (and for the rare unmeasured gap that spans it). Outputs: `outputs/fig1_BA_stock.pdf`, `outputs/fig2_BA_fluxes.pdf`, `outputs/fig3_BA_trajectories.pdf`, plus feather tables of the exported reconstruction (`ba_map_*`) and of the realizations (`ba_mc_*`).

**`general_plot_information.R`** — descriptive summaries of forest structure and composition: individuals, stems, basal area, species, genera, families, lifeforms and diversity in the most recent census, bootstrapped per-hectare estimates, temporal trends and changes between censuses, and the DBH size-class distribution.

---

## Engine Output Reference

After the engine and all post-processing helpers run, `ReconstructedStemID` values are renumbered sequentially from 1 to N within each tag, ordered by the earliest census in which each stem appears (ties broken by largest DBH, then original ID). The BCI merge step's measurement rejoin (`apply_measurement_rejoin()`) runs after this renumbering: a joined stem keeps one of its two IDs, so a tree it changes skips one number per join.

---

## Key Documentation

| Document | Contents |
| ---------- | ---------- |
| `dp_global/README.md` | Algorithm details, cost model, data requirements, parameter estimation, fallback mechanisms |
| `dp_global/scripts/README.md` | CLI flags, chunking, resume, example invocations, basal area uncertainty |
| `dp_global/src/README.md` | C++ acceleration API and validation |
| `dp_global/R/dpglobal_bundle/README.md` | Portable bundle builder for deploying the algorithm on other machines |
| `data_simulation/README.md` | Simulation parameters, biological models, output format |
| `BCI_stem_reconstruction/README.md` | The four stages, data flow and output locations |
| `BCI_stem_reconstruction/1_DATA_PREPARATION/README.md` | Build species tables and cleaned ViewFullTable from ForestGEO exports |
| `BCI_stem_reconstruction/2_STEM_IDENTIFICATION/README.md` | Chunked BCI runner and chunk merger (`run_chunk_bci.md`: run and resume commands) |
| `BCI_stem_reconstruction/3_PREPARE_R_TABLES/README.md` | Consolidate posteriors and build ForestGEO-format R tables |
| `BCI_stem_reconstruction/4_EXAMPLE_STRUCTURE_ASSESSMENT/README.md` | Plot summaries, AGB stocks/fluxes and basal-area uncertainty for the BCI 50-ha plot |
