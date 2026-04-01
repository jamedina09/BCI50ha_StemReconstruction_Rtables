# STEM_IDENTIFICATION_TEST

## Overview

Biologically informed dynamic-programming (DP) solver that reconstructs stem identities across forest censuses backward in time from a known anchor census.  Given multi-stem tree measurements and a late-census anchor with trusted `TrueStemID`, the algorithm assigns each earlier observation to a latent identity track by minimising negative log-likelihood costs that encode growth, mortality, and recruitment biology.  Uncertainty is quantified via forward-backward marginals and optional posterior sampling.

## Directory Layout

```
├── dp_global/                 # Core algorithm, drivers, and C++ acceleration
│   ├── R/                     # R modules (sourced by dp_global_main.R)
│   │   ├── complexity/        # DP complexity estimator (predict runtime)
│   │   └── dpglobal_bundle/   # Portable deployment bundle builder
│   ├── scripts/               # CLI driver scripts (main_cpp*.R)
│   └── src/                   # C++ transition cost (Rcpp)
├── data_simulation/           # Simulated forest-census data generator
│   └── data/                  # Generated test datasets (CSV)
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

See `dp_global/scripts/README.md` for the full CLI flag reference.

## Key Documentation

| Document | Contents |
|----------|----------|
| `dp_global/README.md` | Algorithm details, cost model, data requirements, parameter estimation |
| `dp_global/scripts/README.md` | CLI flags, chunking, resume, example invocations |
| `dp_global/src/README.md` | C++ acceleration API and validation |
| `data_simulation/README.md` | Simulation parameters, biological models, output format |

## Conventions

- Driver scripts: `dp_global/scripts/` — `main_cpp.R` (single-tag/small), `main_cpp_chunk.R` (large chunked), `main_cpp_bci.R` (BCI-specific with `withr` bundle sourcing).
- Module internals: `dp_global/R/` — sourced in order by `dp_global_main.R`.
- Output directories: auto-created under `dp_global/output/` (not tracked by git).
