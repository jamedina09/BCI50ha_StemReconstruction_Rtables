# dp_global/R/dpglobal_bundle

Packages the `dp_global` engine into one archive, to run it on another
machine. The archive holds source code; the C++ functions are compiled on the
target machine.

## Table of contents

- [Contents](#contents)
- [How to run](#how-to-run)
  - [Build (source machine)](#build-source-machine)
  - [Deploy (target machine)](#deploy-target-machine)
- [Outputs](#outputs)
- [Checks](#checks)
- [Requirements](#requirements)

## Contents

| File | What it does |
| --- | --- |
| `package_bundle.sh` | Builds `dist/dpglobal_bundle_full_<timestamp>.tar.gz` and its SHA256 checksum (bash) |
| `dpglobal_bundle_loader.R` | Sources the engine and writes, next to it, `dpglobal_bundle.RData` (the engine's R objects) and `dpglobal_bundle_manifest.rds` (the R packages each file needs) |
| `verify_bundle.R` | Smoke test: the R modules load and the core functions exist |
| `dist/` | Archives built by `package_bundle.sh`. An archive holds the code of the day it was built: rebuild before deploying |

`dist/`, `dpglobal_bundle.RData` and `dpglobal_bundle_manifest.rds` are build
products and are not tracked by git (`.gitignore`).

## How to run

### Build (source machine)

From the project root:

```sh
bash dp_global/R/dpglobal_bundle/package_bundle.sh --build-bundle   # rebuilds the RData and the manifest, then packages
bash dp_global/R/dpglobal_bundle/package_bundle.sh                  # packages only (RData and manifest as they are)
```

`--build-bundle` runs `dpglobal_bundle_loader.R`, which sources
`dp_global_main.R` and therefore needs the R packages of the engine and a C++
toolchain. Output:

```text
dp_global/R/dpglobal_bundle/dist/dpglobal_bundle_full_YYYYMMDD_HHMMSS.tar.gz
dp_global/R/dpglobal_bundle/dist/dpglobal_bundle_full_YYYYMMDD_HHMMSS.tar.gz.sha256
```

### Deploy (target machine)

1. Unpack:

   ```sh
   mkdir -p ~/dp_global_bundle
   tar -xzf dpglobal_bundle_full_<timestamp>.tar.gz -C ~/dp_global_bundle
   cd ~/dp_global_bundle
   ```

2. Install the missing R packages:

   ```r
   manifest <- readRDS("dp_global/R/dpglobal_bundle/dpglobal_bundle_manifest.rds")
   pkgs <- unique(unlist(manifest$required_pkgs_by_file))
   pkgs <- pkgs[!is.na(pkgs) & nzchar(pkgs) & !pkgs %in% installed.packages()[, "Package"]]
   if (length(pkgs)) install.packages(pkgs)
   install.packages("lpSolve")   # used by the matcher; not in the manifest
   ```

3. Compile the C++ functions (sourcing `dp_global_main.R` does it too; the DP
   cannot run without them):

   ```r
   Rcpp::sourceCpp("dp_global/src/transition_cost_rcpp.cpp")
   ```

4. Check the bundle (optional):

   ```sh
   Rscript dp_global/R/dpglobal_bundle/verify_bundle.R
   ```

5. Load the engine. `dp_global_main.R` finds its files from the working
   directory; `withr::with_dir()` sets it for the call and restores it:

   ```r
   withr::with_dir("~/dp_global_bundle", source(file.path("dp_global", "R", "dp_global_main.R")))
   ```

   For a full reconstruction run a driver, e.g.
   `Rscript dp_global/scripts/main_cpp.R` (`dp_global/scripts/README.md`).

## Outputs

The archive holds:

- the whole `dp_global/` folder, without `dp_global/output/`, the bundle's
  `dist/` and `.git/` (copied with `rsync`, else `cp -r`), including the
  drivers, `basal_area_uncertainty.R`, the probabilistic matcher and the tests;
- `data_simulation/data/` as example input;
- at the archive root: `INSTALL.txt` (generated: how to load the engine and
  compile the C++ file), this README, `verify_bundle.R`,
  `dpglobal_bundle.RData` and `dpglobal_bundle_manifest.rds` (when present);
- extra copies of `src/transition_cost_rcpp.R` and
  `src/transition_cost_rcpp.cpp` under `dp_global/R/dpglobal_bundle/` (the
  path `INSTALL.txt` gives for the manual compilation).

`dpglobal_bundle.RData` is a convenience snapshot and not required: sourcing
`dp_global_main.R` loads everything from the R files. It does not hold the
compiled C++ functions. Its content can be checked directly:

```r
myenv <- new.env()
load("dp_global/R/dpglobal_bundle/dpglobal_bundle.RData", envir = myenv)
exists("estimate_bio_pars", envir = myenv, mode = "function")
exists("match_stems_probabilistic", envir = myenv, mode = "function")
```

## Checks

`package_bundle.sh` checks that the staged copy holds the drivers
(`main_cpp.R`, `main_cpp_chunk.R`, `basal_area_uncertainty.R`), the engine
modules (`dp_global_main.R`, `dp_global_dp.R`, `dp_probabilistic_matching.R`),
the C++ file and `verify_bundle.R`, and warns when one is missing (the archive
is still written). `verify_bundle.R` sources the modules and checks the core
functions (`estimate_bio_pars`, `match_stems_dp_global_backward_marginals_batch`,
`match_stems_optimal_backward`, `match_stems_probabilistic`,
`enumerate_states_injective`, `add_dp_posterior_bins`) and the RData file.

Troubleshooting:

- **Functions missing after loading the RData**: rebuild the bundle on the
  source machine and read the output of `dpglobal_bundle_loader.R` for
  sourcing errors.
- **Missing compiled symbol**: compile with
  `Rcpp::sourceCpp("dp_global/src/transition_cost_rcpp.cpp")` (needs `Rcpp`
  and a C++ toolchain).
- **Shared objects**: do not copy compiled `.so` or `.dll` files between
  machines; compile on the target.
- **Reproducibility**: record the R and package versions
  (`sessioninfo::session_info()`).

## Requirements

R (4.0 or later; preferably the same minor version on both machines), the
packages `Rcpp`, `here`, `data.table`, `MASS`, `igraph`, `lpSolve` (the
drivers also use `arrow`, `ggplot2`, `cowplot`, `reshape2`, `scales`,
`patchwork`; the commands above use `withr`), and a C++ compiler (macOS:
`xcode-select --install`). `package_bundle.sh` needs bash and `rsync` (or
`cp`), and `shasum` or `sha256sum` for the checksum.
