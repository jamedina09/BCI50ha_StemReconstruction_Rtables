# dp_global/R/dpglobal_bundle

Packages the `dp_global` engine into one archive, to run it on another machine.

## Contents

| File | What it does |
| --- | --- |
| `package_bundle.sh` | Builds `dist/dpglobal_bundle_full_<timestamp>.tar.gz` and its SHA256 checksum |
| `dpglobal_bundle_loader.R` | Writes `dpglobal_bundle.RData` (the engine's R objects, sourced) and `dpglobal_bundle_manifest.rds` (the R packages each file needs), next to it |
| `verify_bundle.R` | Checks that the modules load and the core functions exist |
| `dist/` | Archives built by `package_bundle.sh`. An archive holds the code of the day it was built: rebuild before deploying |

`dist/`, `dpglobal_bundle.RData` and `dpglobal_bundle_manifest.rds` are build
products and are not tracked by git (`.gitignore`).

## How to run

Build, from the project root of the source machine:

```sh
bash dp_global/R/dpglobal_bundle/package_bundle.sh --build-bundle   # rebuilds the RData and the manifest first
bash dp_global/R/dpglobal_bundle/package_bundle.sh                  # packages only
```

Deploy, on the target machine:

```sh
mkdir -p ~/dp_global_bundle
tar -xzf dpglobal_bundle_full_<timestamp>.tar.gz -C ~/dp_global_bundle
cd ~/dp_global_bundle
```

```r
manifest <- readRDS("dp_global/R/dpglobal_bundle/dpglobal_bundle_manifest.rds")
pkgs <- unique(unlist(manifest$required_pkgs_by_file))
pkgs <- pkgs[!is.na(pkgs) & nzchar(pkgs) & !pkgs %in% installed.packages()[, "Package"]]
if (length(pkgs)) install.packages(pkgs)
install.packages("lpSolve")                      # used by the matcher; not in the manifest

source(file.path("dp_global", "R", "dp_global_main.R"))   # compiles the C++ code and loads the engine
```

```sh
Rscript dp_global/R/dpglobal_bundle/verify_bundle.R      # optional check
```

Run a reconstruction with a driver, e.g. `Rscript dp_global/scripts/main_cpp.R`
(`dp_global/scripts/README.md`).

## Outputs

The archive holds the whole `dp_global/` folder (without `dp_global/output/`,
`dist/` and `.git/`), `data_simulation/data/` as example input, `INSTALL.txt`,
and this README, `verify_bundle.R`, `dpglobal_bundle.RData` and the manifest
at its root. It holds source code, not compiled code.

`dpglobal_bundle.RData` is optional: sourcing `dp_global_main.R` loads
everything from the R files. It does not hold the compiled C++ functions.

## Requirements

R (4.0 or later), the packages `Rcpp`, `here`, `data.table`, `MASS`, `igraph`,
`lpSolve` (the drivers also use `arrow`, `ggplot2`, `cowplot`, `reshape2`,
`scales`, `patchwork`), and a C++ compiler (macOS: `xcode-select --install`).
Compile on the target machine; compiled `.so` or `.dll` files do not move
between machines. `package_bundle.sh` uses `rsync` (or `cp -r`).
