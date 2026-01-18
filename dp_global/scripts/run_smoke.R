#!/usr/bin/env Rscript
# Simple smoke test: verify core dp_global functions load

suppressPackageStartupMessages({
  library(here)
})

message("[run_smoke] here root: ", here::here())

core_files <- list(
  here::here("dp_global","R","dp_global_biol.R"),
  here::here("dp_global","R","check_functions.r")
)

for (f in core_files) {
  if (file.exists(f)) {
    message("[run_smoke] sourcing: ", f)
    source(f)
  } else {
    stop("Missing core file: ", f)
  }
}

# Quick checks for essential functions
required <- c("estimate_bio_pars","match_stems_dp_global_backward","transition_cost_tracks_bio")
for (fn in required) {
  if (exists(fn, mode = "function")) {
    message("[run_smoke] loaded function: ", fn)
  } else {
    stop("Missing function after sourcing: ", fn)
  }
}

message("[run_smoke] Smoke test passed: core functions loaded successfully.")
