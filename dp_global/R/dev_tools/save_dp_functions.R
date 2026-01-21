#!/usr/bin/env Rscript
# Dev tool: create or refresh dp_functions.rdata from source modules
options(error = function() { quit(status = 1) })

cat("[dev/save] starting\n")

tryCatch({
  source(here::here("dp_global", "R", "load_dp_functions.R"))
  cat("[dev/save] sourced loader script OK\n")
}, error = function(e) {
  cat("[dev/save] ERROR sourcing loader script:\n", e$message, "\n")
  quit(status = 1)
})

tryCatch({
  loader <- load_dp_functions(try_compile_cpp = TRUE)
  out_file <- save_dp_functions_rdata(loader)
  cat(sprintf("[dev/save] Wrote: %s\n", out_file))
}, error = function(e) {
  cat("[dev/save] ERROR saving dp_functions.rdata:\n", e$message, "\n")
  quit(status = 1)
})

cat("[dev/save] completed OK\n")
