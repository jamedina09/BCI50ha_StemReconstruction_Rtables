#!/usr/bin/env Rscript
# Dev tool: verify dp_functions.rdata loads and C++ compiles
options(error = function() { quit(status = 1) })

cat("[dev/test] starting\n")

tryCatch({
  source(here::here("dp_global", "R", "load_dp_functions.R"))
  cat("[dev/test] sourced loader script OK\n")
}, error = function(e) {
  cat("[dev/test] ERROR sourcing loader script:\n", e$message, "\n")
  quit(status = 1)
})

rdata_file <- here::here("dp_global", "R", "dp_functions.rdata")
if (!file.exists(rdata_file)) {
  cat("[dev/test] ERROR: rdata file not found:", rdata_file, "\n")
  quit(status = 1)
}

# Try loading without compiling
tryCatch({
  obj <- load_dp_functions_rdata(rdata_file, compile_cpp = FALSE)
  cat(sprintf("[dev/test] Loaded rdata: functions=%d, cpp_present=%s, compiled_ok=%s\n",
              length(obj$functions), !is.null(obj$cpp_source), as.character(obj$compiled_ok)))
}, error = function(e) {
  cat("[dev/test] ERROR loading rdata:\n", e$message, "\n")
  quit(status = 1)
})

# Attach into a fresh environment
tryCatch({
  target_env <- new.env()
  attach_dp_functions_from_rdata(rdata_file, envir = target_env, compile_cpp = FALSE)
  cat(sprintf("[dev/test] Attached into new env: objects=%d\n", length(ls(envir = target_env))))
}, error = function(e) {
  cat("[dev/test] ERROR attaching functions:\n", e$message, "\n")
  quit(status = 1)
})

# If there's C++ source, try to compile it now
if (!is.null(obj$cpp_source)) {
  cat("[dev/test] C++ source embedded; attempting compilation via Rcpp::sourceCpp()\n")
  if (!requireNamespace("Rcpp", quietly = TRUE)) {
    cat("[dev/test] Rcpp not installed; cannot compile embedded C++\n")
  } else {
    tryCatch({
      res2 <- load_dp_functions_rdata(rdata_file, compile_cpp = TRUE)
      cat(sprintf("[dev/test] Compilation attempted; compiled_ok=%s\n", as.character(res2$compiled_ok)))
    }, error = function(e) {
      cat("[dev/test] ERROR compiling embedded C++:\n", e$message, "\n")
      quit(status = 1)
    })
  }
}

cat("[dev/test] completed OK\n")
