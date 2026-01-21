# -------------------------------------------------------------------------
# dp_global_main.R — Main loader for `dp_global` R modules
# -------------------------------------------------------------------------
# Purpose: central entrypoint that ensures required packages are present,
# attempts to enable C++ acceleration, sources core R modules in a
# deterministic order, and performs quick post-source sanity checks.
# -------------------------------------------------------------------------

## ---- 1) Package preflight checks ----------------------------------------
# Minimal required packages for the DP workflow. We check but do not attach
# everything — only packages that need to be attached for convenient operators
# (e.g., data.table's `:=`) are attached below.
check_pkg <- function(p) {
  if (!requireNamespace(p, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required. Install it with install.packages('%s')", p, p), call. = FALSE)
  }
  invisible(TRUE)
}

check_pkg("data.table")
check_pkg("igraph")
check_pkg("Rcpp")
check_pkg("here")

## ---- 2) C++ acceleration (optional, non-fatal) --------------------------
# Source helper R wrappers (keeps the uncompiled fallback short-circuitable)
source(here::here("dp_global", "src", "transition_cost_rcpp.R"))
tryCatch({
  Rcpp::sourceCpp(here::here("dp_global", "src", "transition_cost_rcpp.cpp"))
  message("[dp_global_main.R] C++ acceleration enabled.")
}, error = function(e) {
  warning(sprintf("C++ compilation (transition_cost_rcpp.cpp) failed or is unavailable: %s. Continuing without compiled acceleration.", e$message))
})

## ---- 3) Module manifest & package requirements --------------------------
# Ordered list of R modules we expect to load. Order chosen to respect
# likely dependencies (utils -> bio -> states -> matchers -> dp -> diag).
r_files <- c(
  "dp_global_utils.R",
  "dp_global_bio.R",
  "dp_global_states.R",
  "dp_global_matchers.R",
  "dp_global_dp.R",
  "dp_global_diag.R"
)

# Per-file package requirements (checked before sourcing each file)
required_pkgs_by_file <- list(
  "dp_global_utils.R" = character(0),
  "dp_global_bio.R" = c("data.table", "MASS"),
  "dp_global_states.R" = character(0),
  "dp_global_matchers.R" = c("data.table", "igraph"),
  "dp_global_dp.R" = c("data.table"),
  "dp_global_diag.R" = c("data.table")
)
# Packages that we prefer to attach (library) because the code uses operators
# or unqualified calls that are convenient when attached.
attach_pkgs <- c("data.table", "MASS")

## ---- 4) Source modules with checks -------------------------------------
for (f in r_files) {
  fp <- here::here("dp_global", "R", f)
  if (!file.exists(fp)) {
    stop(sprintf("Required file not found: %s", fp), call. = FALSE)
  }

  # Ensure packages required by this file are installed and attach if requested
  pkgs <- required_pkgs_by_file[[f]]
  if (!is.null(pkgs) && length(pkgs) > 0) {
    for (p in pkgs) {
      check_pkg(p)
      if (p %in% attach_pkgs && !(p %in% loadedNamespaces())) {
        message(sprintf("[dp_global_main.R] attaching: %s", p))
        suppressMessages(suppressPackageStartupMessages(library(p, character.only = TRUE)))
      }
    }
  }

  message(sprintf("[dp_global_main.R] sourcing: %s", f))
  tryCatch(
    source(fp),
    error = function(e) stop(sprintf("Error sourcing %s: %s", f, e$message), call. = FALSE)
  )
}

## ---- 5) Post-source sanity checks -------------------------------------
# Quick verification that a handful of core functions exist after sourcing.
required_symbols <- c(
  "estimate_bio_pars",
  "transition_cost_tracks_bio_components",
  "match_stems_dp_global_backward_marginals_batch",
  "match_stems_optimal_backward",
  "enumerate_states_injective",
  "add_dp_posterior_bins"
)
missing_symbols <- required_symbols[!vapply(required_symbols, function(s) exists(s, mode = "function", inherits = TRUE), logical(1L))]
if (length(missing_symbols) > 0L) {
  stop(sprintf("After sourcing modules, the following expected functions were missing: %s", paste(missing_symbols, collapse = ", ")), call. = FALSE)
}
message(sprintf("[dp_global_main.R] Sourcing complete; core functions present: %s", paste(required_symbols, collapse = ", ")))

## ---- 6) Namespace housekeeping -----------------------------------------
# Avoid R CMD check notes by declaring commonly used global variables as NULL
Tag <- CensusID <- DBH <- TrueStemID <- ReconstructedStemID <- ConstraintViolation <- ReconstructionMethod <- NULL
ReferenceStemID <- ConstraintViolationFlag <- NULL
DP_MaxStatesPerCensus <- DP_MaxStatesCensusID <- DP_KUsed <- NULL
species <- NULL

# End of dp_global_main.R
# -------------------------------------------------------------------------
