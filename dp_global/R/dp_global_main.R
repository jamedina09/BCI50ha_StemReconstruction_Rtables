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

## ---- 2) C++ acceleration (optional, non-fatal) --------------------------
# Source helper R wrappers (keeps the uncompiled fallback short-circuitable)
# Ensure the R wrapper is defined in the top-level environment so callers can
# access the wrapper regardless of how this loader is invoked.
# Use project root rather than here() to construct file paths (CRAN-friendly style)
root_dir <- getwd()
sys.source(file.path(root_dir, "dp_global", "src", "transition_cost_rcpp.R"), envir = globalenv())
tryCatch({
  Rcpp::sourceCpp(file.path(root_dir, "dp_global", "src", "transition_cost_rcpp.cpp"))
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
  "dp_probabilistic_matching.R",
  "dp_global_dp.R",
  "dp_global_diag.R"
)

# Per-file package requirements (checked before sourcing each file)
required_pkgs_by_file <- list(
  "dp_global_utils.R" = character(0),
  "dp_global_bio.R" = c("data.table", "MASS"),
  "dp_global_states.R" = character(0),
  "dp_global_matchers.R" = c("data.table", "igraph"),
  "dp_probabilistic_matching.R" = c("data.table"),
  "dp_global_dp.R" = c("data.table"),
  "dp_global_diag.R" = c("data.table")
)
# Packages that we prefer to attach (library) because the code uses operators
# or unqualified calls that are convenient when attached.
attach_pkgs <- c("data.table", "MASS")

## ---- 4) Source modules with checks -------------------------------------
for (f in r_files) {
  fp <- file.path(root_dir, "dp_global", "R", f)
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
    sys.source(fp, envir = globalenv()),
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
  "match_stems_probabilistic",
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

## ---- 7) Shared post-engine helper: carried_terminal backfill -----------
# Backfill orphan end-of-trajectory rows whose engine output is NA.
#
# A row is treated as an "orphan terminal" when ALL of the following hold:
#   - ReconstructedStemID is NA after the engine has run
#   - DBH is NA (death/break events typically have no measurement)
#   - Status is one of "dead", "stem dead", "broken below"
#
# For each such row we copy the most recent prior non-NA ReconstructedStemID
# from the same (Tag, OriginalStemID) group (LOCF). Biologically, a terminal
# event ends the trajectory of the most recent prior identity carrying the
# same OriginalStemID; without this fill these rows would be dropped from
# any downstream trajectory.
#
# Returns the (potentially modified) data.table with `ReconstructionMethod`
# set to "carried_terminal" on rows that were filled.
apply_carried_terminal_backfill <- function(out, verbose = TRUE) {
    if (is.null(out) || nrow(out) == 0L) {
        return(out)
    }
    needed <- c("Status", "DBH", "OriginalStemID", "Tag",
                "CensusID", "ReconstructedStemID")
    if (!all(needed %in% names(out))) {
        return(out)
    }
    data.table::setorder(out, Tag, OriginalStemID, CensusID)
    .term_mask <- is.na(out$ReconstructedStemID) &
        is.na(out$DBH) &
        !is.na(out$Status) &
        out$Status %in% c("dead", "stem dead", "broken below")
    if (!any(.term_mask)) {
        return(out)
    }
    out[, .carried := data.table::nafill(ReconstructedStemID, type = "locf"),
        by = .(Tag, OriginalStemID)]
    .fill_mask <- .term_mask & !is.na(out$.carried)
    .n_filled <- sum(.fill_mask)
    if (.n_filled > 0L) {
        if (!("ReconstructionMethod" %in% names(out))) {
            out[, ReconstructionMethod := NA_character_]
        }
        out[.fill_mask, `:=`(
            ReconstructedStemID  = .carried,
            ReconstructionMethod = "carried_terminal"
        )]
        if (isTRUE(verbose)) {
            message(sprintf(
                "[apply_carried_terminal_backfill] backfilled %d orphan terminal-event row(s) from prior same-OriginalStemID Recon.",
                .n_filled
            ))
        }
    }
    out[, .carried := NULL]
    out
}

## ---- 8) Shared post-engine helper: orphan-stem backfill ----------------
# Backfill rows for "born-orphan" stems that the engine cannot reach.
#
# A stem is "born orphan" when its source identifier (StemID in production,
# OriginalStemID in this test repo) first appears with no DBH and no
# upstream TrueStemID anchor (e.g. a brand-new StemID first recorded as
# broken-below at C7+). Without DBH the DP has no signal to disambiguate
# and TrueStemID is never assigned, so the engine leaves
# ReconstructedStemID = NA on every census of that stem.
#
# Rule (post-engine, after carried_terminal backfill):
#   where  is.na(ReconstructedStemID)
#     AND  is.na(TrueStemID)
#     AND  is.na(DBH)
#     AND  source-id column is non-NA
#   then  ReconstructedStemID  := <source-id>
#         ReconstructionMethod := "given_orphan"
#
# Justified because (a) the source ID is unambiguous, (b) DBH=NA means DP
# has no signal, (c) the new method label keeps the trail auditable, and
# (d) collision risk is zero — DP cannot have reached these rows.
apply_orphan_stem_backfill <- function(out, verbose = TRUE) {
    if (is.null(out) || nrow(out) == 0L) {
        return(out)
    }
    src_col <- if ("StemID" %in% names(out)) {
        "StemID"
    } else if ("OriginalStemID" %in% names(out)) {
        "OriginalStemID"
    } else {
        return(out)
    }
    needed <- c(src_col, "TrueStemID", "DBH", "ReconstructedStemID")
    if (!all(needed %in% names(out))) {
        return(out)
    }
    .src <- out[[src_col]]
    .orphan_mask <- is.na(out$ReconstructedStemID) &
        is.na(out$TrueStemID) &
        is.na(out$DBH) &
        !is.na(.src)
    n_orphan <- sum(.orphan_mask)
    if (n_orphan == 0L) {
        return(out)
    }
    if (!("ReconstructionMethod" %in% names(out))) {
        out[, ReconstructionMethod := NA_character_]
    }
    out[.orphan_mask, `:=`(
        ReconstructedStemID  = .src[.orphan_mask],
        ReconstructionMethod = "given_orphan"
    )]
    if (isTRUE(verbose)) {
        message(sprintf(
            "[apply_orphan_stem_backfill] backfilled %d orphan stem row(s) using %s.",
            n_orphan, src_col
        ))
    }
    out
}

# End of dp_global_main.R
# -------------------------------------------------------------------------
