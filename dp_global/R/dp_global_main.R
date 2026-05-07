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
    src_col <- if ("StemID" %in% names(out)) {
        "StemID"
    } else if ("OriginalStemID" %in% names(out)) {
        "OriginalStemID"
    } else {
        return(out)
    }
    needed <- c("Status", "DBH", src_col, "Tag",
                "CensusID", "ReconstructedStemID")
    if (!all(needed %in% names(out))) {
        return(out)
    }
    data.table::setorderv(out, c("Tag", src_col, "CensusID"))
    .term_mask <- is.na(out$ReconstructedStemID) &
        is.na(out$DBH) &
        !is.na(out$Status) &
        out$Status %in% c("dead", "stem dead", "broken below")
    if (!any(.term_mask)) {
        return(out)
    }
    out[, .carried := data.table::nafill(ReconstructedStemID, type = "locf"),
        by = c("Tag", src_col)]
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
                "[apply_carried_terminal_backfill] backfilled %d orphan terminal-event row(s) from prior same-%s Recon.",
                .n_filled, src_col
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

## ---- 9) Shared post-engine helper: broken-below invariant pass --------
# Enforce two invariants on the final reconstruction:
#
#   R1 (split-on-break): within a (Tag, series) group ordered by CensusID,
#       a row with Status == "broken below" and !is.na(DBH) MUST have a
#       ReconstructedStemID that does not equal any prior row's
#       ReconstructedStemID in the same series.  When violated, mint a
#       fresh ID (max(existing)+1, monotonic per-call) and propagate it
#       forward through subsequent rows of the series that currently share
#       the pre-split ID, until the next break event or the next pinned
#       TrueStemID row.  Tag method = "bb_split" (first row) /
#       "bb_split_carry" (carried forward).
#
#   R2 (terminate-on-stump): a row with Status == "broken below" and
#       is.na(DBH) terminates the trajectory.  Any later row in the same
#       series that has !is.na(DBH) and currently shares the terminator's
#       ReconstructedStemID MUST be re-IDed.  Mint a fresh ID and propagate
#       analogously.  Tag method = "bb_post_terminator_split" /
#       "bb_post_terminator_split_carry".  NA-DBH dead/BB/stem-dead corpse
#       rows already labelled by `apply_carried_terminal_backfill` keep
#       their ID (LOCF on terminator is allowed).
#
# Pinned rows (non-NA TrueStemID) are normally respected, BUT when a pin
# conflicts with the contract (typical case: BCI driver pre-stamps
# TrueStemID = OriginalStemID on every BB+DBH row under the assumption
# they begin a new trajectory; ~28% of the time the OriginalStemID is
# reused from a prior alive row in the same StemTag and the pin therefore
# locks in a contract violation), the pass overrides the pin and updates
# both ReconstructedStemID and TrueStemID to the newly-minted ID.  Each
# such override is counted and reported.  The pass is deterministic and
# idempotent: running it twice on the same input yields the same output
# as one run.
#
# Series key: prefer StemTag (mirrors `broken_below_tags.csv` diagnostic)
# then OriginalStemID then StemID.
#
# This function operates only on the MAP-level table.  Posterior CSVs need
# the same operator applied per-sample with `path_sig` recomputation; that
# is handled separately at the posterior writer site.
apply_broken_below_invariants <- function(out, verbose = TRUE) {
    if (is.null(out) || nrow(out) == 0L) {
        return(out)
    }
    needed <- c("Tag", "CensusID", "Status", "DBH", "ReconstructedStemID")
    if (!all(needed %in% names(out))) {
        return(out)
    }
    series_col <- if ("StemTag" %in% names(out)) {
        "StemTag"
    } else if ("OriginalStemID" %in% names(out)) {
        "OriginalStemID"
    } else if ("StemID" %in% names(out)) {
        "StemID"
    } else {
        return(out)
    }
    data.table::setorderv(out, c("Tag", series_col, "CensusID"))
    if (!("ReconstructionMethod" %in% names(out))) {
        out[, ReconstructionMethod := NA_character_]
    }
    cur_max <- suppressWarnings(max(out$ReconstructedStemID, na.rm = TRUE))
    if (!is.finite(cur_max)) cur_max <- 0L
    next_id <- as.integer(cur_max) + 1L

    has_true <- "TrueStemID" %in% names(out)
    rid    <- as.integer(out$ReconstructedStemID)
    mth    <- as.character(out$ReconstructionMethod)
    sts    <- out$Status
    dbh    <- out$DBH
    trueid <- if (has_true) as.integer(out$TrueStemID) else rep(NA_integer_, nrow(out))

    # Group row indices by (Tag, series) — preserves the setorderv ordering.
    groups <- out[, .(.idxs = list(.I)), by = c("Tag", series_col)]

    n_r1 <- 0L; n_r2 <- 0L; n_pin_override <- 0L

    for (g in seq_len(nrow(groups))) {
        ix <- groups$.idxs[[g]]
        n_g <- length(ix)
        if (n_g <= 1L) next

        ## ---- R1: split-on-break ----
        for (a in seq_len(n_g)) {
            i <- ix[a]
            if (is.na(sts[i]) || sts[i] != "broken below" || is.na(dbh[i])) next
            if (a == 1L) next
            prior <- rid[ix[seq_len(a - 1L)]]
            prior <- prior[!is.na(prior)]
            if (length(prior) == 0L || is.na(rid[i]) || !(rid[i] %in% prior)) next
            # Pin override: the BCI driver pre-stamps TrueStemID = OriginalStemID
            # on BB+DBH rows under the assumption that they start a new
            # trajectory.  When the contract is violated despite the pin, the
            # pin was applied on a wrong assumption and we override it.  We
            # also update TrueStemID to keep the two columns consistent.
            if (!is.na(trueid[i])) n_pin_override <- n_pin_override + 1L
            old_id <- rid[i]
            new_id <- next_id; next_id <- next_id + 1L
            for (b in a:n_g) {
                j <- ix[b]
                if (is.na(rid[j]) || rid[j] != old_id) break
                if (b > a && !is.na(sts[j]) && sts[j] == "broken below" && !is.na(dbh[j])) break
                rid[j] <- new_id
                if (!is.na(trueid[j]) && trueid[j] == old_id) trueid[j] <- new_id
                mth[j] <- if (b == a) "bb_split" else "bb_split_carry"
            }
            n_r1 <- n_r1 + 1L
        }

        ## ---- R2: terminate-on-stump ----
        for (a in seq_len(n_g)) {
            i <- ix[a]
            if (is.na(sts[i]) || sts[i] != "broken below" || !is.na(dbh[i])) next
            term_id <- rid[i]
            if (is.na(term_id)) next
            if (a + 1L > n_g) next
            for (b in (a + 1L):n_g) {
                j <- ix[b]
                if (is.na(rid[j]) || rid[j] != term_id || is.na(dbh[j])) next
                # alive (non-NA DBH) row reusing terminator id — split
                if (!is.na(trueid[j])) n_pin_override <- n_pin_override + 1L
                old2 <- rid[j]
                new2 <- next_id; next_id <- next_id + 1L
                for (k_ in b:n_g) {
                    k <- ix[k_]
                    if (is.na(rid[k]) || rid[k] != old2) break
                    if (k_ > b && !is.na(sts[k]) && sts[k] == "broken below" && !is.na(dbh[k])) break
                    rid[k] <- new2
                    if (!is.na(trueid[k]) && trueid[k] == old2) trueid[k] <- new2
                    mth[k] <- if (k_ == b) "bb_post_terminator_split" else "bb_post_terminator_split_carry"
                }
                n_r2 <- n_r2 + 1L
                break
            }
        }
    }

    out[, ReconstructedStemID := rid]
    out[, ReconstructionMethod := mth]
    if (has_true) out[, TrueStemID := trueid]
    if (isTRUE(verbose)) {
        message(sprintf(
            "[apply_broken_below_invariants] R1 splits: %d, R2 splits: %d, pin overrides: %d.",
            n_r1, n_r2, n_pin_override
        ))
    }
    out
}

# -------------------------------------------------------------------------
# apply_bb_invariants_to_samples()
#
# Per-sample relabel of broken-below contract violations on a posterior
# `samples_dt` produced by `dp_global_dp.R` or `dp_probabilistic_matching.R`.
# Strategy A: apply the same deterministic R1/R2 operator implemented in
# `apply_broken_below_invariants()` to each posterior sample independently.
#
# Inputs:
#   samples_dt : data.table with columns at least Sample, CensusID,
#                ReconstructedStemID, ObsRowID. Rows ordered by Sample, CensusID.
#   tree_data  : the engine's per-row data.table; must contain `obs_row_id`,
#                `Status`, `DBH`, and a series column (StemTag /
#                OriginalStemID / StemID). `Tag` is optional but used if present.
#   verbose    : log a single summary message.
#
# Returns: samples_dt with `ReconstructedStemID` rewritten so each sample
# satisfies R1 (BB+DBH must not reuse a prior live ID in the same series)
# and R2 (BB+NA-DBH terminator must not be followed by a live row reusing
# its ID). The relabel is independent per sample, so freshly-minted IDs do
# not need to be globally unique across samples — `path_sig` will continue
# to collapse identical reconstructions correctly because `paste0` over the
# rewritten ids is deterministic given the (sample, ordering) inputs.
#
# Use sites: posterior writers in `dp_global_dp.R` and
# `dp_probabilistic_matching.R`, immediately AFTER samples are sorted and
# BEFORE `path_sig` / `path_count` aggregation.
apply_bb_invariants_to_samples <- function(samples_dt, tree_data, verbose = TRUE) {
    if (is.null(samples_dt) || nrow(samples_dt) == 0L) return(samples_dt)
    if (!all(c("Sample", "CensusID", "ReconstructedStemID", "ObsRowID") %in% names(samples_dt))) {
        return(samples_dt)
    }
    if (is.null(tree_data) || !("obs_row_id" %in% names(tree_data))) {
        return(samples_dt)
    }
    series_col <- if ("StemTag" %in% names(tree_data)) {
        "StemTag"
    } else if ("OriginalStemID" %in% names(tree_data)) {
        "OriginalStemID"
    } else if ("StemID" %in% names(tree_data)) {
        "StemID"
    } else {
        return(samples_dt)
    }
    if (!all(c("Status", "DBH") %in% names(tree_data))) return(samples_dt)

    meta <- as.data.table(tree_data)[, c("obs_row_id", "Status", "DBH", series_col), with = FALSE]
    data.table::setnames(meta, "obs_row_id", "ObsRowID")

    tag_local <- if ("Tag" %in% names(samples_dt)) samples_dt$Tag[1L] else NA_integer_
    samples <- sort(unique(samples_dt$Sample))
    new_rid_all <- integer(nrow(samples_dt))
    n_changed_samples <- 0L

    for (s in samples) {
        sel <- which(samples_dt$Sample == s)
        sdt <- samples_dt[sel, .(ObsRowID, CensusID, ReconstructedStemID)]
        sdt[, .ord := seq_len(.N)]
        sdt <- merge(sdt, meta, by = "ObsRowID", all.x = TRUE, sort = FALSE)
        if (!("Tag" %in% names(sdt))) sdt[, Tag := tag_local]
        before <- sdt$ReconstructedStemID
        sdt2 <- apply_broken_below_invariants(sdt, verbose = FALSE)
        data.table::setorder(sdt2, .ord)
        if (!identical(before, sdt2$ReconstructedStemID)) n_changed_samples <- n_changed_samples + 1L
        new_rid_all[sel] <- sdt2$ReconstructedStemID
    }
    samples_dt[, ReconstructedStemID := new_rid_all]
    if (isTRUE(verbose)) {
        message(sprintf(
            "[apply_bb_invariants_to_samples] Relabeled %d / %d posterior sample(s) for BB invariants.",
            n_changed_samples, length(samples)
        ))
    }
    samples_dt
}

# End of dp_global_main.R
# -------------------------------------------------------------------------
