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

# -------------------------------------------------------------------------
# renumber_engine_minted_ids()
#
# Post-engine renumbering pass that reverses the direction of
# ReconstructedStemID numbering so smaller integers correspond to earlier
# stem appearances (matching the OriginalStemID convention used in BCI /
# ForestGEO databases).
#
# Per-tag algorithm (see dp_global/improvements.md for full design):
#   1. Identify known_ids = real BCI database IDs that must never be renamed:
#       - db_ids: unique non-NA TrueStemID on rows whose ReconstructionMethod
#         is NOT one of ENGINE_MINTED_INTO_TRUESTEMID (provisional_dp,
#         bb_split, bb_split_carry, bb_post_terminator_split,
#         bb_post_terminator_split_carry). Those rows hold engine-minted
#         values in TrueStemID (DP minting site 1; bb pin-overrides in
#         apply_broken_below_invariants).
#       - orphan_ids: ReconstructedStemID on rows tagged "given_orphan"
#         (assigned from StemID by apply_orphan_stem_backfill; real DB IDs
#         even though TrueStemID is NA there).
#   2. engine_ids = setdiff(all_recon, known_ids)  -- black-box partition
#   3. Sort engine_ids by first_census ascending (tie-break by original
#      ReconstructedStemID value ascending for determinism).
#   4. Mint new IDs descending from base = min(known_ids) - length(engine_ids)
#      so engine_ids_sorted[1] (earliest) gets base, and all engine IDs are
#      strictly < min(known_ids).
#   5. Apply mapping to: ReconstructedStemID, ReconstructedStemID_PreSweep,
#      DP_PosteriorTop{k}ID columns (always), and TrueStemID ONLY on rows
#      where ReconstructionMethod %in% ENGINE_MINTED_INTO_TRUESTEMID.
#   6. (Optional) write companion mapping file alongside posterior path
#      files so downstream consumers can translate the `recon` / `path_sig`
#      fields. See improvements.md "Fallback" subsection for the limitation
#      regarding per-sample bb-invariant IDs.
#
# Inputs:
#   out                      : data.table (per-chunk or per-tag) with at
#                              least Tag, CensusID, ReconstructedStemID,
#                              TrueStemID, ReconstructionMethod.
#   posterior_top_k          : (optional) explicit number of
#                              DP_PosteriorTop{k}ID columns. If NULL,
#                              auto-detected from column names.
#   posterior_samples_path   : (optional) directory containing the
#                              `posteriors/` subdirectory. If non-NULL,
#                              writes per-tag companion mapping files there.
#   mapping_format           : "feather" | "rds" | "csv". Default feather
#                              with rds fallback if arrow is unavailable.
#   verbose                  : log a per-call summary.
#
# Returns: list(out = renamed data.table, mapping = combined mapping
# data.table with columns Tag, old_id, new_id, first_census).
renumber_engine_minted_ids <- function(out,
                                       posterior_top_k = NULL,
                                       posterior_samples_path = NULL,
                                       mapping_format = c("feather", "rds", "csv"),
                                       verbose = TRUE) {
    if (is.null(out) || nrow(out) == 0L) {
        return(list(out = out, mapping = data.table::data.table(
            Tag = integer(0), old_id = integer(0), new_id = integer(0),
            first_census = integer(0)
        )))
    }
    needed <- c("Tag", "CensusID", "ReconstructedStemID", "ReconstructionMethod")
    if (!all(needed %in% names(out))) {
        if (isTRUE(verbose)) {
            message("[renumber_engine_minted_ids] Required columns missing; skipping.")
        }
        return(list(out = out, mapping = data.table::data.table(
            Tag = integer(0), old_id = integer(0), new_id = integer(0),
            first_census = integer(0)
        )))
    }
    mapping_format <- match.arg(mapping_format)

    # Methods that write engine-minted IDs into TrueStemID (see
    # improvements.md, Pass 7 fix).
    ENGINE_MINTED_INTO_TRUESTEMID <- c(
        "provisional_dp",
        "bb_split", "bb_split_carry",
        "bb_post_terminator_split", "bb_post_terminator_split_carry"
    )

    # Auto-detect DP_PosteriorTop{k}ID columns if not specified
    posterior_id_cols <- grep("^DP_PosteriorTop\\d+ID$", names(out), value = TRUE)
    if (is.null(posterior_top_k)) {
        posterior_top_k <- length(posterior_id_cols)
    } else {
        posterior_top_k <- as.integer(posterior_top_k)
        wanted <- paste0("DP_PosteriorTop", seq_len(posterior_top_k), "ID")
        posterior_id_cols <- intersect(wanted, names(out))
    }

    has_pre_sweep <- "ReconstructedStemID_PreSweep" %in% names(out)
    has_true <- "TrueStemID" %in% names(out)

    # Local copies (single-pass updates avoid repeated data.table writes)
    rid    <- as.integer(out$ReconstructedStemID)
    trueid <- if (has_true) as.integer(out$TrueStemID) else NULL
    presweep <- if (has_pre_sweep) as.integer(out$ReconstructedStemID_PreSweep) else NULL
    pst_cols <- lapply(posterior_id_cols, function(cn) as.integer(out[[cn]]))
    names(pst_cols) <- posterior_id_cols

    tag_vec <- out$Tag
    census_vec <- as.integer(out$CensusID)
    method_vec <- as.character(out$ReconstructionMethod)

    # Group row indices by Tag (preserves order)
    tag_groups <- out[, .(.idxs = list(.I)), by = Tag]

    mapping_list <- vector("list", nrow(tag_groups))
    n_tags_renumbered <- 0L
    n_tags_no_engine <- 0L
    n_tags_no_known <- 0L
    n_engine_ids_total <- 0L

    for (g in seq_len(nrow(tag_groups))) {
        ix <- tag_groups$.idxs[[g]]
        if (length(ix) == 0L) next
        tag_id <- tag_groups$Tag[g]

        rid_g <- rid[ix]
        method_g <- method_vec[ix]
        census_g <- census_vec[ix]
        trueid_g <- if (has_true) trueid[ix] else rep(NA_integer_, length(ix))

        # ---- Step 1a: db_ids (real DB IDs in TrueStemID) ----
        excl_db <- !is.na(method_g) & (method_g %in% ENGINE_MINTED_INTO_TRUESTEMID)
        db_mask <- !is.na(trueid_g) & !excl_db
        db_ids <- unique(trueid_g[db_mask])

        # ---- Step 1b: orphan_ids (real DB IDs from StemID via orphan backfill) ----
        orph_mask <- !is.na(method_g) & method_g == "given_orphan" & !is.na(rid_g)
        orphan_ids <- unique(rid_g[orph_mask])

        known_ids <- unique(c(db_ids, orphan_ids))

        # ---- Step 2-3: partition ----
        all_recon <- unique(rid_g[!is.na(rid_g)])
        engine_ids <- setdiff(all_recon, known_ids)

        if (length(engine_ids) == 0L) {
            n_tags_no_engine <- n_tags_no_engine + 1L
            next
        }

        # ---- Step 4: first_census per engine_id (within this tag) ----
        # Build a small data.table for efficient grouping
        eng_dt <- data.table::data.table(
            id = rid_g,
            census = census_g
        )[!is.na(id) & id %in% engine_ids,
          .(first_census = min(census, na.rm = TRUE)),
          by = id]

        # ---- Step 5: sort by first_census ascending; tie-break by id ascending ----
        data.table::setorder(eng_dt, first_census, id)
        n_eng <- nrow(eng_dt)

        # ---- Step 6: compute new IDs ----
        if (length(known_ids) == 0L) {
            # Edge case: no real DB IDs anchoring the renumber space; use 0L
            # as base so engine IDs become non-positive sequential integers.
            # Chronological ordering within the tag is still preserved.
            n_tags_no_known <- n_tags_no_known + 1L
            base_id <- 0L - n_eng
            if (isTRUE(verbose)) {
                message(sprintf(
                    "[renumber_engine_minted_ids] Tag %s: known_ids is empty; using base=%d.",
                    as.character(tag_id), base_id
                ))
            }
        } else {
            base_id <- as.integer(min(known_ids)) - n_eng
        }
        eng_dt[, new_id := base_id + seq_len(.N) - 1L]

        # ---- Step 7: build mapping table ----
        map_lookup <- eng_dt$new_id
        names(map_lookup) <- as.character(eng_dt$id)

        translate <- function(v) {
            if (is.null(v)) return(NULL)
            out_v <- v
            hit <- !is.na(v) & as.character(v) %in% names(map_lookup)
            if (any(hit)) {
                out_v[hit] <- map_lookup[as.character(v[hit])]
            }
            out_v
        }

        # ---- Step 8: apply mapping to columns ----
        rid[ix] <- translate(rid_g)
        if (has_pre_sweep) {
            presweep[ix] <- translate(presweep[ix])
        }
        for (cn in posterior_id_cols) {
            pst_cols[[cn]][ix] <- translate(pst_cols[[cn]][ix])
        }
        # TrueStemID: only on rows whose method indicates engine-minted-into-TrueStemID
        if (has_true) {
            engine_minted_rows <- !is.na(method_g) & (method_g %in% ENGINE_MINTED_INTO_TRUESTEMID)
            if (any(engine_minted_rows)) {
                trueid_sub <- trueid_g
                trueid_sub[engine_minted_rows] <- translate(trueid_sub[engine_minted_rows])
                # Only commit changes on engine_minted_rows (leave other rows alone)
                trueid[ix[engine_minted_rows]] <- trueid_sub[engine_minted_rows]
            }
        }

        mapping_list[[g]] <- data.table::data.table(
            Tag = rep(tag_id, n_eng),
            old_id = as.integer(eng_dt$id),
            new_id = as.integer(eng_dt$new_id),
            first_census = as.integer(eng_dt$first_census)
        )
        n_tags_renumbered <- n_tags_renumbered + 1L
        n_engine_ids_total <- n_engine_ids_total + n_eng
    }

    # Commit column changes
    out[, ReconstructedStemID := rid]
    if (has_pre_sweep) out[, ReconstructedStemID_PreSweep := presweep]
    if (has_true) out[, TrueStemID := trueid]
    for (cn in posterior_id_cols) {
        out[, (cn) := pst_cols[[cn]]]
    }

    mapping <- if (n_tags_renumbered > 0L) {
        data.table::rbindlist(Filter(Negate(is.null), mapping_list),
                              use.names = TRUE, fill = TRUE)
    } else {
        data.table::data.table(Tag = integer(0), old_id = integer(0),
                               new_id = integer(0), first_census = integer(0))
    }

    # ---- Step 9: (removed) companion mapping files are no longer written.
    # The recommended architecture rewrites the full paths file in the
    # renumbered ID space via finalize_posterior_paths(); see
    # dp_global/improvements.md. The `posterior_samples_path` and
    # `mapping_format` arguments are retained for backward compatibility
    # but are now no-ops.
    if (!is.null(posterior_samples_path) && isTRUE(verbose)) {
        # silent no-op; argument kept for backward compatibility
    }

    if (isTRUE(verbose)) {
        message(sprintf(
            "[renumber_engine_minted_ids] Renumbered %d tag(s), %d engine ID(s); %d tag(s) had no engine IDs%s.",
            n_tags_renumbered, n_engine_ids_total, n_tags_no_engine,
            if (n_tags_no_known > 0L) sprintf(", %d tag(s) had empty known_ids", n_tags_no_known) else ""
        ))
    }

    list(out = out, mapping = mapping)
}

# -------------------------------------------------------------------------
# finalize_posterior_paths()
#
# Recommended-architecture post-engine step (see dp_global/improvements.md).
# Reads per-tag raw posterior staging files written by the engines (DP and
# probabilistic), translates ReconstructedStemID via the renumber mapping,
# re-runs apply_bb_invariants_to_samples() so per-sample bb-minted IDs are
# derived from the renumbered track IDs, computes path signatures and
# probabilities, and writes the final `tag_{Tag}_posterior_samples_{ts}_paths.{ext}`
# files in the user-requested format. Staging files are deleted on success.
#
# Inputs:
#   out                    : data.table AFTER renumber_engine_minted_ids().
#                            Must contain at least: Tag, CensusID,
#                            ReconstructedStemID, Status, DBH, obs_row_id,
#                            and one of (StemTag | OriginalStemID | StemID).
#   posterior_samples_path : directory whose `posteriors/.staging/` subdir
#                            holds raw staging files produced by the engines.
#                            Must be the SAME path passed to the engines.
#   mapping                : data.table with cols (Tag, old_id, new_id) as
#                            returned by renumber_engine_minted_ids(). If
#                            NULL or empty, ReconstructedStemID values are
#                            assumed to already be in renumbered space (the
#                            no-op identity translation).
#   verbose                : log a per-tag summary message.
#
# Returns: invisible list(n_tags, n_written, n_failed, written_paths).
finalize_posterior_paths <- function(out,
                                     posterior_samples_path,
                                     mapping = NULL,
                                     verbose = TRUE) {
    null_result <- function() invisible(list(
        n_tags = 0L, n_written = 0L, n_failed = 0L, written_paths = character(0)
    ))
    if (is.null(posterior_samples_path) || !nzchar(posterior_samples_path)) {
        return(null_result())
    }
    staging_dir <- file.path(posterior_samples_path, "posteriors", ".staging")
    if (!dir.exists(staging_dir)) return(null_result())
    staging_files <- list.files(staging_dir,
                                pattern = "^tag_.*_samples_raw_.*\\.rds$",
                                full.names = TRUE)
    if (length(staging_files) == 0L) return(null_result())
    if (is.null(out) || nrow(out) == 0L) {
        warning("[finalize_posterior_paths] `out` is empty; cannot finalize ",
                length(staging_files), " staging file(s).")
        return(null_result())
    }

    # Build per-tag mapping lookup (tag → named vector old_id_chr → new_id)
    map_by_tag <- list()
    if (!is.null(mapping) && nrow(mapping) > 0L) {
        for (tg in unique(mapping$Tag)) {
            sub <- mapping[Tag == tg]
            v <- as.integer(sub$new_id)
            names(v) <- as.character(sub$old_id)
            map_by_tag[[as.character(tg)]] <- v
        }
    }

    out_dt <- if (data.table::is.data.table(out)) out else data.table::as.data.table(out)
    n_written <- 0L
    n_failed <- 0L
    written_paths <- character(0)

    for (sf in staging_files) {
        info <- tryCatch(readRDS(sf), error = function(e) NULL)
        if (is.null(info) || is.null(info$samples_dt) || nrow(info$samples_dt) == 0L) {
            n_failed <- n_failed + 1L
            next
        }
        tag_val <- info$tag_val
        samples_dt <- data.table::copy(info$samples_dt)
        fmt <- info$posterior_samples_format
        ts_local <- info$batch_ts

        # Subset out to this tag for bb-invariant metadata
        tree_data_for_tag <- if (is.na(tag_val)) {
            out_dt[is.na(Tag)]
        } else {
            out_dt[Tag == tag_val]
        }
        if (nrow(tree_data_for_tag) == 0L) {
            warning(sprintf(
                "[finalize_posterior_paths] No rows in `out` for tag %s; skipping %s.",
                as.character(tag_val), basename(sf)
            ))
            n_failed <- n_failed + 1L
            next
        }

        # ---- (1) Translate ReconstructedStemID via mapping ----
        tag_key <- as.character(tag_val)
        if (!is.null(map_by_tag[[tag_key]])) {
            lk <- map_by_tag[[tag_key]]
            v <- samples_dt$ReconstructedStemID
            hit <- !is.na(v) & as.character(v) %in% names(lk)
            if (any(hit)) {
                v[hit] <- lk[as.character(v[hit])]
                samples_dt[, ReconstructedStemID := as.integer(v)]
            }
        }

        # ---- (2) Re-run bb invariants in renumbered ID space ----
        samples_dt <- apply_bb_invariants_to_samples(
            samples_dt, tree_data_for_tag, verbose = FALSE
        )

        # ---- (3) path_sig + paths_summary ----
        sample_sigs <- samples_dt[, .(path_sig = paste0(ReconstructedStemID, collapse = "-")),
                                  by = Sample]
        path_counts <- sample_sigs[, .N, by = path_sig]
        data.table::setnames(path_counts, "N", "path_count")
        n_samp <- data.table::uniqueN(sample_sigs$Sample)

        has_logp <- "logp" %in% names(samples_dt)
        if (has_logp) {
            sample_logp <- unique(samples_dt[, .(Sample, logp)])
            maxlp <- max(sample_logp$logp, na.rm = TRUE)
            sample_logp[, sample_weight := exp(logp - maxlp)]
            sample_logp[, sample_prob := sample_weight / sum(sample_weight)]
            sample_sigs <- merge(sample_sigs, sample_logp[, .(Sample, sample_prob)],
                                 by = "Sample", all.x = TRUE)
            path_probs <- sample_sigs[, .(path_prob = sum(sample_prob, na.rm = TRUE)),
                                      by = path_sig]
            paths_summary <- merge(path_counts, path_probs, by = "path_sig", all.x = TRUE)
        } else {
            paths_summary <- path_counts[, .(path_count,
                                             path_prob = path_count / n_samp),
                                         by = path_sig]
        }

        # Compact reconstruction mapping (one row per path_sig)
        sigs_by_sample <- sample_sigs[, .(Sample, path_sig)]
        samples_with_sig <- merge(samples_dt[, .(Sample, ObsRowID, ReconstructedStemID)],
                                  sigs_by_sample, by = "Sample", all.x = TRUE)
        recon_by_path <- samples_with_sig[, .(recon = paste0(ObsRowID, ":",
                                                             ReconstructedStemID,
                                                             collapse = ";")),
                                          by = .(path_sig, Sample)]
        recon_compact <- recon_by_path[, .SD[1], by = path_sig, .SDcols = "recon"]
        paths_summary <- merge(paths_summary, recon_compact, by = "path_sig",
                               all.x = TRUE)

        # ---- (4) Write final paths file ----
        out_dir_post <- file.path(info$posterior_samples_path, "posteriors")
        if (!dir.exists(out_dir_post)) {
            dir.create(out_dir_post, recursive = TRUE, showWarnings = FALSE)
        }
        out_path_base <- file.path(out_dir_post, paste0(
            "tag_", ifelse(is.na(tag_val), "NA", tag_val),
            "_posterior_samples_", ts_local
        ))
        wrote <- tryCatch({
            if (fmt == "feather" && requireNamespace("arrow", quietly = TRUE)) {
                p <- paste0(out_path_base, "_paths.feather")
                arrow::write_feather(paths_summary, p)
            } else if (fmt == "csv") {
                p <- paste0(out_path_base, "_paths.csv")
                data.table::fwrite(paths_summary, p)
            } else {
                p <- paste0(out_path_base, "_paths.rds")
                saveRDS(paths_summary, file = p)
            }
            p
        }, error = function(e) {
            warning(sprintf(
                "[finalize_posterior_paths] Tag %s: write failed: %s",
                as.character(tag_val), e$message
            ))
            NA_character_
        })
        if (!is.na(wrote)) {
            n_written <- n_written + 1L
            written_paths <- c(written_paths, wrote)
            file.remove(sf)
            if (isTRUE(verbose)) {
                message(sprintf(
                    "[finalize_posterior_paths] Tag %s: wrote %d unique path(s) to %s",
                    as.character(tag_val), nrow(paths_summary), wrote
                ))
            }
        } else {
            n_failed <- n_failed + 1L
        }
    }

    if (isTRUE(verbose)) {
        message(sprintf(
            "[finalize_posterior_paths] Finalized %d tag(s); %d written, %d failed.",
            length(staging_files), n_written, n_failed
        ))
    }
    invisible(list(
        n_tags = length(staging_files),
        n_written = n_written,
        n_failed = n_failed,
        written_paths = written_paths
    ))
}

# End of dp_global_main.R
# -------------------------------------------------------------------------
