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
tryCatch(
    {
        Rcpp::sourceCpp(file.path(root_dir, "dp_global", "src", "transition_cost_rcpp.cpp"))
        message("[dp_global_main.R] C++ acceleration enabled.")
    },
    error = function(e) {
        warning(sprintf("C++ compilation (transition_cost_rcpp.cpp) failed or is unavailable: %s. Continuing without compiled acceleration.", e$message))
    }
)

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
#   - Status is one of "dead", "stem dead", "broken below", "missing"
#
# For each such row we copy the most recent prior non-NA ReconstructedStemID
# from the same (Tag, OriginalStemID) group (LOCF). Biologically, a terminal
# event ends the trajectory of the most recent prior identity carrying the
# same OriginalStemID; without this fill these rows would be dropped from
# any downstream trajectory. A "missing" record after the stem's earlier
# records is likewise that stem's record (99.9% of BCI "missing" rows follow
# a record of the same StemID); without the fill it would become a stem of its
# own that never lives. Stage 3 reads "missing" as no data, so no status
# changes. Rows with no earlier record of their StemID are left NA (LOCF only
# looks backward).
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
    needed <- c(
        "Status", "DBH", src_col, "Tag",
        "CensusID", "ReconstructedStemID"
    )
    if (!all(needed %in% names(out))) {
        return(out)
    }
    data.table::setorderv(out, c("Tag", src_col, "CensusID"))
    .term_mask <- is.na(out$ReconstructedStemID) &
        is.na(out$DBH) &
        !is.na(out$Status) &
        out$Status %in% c("dead", "stem dead", "broken below", "missing")
    if (!any(.term_mask)) {
        return(out)
    }
    out[, .carried := data.table::nafill(ReconstructedStemID, type = "locf"),
        by = c("Tag", src_col)
    ]
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

## ---- 8a) Shared post-engine helper: rows left behind by the pin sweep ----
# The TrueStemID sweep sets ReconstructedStemID := TrueStemID on pinned rows
# only. When the engine linked pinned rows to unpinned measurements of the
# same stem, those unpinned rows keep another ID and the stem is cut in two
# (a death plus a recruit that the engine never made). Typical sources:
#   - Step 3a.5 pins the alive rows of a StemID that has a death record; the
#     stem's earlier rows sit under an older StemID (the database renumbered
#     21-28% of stems at every census before 2010) and have no pin (tree
#     100015: StemID 100015 in 1982 -> 638594 in 1985, dead 1990);
#   - the DP resprout segment split offsets the pre-segment IDs, including
#     tracks it had labelled with a pin, so after the sweep the pinned rows
#     hold the pin and their unpinned track-mates the offset ID.
# The posterior samples keep the engine's links, so these cuts inflated the
# exported BA Loss and Gain and put the exported line outside the identity
# Monte Carlo ribbon.
#
# Rule, per engine track (measured rows of one Tag and
# ReconstructedStemID_PreSweep), cut into segments at every broken-below row
# with a DBH (R1: a new stem starts there) and after every NA-R barrier census
# (no measured stem but a stump / R-coded record: the whole tree died back):
#   - the segment's pinned rows all hold their pin and share ONE value P;
#   - an unpinned measured row of the segment holding another ID T joins P,
#     with every row of T (measured or not), when all measured rows of T are
#     in this segment, no row of T carries a pin, and P has no row in any
#     census of T (no collision).
# Keys that hold two measurements in one census are not engine tracks and
# are skipped. Pinned rows never move, so ReconstructedStemID == TrueStemID
# still holds; tracks with two pin values or a pin that is not held (Fix-2
# rollback) are left as they are; the links restored are the engine's own.
# Moved rows get ReconstructionMethod = "pin_track"; TrueStemID is unchanged.
# Runs first in the post-engine chain, after the sweeps, so the backfills,
# apply_terminal_to_host() and the broken-below invariants see whole stems.
# Idempotent: a moved ID T no longer exists.
apply_pin_track_rejoin <- function(out, verbose = TRUE) {
    if (is.null(out) || nrow(out) == 0L) {
        return(out)
    }
    needed <- c("Tag", "CensusID", "DBH", "Status", "TrueStemID", "ReconstructedStemID", "ReconstructedStemID_PreSweep")
    if (!all(needed %in% names(out))) {
        return(out)
    }
    tg <- as.character(out$Tag)
    cen <- as.integer(out$CensusID)
    rid <- as.integer(out$ReconstructedStemID)
    pre <- as.integer(out$ReconstructedStemID_PreSweep)
    pin <- as.integer(out$TrueStemID)
    sts <- as.character(out$Status)
    tsm <- if ("ListOfTSM" %in% names(out)) as.character(out$ListOfTSM) else rep(NA_character_, nrow(out))
    meas <- !is.na(out$DBH)

    # NA-R barrier censuses: no measured stem, but a stump / R-coded record.
    r_na <- !meas & ((!is.na(sts) & sts == "broken below") |
        (!is.na(tsm) & grepl("\\b(R|RP|RF|RT|QR|OR)\\b", tsm, perl = TRUE)))
    bar <- data.table::data.table(tag = tg, cen = cen, meas = meas, r_na = r_na)[,
        .(n_meas = sum(meas), n_rna = sum(r_na)),
        by = .(tag, cen)
    ][n_meas == 0L & n_rna > 0L]

    # Engine tracks and their segments.
    m <- data.table::data.table(.row = which(meas & !is.na(pre) & !is.na(rid)))
    if (nrow(m) == 0L) {
        return(out)
    }
    m[, `:=`(
        tag = tg[.row], cen = cen[.row], rid = rid[.row], pre = pre[.row], pin = pin[.row],
        bb = !is.na(sts[.row]) & sts[.row] == "broken below"
    )]
    m[, nb := 0L]
    if (nrow(bar) > 0L) {
        bl <- split(bar$cen, bar$tag)
        m[tag %in% names(bl), nb := findInterval(cen - 0.5, sort(bl[[tag[1L]]])), by = tag]
    }
    data.table::setorder(m, tag, pre, cen)
    m[, seg := cumsum(bb) + nb, by = .(tag, pre)]
    mixed <- m[, .N, by = .(tag, pre, cen)][N > 1L, unique(paste(tag, pre))]
    m <- m[!paste(tag, pre) %in% mixed]

    # Segments whose pinned rows share one held pin value P.
    sg <- m[, {
        p <- pin[!is.na(pin)]
        ok <- length(p) > 0L && all(rid[!is.na(pin)] == p) && data.table::uniqueN(p) == 1L
        .(P = if (ok) p[1L] else NA_integer_)
    }, by = .(tag, pre, seg)][!is.na(P)]
    cand <- m[sg, on = .(tag, pre, seg), nomatch = 0L][is.na(pin) & rid != P]
    n_moved_ids <- 0L
    n_moved_rows <- 0L
    n_collide <- 0L
    if (nrow(cand) > 0L) {
        # An ID T moves only if all its measured rows are candidates of one segment
        # and none of its rows carries a pin.
        all_rows <- data.table::data.table(.row = seq_along(rid), tag = tg, rid = rid, cen = cen, pin = pin, meas = meas)[!is.na(rid)]
        tstat <- all_rows[, .(n_meas = sum(meas), n_pin = sum(!is.na(pin))), by = .(tag, rid)]
        ct <- cand[, .(n_c = .N, n_seg = data.table::uniqueN(paste(pre, seg)), P = P[1L], n_P = data.table::uniqueN(P)), by = .(tag, rid)]
        ct <- tstat[ct, on = .(tag, rid)][n_c == n_meas & n_seg == 1L & n_P == 1L & n_pin == 0L]
        if (nrow(ct) > 0L) {
            rows_T <- all_rows[ct[, .(tag, rid, P)], on = .(tag, rid), nomatch = 0L]
            # collision with rows that already hold P, or with another ID joining P
            occ <- unique(all_rows[, .(tag, P = rid, cen)])
            rows_T[, hit := FALSE]
            rows_T[occ, on = .(tag, P, cen), hit := TRUE]
            rows_T[, dup := .N > 1L, by = .(tag, P, cen)]
            bad_T <- unique(rows_T[hit | dup, .(tag, rid)])
            n_collide <- nrow(bad_T)
            rows_T <- rows_T[!bad_T, on = .(tag, rid)]
            if (nrow(rows_T) > 0L) {
                val <- if (is.integer(out$ReconstructedStemID)) as.integer(rows_T$P) else as.numeric(rows_T$P)
                data.table::set(out, i = rows_T$.row, j = "ReconstructedStemID", value = val)
                if (!("ReconstructionMethod" %in% names(out))) out[, ReconstructionMethod := NA_character_]
                data.table::set(out, i = rows_T$.row, j = "ReconstructionMethod", value = "pin_track")
                n_moved_ids <- data.table::uniqueN(rows_T[, .(tag, rid)])
                n_moved_rows <- nrow(rows_T)
            }
        }
    }
    if (isTRUE(verbose)) {
        message(sprintf(
            "[apply_pin_track_rejoin] %d left-behind ID(s) rejoined their track's pin (%d row(s)); %d skipped for a collision.",
            n_moved_ids, n_moved_rows, n_collide
        ))
    }
    out
}

## ---- 8b) Shared post-engine helper: terminal records back to their stem --
# An unmeasured dead / stem dead / broken-below record that comes BEFORE the
# first evidence of life of the identity it carries (or sits on an identity
# that never lives) is not the start of a stem: it is the death / break record
# of a stem that ended just before. The BCI database usually stores such a
# record under a NEW StemID (the StemID of the resprout measured later), so
# the pre-DP pins (Step 1(a) StemTag pins, Step 3b propagation) and the
# TrueStemID sweep put it on the resprout's identity, where R2 then cuts it
# off into a stem that never lives. The engine itself had usually placed it
# on the stem that broke (ReconstructedStemID_PreSweep).
#
# Rule, per group of such records (same Tag, identity and database StemID):
#   (a) the engine's choice: the group's ReconstructedStemID_PreSweep, if that
#       identity has evidence of life (DBH or raw "alive") and its last
#       evidence of life is before the group's first record. If that engine
#       ID no longer lives (its track was relabelled, e.g. by
#       apply_pin_track_rejoin()), the current ID of the engine track's last
#       measured row before the record is used instead;
#   (b) otherwise the one stem that ended just before: the only identity of
#       the tag whose last evidence of life is in census c1 - 1 (c1 = the
#       group's first record), when exactly one unresolved group starts at c1;
#   (c) otherwise the records stay where they are (a never-alive stem).
# A target must have no row in any census of the group (no collision).
# Only unmeasured rows move, so measured trajectories, DBH and statuses are
# unchanged, and each record keeps its own Status (DFstatus downstream).
# Moved rows get ReconstructionMethod = "terminal_to_host"; under (a)
# SweepRollbackToPreSweep is set to TRUE because the engine's own choice is
# restored. TrueStemID is not changed (the sweep is not re-run after this
# pass, as for the broken-below invariants).
# Runs after apply_orphan_stem_backfill() and before
# apply_broken_below_invariants(). Idempotent: a moved record sits after its
# new stem's life, so it is no longer a candidate.
apply_terminal_to_host <- function(out, verbose = TRUE) {
    if (is.null(out) || nrow(out) == 0L) {
        return(out)
    }
    needed <- c("Tag", "CensusID", "Status", "DBH", "ReconstructedStemID")
    if (!all(needed %in% names(out))) {
        return(out)
    }
    src_col <- if ("StemID" %in% names(out)) {
        "StemID"
    } else if ("OriginalStemID" %in% names(out)) {
        "OriginalStemID"
    } else {
        NULL
    }
    has_pre <- "ReconstructedStemID_PreSweep" %in% names(out)
    rid <- as.integer(out$ReconstructedStemID)
    cen <- as.integer(out$CensusID)
    tg <- as.character(out$Tag)
    sts <- as.character(out$Status)
    alive_ev <- !is.na(out$DBH) | (!is.na(sts) & sts == "alive")
    is_term <- is.na(out$DBH) & !is.na(rid) & !is.na(sts) & sts %in% c("dead", "stem dead", "broken below")
    if (!any(is_term)) {
        return(out)
    }

    # Evidence-of-life span of every identity that lives.
    live <- alive_ev & !is.na(rid)
    span <- if (any(live)) {
        data.table::data.table(tag = tg[live], host = rid[live], cen = cen[live])[,
            .(h_first = min(cen), h_last = max(cen)),
            by = .(tag, host)
        ]
    } else {
        data.table::data.table(tag = character(0), host = integer(0), h_first = integer(0), h_last = integer(0))
    }
    # Candidate records: before their identity's first life, or on an identity that never lives.
    cand <- data.table::data.table(
        .row = which(is_term), tag = tg[is_term], host = rid[is_term], cen = cen[is_term],
        src = if (is.null(src_col)) NA_character_ else as.character(out[[src_col]][is_term]),
        pre = if (has_pre) as.integer(out$ReconstructedStemID_PreSweep[is_term]) else NA_integer_
    )
    cand <- span[cand, on = .(tag, host)]
    cand <- cand[is.na(h_first) | cen < h_first]
    if (nrow(cand) == 0L) {
        return(out)
    }
    data.table::setnames(cand, "host", "rid")
    cand[, grp := .GRP, by = .(tag, rid, src)]
    G <- cand[, .(
        tag = tag[1L], rid = rid[1L], c1 = min(cen), cens = list(sort(unique(cen))),
        pre = if (all(!is.na(pre)) && data.table::uniqueN(pre) == 1L) pre[1L] else NA_integer_
    ), by = grp]

    # Rows already present per (tag, identity, census), plus claims made here.
    occ <- unique(data.table::data.table(tag = tg, host = rid, cen = cen)[!is.na(host)])
    data.table::setkey(occ, tag, host, cen)
    claimed <- character(0)
    collides <- function(t, h, cs) {
        nrow(occ[data.table::CJ(tag = t, host = h, cen = cs), nomatch = 0L]) > 0L ||
            any(paste(t, h, cs, sep = "\r") %in% claimed)
    }
    target <- rep(NA_integer_, nrow(G))
    rule <- rep(NA_character_, nrow(G))

    # (a) the engine's own choice. When that engine ID no longer lives (a later
    # step relabelled the track, e.g. apply_pin_track_rejoin()), follow the
    # engine track: the current ID of its last measured row before c1.
    ga_in <- G[!is.na(pre), .(grp, tag, host = pre, c1, rid)]
    if (has_pre && nrow(ga_in) > 0L) {
        trk <- data.table::data.table(
            tag = tg, pre = as.integer(out$ReconstructedStemID_PreSweep), cen = cen, cur = rid
        )[!is.na(out$DBH) & !is.na(pre) & !is.na(cur)]
        dead_host <- ga_in[!span, on = .(tag, host), which = TRUE]
        if (length(dead_host) > 0L && nrow(trk) > 0L) {
            q <- ga_in[dead_host, .(k = dead_host, tag, pre = host, c1)]
            last <- trk[q, on = .(tag, pre, cen < c1), .(k, cur = x.cur, cen = x.cen), nomatch = 0L][order(k, -cen)][, .SD[1L], by = k]
            ga_in[last$k, host := last$cur]
        }
    }
    ga <- span[ga_in[host != rid, .(grp, tag, host, c1)], on = .(tag, host), nomatch = 0L][h_last < c1]
    for (k in seq_len(nrow(ga))) {
        gi <- match(ga$grp[k], G$grp)
        cs <- G$cens[[gi]]
        if (!collides(ga$tag[k], ga$host[k], cs)) {
            target[gi] <- ga$host[k]
            rule[gi] <- "a"
            claimed <- c(claimed, paste(ga$tag[k], ga$host[k], cs, sep = "\r"))
        }
    }

    # (b) the one stem that ended just before, when exactly one group starts then
    un <- G[is.na(target)]
    if (nrow(un) > 0L) {
        starts <- un[, .(n_groups = .N), by = .(tag, c1)]
        ends <- span[, .(tag, host, c1 = h_last + 1L)]
        cb <- ends[un[, .(grp, tag, c1, rid)], on = .(tag, c1), nomatch = 0L, allow.cartesian = TRUE][host != rid]
        if (nrow(cb) > 0L) {
            cb[, free := !vapply(seq_len(.N), function(k) collides(tag[k], host[k], G$cens[[match(grp[k], G$grp)]]), logical(1))]
            cb <- starts[cb[free == TRUE], on = .(tag, c1)]
            one <- cb[, .(n_host = .N, host = host[1L], n_groups = n_groups[1L]), by = grp][n_host == 1L & n_groups == 1L]
            for (k in seq_len(nrow(one))) {
                gi <- match(one$grp[k], G$grp)
                cs <- G$cens[[gi]]
                if (!collides(G$tag[gi], one$host[k], cs)) {
                    target[gi] <- one$host[k]
                    rule[gi] <- "b"
                    claimed <- c(claimed, paste(G$tag[gi], one$host[k], cs, sep = "\r"))
                }
            }
        }
    }

    moved <- data.table::data.table(grp = G$grp, target = target, rule = rule)[!is.na(target)]
    if (nrow(moved) > 0L) {
        rows <- cand[moved, on = "grp"]
        val <- if (is.integer(out$ReconstructedStemID)) as.integer(rows$target) else as.numeric(rows$target)
        data.table::set(out, i = rows$.row, j = "ReconstructedStemID", value = val)
        if (!("ReconstructionMethod" %in% names(out))) out[, ReconstructionMethod := NA_character_]
        data.table::set(out, i = rows$.row, j = "ReconstructionMethod", value = "terminal_to_host")
        if ("SweepRollbackToPreSweep" %in% names(out) && any(rows$rule == "a")) {
            data.table::set(out, i = rows$.row[rows$rule == "a"], j = "SweepRollbackToPreSweep", value = TRUE)
        }
    }
    if (isTRUE(verbose)) {
        message(sprintf(
            "[apply_terminal_to_host] %d group(s) / %d row(s): engine choice %d / %d, only ending stem %d / %d, left as never-alive %d / %d.",
            nrow(G), nrow(cand),
            sum(rule %in% "a"), cand[grp %in% moved[rule == "a", grp], .N],
            sum(rule %in% "b"), cand[grp %in% moved[rule == "b", grp], .N],
            sum(is.na(rule)), cand[!grp %in% moved$grp, .N]
        ))
    }
    out
}

## ---- 9) Shared post-engine helper: broken-below invariant pass --------
# Enforce two invariants on the final reconstruction, per reconstructed
# trajectory: the rows of one (Tag, ReconstructedStemID), ordered by CensusID.
#
#   R1 (split-on-break): a row with Status == "broken below" and
#       !is.na(DBH) that has an earlier row on its trajectory MUST start a
#       new identity.  Mint a fresh ID (max(existing)+1, monotonic per-call)
#       and carry it forward through the later rows of the trajectory until
#       the next such break row.  Tag method = "bb_split" (first row) /
#       "bb_split_carry" (carried forward).
#
#   R2 (terminate-on-stump): a row with Status == "broken below" and
#       is.na(DBH) terminates the trajectory.  Any later row of the same
#       trajectory that has !is.na(DBH) MUST be re-IDed.  Mint a fresh ID
#       and propagate analogously.  Tag method = "bb_post_terminator_split" /
#       "bb_post_terminator_split_carry".  NA-DBH dead/BB/stem-dead corpse
#       rows already labelled by `apply_carried_terminal_backfill` keep
#       their ID (LOCF on terminator is allowed).
#
# Why the trajectory (and not StemTag) is the group: StemTags only exist from
# the 2010 census on (C7+). Grouping by StemTag put a stem's pre-2010 rows
# (StemTag NA) and its 2010+ rows (StemTag set) in different groups, so a
# break measured in 2010 was never compared with the trunk it continued
# (e.g. tag 082883: a 155 mm trunk and its 30 mm resprout stayed one stem),
# and in the all-NA pre-2010 group the rows of different stems interleaved,
# so a split could stop carrying forward early (e.g. tag 053376).
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
    data.table::setorderv(out, c("Tag", "ReconstructedStemID", "CensusID"))
    if (!("ReconstructionMethod" %in% names(out))) {
        out[, ReconstructionMethod := NA_character_]
    }
    cur_max <- suppressWarnings(max(out$ReconstructedStemID, na.rm = TRUE))
    if (!is.finite(cur_max)) cur_max <- 0L
    next_id <- as.integer(cur_max) + 1L

    has_true <- "TrueStemID" %in% names(out)
    rid <- as.integer(out$ReconstructedStemID)
    mth <- as.character(out$ReconstructionMethod)
    sts <- out$Status
    dbh <- out$DBH
    trueid <- if (has_true) as.integer(out$TrueStemID) else rep(NA_integer_, nrow(out))

    # Group row indices by trajectory (Tag, ReconstructedStemID), in census
    # order (the setorderv above). Rows without an identity are not grouped.
    groups <- out[!is.na(ReconstructedStemID), .(.idxs = list(.I)), by = .(Tag, ReconstructedStemID)]

    n_r1 <- 0L
    n_r2 <- 0L
    n_pin_override <- 0L

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
            new_id <- next_id
            next_id <- next_id + 1L
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
                new2 <- next_id
                next_id <- next_id + 1L
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
# apply_pins_to_samples()
#
# Per-sample version of the TrueStemID pin treatment the exported table gets
# after the engine: the pin sweep (pinned rows take their pin; a pin that would
# put two rows on one ID in a census is rolled back, as Fix 2 does) followed by
# apply_pin_track_rejoin() (the rest of a sampled track joins its single pin).
# Engines that do not honour every pin while sampling (the DP with a
# provisional anchor, the pre segment of a resprout split, pins joining a
# broken-below row to its continuation across the split) otherwise produce
# samples that group pinned rows differently from the export, so the export is
# not among the sampled paths. Only pinned links change: unpinned links, where
# the identity uncertainty lies, keep their per-sample variation.
#
# Pins are database identities. A row is pinned when TrueStemID is non-NA; its
# pin is TrueStemID, except on rows whose TrueStemID was rewritten by the
# broken-below pass or the renumbering (engine-minted methods), where it is the
# row's own StemID (every database pin is the row's own StemID).
#
# Rejoin, per sample: each sampled track is cut into segments at every
# broken-below row with a DBH and after every NA-R barrier census; when a
# segment's pinned rows all hold one pin P, the track's remaining rows join P
# if they all lie in that segment, none carries a pin, and P has no row in
# their censuses. Runs before apply_bb_invariants_to_samples() (R1/R2), as in
# the export. Labels are re-encoded per tag, consistently across samples.
#
# Inputs: samples_dt (Sample, CensusID, ReconstructedStemID, ObsRowID) and the
# tag's rows of the final table (obs_row_id, CensusID, DBH, Status, TrueStemID,
# ReconstructionMethod, and StemID or OriginalStemID).
# Returns samples_dt with ReconstructedStemID re-encoded.
apply_pins_to_samples <- function(samples_dt, tree_data) {
    if (is.null(samples_dt) || nrow(samples_dt) == 0L || is.null(tree_data) || nrow(tree_data) == 0L) {
        return(samples_dt)
    }
    src_col <- if ("StemID" %in% names(tree_data)) "StemID" else if ("OriginalStemID" %in% names(tree_data)) "OriginalStemID" else NULL
    need <- c("obs_row_id", "CensusID", "DBH", "Status", "TrueStemID")
    if (is.null(src_col) || !all(need %in% names(tree_data)) ||
        !all(c("Sample", "CensusID", "ReconstructedStemID", "ObsRowID") %in% names(samples_dt))) {
        return(samples_dt)
    }
    minted <- c("provisional_dp", "bb_split", "bb_split_carry", "bb_post_terminator_split", "bb_post_terminator_split_carry")
    td <- data.table::as.data.table(tree_data)
    mth <- if ("ReconstructionMethod" %in% names(td)) as.character(td$ReconstructionMethod) else rep(NA_character_, nrow(td))
    sts <- as.character(td$Status)
    tsm <- if ("ListOfTSM" %in% names(td)) as.character(td$ListOfTSM) else rep(NA_character_, nrow(td))
    pin_all <- ifelse(is.na(td$TrueStemID), NA_integer_,
        ifelse(mth %in% minted, suppressWarnings(as.integer(as.character(td[[src_col]]))), as.integer(td$TrueStemID))
    )
    meta <- unique(data.table::data.table(
        ObsRowID = td$obs_row_id, pin = pin_all,
        bb = !is.na(td$DBH) & !is.na(sts) & sts == "broken below"
    )[!is.na(ObsRowID)], by = "ObsRowID")
    if (!any(!is.na(meta$pin))) {
        return(samples_dt)
    }
    # NA-R barrier censuses: no measured row, but a stump / R-coded record
    r_na <- is.na(td$DBH) & ((!is.na(sts) & sts == "broken below") |
        (!is.na(tsm) & grepl("\\b(R|RP|RF|RT|QR|OR)\\b", tsm, perl = TRUE)))
    bar <- data.table::data.table(c = as.integer(td$CensusID), meas = !is.na(td$DBH), r_na = r_na)[,
        .(n_meas = sum(meas), n_rna = sum(r_na)),
        by = c
    ][n_meas == 0L & n_rna > 0L, sort(c)]

    s <- data.table::data.table(.row = seq_len(nrow(samples_dt)), Sample = samples_dt$Sample, c = as.integer(samples_dt$CensusID),
        ObsRowID = samples_dt$ObsRowID, lab = as.integer(samples_dt$ReconstructedStemID))
    s <- meta[s, on = "ObsRowID"]
    s[is.na(bb), bb := FALSE]
    # 1. sweep: pinned rows take their pin; duplicates of a pin within a census keep the sampled label
    s[, dup := !is.na(pin) & (duplicated(data.table::data.table(Sample, c, pin)) | duplicated(data.table::data.table(Sample, c, pin), fromLast = TRUE))]
    s[, holds := !is.na(pin) & !dup]
    s[, new := fifelse(holds, paste0("P", pin), paste0("S", lab))]
    # 2. rejoin: segments of each sampled track (cut at broken-below + DBH rows and after barriers)
    data.table::setorder(s, Sample, lab, c)
    s[, nb := if (length(bar) > 0L) findInterval(c - 0.5, bar) else 0L]
    s[, seg := cumsum(bb) + nb, by = .(Sample, lab)]
    sg <- s[, .(P = {
        p <- unique(pin[!is.na(pin)])
        if (length(p) == 1L && all(holds[!is.na(pin)])) paste0("P", p) else NA_character_
    }), by = .(Sample, lab, seg)]
    s <- sg[s, on = .(Sample, lab, seg)]
    rest <- s[new == paste0("S", lab)] # rows still carrying the sampled label after the sweep
    if (nrow(rest) > 0L) {
        mv <- rest[, .(ok = all(is.na(pin)) && data.table::uniqueN(seg) == 1L && !is.na(P[1L]), P = P[1L], cs = list(c)), by = .(Sample, lab)][ok == TRUE]
        if (nrow(mv) > 0L) {
            occ <- unique(s[holds == TRUE, .(Sample, P = new, c)])[, occupied := TRUE]
            cand <- mv[, .(c = unlist(cs)), by = .(Sample, lab, P)]
            cand[occ, on = .(Sample, P, c), occupied := i.occupied]
            cand[, dup2 := .N > 1L, by = .(Sample, P, c)]
            bad <- unique(cand[occupied %in% TRUE | dup2, .(Sample, lab)])
            mv <- mv[!bad, on = .(Sample, lab)]
            if (nrow(mv) > 0L) s[mv, on = .(Sample, lab), new := fifelse(new == paste0("S", lab), i.P, new)]
        }
    }
    # 3. re-encode: sampled labels keep their value; pins get values above every sampled label
    top <- suppressWarnings(max(as.integer(samples_dt$ReconstructedStemID), na.rm = TRUE))
    if (!is.finite(top)) top <- 0L
    pl <- sort(unique(s$new[startsWith(s$new, "P")]))
    s[, enc := fifelse(startsWith(new, "P"), top + match(new, pl), suppressWarnings(as.integer(substring(new, 2L))))]
    data.table::setorder(s, .row)
    out <- data.table::copy(samples_dt)
    out[, ReconstructedStemID := as.integer(s$enc)]
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
    if (is.null(samples_dt) || nrow(samples_dt) == 0L) {
        return(samples_dt)
    }
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
    if (!all(c("Status", "DBH") %in% names(tree_data))) {
        return(samples_dt)
    }

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
# Post-engine renumbering pass that assigns ReconstructedStemID values
# sequentially from 1 within each Tag, in chronological order of first
# appearance. Within a tie on first_census, the stem with the largest DBH
# at that census gets the lower ID. This yields a clean per-tag space
# 1..N where ID=1 is the earliest-appearing (and, on tie, largest) stem.
#
# Per-tag algorithm:
#   1. Collect all unique non-NA ReconstructedStemID values in the tag.
#   2. For each stem id, compute:
#        - first_census = min(CensusID) where ReconstructedStemID == id
#        - dbh_at_first = max DBH at that (id, first_census) cell
#                        (NA if all rows there have NA DBH)
#   3. Sort by (first_census ASC, dbh_at_first DESC, id ASC) so the
#      earliest-appearing stem with the largest first-census DBH gets
#      new_id = 1, the next gets 2, and so on through .N.
#   4. Apply the mapping to: ReconstructedStemID,
#      ReconstructedStemID_PreSweep, DP_PosteriorTop{k}ID columns
#      (always; in the last two, engine labels that are no longer a final ID
#      get their own numbers after .N, so they never read as another stem),
#      and TrueStemID ONLY on rows whose ReconstructionMethod
#      %in% ENGINE_MINTED_INTO_TRUESTEMID (provisional_dp, bb_split,
#      bb_split_carry, bb_post_terminator_split,
#      bb_post_terminator_split_carry — those rows hold engine-minted
#      values in TrueStemID and must follow the rename; rows holding
#      real DB IDs in TrueStemID are left untouched as ground truth).
#
# Posterior path files are no longer written by this function. The post-
# engine driver should call finalize_posterior_paths() afterwards, passing
# the returned `mapping`, to translate staged per-sample samples_dt into
# the renumbered ID space and write the final paths file. Because the
# mapping interface is unchanged (Tag/old_id/new_id), posterior paths
# automatically pick up the new 1..N numbering.
#
# Inputs:
#   out                      : data.table (per-chunk or per-tag) with at
#                              least Tag, CensusID, DBH,
#                              ReconstructedStemID, ReconstructionMethod.
#   posterior_top_k          : (optional) explicit number of
#                              DP_PosteriorTop{k}ID columns. If NULL,
#                              auto-detected from column names.
#   posterior_samples_path   : retained for backward compatibility; no-op.
#   mapping_format           : retained for backward compatibility; no-op.
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

    # Methods that write engine-minted IDs into TrueStemID. For these
    # rows TrueStemID is renumbered alongside ReconstructedStemID;
    # rows holding real DB IDs in TrueStemID are left untouched.
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
    has_dbh <- "DBH" %in% names(out)

    # Local copies (single-pass updates avoid repeated data.table writes)
    rid <- as.integer(out$ReconstructedStemID)
    trueid <- if (has_true) as.integer(out$TrueStemID) else NULL
    presweep <- if (has_pre_sweep) as.integer(out$ReconstructedStemID_PreSweep) else NULL
    pst_cols <- lapply(posterior_id_cols, function(cn) as.integer(out[[cn]]))
    names(pst_cols) <- posterior_id_cols

    census_vec <- as.integer(out$CensusID)
    method_vec <- as.character(out$ReconstructionMethod)
    dbh_vec <- if (has_dbh) as.numeric(out$DBH) else rep(NA_real_, nrow(out))

    # Group row indices by Tag (preserves order)
    tag_groups <- out[, .(.idxs = list(.I)), by = Tag]

    mapping_list <- vector("list", nrow(tag_groups))
    n_tags_renumbered <- 0L
    n_tags_no_stems <- 0L
    n_ids_total <- 0L

    for (g in seq_len(nrow(tag_groups))) {
        ix <- tag_groups$.idxs[[g]]
        if (length(ix) == 0L) next
        tag_id <- tag_groups$Tag[g]

        rid_g <- rid[ix]
        method_g <- method_vec[ix]
        census_g <- census_vec[ix]
        dbh_g <- dbh_vec[ix]
        trueid_g <- if (has_true) trueid[ix] else rep(NA_integer_, length(ix))

        all_recon <- unique(rid_g[!is.na(rid_g)])
        if (length(all_recon) == 0L) {
            n_tags_no_stems <- n_tags_no_stems + 1L
            next
        }

        # Per-stem first_census and DBH at that first census (max DBH if
        # multiple rows share the (id, first_census) cell).
        stem_dt <- data.table::data.table(
            id = rid_g, census = census_g, dbh = dbh_g
        )[!is.na(id)]
        stem_meta <- stem_dt[, .(first_census = min(census, na.rm = TRUE)),
            by = id
        ]
        stem_meta <- merge(
            stem_meta,
            stem_dt[, .(id, census, dbh)],
            by.x = c("id", "first_census"),
            by.y = c("id", "census"),
            all.x = TRUE
        )[, .(dbh_at_first = suppressWarnings(max(dbh, na.rm = TRUE))),
            by = .(id, first_census)
        ]
        # Replace -Inf from all-NA DBH groups with NA so ordering puts them last
        stem_meta[!is.finite(dbh_at_first), dbh_at_first := NA_real_]

        # Sort: chronological first, then larger DBH first within same census,
        # then by original id for determinism. NA dbh sorts last via na.last=TRUE.
        data.table::setorder(stem_meta, first_census, -dbh_at_first, id, na.last = TRUE)
        stem_meta[, new_id := seq_len(.N)]

        map_lookup <- as.integer(stem_meta$new_id)
        names(map_lookup) <- as.character(stem_meta$id)

        translate <- function(v) {
            if (is.null(v)) {
                return(NULL)
            }
            out_v <- v
            hit <- !is.na(v) & as.character(v) %in% names(map_lookup)
            if (any(hit)) {
                out_v[hit] <- map_lookup[as.character(v[hit])]
            }
            out_v
        }

        # Apply mapping
        rid[ix] <- translate(rid_g)
        # Audit columns (engine labels) can hold labels that are no longer a
        # final ID (post-engine steps relabelled the track, e.g.
        # apply_pin_track_rejoin()). Left raw they could equal another stem's
        # new 1..N number, so they get their own numbers after N, shared by
        # ReconstructedStemID_PreSweep and the Top-K ID columns. The returned
        # mapping (used for the posterior paths) is unchanged.
        audit_vals <- c(if (has_pre_sweep) presweep[ix], unlist(lapply(pst_cols, `[`, ix), use.names = FALSE))
        unmapped <- sort(unique(audit_vals[!is.na(audit_vals) & !(as.character(audit_vals) %in% names(map_lookup))]))
        audit_lookup <- map_lookup
        if (length(unmapped) > 0L) {
            extra <- length(map_lookup) + seq_along(unmapped)
            names(extra) <- as.character(unmapped)
            audit_lookup <- c(map_lookup, extra)
        }
        translate_audit <- function(v) {
            out_v <- v
            hit <- !is.na(v)
            if (any(hit)) out_v[hit] <- audit_lookup[as.character(v[hit])]
            out_v
        }
        if (has_pre_sweep) {
            presweep[ix] <- translate_audit(presweep[ix])
        }
        for (cn in posterior_id_cols) {
            pst_cols[[cn]][ix] <- translate_audit(pst_cols[[cn]][ix])
        }
        # TrueStemID: only on rows whose method indicates engine-minted-into-TrueStemID
        if (has_true) {
            engine_minted_rows <- !is.na(method_g) & (method_g %in% ENGINE_MINTED_INTO_TRUESTEMID)
            if (any(engine_minted_rows)) {
                trueid_sub <- trueid_g
                trueid_sub[engine_minted_rows] <- translate(trueid_sub[engine_minted_rows])
                trueid[ix[engine_minted_rows]] <- trueid_sub[engine_minted_rows]
            }
        }

        mapping_list[[g]] <- data.table::data.table(
            Tag = rep(tag_id, nrow(stem_meta)),
            old_id = as.integer(stem_meta$id),
            new_id = as.integer(stem_meta$new_id),
            first_census = as.integer(stem_meta$first_census)
        )
        n_tags_renumbered <- n_tags_renumbered + 1L
        n_ids_total <- n_ids_total + nrow(stem_meta)
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
            use.names = TRUE, fill = TRUE
        )
    } else {
        data.table::data.table(
            Tag = integer(0), old_id = integer(0),
            new_id = integer(0), first_census = integer(0)
        )
    }

    # Posterior companion mapping files are no longer written here (see
    # finalize_posterior_paths()). The `posterior_samples_path` and
    # `mapping_format` arguments are retained for backward compatibility
    # but are now no-ops.

    if (isTRUE(verbose)) {
        message(sprintf(
            "[renumber_engine_minted_ids] Renumbered %d tag(s) to 1..N (%d total IDs); %d tag(s) had no stems.",
            n_tags_renumbered, n_ids_total, n_tags_no_stems
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
    null_result <- function() {
        invisible(list(
            n_tags = 0L, n_written = 0L, n_failed = 0L, written_paths = character(0)
        ))
    }
    if (is.null(posterior_samples_path) || !nzchar(posterior_samples_path)) {
        return(null_result())
    }
    staging_dir <- file.path(posterior_samples_path, "posteriors", ".staging")
    if (!dir.exists(staging_dir)) {
        return(null_result())
    }
    staging_files <- list.files(staging_dir,
        pattern = "^tag_.*_samples_raw_.*\\.rds$",
        full.names = TRUE
    )
    if (length(staging_files) == 0L) {
        return(null_result())
    }
    if (is.null(out) || nrow(out) == 0L) {
        warning(
            "[finalize_posterior_paths] `out` is empty; cannot finalize ",
            length(staging_files), " staging file(s)."
        )
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
            # IDs missing from the map (tracks used only in some samples,
            # sample-level break IDs) get fresh IDs above the renumbered range,
            # the same ID for the same raw value in every sample. Left raw, a
            # small engine ID could equal a renumbered ID and merge two
            # different stems into one label.
            unmapped <- !is.na(v) & !hit
            new_v <- v
            new_v[hit] <- lk[as.character(v[hit])]
            if (any(unmapped)) {
                raw_ids <- sort(unique(v[unmapped]))
                new_v[unmapped] <- max(lk) + match(v[unmapped], raw_ids)
            }
            if (any(hit) || any(unmapped)) {
                samples_dt[, ReconstructedStemID := as.integer(new_v)]
            }
        }

        # ---- (2a) Database pins, as the export applies them after the engine ----
        samples_dt <- apply_pins_to_samples(samples_dt, tree_data_for_tag)

        # ---- (2) Re-run bb invariants in renumbered ID space ----
        samples_dt <- apply_bb_invariants_to_samples(
            samples_dt, tree_data_for_tag,
            verbose = FALSE
        )

        # ---- (3) path_sig + paths_summary ----
        sample_sigs <- samples_dt[, .(path_sig = paste0(ReconstructedStemID, collapse = "-")),
            by = Sample
        ]
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
                by = "Sample", all.x = TRUE
            )
            path_probs <- sample_sigs[, .(path_prob = sum(sample_prob, na.rm = TRUE)),
                by = path_sig
            ]
            paths_summary <- merge(path_counts, path_probs, by = "path_sig", all.x = TRUE)
        } else {
            paths_summary <- path_counts[, .(path_count,
                path_prob = path_count / n_samp
            ),
            by = path_sig
            ]
        }

        # Compact reconstruction mapping (one row per path_sig)
        sigs_by_sample <- sample_sigs[, .(Sample, path_sig)]
        samples_with_sig <- merge(samples_dt[, .(Sample, ObsRowID, ReconstructedStemID)],
            sigs_by_sample,
            by = "Sample", all.x = TRUE
        )
        recon_by_path <- samples_with_sig[, .(recon = paste0(ObsRowID, ":",
            ReconstructedStemID,
            collapse = ";"
        )),
        by = .(path_sig, Sample)
        ]
        recon_compact <- recon_by_path[, .SD[1], by = path_sig, .SDcols = "recon"]
        paths_summary <- merge(paths_summary, recon_compact,
            by = "path_sig",
            all.x = TRUE
        )

        # ---- (4) Write final paths file ----
        out_dir_post <- file.path(info$posterior_samples_path, "posteriors")
        if (!dir.exists(out_dir_post)) {
            dir.create(out_dir_post, recursive = TRUE, showWarnings = FALSE)
        }
        out_path_base <- file.path(out_dir_post, paste0(
            "tag_", ifelse(is.na(tag_val), "NA", tag_val),
            "_posterior_samples_", ts_local
        ))
        wrote <- tryCatch(
            {
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
            },
            error = function(e) {
                warning(sprintf(
                    "[finalize_posterior_paths] Tag %s: write failed: %s",
                    as.character(tag_val), e$message
                ))
                NA_character_
            }
        )
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
