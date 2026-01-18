############################################################
### Global DP stem-ID reconstruction (backward)
############################################################
#
# What this script solves
# - Each tree (Tag) can have multiple stems.
# - In a later “anchor” census, each observed stem has a known `TrueStemID`.
# - In earlier censuses, stems are observed (`DBH`) but their identity (which anchor
#   stem they correspond to) is unknown.
# - Goal: assign a stable `ReconstructedStemID` to each earlier observation so each
#   stem forms a plausible DBH trajectory through time.
#
# Why global dynamic programming (DP)
# - Stepwise matching can make locally good links that cause global ID swaps later.
# - DP chooses a *whole multi-census assignment* that minimizes total cost across
#   all adjacent census transitions.
#
# Core representation
# - Tracks: $K$ latent identity “slots”. Within each census, each observed stem maps
#   injectively to one track.
# - State: the injective mapping (observations -> tracks) for a single census.
#
# Important modeling limitation
# - Missing DBH is treated as “not observed” (NA). There is no explicit latent
#   “alive but unobserved” state.
#
# Transition cost (biological)
# - Per track, per adjacent censuses $t \to t+1$:
#   - NA -> DBH : recruitment penalty (lognormal on recruited DBH)
#   - DBH -> NA : disappearance penalty (hazard function of DBH and interval)
#   - DBH -> DBH: growth penalty (Normal likelihood on annual growth)
# - Small tie-break discourages unnecessary rank crossings when costs tie.
#
# Computational guardrails
# - State count grows as $P(K, n_{obs}) = K (K-1) \cdots (K-n_{obs}+1)$.
# - If the state space is too large, we fall back to a safe stepwise igraph matcher.
#
# Side effects
# - This file defines functions when sourced.
# - The example “main” block at the bottom is opt-in so sourcing won’t run IO.
#
# If you really want this script to clear your workspace when run, set:
#   options(dp_global_biol.clear_env = TRUE)
if (isTRUE(getOption("dp_global_biol.clear_env", FALSE))) {
    rm(list = ls())
}

if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required. Install it with install.packages('data.table')")
}

if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required. Install it with install.packages('igraph')")
}

Tag <- CensusID <- DBH <- TrueStemID <- ReconstructedStemID <- ConstraintViolation <- ReconstructionMethod <- NULL
ReferenceStemID <- ConstraintViolationFlag <- NULL
DP_MaxStatesPerCensus <- DP_MaxStatesCensusID <- DP_KUsed <- NULL
species <- NULL

count_injective_states <- function(K, n_obs) {
    # PURPOSE
    # - Upper bound helper for DP state-space size.
    # - Counts the number of injective assignments of `n_obs` labeled observations
    #   into `K` labeled tracks (identity slots).
    #
    # INPUTS
    # - K: integer-ish; number of available tracks.
    # - n_obs: integer-ish; number of observed stems in a census.
    #
    # OUTPUT
    # - Numeric scalar (double): P(K, n_obs) = K*(K-1)*...*(K-n_obs+1) = K!/(K-n_obs)!
    # - Returns 0 when n_obs > K.
    # - Allows Inf on overflow.
    #
    # NOTES
    # - This grows factorially and is the main reason global DP can become expensive.
    if (!is.finite(K) || !is.finite(n_obs)) {
        return(NA_real_)
    }
    K <- as.integer(K)
    n_obs <- as.integer(n_obs)
    if (n_obs < 0L || K < 0L) {
        return(NA_real_)
    }
    if (n_obs == 0L) {
        return(1)
    }
    if (n_obs > K) {
        return(0)
    }
    prod(seq.int(from = K, to = K - n_obs + 1L, by = -1L))
}

resolve_interval_years <- function(tree_data,
                                   interval_years = NULL,
                                   interval_col_candidates = c(
                                       "Bio_IntervalYears",
                                       "IntervalYears",
                                       "interval_years",
                                       "census_interval_years",
                                       "CensusIntervalYears"
                                   )) {
    # PURPOSE
    # - Determine the census interval (years) for the DP.
    # - Supports the newer workflow where interval is stored in `tree_data` as a
    #   column (similar to Bio_* columns), so callers don't have to pass it around.
    #
    # CONTRACT
    # - The DP currently assumes a single constant interval for all adjacent census
    #   transitions within a (Tag, species) group.
    # - If a candidate interval column contains multiple distinct finite values,
    #   we error rather than silently picking one.

    if (!is.null(interval_years)) {
        interval_years <- suppressWarnings(as.numeric(interval_years))
        if (!is.finite(interval_years) || is.na(interval_years) || interval_years <= 0) {
            stop("interval_years must be a positive finite number.", call. = FALSE)
        }
        return(interval_years)
    }

    if (is.null(tree_data) || !is.data.frame(tree_data)) {
        stop(
            "interval_years not provided and tree_data is not a data.frame/data.table; ",
            "provide interval_years or add a column like 'Bio_IntervalYears'.",
            call. = FALSE
        )
    }

    for (col in interval_col_candidates) {
        if (!(col %in% names(tree_data))) next
        v <- suppressWarnings(as.numeric(tree_data[[col]]))
        v <- v[is.finite(v) & !is.na(v)]
        if (length(v) == 0L) next
        u <- unique(v)
        if (length(u) == 1L) {
            if (u[[1L]] <= 0) {
                stop("Interval column '", col, "' must be positive.", call. = FALSE)
            }
            return(u[[1L]])
        }
        stop(
            "Interval column '", col, "' has multiple distinct values (", paste(head(u, 10L), collapse = ", "),
            if (length(u) > 10L) ", ..." else "",
            "). The DP currently requires a single constant interval_years per group.",
            call. = FALSE
        )
    }

    stop(
        "interval_years not provided and no interval column found in tree_data. ",
        "Add a constant column like 'Bio_IntervalYears' (recommended) or pass interval_years explicitly.",
        call. = FALSE
    )
}

add_constraint_violation <- function(x, id_col = "ReconstructedStemID", min_growth, max_growth, interval_years = NULL) {
    # PURPOSE
    # - Post-hoc diagnostic: flag potentially implausible links along each reconstructed
    #   track when the *implied* per-year growth between adjacent censuses falls outside
    #   [min_growth, max_growth].
    #
    # INPUTS
    # - x: data.table with at least id_col, DBH, CensusID.
    # - id_col: which ID column defines a track (defaults to ReconstructedStemID).
    # - min_growth/max_growth: allowable annual growth bounds (cm/year).
    # - interval_years: years between consecutive censuses (assumed constant here).
    #
    # OUTPUT
    # - `x` with/updated `ConstraintViolation` logical column (TRUE for flagged rows).
    #
    # NOTES
    # - Only evaluates consecutive censuses (CensusID increases by 1).
    # - Flags the earlier observation in the violating pair (so the “bad link” is easy
    #   to see in time series plots).
    if (!("ConstraintViolation" %in% names(x))) {
        x[, ConstraintViolation := NA]
    }
    if (!all(c(id_col, "DBH", "CensusID") %in% names(x))) {
        return(x)
    }

    interval_years <- resolve_interval_years(x, interval_years = interval_years)
    data.table::setorder(x, CensusID)
    ids <- unique(x[[id_col]])
    ids <- ids[!is.na(ids)]
    if (length(ids) == 0L) {
        return(x)
    }
    for (sid in ids) {
        ii <- which(x[[id_col]] == sid & !is.na(x$DBH))
        if (length(ii) < 2L) next
        ii <- ii[order(x$CensusID[ii])]
        for (k in seq_len(length(ii) - 1L)) {
            i0 <- ii[k]
            i1 <- ii[k + 1L]
            if (x$CensusID[i1] != x$CensusID[i0] + 1L) next
            g <- (x$DBH[i1] - x$DBH[i0]) / interval_years
            cond <- isTRUE((g < min_growth) | (g > max_growth))
            if (cond || isTRUE(x$ConstraintViolation[i0])) {
                x$ConstraintViolation[i0] <- TRUE
            }
        }
    }
    x
}

add_dp_posterior_bins <- function(
  x,
  confident_prob = 0.95,
  unlinked_prob = 0.50,
  use_reconstructed_prob = TRUE,
  out_col = "DP_PosteriorBin"
) {
    # PURPOSE
    # - Convenience label per observation based on the DP marginal posterior.
    #
    # BINS
    # - "unlinked-likely": posterior probability of being unlinked >= unlinked_prob
    # - "confident": posterior probability of the chosen reconstructed ID >= confident_prob
    #   (or Top1Prob if use_reconstructed_prob=FALSE)
    # - "ambiguous": everything else (posterior spread across multiple IDs)
    #
    # NOTES
    # - Rows with missing posterior columns (e.g., anchor / given IDs) get NA.
    # - This is a summary; it does not change `ReconstructedStemID`.
    if (!data.table::is.data.table(x)) {
        x <- data.table::as.data.table(x)
    }

    if (!is.character(out_col) || length(out_col) != 1L || !nzchar(out_col)) {
        stop("out_col must be a single non-empty column name")
    }

    # Minimal required posterior columns
    if (!("DP_PosteriorUnlinkedProb" %in% names(x))) {
        x[, (out_col) := NA_character_]
        x[, (out_col) := factor(get(out_col), levels = c("confident", "ambiguous", "unlinked-likely"))]
        return(x)
    }

    score_col <- if (isTRUE(use_reconstructed_prob) && ("DP_PosteriorReconstructedProb" %in% names(x))) {
        "DP_PosteriorReconstructedProb"
    } else if ("DP_PosteriorTop1Prob" %in% names(x)) {
        "DP_PosteriorTop1Prob"
    } else {
        NA_character_
    }

    if (!is.character(score_col) || is.na(score_col)) {
        x[, (out_col) := NA_character_]
        x[, (out_col) := factor(get(out_col), levels = c("confident", "ambiguous", "unlinked-likely"))]
        return(x)
    }

    # Classify
    x[, (out_col) := {
        p_unlinked <- DP_PosteriorUnlinkedProb
        p_score <- get(score_col)

        # Only label when we actually have posterior numbers
        has_post <- is.finite(p_unlinked) | is.finite(p_score)

        out <- rep(NA_character_, .N)
        out[has_post] <- "ambiguous"
        out[has_post & is.finite(p_unlinked) & (p_unlinked >= unlinked_prob)] <- "unlinked-likely"
        out[has_post & is.finite(p_score) & (p_score >= confident_prob) &
            !(is.finite(p_unlinked) & (p_unlinked >= unlinked_prob))] <- "confident"
        out
    }]

    x[, (out_col) := factor(get(out_col), levels = c("confident", "ambiguous", "unlinked-likely"))]
    x
}

match_stems_optimal_backward <- function(tree_data, min_growth, max_growth, interval_years = NULL, anchor_start) {
    # PURPOSE
    # - "Safe" fallback when global DP is too expensive: perform stepwise matching
    #   backward in time using a deterministic bipartite matching.
    # - This is intentionally conservative and fast; it is not guaranteed to be
    #   globally optimal across multiple censuses (unlike DP).
    #
    # HOW IT WORKS (high level)
    # - Starting from an anchor census (where TrueStemID is known), walk backward:
    #   for each adjacent pair c <- c+1, solve a one-step assignment from stems at
    #   census c to stems at census c+1.
    # - Only allow edges where implied annual growth is within [min_growth, max_growth].
    # - If a stem cannot be matched, assign a new ID.
    #
    # INPUTS
    # - tree_data: data.table for ONE (Tag, species) group.
    # - min_growth/max_growth: hard bounds on annual growth (cm/year).
    # - interval_years: years between adjacent censuses.
    # - anchor_start: integer CensusID to start backward matching from.
    #
    # OUTPUT
    # - tree_data with `ReconstructedStemID` filled and diagnostic columns:
    #   - `ReconstructionMethod`: "given" (from TrueStemID) or "igraph" (fallback)
    #   - `ConstraintViolation`: TRUE for links that required force/creation.
    #
    # NOTES
    # - This is conceptually similar to STEM_IDENTIFICATION_TEST/2_igraph_stepwise.R,
    #   but defined locally so DP code can fall back without sourcing other scripts.
    tree_data <- tree_data[order(CensusID)]
    interval_years <- resolve_interval_years(tree_data, interval_years = interval_years)
    if (!("ReconstructionMethod" %in% names(tree_data))) {
        tree_data[, ReconstructionMethod := NA_character_]
    }
    if (!("ConstraintViolation" %in% names(tree_data))) {
        tree_data[, ConstraintViolation := NA]
    }
    tree_data[!is.na(TrueStemID), `:=`(
        ReconstructedStemID = as.integer(TrueStemID),
        ReconstructionMethod = "given"
    )]

    current_max <- suppressWarnings(max(tree_data$ReconstructedStemID, na.rm = TRUE))
    if (!is.finite(current_max)) current_max <- 0L
    next_new_id <- as.integer(current_max + 1L)

    ensure_future_ids <- function(idx) {
        # Helper: ensure all "future" stems (at c+1) have a concrete ID so the
        # one-step match at census c can propagate that ID backward.
        if (length(idx) == 0L) {
            return()
        }
        missing <- is.na(tree_data$ReconstructedStemID[idx])
        if (any(missing)) {
            n <- sum(missing)
            tree_data$ReconstructedStemID[idx[missing]] <<- seq.int(from = next_new_id, length.out = n)
            idx_to_mark <- idx[missing]
            idx_to_mark <- idx_to_mark[is.na(tree_data$TrueStemID[idx_to_mark])]
            if (length(idx_to_mark) > 0L) {
                tree_data$ReconstructionMethod[idx_to_mark] <<- "igraph"
            }
            next_new_id <<- next_new_id + n
        }
    }

    compute_edge_weights <- function(dist_mat, curr_dbh, fut_dbh, valid_ij, n_curr, n_fut) {
        # Helper: deterministic edge weights for max_bipartite_match().
        #
        # Objective: prefer small |DBH_c - DBH_{c+1}| (distance) while gently
        # breaking ties to reduce unnecessary rank crossings.
        #
        # Implementation detail
        # - igraph matches by maximizing total weight.
        # - Convert distance to weight by: weight = (M - distance) with M > max(distance).
        M <- max(dist_mat[is.finite(dist_mat)], 0) + 1

        i_vec <- valid_ij[, 1]
        j_vec <- valid_ij[, 2]

        uniq_d <- sort(unique(as.vector(dist_mat[valid_ij])))
        min_gap <- suppressWarnings(min(diff(uniq_d), na.rm = TRUE))
        eps <- if (is.finite(min_gap) && min_gap > 0) min_gap / 1000 else 1e-9

        curr_rank <- rank(curr_dbh, ties.method = "first")
        fut_rank <- rank(fut_dbh, ties.method = "first")
        rank_gap <- abs(curr_rank[i_vec] - fut_rank[j_vec])

        (M - dist_mat[valid_ij]) -
            eps * rank_gap +
            (eps * 1e-3) * (n_fut - j_vec) +
            (eps * 1e-6) * (n_curr - i_vec)
    }

    build_bipartite_graph_valid_edges <- function(valid, dist_mat, curr_dbh, fut_dbh) {
        # Helper: build the bipartite graph with only growth-feasible edges.
        # - Left vertices: stems observed at census c.
        # - Right vertices: stems observed at census c+1.
        # - Add edge (i,j) only if growth(i->j) is within [min_growth, max_growth].
        n_curr <- length(curr_dbh)
        n_fut <- length(fut_dbh)
        g <- igraph::make_empty_graph(n = n_curr + n_fut, directed = FALSE)
        igraph::V(g)$type <- c(rep(TRUE, n_curr), rep(FALSE, n_fut))

        valid_ij <- which(valid, arr.ind = TRUE)
        if (nrow(valid_ij) == 0L) {
            return(list(g = g, n_curr = n_curr, n_fut = n_fut))
        }

        from <- valid_ij[, 1]
        to <- n_curr + valid_ij[, 2]
        g <- igraph::add_edges(g, as.vector(rbind(from, to)))
        igraph::E(g)$weight <- compute_edge_weights(dist_mat, curr_dbh, fut_dbh, valid_ij, n_curr, n_fut)
        list(g = g, n_curr = n_curr, n_fut = n_fut)
    }

    for (c in seq.int(from = anchor_start - 1L, to = 1L, by = -1L)) {
        # Work backward one census at a time (c <- c+1).
        # - curr_idx: observed stems at census c.
        # - fut_idx: observed stems at census c+1.
        curr_idx <- which(tree_data$CensusID == c & !is.na(tree_data$DBH))
        if (length(curr_idx) == 0L) next
        fut_idx <- which(tree_data$CensusID == c + 1L & !is.na(tree_data$DBH))
        ensure_future_ids(fut_idx)

        n_curr <- length(curr_idx)
        n_fut <- length(fut_idx)

        if (n_fut == 0L) {
            # If nothing is observed at c+1, we cannot link; assign new IDs.
            tree_data$ReconstructedStemID[curr_idx] <- seq.int(from = next_new_id, length.out = n_curr)
            tree_data$ConstraintViolation[curr_idx] <- TRUE
            curr_to_mark <- curr_idx[is.na(tree_data$TrueStemID[curr_idx])]
            if (length(curr_to_mark) > 0L) {
                tree_data$ReconstructionMethod[curr_to_mark] <- "igraph"
            }
            next_new_id <- next_new_id + n_curr
            next
        }

        curr_dbh <- tree_data$DBH[curr_idx]
        fut_dbh <- tree_data$DBH[fut_idx]

        dist_mat <- matrix(
            abs(outer(curr_dbh, fut_dbh, FUN = "-")),
            nrow = n_curr,
            ncol = n_fut
        )
        if (is.null(dim(dist_mat))) {
            dist_mat <- matrix(dist_mat, nrow = n_curr, ncol = n_fut)
        }

        growth_mat <- matrix(
            outer(curr_dbh, fut_dbh, FUN = function(d0, d1) (d1 - d0) / interval_years),
            nrow = n_curr,
            ncol = n_fut
        )

        valid <- (growth_mat >= min_growth) & (growth_mat <= max_growth)

        # Recruitment handling (same rule as 2_igraph_stepwise.R):
        # if there are more stems at c+1 than at c, treat the smallest DBH future stems
        # as recruits by forbidding any matches into them.
        if (n_fut > n_curr) {
            k <- n_fut - n_curr
            recruit_cols <- order(fut_dbh, decreasing = FALSE)[seq_len(k)]
            valid[, recruit_cols] <- FALSE
        }

        built <- build_bipartite_graph_valid_edges(valid, dist_mat, curr_dbh, fut_dbh)
        g <- built$g

        mbm <- igraph::max_bipartite_match(
            g,
            types = igraph::V(g)$type,
            weights = igraph::E(g)$weight
        )

        ids_assigned <- rep(NA_integer_, n_curr)
        violation_flag <- rep(NA, n_curr)

        m <- mbm$matching[seq_len(n_curr)]
        for (i in seq_len(n_curr)) {
            v_fut <- m[[i]]
            if (!is.na(v_fut)) {
                j_right <- as.integer(v_fut) - n_curr
                if (j_right >= 1L && j_right <= n_fut) {
                    ids_assigned[i] <- tree_data$ReconstructedStemID[fut_idx[j_right]]
                    violation_flag[i] <- !valid[i, j_right]
                }
            }
        }

        if (any(is.na(ids_assigned))) {
            n_new <- sum(is.na(ids_assigned))
            ids_assigned[is.na(ids_assigned)] <- seq.int(from = next_new_id, length.out = n_new)
            violation_flag[is.na(violation_flag)] <- TRUE
            next_new_id <- next_new_id + n_new
        }

        tree_data$ReconstructedStemID[curr_idx] <- ids_assigned
        tree_data$ConstraintViolation[curr_idx] <- violation_flag
        curr_to_mark <- curr_idx[is.na(tree_data$TrueStemID[curr_idx])]
        if (length(curr_to_mark) > 0L) {
            tree_data$ReconstructionMethod[curr_to_mark] <- "igraph"
        }
    }

    tree_data
}

enumerate_states_injective <- function(K, n_obs, max_states) {
    # PURPOSE
    # - Enumerate the DP "state space" for a single census: all injective mappings
    #   of n_obs labeled observations into K labeled tracks.
    #
    # INPUTS
    # - K: integer; number of tracks (identity slots).
    # - n_obs: integer; number of observed stems in this census.
    # - max_states: hard cap on number of states to enumerate.
    #
    # OUTPUT
    # - Matrix with one row per state and n_obs columns.
    #   Each row is a length-n_obs integer vector of track indices; e.g. c(2,5,1)
    #   means obs1->track2, obs2->track5, obs3->track1.
    # - Returns NULL when:
    #   - n_obs > K (impossible to assign injectively), or
    #   - the estimated number of states P(K, n_obs) exceeds max_states.
    #
    # COMPLEXITY
    # - Time and memory scale as O(P(K, n_obs) * n_obs).
    # - This can explode quickly; callers must provide a small max_states and fall
    #   back when enumeration is refused.
    if (n_obs == 0L) {
        return(matrix(integer(0), nrow = 1L, ncol = 0L))
    }

    # If we have more observations than tracks, it is impossible to assign
    # each observation to a unique track.
    if (n_obs > K) {
        return(NULL)
    }

    # Before we try to build the list, estimate how big it will be.
    # Number of injective mappings: P(K, n_obs) = K*(K-1)*...*(K-n_obs+1)
    n_states <- 1
    for (i in 0:(n_obs - 1L)) n_states <- n_states * (K - i)

    # Safety guard: if it's too big, we refuse and make the caller fall back.
    if (is.finite(n_states) && n_states > max_states) {
        return(NULL)
    }
    tracks <- seq_len(K)

    build <- function(prefix, remaining) {
        # Recursive constructor.
        #
        # prefix    = the track choices we have already made for obs 1..length(prefix)
        # remaining = which tracks are still available to use
        #
        # Base case:
        # - If we've chosen n_obs tracks, we have completed one full state.
        if (length(prefix) == n_obs) {
            return(list(prefix))
        }

        # Otherwise, we need to pick a track for the next observation.
        # We'll try each available track and recurse.
        out <- vector("list", 0L)
        for (t in remaining) {
            # Choose track t for the next observation.
            # Then remove t from remaining (cannot reuse it).
            out <- c(out, build(c(prefix, t), remaining[remaining != t]))
        }
        out
    }
    combos <- build(integer(0), tracks)
    do.call(rbind, lapply(combos, function(v) matrix(v, nrow = 1L)))
}

state_key <- function(state_vec) {
    # PURPOSE
    # - Canonical string key for a census "state" (injective observation->track map).
    #
    # INPUT
    # - state_vec: integer vector of length n_obs; entries are track indices.
    #
    # OUTPUT
    # - Single string like "2,5,1" used to index DP tables (named vectors/lists).
    if (length(state_vec) == 0L) {
        return("")
    }
    paste(state_vec, collapse = ",")
}

state_to_track_dbh <- function(state_vec, obs_dbh, K) {
    # PURPOSE
    # - Convert a census "state" into a track-indexed DBH vector so transition costs
    #   can be computed track-by-track.
    #
    # INPUTS
    # - state_vec: integer vector of length n_obs mapping obs index -> track index.
    # - obs_dbh: numeric vector of length n_obs; DBH values for each observation.
    # - K: total number of tracks.
    #
    # OUTPUT
    # - Numeric vector of length K where unused tracks are NA and used tracks receive
    #   the DBH of the observation assigned to that track.
    out <- rep(NA_real_, K)
    if (length(state_vec) == 0L) {
        return(out)
    }
    out[state_vec] <- obs_dbh
    out
}

############################################################
### BIO PARAMETER ESTIMATION (WHOLE DATASET)
############################################################
estimate_bio_pars <- function(
  x,
  interval_years,
  census_ids = NULL,
  mortality_start = c(log(0.01), 0),
  # -----------------------------------------------------------------
  # Measurement error model (DBH remeasurement)
  # -----------------------------------------------------------------
  # This is used in TWO places:
  #   (A) to correct growth-variance estimation (separating process vs measurement)
  #   (B) to set conservative guardrails (e.g., max_shrink)
  #
  # Reference (used as the basis for the measurement-error scaling):
  #   https://royalsocietypublishing.org/rstb/article/359/1443/409/20356/Error-propagation-and-scaling-for-tropical-forest
  #
  # Model assumption (mixture, per census measurement):
  #   With probability (1 - p_big): small measurement error
  #   With probability p_big: large measurement error (blunders)
  #
  # Small-error SD (cm) is diameter-dependent:
  #   SD1(D) = a * D + b
  # where D is DBH in cm.
  # (In our workflow, the fitted values are: a=0.0062, b=0.0904.)
  meas_sd1_a = 0.0062,
  meas_sd1_b = 0.0904,
  # Large-error SD (cm) and mixture weight
  meas_sd2 = 4.64,
  meas_p_big = 0.05,
  # Whether to correct growth-variance estimation for measurement error
  use_measurement_error = TRUE,
  # Quantiles used to set conservative guardrails
  shrink_hard_prob = 1e-4,
  shrink_data_quantile = 0.001,
  # Hard shrink guardrail (max_shrink)
  # - "data": estimated from observed shrink tail (with measurement-error support)
  # - "fixed": use a fixed constant bound (cm/year)
  max_shrink_source = c("data", "fixed"),
  max_shrink_fixed = -2,
  # Soft shrinkage penalty strength (k_shrink)
  # - "data": estimate from measurement-error scale (preferred) or from data variance
  # - "fixed": use a fixed constant (units: 1/cm^2)
  k_shrink_source = c("data", "fixed"),
  k_shrink_fixed = 50,
  # Soft extreme-growth penalty strength (k_growth), analogous to k_shrink
  # - "data": estimate from measurement-error scale (preferred) or from data variance
  # - "fixed": use a fixed constant (units: 1/cm^2); set 0 to disable soft penalty
  k_growth_source = c("data", "fixed"),
  k_growth_fixed = 50,
  # Extreme-growth guardrails (upper tail)
  # - growth_hard_prob is the *upper-tail* probability (e.g., 1e-4 means 99.99th percentile)
  # - growth_data_quantile is the empirical upper quantile used as a guardrail
  # - growth_soft_quantile sets a softer threshold used for a quadratic penalty
  growth_hard_prob = 1e-4,
  growth_data_quantile = 0.999,
  growth_soft_quantile = 0.99,
  # Hard growth guardrail (max_growth)
  # - "data": estimated from observed extreme-growth tail
  # - "fixed": use a fixed constant bound (cm/year)
  max_growth_source = c("data", "fixed"),
  max_growth_fixed = 7.5,
  # Recruitment max DBH (upper bound for recruits)
  recruit_max_quantile = 0.999
) {
    # =====================================================================
    # estimate_bio_pars()
    # =====================================================================
    # Goal
    #   Estimate a set of "biologically plausible" parameters that make the
    #   DP stem-tracking likelihood behave sensibly on a given dataset.
    #
    # What this function returns
    #   A nested list with:
    #   - growth:
    #       mu     : mean annual diameter increment (cm / year)
    #       sigma0 : baseline SD of annual increment (cm / year)
    #       sigma1 : slope for SD vs DBH (cm / year per cm DBH)
    #   - mortality:
    #       h0, beta : parameters of a hazard model for DBH -> NA transitions
    #   - recruitment:
    #       meanlog, sdlog : lognormal parameters for recruit DBH (cm)
    #       recruit_max_dbh: hard-ish upper guardrail for recruit DBH (cm)
    #       lambda         : recruit rate per available slot per year
    #   - shrinkage:
    #       k_shrink  : soft penalty weight for shrinkage (1/cm^2)
    #       max_shrink: conservative lower bound on annual growth (cm / year)
    #   - measurement_error:
    #       echo of SD1/SD2/p_big settings used
    #   - settings:
    #       echo of use_measurement_error
    #
    # How these parameters are used later
    #   They are passed into transition_cost_tracks_bio(), which scores a single
    #   transition between two adjacent censuses. The DP solver then finds the
    #   lowest-cost set of trajectories consistent with an anchor census.
    #
    # Key modeling assumptions
    #   1) Growth increments (annualized) are approximately Normal with
    #      heteroskedastic SD and a DBH-dependent mean:
    #        g = (DBH_{t+1} - DBH_t) / T
    #        g | DBH_t ~ Normal(mu(DBH_t), sigma(DBH_t)^2)
    #        mu(DBH)    = alpha + gamma * log(DBH)
    #        sigma(DBH) = sigma0 + sigma1 * DBH
    #   2) DBH measurement has nontrivial noise. If enabled, we treat measured
    #      DBH as:
    #        DBH_obs = DBH_true + epsilon
    #      where epsilon is a mixture of a "small" Normal error and a "large"
    #      Normal error (blunders). The small SD increases with DBH.
    #   3) Mortality is modeled as a hazard over the census interval:
    #        hazard(DBH) = h0 * exp(beta * DBH)
    #        P(death over interval T) = 1 - exp(-hazard(DBH) * T)
    #   4) Recruitment is modeled as NA -> DBH with:
    #        - recruit sizes ~ LogNormal(meanlog, sdlog)
    #        - recruit rate lambda per available NA slot per year
    #
    # Practical intent
    #   This is NOT meant to be a perfect ecological model. It is meant to:
    #   - make the cost function realistic enough that the DP doesn't "cheat"
    #     by preferring ID swaps / forced deaths / forced recruits.
    #   - provide reasonable defaults for sensitivity analysis.

    # ---------------------------------------------------------------------
    # Inputs / arguments (detailed)
    # ---------------------------------------------------------------------
    # x
    #   A data.frame/data.table with at least:
    #     - Tag        : group identifier (integer-like)
    #     - CensusID   : census index (integer-like, increasing)
    #     - DBH        : diameter at breast height (cm)
    #     - TrueStemID : a "ground truth" stem identity used ONLY for parameter estimation
    #     - species    : species code (optional; if missing, caller should add)
    #
    # interval_years
    #   Numeric; time between adjacent censuses used for annualization.
    #   Example: if censuses are 5 years apart, interval_years=5.
    #
    # census_ids
    #   Optional integer vector of CensusID values to use. If NULL, uses all
    #   CensusIDs present after filtering.
    #
    # mortality_start
    #   Starting values for optim() in the mortality fit:
    #     c(log(h0_start), beta_start)
    #
    # meas_sd1_a, meas_sd1_b, meas_sd2, meas_p_big
    #   Parameters of the DBH measurement error model.
    #   We use the (linear-in-DBH) small-error SD:
    #     SD1(D) = a*D + b
    #   and a "large error" SD2 (cm) with mixture probability p_big.
    #
    #   These values follow the remeasurement-based error propagation discussion
    #   in the paper linked above (Royal Society Phil. Trans. B article).
    #   We keep the model in code form only (not reproducing paper text).
    #
    # use_measurement_error
    #   If TRUE:
    #     - subtract expected measurement variance when estimating growth SD
    #     - set shrinkage guardrails informed by measurement error
    #   If FALSE:
    #     - treat all growth variability as process variability
    #     - set shrinkage guardrails purely from data quantiles
    #
    # shrink_hard_prob
    #   Lower-tail probability used for the measurement-noise-only shrink quantile.
    #   Smaller values make max_shrink more permissive (more negative).
    #
    # shrink_data_quantile
    #   Lower quantile of observed annual increments used as a data-driven guardrail.
    #   Smaller values make max_shrink more permissive (more negative).
    #
    # recruit_max_quantile
    #   Upper quantile of observed recruits used as a guardrail for recruit size.
    #   Larger values make recruit_max_dbh more permissive.

    library(data.table)
    library(MASS)

    interval_years <- as.numeric(interval_years)
    if (!is.finite(interval_years) || interval_years <= 0) {
        stop("interval_years must be positive.", call. = FALSE)
    }

    # ---------------------------------------------------------------------
    # 1) Clean input
    # ---------------------------------------------------------------------
    # We only use rows where:
    # - DBH is observed and positive
    # - TrueStemID is present (this function is estimating parameters from
    #   "known tracks" based on TrueStemID)
    #
    # This is a key point: estimate_bio_pars() uses TrueStemID to assemble
    # empirical growth/mortality/recruitment signals. If TrueStemID is missing
    # or unreliable, you should not trust the resulting parameters.
    dt <- as.data.table(x)[
        !is.na(DBH) & DBH > 0 & !is.na(TrueStemID),
        .(
            Tag,
            TrueStemID = as.integer(TrueStemID),
            CensusID   = as.integer(CensusID),
            DBH        = as.numeric(DBH),
            species    = as.character(species)
        )
    ]

    if (nrow(dt) == 0) {
        stop("No usable rows after filtering.", call. = FALSE)
    }

    if (is.null(census_ids)) {
        census_ids <- sort(unique(dt$CensusID))
    }

    census_ids <- census_ids[is.finite(census_ids)]
    if (length(census_ids) < 2) {
        stop("Need at least two censuses.", call. = FALSE)
    }

    # ---------------------------------------------------------------------
    # 2) Wide format
    # ---------------------------------------------------------------------
    # Convert to one row per (Tag, TrueStemID, species) with one column per census.
    # This makes it easy to compute adjacent-census transitions.
    dw <- dcast(
        dt,
        Tag + TrueStemID + species ~ CensusID,
        value.var = "DBH"
    )

    # ---------------------------------------------------------------------
    # 3) Growth increments
    # ---------------------------------------------------------------------
    # For each adjacent census pair (t0, t1):
    #   d0 = DBH at t0 (cm)
    #   d1 = DBH at t1 (cm)
    #   g  = (d1 - d0) / interval_years  (cm/year)
    #
    # We collect:
    #   g_all           : all observed annualized increments
    #   d0_all, d1_all  : the corresponding DBHs for regression/diagnostics
    #   var_meas_g_all  : expected measurement-variance contribution to g
    #
    # Measurement error model (mixture)
    #   We treat per-census DBH measurement error epsilon as:
    #     epsilon ~ (1-p) Normal(0, SD1(DBH)^2) + p Normal(0, SD2^2)
    #   where SD1(DBH) is linear in DBH.
    #
    #   For annualized increments g = (d1 - d0)/T, if we assume independent
    #   measurement errors at t0 and t1, then the variance of the annualized
    #   measurement difference is approximately:
    #     Var( (e1 - e0)/T ) = (Var(e1) + Var(e0)) / T^2
    #   where Var(e) for the mixture is:
    #     Var(e) = (1-p)*SD1^2 + p*SD2^2
    g_all <- c()
    d0_all <- c()
    d1_all <- c()
    var_meas_g_all <- c()

    sd1 <- function(d) {
        d <- as.numeric(d)
        pmax(meas_sd1_a * d + meas_sd1_b, 1e-6)
    }

    # Var(epsilon) under the mixture model.
    meas_var_eps <- function(d) {
        s1 <- sd1(d)
        (1 - meas_p_big) * (s1^2) + meas_p_big * (meas_sd2^2)
    }

    for (i in seq_len(length(census_ids) - 1)) {
        t0 <- as.character(census_ids[i])
        t1 <- as.character(census_ids[i + 1])

        if (!all(c(t0, t1) %in% names(dw))) next

        ok <- !is.na(dw[[t0]]) & !is.na(dw[[t1]])

        d0 <- dw[[t0]][ok]
        d1 <- dw[[t1]][ok]
        g <- (d1 - d0) / interval_years

        v_meas_g <- (meas_var_eps(d0) + meas_var_eps(d1)) / (interval_years^2)

        g_all <- c(g_all, g)
        d0_all <- c(d0_all, d0)
        d1_all <- c(d1_all, d1)
        var_meas_g_all <- c(var_meas_g_all, v_meas_g)
    }

    if (length(g_all) < 5) {
        stop("Not enough growth observations.", call. = FALSE)
    }

    # ---------------------------------------------------------------------
    # 3a) Estimate mean growth: mu(DBH) = alpha + gamma*log(DBH)
    # ---------------------------------------------------------------------
    # We fit a simple size-dependent mean model using the starting DBH (d0).
    # For robustness on small datasets, we fall back to a constant mean.
    mu_hat <- mean(g_all)

    ok_mu <- is.finite(g_all) & is.finite(d0_all) & (d0_all > 0)
    if (sum(ok_mu) >= 10L && stats::var(log(d0_all[ok_mu])) > 0) {
        fit_mu <- stats::lm(g_all[ok_mu] ~ log(d0_all[ok_mu]))
        alpha_hat <- as.numeric(stats::coef(fit_mu)[1])
        gamma_hat <- as.numeric(stats::coef(fit_mu)[2])
        if (!is.finite(alpha_hat)) alpha_hat <- mu_hat
        if (!is.finite(gamma_hat)) gamma_hat <- 0
    } else {
        alpha_hat <- mu_hat
        gamma_hat <- 0
    }

    mu_pred <- rep(alpha_hat, length(g_all))
    mu_pred[ok_mu] <- alpha_hat + gamma_hat * log(d0_all[ok_mu])

    # ---------------------------------------------------------------------
    # 3b) Estimate growth variance (sigma0, sigma1)
    # ---------------------------------------------------------------------
    # We want the *process* SD of annual increments, not the observed SD which
    # includes measurement error.
    #
    # Target model:
    #   SD_process(g | d0) = sigma0 + sigma1*d0
    #
    # Observed increments include measurement noise, so we estimate a total SD
    # and then (optionally) subtract expected measurement variance in quadrature:
    #   SD_process^2 ≈ max( SD_total^2 - Var_meas(g), 0 )
    #
    # Robust SD estimation trick
    #   For X ~ Normal(0, sd^2), E|X| = sd*sqrt(2/pi).
    #   Rearranging yields an SD proxy:
    #     sd ≈ |X| * sqrt(pi/2)
    #   We apply this to residuals (g - mu_hat) to reduce sensitivity to outliers.
    resid_abs <- abs(g_all - mu_pred)
    # For Normal(0, sd^2), E|X| = sd*sqrt(2/pi), so sd_hat ≈ |X|*sqrt(pi/2)
    sd_total_hat <- resid_abs * sqrt(pi / 2)

    if (isTRUE(use_measurement_error)) {
        sd_proc_hat <- sqrt(pmax(sd_total_hat^2 - var_meas_g_all, 1e-8))
    } else {
        sd_proc_hat <- pmax(sd_total_hat, 1e-6)
    }

    # Fit a simple linear model for SD vs DBH.
    # Notes:
    # - sigma0_hat is constrained to be positive
    # - sigma1_hat is constrained to be non-negative
    fit_sd <- lm(sd_proc_hat ~ d0_all)
    sigma0_hat <- max(coef(fit_sd)[1], 1e-4)
    sigma1_hat <- max(coef(fit_sd)[2], 0)

    # ---------------------------------------------------------------------
    # 4) Shrinkage penalty (k_shrink)
    # ---------------------------------------------------------------------
    # Shrinkage here means d1 < d0 (negative increment), which can occur due to:
    # - real biological shrinkage (limited)
    # - measurement error (common)
    # - ID swaps / bad matches (the thing we want to discourage)
    #
    # In transition_cost_tracks_bio(), shrinkage is penalized softly as:
    #   cost_shrink = k_shrink * (d0 - d1)^2
    # Units:
    #   (d0 - d1) is cm; to make cost dimensionless, k_shrink has units 1/cm^2.
    #
    # Heuristic used here
    #   If measurement error is enabled, set k_shrink so that shrinkage of about
    #   one typical measurement SD costs O(1).
    #   If measurement error is disabled, fall back to a crude estimate based on
    #   the variance of negative increments.
    # Estimate a soft shrink penalty from measurement error scale.
    # k_shrink is applied as:  k_shrink * (d0-d1)^2  (units: 1/cm^2).
    # With measurement noise, small shrinkage can be expected; we therefore scale
    # k_shrink so that shrinkage of ~1 SD costs O(1).
    k_shrink_source <- match.arg(k_shrink_source)
    k_shrink_fixed <- as.numeric(k_shrink_fixed)
    if (identical(k_shrink_source, "fixed")) {
        if (!is.finite(k_shrink_fixed) || k_shrink_fixed < 0) {
            stop("k_shrink_fixed must be a finite non-negative number when k_shrink_source='fixed'.", call. = FALSE)
        }
        k_shrink_hat_est <- NA_real_
        k_shrink_hat <- k_shrink_fixed
    } else {
        sd_meas_diff <- sqrt((meas_var_eps(d0_all) + meas_var_eps(d1_all)))
        sd_meas_diff <- sd_meas_diff[is.finite(sd_meas_diff) & sd_meas_diff > 0]
        if (isTRUE(use_measurement_error) && length(sd_meas_diff) >= 10) {
            s_typ <- stats::median(sd_meas_diff, na.rm = TRUE)
            k_shrink_hat_est <- 1 / (2 * (s_typ^2))
            k_shrink_hat_est <- min(max(k_shrink_hat_est, 1e-6), 1e6)
        } else {
            # Fallback (no measurement-error model): estimate from observed shrink magnitudes
            # on the *DBH difference* scale (cm), so k_shrink retains units 1/cm^2.
            delta_shrink <- (d0_all - d1_all)[is.finite(d0_all) & is.finite(d1_all) & (d1_all < d0_all)]
            if (length(delta_shrink) >= 5 && is.finite(stats::var(delta_shrink)) && stats::var(delta_shrink) > 0) {
                k_shrink_hat_est <- 1 / (2 * stats::var(delta_shrink))
            } else {
                k_shrink_hat_est <- 50
            }
        }
        k_shrink_hat <- k_shrink_hat_est
    }

    # ---------------------------------------------------------------------
    # 5) Mortality model
    # ---------------------------------------------------------------------
    # Mortality is inferred from TrueStemID tracks as:
    #   alive at t0 (DBH observed) and missing at t1 (DBH NA)
    #
    # We fit a simple hazard model:
    #   hazard(DBH) = h0 * exp(beta * DBH)
    #   P(death over interval T) = 1 - exp(-hazard(DBH)*T)
    #
    # mortality_start is on the unconstrained scale used by optim:
    #   par[1] = log(h0), par[2] = beta
    d0_m <- c()
    died <- c()

    for (i in seq_len(length(census_ids) - 1)) {
        t0 <- as.character(census_ids[i])
        t1 <- as.character(census_ids[i + 1])

        if (!all(c(t0, t1) %in% names(dw))) next

        at_risk <- !is.na(dw[[t0]])

        d0_m <- c(d0_m, dw[[t0]][at_risk])
        died <- c(died, as.integer(is.na(dw[[t1]][at_risk])))
    }

    negloglik_mort <- function(par, d0, died) {
        h0 <- exp(par[1])
        beta <- par[2]

        hazard <- h0 * exp(beta * d0)
        p <- 1 - exp(-hazard * interval_years)
        p <- pmin(pmax(p, 1e-12), 1 - 1e-12)

        -sum(died * log(p) + (1 - died) * log(1 - p))
    }

    fit_m <- optim(
        mortality_start,
        negloglik_mort,
        d0 = d0_m,
        died = died,
        method = "BFGS"
    )

    h0_hat <- exp(fit_m$par[1])
    beta_hat <- fit_m$par[2]

    # ---------------------------------------------------------------------
    # 6) Recruitment model
    # ---------------------------------------------------------------------
    # We identify "recruitment events" from TrueStemID tracks as:
    #   missing at t0 (DBH NA) and observed at t1 (DBH > 0)
    #
    # We estimate:
    #   - recruit sizes via a lognormal fit
    #   - recruit_max_dbh as a high quantile guardrail
    #   - recruit rate lambda as recruits per available NA slot per year
    recruit_dbh <- c()
    n_risk <- 0
    n_rec <- 0

    for (i in seq_len(length(census_ids) - 1)) {
        t0 <- as.character(census_ids[i])
        t1 <- as.character(census_ids[i + 1])

        if (!all(c(t0, t1) %in% names(dw))) next

        at_risk <- is.na(dw[[t0]])
        d1_at_risk <- dw[[t1]][at_risk]
        # Recruits must have a positive observed size. Non-positive values will
        # break the lognormal fit and are not meaningful DBH measurements.
        rec <- at_risk
        rec[at_risk] <- !is.na(d1_at_risk) & is.finite(d1_at_risk) & (d1_at_risk > 0)

        recruit_dbh <- c(recruit_dbh, dw[[t1]][rec])

        n_risk <- n_risk + sum(at_risk)
        n_rec <- n_rec + sum(rec)
    }

    recruit_dbh <- recruit_dbh[is.finite(recruit_dbh) & recruit_dbh > 0]

    if (length(recruit_dbh) >= 2) {
        fit_r <- fitdistr(recruit_dbh, "lognormal")
        mu_r <- fit_r$estimate["meanlog"]
        sd_r <- fit_r$estimate["sdlog"]
    } else {
        mu_r <- log(2)
        sd_r <- 0.5
    }

    # Guardrail: prevent the DP from treating very large stems as recruits.
    recruit_max_dbh <- if (length(recruit_dbh) > 0) {
        as.numeric(stats::quantile(recruit_dbh, recruit_max_quantile, na.rm = TRUE))
    } else {
        5
    }

    # Recruitment rate (Poisson)
    lambda_hat <- if (n_risk > 0) {
        n_rec / (n_risk * interval_years)
    } else {
        0
    }

    # ---------------------------------------------------------------------
    # 7) Shrink hard bound (max_shrink) from measurement + data
    # ---------------------------------------------------------------------
    # max_shrink is a conservative lower bound on annual growth (cm/year).
    # It is used as a guardrail to reject *absurd* shrinkage that is very
    # unlikely to be real or due to measurement error.
    #
    # Construction
    #   max_shrink_data : a small quantile of observed growth increments
    #   max_shrink_meas : a lower quantile of the measurement-noise-only
    #                    annualized difference distribution
    #   max_shrink_hat  : the minimum of these (more conservative)
    #
    # Measurement-noise-only lower quantile
    #   We compute a lower quantile of (e1 - e0)/T where e0 and e1 follow the
    #   mixture model. This yields a realistic lower tail for shrinkage driven
    #   purely by measurement error.
    #
    # Important: this is a guardrail, not a hard ecological law.
    # Measurement-informed lower bound on annual growth (mostly to prevent absurd matches).
    # We compute a conservative lower quantile of the *measurement-noise-only* annualized
    # DBH difference, using a typical diameter.
    meas_lower_quantile_g <- function(p, d_typ) {
        p <- as.numeric(p)
        if (!is.finite(p) || p <= 0 || p >= 1) {
            return(NA_real_)
        }
        d_typ <- as.numeric(d_typ)
        if (!is.finite(d_typ) || d_typ <= 0) {
            return(NA_real_)
        }

        s_small <- sd1(d_typ)
        s_big <- meas_sd2
        w_small <- 1 - meas_p_big
        w_big <- meas_p_big

        # Four-component mixture for (e1-e0)/T
        sds <- c(
            sqrt(s_small^2 + s_small^2) / interval_years,
            sqrt(s_small^2 + s_big^2) / interval_years,
            sqrt(s_big^2 + s_small^2) / interval_years,
            sqrt(s_big^2 + s_big^2) / interval_years
        )
        wts <- c(w_small * w_small, w_small * w_big, w_big * w_small, w_big * w_big)

        cdf <- function(x) sum(wts * pnorm(x, mean = 0, sd = sds))

        lo <- -10
        hi <- 0
        # Expand if needed
        if (cdf(lo) > p) {
            lo2 <- -50
            if (cdf(lo2) > p) {
                return(lo2)
            }
            lo <- lo2
        }
        if (cdf(hi) < p) {
            hi <- 10
            if (cdf(hi) < p) {
                return(hi)
            }
        }

        out <- tryCatch(
            uniroot(function(x) cdf(x) - p, lower = lo, upper = hi, tol = 1e-6)$root,
            error = function(e) NA_real_
        )
        out
    }

    d_typ <- stats::median(d0_all[is.finite(d0_all) & d0_all > 0], na.rm = TRUE)
    max_shrink_meas <- if (isTRUE(use_measurement_error)) meas_lower_quantile_g(shrink_hard_prob, d_typ) else NA_real_
    max_shrink_data <- as.numeric(stats::quantile(g_all, shrink_data_quantile, na.rm = TRUE))
    max_shrink_hat_est <- if (isTRUE(use_measurement_error) && is.finite(max_shrink_meas)) {
        min(max_shrink_data, max_shrink_meas)
    } else {
        max_shrink_data
    }

    max_shrink_source <- match.arg(max_shrink_source)
    max_shrink_fixed <- as.numeric(max_shrink_fixed)
    if (identical(max_shrink_source, "fixed")) {
        if (!is.finite(max_shrink_fixed)) {
            stop("max_shrink_fixed must be a finite number when max_shrink_source='fixed'.", call. = FALSE)
        }
        max_shrink_hat <- max_shrink_fixed
    } else {
        max_shrink_hat <- max_shrink_hat_est
    }

    # ---------------------------------------------------------------------
    # 7b) Extreme-growth guardrails (upper tail)
    # ---------------------------------------------------------------------
    # We treat very large positive increments as a likely sign of an incorrect
    # match (ID swap, mis-measurement, etc.). To avoid the DP taking such edges,
    # we set a conservative (permissive) hard upper bound, and also provide a
    # softer threshold for a quadratic penalty.
    #
    # Construction mirrors the shrinkage guardrail, but for the upper tail:
    #   max_growth_data : upper quantile of observed annualized increments
    #   max_growth_meas : upper quantile of the measurement-noise-only mixture for (e1-e0)/T
    #   max_growth_hat  : max(max_growth_data, max_growth_meas)  (more permissive)
    meas_upper_quantile_g <- function(p, d_typ) {
        p <- as.numeric(p)
        if (!is.finite(p) || p <= 0 || p >= 1) {
            return(NA_real_)
        }
        d_typ <- as.numeric(d_typ)
        if (!is.finite(d_typ) || d_typ <= 0) {
            return(NA_real_)
        }

        s_small <- sd1(d_typ)
        s_big <- meas_sd2
        w_small <- 1 - meas_p_big
        w_big <- meas_p_big

        # Four-component mixture for (e1-e0)/T
        sds <- c(
            sqrt(s_small^2 + s_small^2) / interval_years,
            sqrt(s_small^2 + s_big^2) / interval_years,
            sqrt(s_big^2 + s_small^2) / interval_years,
            sqrt(s_big^2 + s_big^2) / interval_years
        )
        wts <- c(w_small * w_small, w_small * w_big, w_big * w_small, w_big * w_big)

        cdf <- function(x) sum(wts * pnorm(x, mean = 0, sd = sds))

        lo <- 0
        hi <- 10
        if (cdf(hi) < p) {
            hi2 <- 50
            if (cdf(hi2) < p) {
                return(hi2)
            }
            hi <- hi2
        }
        if (cdf(lo) > p) {
            lo2 <- -10
            if (cdf(lo2) > p) {
                return(lo2)
            }
            lo <- lo2
        }

        out <- tryCatch(
            uniroot(function(x) cdf(x) - p, lower = lo, upper = hi, tol = 1e-6)$root,
            error = function(e) NA_real_
        )
        out
    }

    growth_hard_prob <- as.numeric(growth_hard_prob)
    growth_data_quantile <- as.numeric(growth_data_quantile)
    growth_soft_quantile <- as.numeric(growth_soft_quantile)

    max_growth_data <- as.numeric(stats::quantile(g_all, growth_data_quantile, na.rm = TRUE))
    max_growth_meas <- if (isTRUE(use_measurement_error) && is.finite(d_typ)) {
        p_hi <- 1 - growth_hard_prob
        meas_upper_quantile_g(p_hi, d_typ)
    } else {
        NA_real_
    }
    max_growth_hat_est <- if (isTRUE(use_measurement_error) && is.finite(max_growth_meas)) {
        max(max_growth_data, max_growth_meas)
    } else {
        max_growth_data
    }

    max_growth_source <- match.arg(max_growth_source)
    max_growth_fixed <- as.numeric(max_growth_fixed)
    if (identical(max_growth_source, "fixed")) {
        if (!is.finite(max_growth_fixed) || max_growth_fixed <= 0) {
            stop("max_growth_fixed must be a finite positive number when max_growth_source='fixed'.", call. = FALSE)
        }
        max_growth_hat <- max_growth_fixed
    } else {
        max_growth_hat <- max_growth_hat_est
    }

    max_growth_soft_data <- as.numeric(stats::quantile(g_all, growth_soft_quantile, na.rm = TRUE))
    max_growth_soft_hat <- if (is.finite(max_growth_hat) && is.finite(max_growth_soft_data)) {
        min(max_growth_hat, max_growth_soft_data)
    } else {
        max_growth_soft_data
    }

    # Soft extreme-growth penalty strength (units: 1/cm^2), analogous to k_shrink.
    # In transition_cost_tracks_bio(), the penalty is applied to the *excess DBH* above
    # the soft cap in cm:
    #   excess = d1 - (d0 + max_growth_soft*T)
    #   cost_growth_soft = k_growth * excess^2
    # We therefore scale k_growth so that an excess of ~1 typical measurement SD costs O(1).
    sd_meas_diff2 <- sqrt((meas_var_eps(d0_all) + meas_var_eps(d1_all)))
    sd_meas_diff2 <- sd_meas_diff2[is.finite(sd_meas_diff2) & sd_meas_diff2 > 0]
    if (isTRUE(use_measurement_error) && length(sd_meas_diff2) >= 10) {
        s_typ2 <- stats::median(sd_meas_diff2, na.rm = TRUE)
        k_growth_hat_est <- 1 / (2 * (s_typ2^2))
        k_growth_hat_est <- min(max(k_growth_hat_est, 1e-6), 1e6)
    } else {
        # Fallback (no measurement-error model): estimate from observed DBH increments (cm)
        # so k_growth retains units 1/cm^2.
        delta_pos <- (d1_all - d0_all)[is.finite(d0_all) & is.finite(d1_all) & (d1_all > d0_all)]
        if (length(delta_pos) >= 5 && is.finite(stats::var(delta_pos)) && stats::var(delta_pos) > 0) {
            k_growth_hat_est <- 1 / (2 * stats::var(delta_pos))
        } else {
            k_growth_hat_est <- 50
        }
    }

    k_growth_source <- match.arg(k_growth_source)
    k_growth_fixed <- as.numeric(k_growth_fixed)
    if (identical(k_growth_source, "fixed")) {
        if (!is.finite(k_growth_fixed) || k_growth_fixed < 0) {
            stop("k_growth_fixed must be a finite non-negative number when k_growth_source='fixed'.", call. = FALSE)
        }
        k_growth_hat <- k_growth_fixed
    } else {
        k_growth_hat <- k_growth_hat_est
    }

    # ---------------------------------------------------------------------
    # 8) Return parameters
    # ---------------------------------------------------------------------
    # The returned structure is intentionally aligned with:
    # - bio_pars_to_transition_args()
    # - transition_cost_tracks_bio()
    # - realism_calibration.R diagnostics
    list(
        growth = list(
            # Mean annual growth model: mu(DBH) = alpha + gamma*log(DBH)
            alpha = alpha_hat,
            gamma = gamma_hat,
            # Backward-compatible summary (empirical mean of g)
            mu = mu_hat,
            sigma0 = sigma0_hat,
            sigma1 = sigma1_hat,
            # Extreme-growth guardrails
            max_growth_soft = max_growth_soft_hat,
            max_growth = max_growth_hat,
            k_growth = k_growth_hat,
            k_growth_source = k_growth_source,
            k_growth_fixed = k_growth_fixed,
            k_growth_estimated_value = k_growth_hat_est,
            # Diagnostics / provenance (analogous to shrinkage)
            max_growth_data_quantile = growth_data_quantile,
            max_growth_data_value = max_growth_data,
            max_growth_meas_prob = growth_hard_prob,
            max_growth_meas_value = max_growth_meas,
            max_growth_soft_quantile = growth_soft_quantile,
            max_growth_soft_value = max_growth_soft_data,
            max_growth_source = max_growth_source,
            max_growth_fixed = max_growth_fixed,
            max_growth_estimated_value = max_growth_hat_est,
            # -----------------------------------------------------------------
            # Uniform nested layout (organizational; flat fields above remain)
            # -----------------------------------------------------------------
            guardrails = list(
                hard = list(
                    value = max_growth_hat,
                    source = max_growth_source,
                    fixed = max_growth_fixed,
                    estimated_value = max_growth_hat_est,
                    data_quantile = growth_data_quantile,
                    data_value = max_growth_data,
                    meas_prob = growth_hard_prob,
                    meas_value = max_growth_meas
                ),
                soft = list(
                    value = max_growth_soft_hat,
                    quantile = growth_soft_quantile,
                    data_value = max_growth_soft_data
                )
            ),
            penalties = list(
                soft = list(
                    k = k_growth_hat,
                    source = k_growth_source,
                    fixed = k_growth_fixed,
                    estimated_value = k_growth_hat_est
                )
            )
        ),
        mortality = list(
            h0   = h0_hat,
            beta = beta_hat
        ),
        recruitment = list(
            meanlog = mu_r,
            sdlog = sd_r,
            recruit_max_dbh = recruit_max_dbh,
            lambda = lambda_hat
        ),
        shrinkage = list(
            k_shrink = k_shrink_hat,
            max_shrink = max_shrink_hat,
            k_shrink_source = k_shrink_source,
            k_shrink_fixed = k_shrink_fixed,
            k_shrink_estimated_value = k_shrink_hat_est,
            max_shrink_source = max_shrink_source,
            max_shrink_fixed = max_shrink_fixed,
            max_shrink_estimated_value = max_shrink_hat_est,
            max_shrink_data_quantile = shrink_data_quantile,
            max_shrink_data_value = max_shrink_data,
            max_shrink_meas_prob = shrink_hard_prob,
            max_shrink_meas_value = max_shrink_meas,
            # -----------------------------------------------------------------
            # Uniform nested layout (organizational; flat fields above remain)
            # -----------------------------------------------------------------
            guardrails = list(
                hard = list(
                    value = max_shrink_hat,
                    source = max_shrink_source,
                    fixed = max_shrink_fixed,
                    estimated_value = max_shrink_hat_est,
                    data_quantile = shrink_data_quantile,
                    data_value = max_shrink_data,
                    meas_prob = shrink_hard_prob,
                    meas_value = max_shrink_meas
                ),
                # For shrinkage, the soft penalty begins as soon as DBH decreases (d1 < d0),
                # i.e., the soft "threshold" is effectively 0 cm shrink.
                soft = list(
                    value = 0
                )
            ),
            penalties = list(
                soft = list(
                    k = k_shrink_hat,
                    source = k_shrink_source,
                    fixed = k_shrink_fixed,
                    estimated_value = k_shrink_hat_est
                )
            )
        ),
        measurement_error = list(
            sd1_a = meas_sd1_a,
            sd1_b = meas_sd1_b,
            sd2 = meas_sd2,
            p_big = meas_p_big
        ),
        settings = list(
            use_measurement_error = use_measurement_error
        )
    )
}

############################################################
### BIO TRANSITION COST
############################################################

transition_cost_tracks_bio <- function(
  track_dbh_t,
  track_dbh_tp1,
  interval_years,
  # -----------------------
  # GROWTH MODEL PARAMETERS
  # -----------------------
  # Mean annual diameter increment (cm / year)
  # mu(DBH) = mu_const + mu_gamma * log(DBH)
  # If mu_gamma == 0, this reduces to a constant mean.
  mu_const = Bio_Mu_Growth_unit,
  mu_gamma = 0,
  # Size-dependent growth variance:
  #   sigma(d) = sigma0 + sigma1 * d
  sigma0 = Bio_Sigma0_unit,
  sigma1 = Bio_Sigma1_unit,
  # Maximum biologically plausible shrinkage (cm / year)
  max_shrink = Bio_max_shrink_unit, #-0.2,
  # Soft penalty strength for shrinkage
  k_shrink = Bio_k_shrink_unit, # 50,
  # Maximum biologically plausible *positive* growth (cm / year)
  # (hard guardrail; defaults to Inf = disabled)
  max_growth = Inf,
  # Soft penalty threshold for extreme growth (cm / year)
  # If finite, excess growth above this threshold is penalized as a quadratic
  # in *DBH units* (cm^2), analogous to shrinkage.
  max_growth_soft = Inf,
  # Soft penalty strength for extreme growth (units: 1/cm^2)
  k_growth = 0,
  # -----------------
  # MEASUREMENT ERROR (optional)
  # -----------------
  use_measurement_error = FALSE,
  meas_sd1_a = 0.0062,
  meas_sd1_b = 0.0904,
  meas_sd2 = 4.64,
  meas_p_big = 0.05,
  # -------------------------
  # MORTALITY MODEL PARAMETERS
  # -------------------------
  h0 = Bio_H0_Mortality,
  beta = Bio_Beta_Mortality,
  # ----------------------------
  # RECRUITMENT MODEL PARAMETERS
  # ----------------------------
  recruit_meanlog = Bio_Recruit_Meanlog_unit,
  recruit_sdlog = Bio_Recruit_Sdlog_unit,
  recruit_max_dbh = Bio_Recruit_MaxDBH_unit,
  recruit_lambda = Bio_Recruitment_lambda_unit,
  # -----------------
  # DETERMINISTIC TIE-BREAK
  # -----------------
  eps_tiebreak = 1e-6
) {
    # PURPOSE
    # - Compute the negative log-likelihood ("cost") of transitioning from census t
    #   to census t+1, given a *track-indexed* DBH representation.
    # - This is the main scoring function used by the global DP solver.
    #
    # INPUTS (core)
    # - track_dbh_t: numeric vector length K; DBH for each track at census t, or NA.
    # - track_dbh_tp1: numeric vector length K; DBH for each track at census t+1, or NA.
    # - interval_years: years between censuses (T).
    #
    # MODEL (per track, per transition)
    # - There are 4 mutually exclusive cases:
    #   1) NA -> NA: no recruit during interval, penalty = -log(1 - p_recruit)
    #   2) NA -> DBH: recruit, penalty = -log(p_recruit) - log f_recruit(DBH)
    #   3) DBH -> NA: mortality/disappearance, penalty = -log(p_death(DBH, T))
    #   4) DBH -> DBH: growth, penalty from growth likelihood + shrinkage penalties
    #
    # SHRINKAGE HANDLING
    # - Hard guardrail: if g = (d1-d0)/T < max_shrink, add a large penalty and stop
    #   evaluating this track.
    # - Soft penalty: if d1 < d0, add k_shrink * (d0 - d1)^2 (discourages shrinkage
    #   without strictly forbidding small decreases).
    #
    # EXTREME GROWTH HANDLING (analogous to shrinkage)
    # - Hard guardrail: if g = (d1-d0)/T > max_growth, add a large penalty and stop.
    # - Soft penalty: if g > max_growth_soft, add k_growth * (d1 - (d0 + max_growth_soft*T))^2.
    #   This penalizes only the *excess* beyond the soft threshold.
    #
    # MEASUREMENT ERROR (optional)
    # - If use_measurement_error=TRUE, the DBH->DBH likelihood becomes a mixture:
    #   observed annual growth g_obs = g_true + (e1 - e0)/T.
    # - Each measurement error e_t is a 2-component Normal mixture; differencing
    #   yields a 4-component mixture for (e1-e0)/T.
    # - Convolving that with the process growth Normal yields a 4-component Normal
    #   mixture on g_obs. We evaluate that mixture with a stable log-sum-exp.
    #
    # TIE-BREAK (non-biological)
    # - eps_tiebreak adds a tiny penalty for rank crossings to make DP decoding more
    #   stable when two paths have nearly identical biological cost.
    #
    # OUTPUT
    # - Numeric scalar: total cost summed over tracks + (optional) tie-break.
    # ---------------------------------------------------------------------
    # Setup
    # ---------------------------------------------------------------------
    K <- length(track_dbh_t)
    cost <- 0

    # Recruitment probability over interval
    p_recruit <- 1 - exp(-recruit_lambda * interval_years)
    p_recruit <- pmin(pmax(p_recruit, 1e-12), 1 - 1e-12)

    # Mean growth function: mu(DBH) = alpha + gamma*log(DBH)
    mu_growth <- function(d) {
        if (!is.finite(mu_gamma) || mu_gamma == 0 || !is.finite(d) || d <= 0) {
            return(mu_const)
        }
        mu_const + mu_gamma * log(d)
    }

    log_sum_exp <- function(x) {
        m <- max(x)
        if (!is.finite(m)) {
            return(m)
        }
        m + log(sum(exp(x - m)))
    }

    meas_sd1 <- function(d) pmax(meas_sd1_a * d + meas_sd1_b, 1e-6)

    # ---------------------------------------------------------------------
    # Loop over tracks
    # ---------------------------------------------------------------------
    for (k in seq_len(K)) {
        d0 <- track_dbh_t[k]
        d1 <- track_dbh_tp1[k]

        # -------------------------------------------------------------
        # CASE 1: NA -> NA  (no tree, no recruit)
        # -------------------------------------------------------------
        if (is.na(d0) && is.na(d1)) {
            cost <- cost - log(1 - p_recruit)
            next
        }

        # -------------------------------------------------------------
        # CASE 2: NA -> DBH  (RECRUITMENT)
        # -------------------------------------------------------------
        # “Model whether a new individual appears, and if so, what size it has at census 2.”
        # This is a two-part (hurdle) model:
        # Whether recruitment occurs
        # What size the recruit has.

        # This recruitment formulation assumes:
        # Recruitment timing is irrelevant
        # Only presence/absence at census 2 matters
        # All recruits are detected
        # Recruit size distribution is stationary
        if (is.na(d0) && !is.na(d1)) {
            if (!is.finite(d1) || d1 <= 0 || d1 > recruit_max_dbh) {
                # Biologically impossible recruit
                cost <- cost + 1e6
            } else {
                cost <- cost -
                    log(p_recruit) -
                    dlnorm(d1, recruit_meanlog, recruit_sdlog, log = TRUE)
            }
            next
        }

        # -------------------------------------------------------------
        # CASE 3: DBH -> NA  (MORTALITY)
        # -------------------------------------------------------------
        # Assumptions
        # Constant hazard within the interval
        # DBH effect does not change between censuses
        # DBH measured at census 1 drives mortality
        # Growth after census 1 does not affect hazard
        # Independent individuals
        # No unobserved heterogeneity (frailty)
        if (!is.na(d0) && is.na(d1)) {
            # maybe: hazard <- h0 * exp(beta * (d0 - mean_dbh))
            hazard <- h0 * exp(beta * d0)
            p_death <- 1 - exp(-hazard * interval_years)

            p_death <- pmin(pmax(p_death, 1e-12), 1 - 1e-12)
            cost <- cost - log(p_death)
            next
        }

        # -------------------------------------------------------------
        # CASE 4: DBH -> DBH  (GROWTH)
        # -------------------------------------------------------------
        g <- (d1 - d0) / interval_years

        # ---------------- FIX A ----------------
        # Hard biological constraint on shrinkage
        if (is.finite(max_shrink) && (g < max_shrink)) {
            cost <- cost + 1e6
            next
        }

        # Hard biological constraint on extreme positive growth
        if (is.finite(max_growth) && (g > max_growth)) {
            cost <- cost + 1e6
            next
        }

        # ---------------- FIX B ----------------
        # Size-dependent growth variance
        sigma_d <- sigma0 + sigma1 * d0
        sigma_d <- pmax(sigma_d, 1e-6)

        mu <- mu_growth(d0)

        if (isTRUE(use_measurement_error)) {
            # Measurement-error-aware likelihood:
            # g_obs = g_true + (e1-e0)/T.
            # With e_t being a 2-component Normal mixture, the differenced error is a
            # 4-component Normal mixture. Convolving with the process Normal gives a
            # 4-component Normal mixture on g_obs with the same mean mu.
            s_small0 <- meas_sd1(d0)
            s_small1 <- meas_sd1(d1)
            s_big <- meas_sd2
            w_small <- 1 - meas_p_big
            w_big <- meas_p_big

            sd_meas_mix <- c(
                sqrt(s_small0^2 + s_small1^2) / interval_years,
                sqrt(s_small0^2 + s_big^2) / interval_years,
                sqrt(s_big^2 + s_small1^2) / interval_years,
                sqrt(s_big^2 + s_big^2) / interval_years
            )
            wt_meas_mix <- c(w_small * w_small, w_small * w_big, w_big * w_small, w_big * w_big)
            sd_tot <- sqrt(sigma_d^2 + sd_meas_mix^2)

            ll <- log(wt_meas_mix) + stats::dnorm(g, mean = mu, sd = sd_tot, log = TRUE)
            cost <- cost - log_sum_exp(ll)
        } else {
            # Gaussian growth likelihood
            cost <- cost +
                (g - mu)^2 / (2 * sigma_d^2) +
                log(sigma_d) +
                0.5 * log(2 * pi)
        }

        # ---------------- FIX C ----------------
        # Soft penalty for any shrinkage
        if (d1 < d0) {
            cost <- cost + k_shrink * (d0 - d1)^2
        }

        # Soft penalty for extreme positive growth (only above max_growth_soft)
        if (is.finite(max_growth_soft) && is.finite(k_growth) && k_growth > 0) {
            d1_soft_cap <- d0 + max_growth_soft * interval_years
            if (is.finite(d1_soft_cap) && d1 > d1_soft_cap) {
                cost <- cost + k_growth * (d1 - d1_soft_cap)^2
            }
        }
    }

    # ---------------------------------------------------------------------
    # Deterministic tie-break (non-biological)
    # ---------------------------------------------------------------------
    if (eps_tiebreak > 0) {
        r0 <- rank(track_dbh_t, ties.method = "first")
        r1 <- rank(track_dbh_tp1, ties.method = "first")
        both_obs <- !is.na(track_dbh_t) & !is.na(track_dbh_tp1)

        if (any(both_obs)) {
            cost <- cost +
                eps_tiebreak * sum(abs(r0[both_obs] - r1[both_obs]))
        }
    }

    return(cost)
}

transition_cost_tracks_bio_batch <- function(
  track_dbh_t,
  track_dbh_tp1,
  interval_years,
  # -----------------------
  # GROWTH MODEL PARAMETERS
  # -----------------------
  mu_const = Bio_Mu_Growth_unit,
  mu_gamma = 0,
  sigma0 = Bio_Sigma0_unit,
  sigma1 = Bio_Sigma1_unit,
  max_shrink = Bio_max_shrink_unit,
  k_shrink = Bio_k_shrink_unit,
  max_growth = Inf,
  max_growth_soft = Inf,
  k_growth = 0,
  # -----------------
  # MEASUREMENT ERROR (optional)
  # -----------------
  use_measurement_error = FALSE,
  meas_sd1_a = 0.0062,
  meas_sd1_b = 0.0904,
  meas_sd2 = 4.64,
  meas_p_big = 0.05,
  # -------------------------
  # MORTALITY MODEL PARAMETERS
  # -------------------------
  h0 = Bio_H0_Mortality,
  beta = Bio_Beta_Mortality,
  # ----------------------------
  # RECRUITMENT MODEL PARAMETERS
  # ----------------------------
  recruit_meanlog = Bio_Recruit_Meanlog_unit,
  recruit_sdlog = Bio_Recruit_Sdlog_unit,
  recruit_max_dbh = Bio_Recruit_MaxDBH_unit,
  recruit_lambda = Bio_Recruitment_lambda_unit,
  # -----------------
  # DETERMINISTIC TIE-BREAK
  # -----------------
  eps_tiebreak = 1e-6,
  hard_penalty = 1e6
) {
    # PURPOSE
    # - Batched version of transition_cost_tracks_bio().
    # - Computes the transition cost for many candidate next-track DBH vectors
    #   given a fixed current-track DBH vector.
    #
    # INPUTS
    # - track_dbh_t: numeric vector length K
    # - track_dbh_tp1: either
    #     * a numeric matrix with nrow = n_batch and ncol = K (each row is a candidate), OR
    #     * a list of numeric vectors length K
    #
    # OUTPUT
    # - numeric vector length n_batch

    K <- length(track_dbh_t)
    interval_years <- as.numeric(interval_years)
    if (!is.finite(interval_years) || interval_years <= 0) {
        stop("interval_years must be positive.", call. = FALSE)
    }

    if (is.list(track_dbh_tp1)) {
        if (length(track_dbh_tp1) < 1L) {
            return(numeric(0))
        }
        if (any(vapply(track_dbh_tp1, length, integer(1L)) != K)) {
            stop("All elements of track_dbh_tp1 list must have length K.", call. = FALSE)
        }
        mat_tp1 <- do.call(rbind, track_dbh_tp1)
    } else {
        mat_tp1 <- as.matrix(track_dbh_tp1)
    }

    if (ncol(mat_tp1) != K) {
        stop("track_dbh_tp1 must have K columns (same length as track_dbh_t).", call. = FALSE)
    }

    n_batch <- nrow(mat_tp1)
    if (n_batch < 1L) {
        return(numeric(0))
    }

    cost <- rep(0, n_batch)

    # Recruitment probability over interval
    p_recruit <- 1 - exp(-recruit_lambda * interval_years)
    p_recruit <- pmin(pmax(p_recruit, 1e-12), 1 - 1e-12)

    mu_growth_scalar <- function(d) {
        if (!is.finite(mu_gamma) || mu_gamma == 0 || !is.finite(d) || d <= 0) {
            return(mu_const)
        }
        mu_const + mu_gamma * log(d)
    }

    meas_sd1 <- function(d) pmax(meas_sd1_a * d + meas_sd1_b, 1e-6)

    # ---------------------------------------------------------------------
    # Loop over tracks (vectorized across candidates)
    # ---------------------------------------------------------------------
    for (k in seq_len(K)) {
        d0 <- track_dbh_t[[k]]
        d1 <- mat_tp1[, k]

        # CASE 1 + 2: NA -> *
        if (is.na(d0)) {
            mask_na_na <- is.na(d1)
            if (any(mask_na_na)) {
                cost[mask_na_na] <- cost[mask_na_na] - log(1 - p_recruit)
            }

            mask_na_dbh <- !mask_na_na
            if (any(mask_na_dbh)) {
                d1v <- d1[mask_na_dbh]
                hard <- (!is.finite(d1v)) | (d1v <= 0) | (d1v > recruit_max_dbh)
                if (any(hard)) {
                    cost[which(mask_na_dbh)[hard]] <- cost[which(mask_na_dbh)[hard]] + hard_penalty
                }
                ok <- !hard
                if (any(ok)) {
                    idx_ok <- which(mask_na_dbh)[ok]
                    d1_ok <- d1v[ok]
                    cost[idx_ok] <- cost[idx_ok] - log(p_recruit) - stats::dlnorm(d1_ok, recruit_meanlog, recruit_sdlog, log = TRUE)
                }
            }
            next
        }

        # CASE 3: DBH -> NA
        mask_dbh_na <- is.na(d1)
        if (any(mask_dbh_na)) {
            hazard <- h0 * exp(beta * d0)
            p_death <- 1 - exp(-hazard * interval_years)
            p_death <- pmin(pmax(p_death, 1e-12), 1 - 1e-12)
            cost[mask_dbh_na] <- cost[mask_dbh_na] - log(p_death)
        }

        # CASE 4: DBH -> DBH
        mask_dbh_dbh <- !mask_dbh_na
        if (!any(mask_dbh_dbh)) next

        idx <- which(mask_dbh_dbh)
        d1v <- d1[idx]
        g <- (d1v - d0) / interval_years

        hard <- rep(FALSE, length(g))
        if (is.finite(max_shrink)) {
            hard <- hard | (g < max_shrink)
        }
        if (is.finite(max_growth)) {
            hard <- hard | (g > max_growth)
        }

        if (any(hard)) {
            cost[idx[hard]] <- cost[idx[hard]] + hard_penalty
        }

        ok <- !hard
        if (any(ok)) {
            idx_ok <- idx[ok]
            d1_ok <- d1v[ok]
            g_ok <- g[ok]

            sigma_d <- sigma0 + sigma1 * d0
            sigma_d <- pmax(sigma_d, 1e-6)
            mu <- mu_growth_scalar(d0)

            if (isTRUE(use_measurement_error)) {
                s_small0 <- meas_sd1(d0)
                s_small1 <- meas_sd1(d1_ok)
                s_big <- meas_sd2
                w_small <- 1 - meas_p_big
                w_big <- meas_p_big

                # Build m x 4 matrices of component SDs and log-weights
                sd_meas_1 <- sqrt(s_small0^2 + s_small1^2) / interval_years
                sd_meas_2 <- sqrt(s_small0^2 + s_big^2) / interval_years
                sd_meas_3 <- sqrt(s_big^2 + s_small1^2) / interval_years
                sd_meas_4 <- sqrt(s_big^2 + s_big^2) / interval_years

                wt1 <- w_small * w_small
                wt2 <- w_small * w_big
                wt3 <- w_big * w_small
                wt4 <- w_big * w_big

                sd_tot_1 <- sqrt(sigma_d^2 + sd_meas_1^2)
                sd_tot_2 <- sqrt(sigma_d^2 + sd_meas_2^2)
                sd_tot_3 <- sqrt(sigma_d^2 + sd_meas_3^2)
                sd_tot_4 <- sqrt(sigma_d^2 + sd_meas_4^2)

                ll1 <- log(wt1) + stats::dnorm(g_ok, mean = mu, sd = sd_tot_1, log = TRUE)
                ll2 <- log(wt2) + stats::dnorm(g_ok, mean = mu, sd = sd_tot_2, log = TRUE)
                ll3 <- log(wt3) + stats::dnorm(g_ok, mean = mu, sd = sd_tot_3, log = TRUE)
                ll4 <- log(wt4) + stats::dnorm(g_ok, mean = mu, sd = sd_tot_4, log = TRUE)

                mmax <- pmax(ll1, ll2, ll3, ll4)
                lse <- mmax + log(exp(ll1 - mmax) + exp(ll2 - mmax) + exp(ll3 - mmax) + exp(ll4 - mmax))
                cost[idx_ok] <- cost[idx_ok] - lse
            } else {
                cost[idx_ok] <- cost[idx_ok] +
                    (g_ok - mu)^2 / (2 * sigma_d^2) +
                    log(sigma_d) +
                    0.5 * log(2 * pi)
            }

            # Soft penalty for shrinkage
            if (is.finite(k_shrink) && k_shrink > 0) {
                shrink <- d1_ok < d0
                if (any(shrink)) {
                    dd <- d0 - d1_ok[shrink]
                    cost[idx_ok[shrink]] <- cost[idx_ok[shrink]] + k_shrink * (dd^2)
                }
            }

            # Soft penalty for extreme positive growth
            if (is.finite(max_growth_soft) && is.finite(k_growth) && k_growth > 0) {
                d1_soft_cap <- d0 + max_growth_soft * interval_years
                if (is.finite(d1_soft_cap)) {
                    exceed <- d1_ok > d1_soft_cap
                    if (any(exceed)) {
                        dd <- d1_ok[exceed] - d1_soft_cap
                        cost[idx_ok[exceed]] <- cost[idx_ok[exceed]] + k_growth * (dd^2)
                    }
                }
            }
        }
    }

    # ---------------------------------------------------------------------
    # Deterministic tie-break (non-biological)
    # ---------------------------------------------------------------------
    if (eps_tiebreak > 0) {
        r0 <- rank(track_dbh_t, ties.method = "first")
        for (i in seq_len(n_batch)) {
            row <- mat_tp1[i, ]
            both_obs <- !is.na(track_dbh_t) & !is.na(row)
            if (any(both_obs)) {
                r1 <- rank(row, ties.method = "first")
                cost[i] <- cost[i] + eps_tiebreak * sum(abs(r0[both_obs] - r1[both_obs]))
            }
        }
    }

    cost
}

transition_cost_tracks_bio_components <- function(
  track_dbh_t,
  track_dbh_tp1,
  interval_years,
  # -----------------------
  # GROWTH MODEL PARAMETERS
  # -----------------------
  mu_const,
  mu_gamma = 0,
  sigma0,
  sigma1,
  max_shrink,
  k_shrink,
  max_growth = Inf,
  max_growth_soft = Inf,
  k_growth = 0,
  # -------------------------
  # MORTALITY MODEL PARAMETERS
  # -------------------------
  h0,
  beta,
  # ----------------------------
  # RECRUITMENT MODEL PARAMETERS
  # ----------------------------
  recruit_meanlog,
  recruit_sdlog,
  recruit_max_dbh,
  recruit_lambda,
  # -----------------
  # MEASUREMENT ERROR (optional)
  # -----------------
  use_measurement_error = FALSE,
  meas_sd1_a = 0.0062,
  meas_sd1_b = 0.0904,
  meas_sd2 = 4.64,
  meas_p_big = 0.05,
  # -----------------
  # DETERMINISTIC TIE-BREAK
  # -----------------
  eps_tiebreak = 1e-6,
  hard_penalty = 1e6
) {
    # PURPOSE
    # - Diagnostics companion to transition_cost_tracks_bio().
    # - Returns a per-track breakdown of the *same* terms used in the scalar cost,
    #   plus the tie-break contribution.
    #
    # OUTPUT
    # - list(per_track=..., tiebreak=..., total=..., p_recruit=...)
    #   where per_track is a data.table with case labels and component costs.

    if (length(track_dbh_tp1) != length(track_dbh_t)) {
        stop("track_dbh_t and track_dbh_tp1 must have the same length.", call. = FALSE)
    }

    K <- length(track_dbh_t)
    interval_years <- as.numeric(interval_years)
    if (!is.finite(interval_years) || interval_years <= 0) {
        stop("interval_years must be positive.", call. = FALSE)
    }

    # Recruitment probability over interval
    p_recruit <- 1 - exp(-recruit_lambda * interval_years)
    p_recruit <- pmin(pmax(p_recruit, 1e-12), 1 - 1e-12)

    # Pre-allocate breakdown table
    if (!requireNamespace("data.table", quietly = TRUE)) {
        stop("Package 'data.table' is required for transition_cost_tracks_bio_components().")
    }

    dt <- data.table::data.table(
        track = seq_len(K),
        d0 = as.numeric(track_dbh_t),
        d1 = as.numeric(track_dbh_tp1),
        case = NA_character_,
        # Components (all non-negative; some may be 0)
        cost_recruit = 0,
        cost_no_recruit = 0,
        cost_mortality = 0,
        cost_growth_lik = 0,
        cost_shrink_soft = 0,
        cost_growth_soft = 0,
        cost_hard = 0,
        stringsAsFactors = FALSE
    )

    mu_growth <- function(d) {
        if (!is.finite(mu_gamma) || mu_gamma == 0 || !is.finite(d) || d <= 0) {
            return(mu_const)
        }
        mu_const + mu_gamma * log(d)
    }

    log_sum_exp <- function(x) {
        m <- max(x)
        if (!is.finite(m)) {
            return(m)
        }
        m + log(sum(exp(x - m)))
    }

    meas_sd1 <- function(d) pmax(meas_sd1_a * d + meas_sd1_b, 1e-6)

    for (k in seq_len(K)) {
        d0 <- dt$d0[k]
        d1 <- dt$d1[k]

        # CASE 1: NA -> NA
        if (is.na(d0) && is.na(d1)) {
            dt$case[k] <- "NA->NA"
            dt$cost_no_recruit[k] <- -log(1 - p_recruit)
            next
        }

        # CASE 2: NA -> DBH (recruitment)
        if (is.na(d0) && !is.na(d1)) {
            dt$case[k] <- "NA->DBH"
            if (d1 > recruit_max_dbh) {
                dt$cost_hard[k] <- hard_penalty
            } else {
                dt$cost_recruit[k] <- -log(p_recruit) - dlnorm(d1, recruit_meanlog, recruit_sdlog, log = TRUE)
            }
            next
        }

        # CASE 3: DBH -> NA (mortality)
        if (!is.na(d0) && is.na(d1)) {
            dt$case[k] <- "DBH->NA"
            hazard <- h0 * exp(beta * d0)
            p_death <- 1 - exp(-hazard * interval_years)
            p_death <- pmin(pmax(p_death, 1e-12), 1 - 1e-12)
            dt$cost_mortality[k] <- -log(p_death)
            next
        }

        # CASE 4: DBH -> DBH (growth)
        dt$case[k] <- "DBH->DBH"
        g <- (d1 - d0) / interval_years

        # Hard constraint on shrinkage
        if (is.finite(max_shrink) && (g < max_shrink)) {
            dt$cost_hard[k] <- hard_penalty
            next
        }

        # Hard constraint on extreme positive growth
        if (is.finite(max_growth) && (g > max_growth)) {
            dt$cost_hard[k] <- hard_penalty
            next
        }

        sigma_d <- sigma0 + sigma1 * d0
        sigma_d <- pmax(sigma_d, 1e-6)
        mu <- mu_growth(d0)

        if (isTRUE(use_measurement_error)) {
            s_small0 <- meas_sd1(d0)
            s_small1 <- meas_sd1(d1)
            s_big <- meas_sd2
            w_small <- 1 - meas_p_big
            w_big <- meas_p_big

            sd_meas_mix <- c(
                sqrt(s_small0^2 + s_small1^2) / interval_years,
                sqrt(s_small0^2 + s_big^2) / interval_years,
                sqrt(s_big^2 + s_small1^2) / interval_years,
                sqrt(s_big^2 + s_big^2) / interval_years
            )
            wt_meas_mix <- c(w_small * w_small, w_small * w_big, w_big * w_small, w_big * w_big)
            sd_tot <- sqrt(sigma_d^2 + sd_meas_mix^2)

            ll <- log(wt_meas_mix) + stats::dnorm(g, mean = mu, sd = sd_tot, log = TRUE)
            dt$cost_growth_lik[k] <- -log_sum_exp(ll)
        } else {
            dt$cost_growth_lik[k] <-
                (g - mu)^2 / (2 * sigma_d^2) +
                log(sigma_d) +
                0.5 * log(2 * pi)
        }

        if (d1 < d0) {
            dt$cost_shrink_soft[k] <- k_shrink * (d0 - d1)^2
        }

        if (is.finite(max_growth_soft) && is.finite(k_growth) && k_growth > 0) {
            d1_soft_cap <- d0 + max_growth_soft * interval_years
            if (is.finite(d1_soft_cap) && d1 > d1_soft_cap) {
                dt$cost_growth_soft[k] <- k_growth * (d1 - d1_soft_cap)^2
            }
        }
    }

    tie_cost <- 0
    if (eps_tiebreak > 0) {
        r0 <- rank(track_dbh_t, ties.method = "first")
        r1 <- rank(track_dbh_tp1, ties.method = "first")
        both_obs <- !is.na(track_dbh_t) & !is.na(track_dbh_tp1)
        if (any(both_obs)) {
            tie_cost <- eps_tiebreak * sum(abs(r0[both_obs] - r1[both_obs]))
        }
    }

    dt[, total_track := cost_recruit + cost_no_recruit + cost_mortality + cost_growth_lik + cost_shrink_soft + cost_growth_soft + cost_hard]

    list(
        per_track = dt,
        tiebreak = tie_cost,
        total = sum(dt$total_track) + tie_cost,
        p_recruit = p_recruit
    )
}

# ## validation via simulation
# # Simulate transitions from the model
# simulate_transition <- function(n = 1000) {
#     d0 <- runif(n, 10, 60)
#     # Survival
#     hazard <- Bio_H0_unit * exp(Bio_Beta_unit * d0)
#     alive <- runif(n) > (1 - exp(-hazard))
#     d1 <- rep(NA_real_, n)
#     # Growth
#     g <- rnorm(n, Bio_Mu_Growth_unit, Bio_Sigma_Growth_unit)
#     d1[alive] <- d0[alive] + g[alive]
#     list(d0 = d0, d1 = d1)
# }
# sim <- simulate_transition()

# # Gold-standard test (the one reviewers love)
# cost_true <- transition_cost_tracks_bio(sim$d0, sim$d1, 1)
# cost_random <- replicate(1000, {
#     transition_cost_tracks_bio(sim$d0, sample(sim$d1), 1)
# })
# mean(cost_random > cost_true)
# # Likelihood ratio test

match_stems_dp_global_backward <- function(tree_data,
                                           min_growth,
                                           max_growth,
                                           interval_years = NULL,
                                           anchor_start,
                                           max_tracks,
                                           max_states,
                                           slack_tracks = 1L,
                                           # Optional: measurement error (Condit remeasurement study)
                                           use_measurement_error = TRUE,
                                           meas_sd1_a = 0.0062,
                                           meas_sd1_b = 0.0904,
                                           meas_sd2 = 4.64,
                                           meas_p_big = 0.05,
                                           verbose = FALSE) {
    # PURPOSE
    # - Main global dynamic-programming (DP) stem-ID reconstruction solver.
    # - Run per (Tag, species) group: it reconstructs `ReconstructedStemID` for
    #   censuses 1..anchor_start by choosing the lowest-cost assignment path.
    #
    # REQUIRED INPUT DATA
    # - `tree_data` should contain ONE Tag (and ideally one species) with columns:
    #   - CensusID (1..), DBH (cm; NA means not observed), TrueStemID (known at anchor)
    #   - species (used only for grouping / bookkeeping)
    #   - Bio_* columns used by transition_cost_tracks_bio() (growth/mortality/
    #     recruitment/shrinkage parameters).
    #
    # KEY IDEAS
    # - Tracks: K latent identity slots; each census maps observed stems injectively
    #   onto tracks (a "state").
    # - DP chooses a globally consistent sequence of states across censuses.
    # - Life-cycle constraint is enforced via an extra per-track phase (prebirth/
    #   alive/dead) carried in the DP key to forbid resurrection (OBS->NA->OBS).
    #
    # PARAMETERS / GUARDRAILS
    # - max_tracks: hard cap for K.
    # - max_states: hard cap for per-census injective state enumeration.
    # - slack_tracks: optionally add extra unused tracks to represent simultaneous
    #   death+birth in constant-count transitions.
    # - min_growth/max_growth: used by the stepwise fallback + post-hoc diagnostics
    #   (they are not the DP transition model itself).
    # - use_measurement_error and meas_*: passed into growth likelihood inside
    #   transition_cost_tracks_bio().
    #
    # OUTPUT
    # - Returns `tree_data` with:
    #   - ReconstructedStemID filled
    #   - ReconstructionMethod marked as:
    #       "given"  (TrueStemID rows)
    #       "dp"     (filled by DP decoding)
    #       "igraph"  (filled by fallback when DP is refused)
    #   - DP_KUsed / DP_MaxStatesPerCensus / DP_MaxStatesCensusID diagnostics.
    #
    # FALLBACK
    # - If anchoring fails, state-space is too large, or DP becomes inconsistent,
    #   falls back to match_stems_optimal_backward() (fast, stepwise, not global).

    # tree_data <- xraw[Tag == 1 & species == "sp1", ]

    tree_data <- tree_data[order(CensusID)]
    interval_years <- resolve_interval_years(tree_data, interval_years = interval_years)

    verbose <- isTRUE(verbose) || isTRUE(getOption("dp_global_biol.verbose", FALSE))
    vcat <- function(...) {
        if (!isTRUE(verbose)) {
            return(invisible(NULL))
        }
        cat(..., "\n")
        flush.console()
        invisible(NULL)
    }
    tic <- function() as.numeric(proc.time()[[3L]])
    t_start <- tic()

    tag_val <- tryCatch(
        {
            u <- unique(tree_data$Tag)
            u <- u[!is.na(u)]
            if (length(u) == 1L) u[[1L]] else NA
        },
        error = function(e) NA
    )
    sp_val <- tryCatch(
        {
            u <- unique(tree_data$species)
            u <- u[!is.na(u) & nzchar(as.character(u))]
            if (length(u) == 1L) as.character(u[[1L]]) else NA_character_
        },
        error = function(e) NA_character_
    )
    prefix <- paste0(
        "[DP_GLOBAL",
        if (!is.na(tag_val)) paste0(" Tag=", tag_val) else "",
        if (!is.na(sp_val)) paste0(" species=", sp_val) else "",
        "] "
    )
    vcat(prefix, "Starting MAP DP (Viterbi): anchor_start=", anchor_start, ", interval_years=", interval_years)
    if (!("ReconstructionMethod" %in% names(tree_data))) {
        tree_data[, ReconstructionMethod := NA_character_]
    }
    if (!("DP_MaxStatesPerCensus" %in% names(tree_data))) {
        tree_data[, `:=`(
            DP_MaxStatesPerCensus = NA_real_,
            DP_MaxStatesCensusID = NA_integer_,
            DP_KUsed = NA_integer_
        )]
    }
    tree_data[!is.na(TrueStemID), `:=`(
        ReconstructedStemID = as.integer(TrueStemID),
        ReconstructionMethod = "given"
    )]

    # Precompute observed-stem counts per census up to the anchor.
    # (These drive the DP state-space size.)
    obs_counts <- vapply(
        seq_len(anchor_start),
        function(cc) nrow(tree_data[CensusID == cc & !is.na(DBH)]),
        integer(1L)
    )
    max_obs <- if (length(obs_counts) > 0L) max(obs_counts) else 0L

    # We always preserve any provided TrueStemID (they act as hard anchors).
    # The solver only fills in missing IDs in earlier censuses.

    anchor_obs <- tree_data[CensusID == anchor_start & !is.na(DBH)]
    # If the anchor census has no observed stems or missing TrueStemID, DP cannot anchor.
    # Return a safe stepwise method.
    if (nrow(anchor_obs) == 0L || any(is.na(anchor_obs$TrueStemID))) {
        vcat(prefix, "Cannot anchor DP (missing anchor observations or TrueStemID). Falling back to igraph.")
        # Can't anchor DP. Still report the implied worst-case DP state count using
        # an observed-only K (capped).
        K_used <- as.integer(min(max_obs, max_tracks))
        n_states_by_census <- vapply(obs_counts, function(n_obs) count_injective_states(K_used, n_obs), numeric(1L))
        tree_data[, `:=`(
            DP_KUsed = K_used,
            DP_MaxStatesPerCensus = max(n_states_by_census, na.rm = TRUE),
            DP_MaxStatesCensusID = as.integer(which.max(n_states_by_census))
        )]
        return(match_stems_optimal_backward(tree_data, min_growth, max_growth, interval_years, anchor_start))
    }

    anchor_ids <- sort(unique(anchor_obs$TrueStemID))
    anchor_ids <- anchor_ids[!is.na(anchor_ids)]
    if (length(anchor_ids) == 0L) {
        vcat(prefix, "Cannot anchor DP (no anchor IDs). Falling back to igraph.")
        K_used <- as.integer(min(max_obs, max_tracks))
        n_states_by_census <- vapply(obs_counts, function(n_obs) count_injective_states(K_used, n_obs), numeric(1L))
        tree_data[, `:=`(
            DP_KUsed = K_used,
            DP_MaxStatesPerCensus = max(n_states_by_census, na.rm = TRUE),
            DP_MaxStatesCensusID = as.integer(which.max(n_states_by_census))
        )]
        return(match_stems_optimal_backward(tree_data, min_growth, max_growth, interval_years, anchor_start))
    }

    # Determine K = number of tracks.
    #
    # With the life-cycle constraint (one birth per track; no resurrection), K must
    # be large enough to represent the number of *unique* stems that can appear over
    # time, not just the maximum number observed in a single census.
    #
    # Minimal unique-stem lower bound from observed counts:
    #   unique >= obs_counts[1] + sum_{t} max(0, obs_counts[t+1] - obs_counts[t])
    # (i.e., initial alive stems + total births implied by count increases).
    births_needed <- if (length(obs_counts) >= 2L) sum(pmax(0L, diff(obs_counts))) else 0L
    K_from_counts <- as.integer(if (length(obs_counts) > 0L) obs_counts[1L] + births_needed else 0L)

    # K_base is the minimum number of tracks needed to represent the anchor stems and
    # any births implied by count increases.
    K_base <- max(length(anchor_ids), max_obs, K_from_counts)

    # IMPORTANT (fix for constant-count turnover):
    # If K == max_obs, then in the densest censuses every track is occupied, so the DP
    # cannot represent a "death" (obs->NA) plus a "birth" (NA->obs) in the same
    # transition when the observed count stays constant.
    # In that case, the DP is forced to do obs->obs relabeling (ID swaps).
    # Adding one slack track (K = max_obs + 1) allows simultaneous death+birth.
    slack_tracks <- suppressWarnings(as.integer(slack_tracks))
    if (!is.finite(slack_tracks) || is.na(slack_tracks) || slack_tracks < 0L) {
        slack_tracks <- 0L
    }
    K_target <- K_base
    if (slack_tracks > 0L && K_base == max_obs) {
        K_target <- K_base + slack_tracks
    }
    K <- min(K_target, max_tracks)

    # Report the true worst-case theoretical number of DP states in any census.
    n_states_by_census <- vapply(obs_counts, function(n_obs) count_injective_states(K, n_obs), numeric(1L))
    tree_data[, `:=`(
        DP_KUsed = as.integer(K),
        DP_MaxStatesPerCensus = max(n_states_by_census, na.rm = TRUE),
        DP_MaxStatesCensusID = as.integer(which.max(n_states_by_census))
    )]
    if (K < max(obs_counts)) {
        vcat(prefix, "K too small for observed counts (K=", K, ", max_obs=", max(obs_counts), "). Falling back to igraph.")
        return(match_stems_optimal_backward(tree_data, min_growth, max_growth, interval_years, anchor_start))
    }

    vcat(prefix, "Chosen K=", K, " tracks; max theoretical states=", format(max(n_states_by_census, na.rm = TRUE), scientific = TRUE))

    n_extra <- K - length(anchor_ids)
    current_max <- suppressWarnings(max(tree_data$TrueStemID, na.rm = TRUE))
    if (!is.finite(current_max)) current_max <- 0
    # Track IDs correspond to the numeric IDs we will write into ReconstructedStemID.
    # - Anchor IDs are preserved.
    # - Extra tracks get fresh IDs above the current maximum TrueStemID.
    track_ids <- c(anchor_ids, if (n_extra > 0L) seq.int(from = current_max + 1L, length.out = n_extra) else integer(0))

    # track_ids is the numeric label of each track.
    # When DP assigns an observation to track k, we write track_ids[k] into
    # ReconstructedStemID.

    growth_penalty <- 1e4

    obs_dbh <- vector("list", anchor_start)
    obs_species <- vector("list", anchor_start)
    state_mats <- vector("list", anchor_start)
    state_keys <- vector("list", anchor_start)

    for (cc in seq_len(anchor_start)) {
        # Enumerate all states for census cc.
        # If enumeration is too large (returns NULL), fall back.
        obs <- tree_data[CensusID == cc & !is.na(DBH)]
        obs_dbh[[cc]] <- obs$DBH
        # Species is treated as constant within a (Tag, species) group, but we still
        # store it per-observation for potential future extensions.
        obs_species[[cc]] <- if (nrow(obs) > 0L) as.character(obs$species) else character(0)
        n_obs <- length(obs_dbh[[cc]])

        # At census cc, we consider only stems with a measured DBH.
        # `mat` contains every possible injective assignment of those observations to tracks.
        mat <- enumerate_states_injective(K, n_obs, max_states = max_states)
        if (is.null(mat)) {
            vcat(prefix, "State enumeration exceeded max_states at CensusID=", cc, " (n_obs=", n_obs, "). Falling back to igraph.")
            return(match_stems_optimal_backward(tree_data, min_growth, max_growth, interval_years, anchor_start))
        }
        state_mats[[cc]] <- mat

        # We need a stable key string so we can store DP values in a named vector.
        # Example: state c(2,5,1) becomes the key "2,5,1".
        state_keys[[cc]] <- apply(mat, 1L, state_key)

        vcat(prefix, "Enumerated CensusID=", cc, ": n_obs=", n_obs, ", n_states=", nrow(mat))
    }

    anchor_obs_ordered <- tree_data[CensusID == anchor_start & !is.na(DBH)]
    # Anchor state: map each anchor observation to the track with the same TrueStemID.
    # This pins the DP endpoint.
    anchor_track_idx <- match(anchor_obs_ordered$TrueStemID, track_ids)
    if (any(is.na(anchor_track_idx))) {
        vcat(prefix, "Anchor TrueStemID not found in track_ids. Falling back to igraph.")
        return(match_stems_optimal_backward(tree_data, min_growth, max_growth, interval_years, anchor_start))
    }
    anchor_state_key <- state_key(anchor_track_idx)

    # ---------------------------------------------------------------------
    # Life-cycle constraint (one birth per track; no resurrection)
    # ---------------------------------------------------------------------
    # A track's observation status over time must be contiguous:
    #   NA*  -> (birth) -> OBS+ -> (death) -> NA*
    # i.e., we forbid OBS -> NA -> OBS for the same track.
    #
    # We implement this by carrying a per-track phase in the DP state:
    #   0 = prebirth (not yet born), 1 = alive (observed), 2 = dead (cannot return)
    #
    # DP keys become: "<assignment_key>|<phase0,phase1,...,phaseK>".
    # The assignment part is exactly the old state_key() output.

    encode_full_key <- function(assign_vec, phase_vec) {
        paste0(state_key(assign_vec), "|", paste(as.integer(phase_vec), collapse = ","))
    }

    decode_full_key <- function(k) {
        parts <- strsplit(k, "|", fixed = TRUE)[[1L]]
        a <- parts[[1L]]
        p <- if (length(parts) >= 2L) parts[[2L]] else ""
        assign_vec <- if (a == "") integer(0) else as.integer(strsplit(a, ",", fixed = TRUE)[[1L]])
        phase_vec <- if (p == "") integer(0) else as.integer(strsplit(p, ",", fixed = TRUE)[[1L]])
        list(assign = assign_vec, phase = phase_vec)
    }

    derive_phase_prev <- function(phase_tp1, tdbh_t, tdbh_tp1) {
        # Given phase at t+1 and observation status at t and t+1, derive phase at t.
        # Returns integer(K) phase vector or NULL if inconsistent.
        K_loc <- length(tdbh_t)
        if (length(phase_tp1) != K_loc) {
            return(NULL)
        }
        alive_t <- !is.na(tdbh_t)
        alive_tp1 <- !is.na(tdbh_tp1)

        # Phase must match observed status at t+1.
        if (any(alive_tp1 & phase_tp1 != 1L)) {
            return(NULL)
        }
        if (any((!alive_tp1) & phase_tp1 == 1L)) {
            return(NULL)
        }

        phase_t <- rep.int(NA_integer_, K_loc)
        for (k in seq_len(K_loc)) {
            if (alive_tp1[k]) {
                # t+1 alive: at t either alive (alive->alive) or prebirth (birth).
                if (alive_t[k]) {
                    phase_t[k] <- 1L
                } else {
                    phase_t[k] <- 0L
                }
            } else {
                # t+1 not alive: either still prebirth, or dead.
                if (phase_tp1[k] == 0L) {
                    # Still prebirth: must also be unobserved at t.
                    if (alive_t[k]) {
                        return(NULL)
                    }
                    phase_t[k] <- 0L
                } else if (phase_tp1[k] == 2L) {
                    # Dead: at t either alive (death) or dead.
                    if (alive_t[k]) {
                        phase_t[k] <- 1L
                    } else {
                        phase_t[k] <- 2L
                    }
                } else {
                    return(NULL)
                }
            }
        }

        # Sanity: alive <-> phase==1
        if (any(alive_t & phase_t != 1L)) {
            return(NULL)
        }
        if (any((!alive_t) & phase_t == 1L)) {
            return(NULL)
        }
        phase_t
    }

    # DP tables (now phase-aware):
    # - dp_next: named vector mapping full_key at census cc+1 -> best future cost
    # - backptr[[cc]]: named vector mapping full_key at census cc -> chosen next full_key
    #
    # Initialize at anchor census with cost 0 for the fixed anchor state.
    phase_anchor <- rep.int(2L, K)
    phase_anchor[anchor_track_idx] <- 1L
    anchor_full_key <- encode_full_key(anchor_track_idx, phase_anchor)
    dp_next <- setNames(0, anchor_full_key)
    backptr <- vector("list", anchor_start)

    ## extract bio params
    # Bio_Mu_Growth_unit <- unique(tree_data$Bio_Mu_Growth)
    # Bio_Sigma_Growth_unit <- unique(tree_data$Bio_Sigma_Growth)
    # Bio_H0_unit <- unique(tree_data$Bio_H0)
    # Bio_Beta_unit <- unique(tree_data$Bio_Beta)
    # Bio_Recruit_Meanlog_unit <- unique(tree_data$Bio_Recruit_Meanlog)
    # Bio_Recruit_Sdlog_unit <- unique(tree_data$Bio_Recruit_Sdlog)
    # Bio_Recruit_MaxDBH_unit <- unique(tree_data$Bio_Recruit_MaxDBH_unit)

    ## extract bio params
    Bio_Mu_Growth_unit <- unique(tree_data$Bio_Mu_Growth)
    Bio_Gamma_Growth_unit <- if ("Bio_Gamma_Growth" %in% names(tree_data)) {
        unique(tree_data$Bio_Gamma_Growth)
    } else {
        0
    }
    # Bio_Sigma_Growth_unit <- unique(tree_data$Bio_Sigma_Growth)
    Bio_Sigma0_unit <- unique(tree_data$Bio_Sigma0_Growth)
    Bio_Sigma1_unit <- unique(tree_data$Bio_Sigma1_Growth)
    Bio_max_shrink_unit <- unique(tree_data$Bio_Max_Shrink)
    Bio_k_shrink_unit <- unique(tree_data$Bio_K_Shrink)
    Bio_max_growth_unit <- if ("Bio_Max_Growth" %in% names(tree_data)) {
        unique(tree_data$Bio_Max_Growth)
    } else {
        Inf
    }
    Bio_max_growth_soft_unit <- if ("Bio_Max_Growth_Soft" %in% names(tree_data)) {
        unique(tree_data$Bio_Max_Growth_Soft)
    } else {
        Inf
    }
    Bio_k_growth_unit <- if ("Bio_K_Growth" %in% names(tree_data)) {
        unique(tree_data$Bio_K_Growth)
    } else {
        0
    }
    Bio_max_growth_unit <- if ("Bio_Max_Growth" %in% names(tree_data)) {
        unique(tree_data$Bio_Max_Growth)
    } else {
        Inf
    }
    Bio_max_growth_soft_unit <- if ("Bio_Max_Growth_Soft" %in% names(tree_data)) {
        unique(tree_data$Bio_Max_Growth_Soft)
    } else {
        Inf
    }
    Bio_k_growth_unit <- if ("Bio_K_Growth" %in% names(tree_data)) {
        unique(tree_data$Bio_K_Growth)
    } else {
        0
    }
    # Mortality parameter column names vary across datasets/scripts.
    # Use only Bio_H0_Mortality and Bio_Beta_Mortality for mortality parameters.
    Bio_H0_Mortality <- unique(tree_data$Bio_H0_Mortality)
    Bio_Beta_Mortality <- unique(tree_data$Bio_Beta_Mortality)
    Bio_Recruit_Meanlog_unit <- unique(tree_data$Bio_Recruit_Meanlog)
    Bio_Recruit_Sdlog_unit <- unique(tree_data$Bio_Recruit_Sdlog)
    Bio_Recruit_MaxDBH_unit <- unique(tree_data$Bio_Recruit_MaxDBH_unit)
    Bio_Recruitment_lambda <- unique(tree_data$Bio_Recruitment_lambda)
    # eps_tie <- 1e-06

    vcat(prefix, "DP recursion (backward) starting ...")
    for (cc in seq.int(anchor_start - 1L, 1L, by = -1L)) {
        # cc <- 6
        # DP recursion step (backward in time):
        # For each state at census cc, choose the best next state at census cc+1.

        # DP bookkeeping
        # - dp_next contains the best achievable cost from census (cc+1) onward,
        #   keyed by the chosen (assignment, phase) at census (cc+1).
        # - For each candidate state at cc, we try all reachable next states and
        #   keep the transition that minimizes (transition_cost + dp_next).
        mat_cc <- state_mats[[cc]]
        n_states_cc <- nrow(mat_cc)
        dp_curr <- numeric(0)
        ptr_curr <- character(0)

        # Precompute track-wise DBH/species vectors for each state at census cc.
        # This makes the transition cost computation simpler.
        track_dbh_cc <- vector("list", n_states_cc)
        for (i in seq_len(n_states_cc)) {
            # Convert a state (a mapping obs->track) into track DBHs.
            # This turns the state into something we can score track-by-track.
            track_dbh_cc[[i]] <- state_to_track_dbh(mat_cc[i, ], obs_dbh[[cc]], K)
        }

        keys_next <- names(dp_next)
        keys_next <- keys_next[!is.na(keys_next)]
        if (length(keys_next) == 0L) {
            vcat(prefix, "No reachable next states at CensusID=", cc + 1L, ". Falling back to igraph.")
            return(match_stems_optimal_backward(tree_data, min_growth, max_growth, interval_years, anchor_start))
        }

        t_cc0 <- tic()
        vcat(
            prefix, "Backward step CensusID=", cc, ": n_states=", n_states_cc, ", n_next=", length(keys_next),
            ", approx transitions=", format(n_states_cc * length(keys_next), scientific = TRUE)
        )

        for (i in seq_len(n_states_cc)) {
            # i <- 1
            tdbh0 <- track_dbh_cc[[i]]

            # Try every reachable next full_key.
            for (kn in keys_next) {
                # kn <- keys_next[1L]
                dec <- decode_full_key(kn)
                s1 <- dec$assign
                phase_tp1 <- dec$phase
                if (length(phase_tp1) != K) next
                tdbh1 <- state_to_track_dbh(s1, obs_dbh[[cc + 1L]], K)
                fut_cost <- dp_next[kn]
                if (!is.finite(fut_cost)) next

                phase_t <- derive_phase_prev(phase_tp1, tdbh0, tdbh1)
                if (is.null(phase_t)) next

                curr_full_key <- encode_full_key(mat_cc[i, ], phase_t)

                # cst <- transition_cost_tracks_bio(
                #     track_dbh_t = tdbh0,
                #     track_dbh_tp1 = tdbh1,
                #     interval_years = interval_years,
                #     # --- growth model ---
                #     mu_const = Bio_Mu_Growth_unit,
                #     sigma_growth = Bio_Sigma_Growth_unit,
                #     # --- mortality model ---
                #     h0 = Bio_H0_Mortality,
                #     beta = Bio_Beta_Mortality,
                # max_growth = Bio_max_growth_unit,
                # max_growth_soft = Bio_max_growth_soft_unit,
                # k_growth = Bio_k_growth_unit,
                #     # --- recruitment model ---
                #     recruit_meanlog = Bio_Recruit_Meanlog_unit,
                #     recruit_sdlog = Bio_Recruit_Sdlog_unit,
                #     eps_tiebreak = 1e-6
                # ) + fut_cost

                cst <- transition_cost_tracks_bio(
                    track_dbh_t = tdbh0,
                    track_dbh_tp1 = tdbh1,
                    interval_years = interval_years,
                    # --- growth model ---
                    mu_const = Bio_Mu_Growth_unit,
                    mu_gamma = Bio_Gamma_Growth_unit,
                    # sigma_growth = Bio_Sigma_Growth_unit,
                    sigma0 = Bio_Sigma0_unit,
                    sigma1 = Bio_Sigma1_unit,
                    max_shrink = Bio_max_shrink_unit,
                    k_shrink = Bio_k_shrink_unit,
                    max_growth = Bio_max_growth_unit,
                    max_growth_soft = Bio_max_growth_soft_unit,
                    k_growth = Bio_k_growth_unit,
                    # --- measurement error (optional) ---
                    use_measurement_error = use_measurement_error,
                    meas_sd1_a = meas_sd1_a,
                    meas_sd1_b = meas_sd1_b,
                    meas_sd2 = meas_sd2,
                    meas_p_big = meas_p_big,
                    # --- mortality model ---
                    h0 = Bio_H0_Mortality,
                    beta = Bio_Beta_Mortality,
                    # --- recruitment model ---
                    recruit_meanlog = Bio_Recruit_Meanlog_unit,
                    recruit_sdlog = Bio_Recruit_Sdlog_unit,
                    recruit_max_dbh = Bio_Recruit_MaxDBH_unit,
                    recruit_lambda = Bio_Recruitment_lambda,
                    eps_tiebreak = 1e-6
                ) + fut_cost

                prev_best <- dp_curr[curr_full_key]
                prev_best_val <- if (length(prev_best) == 0L) NA_real_ else as.numeric(prev_best)
                if (!is.finite(prev_best_val) || cst < prev_best_val) {
                    dp_curr[curr_full_key] <- cst
                    ptr_curr[curr_full_key] <- kn
                }
            }
        }

        if (length(dp_curr) == 0L) {
            vcat(prefix, "DP produced no states at CensusID=", cc, ". Falling back to igraph.")
            return(match_stems_optimal_backward(tree_data, min_growth, max_growth, interval_years, anchor_start))
        }

        vcat(prefix, "Finished CensusID=", cc, ": kept ", length(dp_curr), " states; dt=", sprintf("%.2fs", tic() - t_cc0))

        backptr[[cc]] <- ptr_curr
        dp_next <- dp_curr
    }

    if (all(!is.finite(dp_next))) {
        vcat(prefix, "All start-state costs non-finite. Falling back to igraph.")
        return(match_stems_optimal_backward(tree_data, min_growth, max_growth, interval_years, anchor_start))
    }
    # Choose the best starting state at census 1 (minimum total cost).
    start_key <- names(dp_next)[which.min(dp_next)]

    # Decode best path of state keys from census 1 to anchor.
    chosen_keys <- rep("", anchor_start)
    chosen_keys[1L] <- start_key
    for (cc in seq_len(anchor_start - 1L)) {
        # Follow the backpointer chain:
        # - chosen_keys[cc] tells us the chosen state at census cc.
        # - backptr[[cc]][[chosen_keys[cc]]] tells us which state to use at cc+1.
        chosen_keys[cc + 1L] <- backptr[[cc]][[chosen_keys[cc]]]
    }
    chosen_keys[anchor_start] <- anchor_full_key

    vcat(prefix, "Decoding MAP path and writing ReconstructedStemID ...")
    for (cc in seq_len(anchor_start)) {
        # Convert the chosen state's track indices into actual IDs, one census at a time.
        # obs_idx are the row indices in tree_data where DBH is observed.
        obs_idx <- tree_data[CensusID == cc & !is.na(DBH), which = TRUE]
        if (length(obs_idx) == 0L) next

        # chosen_keys[cc] is a string like "2,5,1".
        # With the life-cycle constraint enabled, chosen_keys[cc] is "<assign>|<phase>".
        dec <- decode_full_key(chosen_keys[cc])
        sv <- dec$assign
        if (length(sv) != length(obs_idx)) {
            vcat(prefix, "Path decode inconsistency at CensusID=", cc, ". Falling back to igraph.")
            return(match_stems_optimal_backward(tree_data, min_growth, max_growth, interval_years, anchor_start))
        }

        # Now the key idea:
        # - sv[j] is the track index used by the j-th observed stem at this census.
        # - track_ids[sv[j]] is the real numeric ID we want to write.
        tree_data[obs_idx, ReconstructedStemID := track_ids[sv]]
        obs_to_mark <- obs_idx[is.na(tree_data$TrueStemID[obs_idx])]
        if (length(obs_to_mark) > 0L) {
            tree_data[obs_to_mark, ReconstructionMethod := "dp"]
        }
    }

    tree_data <- add_constraint_violation(
        tree_data,
        id_col = "ReconstructedStemID",
        min_growth = min_growth,
        max_growth = max_growth,
        interval_years = interval_years
    )

    vcat(prefix, "Done. Total elapsed ", sprintf("%.2fs", tic() - t_start))
    return(tree_data)
}

# ---------------------------------------------------------------------
# Global DP WITH TRUE MARGINAL PROBABILITIES (sum-product)
# ---------------------------------------------------------------------
# This is a probabilistic counterpart to match_stems_dp_global_backward().
#
# Key difference:
# - The original DP is Viterbi/MAP: it keeps only the minimum-cost path.
# - This function computes *exact marginal probabilities* under the model
#   P(path) ∝ exp(- total_cost / temperature)
#   by running a backward log-sum-exp recursion and a forward recursion.
#
# Important:
# - It reuses the expensive transition-cost evaluations: during the backward pass
#   we compute transition_cost_tracks_bio() once per (curr_state, next_state)
#   and store the resulting log-weights for the forward pass.
# - If the algorithm falls back to igraph (state space too large / cannot anchor),
#   marginals are returned as NA.

match_stems_dp_global_backward_marginals <- function(tree_data,
                                                     min_growth = -Inf,
                                                     max_growth = Inf,
                                                     interval_years = NULL,
                                                     anchor_start,
                                                     max_tracks = 30L,
                                                     slack_tracks = 1L,
                                                     max_states = 50000L,
                                                     temperature = 1.0,
                                                     posterior_top_k = 2L,
                                                     eps_tiebreak = 1e-6,
                                                     # --- measurement error (optional) ---
                                                     use_measurement_error = FALSE,
                                                     meas_sd1_a = 0.0062,
                                                     meas_sd1_b = 0.0904,
                                                     meas_sd2 = 4.64,
                                                     meas_p_big = 0.05,
                                                     verbose = FALSE) {
    # Safety
    posterior_top_k <- as.integer(posterior_top_k)
    if (!is.finite(posterior_top_k) || is.na(posterior_top_k) || posterior_top_k < 1L) {
        posterior_top_k <- 1L
    }
    temperature <- suppressWarnings(as.numeric(temperature))
    if (!is.finite(temperature) || is.na(temperature) || temperature <= 0) {
        stop("temperature must be a positive finite number")
    }

    verbose <- isTRUE(verbose) || isTRUE(getOption("dp_global_biol.verbose", FALSE))
    vcat <- function(...) {
        if (!isTRUE(verbose)) {
            return(invisible(NULL))
        }
        cat(..., "\n")
        flush.console()
        invisible(NULL)
    }
    tic <- function() as.numeric(proc.time()[[3L]])
    t_start <- tic()

    tag_val <- tryCatch(
        {
            u <- unique(tree_data$Tag)
            u <- u[!is.na(u)]
            if (length(u) == 1L) u[[1L]] else NA
        },
        error = function(e) NA
    )
    sp_val <- tryCatch(
        {
            u <- unique(tree_data$species)
            u <- u[!is.na(u) & nzchar(as.character(u))]
            if (length(u) == 1L) as.character(u[[1L]]) else NA_character_
        },
        error = function(e) NA_character_
    )
    prefix <- paste0(
        "[DP_GLOBAL",
        if (!is.na(tag_val)) paste0(" Tag=", tag_val) else "",
        if (!is.na(sp_val)) paste0(" species=", sp_val) else "",
        "] "
    )

    # Ensure deterministic ordering + resolve interval before logging.
    tree_data <- tree_data[order(CensusID)]
    interval_years <- resolve_interval_years(tree_data, interval_years = interval_years)

    vcat(
        prefix, "Starting marginal DP (sum-product): anchor_start=", anchor_start, ", interval_years=", interval_years,
        ", temperature=", temperature, ", top_k=", posterior_top_k
    )

    # Helpers
    log_add_exp <- function(a, b) {
        # stable log(exp(a) + exp(b)) for scalars
        if (!is.finite(a)) {
            return(b)
        }
        if (!is.finite(b)) {
            return(a)
        }
        m <- max(a, b)
        m + log(exp(a - m) + exp(b - m))
    }
    log_sum_exp <- function(x) {
        x <- x[is.finite(x)]
        if (length(x) == 0L) {
            return(-Inf)
        }
        m <- max(x)
        m + log(sum(exp(x - m)))
    }

    # Ensure columns exist
    if (!("ReconstructionMethod" %in% names(tree_data))) {
        tree_data[, ReconstructionMethod := NA_character_]
    }
    if (!("DP_MaxStatesPerCensus" %in% names(tree_data))) {
        tree_data[, `:=`(
            DP_MaxStatesPerCensus = NA_real_,
            DP_MaxStatesCensusID = NA_integer_,
            DP_KUsed = NA_integer_
        )]
    }
    # Posterior output columns
    post_cols <- c(
        "DP_PosteriorTop1ID", "DP_PosteriorTop1Prob",
        "DP_PosteriorTop2ID", "DP_PosteriorTop2Prob",
        "DP_PosteriorEntropy", "DP_PosteriorReconstructedProb",
        "DP_PosteriorUnlinkedProb"
    )

    ensure_posterior_columns <- function(dt) {
        # data.table uses the type of the RHS to infer the column type.
        # Bare NA creates a logical column, which later triggers warnings when
        # assigning integer/numeric posteriors.
        if (!("DP_PosteriorTop1ID" %in% names(dt))) dt[, DP_PosteriorTop1ID := NA_integer_]
        if (!("DP_PosteriorTop2ID" %in% names(dt))) dt[, DP_PosteriorTop2ID := NA_integer_]
        if (!("DP_PosteriorTop1Prob" %in% names(dt))) dt[, DP_PosteriorTop1Prob := NA_real_]
        if (!("DP_PosteriorTop2Prob" %in% names(dt))) dt[, DP_PosteriorTop2Prob := NA_real_]
        if (!("DP_PosteriorEntropy" %in% names(dt))) dt[, DP_PosteriorEntropy := NA_real_]
        if (!("DP_PosteriorReconstructedProb" %in% names(dt))) dt[, DP_PosteriorReconstructedProb := NA_real_]
        if (!("DP_PosteriorUnlinkedProb" %in% names(dt))) dt[, DP_PosteriorUnlinkedProb := NA_real_]

        # If columns already exist but have the wrong type, coerce once.
        if (!(is.integer(dt$DP_PosteriorTop1ID))) dt[, DP_PosteriorTop1ID := as.integer(DP_PosteriorTop1ID)]
        if (!(is.integer(dt$DP_PosteriorTop2ID))) dt[, DP_PosteriorTop2ID := as.integer(DP_PosteriorTop2ID)]

        if (!(is.numeric(dt$DP_PosteriorTop1Prob))) dt[, DP_PosteriorTop1Prob := as.numeric(DP_PosteriorTop1Prob)]
        if (!(is.numeric(dt$DP_PosteriorTop2Prob))) dt[, DP_PosteriorTop2Prob := as.numeric(DP_PosteriorTop2Prob)]
        if (!(is.numeric(dt$DP_PosteriorEntropy))) dt[, DP_PosteriorEntropy := as.numeric(DP_PosteriorEntropy)]
        if (!(is.numeric(dt$DP_PosteriorReconstructedProb))) dt[, DP_PosteriorReconstructedProb := as.numeric(DP_PosteriorReconstructedProb)]
        if (!(is.numeric(dt$DP_PosteriorUnlinkedProb))) dt[, DP_PosteriorUnlinkedProb := as.numeric(DP_PosteriorUnlinkedProb)]

        dt
    }

    tree_data <- ensure_posterior_columns(tree_data)

    # Preserve any provided TrueStemID as hard values in output.
    tree_data[!is.na(TrueStemID), `:=`(
        ReconstructedStemID = as.integer(TrueStemID),
        ReconstructionMethod = "given"
    )]

    # Observed-stem counts up to anchor (for state-space diagnostics)
    obs_counts <- vapply(
        seq_len(anchor_start),
        function(cc) nrow(tree_data[CensusID == cc & !is.na(DBH)]),
        integer(1L)
    )
    max_obs <- if (length(obs_counts) > 0L) max(obs_counts) else 0L

    # Need a fully-anchored endpoint
    anchor_obs <- tree_data[CensusID == anchor_start & !is.na(DBH)]
    if (nrow(anchor_obs) == 0L || any(is.na(anchor_obs$TrueStemID))) {
        vcat(prefix, "Cannot anchor DP (missing anchor observations or TrueStemID). Falling back to igraph.")
        K_used <- as.integer(min(max_obs, max_tracks))
        n_states_by_census <- vapply(obs_counts, function(n_obs) count_injective_states(K_used, n_obs), numeric(1L))
        tree_data[, `:=`(
            DP_KUsed = K_used,
            DP_MaxStatesPerCensus = max(n_states_by_census, na.rm = TRUE),
            DP_MaxStatesCensusID = as.integer(which.max(n_states_by_census))
        )]
        out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, interval_years, anchor_start)
        out <- ensure_posterior_columns(out)
        return(out)
    }
    anchor_ids <- sort(unique(anchor_obs$TrueStemID))
    anchor_ids <- anchor_ids[!is.na(anchor_ids)]
    if (length(anchor_ids) == 0L) {
        vcat(prefix, "Cannot anchor DP (no anchor IDs). Falling back to igraph.")
        K_used <- as.integer(min(max_obs, max_tracks))
        n_states_by_census <- vapply(obs_counts, function(n_obs) count_injective_states(K_used, n_obs), numeric(1L))
        tree_data[, `:=`(
            DP_KUsed = K_used,
            DP_MaxStatesPerCensus = max(n_states_by_census, na.rm = TRUE),
            DP_MaxStatesCensusID = as.integer(which.max(n_states_by_census))
        )]
        out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, interval_years, anchor_start)
        out <- ensure_posterior_columns(out)
        return(out)
    }

    # Choose K (tracks) same logic as the MAP DP
    births_needed <- if (length(obs_counts) >= 2L) sum(pmax(0L, diff(obs_counts))) else 0L
    K_from_counts <- as.integer(if (length(obs_counts) > 0L) obs_counts[1L] + births_needed else 0L)
    K_base <- max(length(anchor_ids), max_obs, K_from_counts)

    slack_tracks <- suppressWarnings(as.integer(slack_tracks))
    if (!is.finite(slack_tracks) || is.na(slack_tracks) || slack_tracks < 0L) {
        slack_tracks <- 0L
    }
    K_target <- K_base
    if (slack_tracks > 0L && K_base == max_obs) {
        K_target <- K_base + slack_tracks
    }
    K <- min(K_target, max_tracks)

    # Report theoretical worst-case state count
    n_states_by_census <- vapply(obs_counts, function(n_obs) count_injective_states(K, n_obs), numeric(1L))
    tree_data[, `:=`(
        DP_KUsed = as.integer(K),
        DP_MaxStatesPerCensus = max(n_states_by_census, na.rm = TRUE),
        DP_MaxStatesCensusID = as.integer(which.max(n_states_by_census))
    )]
    if (K < max(obs_counts)) {
        vcat(prefix, "K too small for observed counts (K=", K, ", max_obs=", max(obs_counts), "). Falling back to igraph.")
        out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, interval_years, anchor_start)
        out <- ensure_posterior_columns(out)
        return(out)
    }

    vcat(prefix, "Chosen K=", K, " tracks; max theoretical states=", format(max(n_states_by_census, na.rm = TRUE), scientific = TRUE))

    n_extra <- K - length(anchor_ids)
    current_max <- suppressWarnings(max(tree_data$TrueStemID, na.rm = TRUE))
    if (!is.finite(current_max)) current_max <- 0
    track_ids <- c(anchor_ids, if (n_extra > 0L) seq.int(from = current_max + 1L, length.out = n_extra) else integer(0))

    # Pre-enumerate assignment states (injective obs->track) for each census
    obs_dbh <- vector("list", anchor_start)
    state_mats <- vector("list", anchor_start)
    state_keys <- vector("list", anchor_start)
    for (cc in seq_len(anchor_start)) {
        obs <- tree_data[CensusID == cc & !is.na(DBH)]
        obs_dbh[[cc]] <- obs$DBH
        n_obs <- length(obs_dbh[[cc]])
        mat <- enumerate_states_injective(K, n_obs, max_states = max_states)
        if (is.null(mat)) {
            vcat(prefix, "State enumeration exceeded max_states at CensusID=", cc, " (n_obs=", n_obs, "). Falling back to igraph.")
            out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, interval_years, anchor_start)
            out <- ensure_posterior_columns(out)
            return(out)
        }
        state_mats[[cc]] <- mat
        state_keys[[cc]] <- apply(mat, 1L, state_key)

        vcat(prefix, "Enumerated CensusID=", cc, ": n_obs=", n_obs, ", n_states=", nrow(mat))
    }

    # Anchor state assignment (pins endpoint)
    anchor_obs_ordered <- tree_data[CensusID == anchor_start & !is.na(DBH)]
    anchor_track_idx <- match(anchor_obs_ordered$TrueStemID, track_ids)
    if (any(is.na(anchor_track_idx))) {
        vcat(prefix, "Anchor TrueStemID not found in track_ids. Falling back to igraph.")
        out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, interval_years, anchor_start)
        out <- ensure_posterior_columns(out)
        return(out)
    }

    # Life-cycle phase constraint helpers (internal representation)
    # Phase values are in {0,1,2}. For speed, we encode the phase vector as a compact
    # byte-string of digits (e.g., "201102...") instead of comma-separated integers.
    encode_phase_key <- function(phase_vec) {
        # phase_vec is integer vector with values 0/1/2
        rawToChar(as.raw(as.integer(phase_vec) + 48L))
    }
    decode_phase_key <- function(p) {
        if (!nzchar(p)) {
            return(integer(0))
        }
        as.integer(utf8ToInt(p) - 48L)
    }
    encode_full_key <- function(assign_vec, phase_vec) {
        paste0(state_key(assign_vec), "|", encode_phase_key(phase_vec))
    }
    decode_full_key <- function(k) {
        parts <- strsplit(k, "|", fixed = TRUE)[[1L]]
        a <- parts[[1L]]
        p <- if (length(parts) >= 2L) parts[[2L]] else ""
        assign_vec <- if (a == "") integer(0) else as.integer(strsplit(a, ",", fixed = TRUE)[[1L]])
        phase_vec <- decode_phase_key(p)
        list(assign = assign_vec, phase = phase_vec)
    }
    derive_phase_prev <- function(phase_tp1, tdbh_t, tdbh_tp1) {
        K_loc <- length(tdbh_t)
        if (length(phase_tp1) != K_loc) {
            return(NULL)
        }
        alive_t <- !is.na(tdbh_t)
        alive_tp1 <- !is.na(tdbh_tp1)
        if (any(alive_tp1 & phase_tp1 != 1L)) {
            return(NULL)
        }
        if (any((!alive_tp1) & phase_tp1 == 1L)) {
            return(NULL)
        }

        phase_t <- integer(K_loc)
        for (k in seq_len(K_loc)) {
            if (alive_tp1[k]) {
                phase_t[k] <- if (alive_t[k]) 1L else 0L
            } else {
                if (phase_tp1[k] == 0L) {
                    if (alive_t[k]) {
                        return(NULL)
                    }
                    phase_t[k] <- 0L
                } else if (phase_tp1[k] == 2L) {
                    phase_t[k] <- if (alive_t[k]) 1L else 2L
                } else {
                    return(NULL)
                }
            }
        }
        if (any(alive_t & phase_t != 1L)) {
            return(NULL)
        }
        if (any((!alive_t) & phase_t == 1L)) {
            return(NULL)
        }
        phase_t
    }

    # Bio params (same extraction logic as MAP DP)
    Bio_Mu_Growth_unit <- unique(tree_data$Bio_Mu_Growth)
    Bio_Gamma_Growth_unit <- if ("Bio_Gamma_Growth" %in% names(tree_data)) {
        unique(tree_data$Bio_Gamma_Growth)
    } else {
        0
    }
    Bio_Sigma0_unit <- unique(tree_data$Bio_Sigma0_Growth)
    Bio_Sigma1_unit <- unique(tree_data$Bio_Sigma1_Growth)
    Bio_max_shrink_unit <- unique(tree_data$Bio_Max_Shrink)
    Bio_k_shrink_unit <- unique(tree_data$Bio_K_Shrink)
    Bio_max_growth_unit <- if ("Bio_Max_Growth" %in% names(tree_data)) {
        unique(tree_data$Bio_Max_Growth)
    } else {
        Inf
    }
    Bio_max_growth_soft_unit <- if ("Bio_Max_Growth_Soft" %in% names(tree_data)) {
        unique(tree_data$Bio_Max_Growth_Soft)
    } else {
        Inf
    }
    Bio_k_growth_unit <- if ("Bio_K_Growth" %in% names(tree_data)) {
        unique(tree_data$Bio_K_Growth)
    } else {
        0
    }
    Bio_H0_Mortality <- unique(tree_data$Bio_H0_Mortality)
    Bio_Beta_Mortality <- unique(tree_data$Bio_Beta_Mortality)
    Bio_Recruit_Meanlog_unit <- unique(tree_data$Bio_Recruit_Meanlog)
    Bio_Recruit_Sdlog_unit <- unique(tree_data$Bio_Recruit_Sdlog)
    Bio_Recruit_MaxDBH_unit <- unique(tree_data$Bio_Recruit_MaxDBH_unit)
    Bio_Recruitment_lambda <- unique(tree_data$Bio_Recruitment_lambda)

    # Precompute track-wise DBH vectors for each *assignment* state (phase doesn't matter for costs)
    track_dbh_by_state <- vector("list", anchor_start)
    for (cc in seq_len(anchor_start)) {
        mat <- state_mats[[cc]]
        n_states <- nrow(mat)
        tdbh_list <- vector("list", n_states)
        for (i in seq_len(n_states)) {
            tdbh_list[[i]] <- state_to_track_dbh(mat[i, ], obs_dbh[[cc]], K)
        }
        track_dbh_by_state[[cc]] <- tdbh_list
    }

    # Backward tables (sum-product) and Viterbi backpointers
    # For each census cc we store:
    # - keys_full[[cc]]: vector of full-state keys (phase-aware)
    # - assign_full[[cc]]: list of assignment vectors aligned with keys_full
    # - logB[[cc]]: backward log-sum weights aligned with keys_full
    # - vit_cost[[cc]]: Viterbi cost-to-go aligned with keys_full
    # - vit_ptr[[cc]]: integer index of chosen next state (into keys_full[[cc+1]])
    # - edges[[cc]]: data.table(from_idx, to_idx, logw) for forward pass

    keys_full <- vector("list", anchor_start)
    assign_full <- vector("list", anchor_start)
    logB <- vector("list", anchor_start)
    vit_cost <- vector("list", anchor_start)
    vit_ptr <- vector("list", anchor_start)
    edges <- vector("list", anchor_start)

    # Initialize at anchor census: only one allowed full-state
    phase_anchor <- rep.int(2L, K)
    phase_anchor[anchor_track_idx] <- 1L
    anchor_full_key <- encode_full_key(anchor_track_idx, phase_anchor)
    keys_full[[anchor_start]] <- anchor_full_key
    assign_full[[anchor_start]] <- list(as.integer(anchor_track_idx))
    logB[[anchor_start]] <- 0
    vit_cost[[anchor_start]] <- 0
    vit_ptr[[anchor_start]] <- integer(0)
    edges[[anchor_start]] <- NULL

    # Backward recursion cc = anchor_start-1 .. 1
    vcat(prefix, "Backward pass (log-sum-exp + Viterbi) starting ...")
    for (cc in seq.int(anchor_start - 1L, 1L, by = -1L)) {
        mat_cc <- state_mats[[cc]]
        n_states_cc <- nrow(mat_cc)

        t_cc0 <- tic()
        vcat(prefix, "Backward step CensusID=", cc, ": n_assignment_states=", n_states_cc, ", n_next_full_states=", length(keys_full[[cc + 1L]]))

        next_keys <- keys_full[[cc + 1L]]
        n_next <- length(next_keys)
        if (n_next == 0L) {
            vcat(prefix, "No reachable next full-states at CensusID=", cc + 1L, ". Falling back to igraph.")
            out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, interval_years, anchor_start)
            out <- ensure_posterior_columns(out)
            return(out)
        }
        next_index <- seq_len(n_next)
        names(next_index) <- next_keys
        logB_next <- as.numeric(logB[[cc + 1L]])
        vit_next <- as.numeric(vit_cost[[cc + 1L]])

        # Fast lookup for next assignment DBHs via assignment key
        next_assign_list <- assign_full[[cc + 1L]]
        # Build assignment-key -> state index for next census (since phase differs but assignment cost uses assignment)
        next_assign_key <- vapply(next_assign_list, state_key, character(1L))
        next_assign_row_idx <- match(next_assign_key, state_keys[[cc + 1L]])
        if (any(is.na(next_assign_row_idx))) {
            # Should not happen; indicates mismatch in state enumeration.
            out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, interval_years, anchor_start)
            out <- ensure_posterior_columns(out)
            return(out)
        }

        # Preallocate edge arrays (worst-case: every (assignment_state, next_full_state) is feasible)
        upper_edges <- n_states_cc * n_next
        from_idx <- integer(upper_edges)
        to_idx <- integer(upper_edges)
        logw <- numeric(upper_edges)
        used_edges <- 0L

        # Dynamic creation of current full-states
        key_to_idx <- new.env(parent = emptyenv())
        curr_keys_list <- list()
        curr_assign_list <- list()
        curr_logB <- numeric(0)
        curr_vit <- numeric(0)
        curr_ptr <- integer(0)

        for (i in seq_len(n_states_cc)) {
            tdbh0 <- track_dbh_by_state[[cc]][[i]]
            assign0 <- mat_cc[i, ]

            for (j in seq_len(n_next)) {
                # Next full state j
                assign1 <- next_assign_list[[j]]
                # track DBH at cc+1 depends only on assignment
                tdbh1 <- track_dbh_by_state[[cc + 1L]][[next_assign_row_idx[j]]]
                phase_tp1 <- decode_full_key(next_keys[j])$phase
                if (length(phase_tp1) != K) next
                phase_t <- derive_phase_prev(phase_tp1, tdbh0, tdbh1)
                if (is.null(phase_t)) next
                curr_key <- encode_full_key(assign0, phase_t)

                idx <- key_to_idx[[curr_key]]
                if (is.null(idx)) {
                    idx <- length(curr_keys_list) + 1L
                    key_to_idx[[curr_key]] <- idx
                    curr_keys_list[[idx]] <- curr_key
                    curr_assign_list[[idx]] <- as.integer(assign0)
                    curr_logB[idx] <- -Inf
                    curr_vit[idx] <- Inf
                    curr_ptr[idx] <- NA_integer_
                }

                # Transition cost computed ONCE; reused for both Viterbi and marginals
                # print a message that shows the percentage missing to compute
                # vcat(prefix, "Computing transition cost for CensusID=", cc, ", curr_state=", i, "/", n_states_cc, ", next_state=", j, "/", n_next, " ...")
                c_trans <- transition_cost_tracks_bio(
                    track_dbh_t = tdbh0,
                    track_dbh_tp1 = tdbh1,
                    interval_years = interval_years,
                    # --- growth model ---
                    mu_const = Bio_Mu_Growth_unit,
                    mu_gamma = Bio_Gamma_Growth_unit,
                    sigma0 = Bio_Sigma0_unit,
                    sigma1 = Bio_Sigma1_unit,
                    max_shrink = Bio_max_shrink_unit,
                    k_shrink = Bio_k_shrink_unit,
                    max_growth = Bio_max_growth_unit,
                    max_growth_soft = Bio_max_growth_soft_unit,
                    k_growth = Bio_k_growth_unit,
                    # --- measurement error (optional) ---
                    use_measurement_error = use_measurement_error,
                    meas_sd1_a = meas_sd1_a,
                    meas_sd1_b = meas_sd1_b,
                    meas_sd2 = meas_sd2,
                    meas_p_big = meas_p_big,
                    # --- mortality model ---
                    h0 = Bio_H0_Mortality,
                    beta = Bio_Beta_Mortality,
                    # --- recruitment model ---
                    recruit_meanlog = Bio_Recruit_Meanlog_unit,
                    recruit_sdlog = Bio_Recruit_Sdlog_unit,
                    recruit_max_dbh = Bio_Recruit_MaxDBH_unit,
                    recruit_lambda = Bio_Recruitment_lambda,
                    eps_tiebreak = eps_tiebreak
                )

                # Viterbi update (MAP)
                cand_vit <- c_trans + vit_next[j]
                if (!is.finite(curr_vit[idx]) || cand_vit < curr_vit[idx]) {
                    curr_vit[idx] <- cand_vit
                    curr_ptr[idx] <- j
                }

                # Backward log-sum-exp update
                cand_log <- (-c_trans / temperature) + logB_next[j]
                curr_logB[idx] <- log_add_exp(curr_logB[idx], cand_log)

                # Store edge for forward pass
                used_edges <- used_edges + 1L
                from_idx[used_edges] <- idx
                to_idx[used_edges] <- j
                logw[used_edges] <- (-c_trans / temperature)
            }
        }

        if (length(curr_keys_list) == 0L) {
            vcat(prefix, "DP produced no states at CensusID=", cc, ". Falling back to igraph.")
            out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, interval_years, anchor_start)
            out <- ensure_posterior_columns(out)
            return(out)
        }

        keys_full[[cc]] <- unlist(curr_keys_list, use.names = FALSE)
        assign_full[[cc]] <- curr_assign_list
        logB[[cc]] <- curr_logB
        vit_cost[[cc]] <- curr_vit
        vit_ptr[[cc]] <- curr_ptr

        if (used_edges == 0L) {
            vcat(prefix, "No feasible edges at CensusID=", cc, ". Falling back to igraph.")
            out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, interval_years, anchor_start)
            out <- ensure_posterior_columns(out)
            return(out)
        }
        edges[[cc]] <- data.table::data.table(
            from_idx = from_idx[seq_len(used_edges)],
            to_idx = to_idx[seq_len(used_edges)],
            logw = logw[seq_len(used_edges)]
        )

        vcat(prefix, "Finished CensusID=", cc, ": full_states=", length(keys_full[[cc]]), ", edges=", used_edges, ", dt=", sprintf("%.2fs", tic() - t_cc0))
    }

    # -----------------
    # Decode MAP path
    # -----------------
    vcat(prefix, "Decoding MAP path and writing ReconstructedStemID ...")
    start_idx <- which.min(vit_cost[[1L]])
    if (length(start_idx) == 0L || !is.finite(vit_cost[[1L]][start_idx])) {
        out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, interval_years, anchor_start)
        out <- ensure_posterior_columns(out)
        return(out)
    }
    map_idx <- integer(anchor_start)
    map_idx[1L] <- start_idx
    for (cc in seq_len(anchor_start - 1L)) {
        nxt <- vit_ptr[[cc]][map_idx[cc]]
        if (!is.finite(nxt) || is.na(nxt) || nxt < 1L) {
            out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, interval_years, anchor_start)
            out <- ensure_posterior_columns(out)
            return(out)
        }
        map_idx[cc + 1L] <- nxt
    }

    for (cc in seq_len(anchor_start)) {
        obs_idx <- tree_data[CensusID == cc & !is.na(DBH), which = TRUE]
        if (length(obs_idx) == 0L) next
        sv <- assign_full[[cc]][[map_idx[cc]]]
        if (length(sv) != length(obs_idx)) {
            out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, interval_years, anchor_start)
            out <- ensure_posterior_columns(out)
            return(out)
        }
        tree_data[obs_idx, ReconstructedStemID := track_ids[sv]]
        obs_to_mark <- obs_idx[is.na(tree_data$TrueStemID[obs_idx])]
        if (length(obs_to_mark) > 0L) {
            tree_data[obs_to_mark, ReconstructionMethod := "dp"]
        }
    }

    # -----------------
    # Forward pass for marginals
    # -----------------
    vcat(prefix, "Forward pass starting ...")
    # Start distribution: uniform over all reachable states at census 1.
    logalpha <- vector("list", anchor_start)
    logalpha[[1L]] <- rep.int(0, length(keys_full[[1L]]))

    for (cc in seq_len(anchor_start - 1L)) {
        ed <- edges[[cc]]
        if (is.null(ed) || nrow(ed) == 0L) {
            out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, interval_years, anchor_start)
            out <- ensure_posterior_columns(out)
            return(out)
        }
        la_from <- logalpha[[cc]][ed$from_idx]
        vals <- la_from + ed$logw
        dt <- data.table::data.table(to_idx = ed$to_idx, v = vals)
        dt <- dt[is.finite(v)]
        if (nrow(dt) == 0L) {
            out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, interval_years, anchor_start)
            out <- ensure_posterior_columns(out)
            return(out)
        }
        la_next_dt <- dt[, .(logalpha = log_sum_exp(v)), by = to_idx]
        la_next <- rep.int(-Inf, length(keys_full[[cc + 1L]]))
        la_next[la_next_dt$to_idx] <- la_next_dt$logalpha
        logalpha[[cc + 1L]] <- la_next

        vcat(prefix, "Forward step CensusID=", cc + 1L, ": reached ", sum(is.finite(la_next)), " / ", length(la_next), " states")
    }

    # Partition function Z = total weight of all paths ending at the fixed anchor state.
    # At anchor_start there is exactly one state.
    logZ <- logalpha[[anchor_start]][1L]
    if (!is.finite(logZ)) {
        # Fallback: compute from backward at census 1.
        logZ <- log_sum_exp(logB[[1L]])
    }

    vcat(prefix, "Computed logZ=", sprintf("%.3f", logZ))

    # -----------------
    # Observation-level marginals
    # -----------------
    anchor_set <- anchor_ids
    is_anchor_track <- track_ids %in% anchor_set

    for (cc in seq_len(anchor_start)) {
        obs_idx <- tree_data[CensusID == cc & !is.na(DBH), which = TRUE]
        if (length(obs_idx) == 0L) next

        # State posterior weights at this census
        lg <- logalpha[[cc]] + logB[[cc]] - logZ
        # Normalize defensively
        lg <- lg - log_sum_exp(lg)
        w <- exp(lg)

        n_obs <- length(obs_idx)
        prob_mat <- matrix(0, nrow = n_obs, ncol = K)
        st_assign <- assign_full[[cc]]
        for (s in seq_along(st_assign)) {
            ww <- w[s]
            if (!is.finite(ww) || ww <= 0) next
            a <- st_assign[[s]]
            # a is length n_obs, entries in 1..K
            for (j in seq_len(n_obs)) {
                prob_mat[j, a[j]] <- prob_mat[j, a[j]] + ww
            }
        }
        # Re-normalize each row (floating error)
        row_sums <- rowSums(prob_mat)
        for (j in seq_len(n_obs)) {
            if (row_sums[j] > 0) {
                prob_mat[j, ] <- prob_mat[j, ] / row_sums[j]
            }
        }

        # MAP assignment for DP_PosteriorReconstructedProb
        map_assign <- assign_full[[cc]][[map_idx[cc]]]

        top1_id <- integer(n_obs)
        top2_id <- integer(n_obs)
        top1_p <- numeric(n_obs)
        top2_p <- numeric(n_obs)
        entropy <- numeric(n_obs)
        p_map <- numeric(n_obs)
        p_unlinked <- numeric(n_obs)

        for (j in seq_len(n_obs)) {
            p <- prob_mat[j, ]
            # numeric stability
            p[p < 0] <- 0
            sp <- sum(p)
            if (sp > 0) p <- p / sp
            ord <- order(p, decreasing = TRUE)
            i1 <- ord[1L]
            top1_id[j] <- track_ids[i1]
            top1_p[j] <- p[i1]
            if (posterior_top_k >= 2L && length(ord) >= 2L) {
                i2 <- ord[2L]
                top2_id[j] <- track_ids[i2]
                top2_p[j] <- p[i2]
            } else {
                top2_id[j] <- NA_integer_
                top2_p[j] <- NA_real_
            }
            entropy[j] <- -sum(ifelse(p > 0, p * log(p), 0))
            p_map[j] <- p[map_assign[j]]
            p_unlinked[j] <- sum(p[!is_anchor_track])
        }

        tree_data[obs_idx, `:=`(
            DP_PosteriorTop1ID = top1_id,
            DP_PosteriorTop1Prob = top1_p,
            DP_PosteriorTop2ID = top2_id,
            DP_PosteriorTop2Prob = top2_p,
            DP_PosteriorEntropy = entropy,
            DP_PosteriorReconstructedProb = p_map,
            DP_PosteriorUnlinkedProb = p_unlinked
        )]
    }

    tree_data <- add_constraint_violation(
        tree_data,
        id_col = "ReconstructedStemID",
        min_growth = min_growth,
        max_growth = max_growth,
        interval_years = interval_years
    )

    vcat(prefix, "Done. Total elapsed ", sprintf("%.2fs", tic() - t_start))
    return(tree_data)
}

match_stems_dp_global_backward_marginals_batch <- function(tree_data,
                                                           min_growth = -Inf,
                                                           max_growth = Inf,
                                                           interval_years = NULL,
                                                           anchor_start,
                                                           max_tracks = 30L,
                                                           slack_tracks = 1L,
                                                           max_states = 50000L,
                                                           temperature = 1.0,
                                                           posterior_top_k = 2L,
                                                           eps_tiebreak = 1e-6,
                                                           # --- measurement error (optional) ---
                                                           use_measurement_error = FALSE,
                                                           meas_sd1_a = 0.0062,
                                                           meas_sd1_b = 0.0904,
                                                           meas_sd2 = 4.64,
                                                           meas_p_big = 0.05,
                                                           verbose = FALSE) {
    # Safety
    posterior_top_k <- as.integer(posterior_top_k)
    if (!is.finite(posterior_top_k) || is.na(posterior_top_k) || posterior_top_k < 1L) {
        posterior_top_k <- 1L
    }
    temperature <- suppressWarnings(as.numeric(temperature))
    if (!is.finite(temperature) || is.na(temperature) || temperature <= 0) {
        stop("temperature must be a positive finite number")
    }

    verbose <- isTRUE(verbose) || isTRUE(getOption("dp_global_biol.verbose", FALSE))
    vcat <- function(...) {
        if (!isTRUE(verbose)) {
            return(invisible(NULL))
        }
        cat(..., "\n")
        flush.console()
        invisible(NULL)
    }
    tic <- function() as.numeric(proc.time()[[3L]])
    t_start <- tic()

    tag_val <- tryCatch(
        {
            u <- unique(tree_data$Tag)
            u <- u[!is.na(u)]
            if (length(u) == 1L) u[[1L]] else NA
        },
        error = function(e) NA
    )
    sp_val <- tryCatch(
        {
            u <- unique(tree_data$species)
            u <- u[!is.na(u) & nzchar(as.character(u))]
            if (length(u) == 1L) as.character(u[[1L]]) else NA_character_
        },
        error = function(e) NA_character_
    )
    prefix <- paste0(
        "[dp_global_batch",
        if (!is.na(tag_val)) paste0(" Tag=", tag_val) else "",
        if (!is.na(sp_val)) paste0(" species=", sp_val) else "",
        "] "
    )

    # Ensure deterministic ordering + resolve interval before logging.
    tree_data <- tree_data[order(CensusID)]
    interval_years <- resolve_interval_years(tree_data, interval_years = interval_years)

    vcat(
        prefix, "Starting marginal DP (sum-product) [batch costs]: anchor_start=", anchor_start, ", interval_years=", interval_years,
        ", temperature=", temperature, ", top_k=", posterior_top_k
    )

    # Helpers
    log_add_exp <- function(a, b) {
        # stable log(exp(a) + exp(b)) for scalars
        if (!is.finite(a)) {
            return(b)
        }
        if (!is.finite(b)) {
            return(a)
        }
        m <- max(a, b)
        m + log(exp(a - m) + exp(b - m))
    }
    log_sum_exp <- function(x) {
        x <- x[is.finite(x)]
        if (length(x) == 0L) {
            return(-Inf)
        }
        m <- max(x)
        m + log(sum(exp(x - m)))
    }

    # Ensure columns exist
    if (!("ReconstructionMethod" %in% names(tree_data))) {
        tree_data[, ReconstructionMethod := NA_character_]
    }
    if (!("DP_MaxStatesPerCensus" %in% names(tree_data))) {
        tree_data[, `:=`(
            DP_MaxStatesPerCensus = NA_real_,
            DP_MaxStatesCensusID = NA_integer_,
            DP_KUsed = NA_integer_
        )]
    }
    # Posterior output columns
    post_cols <- c(
        "DP_PosteriorTop1ID", "DP_PosteriorTop1Prob",
        "DP_PosteriorTop2ID", "DP_PosteriorTop2Prob",
        "DP_PosteriorEntropy", "DP_PosteriorReconstructedProb",
        "DP_PosteriorUnlinkedProb"
    )

    ensure_posterior_columns <- function(dt) {
        # data.table uses the type of the RHS to infer the column type.
        # Bare NA creates a logical column, which later triggers warnings when
        # assigning integer/numeric posteriors.
        if (!("DP_PosteriorTop1ID" %in% names(dt))) dt[, DP_PosteriorTop1ID := NA_integer_]
        if (!("DP_PosteriorTop2ID" %in% names(dt))) dt[, DP_PosteriorTop2ID := NA_integer_]
        if (!("DP_PosteriorTop1Prob" %in% names(dt))) dt[, DP_PosteriorTop1Prob := NA_real_]
        if (!("DP_PosteriorTop2Prob" %in% names(dt))) dt[, DP_PosteriorTop2Prob := NA_real_]
        if (!("DP_PosteriorEntropy" %in% names(dt))) dt[, DP_PosteriorEntropy := NA_real_]
        if (!("DP_PosteriorReconstructedProb" %in% names(dt))) dt[, DP_PosteriorReconstructedProb := NA_real_]
        if (!("DP_PosteriorUnlinkedProb" %in% names(dt))) dt[, DP_PosteriorUnlinkedProb := NA_real_]

        # If columns already exist but have the wrong type, coerce once.
        if (!(is.integer(dt$DP_PosteriorTop1ID))) dt[, DP_PosteriorTop1ID := as.integer(DP_PosteriorTop1ID)]
        if (!(is.integer(dt$DP_PosteriorTop2ID))) dt[, DP_PosteriorTop2ID := as.integer(DP_PosteriorTop2ID)]

        if (!(is.numeric(dt$DP_PosteriorTop1Prob))) dt[, DP_PosteriorTop1Prob := as.numeric(DP_PosteriorTop1Prob)]
        if (!(is.numeric(dt$DP_PosteriorTop2Prob))) dt[, DP_PosteriorTop2Prob := as.numeric(DP_PosteriorTop2Prob)]
        if (!(is.numeric(dt$DP_PosteriorEntropy))) dt[, DP_PosteriorEntropy := as.numeric(DP_PosteriorEntropy)]
        if (!(is.numeric(dt$DP_PosteriorReconstructedProb))) dt[, DP_PosteriorReconstructedProb := as.numeric(DP_PosteriorReconstructedProb)]
        if (!(is.numeric(dt$DP_PosteriorUnlinkedProb))) dt[, DP_PosteriorUnlinkedProb := as.numeric(DP_PosteriorUnlinkedProb)]

        dt
    }

    tree_data <- ensure_posterior_columns(tree_data)

    # Preserve any provided TrueStemID as hard values in output.
    tree_data[!is.na(TrueStemID), `:=`(
        ReconstructedStemID = as.integer(TrueStemID),
        ReconstructionMethod = "given"
    )]

    # Observed-stem counts up to anchor (for state-space diagnostics)
    obs_counts <- vapply(
        seq_len(anchor_start),
        function(cc) nrow(tree_data[CensusID == cc & !is.na(DBH)]),
        integer(1L)
    )
    max_obs <- if (length(obs_counts) > 0L) max(obs_counts) else 0L

    # Need a fully-anchored endpoint
    anchor_obs <- tree_data[CensusID == anchor_start & !is.na(DBH)]
    if (nrow(anchor_obs) == 0L || any(is.na(anchor_obs$TrueStemID))) {
        vcat(prefix, "Cannot anchor DP (missing anchor observations or TrueStemID). Falling back to igraph.")
        K_used <- as.integer(min(max_obs, max_tracks))
        n_states_by_census <- vapply(obs_counts, function(n_obs) count_injective_states(K_used, n_obs), numeric(1L))
        tree_data[, `:=`(
            DP_KUsed = K_used,
            DP_MaxStatesPerCensus = max(n_states_by_census, na.rm = TRUE),
            DP_MaxStatesCensusID = as.integer(which.max(n_states_by_census))
        )]
        out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, interval_years, anchor_start)
        out <- ensure_posterior_columns(out)
        return(out)
    }
    anchor_ids <- sort(unique(anchor_obs$TrueStemID))
    anchor_ids <- anchor_ids[!is.na(anchor_ids)]
    if (length(anchor_ids) == 0L) {
        vcat(prefix, "Cannot anchor DP (no anchor IDs). Falling back to igraph.")
        K_used <- as.integer(min(max_obs, max_tracks))
        n_states_by_census <- vapply(obs_counts, function(n_obs) count_injective_states(K_used, n_obs), numeric(1L))
        tree_data[, `:=`(
            DP_KUsed = K_used,
            DP_MaxStatesPerCensus = max(n_states_by_census, na.rm = TRUE),
            DP_MaxStatesCensusID = as.integer(which.max(n_states_by_census))
        )]
        out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, interval_years, anchor_start)
        out <- ensure_posterior_columns(out)
        return(out)
    }

    # Choose K (tracks) same logic as the MAP DP
    births_needed <- if (length(obs_counts) >= 2L) sum(pmax(0L, diff(obs_counts))) else 0L
    K_from_counts <- as.integer(if (length(obs_counts) > 0L) obs_counts[1L] + births_needed else 0L)
    K_base <- max(length(anchor_ids), max_obs, K_from_counts)

    slack_tracks <- suppressWarnings(as.integer(slack_tracks))
    if (!is.finite(slack_tracks) || is.na(slack_tracks) || slack_tracks < 0L) {
        slack_tracks <- 0L
    }
    K_target <- K_base
    if (slack_tracks > 0L && K_base == max_obs) {
        K_target <- K_base + slack_tracks
    }
    K <- min(K_target, max_tracks)

    # Report theoretical worst-case state count
    n_states_by_census <- vapply(obs_counts, function(n_obs) count_injective_states(K, n_obs), numeric(1L))
    tree_data[, `:=`(
        DP_KUsed = as.integer(K),
        DP_MaxStatesPerCensus = max(n_states_by_census, na.rm = TRUE),
        DP_MaxStatesCensusID = as.integer(which.max(n_states_by_census))
    )]
    if (K < max(obs_counts)) {
        vcat(prefix, "K too small for observed counts (K=", K, ", max_obs=", max(obs_counts), "). Falling back to igraph.")
        out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, interval_years, anchor_start)
        out <- ensure_posterior_columns(out)
        return(out)
    }

    vcat(prefix, "Chosen K=", K, " tracks; max theoretical states=", format(max(n_states_by_census, na.rm = TRUE), scientific = TRUE))

    n_extra <- K - length(anchor_ids)
    current_max <- suppressWarnings(max(tree_data$TrueStemID, na.rm = TRUE))
    if (!is.finite(current_max)) current_max <- 0
    track_ids <- c(anchor_ids, if (n_extra > 0L) seq.int(from = current_max + 1L, length.out = n_extra) else integer(0))

    # Pre-enumerate assignment states (injective obs->track) for each census
    obs_dbh <- vector("list", anchor_start)
    state_mats <- vector("list", anchor_start)
    state_keys <- vector("list", anchor_start)
    for (cc in seq_len(anchor_start)) {
        obs <- tree_data[CensusID == cc & !is.na(DBH)]
        obs_dbh[[cc]] <- obs$DBH
        n_obs <- length(obs_dbh[[cc]])
        mat <- enumerate_states_injective(K, n_obs, max_states = max_states)
        if (is.null(mat)) {
            vcat(prefix, "State enumeration exceeded max_states at CensusID=", cc, " (n_obs=", n_obs, "). Falling back to igraph.")
            out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, interval_years, anchor_start)
            out <- ensure_posterior_columns(out)
            return(out)
        }
        state_mats[[cc]] <- mat
        state_keys[[cc]] <- apply(mat, 1L, state_key)

        vcat(prefix, "Enumerated CensusID=", cc, ": n_obs=", n_obs, ", n_states=", nrow(mat))
    }

    # Anchor state assignment (pins endpoint)
    anchor_obs_ordered <- tree_data[CensusID == anchor_start & !is.na(DBH)]
    anchor_track_idx <- match(anchor_obs_ordered$TrueStemID, track_ids)
    if (any(is.na(anchor_track_idx))) {
        vcat(prefix, "Anchor TrueStemID not found in track_ids. Falling back to igraph.")
        out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, interval_years, anchor_start)
        out <- ensure_posterior_columns(out)
        return(out)
    }

    # Life-cycle phase constraint helpers (internal representation)
    # Phase values are in {0,1,2}. For speed, we encode the phase vector as a compact
    # byte-string of digits (e.g., "201102...") instead of comma-separated integers.
    encode_phase_key <- function(phase_vec) {
        # phase_vec is integer vector with values 0/1/2
        rawToChar(as.raw(as.integer(phase_vec) + 48L))
    }
    decode_phase_key <- function(p) {
        if (!nzchar(p)) {
            return(integer(0))
        }
        as.integer(utf8ToInt(p) - 48L)
    }
    encode_full_key <- function(assign_vec, phase_vec) {
        paste0(state_key(assign_vec), "|", encode_phase_key(phase_vec))
    }
    decode_full_key <- function(k) {
        parts <- strsplit(k, "|", fixed = TRUE)[[1L]]
        a <- parts[[1L]]
        p <- if (length(parts) >= 2L) parts[[2L]] else ""
        assign_vec <- if (a == "") integer(0) else as.integer(strsplit(a, ",", fixed = TRUE)[[1L]])
        phase_vec <- decode_phase_key(p)
        list(assign = assign_vec, phase = phase_vec)
    }
    derive_phase_prev <- function(phase_tp1, tdbh_t, tdbh_tp1) {
        K_loc <- length(tdbh_t)
        if (length(phase_tp1) != K_loc) {
            return(NULL)
        }
        alive_t <- !is.na(tdbh_t)
        alive_tp1 <- !is.na(tdbh_tp1)
        if (any(alive_tp1 & phase_tp1 != 1L)) {
            return(NULL)
        }
        if (any((!alive_tp1) & phase_tp1 == 1L)) {
            return(NULL)
        }

        phase_t <- integer(K_loc)
        for (k in seq_len(K_loc)) {
            if (alive_tp1[k]) {
                phase_t[k] <- if (alive_t[k]) 1L else 0L
            } else {
                if (phase_tp1[k] == 0L) {
                    if (alive_t[k]) {
                        return(NULL)
                    }
                    phase_t[k] <- 0L
                } else if (phase_tp1[k] == 2L) {
                    phase_t[k] <- if (alive_t[k]) 1L else 2L
                } else {
                    return(NULL)
                }
            }
        }
        if (any(alive_t & phase_t != 1L)) {
            return(NULL)
        }
        if (any((!alive_t) & phase_t == 1L)) {
            return(NULL)
        }
        phase_t
    }

    # Bio params (same extraction logic as MAP DP)
    Bio_Mu_Growth_unit <- unique(tree_data$Bio_Mu_Growth)
    Bio_Gamma_Growth_unit <- if ("Bio_Gamma_Growth" %in% names(tree_data)) {
        unique(tree_data$Bio_Gamma_Growth)
    } else {
        0
    }
    Bio_Sigma0_unit <- unique(tree_data$Bio_Sigma0_Growth)
    Bio_Sigma1_unit <- unique(tree_data$Bio_Sigma1_Growth)
    Bio_max_shrink_unit <- unique(tree_data$Bio_Max_Shrink)
    Bio_k_shrink_unit <- unique(tree_data$Bio_K_Shrink)
    Bio_max_growth_unit <- if ("Bio_Max_Growth" %in% names(tree_data)) {
        unique(tree_data$Bio_Max_Growth)
    } else {
        Inf
    }
    Bio_max_growth_soft_unit <- if ("Bio_Max_Growth_Soft" %in% names(tree_data)) {
        unique(tree_data$Bio_Max_Growth_Soft)
    } else {
        Inf
    }
    Bio_k_growth_unit <- if ("Bio_K_Growth" %in% names(tree_data)) {
        unique(tree_data$Bio_K_Growth)
    } else {
        0
    }
    Bio_H0_Mortality <- unique(tree_data$Bio_H0_Mortality)
    Bio_Beta_Mortality <- unique(tree_data$Bio_Beta_Mortality)
    Bio_Recruit_Meanlog_unit <- unique(tree_data$Bio_Recruit_Meanlog)
    Bio_Recruit_Sdlog_unit <- unique(tree_data$Bio_Recruit_Sdlog)
    Bio_Recruit_MaxDBH_unit <- unique(tree_data$Bio_Recruit_MaxDBH_unit)
    Bio_Recruitment_lambda <- unique(tree_data$Bio_Recruitment_lambda)

    # Precompute track-wise DBH vectors for each *assignment* state (phase doesn't matter for costs)
    track_dbh_by_state <- vector("list", anchor_start)
    for (cc in seq_len(anchor_start)) {
        mat <- state_mats[[cc]]
        n_states <- nrow(mat)
        tdbh_list <- vector("list", n_states)
        for (i in seq_len(n_states)) {
            tdbh_list[[i]] <- state_to_track_dbh(mat[i, ], obs_dbh[[cc]], K)
        }
        track_dbh_by_state[[cc]] <- tdbh_list
    }

    # Backward tables (sum-product) and Viterbi backpointers
    keys_full <- vector("list", anchor_start)
    assign_full <- vector("list", anchor_start)
    logB <- vector("list", anchor_start)
    vit_cost <- vector("list", anchor_start)
    vit_ptr <- vector("list", anchor_start)
    edges <- vector("list", anchor_start)

    # Initialize at anchor census: only one allowed full-state
    phase_anchor <- rep.int(2L, K)
    phase_anchor[anchor_track_idx] <- 1L
    anchor_full_key <- encode_full_key(anchor_track_idx, phase_anchor)
    keys_full[[anchor_start]] <- anchor_full_key
    assign_full[[anchor_start]] <- list(as.integer(anchor_track_idx))
    logB[[anchor_start]] <- 0
    vit_cost[[anchor_start]] <- 0
    vit_ptr[[anchor_start]] <- integer(0)
    edges[[anchor_start]] <- NULL

    # Backward recursion cc = anchor_start-1 .. 1
    vcat(prefix, "Backward pass (log-sum-exp + Viterbi) starting ...")
    for (cc in seq.int(anchor_start - 1L, 1L, by = -1L)) {
        mat_cc <- state_mats[[cc]]
        n_states_cc <- nrow(mat_cc)

        t_cc0 <- tic()
        vcat(prefix, "Backward step CensusID=", cc, ": n_assignment_states=", n_states_cc, ", n_next_full_states=", length(keys_full[[cc + 1L]]))

        next_keys <- keys_full[[cc + 1L]]
        n_next <- length(next_keys)
        if (n_next == 0L) {
            vcat(prefix, "No reachable next full-states at CensusID=", cc + 1L, ". Falling back to igraph.")
            out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, interval_years, anchor_start)
            out <- ensure_posterior_columns(out)
            return(out)
        }
        next_index <- seq_len(n_next)
        names(next_index) <- next_keys
        logB_next <- as.numeric(logB[[cc + 1L]])
        vit_next <- as.numeric(vit_cost[[cc + 1L]])

        # Fast lookup for next assignment DBHs via assignment key
        next_assign_list <- assign_full[[cc + 1L]]
        # Build assignment-key -> state index for next census (since phase differs but assignment cost uses assignment)
        next_assign_key <- vapply(next_assign_list, state_key, character(1L))
        next_assign_row_idx <- match(next_assign_key, state_keys[[cc + 1L]])
        if (any(is.na(next_assign_row_idx))) {
            # Should not happen; indicates mismatch in state enumeration.
            out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, interval_years, anchor_start)
            out <- ensure_posterior_columns(out)
            return(out)
        }

        # Pre-decode next phases and next track-DBH vectors (shared across all current states)
        phase_tp1_by_next <- vector("list", n_next)
        for (j in seq_len(n_next)) {
            phase_tp1_by_next[[j]] <- decode_full_key(next_keys[[j]])$phase
        }
        tdbh1_by_next <- vector("list", n_next)
        for (j in seq_len(n_next)) {
            tdbh1_by_next[[j]] <- track_dbh_by_state[[cc + 1L]][[next_assign_row_idx[[j]]]]
        }

        # Preallocate edge arrays (worst-case: every (assignment_state, next_full_state) is feasible)
        upper_edges <- n_states_cc * n_next
        from_idx <- integer(upper_edges)
        to_idx <- integer(upper_edges)
        logw <- numeric(upper_edges)
        used_edges <- 0L

        # Dynamic creation of current full-states
        key_to_idx <- new.env(parent = emptyenv())
        curr_keys_list <- list()
        curr_assign_list <- list()
        curr_logB <- numeric(0)
        curr_vit <- numeric(0)
        curr_ptr <- integer(0)

        for (i in seq_len(n_states_cc)) {
            tdbh0 <- track_dbh_by_state[[cc]][[i]]
            assign0 <- mat_cc[i, ]

            # Collect all feasible next states (based on phase constraints) for this current assignment
            feasible_n <- 0L
            feasible_j <- integer(n_next)
            feasible_key <- character(n_next)
            feasible_tdbh1 <- vector("list", n_next)

            for (j in seq_len(n_next)) {
                phase_tp1 <- phase_tp1_by_next[[j]]
                if (length(phase_tp1) != K) next
                tdbh1 <- tdbh1_by_next[[j]]
                phase_t <- derive_phase_prev(phase_tp1, tdbh0, tdbh1)
                if (is.null(phase_t)) next

                feasible_n <- feasible_n + 1L
                feasible_j[[feasible_n]] <- j
                feasible_key[[feasible_n]] <- encode_full_key(assign0, phase_t)
                feasible_tdbh1[[feasible_n]] <- tdbh1
            }

            if (feasible_n == 0L) next

            feasible_j <- feasible_j[seq_len(feasible_n)]
            feasible_key <- feasible_key[seq_len(feasible_n)]
            feasible_tdbh1 <- feasible_tdbh1[seq_len(feasible_n)]

            # Batch compute all transition costs from this current assignment
            c_trans_vec <- transition_cost_tracks_bio_batch(
                track_dbh_t = tdbh0,
                track_dbh_tp1 = feasible_tdbh1,
                interval_years = interval_years,
                # --- growth model ---
                mu_const = Bio_Mu_Growth_unit,
                mu_gamma = Bio_Gamma_Growth_unit,
                sigma0 = Bio_Sigma0_unit,
                sigma1 = Bio_Sigma1_unit,
                max_shrink = Bio_max_shrink_unit,
                k_shrink = Bio_k_shrink_unit,
                max_growth = Bio_max_growth_unit,
                max_growth_soft = Bio_max_growth_soft_unit,
                k_growth = Bio_k_growth_unit,
                # --- measurement error (optional) ---
                use_measurement_error = use_measurement_error,
                meas_sd1_a = meas_sd1_a,
                meas_sd1_b = meas_sd1_b,
                meas_sd2 = meas_sd2,
                meas_p_big = meas_p_big,
                # --- mortality model ---
                h0 = Bio_H0_Mortality,
                beta = Bio_Beta_Mortality,
                # --- recruitment model ---
                recruit_meanlog = Bio_Recruit_Meanlog_unit,
                recruit_sdlog = Bio_Recruit_Sdlog_unit,
                recruit_max_dbh = Bio_Recruit_MaxDBH_unit,
                recruit_lambda = Bio_Recruitment_lambda,
                eps_tiebreak = eps_tiebreak
            )

            for (e in seq_len(feasible_n)) {
                j <- feasible_j[[e]]
                curr_key <- feasible_key[[e]]

                idx <- key_to_idx[[curr_key]]
                if (is.null(idx)) {
                    idx <- length(curr_keys_list) + 1L
                    key_to_idx[[curr_key]] <- idx
                    curr_keys_list[[idx]] <- curr_key
                    curr_assign_list[[idx]] <- as.integer(assign0)
                    curr_logB[idx] <- -Inf
                    curr_vit[idx] <- Inf
                    curr_ptr[idx] <- NA_integer_
                }

                c_trans <- c_trans_vec[[e]]

                # Viterbi update (MAP)
                cand_vit <- c_trans + vit_next[[j]]
                if (!is.finite(curr_vit[idx]) || cand_vit < curr_vit[idx]) {
                    curr_vit[idx] <- cand_vit
                    curr_ptr[idx] <- j
                }

                # Backward log-sum-exp update
                cand_log <- (-c_trans / temperature) + logB_next[[j]]
                curr_logB[idx] <- log_add_exp(curr_logB[idx], cand_log)

                # Store edge for forward pass
                used_edges <- used_edges + 1L
                from_idx[[used_edges]] <- idx
                to_idx[[used_edges]] <- j
                logw[[used_edges]] <- (-c_trans / temperature)
            }
        }

        if (length(curr_keys_list) == 0L) {
            vcat(prefix, "DP produced no states at CensusID=", cc, ". Falling back to igraph.")
            out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, interval_years, anchor_start)
            out <- ensure_posterior_columns(out)
            return(out)
        }

        keys_full[[cc]] <- unlist(curr_keys_list, use.names = FALSE)
        assign_full[[cc]] <- curr_assign_list
        logB[[cc]] <- curr_logB
        vit_cost[[cc]] <- curr_vit
        vit_ptr[[cc]] <- curr_ptr

        if (used_edges == 0L) {
            vcat(prefix, "No feasible edges at CensusID=", cc, ". Falling back to igraph.")
            out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, interval_years, anchor_start)
            out <- ensure_posterior_columns(out)
            return(out)
        }
        edges[[cc]] <- data.table::data.table(
            from_idx = from_idx[seq_len(used_edges)],
            to_idx = to_idx[seq_len(used_edges)],
            logw = logw[seq_len(used_edges)]
        )

        vcat(prefix, "Finished CensusID=", cc, ": full_states=", length(keys_full[[cc]]), ", edges=", used_edges, ", dt=", sprintf("%.2fs", tic() - t_cc0))
    }

    # -----------------
    # Decode MAP path
    # -----------------
    vcat(prefix, "Decoding MAP path and writing ReconstructedStemID ...")
    start_idx <- which.min(vit_cost[[1L]])
    if (length(start_idx) == 0L || !is.finite(vit_cost[[1L]][start_idx])) {
        out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, interval_years, anchor_start)
        out <- ensure_posterior_columns(out)
        return(out)
    }
    map_idx <- integer(anchor_start)
    map_idx[1L] <- start_idx
    for (cc in seq_len(anchor_start - 1L)) {
        nxt <- vit_ptr[[cc]][map_idx[cc]]
        if (!is.finite(nxt) || is.na(nxt) || nxt < 1L) {
            out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, interval_years, anchor_start)
            out <- ensure_posterior_columns(out)
            return(out)
        }
        map_idx[cc + 1L] <- nxt
    }

    for (cc in seq_len(anchor_start)) {
        obs_idx <- tree_data[CensusID == cc & !is.na(DBH), which = TRUE]
        if (length(obs_idx) == 0L) next
        sv <- assign_full[[cc]][[map_idx[cc]]]
        if (length(sv) != length(obs_idx)) {
            out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, interval_years, anchor_start)
            out <- ensure_posterior_columns(out)
            return(out)
        }
        tree_data[obs_idx, ReconstructedStemID := track_ids[sv]]
        obs_to_mark <- obs_idx[is.na(tree_data$TrueStemID[obs_idx])]
        if (length(obs_to_mark) > 0L) {
            tree_data[obs_to_mark, ReconstructionMethod := "dp"]
        }
    }

    # -----------------
    # Forward pass for marginals
    # -----------------
    vcat(prefix, "Forward pass starting ...")
    # Start distribution: uniform over all reachable states at census 1.
    logalpha <- vector("list", anchor_start)
    logalpha[[1L]] <- rep.int(0, length(keys_full[[1L]]))

    for (cc in seq_len(anchor_start - 1L)) {
        ed <- edges[[cc]]
        if (is.null(ed) || nrow(ed) == 0L) {
            out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, interval_years, anchor_start)
            out <- ensure_posterior_columns(out)
            return(out)
        }
        la_from <- logalpha[[cc]][ed$from_idx]
        vals <- la_from + ed$logw
        dt <- data.table::data.table(to_idx = ed$to_idx, v = vals)
        dt <- dt[is.finite(v)]
        if (nrow(dt) == 0L) {
            out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, interval_years, anchor_start)
            out <- ensure_posterior_columns(out)
            return(out)
        }
        la_next_dt <- dt[, .(logalpha = log_sum_exp(v)), by = to_idx]
        la_next <- rep.int(-Inf, length(keys_full[[cc + 1L]]))
        la_next[la_next_dt$to_idx] <- la_next_dt$logalpha
        logalpha[[cc + 1L]] <- la_next

        vcat(prefix, "Forward step CensusID=", cc + 1L, ": reached ", sum(is.finite(la_next)), " / ", length(la_next), " states")
    }

    # Partition function Z = total weight of all paths ending at the fixed anchor state.
    # At anchor_start there is exactly one state.
    logZ <- logalpha[[anchor_start]][1L]
    if (!is.finite(logZ)) {
        # Fallback: compute from backward at census 1.
        logZ <- log_sum_exp(logB[[1L]])
    }

    vcat(prefix, "Computed logZ=", sprintf("%.3f", logZ))

    # -----------------
    # Observation-level marginals
    # -----------------
    anchor_set <- anchor_ids
    is_anchor_track <- track_ids %in% anchor_set

    for (cc in seq_len(anchor_start)) {
        obs_idx <- tree_data[CensusID == cc & !is.na(DBH), which = TRUE]
        if (length(obs_idx) == 0L) next

        # State posterior weights at this census
        lg <- logalpha[[cc]] + logB[[cc]] - logZ
        # Normalize defensively
        lg <- lg - log_sum_exp(lg)
        w <- exp(lg)

        n_obs <- length(obs_idx)
        prob_mat <- matrix(0, nrow = n_obs, ncol = K)
        st_assign <- assign_full[[cc]]
        for (s in seq_along(st_assign)) {
            ww <- w[s]
            if (!is.finite(ww) || ww <= 0) next
            a <- st_assign[[s]]
            # a is length n_obs, entries in 1..K
            for (j in seq_len(n_obs)) {
                prob_mat[j, a[j]] <- prob_mat[j, a[j]] + ww
            }
        }
        # Re-normalize each row (floating error)
        row_sums <- rowSums(prob_mat)
        for (j in seq_len(n_obs)) {
            if (row_sums[j] > 0) {
                prob_mat[j, ] <- prob_mat[j, ] / row_sums[j]
            }
        }

        # MAP assignment for DP_PosteriorReconstructedProb
        map_assign <- assign_full[[cc]][[map_idx[cc]]]

        top1_id <- integer(n_obs)
        top2_id <- integer(n_obs)
        top1_p <- numeric(n_obs)
        top2_p <- numeric(n_obs)
        entropy <- numeric(n_obs)
        p_map <- numeric(n_obs)
        p_unlinked <- numeric(n_obs)

        for (j in seq_len(n_obs)) {
            p <- prob_mat[j, ]
            # numeric stability
            p[p < 0] <- 0
            sp <- sum(p)
            if (sp > 0) p <- p / sp
            ord <- order(p, decreasing = TRUE)
            i1 <- ord[1L]
            top1_id[j] <- track_ids[i1]
            top1_p[j] <- p[i1]
            if (posterior_top_k >= 2L && length(ord) >= 2L) {
                i2 <- ord[2L]
                top2_id[j] <- track_ids[i2]
                top2_p[j] <- p[i2]
            } else {
                top2_id[j] <- NA_integer_
                top2_p[j] <- NA_real_
            }
            entropy[j] <- -sum(ifelse(p > 0, p * log(p), 0))
            p_map[j] <- p[map_assign[j]]
            p_unlinked[j] <- sum(p[!is_anchor_track])
        }

        tree_data[obs_idx, `:=`(
            DP_PosteriorTop1ID = top1_id,
            DP_PosteriorTop1Prob = top1_p,
            DP_PosteriorTop2ID = top2_id,
            DP_PosteriorTop2Prob = top2_p,
            DP_PosteriorEntropy = entropy,
            DP_PosteriorReconstructedProb = p_map,
            DP_PosteriorUnlinkedProb = p_unlinked
        )]
    }

    tree_data <- add_constraint_violation(
        tree_data,
        id_col = "ReconstructedStemID",
        min_growth = min_growth,
        max_growth = max_growth,
        interval_years = interval_years
    )

    vcat(prefix, "Done. Total elapsed ", sprintf("%.2fs", tic() - t_start))
    return(tree_data)
}

# ----------------------------
# Plotting: one PDF page per Tag
# ----------------------------
# This helper is NOT part of the matching algorithm.
# It is only a visualization tool to help you see whether reconstructed IDs
# create smooth DBH trajectories across censuses.
#
# For each Tag, we optionally create two side-by-side plots:
# 1) DBH trajectories grouped by `ReferenceStemID` (only if present and `include_reference=TRUE`)
# 2) DBH trajectories grouped by `ReconstructedStemID`
#
# Each line is a "stem ID group" and points are DBH measurements by census.
# If DP posterior bins are present (DP_PosteriorBin), point shapes show
# confident/ambiguous/unlinked-likely, and constraint violations are overlaid
# as an "X" symbol.

plot_tag_to_pdf <- function(out, pdf_file, include_reference = FALSE, tag = NULL) {
    ## ---- Package checks (no attach unless needed) ----
    # We use requireNamespace() so we don't attach these packages globally.
    pkgs <- c("ggplot2", "cowplot")
    invisible(lapply(pkgs, function(p) {
        if (!requireNamespace(p, quietly = TRUE)) {
            stop("Package not installed: ", p, call. = FALSE)
        }
    }))

    ## ---- Required columns ----
    # Plot 2 (reconstructed IDs) always requires these.
    required_cols <- c(
        "Tag", "CensusID", "DBH", "species",
        "ReconstructedStemID", "ReconstructionMethod", "ConstraintViolation"
    )
    # Plot 1 (reference/original IDs) is optional.
    if (isTRUE(include_reference)) {
        required_cols <- c(required_cols, "OriginalStemID")
    }

    # Find which required columns are absent.
    missing_cols <- setdiff(required_cols, names(out))
    if (length(missing_cols) > 0) {
        stop("Missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
    }

    # Get out_dir for subtitle if available
    out_dir_value <- if ("out_dir" %in% names(out)) unique(out$out_dir) else NULL

    ## ---- Optional filter: one tag (or subset of tags) ----
    if (!is.null(tag)) {
        tag <- sort(unique(as.integer(tag)))
        out <- out[Tag %in% tag]
        if (nrow(out) == 0L) {
            stop("No rows found for requested tag(s): ", paste(tag, collapse = ", "), call. = FALSE)
        }
    }

    ## ---- Output filename ----
    # If `tag` is specified, we automatically write to a different filename than
    # the "all tags" PDF so you don't overwrite the full report.
    resolve_pdf_file <- function(path, tag) {
        if (!is.character(path) || length(path) != 1L || !nzchar(path)) {
            stop("pdf_file must be a single non-empty path", call. = FALSE)
        }

        if (dir.exists(path)) {
            base <- if (is.null(tag)) {
                "stem_reconstruction_all_tags.pdf"
            } else if (length(tag) == 1L) {
                paste0("stem_reconstruction_tag_", tag, ".pdf")
            } else {
                paste0("stem_reconstruction_tags_", paste(tag, collapse = "_"), ".pdf")
            }
            return(file.path(path, base))
        }

        ext <- tolower(tools::file_ext(path))
        if (ext != "pdf") {
            stop("pdf_file must be a .pdf file or an existing directory", call. = FALSE)
        }
        dirn <- dirname(path)
        base0 <- tools::file_path_sans_ext(basename(path))

        if (is.null(tag)) {
            return(path)
        }
        if (length(tag) == 1L) {
            return(file.path(dirn, paste0(base0, "_tag_", tag, ".pdf")))
        }
        file.path(dirn, paste0(base0, "_tags_", paste(tag, collapse = "_"), ".pdf"))
    }

    pdf_file_final <- resolve_pdf_file(pdf_file, tag)

    ## ---- Helper: common theme & scales ----
    # A common theme makes the two panels comparable and consistent.
    common_theme <- ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
            plot.margin = ggplot2::margin(10, 10, 10, 10),
            # legend to bottom and horizontal
            # legend.position = "bottom"#,
            # legend.direction = "horizontal"
        )

    # Shape scales
    # - If DP posterior bins exist, shapes show the bin.
    # - Otherwise, shapes show constraint violations.
    bin_shape_scale <- ggplot2::scale_shape_manual(
        values = c(
            `confident` = 16,
            `ambiguous` = 17,
            `unlinked-likely` = 1,
            `NA` = 4
        ),
        drop = FALSE
    )

    violation_shape_scale <- ggplot2::scale_shape_manual(
        values = c(`FALSE` = 16, `TRUE` = 17),
        labels = c(`FALSE` = "ok", `TRUE` = "violation")
    )

    ## ---- Open PDF ----
    # onefile=TRUE means: one PDF file with multiple pages.
    grDevices::pdf(pdf_file_final, width = 13, height = 8, onefile = TRUE)

    ## ---- Loop over Tag ----
    # Each iteration creates ONE page in the PDF.
    keys <- data.table::data.table(Tag = sort(unique(out$Tag)))

    for (ii in seq_len(nrow(keys))) {
        # ii <- 1
        tag <- keys$Tag[[ii]]
        species_code <- unique(out[Tag == tag]$species)
        # Filter to the current Tag (and species, if applicable).
        tag_data <- data.table::copy(out[Tag == tag])

        tag_data[, ConstraintViolationFlag := (!is.na(ConstraintViolation) & as.logical(ConstraintViolation))]

        has_bins <- ("DP_PosteriorBin" %in% names(tag_data)) && any(!is.na(tag_data$DP_PosteriorBin))

        p1 <- NULL
        if (isTRUE(include_reference)) {
            # Plot 1: Original/reference IDs (if present in `out`)
            if (isTRUE(has_bins)) {
                p1 <- ggplot2::ggplot(
                    tag_data,
                    ggplot2::aes(
                        x = as.factor(CensusID),
                        y = DBH,
                        group = OriginalStemID,
                        color = factor(OriginalStemID),
                        shape = DP_PosteriorBin
                    )
                ) +
                    ggplot2::geom_line(na.rm = TRUE) +
                    ggplot2::geom_point(size = 3, na.rm = TRUE) +
                    ggplot2::geom_point(
                        data = tag_data[ConstraintViolationFlag == TRUE],
                        ggplot2::aes(x = as.factor(CensusID), y = DBH, group = OriginalStemID),
                        inherit.aes = FALSE,
                        shape = 4,
                        color = "black",
                        size = 3.5,
                        stroke = 1,
                        na.rm = TRUE
                    ) +
                    bin_shape_scale +
                    ggplot2::labs(
                        title = paste("Original Stem IDs - Tag", tag, "(species:", species_code, ")"),
                        subtitle = out_dir_value,
                        color = "OriginalStemID",
                        shape = "Posterior bin"
                    ) +
                    common_theme
            } else {
                p1 <- ggplot2::ggplot(
                    tag_data,
                    ggplot2::aes(
                        x = as.factor(CensusID),
                        y = DBH,
                        group = OriginalStemID,
                        color = factor(OriginalStemID),
                        shape = ConstraintViolationFlag
                    )
                ) +
                    ggplot2::geom_line(na.rm = TRUE) +
                    ggplot2::geom_point(size = 3, na.rm = TRUE) +
                    violation_shape_scale +
                    ggplot2::labs(
                        title = paste("Original Stem IDs - Tag", tag, "(species:", species_code, ")"),
                        subtitle = out_dir_value,
                        color = "OriginalStemID",
                        shape = "ConstraintViolation"
                    ) +
                    common_theme
            }
        }

        ## --- Plot 2: ReconstructedStemID ---
        # Group by ReconstructedStemID: this shows your algorithm's result.
        # Ideally, DBH trajectories should look smooth/biologically plausible.
        if (isTRUE(has_bins)) {
            p2 <- ggplot2::ggplot(
                tag_data,
                ggplot2::aes(
                    x = as.factor(CensusID),
                    y = DBH,
                    group = ReconstructedStemID,
                    color = factor(ReconstructedStemID),
                    shape = DP_PosteriorBin
                )
            ) +
                ggplot2::geom_line(na.rm = TRUE) +
                ggplot2::geom_point(size = 3, na.rm = TRUE) +
                ggplot2::geom_point(
                    data = tag_data[ConstraintViolationFlag == TRUE],
                    ggplot2::aes(x = as.factor(CensusID), y = DBH, group = ReconstructedStemID),
                    inherit.aes = FALSE,
                    shape = 4,
                    color = "black",
                    size = 3.5,
                    stroke = 1,
                    na.rm = TRUE
                ) +
                bin_shape_scale +
                ggplot2::labs(
                    title = paste("Reconstructed Stem IDs - Tag", tag, "(species:", species_code, ")"),
                    subtitle = out_dir_value,
                    color = "ReconstructedStemID",
                    shape = "Posterior bin"
                ) +
                common_theme
        } else {
            p2 <- ggplot2::ggplot(
                tag_data,
                ggplot2::aes(
                    x = as.factor(CensusID),
                    y = DBH,
                    group = ReconstructedStemID,
                    color = factor(ReconstructedStemID),
                    shape = ConstraintViolationFlag
                )
            ) +
                ggplot2::geom_line(na.rm = TRUE) +
                ggplot2::geom_point(size = 3, na.rm = TRUE) +
                violation_shape_scale +
                ggplot2::labs(
                    title = paste("Reconstructed Stem IDs - Tag", tag, "(species:", species_code, ")"),
                    subtitle = out_dir_value,
                    color = "ReconstructedStemID",
                    shape = "ConstraintViolation"
                ) +
                common_theme
        }

        # Combine panels
        if (isTRUE(include_reference)) {
            combined_plot <- cowplot::plot_grid(
                p1, p2,
                ncol = 2,
                align = "h",
                axis = "tb",
                rel_widths = c(1, 1)
            )
            print(combined_plot)
        } else {
            print(p2)
        }
    }

    # Always close the PDF device.
    grDevices::dev.off()
    invisible(pdf_file_final)
}

# To estimate rates of error, we performed a double-blind re-measurement of 1715
# trees in 1995 and 2000 (Condit 1998) and fitted the discrepancies with a sum
# of two normal distributions. The first describes small errors and has an s.d.
# (SD1) proportional to the trunk diameter; the second has a fixed larger s.d.
# (SD2). The 1715 errors were best fit with SD1 = 0.0062 × D + 0.0904, SD2 = 4.64
# (all units in centimetres), with 5% of the trees subject to the larger error.
# For example, the diameter of a 30 cm tree has a typical error of 0.27 cm (95%
# probability) or of 4.63 cm (5% probability).