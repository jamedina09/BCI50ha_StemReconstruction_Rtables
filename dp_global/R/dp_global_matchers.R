############################################################
# dp_global_matchers.R
# Fallback stepwise matching (igraph-based)
############################################################

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
    # interval_years are resolved per-census-pair when needed via resolve_interval_years_pair()
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

        pair_interval <- resolve_interval_years_pair(tree_data, t0 = c, t1 = c + 1L, interval_years = interval_years)

        growth_mat <- matrix(
            outer(curr_dbh, fut_dbh, FUN = function(d0, d1) (d1 - d0) / pair_interval),
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
