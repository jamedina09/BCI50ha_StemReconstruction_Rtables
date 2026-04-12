############################################################
# dp_probabilistic_matching.R
# Probabilistic greedy matching fallback for large state spaces
############################################################
# When the DP state space is too large:
#   1. Pairwise log-likelihoods (same bio model as DP)
#   2. Augment cost matrix with mortality/recruitment slots
#   3. Draw n_samples stochastic assignments via Gumbel-noise greedy
#   4. Stitch per-pair assignments backward from anchor
#   5. Repair growth violations at SAMPLE level (hard-rate + ME cumulative-shrinkage)
#   6. Compute marginal posterior probabilities
#   7. Growth-aware greedy conflict resolution (anchor-outward census ordering,
#      rejects candidate IDs that would violate growth bounds against already-
#      resolved adjacent censuses)
#   8. Re-stamp anchor TrueStemID rows
#
# All Bio_* parameters are read directly from tree_data columns (no new
# estimation needed — they are already computed by dp_global_bio.R).

# ---- Main entry point ----------------------------------------------------

match_stems_probabilistic <- function(tree_data,
                                      min_growth,
                                      max_growth,
                                      anchor_start,
                                      n_samples     = 200L,
                                      temperature   = 1.0,
                                      posterior_top_k = 2L,
                                      posterior_samples_path = NULL,
                                      posterior_samples_format = "csv",
                                      posterior_sample_seed = NULL,
                                      prune_min_growth    = NULL,
                                      prune_max_growth    = NULL,
                                      prune_recruit_max_dbh = NULL,
                                      prob_lookahead_weight = 0.5,  # backward conditioning weight [0,1]; 0 = independent
                                      pin_truestemid = TRUE,        # pin obs with known TrueStemID to their track
                                      verbose       = FALSE) {
    tree_data <- tree_data[order(CensusID)]
    n_samples <- as.integer(n_samples)
    posterior_top_k <- max(1L, as.integer(posterior_top_k))

    vcat <- function(...) {
        if (!isTRUE(verbose)) return(invisible(NULL))
        cat(..., "\n"); flush.console(); invisible(NULL)
    }

    tag_val <- tryCatch({
        u <- unique(tree_data$Tag); u <- u[!is.na(u)]
        if (length(u) == 1L) u[[1L]] else NA
    }, error = function(e) NA)
    prefix <- paste0("[prob_match Tag=", if (!is.na(tag_val)) tag_val else "?", "] ")

    # --- Resolve effective prune bounds (mirrors DP logic) ----------------
    eff_min_growth <- if (!is.null(prune_min_growth)) prune_min_growth else min_growth
    eff_max_growth <- if (!is.null(prune_max_growth)) prune_max_growth else max_growth

    # --- Extract bio parameters from tree_data (same as dp_global_dp.R) ---
    bio <- list(
        mu_const   = unique(tree_data$Bio_Mu_Growth)[1],
        mu_gamma   = if ("Bio_Gamma_Growth" %in% names(tree_data)) unique(tree_data$Bio_Gamma_Growth)[1] else 0,
        sigma0     = unique(tree_data$Bio_Sigma0_Growth)[1],
        sigma1     = unique(tree_data$Bio_Sigma1_Growth)[1],
        max_shrink = unique(tree_data$Bio_Max_Shrink)[1],
        k_shrink   = unique(tree_data$Bio_K_Shrink)[1],
        max_growth_bio = if ("Bio_Max_Growth" %in% names(tree_data)) unique(tree_data$Bio_Max_Growth)[1] else Inf,
        k_growth   = if ("Bio_K_Growth" %in% names(tree_data)) unique(tree_data$Bio_K_Growth)[1] else 0,
        h0         = unique(tree_data$Bio_H0_Mortality)[1],
        beta_mort  = unique(tree_data$Bio_Beta_Mortality)[1],
        recruit_meanlog = unique(tree_data$Bio_Recruit_Meanlog)[1],
        recruit_sdlog   = unique(tree_data$Bio_Recruit_Sdlog)[1],
        recruit_max_dbh = if (!is.null(prune_recruit_max_dbh)) prune_recruit_max_dbh else unique(tree_data$Bio_Recruit_MaxDBH_unit)[1],
        recruit_lambda  = unique(tree_data$Bio_Recruitment_lambda)[1]
    )

    # --- Identify censuses and observations --------------------------------
    if (!("ReconstructionMethod" %in% names(tree_data)))
        tree_data[, ReconstructionMethod := NA_character_]
    if (!("ConstraintViolation" %in% names(tree_data)))
        tree_data[, ConstraintViolation := NA]
    if (!("obs_row_id" %in% names(tree_data)))
        tree_data[, obs_row_id := seq_len(.N)]

    # Anchor census: use TrueStemID if available, else provisional
    anchor_obs <- tree_data[CensusID == anchor_start & !is.na(DBH)]
    has_anchor <- nrow(anchor_obs) > 0L && any(!is.na(anchor_obs$TrueStemID))
    if (!has_anchor) {
        # Try last census with observed DBH as provisional anchor
        obs_census <- sort(unique(tree_data$CensusID[!is.na(tree_data$DBH)]))
        if (length(obs_census) == 0L) {
            # No DBH anywhere — cannot match; leave ReconstructedStemID as NA
            tree_data[, ReconstructedStemID := NA_integer_]
            tree_data[, ReconstructionMethod := "probabilistic"]
            return(tree_data)
        }
        anchor_start <- max(obs_census)
        anchor_obs <- tree_data[CensusID == anchor_start & !is.na(DBH)]
    }

    # Census range (only censuses with observed DBH, up to anchor)
    obs_census <- sort(unique(tree_data$CensusID[!is.na(tree_data$DBH) & tree_data$CensusID <= anchor_start]))
    n_census <- length(obs_census)

    if (n_census <= 1L) {
        # Single census — assign sequential IDs
        tree_data[CensusID == anchor_start & !is.na(DBH),
                  ReconstructedStemID := seq_len(.N)]
        tree_data[is.na(ReconstructedStemID), ReconstructedStemID := NA_integer_]
        tree_data[, ReconstructionMethod := "probabilistic"]
        return(tree_data)
    }

    # Build per-census observation data
    obs_data <- vector("list", n_census)
    for (i in seq_len(n_census)) {
        cc <- obs_census[i]
        idx <- which(tree_data$CensusID == cc & !is.na(tree_data$DBH))
        obs_data[[i]] <- list(
            census_id = cc,
            idx       = idx,                       # row indices in tree_data
            dbh       = tree_data$DBH[idx],
            row_id    = tree_data$obs_row_id[idx],
            n         = length(idx)
        )
    }

    # Compute interval years between consecutive observed censuses
    .dt_pi <- tree_data[, .(MeanDate = as.numeric(mean(as.numeric(ExactDate), na.rm = TRUE))), by = CensusID]
    intervals <- numeric(n_census - 1L)
    for (i in seq_len(n_census - 1L)) {
        d0 <- .dt_pi$MeanDate[.dt_pi$CensusID == obs_census[i]]
        d1 <- .dt_pi$MeanDate[.dt_pi$CensusID == obs_census[i + 1L]]
        intervals[i] <- (d1 - d0) / 365.25
        if (!is.finite(intervals[i]) || intervals[i] <= 0) intervals[i] <- 5.0  # safe fallback
    }

    # --- Anchor IDs --------------------------------------------------------
    # Build IDs for ALL anchor observations (one per row).  Where TrueStemID
    # is available, use it; where NA, assign new sequential IDs starting
    # above the max known TrueStemID so they don't collide.
    anchor_pos <- n_census  # anchor is the last observed census
    anchor_ids <- integer(nrow(anchor_obs))
    has_true <- !is.na(anchor_obs$TrueStemID)
    if (any(has_true)) {
        anchor_ids[has_true] <- as.integer(anchor_obs$TrueStemID[has_true])
        next_id <- max(anchor_ids[has_true]) + 1L
    } else {
        next_id <- 1L
    }
    if (any(!has_true)) {
        n_missing <- sum(!has_true)
        anchor_ids[!has_true] <- seq.int(next_id, length.out = n_missing)
    }
    # Ensure all anchor obs have IDs
    # Preserve "provisional_dp" for rows where TrueStemID was fabricated;
    # mark real TrueStemID as "given"; the rest are "probabilistic".
    .anchor_rows <- which(tree_data$CensusID == anchor_start & !is.na(tree_data$DBH))
    tree_data[.anchor_rows, ReconstructedStemID := anchor_ids]
    .is_provisional <- tree_data$ReconstructionMethod[.anchor_rows] %in% "provisional_dp"
    .has_tsid       <- !is.na(tree_data$TrueStemID[.anchor_rows])
    .method_vec     <- ifelse(.is_provisional, "provisional_dp",
                       ifelse(.has_tsid, "given", "probabilistic"))
    tree_data[.anchor_rows, ReconstructionMethod := .method_vec]

    # K = number of tracks (at least max obs across any census)
    max_obs <- max(vapply(obs_data, function(x) x$n, integer(1)))
    K <- max(length(anchor_ids), max_obs)

    vcat(prefix, "Probabilistic matching: ", n_census, " censuses, K=", K,
         ", max_obs=", max_obs, ", n_samples=", n_samples)

    # --- Pre-compute TrueStemID pin map for non-anchor censuses -------------
    # pin_info[[i]][j] = anchor-position index (1..n_anchor) for obs j, or NA
    pin_info <- vector("list", n_census)
    .any_pins <- FALSE
    if (isTRUE(pin_truestemid)) {
        n_anchor <- length(anchor_ids)
        for (i in seq_len(n_census)) {
            if (i == n_census) next  # anchor pinned via anchor_ids directly
            n_obs_i <- obs_data[[i]]$n
            if (n_obs_i == 0L) next
            tsid <- tree_data$TrueStemID[obs_data[[i]]$idx]
            tidx <- match(as.integer(tsid), anchor_ids)
            tidx[is.na(tsid)] <- NA_integer_
            # Duplicate-pin guard: if two obs claim the same anchor position, keep first
            .seen <- integer(0)
            for (.j in seq_along(tidx)) {
                if (is.na(tidx[.j])) next
                if (tidx[.j] %in% .seen) {
                    vcat(prefix, "WARNING: duplicate TrueStemID pin at C",
                         obs_data[[i]]$census_id, " for anchor ID ", anchor_ids[tidx[.j]],
                         "; keeping first, releasing obs ", .j)
                    tidx[.j] <- NA_integer_
                } else {
                    .seen <- c(.seen, tidx[.j])
                }
            }
            if (any(!is.na(tidx))) .any_pins <- TRUE
            pin_info[[i]] <- tidx
        }
    }

    # --- Per-pair log-likelihood matrices (backward) -----------------------
    # For each pair (c, c+1) compute the pairwise + augmented cost matrix
    pair_data <- vector("list", n_census - 1L)
    for (i in seq_len(n_census - 1L)) {
        dbh_curr <- obs_data[[i]]$dbh
        dbh_next <- obs_data[[i + 1L]]$dbh
        iv <- intervals[i]

        L <- compute_pairwise_log_likelihood(dbh_curr, dbh_next, iv, bio,
                                             eff_min_growth, eff_max_growth)
        aug <- augment_cost_matrix(L, dbh_curr, dbh_next, iv, bio)

        pair_data[[i]] <- list(
            log_cost = aug,
            n_curr   = length(dbh_curr),
            n_next   = length(dbh_next)
        )
    }

    # --- Draw n_samples stochastic assignments (backward from anchor) ------
    if (!is.null(posterior_sample_seed)) set.seed(as.integer(posterior_sample_seed))

    all_samples <- vector("list", n_samples)
    use_lookahead <- is.finite(prob_lookahead_weight) && prob_lookahead_weight > 0 &&
                     n_census >= 3L && K >= 4L

    for (s in seq_len(n_samples)) {
        # For each census pair (working backward from anchor-1 to 1),
        # sample an assignment.  When lookahead is enabled, after sampling
        # pair (i+1), condition pair (i)'s cost matrix on the result.
        # When pinning is active, maintain per-sample track-to-anchor mapping
        # so pinned obs can be forced to the correct column.
        per_pair_assignments <- vector("list", n_census - 1L)

        # Initialize anchor-position mapping: at anchor, obs j IS position j
        .n_anchor_obs <- obs_data[[n_census]]$n
        .next_obs_to_anchor_pos <- seq_len(.n_anchor_obs)

        # Last pair (closest to anchor): no conditioning available
        last_pair <- n_census - 1L
        .cost_last <- pair_data[[last_pair]]$log_cost
        if (.any_pins && !is.null(pin_info[[last_pair]])) {
            .cost_last <- apply_pin_mask(
                .cost_last, pin_info[[last_pair]], .next_obs_to_anchor_pos,
                pair_data[[last_pair]]$n_curr, pair_data[[last_pair]]$n_next
            )
        }
        per_pair_assignments[[last_pair]] <- greedy_assignment_gumbel(
            .cost_last,
            temperature = temperature
        )
        if (.any_pins) {
            .next_obs_to_anchor_pos <- propagate_track_backward(
                per_pair_assignments[[last_pair]], .next_obs_to_anchor_pos,
                pair_data[[last_pair]]$n_curr, pair_data[[last_pair]]$n_next
            )
        }

        # Remaining pairs moving backward: condition on next pair's assignment
        if (last_pair >= 2L) {
            for (i in seq.int(last_pair - 1L, 1L, by = -1L)) {
                cost_i <- pair_data[[i]]$log_cost

                if (use_lookahead) {
                    cost_i <- condition_cost_matrix(
                        aug_cost        = cost_i,
                        n_curr          = pair_data[[i]]$n_curr,
                        n_next          = pair_data[[i]]$n_next,
                        dbh_next        = obs_data[[i + 1L]]$dbh,
                        next_assignment = per_pair_assignments[[i + 1L]],
                        n_next_next     = pair_data[[i + 1L]]$n_next,
                        dbh_further     = obs_data[[i + 2L]]$dbh,
                        interval_next   = intervals[i + 1L],
                        bio             = bio,
                        weight          = prob_lookahead_weight
                    )
                }

                if (.any_pins && !is.null(pin_info[[i]])) {
                    cost_i <- apply_pin_mask(
                        cost_i, pin_info[[i]], .next_obs_to_anchor_pos,
                        pair_data[[i]]$n_curr, pair_data[[i]]$n_next
                    )
                }

                per_pair_assignments[[i]] <- greedy_assignment_gumbel(
                    cost_i,
                    temperature = temperature
                )
                if (.any_pins) {
                    .next_obs_to_anchor_pos <- propagate_track_backward(
                        per_pair_assignments[[i]], .next_obs_to_anchor_pos,
                        pair_data[[i]]$n_curr, pair_data[[i]]$n_next
                    )
                }
            }
        }

        all_samples[[s]] <- per_pair_assignments
    }

    # --- Stitch assignments backward from anchor ---------------------------
    stitched <- stitch_assignments_backward(all_samples, obs_data, anchor_ids, K)

    # --- Repair growth violations at the SAMPLE level (before marginals) ---
    # This ensures probabilities only count biologically valid paths.
    # Two layers: hard-rate bounds + ME-informed cumulative shrinkage.
    stitched <- repair_stitched_growth_violations(
        stitched, obs_data, intervals, eff_min_growth, eff_max_growth,
        me_sd1_a = 0.0062, me_sd1_b = 0.0904, n_sigma_me = 3
    )
    .sample_breaks    <- attr(stitched, "sample_level_breaks")
    .sample_me_breaks <- attr(stitched, "sample_level_me_breaks")
    if (!is.null(.sample_breaks) && .sample_breaks > 0L) {
        .msg <- paste0(prefix, "Sample-level repair: ", .sample_breaks,
             " growth violation(s) broken across ", n_samples, " samples",
             if (!is.null(.sample_me_breaks) && .sample_me_breaks > 0L)
                 paste0(" (", .sample_me_breaks, " from ME cumulative-shrinkage check)")
             else "")
        vcat(.msg)
        message(.msg)  # ensure it appears on stderr / captured by log redirection
    }

    # --- Filter to pin-consistent samples (before marginals) ---------------
    # Discard samples where any pinned obs got the wrong track ID, so
    # marginals are conditioned on all pins being correct.
    if (.any_pins) {
        stitched <- filter_pin_consistent_samples(
            stitched, pin_info, anchor_ids, n_census,
            min_keep = 10L, vcat = vcat, prefix = prefix
        )
    }

    # --- Compute marginals and fill tree_data ------------------------------
    # Growth-aware greedy resolver: resolves censuses from anchor outward and
    # rejects candidate IDs whose growth rate against the nearest already-resolved
    # census violates hard bounds.  Posteriors (Top-K, entropy) are unaffected.
    tree_data <- compute_marginals_from_samples(stitched, tree_data, obs_data,
                                                obs_census, anchor_pos,
                                                posterior_top_k,
                                                intervals = intervals,
                                                min_rate = eff_min_growth,
                                                max_rate = eff_max_growth)

    # --- Diagnostic check: count residual growth violations from greedy
    #     conflict resolution.  With growth-aware resolver + pin-consistent
    #     sample filtering + sample-level repair, ideally 0 violations remain.
    #     We do NOT modify tree_data here — posteriors must stay pristine.
    {
        .n_violations <- 0L
        .stem_ids <- unique(tree_data$ReconstructedStemID[!is.na(tree_data$ReconstructedStemID)])
        for (.sid in .stem_ids) {
            .rows <- which(tree_data$ReconstructedStemID == .sid)
            if (length(.rows) < 2L) next
            .rows <- .rows[order(tree_data$CensusID[.rows])]
            for (.ri in seq_len(length(.rows) - 1L)) {
                .d1 <- tree_data$DBH[.rows[.ri]]
                .d2 <- tree_data$DBH[.rows[.ri + 1L]]
                if (is.na(.d1) || is.na(.d2)) next
                .iv <- intervals[match(tree_data$CensusID[.rows[.ri]], obs_census)]
                if (is.na(.iv) || .iv <= 0) next
                .rate <- (.d2 - .d1) / .iv
                if (.rate < eff_min_growth || .rate > eff_max_growth) {
                    .n_violations <- .n_violations + 1L
                }
            }
        }
        if (.n_violations > 0L) {
            .msg <- paste0(prefix, "WARNING: ", .n_violations,
                           " residual growth violation(s) in marginal trajectory ",
                           "(post-marginal repair DISABLED — posteriors preserved)")
            vcat(.msg)
            message(.msg)
        }
    }

    # --- Label assignment (does NOT modify ReconstructedStemID or posteriors)
    # With pin-consistent sample filtering, marginals already reflect the
    # correct IDs at anchor and pinned rows.  Only set method labels here.
    #
    # Label logic:
    #   anchor + real TrueStemID (not provisional)  -> "given"
    #   anchor + provisional TrueStemID             -> "provisional_dp"
    #   pre-anchor + real TrueStemID + pin active   -> "given"
    #   everything else                             -> "probabilistic"
    .is_anchor <- tree_data$CensusID == anchor_start
    .has_tsid  <- !is.na(tree_data$TrueStemID)
    .is_prov   <- tree_data$ReconstructionMethod %in% "provisional_dp"
    .is_pinned <- .has_tsid & !.is_prov  # real TrueStemID = pinned
    tree_data[, ReconstructionMethod := "probabilistic"]
    tree_data[.is_pinned & (.is_anchor | isTRUE(pin_truestemid)), ReconstructionMethod := "given"]
    tree_data[.is_prov, ReconstructionMethod := "provisional_dp"]

    # --- Export posterior samples (same format as DP) -----------------------
    if (n_samples > 0L) {
        export_probabilistic_posteriors(
            stitched, tree_data, obs_data, obs_census,
            tag_val         = tag_val,
            n_samples       = length(stitched),
            posterior_samples_path   = posterior_samples_path,
            posterior_samples_format = posterior_samples_format,
            verbose         = verbose,
            prefix          = prefix,
            vcat            = vcat
        )
    }

    tree_data
}

# ---- Pairwise log-likelihood matrix --------------------------------------

compute_pairwise_log_likelihood <- function(dbh_curr, dbh_next, interval_years,
                                            bio, min_growth, max_growth) {
    n_curr <- length(dbh_curr)
    n_next <- length(dbh_next)
    L <- matrix(-Inf, nrow = n_curr, ncol = n_next)

    mu_growth_fn <- function(d) {
        if (!is.finite(bio$mu_gamma) || bio$mu_gamma == 0 || !is.finite(d) || d <= 0)
            return(bio$mu_const)
        bio$mu_const + bio$mu_gamma * log(d)
    }

    for (i in seq_len(n_curr)) {
        d0 <- dbh_curr[i]
        if (!is.finite(d0)) next
        for (j in seq_len(n_next)) {
            d1 <- dbh_next[j]
            if (!is.finite(d1)) next

            g <- (d1 - d0) / interval_years
            # Hard growth constraints — infeasible edge
            if (is.finite(min_growth) && g < min_growth) next
            if (is.finite(max_growth) && g > max_growth) next
            if (is.finite(bio$max_shrink) && g < bio$max_shrink) next

            # Growth likelihood (Gaussian)
            sigma_d <- max(bio$sigma0 + bio$sigma1 * d0, 1e-6)
            mu <- mu_growth_fn(d0)
            ll_growth <- dnorm(g, mean = mu, sd = sigma_d, log = TRUE)

            # Survival probability
            hazard <- bio$h0 * exp(bio$beta_mort * d0)
            p_surv <- exp(-hazard * interval_years)
            p_surv <- max(1e-12, min(1 - 1e-12, p_surv))
            ll_surv <- log(p_surv)

            # Soft penalties (same as C++)
            ll_soft <- 0
            if (is.finite(bio$k_shrink) && bio$k_shrink > 0 && d1 < d0) {
                ll_soft <- ll_soft - bio$k_shrink * (d0 - d1)^2
            }
            if (is.finite(bio$k_growth) && bio$k_growth > 0 &&
                is.finite(bio$max_growth_bio)) {
                d1_cap <- d0 + bio$max_growth_bio * interval_years
                if (is.finite(d1_cap) && d1 > d1_cap) {
                    ll_soft <- ll_soft - bio$k_growth * (d1 - d1_cap)^2
                }
            }

            L[i, j] <- ll_growth + ll_surv + ll_soft
        }
    }
    L
}

# ---- Augment cost matrix with mortality/recruitment ----------------------

augment_cost_matrix <- function(L, dbh_curr, dbh_next, interval_years, bio) {
    n_curr <- length(dbh_curr)
    n_next <- length(dbh_next)

    # Adaptive K: start from max(n_curr, n_next), then ensure enough
    # death columns for rows where ALL survival entries are -Inf, and
    # enough recruit rows for cols where ALL survival entries are -Inf.
    must_die <- 0L
    if (n_curr > 0L && n_next > 0L) {
        for (i in seq_len(n_curr))
            if (all(L[i, ] == -Inf)) must_die <- must_die + 1L
    } else if (n_curr > 0L) {
        must_die <- n_curr
    }

    must_recruit <- 0L
    if (n_curr > 0L && n_next > 0L) {
        for (j in seq_len(n_next))
            if (all(L[, j] == -Inf)) must_recruit <- must_recruit + 1L
    } else if (n_next > 0L) {
        must_recruit <- n_next
    }

    K <- max(n_curr, n_next)
    # Add extra death columns if needed
    death_avail <- max(0L, K - n_next)
    if (death_avail < must_die) K <- K + (must_die - death_avail)
    # Add extra recruit rows if needed
    recruit_avail <- max(0L, K - n_curr)
    if (recruit_avail < must_recruit) K <- K + (must_recruit - recruit_avail)

    # Augmented matrix: K rows × K cols
    # Rows 1..n_curr are real current stems; rows (n_curr+1)..K are virtual recruit sources
    # Cols 1..n_next are real next stems; cols (n_next+1)..K are virtual death sinks
    A <- matrix(-Inf, nrow = K, ncol = K)

    # Fill the survival sub-matrix
    if (n_curr > 0 && n_next > 0) {
        A[seq_len(n_curr), seq_len(n_next)] <- L
    }

    # Death columns: current stem dies (cols n_next+1 .. K)
    if (n_next < K) {
        for (i in seq_len(n_curr)) {
            d0 <- dbh_curr[i]
            hazard <- bio$h0 * exp(bio$beta_mort * d0)
            p_death <- 1 - exp(-hazard * interval_years)
            p_death <- max(1e-12, min(1 - 1e-12, p_death))
            # Each death column is equivalent — log P(death)
            for (jj in (n_next + 1L):K) {
                A[i, jj] <- log(p_death)
            }
        }
    }

    # Recruitment rows: new stem appears (rows n_curr+1 .. K)
    if (n_curr < K) {
        p_recruit <- 1 - exp(-bio$recruit_lambda * interval_years)
        p_recruit <- max(1e-12, min(1 - 1e-12, p_recruit))
        for (ii in (n_curr + 1L):K) {
            for (j in seq_len(n_next)) {
                d1 <- dbh_next[j]
                if (!is.finite(d1) || d1 <= 0) next
                # Recruitment: lognormal size distribution
                ll_recruit <- log(p_recruit) +
                    dlnorm(d1, meanlog = bio$recruit_meanlog,
                           sdlog = bio$recruit_sdlog, log = TRUE)
                # Respect recruit max DBH
                if (is.finite(bio$recruit_max_dbh) && d1 > bio$recruit_max_dbh) {
                    ll_recruit <- -Inf
                }
                A[ii, j] <- ll_recruit
            }
        }
    }

    # Virtual-to-virtual (recruit source dies): very low probability placeholder
    if (n_curr < K && n_next < K) {
        for (ii in (n_curr + 1L):K) {
            for (jj in (n_next + 1L):K) {
                A[ii, jj] <- -20  # small log-prob: neither recruit nor die
            }
        }
    }

    A
}

# ---- Condition cost matrix with lookahead --------------------------------
# After sampling pair (i+1), adjust the cost matrix for pair (i) so that
# each survival edge (r -> j) gets a bonus for two-step biological
# plausibility: does the trajectory r -> j -> k (where k is j's
# forward assignment) have plausible growth?
#
# Arguments:
#   aug_cost   : K×K augmented log-cost matrix for pair i
#   n_curr     : number of real observations at census i
#   n_next     : number of real observations at census i+1
#   dbh_next   : DBH vector at census i+1 (length n_next)
#   next_assignment : K-length integer vector (assignment[row]=col at pair i+1)
#   n_next_next : number of real observations at census i+2
#   dbh_further : DBH vector at census i+2 (length n_next_next)
#   interval_next : interval in years between census i+1 and i+2
#   bio        : bio parameter list
#   weight     : lookahead weight (0 disables, 0.5 default)
#
# Returns: modified aug_cost matrix (same dimensions)

condition_cost_matrix <- function(aug_cost, n_curr, n_next,
                                  dbh_next, next_assignment,
                                  n_next_next, dbh_further,
                                  interval_next, bio, weight) {
    if (weight <= 0 || n_next == 0L || n_next_next == 0L) return(aug_cost)

    # Growth mean function (same as in compute_pairwise_log_likelihood)
    mu_growth_fn <- function(d) {
        if (!is.finite(bio$mu_gamma) || bio$mu_gamma == 0 ||
            !is.finite(d) || d <= 0)
            return(bio$mu_const)
        bio$mu_const + bio$mu_gamma * log(d)
    }

    # Pass 1: compute raw continuity log-LL for each real column j
    raw_bonus <- rep(0, n_next)  # 0 = neutral (no info)
    has_info  <- logical(n_next)
    for (j in seq_len(n_next)) {
        k <- next_assignment[j]
        if (k > n_next_next) next  # j died — no forward info

        d_j <- dbh_next[j]
        d_k <- dbh_further[k]
        if (!is.finite(d_j) || !is.finite(d_k) || d_j <= 0) next

        g2 <- (d_k - d_j) / interval_next
        sigma_j <- max(bio$sigma0 + bio$sigma1 * d_j, 1e-6)
        mu_j <- mu_growth_fn(d_j)
        raw_bonus[j] <- dnorm(g2, mean = mu_j, sd = sigma_j, log = TRUE)
        has_info[j] <- TRUE
    }

    # Normalize: shift so best column = 0, others get negative bonuses.
    # Columns without forward info (death / missing) get 0 (neutral).
    # Cap the maximum penalty at -2 log units to prevent small-DBH stems
    # (which have tight growth variance) from dominating the cost matrix.
    info_vals <- raw_bonus[has_info]
    if (length(info_vals) == 0L) return(aug_cost)  # nothing to condition on

    max_bonus <- max(info_vals)
    bonus <- rep(0, n_next)
    bonus[has_info] <- pmax(raw_bonus[has_info] - max_bonus, -2)

    # Pass 2: apply weighted normalized bonus to all feasible cells
    K <- nrow(aug_cost)
    for (j in seq_len(n_next)) {
        if (bonus[j] == 0) next  # no adjustment needed
        adj <- weight * bonus[j]
        for (r in seq_len(min(n_curr, K))) {
            if (is.finite(aug_cost[r, j])) {
                aug_cost[r, j] <- aug_cost[r, j] + adj
            }
        }
    }

    aug_cost
}

# ---- TrueStemID pin helpers for probabilistic matcher --------------------

# apply_pin_mask: For each pinned obs r at current census, find the column j
# at the next census that carries the pinned track, and set all other columns
# to -Inf so greedy_assignment_gumbel() is forced to pick column j.
#
# INPUTS
#   cost_matrix            K×K augmented cost matrix (modified in place)
#   pin_for_curr           integer vector length n_curr; pin_for_curr[r] =
#                          anchor position the obs is pinned to, or NA
#   next_obs_to_anchor_pos integer vector length n_next; which anchor position
#                          obs j at next census is currently carrying
#   n_curr, n_next         number of real observations at current / next census
#
# RETURNS  modified cost_matrix
apply_pin_mask <- function(cost_matrix, pin_for_curr, next_obs_to_anchor_pos,
                           n_curr, n_next) {
    for (r in seq_len(n_curr)) {
        target <- pin_for_curr[r]
        if (is.na(target)) next
        # Find column j at next census carrying this anchor position
        j_candidates <- which(next_obs_to_anchor_pos[seq_len(n_next)] == target)
        if (length(j_candidates) != 1L) next  # target died or ambiguous — skip
        j <- j_candidates[1L]
        # Mask all columns except j to -Inf for row r
        cost_matrix[r, -j] <- -Inf
    }
    cost_matrix
}

# propagate_track_backward: After an assignment is drawn, compute which
# anchor position each current-census obs now carries.
#
# INPUTS
#   assignment             K-length integer vector: assignment[r] = col
#   next_obs_to_anchor_pos integer vector: anchor position for each next-census obs
#   n_curr, n_next         number of real observations at current / next census
#
# RETURNS  integer vector length n_curr: anchor position per obs (NA = died)
propagate_track_backward <- function(assignment, next_obs_to_anchor_pos,
                                     n_curr, n_next) {
    curr <- rep(NA_integer_, n_curr)
    for (r in seq_len(n_curr)) {
        col <- assignment[r]
        if (col <= n_next) {
            curr[r] <- next_obs_to_anchor_pos[col]
        }
    }
    curr
}

# ---- Gumbel-noise greedy assignment --------------------------------------
# Draw one stochastic assignment from an augmented log-cost matrix using the
# Gumbel-max trick.  Each row is assigned to the highest-scoring available
# column after adding Gumbel(0, temperature) noise, processed in descending
# order of row maxima.
#
# INPUTS
#   log_cost_matrix  K×K matrix of log-likelihoods (augmented with
#                    mortality/recruitment slots so it is square).
#   temperature      Gumbel noise scale; higher = more random.
#
# RETURNS
#   Integer vector of length K: assignment[row] = assigned column index.

greedy_assignment_gumbel <- function(log_cost_matrix, temperature = 1.0) {
    K <- nrow(log_cost_matrix)
    stopifnot(ncol(log_cost_matrix) == K)

    # Add Gumbel(0, temperature) noise:  -temperature * log(-log(U))
    noise <- matrix(-temperature * log(-log(runif(K * K))), nrow = K, ncol = K)
    noisy <- log_cost_matrix + noise

    # Greedy assignment: for each row in descending max-noisy-score order,
    # assign to the best available column
    assignment <- integer(K)  # assignment[row] = col
    used_cols <- logical(K)

    # Order rows by their maximum noisy value (descending)
    row_max <- apply(noisy, 1, max, na.rm = TRUE)
    row_order <- order(row_max, decreasing = TRUE)

    for (r in row_order) {
        available <- which(!used_cols)
        if (length(available) == 0L) break
        scores <- noisy[r, available]
        best_idx <- available[which.max(scores)]
        assignment[r] <- best_idx
        used_cols[best_idx] <- TRUE
    }

    assignment
}

# ---- Repair stitched samples BEFORE marginal aggregation -----------------
# Walk each sample's trajectories and break links that violate growth
# constraints.  Two layers of defense (mirroring the DP pathway):
#
#   1. Hard-rate check: annualized growth outside [min_rate, max_rate]
#      is severed immediately (same as before).
#
#   2. ME-informed cumulative-shrinkage check: even when each consecutive
#      pair passes the hard rate, a long run of small decreases can
#      accumulate more shrinkage than measurement error can explain.
#      We track cumulative shrinkage along each trajectory and compare
#      it against an n_sigma_me threshold derived from the small-error
#      component of the BCI measurement-error model:
#          SD(D) = me_sd1_a * D + me_sd1_b
#      Threshold = n_sigma_me * sqrt( SD(d_start)^2 + SD(d_curr)^2 )
#      where d_start is the DBH at the beginning of the shrinkage run
#      and d_curr is the current DBH.
#      This mirrors DP's global cost accumulation which naturally penalises
#      consecutive shrinkage through likelihood, but adapted for the
#      per-pair greedy matcher.
#
# When a violation is found, the EARLIER observation is severed by assigning
# it a new unique break-ID.  Up to max_passes iterations per sample (a
# break can shorten trajectories and expose new violations).
#
# Returns the modified stitched list (same structure).

repair_stitched_growth_violations <- function(stitched, obs_data, intervals,
                                              min_rate, max_rate,
                                              me_sd1_a    = 0.0062,
                                              me_sd1_b    = 0.0904,
                                              n_sigma_me  = 3,
                                              max_passes  = 10L) {
    n_samples <- length(stitched)
    n_census  <- length(obs_data)
    if (n_census < 2L) return(stitched)

    # Global break-ID counter: start above the max ID in any sample to
    # avoid collisions when marginals aggregate across samples.
    break_base <- 0L
    for (s in seq_len(n_samples)) {
        for (ci in seq_len(n_census)) {
            ids <- stitched[[s]][[ci]]
            if (length(ids) > 0L) {
                mx <- max(ids, na.rm = TRUE)
                if (is.finite(mx) && mx > break_base) break_base <- mx
            }
        }
    }

    # ME helper: SD for the small-error component
    me_sd <- function(d) me_sd1_a * d + me_sd1_b

    total_breaks      <- 0L
    total_me_breaks   <- 0L

    for (s in seq_len(n_samples)) {
        for (pass in seq_len(max_passes)) {
            breaks_this_pass <- 0L

            # Build reverse map: stem_id -> list of (ci, oi, dbh)
            traj_map <- list()
            for (ci in seq_len(n_census)) {
                ids  <- stitched[[s]][[ci]]
                dbhs <- obs_data[[ci]]$dbh
                n_obs <- obs_data[[ci]]$n
                for (oi in seq_len(n_obs)) {
                    sid <- ids[oi]
                    key <- as.character(sid)
                    traj_map[[key]] <- c(traj_map[[key]], list(list(ci = ci, oi = oi, dbh = dbhs[oi])))
                }
            }

            # Walk each trajectory
            for (key in names(traj_map)) {
                entries <- traj_map[[key]]
                if (length(entries) < 2L) next

                # entries are already in census order (built ci=1..n_census)
                cumul_shrink   <- 0
                shrink_run_start <- 1L  # index into entries where current shrinkage run began
                d_run_start    <- entries[[1L]]$dbh  # DBH at start of shrinkage run

                for (r in 2:length(entries)) {
                    ci_prev <- entries[[r - 1L]]$ci
                    ci_curr <- entries[[r]]$ci

                    # Compute interval: sum of intervals between the two censuses
                    # (they may not be consecutive if obs are missing in between)
                    iv <- 0
                    if (ci_curr > ci_prev && ci_prev < n_census) {
                        for (ii in ci_prev:(ci_curr - 1L)) {
                            if (ii <= length(intervals)) iv <- iv + intervals[ii]
                        }
                    }
                    if (!is.finite(iv) || iv <= 0) iv <- 5.0

                    d_prev <- entries[[r - 1L]]$dbh
                    d_curr <- entries[[r]]$dbh
                    rate   <- (d_curr - d_prev) / iv

                    # --- Layer 1: hard-rate check --------------------------
                    if (rate < min_rate || rate > max_rate) {
                        break_base <- break_base + 1L
                        stitched[[s]][[entries[[r - 1L]]$ci]][entries[[r - 1L]]$oi] <- break_base
                        breaks_this_pass <- breaks_this_pass + 1L
                        break  # re-evaluate shortened trajectory in next pass
                    }

                    # --- Layer 2: ME cumulative-shrinkage check ------------
                    if (d_curr < d_prev) {
                        cumul_shrink <- cumul_shrink + (d_prev - d_curr)
                        thresh <- n_sigma_me * sqrt(me_sd(d_run_start)^2 + me_sd(d_curr)^2)
                        if (cumul_shrink > thresh) {
                            # Sever at the start of the shrinkage run
                            break_base <- break_base + 1L
                            stitched[[s]][[entries[[shrink_run_start]]$ci]][entries[[shrink_run_start]]$oi] <- break_base
                            breaks_this_pass <- breaks_this_pass + 1L
                            total_me_breaks  <- total_me_breaks + 1L
                            break  # re-evaluate shortened trajectory
                        }
                    } else {
                        # Growth step: reset cumulative shrinkage tracker
                        cumul_shrink     <- 0
                        shrink_run_start <- r
                        d_run_start      <- d_curr
                    }
                }
            }

            total_breaks <- total_breaks + breaks_this_pass
            if (breaks_this_pass == 0L) break  # this sample converged
        }
    }

    attr(stitched, "sample_level_breaks")    <- total_breaks
    attr(stitched, "sample_level_me_breaks") <- total_me_breaks
    stitched
}

# ---- Filter samples to pin-consistent ones -------------------------------
# After stitching + sample-level repair, discard any sample where a pinned
# observation ended up on the wrong track.  This guarantees that marginals
# computed from the surviving samples are conditioned on all pins being
# satisfied, so posteriors are naturally correct without post-hoc patches.
#
# INPUTS
#   stitched    list of n_samples; each element is a list of n_census
#               integer vectors (obs → track_id mapping)
#   pin_info    list of n_census; pin_info[[i]][j] = anchor-position
#               index (1..n_anchor) for obs j, or NA if not pinned
#   anchor_ids  integer vector of anchor IDs
#   n_census    number of censuses
#   min_keep    minimum number of surviving samples (safety net)
#   vcat        verbose logger
#   prefix      log prefix
#
# RETURNS  filtered stitched list (possibly unchanged if no pins or
#          too few survive).  attr("n_pin_filtered") records how many
#          were dropped.

filter_pin_consistent_samples <- function(stitched, pin_info, anchor_ids,
                                          n_census, min_keep = 10L,
                                          vcat = function(...) invisible(NULL),
                                          prefix = "") {
    n_samples <- length(stitched)
    if (n_samples == 0L) return(stitched)

    # Gather all (census, obs, expected_track_id) triples
    pin_checks <- list()
    for (i in seq_len(n_census)) {
        pi <- pin_info[[i]]
        if (is.null(pi)) next
        pinned_obs <- which(!is.na(pi))
        if (length(pinned_obs) == 0L) next
        for (j in pinned_obs) {
            pin_checks[[length(pin_checks) + 1L]] <- list(
                census = i, obs = j, expected = anchor_ids[pi[j]]
            )
        }
    }
    if (length(pin_checks) == 0L) {
        attr(stitched, "n_pin_filtered") <- 0L
        return(stitched)
    }

    # Check each sample
    keep <- logical(n_samples)
    for (s in seq_len(n_samples)) {
        ok <- TRUE
        for (pc in pin_checks) {
            actual <- stitched[[s]][[pc$census]][pc$obs]
            if (is.na(actual) || actual != pc$expected) {
                ok <- FALSE
                break
            }
        }
        keep[s] <- ok
    }

    n_kept <- sum(keep)
    n_dropped <- n_samples - n_kept

    if (n_dropped == 0L) {
        attr(stitched, "n_pin_filtered") <- 0L
        return(stitched)
    }

    # Safety net: if too few survive, warn and keep all
    .min_safe <- min(min_keep, max(1L, as.integer(n_samples / 4L)))
    if (n_kept < .min_safe) {
        .msg <- paste0(prefix, "WARNING: pin-consistent filter would keep only ",
                       n_kept, "/", n_samples, " samples (min_safe=", .min_safe,
                       "); keeping ALL samples (degrading to soft-pin behavior)")
        vcat(.msg)
        message(.msg)
        attr(stitched, "n_pin_filtered") <- 0L
        return(stitched)
    }

    .msg <- paste0(prefix, "Pin-consistent filter: kept ", n_kept, "/",
                   n_samples, " samples (", n_dropped, " dropped)")
    vcat(.msg)
    message(.msg)

    out <- stitched[keep]
    attr(out, "n_pin_filtered") <- n_dropped
    out
}

# ---- Stitch assignments backward from anchor ----------------------------

stitch_assignments_backward <- function(all_samples, obs_data, anchor_ids, K) {
    n_samples <- length(all_samples)
    n_census <- length(obs_data)
    anchor_pos <- n_census

    # Pre-compute deterministic death-track IDs for each (pair, obs).
    # When obs r at pair i is assigned to ANY death column, it always gets
    # death_ids[[i]][r] — the same ID across all samples.  This prevents
    # fragmentation caused by the Gumbel sampler picking different death
    # columns in different samples.
    death_base <- max(anchor_ids, na.rm = TRUE) + 1L
    death_ids <- vector("list", n_census - 1L)
    for (i in seq.int(anchor_pos - 1L, 1L, by = -1L)) {
        n_curr_i <- obs_data[[i]]$n
        death_ids[[i]] <- seq.int(death_base, length.out = n_curr_i)
        death_base <- death_base + n_curr_i
    }

    results <- vector("list", n_samples)

    for (s in seq_len(n_samples)) {
        sample_assignments <- all_samples[[s]]
        recon_by_census <- vector("list", n_census)

        # Anchor: obs positions map directly to anchor_ids
        n_anchor <- obs_data[[anchor_pos]]$n
        next_obs_to_track <- anchor_ids  # length n_anchor
        recon_by_census[[anchor_pos]] <- next_obs_to_track

        # Walk backward
        for (i in seq.int(anchor_pos - 1L, 1L, by = -1L)) {
            assignment <- sample_assignments[[i]]  # K_pair-length: assignment[row] = col
            n_curr <- obs_data[[i]]$n
            n_next <- obs_data[[i + 1L]]$n

            # Map current real observations to tracks
            curr_obs_to_track <- integer(n_curr)
            for (r in seq_len(n_curr)) {
                assigned_col <- assignment[r]
                if (assigned_col <= n_next) {
                    # Survival: inherit track from next census
                    curr_obs_to_track[r] <- next_obs_to_track[assigned_col]
                } else {
                    # Death: use deterministic ID for this (pair, obs)
                    curr_obs_to_track[r] <- death_ids[[i]][r]
                }
            }

            recon_by_census[[i]] <- curr_obs_to_track
            next_obs_to_track <- curr_obs_to_track
        }

        results[[s]] <- recon_by_census
    }

    results
}

# ---- Compute marginals from samples -------------------------------------

compute_marginals_from_samples <- function(stitched, tree_data, obs_data,
                                           obs_census, anchor_pos,
                                           posterior_top_k,
                                           intervals = NULL,
                                           min_rate = NULL,
                                           max_rate = NULL) {
    n_samples <- length(stitched)
    n_census <- length(obs_data)

    # Ensure posterior columns exist before writing
    for (k in seq_len(posterior_top_k)) {
        id_col <- paste0("DP_PosteriorTop", k, "ID")
        prob_col <- paste0("DP_PosteriorTop", k, "Prob")
        if (!(id_col %in% names(tree_data))) tree_data[, (id_col) := NA_integer_]
        if (!(prob_col %in% names(tree_data))) tree_data[, (prob_col) := NA_real_]
    }
    if (!("DP_PosteriorEntropy" %in% names(tree_data)))
        tree_data[, DP_PosteriorEntropy := NA_real_]
    if (!("DP_PosteriorReconstructedProb" %in% names(tree_data)))
        tree_data[, DP_PosteriorReconstructedProb := NA_real_]

    # ---- Pass 1: compute per-obs marginal posteriors ----
    all_posteriors <- vector("list", n_census)
    for (ci in seq_len(n_census)) {
        n_obs <- obs_data[[ci]]$n
        if (n_obs == 0L) next
        census_posts <- vector("list", n_obs)

        for (oi in seq_len(n_obs)) {
            assigned_ids <- vapply(stitched, function(s) {
                if (length(s[[ci]]) >= oi) s[[ci]][oi] else NA_integer_
            }, integer(1))

            id_table <- table(assigned_ids[!is.na(assigned_ids)])
            if (length(id_table) == 0L) {
                census_posts[[oi]] <- list(ids = integer(0), probs = numeric(0))
                next
            }
            sorted_ids <- sort(id_table, decreasing = TRUE)
            census_posts[[oi]] <- list(
                ids   = as.integer(names(sorted_ids)),
                probs = as.numeric(sorted_ids) / n_samples
            )
        }
        all_posteriors[[ci]] <- census_posts
    }

    # ---- Pass 2: growth-aware greedy resolver ---------------------------------
    # Resolves censuses from anchor outward so that each candidate ID can be
    # checked for growth-bound compatibility against the nearest already-resolved
    # assignment for that stem.  Marginal posteriors (Top-K, entropy) are pure
    # sample statistics and stay unchanged; only the ReconstructedStemID and its
    # associated DP_PosteriorReconstructedProb may differ from the naive MAP when
    # a growth-violating candidate is skipped.
    growth_aware <- !is.null(intervals) && !is.null(min_rate) && !is.null(max_rate)

    if (growth_aware) {
        # Census order: anchor first, then alternating ±1, ±2, ...
        anchor_out_order <- anchor_pos
        for (.d in seq_len(n_census - 1L)) {
            .before <- anchor_pos - .d
            .after  <- anchor_pos + .d
            if (.before >= 1L) anchor_out_order <- c(anchor_out_order, .before)
            if (.after <= n_census) anchor_out_order <- c(anchor_out_order, .after)
        }
        # Per-stem resolved track: stem_id_string -> list of (ci, dbh) entries
        .stem_tracks <- new.env(hash = TRUE, parent = emptyenv())

        # Interval between two census positions (handles gaps)
        .census_iv <- function(ci_a, ci_b) {
            lo <- min(ci_a, ci_b); hi <- max(ci_a, ci_b)
            if (lo == hi) return(0)
            idx_rng <- lo:(hi - 1L)
            idx_rng <- idx_rng[idx_rng <= length(intervals)]
            if (length(idx_rng) == 0L) return(5.0)
            iv <- sum(intervals[idx_rng])
            if (!is.finite(iv) || iv <= 0) 5.0 else iv
        }

        # Check growth rate against nearest resolved assignment on each side
        .growth_ok <- function(sid_key, ci_new, dbh_new) {
            trk <- .stem_tracks[[sid_key]]
            if (is.null(trk)) return(TRUE)
            trk_cis <- vapply(trk, function(x) x$ci, numeric(1))
            # Closest already-resolved census BEFORE ci_new
            below_idx <- which(trk_cis < ci_new)
            if (length(below_idx) > 0L) {
                j <- below_idx[which.max(trk_cis[below_idx])]
                iv <- .census_iv(trk_cis[j], ci_new)
                rate <- (dbh_new - trk[[j]]$dbh) / iv
                if (rate < min_rate || rate > max_rate) return(FALSE)
            }
            # Closest already-resolved census AFTER ci_new
            above_idx <- which(trk_cis > ci_new)
            if (length(above_idx) > 0L) {
                j <- above_idx[which.min(trk_cis[above_idx])]
                iv <- .census_iv(ci_new, trk_cis[j])
                rate <- (trk[[j]]$dbh - dbh_new) / iv
                if (rate < min_rate || rate > max_rate) return(FALSE)
            }
            TRUE
        }
    } else {
        anchor_out_order <- seq_len(n_census)
    }

    for (ci in anchor_out_order) {
        n_obs <- obs_data[[ci]]$n
        if (n_obs == 0L) next
        idx <- obs_data[[ci]]$idx
        census_posts <- all_posteriors[[ci]]

        # Greedy: sort obs by max posterior prob (descending)
        max_probs <- vapply(census_posts, function(p) {
            if (length(p$probs) > 0) p$probs[1] else 0
        }, numeric(1))
        order_by_conf <- order(max_probs, decreasing = TRUE)

        used_ids    <- integer(0)
        resolved_ids   <- integer(n_obs)
        resolved_probs <- numeric(n_obs)

        for (rank_pos in seq_along(order_by_conf)) {
            oi <- order_by_conf[rank_pos]
            post <- census_posts[[oi]]
            assigned <- FALSE
            for (j in seq_along(post$ids)) {
                cand_id <- post$ids[j]
                if (cand_id %in% used_ids) next
                # Growth-aware check: reject candidate if it violates growth
                # bounds against the nearest already-resolved census for this stem
                if (growth_aware) {
                    dbh_oi <- obs_data[[ci]]$dbh[oi]
                    if (!is.na(dbh_oi) &&
                        !.growth_ok(as.character(cand_id), ci, dbh_oi)) next
                }
                resolved_ids[oi]   <- cand_id
                resolved_probs[oi] <- post$probs[j]
                used_ids <- c(used_ids, cand_id)
                assigned <- TRUE
                break
            }
            if (!assigned) {
                # All posterior alternatives are taken or growth-violated —
                # assign a new unique (break) ID
                new_id <- max(c(tree_data$ReconstructedStemID, used_ids,
                                resolved_ids), na.rm = TRUE) + 1L
                if (!is.finite(new_id)) new_id <- 1L
                resolved_ids[oi]   <- new_id
                resolved_probs[oi] <- 0
                used_ids <- c(used_ids, new_id)
            }
        }

        # Register resolved assignments for growth tracking
        if (growth_aware) {
            for (oi in seq_len(n_obs)) {
                dbh_oi <- obs_data[[ci]]$dbh[oi]
                if (!is.na(dbh_oi)) {
                    sid_key <- as.character(resolved_ids[oi])
                    trk <- .stem_tracks[[sid_key]]
                    if (is.null(trk)) trk <- list()
                    trk[[length(trk) + 1L]] <- list(ci = ci, dbh = dbh_oi)
                    .stem_tracks[[sid_key]] <- trk
                }
            }
        }

        # Write resolved assignments + posteriors into tree_data
        for (oi in seq_len(n_obs)) {
            tree_data_row <- idx[oi]
            data.table::set(tree_data, tree_data_row, "ReconstructedStemID",
                            resolved_ids[oi])
            data.table::set(tree_data, tree_data_row, "DP_PosteriorReconstructedProb",
                            resolved_probs[oi])

            # Top-k posteriors (marginal, may differ from resolved assignment)
            post <- census_posts[[oi]]
            for (k in seq_len(min(posterior_top_k, length(post$ids)))) {
                id_col <- paste0("DP_PosteriorTop", k, "ID")
                prob_col <- paste0("DP_PosteriorTop", k, "Prob")
                data.table::set(tree_data, tree_data_row, id_col, post$ids[k])
                data.table::set(tree_data, tree_data_row, prob_col, post$probs[k])
            }

            # Entropy (from marginal posterior)
            if (length(post$probs) > 0) {
                ent <- -sum(post$probs * log(post$probs + 1e-30))
            } else {
                ent <- NA_real_
            }
            data.table::set(tree_data, tree_data_row, "DP_PosteriorEntropy", ent)
        }
    }

    # NA-DBH rows keep ReconstructedStemID = NA (no observation to match)

    tree_data
}

# ---- Repair growth violations from marginal resolution -------------------
# The per-census greedy marginal resolution can assign the same StemID to
# observations at consecutive censuses that violate growth bounds.  This
# happens because the marginals are computed independently per census.
#
# Repair strategy: walk each StemID's trajectory.  When a violation is
# found between census c and c+1, try to reassign the observation at
# census c to one of its top-k posterior alternatives.  If no alternative
# resolves the violation, break the track by assigning a new unique ID.

repair_marginal_growth_violations <- function(tree_data, obs_data, obs_census,
                                              intervals, min_rate, max_rate,
                                              posterior_top_k, vcat, prefix) {
    n_census <- length(obs_census)
    if (n_census < 2L) return(tree_data)

    # Build date-based interval lookup for arbitrary census pairs
    dt_dates <- tree_data[, .(MeanDate = mean(as.numeric(as.Date(ExactDate)),
                              na.rm = TRUE)), by = CensusID]

    max_id <- suppressWarnings(max(tree_data$ReconstructedStemID, na.rm = TRUE))
    if (!is.finite(max_id)) max_id <- 0L
    total_breaks <- 0L
    max_passes <- 10L  # safety limit

    for (pass in seq_len(max_passes)) {
        # Rebuild observation table each pass to reflect prior breaks
        obs_rows <- which(!is.na(tree_data$ReconstructedStemID) & !is.na(tree_data$DBH))
        obs_dt <- data.table::data.table(
            row_idx = obs_rows,
            CensusID = tree_data$CensusID[obs_rows],
            DBH = tree_data$DBH[obs_rows],
            ReconstructedStemID = tree_data$ReconstructedStemID[obs_rows]
        )
        setorder(obs_dt, ReconstructedStemID, CensusID)

        breaks_this_pass <- 0L
        for (sid in unique(obs_dt$ReconstructedStemID)) {
            dsub <- obs_dt[ReconstructedStemID == sid]
            if (nrow(dsub) < 2L) next

            for (r in 2:nrow(dsub)) {
                c_prev <- dsub$CensusID[r - 1L]
                c_curr <- dsub$CensusID[r]
                md0 <- dt_dates$MeanDate[dt_dates$CensusID == c_prev]
                md1 <- dt_dates$MeanDate[dt_dates$CensusID == c_curr]
                iv <- if (length(md0) > 0L && length(md1) > 0L)
                          (md1[1] - md0[1]) / 365.25
                      else 5.0
                if (!is.finite(iv) || iv <= 0) iv <- 5.0

                rate <- (dsub$DBH[r] - dsub$DBH[r - 1L]) / iv
                if (rate >= min_rate && rate <= max_rate) next

                # Violation: break the earlier census obs to a new unique ID
                row_prev <- dsub$row_idx[r - 1L]
                max_id <- max_id + 1L
                data.table::set(tree_data, as.integer(row_prev),
                                "ReconstructedStemID", as.integer(max_id))
                breaks_this_pass <- breaks_this_pass + 1L
                # After breaking, don't check further pairs for this stem
                # in this pass — re-evaluate in next pass
                break
            }
        }

        total_breaks <- total_breaks + breaks_this_pass
        if (breaks_this_pass == 0L) break  # converged
    }

    if (total_breaks > 0L) {
        .warn_msg <- paste0(prefix, "WARNING: Post-marginal safety-net repair fired ",
             total_breaks, " break(s) — probabilities for these rows are stale ",
             "(greedy conflict resolution created new violations)")
        vcat(.warn_msg)
        message(.warn_msg)  # ensure it appears on stderr / captured by log redirection
    }
    tree_data
}

# ---- Export posterior samples (mirrors DP format) -------------------------

export_probabilistic_posteriors <- function(stitched, tree_data, obs_data,
                                            obs_census, tag_val, n_samples,
                                            posterior_samples_path,
                                            posterior_samples_format,
                                            verbose, prefix, vcat) {
    fmt <- match.arg(posterior_samples_format, c("rds", "feather", "csv"))

    out_dir_local <- if (!is.null(posterior_samples_path)) {
        posterior_samples_path
    } else {
        get0("out_dir", ifnotfound = NULL)
    }
    if (is.null(out_dir_local) || !nzchar(out_dir_local)) out_dir_local <- tempdir()

    out_dir_post <- file.path(out_dir_local, "posteriors")
    if (!dir.exists(out_dir_post)) dir.create(out_dir_post, recursive = TRUE, showWarnings = FALSE)

    ts_local <- get0("BATCH_TS", ifnotfound = format(Sys.time(), "%Y%m%d_%H%M%S"))
    out_path_base <- file.path(out_dir_post,
                               paste0("tag_", ifelse(is.na(tag_val), "NA", tag_val),
                                      "_posterior_samples_", ts_local))

    n_census <- length(obs_data)

    # Build samples data.table matching DP format
    samples_list <- vector("list", n_samples)
    for (s in seq_len(n_samples)) {
        rows <- vector("list", n_census)
        for (ci in seq_len(n_census)) {
            n_obs <- obs_data[[ci]]$n
            if (n_obs == 0L) next
            recon_ids <- stitched[[s]][[ci]]
            rows[[ci]] <- data.table::data.table(
                Tag = tag_val,
                Sample = s,
                CensusID = rep(obs_census[ci], n_obs),
                ReconstructedStemID = recon_ids[seq_len(n_obs)],
                ObsRowID = obs_data[[ci]]$row_id
            )
        }
        samples_list[[s]] <- data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
    }
    samples_dt <- data.table::rbindlist(samples_list, use.names = TRUE, fill = TRUE)
    samples_dt <- samples_dt[order(Sample, CensusID)]

    # Path signatures (same as DP)
    sample_sigs <- samples_dt[, .(path_sig = paste0(ReconstructedStemID, collapse = "-")), by = Sample]
    path_counts <- sample_sigs[, .N, by = path_sig]
    data.table::setnames(path_counts, "N", "path_count")
    sample_sigs <- merge(sample_sigs, path_counts, by = "path_sig")

    n_unique <- data.table::uniqueN(sample_sigs$path_sig)
    vcat(prefix, "Posterior sampling: ", n_samples, " samples (", n_unique, " unique paths)")

    # Probabilities (uniform weights since greedy sampling doesn't have exact logp)
    paths_summary <- path_counts[, .(path_count, path_prob = path_count / n_samples), by = path_sig]

    # Compact reconstruction mapping
    recon_by_path <- samples_dt[, .(recon = paste0(ObsRowID, ":", ReconstructedStemID, collapse = ";")),
                                by = .(path_sig = sample_sigs$path_sig[match(Sample, sample_sigs$Sample)])]
    recon_compact <- recon_by_path[, .SD[1], by = path_sig, .SDcols = "recon"]
    paths_summary <- merge(paths_summary, recon_compact, by = "path_sig", all.x = TRUE)

    # Export
    if (fmt == "feather" && requireNamespace("arrow", quietly = TRUE)) {
        p2 <- paste0(out_path_base, "_paths.feather")
        arrow::write_feather(paths_summary, p2)
        vcat(prefix, "Wrote posterior paths summary to: ", p2)
    } else if (fmt == "csv") {
        p2 <- paste0(out_path_base, "_paths.csv")
        data.table::fwrite(paths_summary, p2)
        vcat(prefix, "Wrote posterior paths summary to: ", p2)
    } else {
        p1 <- paste0(out_path_base, "_paths.rds")
        saveRDS(paths_summary, file = p1)
        vcat(prefix, "Wrote posterior paths summary to: ", p1)
    }
}
