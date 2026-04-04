############################################################
# dp_probabilistic_matching.R
# Probabilistic greedy matching fallback for large state spaces
############################################################
# When the DP state space is too large (enum_exceeded or edge_count_exceeded),
# this module provides a stochastic per-census-pair matching approach that:
#   1. Computes pairwise log-likelihoods using the same biological model as DP
#   2. Augments the cost matrix with mortality/recruitment slots
#   3. Draws n_samples stochastic assignments via Gumbel-noise greedy
#   4. Stitches per-pair assignments backward from the anchor
#   5. Computes marginal posterior probabilities from samples
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
        recruit_max_dbh = unique(tree_data$Bio_Recruit_MaxDBH_unit)[1],
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
            # No DBH anywhere — cannot match
            tree_data[, ReconstructedStemID := seq_len(.N)]
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
    anchor_pos <- n_census  # anchor is the last observed census
    anchor_ids <- if (any(!is.na(anchor_obs$TrueStemID))) {
        as.integer(anchor_obs$TrueStemID[!is.na(anchor_obs$TrueStemID)])
    } else {
        seq_len(nrow(anchor_obs))
    }
    # Ensure all anchor obs have IDs
    tree_data[CensusID == anchor_start & !is.na(DBH), `:=`(
        ReconstructedStemID = anchor_ids,
        ReconstructionMethod = ifelse(!is.na(TrueStemID), "given", "probabilistic")
    )]

    # K = number of tracks (at least max obs across any census)
    max_obs <- max(vapply(obs_data, function(x) x$n, integer(1)))
    K <- max(length(anchor_ids), max_obs)

    vcat(prefix, "Probabilistic matching: ", n_census, " censuses, K=", K,
         ", max_obs=", max_obs, ", n_samples=", n_samples)

    # --- Per-pair log-likelihood matrices (backward) -----------------------
    # For each pair (c, c+1) compute the pairwise + augmented cost matrix
    pair_data <- vector("list", n_census - 1L)
    for (i in seq_len(n_census - 1L)) {
        dbh_curr <- obs_data[[i]]$dbh
        dbh_next <- obs_data[[i + 1L]]$dbh
        iv <- intervals[i]

        L <- compute_pairwise_log_likelihood(dbh_curr, dbh_next, iv, bio,
                                             min_growth, max_growth)
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
    for (s in seq_len(n_samples)) {
        # For each census pair (working backward from anchor-1 to 1),
        # sample an assignment
        per_pair_assignments <- vector("list", n_census - 1L)
        for (i in rev(seq_len(n_census - 1L))) {
            per_pair_assignments[[i]] <- greedy_assignment_gumbel(
                pair_data[[i]]$log_cost,
                temperature = temperature
            )
        }
        all_samples[[s]] <- per_pair_assignments
    }

    # --- Stitch assignments backward from anchor ---------------------------
    stitched <- stitch_assignments_backward(all_samples, obs_data, anchor_ids, K)

    # --- Compute marginals and fill tree_data ------------------------------
    tree_data <- compute_marginals_from_samples(stitched, tree_data, obs_data,
                                                obs_census, anchor_pos,
                                                posterior_top_k)

    tree_data[, ReconstructionMethod := ifelse(
        !is.na(TrueStemID) & ReconstructionMethod == "given",
        "given", "probabilistic"
    )]

    # --- Export posterior samples (same format as DP) -----------------------
    if (n_samples > 0L) {
        export_probabilistic_posteriors(
            stitched, tree_data, obs_data, obs_census,
            tag_val         = tag_val,
            n_samples       = n_samples,
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

    # Total size: (n_curr + n_recruit_slots) rows × (n_next + n_death_slots) cols
    # n_recruit_slots = max(0, n_next - n_curr) extra rows for new recruits
    # n_death_slots = max(0, n_curr - n_next) extra cols for dying tracks
    # But we always allow some flexibility: use K = max(n_curr, n_next)
    K <- max(n_curr, n_next)

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

# ---- Gumbel-noise greedy assignment --------------------------------------

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

# ---- Stitch assignments backward from anchor ----------------------------

stitch_assignments_backward <- function(all_samples, obs_data, anchor_ids, K) {
    n_samples <- length(all_samples)
    n_census <- length(obs_data)
    anchor_pos <- n_census

    results <- vector("list", n_samples)

    for (s in seq_len(n_samples)) {
        sample_assignments <- all_samples[[s]]
        recon_by_census <- vector("list", n_census)

        # Anchor: obs positions map directly to anchor_ids
        n_anchor <- obs_data[[anchor_pos]]$n
        next_obs_to_track <- anchor_ids  # length n_anchor
        recon_by_census[[anchor_pos]] <- next_obs_to_track

        # Global ID counter for new tracks (recruits at earlier censuses)
        next_new_id <- max(anchor_ids, na.rm = TRUE) + 1L

        # Walk backward
        for (i in seq.int(anchor_pos - 1L, 1L, by = -1L)) {
            assignment <- sample_assignments[[i]]  # K_pair-length: assignment[row] = col
            n_curr <- obs_data[[i]]$n
            n_next <- obs_data[[i + 1L]]$n
            K_pair <- length(assignment)

            # Build col_to_track for this pair:
            # cols 1..n_next -> known track IDs from next census
            # cols (n_next+1)..K_pair -> death sinks (stem at curr dies before next)
            col_to_track <- integer(K_pair)
            col_to_track[seq_len(min(n_next, K_pair))] <- next_obs_to_track[seq_len(min(n_next, K_pair))]
            if (K_pair > n_next) {
                for (jj in (n_next + 1L):K_pair) {
                    col_to_track[jj] <- next_new_id
                    next_new_id <- next_new_id + 1L
                }
            }

            # Map current real observations to tracks
            curr_obs_to_track <- integer(n_curr)
            for (r in seq_len(n_curr)) {
                curr_obs_to_track[r] <- col_to_track[assignment[r]]
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
                                           posterior_top_k) {
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

    # For each observation (census x position), count how often each
    # ReconstructedStemID was assigned across samples
    for (ci in seq_len(n_census)) {
        n_obs <- obs_data[[ci]]$n
        if (n_obs == 0L) next
        idx <- obs_data[[ci]]$idx

        for (oi in seq_len(n_obs)) {
            # Collect all assigned IDs for this observation across samples
            assigned_ids <- vapply(stitched, function(s) {
                if (length(s[[ci]]) >= oi) s[[ci]][oi] else NA_integer_
            }, integer(1))

            # Count frequencies
            id_table <- table(assigned_ids[!is.na(assigned_ids)])
            if (length(id_table) == 0L) next

            sorted_ids <- sort(id_table, decreasing = TRUE)
            probs <- as.numeric(sorted_ids) / n_samples
            ids <- as.integer(names(sorted_ids))

            # MAP assignment
            tree_data_row <- idx[oi]
            data.table::set(tree_data, tree_data_row, "ReconstructedStemID", ids[1])

            # Fill posterior columns
            for (k in seq_len(min(posterior_top_k, length(ids)))) {
                id_col <- paste0("DP_PosteriorTop", k, "ID")
                prob_col <- paste0("DP_PosteriorTop", k, "Prob")
                data.table::set(tree_data, tree_data_row, id_col, ids[k])
                data.table::set(tree_data, tree_data_row, prob_col, probs[k])
            }

            # Entropy
            ent <- -sum(probs * log(probs + 1e-30))
            data.table::set(tree_data, tree_data_row, "DP_PosteriorEntropy", ent)

            # Reconstructed prob = probability of the MAP assignment
            data.table::set(tree_data, tree_data_row, "DP_PosteriorReconstructedProb", probs[1])
        }
    }

    # Handle non-observed rows (NA DBH)
    na_rows <- which(is.na(tree_data$DBH) & is.na(tree_data$ReconstructedStemID))
    if (length(na_rows) > 0L) {
        current_max <- max(tree_data$ReconstructedStemID, na.rm = TRUE)
        if (!is.finite(current_max)) current_max <- 0L
        tree_data[na_rows, ReconstructedStemID := seq.int(current_max + 1L,
                                                          current_max + length(na_rows))]
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
