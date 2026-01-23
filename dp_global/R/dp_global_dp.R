############################################################
# dp_global_dp.R
# Core dynamic programming (MAP and marginal DP functions)
############################################################

# # Guard against accidental top-level side-effect when sourcing in tests
match_stems_dp_global_backward_marginals_batch <- function(tree_data,
                                                           min_growth = -Inf,
                                                           max_growth = Inf,
                                                           anchor_start,
                                                           max_tracks = 30L,
                                                           slack_tracks = 1L,
                                                           slack_require_anchor_recruitable = FALSE,
                                                           slack_require_anchor_eps = 1e-6,
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
                                                           # --- posterior sampling options ---
                                                           posterior_samples = 0L,
                                                           posterior_samples_format = c("rds", "feather", "csv"),
                                                           posterior_samples_path = NULL,
                                                           posterior_sample_seed = NULL,
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

    tree_data <- tree_data[order(CensusID)]

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
    # Posterior output columns (dynamically generated from `posterior_top_k`)
    post_cols <- c()
    for (k in seq_len(max(1L, as.integer(posterior_top_k)))) {
        post_cols <- c(post_cols, paste0("DP_PosteriorTop", k, "ID"), paste0("DP_PosteriorTop", k, "Prob"))
    }
    post_cols <- c(post_cols, "DP_PosteriorEntropy", "DP_PosteriorReconstructedProb", "DP_PosteriorUnlinkedProb")

    ensure_posterior_columns <- function(dt) {
        # data.table uses the type of the RHS to infer the column type.
        # Create ID (integer) and Prob (numeric) columns for 1..posterior_top_k
        for (k in seq_len(max(1L, as.integer(posterior_top_k)))) {
            id_col <- paste0("DP_PosteriorTop", k, "ID")
            prob_col <- paste0("DP_PosteriorTop", k, "Prob")
            if (!(id_col %in% names(dt))) dt[, (id_col) := NA_integer_]
            if (!(prob_col %in% names(dt))) dt[, (prob_col) := NA_real_]
            # Coerce types if present but wrong
            if (!(is.integer(dt[[id_col]]))) dt[, (id_col) := as.integer(get(id_col))]
            if (!(is.numeric(dt[[prob_col]]))) dt[, (prob_col) := as.numeric(get(prob_col))]
        }

        if (!("DP_PosteriorEntropy" %in% names(dt))) dt[, DP_PosteriorEntropy := NA_real_]
        if (!("DP_PosteriorReconstructedProb" %in% names(dt))) dt[, DP_PosteriorReconstructedProb := NA_real_]
        if (!("DP_PosteriorUnlinkedProb" %in% names(dt))) dt[, DP_PosteriorUnlinkedProb := NA_real_]

        if (!(is.numeric(dt$DP_PosteriorEntropy))) dt[, DP_PosteriorEntropy := as.numeric(DP_PosteriorEntropy)]
        if (!(is.numeric(dt$DP_PosteriorReconstructedProb))) dt[, DP_PosteriorReconstructedProb := as.numeric(DP_PosteriorReconstructedProb)]
        if (!(is.numeric(dt$DP_PosteriorUnlinkedProb))) dt[, DP_PosteriorUnlinkedProb := as.numeric(DP_PosteriorUnlinkedProb)]

        # Add a stable per-row identifier useful for posterior matching (preserve existing if present)
        if (!("obs_row_id" %in% names(dt))) {
            dt[, obs_row_id := seq_len(.N)]
        } else {
            # ensure integer type
            if (!is.integer(dt$obs_row_id)) dt[, obs_row_id := as.integer(obs_row_id)]
        }
        # Return the data.table (explicit)
        dt
    }

    # Ensure tree_data has posterior columns and obs_row_id assigned
    tree_data <- ensure_posterior_columns(tree_data)
    vcat(prefix, "Ensured posterior columns; obs_row_id present (or created)")

    # Preserve any provided TrueStemID as hard values in output.
    tree_data[!is.na(TrueStemID), `:=`(
        ReconstructedStemID = as.integer(TrueStemID),
        ReconstructionMethod = "given"
    )]

    # Observed-stem counts over census interval from first observed census up to anchor (for state-space diagnostics)
    # Determine the first census with any DBH observation at or before the anchor
    census_ids_up_to_anchor <- sort(unique(tree_data$CensusID[tree_data$CensusID <= anchor_start]))
    first_obs_census <- if (length(census_ids_up_to_anchor) > 0L) {
        tmp <- tree_data[CensusID <= anchor_start & !is.na(DBH), CensusID]
        if (length(tmp) > 0L) min(tmp) else NA_integer_
    } else {
        NA_integer_
    }
    if (is.na(first_obs_census)) {
        vcat(prefix, "Cannot find any observations up to anchor_start. Falling back to igraph.")
        K_used <- as.integer(min(0L, max_tracks))
        tree_data[, `:=`(
            DP_KUsed = K_used,
            DP_MaxStatesPerCensus = 0L,
            DP_MaxStatesCensusID = NA_integer_
        )]
        out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, anchor_start)
        out <- ensure_posterior_columns(out)
        return(out)
    }
    census_range <- seq.int(from = first_obs_census, to = anchor_start)
    n_census <- length(census_range)
    vcat(prefix, "first_obs_census=", first_obs_census, "census_range=", paste(census_range, collapse = ","))
    obs_counts <- vapply(
        census_range,
        function(cc) nrow(tree_data[CensusID == cc & !is.na(DBH)]),
        integer(1L)
    )
    max_obs <- if (length(obs_counts) > 0L) max(obs_counts) else 0L

    # Need a fully-anchored endpoint
    # If the requested anchor census is missing entirely or all rows at that census
    # have both NA DBH and NA TrueStemID, prefer the most recent earlier census
    # that has at least one row with non-NA DBH and non-NA TrueStemID and use that
    # as the anchor instead of immediately falling back to the igraph matcher.
    anchor_rows_all <- tree_data[CensusID == anchor_start]
    if (nrow(anchor_rows_all) == 0L || (all(is.na(anchor_rows_all$DBH)) && all(is.na(anchor_rows_all$TrueStemID)))) {
        cand_census <- sort(unique(tree_data$CensusID[tree_data$CensusID < anchor_start & !is.na(tree_data$DBH) & !is.na(tree_data$TrueStemID)]))
        if (length(cand_census) > 0L) {
            new_anchor <- as.integer(max(cand_census))
            vcat(prefix, "Requested anchor census=", anchor_start, " had no DBH/TrueStemID; using earlier anchor census=", new_anchor)
            anchor_start <- new_anchor
            census_range <- seq.int(from = first_obs_census, to = anchor_start)
            n_census <- length(census_range)
            vcat(prefix, "first_obs_census=", first_obs_census, "census_range=", paste(census_range, collapse = ","))
            obs_counts <- vapply(
                census_range,
                function(cc) nrow(tree_data[CensusID == cc & !is.na(DBH)]),
                integer(1L)
            )
            max_obs <- if (length(obs_counts) > 0L) max(obs_counts) else 0L
        } else {
            vcat(prefix, "Cannot anchor DP (missing anchor observations or TrueStemID). Falling back to igraph.")
            K_used <- as.integer(min(max_obs, max_tracks))
            n_states_by_census <- vapply(obs_counts, function(n_obs) count_injective_states(K_used, n_obs), numeric(1L))
            tree_data[, `:=`(
                DP_KUsed = K_used,
                DP_MaxStatesPerCensus = max(n_states_by_census, na.rm = TRUE),
                DP_MaxStatesCensusID = as.integer(census_range[which.max(n_states_by_census)])
            )]
            out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, anchor_start)
            out <- ensure_posterior_columns(out)
            return(out)
        }
    }

    anchor_obs <- tree_data[CensusID == anchor_start & !is.na(DBH)]
    if (nrow(anchor_obs) == 0L || any(is.na(anchor_obs$TrueStemID))) {
        vcat(prefix, "Cannot anchor DP (missing anchor observations or TrueStemID). Falling back to igraph.")
        K_used <- as.integer(min(max_obs, max_tracks))
        n_states_by_census <- vapply(obs_counts, function(n_obs) count_injective_states(K_used, n_obs), numeric(1L))
        tree_data[, `:=`(
            DP_KUsed = K_used,
            DP_MaxStatesPerCensus = max(n_states_by_census, na.rm = TRUE),
            DP_MaxStatesCensusID = as.integer(census_range[which.max(n_states_by_census)])
        )]
        out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, anchor_start)
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
            DP_MaxStatesCensusID = as.integer(census_range[which.max(n_states_by_census)])
        )]
        out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, anchor_start)
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

    # Optionally require that an anchor DBH be recruitable before granting slack.
    anchor_ok <- TRUE
    if (isTRUE(slack_require_anchor_recruitable)) {
        eps <- suppressWarnings(as.numeric(slack_require_anchor_eps))
        if (!is.finite(eps) || is.na(eps)) eps <- 0
        anchor_ok <- any(
            !is.na(anchor_obs$DBH) &
                !is.na(anchor_obs$Bio_Recruit_MaxDBH_unit) &
                (anchor_obs$DBH <= (anchor_obs$Bio_Recruit_MaxDBH_unit + eps))
        )
        if (!anchor_ok) {
            vcat(prefix, "slack requested but not granted: no anchor DBH <= recruit_max_dbh (eps=", eps, ")")
        }
    }

    # Apply slack only when requested and when anchor_ok (if required).
    K_target <- K_base + ifelse(slack_tracks > 0L && isTRUE(anchor_ok), slack_tracks, 0L)
    K <- min(K_target, max_tracks)

    # Report theoretical worst-case state count
    n_states_by_census <- vapply(obs_counts, function(n_obs) count_injective_states(K, n_obs), numeric(1L))
    tree_data[, `:=`(
        DP_KUsed = as.integer(K),
        DP_MaxStatesPerCensus = max(n_states_by_census, na.rm = TRUE),
        DP_MaxStatesCensusID = as.integer(census_range[which.max(n_states_by_census)])
    )]

    if (K < max(obs_counts)) {
        vcat(prefix, "K too small for observed counts (K=", K, ", max_obs=", max(obs_counts), "). Falling back to igraph.")
        out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, anchor_start)
        out <- ensure_posterior_columns(out)
        return(out)
    }

    vcat(prefix, "Chosen K=", K, " tracks; max theoretical states=", format(max(n_states_by_census, na.rm = TRUE), scientific = TRUE))

    n_extra <- K - length(anchor_ids)
    current_max <- suppressWarnings(max(tree_data$TrueStemID, na.rm = TRUE))
    if (!is.finite(current_max)) current_max <- 0
    track_ids <- c(anchor_ids, if (n_extra > 0L) seq.int(from = current_max + 1L, length.out = n_extra) else integer(0))

    # Pre-enumerate assignment states (injective obs->track) for each census in census_range
    obs_dbh <- vector("list", n_census)
    state_mats <- vector("list", n_census)
    state_keys <- vector("list", n_census)
    for (p in seq_len(n_census)) {
        cc <- census_range[p]
        obs <- tree_data[CensusID == cc & !is.na(DBH)]
        obs_dbh[[p]] <- obs$DBH
        n_obs <- length(obs_dbh[[p]])
        mat <- enumerate_states_injective(K, n_obs, max_states = max_states)

        if (is.null(mat)) {
            vcat(prefix, "State enumeration exceeded max_states at CensusID=", cc, " (n_obs=", n_obs, "). Falling back to igraph.")
            out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, anchor_start)
            out <- ensure_posterior_columns(out)
            return(out)
        }
        state_mats[[p]] <- mat
        state_keys[[p]] <- apply(mat, 1L, state_key)

        vcat(prefix, "Enumerated CensusID=", cc, ": n_obs=", n_obs, ", n_states=", nrow(mat))
    }

    # Anchor state assignment (pins endpoint)
    anchor_obs_ordered <- tree_data[CensusID == anchor_start & !is.na(DBH)]
    anchor_track_idx <- match(anchor_obs_ordered$TrueStemID, track_ids)

    if (any(is.na(anchor_track_idx))) {
        vcat(prefix, "Anchor TrueStemID not found in track_ids. Falling back to igraph.")
        out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, anchor_start)
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
    track_dbh_by_state <- vector("list", n_census)
    for (p in seq_len(n_census)) {
        mat <- state_mats[[p]]
        n_states <- nrow(mat)
        tdbh_list <- vector("list", n_states)
        for (i in seq_len(n_states)) {
            tdbh_list[[i]] <- state_to_track_dbh(mat[i, ], obs_dbh[[p]], K)
        }
        track_dbh_by_state[[p]] <- tdbh_list
    }

    # Backward tables (sum-product) and Viterbi backpointers
    keys_full <- vector("list", n_census)
    assign_full <- vector("list", n_census)
    logB <- vector("list", n_census)
    vit_cost <- vector("list", n_census)
    vit_ptr <- vector("list", n_census)
    edges <- vector("list", n_census)

    # Initialize at anchor census: only one allowed full-state (use position index)
    anchor_pos <- which(census_range == anchor_start)
    phase_anchor <- rep.int(2L, K)
    phase_anchor[anchor_track_idx] <- 1L
    anchor_full_key <- encode_full_key(anchor_track_idx, phase_anchor)
    keys_full[[anchor_pos]] <- anchor_full_key
    assign_full[[anchor_pos]] <- list(as.integer(anchor_track_idx))
    logB[[anchor_pos]] <- 0
    vit_cost[[anchor_pos]] <- 0
    vit_ptr[[anchor_pos]] <- integer(0)
    edges[[anchor_pos]] <- NULL

    # Backward recursion p = anchor_pos-1 .. 1 (maps to CensusID via census_range)
    vcat(prefix, "Backward pass (log-sum-exp + Viterbi) starting ...")
    for (p in seq.int(anchor_pos - 1L, 1L, by = -1L)) {
        # p: position in census_range (for testing only)
        cc <- census_range[p]
        next_cc <- census_range[p + 1L]
        mat_cc <- state_mats[[p]]
        n_states_cc <- nrow(mat_cc)

        t_cc0 <- tic()
        vcat(prefix, "Backward step CensusID=", cc, ": n_assignment_states=", n_states_cc, ", n_next_full_states=", length(keys_full[[p + 1L]]))

        next_keys <- keys_full[[p + 1L]]
        n_next <- length(next_keys)
        if (n_next == 0L) {
            vcat(prefix, "No reachable next full-states at CensusID=", next_cc, ". Falling back to igraph.")
            out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, anchor_start)
            out <- ensure_posterior_columns(out)
            return(out)
        }
        next_index <- seq_len(n_next)
        names(next_index) <- next_keys
        logB_next <- as.numeric(logB[[p + 1L]])
        vit_next <- as.numeric(vit_cost[[p + 1L]])

        # Fast lookup for next assignment DBHs via assignment key
        next_assign_list <- assign_full[[p + 1L]]
        # Build assignment-key -> state index for next census (since phase differs but assignment cost uses assignment)
        next_assign_key <- vapply(next_assign_list, state_key, character(1L))
        next_assign_row_idx <- match(next_assign_key, state_keys[[p + 1L]])
        if (any(is.na(next_assign_row_idx))) {
            # Should not happen; indicates mismatch in state enumeration.
            out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, anchor_start)
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
            tdbh1_by_next[[j]] <- track_dbh_by_state[[p + 1L]][[next_assign_row_idx[[j]]]]
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

        resolve_interval_years_pair <- function(tree_data) {
            dt <- tree_data[, .(CensusID, ExactDate)]
            ## get mean exactdate per census
            dt_mean <- dt[, .(MeanDate = mean(ExactDate, na.rm = TRUE)), by = CensusID]
            setorder(dt_mean, CensusID)

            ## dcast to wide format
            dt_wide <- dcast(dt_mean, 1 ~ CensusID, value.var = "MeanDate")
            ## compute interval between t0 and t1
            return(dt_wide)
        }

        pair_interval <- resolve_interval_years_pair(tree_data)

        for (i in seq_len(n_states_cc)) {
            # i <- 1L # For testing only
            tdbh0 <- track_dbh_by_state[[p]][[i]]
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

            interval_val <- (as.numeric(pair_interval[[as.character(next_cc)]]) - as.numeric(pair_interval[[as.character(cc)]])) / 365.25
            # print(interval_val)
            # Batch compute all transition costs from this current assignment
            c_trans_vec <- transition_cost_tracks_bio_batch_rcpp(
                track_dbh_t = tdbh0,
                track_dbh_tp1 = feasible_tdbh1,
                interval_years = interval_val,
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
            out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, anchor_start)
            out <- ensure_posterior_columns(out)
            return(out)
        }

        keys_full[[p]] <- unlist(curr_keys_list, use.names = FALSE)
        assign_full[[p]] <- curr_assign_list
        logB[[p]] <- curr_logB
        vit_cost[[p]] <- curr_vit
        vit_ptr[[p]] <- curr_ptr

        if (used_edges == 0L) {
            vcat(prefix, "No feasible edges at CensusID=", cc, ". Falling back to igraph.")
            out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, anchor_start)
            out <- ensure_posterior_columns(out)
            return(out)
        }
        edges[[p]] <- data.table::data.table(
            from_idx = from_idx[seq_len(used_edges)],
            to_idx = to_idx[seq_len(used_edges)],
            logw = logw[seq_len(used_edges)]
        )

        vcat(prefix, "Finished CensusID=", cc, ": full_states=", length(keys_full[[p]]), ", edges=", used_edges, ", dt=", sprintf("%.2fs", tic() - t_cc0))
    }

    # -----------------
    # Decode MAP path
    # -----------------
    vcat(prefix, "Decoding MAP path and writing ReconstructedStemID ...")
    start_idx <- which.min(vit_cost[[1L]])
    if (length(start_idx) == 0L || !is.finite(vit_cost[[1L]][start_idx])) {
        out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, anchor_start)
        out <- ensure_posterior_columns(out)
        return(out)
    }
    map_idx <- integer(n_census)
    map_idx[1L] <- start_idx
    for (p in seq_len(n_census - 1L)) {
        nxt <- vit_ptr[[p]][map_idx[p]]
        if (!is.finite(nxt) || is.na(nxt) || nxt < 1L) {
            out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, anchor_start)
            out <- ensure_posterior_columns(out)
            return(out)
        }
        map_idx[p + 1L] <- nxt
    }

    for (p in seq_len(n_census)) {
        cc <- census_range[p]
        obs_idx <- tree_data[CensusID == cc & !is.na(DBH), which = TRUE]
        if (length(obs_idx) == 0L) next
        sv <- assign_full[[p]][[map_idx[p]]]
        if (length(sv) != length(obs_idx)) {
            out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, anchor_start)
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
    logalpha <- vector("list", n_census)
    logalpha[[1L]] <- rep.int(0, length(keys_full[[1L]]))

    for (p in seq_len(n_census - 1L)) {
        ed <- edges[[p]]
        if (is.null(ed) || nrow(ed) == 0L) {
            out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, anchor_start)
            out <- ensure_posterior_columns(out)
            return(out)
        }
        la_from <- logalpha[[p]][ed$from_idx]
        vals <- la_from + ed$logw
        dt <- data.table::data.table(to_idx = ed$to_idx, v = vals)
        dt <- dt[is.finite(v)]
        if (nrow(dt) == 0L) {
            out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, anchor_start)
            out <- ensure_posterior_columns(out)
            return(out)
        }
        la_next_dt <- dt[, .(logalpha = log_sum_exp(v)), by = to_idx]
        la_next <- rep.int(-Inf, length(keys_full[[p + 1L]]))
        la_next[la_next_dt$to_idx] <- la_next_dt$logalpha
        logalpha[[p + 1L]] <- la_next

        vcat(prefix, "Forward step CensusID=", census_range[p + 1L], ": reached ", sum(is.finite(la_next)), " / ", length(la_next), " states")
    }

    # Partition function Z = total weight of all paths ending at the fixed anchor state.
    # At anchor_start there is exactly one state.
    logZ <- logalpha[[anchor_pos]][1L]
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

    for (p in seq_len(n_census)) {
        cc <- census_range[p]
        obs_idx <- tree_data[CensusID == cc & !is.na(DBH), which = TRUE]
        if (length(obs_idx) == 0L) next

        # State posterior weights at this census
        lg <- logalpha[[p]] + logB[[p]] - logZ
        # Normalize defensively
        lg <- lg - log_sum_exp(lg)
        w <- exp(lg)

        n_obs <- length(obs_idx)
        prob_mat <- matrix(0, nrow = n_obs, ncol = K)
        st_assign <- assign_full[[p]]
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
        map_assign <- assign_full[[p]][[map_idx[p]]]

        # Prepare containers for top-K posterior summaries
        Ktop <- max(1L, as.integer(posterior_top_k))
        top_id_mat <- matrix(NA_integer_, nrow = n_obs, ncol = Ktop)
        top_p_mat <- matrix(NA_real_, nrow = n_obs, ncol = Ktop)
        entropy <- numeric(n_obs)
        p_map <- numeric(n_obs)
        p_unlinked <- numeric(n_obs)

        for (j in seq_len(n_obs)) {
            pvec <- prob_mat[j, ]
            # numeric stability
            pvec[pvec < 0] <- 0
            sp <- sum(pvec)
            if (sp > 0) pvec <- pvec / sp
            ord <- order(pvec, decreasing = TRUE)
            m <- min(length(ord), Ktop)
            if (m >= 1L) {
                for (kk in seq_len(m)) {
                    idxk <- ord[kk]
                    top_id_mat[j, kk] <- track_ids[idxk]
                    top_p_mat[j, kk] <- pvec[idxk]
                }
            }
            entropy[j] <- -sum(ifelse(pvec > 0, pvec * log(pvec), 0))
            p_map[j] <- pvec[map_assign[j]]
            p_unlinked[j] <- sum(pvec[!is_anchor_track])
        }

        # Build named list for data.table assignment
        assign_list <- list()
        for (kk in seq_len(Ktop)) {
            id_col <- paste0("DP_PosteriorTop", kk, "ID")
            prob_col <- paste0("DP_PosteriorTop", kk, "Prob")
            assign_list[[id_col]] <- top_id_mat[, kk]
            assign_list[[prob_col]] <- top_p_mat[, kk]
        }
        assign_list[["DP_PosteriorEntropy"]] <- entropy
        assign_list[["DP_PosteriorReconstructedProb"]] <- p_map
        assign_list[["DP_PosteriorUnlinkedProb"]] <- p_unlinked

        tree_data[obs_idx, (names(assign_list)) := assign_list]
    }

    tree_data <- add_constraint_violation(
        tree_data,
        id_col = "ReconstructedStemID",
        min_growth = min_growth,
        max_growth = max_growth,
        pair_interval = pair_interval
    )

    # ------------------------------
    # Optional: posterior sampling of full reconstructions
    # ------------------------------
    posterior_samples <- as.integer(posterior_samples)
    if (!is.null(posterior_samples) && posterior_samples > 0L) {
        fmt <- match.arg(posterior_samples_format)
        out_dir_local <- if (!is.null(posterior_samples_path)) posterior_samples_path else get0("out_dir", ifnotfound = NULL)
        if (is.null(out_dir_local) || (is.character(out_dir_local) && nzchar(out_dir_local) == FALSE)) {
            out_dir_local <- getwd()
        }
        if (is.null(out_dir_local) || !nzchar(out_dir_local)) out_dir_local <- tempdir()

        # Fallback for Tag if not present in tree_data
        tag_local <- if (!is.na(tag_val)) tag_val else get0("which_tag", ifnotfound = NA_integer_)

        # Build adjacency lookup for quick sampling
        adj_by_p <- vector("list", n_census - 1L)
        for (p in seq_len(n_census - 1L)) {
            ed <- edges[[p]]
            # rows grouped by to_idx
            to_idxs <- unique(ed$to_idx)
            m <- vector("list", length(ed$to_idx))
            lookup <- vector("list", max(to_idxs))
            for (r in seq_len(nrow(ed))) {
                j <- ed$to_idx[r]
                if (is.null(lookup[[j]])) lookup[[j]] <- integer(0)
                lookup[[j]] <- c(lookup[[j]], r)
            }
            adj_by_p[[p]] <- lookup
        }

        # sampler: backwards sampling from anchor_pos to census 1
        if (!is.null(posterior_sample_seed)) set.seed(as.integer(posterior_sample_seed))
        samples_list <- vector("list", posterior_samples)
        for (m in seq_len(posterior_samples)) {
            sampled_idx <- integer(n_census)
            sampled_idx[anchor_pos] <- 1L # anchor position has single full-state at index 1
            logp <- 0
            # sample backwards
            for (p in seq.int(anchor_pos - 1L, 1L, by = -1L)) {
                j <- sampled_idx[p + 1L]
                rows <- adj_by_p[[p]][[j]]
                # rows are indices into edges[[p]]
                from_idx <- edges[[p]]$from_idx[rows]
                logw <- edges[[p]]$logw[rows]
                loga <- logalpha[[p]][from_idx]
                L <- loga + logw
                # normalize to probabilities
                Lmax <- max(L)
                probs <- exp(L - Lmax)
                probs <- probs / sum(probs)
                k <- sample.int(length(probs), size = 1L, prob = probs)
                chosen_row <- rows[k]
                sampled_idx[p] <- edges[[p]]$from_idx[chosen_row]
                logp <- logp + log(probs[k])
            }
            # Convert sampled full-state indices to ReconstructedStemID per census
            sample_dt <- data.table::data.table(Tag = tag_local, Sample = m, CensusID = integer(0), ReconstructedStemID = integer(0), ObsRowID = integer(0))
            for (p in seq_len(n_census)) {
                cc <- census_range[p]
                assign_vec <- assign_full[[p]][[sampled_idx[p]]]
                track_ids_loc <- track_ids[assign_vec]
                # attach as multiple rows (one per observed tree)
                obs_idx <- tree_data[CensusID == cc & !is.na(DBH), which = TRUE]
                if (length(obs_idx) > 0) {
                    obs_row_ids <- tree_data$obs_row_id[obs_idx]
                    sample_dt <- rbind(sample_dt, data.table::data.table(Tag = tag_local, Sample = m, CensusID = rep(cc, length(obs_idx)), ReconstructedStemID = track_ids_loc, ObsRowID = obs_row_ids))
                }
            }
            sample_dt[, logp := logp]
            samples_list[[m]] <- sample_dt
        }

        samples_dt <- data.table::rbindlist(samples_list, use.names = TRUE, fill = TRUE)
        # Ensure rows are ordered for signature construction
        samples_dt <- samples_dt[order(Sample, CensusID)]

        # Collate per-sample path signature and counts (useful diagnostics)
        sample_sigs <- samples_dt[, .(path_sig = paste0(ReconstructedStemID, collapse = "-")), by = Sample]
        path_counts <- sample_sigs[, .N, by = path_sig]
        sample_sigs <- merge(sample_sigs, path_counts, by = "path_sig")
        setnames(sample_sigs, "N", "path_count")
        samples_dt <- merge(samples_dt, sample_sigs, by = "Sample", all.x = TRUE)

        # Normalize sample weights from logp using per-Sample unique logp
        sample_logp <- unique(samples_dt[, .(Sample, logp)])
        maxlp <- max(sample_logp$logp, na.rm = TRUE)
        sample_logp[, sample_weight := exp(logp - maxlp)]
        sample_logp[, sample_prob := sample_weight / sum(sample_weight)]
        samples_dt <- merge(samples_dt, sample_logp[, .(Sample, sample_weight, sample_prob)], by = "Sample", all.x = TRUE)

        # Diagnostic summary
        n_unique_paths <- uniqueN(sample_sigs$path_sig)
        vcat(prefix, sprintf("Posterior sampling: generated %d samples (%d unique paths)", posterior_samples, n_unique_paths))
        if (uniqueN(sample_logp$logp) == 1L) {
            vcat(prefix, sprintf("Warning: all sampled logp identical (logp=%f); posterior may be degenerate or sampling explored equivalent paths", sample_logp$logp[1]))
        }

        # Export samples: prefer feather via arrow for speed if available
        ts_local <- get0("BATCH_TS", ifnotfound = format(Sys.time(), "%Y%m%d_%H%M%S"))
        out_dir_post <- file.path(out_dir_local, "posteriors")
        if (!dir.exists(out_dir_post)) dir.create(out_dir_post, recursive = TRUE, showWarnings = FALSE)
        out_path_base <- file.path(out_dir_post, paste0("tag_", ifelse(is.na(tag_local), "NA", tag_local), "_posterior_samples", "_", ts_local))

        # Prepare summary tables useful for downstream uncertainty propagation
        # Attach per-sample ObsRowID list when available (semicolon-separated string of ObsRowIDs in sample order)
        if ("ObsRowID" %in% names(samples_dt)) {
            sample_obsids <- samples_dt[, .(ObsRowIDs = paste0(ObsRowID, collapse = ";")), by = Sample]
        } else {
            sample_obsids <- data.table::data.table(Sample = integer(0), ObsRowIDs = character(0))
        }
        samples_summary <- unique(samples_dt[, .(Sample, Tag, logp, path_sig, path_count, sample_weight, sample_prob)])
        if (nrow(sample_obsids) > 0) samples_summary <- merge(samples_summary, sample_obsids, by = "Sample", all.x = TRUE)
        # per-path aggregated probabilities
        paths_summary <- samples_summary[, .(path_count = sum(path_count, na.rm = TRUE), path_prob = sum(sample_prob, na.rm = TRUE)), by = path_sig]
        # also create a compact per-path reconstruction mapping (one row per path)
        if ("ObsRowID" %in% names(samples_dt)) {
            recon_by_path <- samples_dt[, .(recon = paste0(ObsRowID, ":", ReconstructedStemID, collapse = ";")), by = .(path_sig, Sample)]
        } else {
            recon_by_path <- samples_dt[, .(recon = paste0(CensusID, ":", ReconstructedStemID, collapse = ";")), by = .(path_sig, Sample)]
        }
        # take first recon per path_sig (they are identical across samples with same path_sig)
        recon_compact <- recon_by_path[, .SD[1], by = path_sig, .SDcols = "recon"]
        paths_summary <- merge(paths_summary, recon_compact, by = "path_sig", all.x = TRUE)

        # Export samples and summaries: prefer feather via arrow for speed if available
        if (fmt == "feather" && requireNamespace("arrow", quietly = TRUE)) {
            arrow::write_feather(samples_dt, paste0(out_path_base, ".feather"))
            arrow::write_feather(samples_summary, paste0(out_path_base, "_summary.feather"))
            arrow::write_feather(paths_summary, paste0(out_path_base, "_paths.feather"))
            vcat(prefix, "Wrote posterior samples to: ", paste0(out_path_base, ".feather"))
            vcat(prefix, "Wrote posterior samples summary to: ", paste0(out_path_base, "_summary.feather"))
            vcat(prefix, "Wrote posterior paths summary to: ", paste0(out_path_base, "_paths.feather"))
        } else if (fmt == "csv") {
            data.table::fwrite(samples_dt, paste0(out_path_base, ".csv"))
            data.table::fwrite(samples_summary, paste0(out_path_base, "_summary.csv"))
            data.table::fwrite(paths_summary, paste0(out_path_base, "_paths.csv"))
            vcat(prefix, "Wrote posterior samples to: ", paste0(out_path_base, ".csv"))
            vcat(prefix, "Wrote posterior samples summary to: ", paste0(out_path_base, "_summary.csv"))
            vcat(prefix, "Wrote posterior paths summary to: ", paste0(out_path_base, "_paths.csv"))
        } else {
            saveRDS(list(full = samples_dt, summary = samples_summary, paths = paths_summary), file = paste0(out_path_base, ".rds"))
            vcat(prefix, "Wrote posterior samples to: ", paste0(out_path_base, ".rds"))
        }
    }

    vcat(prefix, "Done. Total elapsed ", sprintf("%.2fs", tic() - t_start))
    return(tree_data)
}

# Summary of what's available in posteriors

# Full long-format samples file (one row per observed reconstructed tree per sample):
# Example: posteriors/tag_20_posterior_samples_<ts>.csv
# Per-sample summary (one row per sampled full-reconstruction):
# posteriors/tag_20_posterior_samples_<ts>_summary.csv
# Columns: Sample, Tag, logp, path_sig, path_count, sample_weight, sample_prob
# sample_prob is normalized over all drawn samples and can be used as sampling weight for downstream propagation.
# Per-path aggregated summary (unique reconstructions):
# posteriors/tag_20_posterior_samples_<ts>_paths.csv
# Columns: path_sig, path_count, path_prob, recon (compact reconstruction mapping like "1:8;2:3;3:4;...")

# How to use these for error propagation (suggestions)
# Use paths.csv directly: each row is a unique reconstruction with probability path_prob (sums to 1 across unique paths) — convenient for expectation of downstream metrics without resampling.
# Or sample reconstructions according to sample_prob in the summary file to create Monte Carlo realizations for error propagation; then expand each sample using the full long file if needed to attach per-census reconstructed IDs.
# The recon column in paths.csv is handy to quickly apply a mapping (parse "CensusID:ReconstructedStemID" pairs).
