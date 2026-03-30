############################################################
### estimate_dp_complexity_function.R
### Accurate estimator of DP transition cost calls per tag
###
### Mirrors exact logic from dp_global_dp.R to estimate —
### without running the DP — how many transition cost evaluations
### will be performed per tag.  Use this to rank tags by expected
### run time before submitting batch jobs.
###
### Key mechanics taken directly from dp_global_dp.R:
###  1. census_range  — from first_obs_census up to anchor_start
###                     (shifted when last observed census < anchor_start)
###  2. K (tracks)    — max(anchor_ids, max_obs, K_from_counts)
###                     + resprout barrier + slack, capped to max_tracks
###  3. State space   — P(K, n_obs) per census;
###                     exceeding max_states flags estimated_fallback
###  4. Pruned edges  — enumerate every (state_t × state_{t+1}) pair
###                     and apply the same hard-prune guard:
###                       DBH->DBH  : (d1-d0)/interval in [eff_min, eff_max]
###                       NA ->DBH  : d1 <= eff_recruit_max
###                       DBH->NA   : always OK (death)
###                       NA ->NA   : always OK (unborn track)
###  5. Output        — rows per tag, sorted by estimated_edges_pruned desc
###                     + predicted_seconds / predicted_hours from calibrated model
############################################################

#' Estimate DP transition-cost evaluations per tag
#'
#' @param data RDS path, CSV path, or already-loaded data.table/data.frame
#' @param anchor_start Integer: DP anchor census (default 7)
#' @param slack_tracks Integer: extra tracks for simultaneous death+birth (default 1)
#' @param max_states Integer: same as DP max_states – tags exceeding this per
#'   census would fall back to igraph and are flagged (default 40000)
#' @param max_tracks Integer: hard cap on K (default 9999)
#' @param slack_require_anchor_recruitable Logical: only grant slack when any
#'   anchor DBH <= recruit_max_dbh (default FALSE)
#' @param min_growth Numeric: DP min_growth (default -Inf)
#' @param max_growth Numeric: DP max_growth (default Inf)
#' @param prune_min_growth Numeric|NULL: overrides min_growth in pruning
#' @param prune_max_growth Numeric|NULL: overrides max_growth in pruning
#' @param prune_use_bio_bounds Logical: intersect user bounds with Bio columns (default TRUE)
#' @param bio_max_shrink Numeric: fallback Bio_Max_Shrink when column absent (default -Inf)
#' @param bio_max_growth Numeric: fallback Bio_Max_Growth when column absent (default Inf)
#' @param recruit_max_dbh Numeric: fallback recruit size cap when column absent (default Inf)
#' @param prune_recruit_max_dbh Numeric|NULL: explicit recruit prune cap
#' @param prune_use_bio_recruit Logical: use bio recruit max DBH (default TRUE)
#' @param fast Logical: when TRUE (default) skip the O(states^2) pruned-edge
#'   enumeration and rank by unpruned edge count instead.  Use fast=FALSE only
#'   when you need accurate pruned counts (e.g., for model validation).
#' @return data.table sorted descending by estimated_edges_unpruned (fast=TRUE)
#'   or estimated_edges_pruned (fast=FALSE), with columns `predicted_seconds`
#'   and `predicted_hours` from a log-log polynomial model calibrated on actual
#'   benchmark runs (see dp_global/dev/run_dp_benchmark.R)
#' @export
estimate_dp_complexity <- function(data,
                                   anchor_start = 7L,
                                   slack_tracks = 1L,
                                   max_states = 40000L,
                                   max_tracks = 9999L,
                                   slack_require_anchor_recruitable = FALSE,
                                   min_growth = -Inf,
                                   max_growth = Inf,
                                   prune_min_growth = NULL,
                                   prune_max_growth = NULL,
                                   prune_use_bio_bounds = TRUE,
                                   bio_max_shrink = -Inf,
                                   bio_max_growth = Inf,
                                   recruit_max_dbh = Inf,
                                   prune_recruit_max_dbh = NULL,
                                   prune_use_bio_recruit = TRUE,
                                   fast = TRUE) {
    if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table required")
    library(data.table)

    # ------------------------------------------------------------------
    # 1. Load data
    # ------------------------------------------------------------------
    if (is.character(data)) {
        if (grepl("\\.rds$", data, ignore.case = TRUE)) {
            data <- readRDS(data)
            setDT(data)
        } else {
            data <- fread(data)
        }
    } else {
        setDT(copy(data))
    }

    # Normalise required columns (case-insensitive)
    for (req in c("Tag", "CensusID", "DBH", "TrueStemID")) {
        if (!(req %in% names(data))) {
            m <- names(data)[tolower(names(data)) == tolower(req)]
            if (length(m) == 1L) setnames(data, m, req) else stop("Missing column: ", req)
        }
    }

    species_col <- if ("species" %in% names(data)) {
        "species"
    } else if ("Species" %in% names(data)) {
        "Species"
    } else {
        data[, species := "unknown"]
        "species"
    }

    # ------------------------------------------------------------------
    # 2. Per-census-pair interval lookup (mirrors resolve_interval_years_pair)
    # ------------------------------------------------------------------
    census_interval_lookup <- NULL
    if ("ExactDate" %in% names(data)) {
        dates <- unique(data[, .(CensusID, ExactDate)])
        dates[, ExactDate := as.Date(as.character(ExactDate))]
        dates <- dates[!is.na(ExactDate)][order(CensusID)]
        dmean <- dates[, .(MeanDate = mean(as.numeric(ExactDate), na.rm = TRUE)), by = CensusID]
        setorder(dmean, CensusID)
        if (nrow(dmean) >= 2L) {
            census_interval_lookup <- data.table(
                c0             = dmean$CensusID[-nrow(dmean)],
                c1             = dmean$CensusID[-1L],
                interval_years = diff(dmean$MeanDate) / 365.25
            )
        }
    }

    get_interval <- function(c0, c1) {
        if (!is.null(census_interval_lookup)) {
            r <- census_interval_lookup[.(c0, c1), on = .(c0, c1)]
            if (nrow(r) == 1L && is.finite(r$interval_years) && r$interval_years > 0) {
                return(r$interval_years)
            }
        }
        NA_real_
    }

    # ------------------------------------------------------------------
    # 3. Helpers mirroring dp_global_states.R
    # ------------------------------------------------------------------

    # P(K, n) = K * (K-1) * ... * (K-n+1)
    count_injective_states <- function(K, n_obs) {
        K <- as.integer(K)
        n_obs <- as.integer(n_obs)
        if (n_obs == 0L) {
            return(1.0)
        }
        if (n_obs > K || K <= 0L) {
            return(0.0)
        }
        prod(seq.int(from = K, to = K - n_obs + 1L, by = -1L))
    }

    # Enumerate all injective assignments; returns NULL when > max_states
    enumerate_assignments <- function(K, n_obs, max_s) {
        if (n_obs == 0L) {
            return(matrix(integer(0), nrow = 1L, ncol = 0L))
        }
        if (n_obs > K) {
            return(NULL)
        }
        if (count_injective_states(K, n_obs) > max_s) {
            return(NULL)
        }
        build <- function(pfx, rem) {
            if (length(pfx) == n_obs) {
                return(list(pfx))
            }
            out <- vector("list", 0L)
            for (t in rem) out <- c(out, build(c(pfx, t), rem[rem != t]))
            out
        }
        combos <- build(integer(0L), seq_len(K))
        do.call(rbind, lapply(combos, function(v) matrix(v, nrow = 1L)))
    }

    # Map assignment row -> K-length track DBH vector (NA = unoccupied)
    state_to_track_dbh <- function(assign_vec, dbh_vec, K) {
        out <- rep(NA_real_, K)
        if (length(assign_vec) > 0L) out[assign_vec] <- dbh_vec
        out
    }

    # ------------------------------------------------------------------
    # 4. Hard-prune check — exact mirror of candidate_ok in dp_global_dp.R
    # ------------------------------------------------------------------
    transition_feasible <- function(tdbh0, tdbh1, interval_val, eff_min, eff_max, eff_rec) {
        if (!is.finite(interval_val)) {
            return(TRUE)
        }
        for (k in seq_along(tdbh0)) {
            v0 <- tdbh0[k]
            v1 <- tdbh1[k]
            if (!is.na(v0) && !is.na(v1)) {
                g <- (v1 - v0) / interval_val
                if (g < eff_min || g > eff_max) {
                    return(FALSE)
                }
            } else if (is.na(v0) && !is.na(v1)) {
                if (isTRUE(is.finite(eff_rec) && v1 > eff_rec)) {
                    return(FALSE)
                }
            }
            # DBH->NA and NA->NA always feasible
        }
        TRUE
    }

    # ------------------------------------------------------------------
    # 4b. Vectorized pruned-edge count for one census transition.
    #     Uses outer() matrix operations instead of nested R for loops.
    #     Memory: O(n0 * n1) per track — feasible up to ~6000 states.
    # ------------------------------------------------------------------
    count_pruned_edges_step_vec <- function(mat0, mat1, dbh0, dbh1, K,
                                           interval_val, eff_min, eff_max, eff_rec) {
        n0 <- nrow(mat0)
        n1 <- nrow(mat1)
        if (n0 == 0L || n1 == 0L) return(0L)

        # Build track-DBH matrices: (n_states x K), NA where track is unoccupied
        tdbh0 <- matrix(NA_real_, nrow = n0, ncol = K)
        tdbh1 <- matrix(NA_real_, nrow = n1, ncol = K)
        for (i in seq_len(n0)) {
            a <- mat0[i, ]
            if (length(a) > 0L) tdbh0[i, a] <- dbh0
        }
        for (j in seq_len(n1)) {
            a <- mat1[j, ]
            if (length(a) > 0L) tdbh1[j, a] <- dbh1
        }

        # feasible[i,j] starts all TRUE; we AND-in per-track constraints
        feasible <- matrix(TRUE, nrow = n0, ncol = n1)

        if (!is.finite(interval_val)) return(sum(feasible))

        for (k in seq_len(K)) {
            d0 <- tdbh0[, k]   # length n0
            d1 <- tdbh1[, k]   # length n1
            have0 <- !is.na(d0)
            have1 <- !is.na(d1)

            # Case 1: both observed -> growth must be in [eff_min, eff_max]
            both <- outer(have1, have0, "&")  # n1 x n0 (R outer: rows=d1, cols=d0)
            if (any(both)) {
                g <- outer(d1, d0, "-") / interval_val  # n1 x n0
                bad <- both & (g < eff_min | g > eff_max)
                feasible <- feasible & !t(bad)   # transpose to n0 x n1
            }

            # Case 2: recruit (NA at t, observed at t+1) -> DBH at t+1 <= eff_rec
            if (is.finite(eff_rec)) {
                recruit <- outer(have1, !have0, "&")  # n1 x n0
                bad_rec <- recruit & outer(d1 > eff_rec, rep(TRUE, length(d0)), "&")
                feasible <- feasible & !t(bad_rec)
            }
            # Cases 3,4 (death, both-NA) are always feasible — no update needed
        }
        sum(feasible)
    }

    # ------------------------------------------------------------------
    # 5. Per-tag estimator — mirrors match_stems_dp_global_backward_marginals_batch
    # ------------------------------------------------------------------
    estimate_one_tag <- function(tag_data, tag_val, species_val) {
        # 5a. Effective anchor (shift down if no obs reach anchor_start)
        obs_census_all <- sort(unique(tag_data$CensusID[!is.na(tag_data$DBH)]))
        last_obs <- if (length(obs_census_all) > 0L) max(obs_census_all) else NA_integer_
        eff_anchor <- if (!is.na(last_obs) && last_obs < anchor_start) {
            as.integer(last_obs)
        } else {
            as.integer(anchor_start)
        }

        # 5b. census_range: first obs census up to eff_anchor
        # Use pre-aggregated summary (keyed) for O(1) lookup per tag.
        tag_census_obs <- obs_summary_pre[.(tag_val)][CensusID <= eff_anchor]
        if (nrow(tag_census_obs) == 0L) {
            return(data.table(
                Tag = tag_val, Species = species_val,
                census_range = NA_character_, n_censuses = 0L, K = 0L,
                max_obs = 0L, max_states_per_census = 0,
                total_states = 0, estimated_edges_unpruned = 0,
                estimated_edges_pruned = NA_real_,
                estimated_fallback = FALSE,
                fallback_reason = "no_obs_up_to_anchor"
            ))
        }
        first_obs <- as.integer(min(tag_census_obs$CensusID))
        # Keep all observed censuses up to anchor; add anchor itself if no obs there
        census_range <- sort(unique(c(
            tag_census_obs$CensusID,
            if (eff_anchor > max(tag_census_obs$CensusID)) eff_anchor else integer(0L)
        )))
        census_range <- census_range[census_range <= eff_anchor]
        n_census <- length(census_range)

        obs_counts <- tag_census_obs$N[match(census_range, tag_census_obs$CensusID)]
        obs_counts[is.na(obs_counts)] <- 0L
        max_obs <- if (length(obs_counts) > 0L) max(obs_counts) else 0L

        # 5c. K — exact replica of dp_global_dp.R K computation
        anchor_obs <- tag_data[CensusID == eff_anchor & !is.na(DBH)]
        anchor_ids <- sort(unique(na.omit(anchor_obs$TrueStemID)))

        if (length(anchor_ids) == 0L) {
            return(data.table(
                Tag = tag_val, Species = species_val,
                census_range = paste(census_range, collapse = ","),
                n_censuses = n_census, K = 0L, max_obs = max_obs,
                max_states_per_census = 0, total_states = 0,
                estimated_edges_unpruned = 0,
                estimated_edges_pruned = NA_real_,
                estimated_fallback = FALSE,
                fallback_reason = "anchor_ids_missing"
            ))
        }

        births_needed <- if (length(obs_counts) >= 2L) sum(pmax(0L, diff(obs_counts))) else 0L
        K_from_counts <- as.integer(if (length(obs_counts) > 0L) obs_counts[1L] + births_needed else 0L)
        K_base <- max(length(anchor_ids), max_obs, K_from_counts)

        # Resprout barrier — use pre-aggregated summary (keyed lookup)
        n_resprout <- 0L
        if (!is.null(resprout_summary_pre)) {
            tag_resp <- resprout_summary_pre[.(tag_val)][CensusID %in% census_range]
            n_resprout <- if (nrow(tag_resp) > 0L) sum(tag_resp$N) else 0L
        }
        K_base <- K_base + n_resprout

        # Recruit max for slack eligibility guard
        tag_recruit_max <- recruit_max_dbh
        for (col in c("Bio_Recruit_MaxDBH_unit", "Bio_recruit_maxdbh_unit")) {
            if (col %in% names(anchor_obs)) {
                v <- unique(na.omit(anchor_obs[[col]]))
                if (length(v) == 1L && is.finite(v)) {
                    tag_recruit_max <- v
                    break
                }
            }
        }
        grant_slack <- slack_tracks > 0L
        if (isTRUE(slack_require_anchor_recruitable) && is.finite(tag_recruit_max)) {
            grant_slack <- grant_slack &&
                any(!is.na(anchor_obs$DBH) & anchor_obs$DBH <= tag_recruit_max)
        }
        K_target <- K_base + ifelse(grant_slack, as.integer(slack_tracks), 0L)
        K <- min(K_target, as.integer(max_tracks))

        if (K < max_obs) {
            return(data.table(
                Tag = tag_val, Species = species_val,
                census_range = paste(census_range, collapse = ","),
                n_censuses = n_census, K = K, max_obs = max_obs,
                max_states_per_census = 0, total_states = 0,
                estimated_edges_unpruned = 0,
                estimated_edges_pruned = NA_real_,
                estimated_fallback = TRUE,
                fallback_reason = "K_too_small"
            ))
        }

        # 5d. State counts per census
        n_states_c <- vapply(obs_counts, function(n) count_injective_states(K, n), numeric(1L))
        max_s_c <- max(n_states_c, na.rm = TRUE)
        total_s <- sum(n_states_c, na.rm = TRUE)
        fallback <- isTRUE(max_s_c > max_states)

        # 5e. Effective pruning bounds — mirrors dp_global_dp.R lines ~978-999
        read_bio_col <- function(candidates, default_val) {
            for (col in candidates) {
                if (col %in% names(tag_data)) {
                    v <- unique(na.omit(tag_data[[col]]))
                    if (length(v) == 1L && is.finite(v)) {
                        return(v)
                    }
                }
            }
            default_val
        }
        tag_bio_shrink <- if (isTRUE(prune_use_bio_bounds)) {
            read_bio_col(
                c("Bio_Max_Shrink_unit", "Bio_max_shrink_unit", "Bio_Max_Shrink"),
                bio_max_shrink
            )
        } else {
            bio_max_shrink
        }
        tag_bio_growth <- if (isTRUE(prune_use_bio_bounds)) {
            read_bio_col(
                c("Bio_Max_Growth_unit", "Bio_max_growth_unit", "Bio_Max_Growth"),
                bio_max_growth
            )
        } else {
            bio_max_growth
        }

        user_min <- if (!is.null(prune_min_growth)) prune_min_growth else min_growth
        user_max <- if (!is.null(prune_max_growth)) prune_max_growth else max_growth
        eff_min <- if (isTRUE(prune_use_bio_bounds)) max(user_min, tag_bio_shrink) else user_min
        eff_max <- if (isTRUE(prune_use_bio_bounds)) min(user_max, tag_bio_growth) else user_max

        eff_rec <- if (!is.null(prune_recruit_max_dbh)) {
            if (isTRUE(prune_use_bio_recruit)) {
                min(tag_recruit_max, prune_recruit_max_dbh)
            } else {
                prune_recruit_max_dbh
            }
        } else {
            tag_recruit_max
        }

        # 5f. Unpruned edge count (product of adjacent state counts)
        edges_unpruned <- 0
        if (n_census >= 2L) {
            edges_unpruned <- sum(n_states_c[-n_census] * n_states_c[-1L])
        }

        # 5g. Pruned edge count: full state enumeration + per-pair feasibility check
        #     Skipped entirely when fast=TRUE (use unpruned count for ranking).
        #     Skipped (NA) when tag would fall back to igraph.
        edges_pruned <- NA_real_
        if (!isTRUE(fast) && !fallback && n_census >= 2L) {
            all_mats <- vector("list", n_census)
            ok <- TRUE
            for (p in seq_len(n_census)) {
                all_mats[[p]] <- enumerate_assignments(K, obs_counts[p], max_states)
                if (is.null(all_mats[[p]])) {
                    ok <- FALSE
                    break
                }
            }
            if (ok) {
                edges_pruned <- 0
                for (p in seq_len(n_census - 1L)) {
                    mat0 <- all_mats[[p]]
                    mat1 <- all_mats[[p + 1L]]
                    dbh0 <- tag_data[CensusID == census_range[p] & !is.na(DBH), DBH]
                    dbh1 <- tag_data[CensusID == census_range[p + 1L] & !is.na(DBH), DBH]
                    iv <- get_interval(census_range[p], census_range[p + 1L])
                    step <- 0L
                    for (i in seq_len(nrow(mat0))) {
                        td0 <- state_to_track_dbh(mat0[i, ], dbh0, K)
                        for (j in seq_len(nrow(mat1))) {
                            td1 <- state_to_track_dbh(mat1[j, ], dbh1, K)
                            if (transition_feasible(td0, td1, iv, eff_min, eff_max, eff_rec)) {
                                step <- step + 1L
                            }
                        }
                    }
                    edges_pruned <- edges_pruned + step
                }
            } else {
                edges_pruned <- edges_unpruned # conservative upper bound
            }
        }

        data.table(
            Tag                      = tag_val,
            Species                  = species_val,
            census_range             = paste(census_range, collapse = ","),
            n_censuses               = n_census,
            K                        = K,
            max_obs                  = max_obs,
            max_states_per_census    = max_s_c,
            total_states             = total_s,
            estimated_edges_unpruned = edges_unpruned,
            estimated_edges_pruned   = edges_pruned,
            estimated_fallback       = fallback,
            fallback_reason          = if (fallback) "enum_exceeded" else NA_character_
        )
    }

    # ------------------------------------------------------------------
    # 6. Runtime predictor
    #    Calibrated polynomial in log-log space from dp_global/dev/run_dp_benchmark.R.
    #    Uses estimated_edges_pruned when available (= actual C++ calls after pruning),
    #    otherwise estimated_edges_unpruned (upper bound, same as old TransitionComputations).
    #    Note: the model was fitted on unpruned counts; pruned predictions are optimistic.
    # ------------------------------------------------------------------
    predict_seconds <- function(N) {
        logN <- log10(pmax(N, 1))  # guard against 0
        10^(-1.45817 + (-0.07313 * logN) + (0.14164 * logN^2))
    }

    # ------------------------------------------------------------------
    # 7. Run over all tags and sort
    # ------------------------------------------------------------------
    # FAST PATH (fast=TRUE, default): fully vectorized data.table operations,
    # no per-tag for loop.  Scales to 100k+ tags in seconds.
    # SLOW PATH (fast=FALSE): per-tag for loop used only when pruned-edge
    # counts are needed (exact complexity measurement).
    # ------------------------------------------------------------------
    if (isTRUE(fast)) {
        cat("[estimate_dp_complexity] Pre-aggregating (vectorized)...\n"); flush.console()
        setkey(data, Tag, CensusID)

        # 7a. Observation counts per (Tag, CensusID)
        obs_sum <- data[!is.na(DBH), .N, by = .(Tag, CensusID)]
        setkey(obs_sum, Tag, CensusID)

        # 7b. Per-tag last observed census -> effective anchor
        tag_meta <- obs_sum[, .(last_obs = max(CensusID)), by = Tag]
        tag_meta[, eff_anchor := as.integer(pmin(last_obs, as.integer(anchor_start)))]
        setkey(tag_meta, Tag)

        # 7c. Filter obs to <= eff_anchor per tag
        obs_pre <- obs_sum[tag_meta[, .(Tag, eff_anchor)], on = "Tag"][CensusID <= eff_anchor]
        setorder(obs_pre, Tag, CensusID)

        # 7d. Per-tag census stats + K_from_counts (births-based)
        tag_stats <- obs_pre[, {
            n <- N
            list(
                n_censuses    = .N,
                max_obs       = max(n),
                K_from_counts = n[1L] + sum(pmax(0L, diff(n)))
            )
        }, by = Tag]

        # 7e. Anchor ids at effective anchor (one join, not per-tag)
        data_anch <- data[!is.na(DBH)][
            tag_meta[, .(Tag, CensusID = eff_anchor)],
            on = .(Tag, CensusID), nomatch = NULL
        ]
        anchor_stats <- data_anch[,
            .(n_anchor_ids = as.integer(length(unique(na.omit(TrueStemID))))),
            by = Tag
        ]

        # 7f. Resprouts (one grep pass on full data, then aggregate per tag)
        resprout_regex <- "\\b(R|RP|RF|RT|QR)\\b"
        if ("ListOfTSM" %in% names(data)) {
            resp_sum <- data[
                !is.na(DBH) & !is.na(ListOfTSM) & grepl(resprout_regex, ListOfTSM, perl = TRUE),
                .N, by = .(Tag, CensusID)
            ]
            resp_pre <- resp_sum[tag_meta[, .(Tag, eff_anchor)], on = "Tag"][CensusID <= eff_anchor]
            resp_per_tag <- resp_pre[, .(n_resprout = sum(N)), by = Tag]
        } else {
            resp_per_tag <- data.table(Tag = integer(0L), n_resprout = integer(0L))
        }

        # 7g. Species per tag
        if (species_col %in% names(data)) {
            sp_tab <- data[, .(Species = as.character(.SD[[1L]][1L])), by = Tag, .SDcols = species_col]
        } else {
            sp_tab <- data.table(Tag = unique(data$Tag), Species = NA_character_)
        }

        # 7h. Assemble master per-tag table
        out <- tag_meta[tag_stats, on = "Tag"]
        out <- anchor_stats[out, on = "Tag"]
        out[is.na(n_anchor_ids), n_anchor_ids := 0L]
        out <- resp_per_tag[out, on = "Tag"]
        out[is.na(n_resprout), n_resprout := 0L]
        out <- sp_tab[out, on = "Tag"]

        # census_range string (cosmetic)
        census_range_str <- obs_pre[, .(census_range = paste(CensusID, collapse = ",")), by = Tag]
        out <- census_range_str[out, on = "Tag"]
        out[is.na(census_range), census_range := NA_character_]

        # 7i. K computation (vectorized)
        out[, K_base := pmax(n_anchor_ids, max_obs, K_from_counts) + n_resprout]

        if (isTRUE(slack_require_anchor_recruitable) && is.finite(recruit_max_dbh)) {
            anch_slack <- data_anch[!is.na(DBH),
                .(grant_slack = any(DBH <= recruit_max_dbh, na.rm = TRUE)), by = Tag]
            out <- anch_slack[out, on = "Tag"]
            out[is.na(grant_slack), grant_slack := FALSE]
        } else {
            out[, grant_slack := (slack_tracks > 0L)]
        }
        out[, K := as.integer(pmin(
            K_base + ifelse(grant_slack, as.integer(slack_tracks), 0L),
            as.integer(max_tracks)
        ))]

        # 7j. State counts: precompute all unique (K, n) pairs once, then join
        obs_for_states <- out[, .(Tag, K)][obs_pre, on = "Tag", nomatch = NULL]
        unique_kn <- unique(obs_for_states[, .(K, N)])
        unique_kn[, n_states := mapply(function(Kv, nv) {
            Kv <- as.integer(Kv); nv <- as.integer(nv)
            if (nv == 0L || Kv <= 0L) return(1.0)
            if (nv > Kv) return(0.0)
            prod(seq.int(from = Kv, to = Kv - nv + 1L, by = -1L))
        }, K, N)]
        obs_for_states <- unique_kn[obs_for_states, on = .(K, N)]
        setorder(obs_for_states, Tag, CensusID)

        # 7k. Per-tag: max states, total states, unpruned edge count
        census_stats <- obs_for_states[, .(
            max_states_per_census    = max(n_states),
            total_states             = sum(n_states),
            estimated_edges_unpruned = if (.N >= 2L) sum(n_states[-.N] * n_states[-1L]) else 0
        ), by = Tag]
        out <- census_stats[out, on = "Tag"]

        # 7l. Fallback flags
        out[, estimated_fallback := (n_anchor_ids == 0L) |
                (K < max_obs) |
                (max_states_per_census > as.numeric(max_states))]
        out[, fallback_reason := fcase(
            n_anchor_ids == 0L,                                    "anchor_ids_missing",
            K < max_obs,                                           "K_too_small",
            max_states_per_census > as.numeric(max_states),        "enum_exceeded",
            default                                                = NA_character_
        )]
        out[is.na(estimated_edges_unpruned), estimated_edges_unpruned := 0]
        out[, estimated_edges_pruned := NA_real_]
        out[, n_censuses := as.integer(n_censuses)]

        # 7m. Pruned-edge computation (vectorized matrix ops per tag)
        #     Only for non-fallback tags with >= 2 censuses.
        #     Uses count_pruned_edges_step_vec() which does outer() products
        #     per track — memory bounded by O(max_states^2).
        prune_candidates <- out[estimated_fallback == FALSE & n_censuses >= 2L, Tag]
        if (length(prune_candidates) > 0L) {
            # Effective bounds (same logic as slow path; bio_bounds not used here when prune_use_bio_bounds=FALSE)
            user_min_g <- if (!is.null(prune_min_growth)) prune_min_growth else min_growth
            user_max_g <- if (!is.null(prune_max_growth)) prune_max_growth else max_growth
            eff_min_g  <- if (isTRUE(prune_use_bio_bounds)) max(user_min_g, bio_max_shrink) else user_min_g
            eff_max_g  <- if (isTRUE(prune_use_bio_bounds)) min(user_max_g, bio_max_growth) else user_max_g
            eff_rec_g  <- if (!is.null(prune_recruit_max_dbh)) {
                if (isTRUE(prune_use_bio_recruit)) min(recruit_max_dbh, prune_recruit_max_dbh)
                else prune_recruit_max_dbh
            } else {
                recruit_max_dbh
            }

            cat(sprintf("[estimate_dp_complexity] Computing pruned edges for %d tags...\n",
                length(prune_candidates))); flush.console()

            for (tg in prune_candidates) {
                tg_K     <- out[Tag == tg, K]
                tg_data  <- data[Tag == tg & !is.na(DBH)]
                tg_meta  <- tag_meta[Tag == tg]
                eff_anch <- tg_meta$eff_anchor
                census_ids <- sort(unique(tg_data[CensusID <= eff_anch, CensusID]))
                n_c <- length(census_ids)
                if (n_c < 2L) next

                # Enumerate states per census (re-uses existing helper)
                all_mats <- vector("list", n_c)
                ok <- TRUE
                for (p in seq_len(n_c)) {
                    obs_in_census <- tg_data[CensusID == census_ids[p]]
                    n_obs_p <- nrow(obs_in_census)
                    all_mats[[p]] <- enumerate_assignments(tg_K, n_obs_p, max_states)
                    if (is.null(all_mats[[p]])) { ok <- FALSE; break }
                }
                if (!ok) next  # shouldn't happen (already non-fallback), but guard

                edges_p <- 0L
                for (p in seq_len(n_c - 1L)) {
                    dbh0 <- tg_data[CensusID == census_ids[p], DBH]
                    dbh1 <- tg_data[CensusID == census_ids[p + 1L], DBH]
                    iv   <- get_interval(census_ids[p], census_ids[p + 1L])
                    edges_p <- edges_p + count_pruned_edges_step_vec(
                        all_mats[[p]], all_mats[[p + 1L]],
                        dbh0, dbh1, tg_K, iv, eff_min_g, eff_max_g, eff_rec_g
                    )
                }
                out[Tag == tg, estimated_edges_pruned := as.numeric(edges_p)]
            }
            cat("[estimate_dp_complexity] Pruned edges done.\n"); flush.console()
        }

        # Tags in data but absent from obs_pre (all-NA DBH)
        tags_no_obs <- setdiff(unique(data$Tag), tag_meta$Tag)
        if (length(tags_no_obs) > 0L) {
            out_noobs <- data.table(
                Tag = tags_no_obs, Species = NA_character_,
                census_range = NA_character_, n_censuses = 0L, K = 0L,
                max_obs = 0L, max_states_per_census = 0, total_states = 0,
                estimated_edges_unpruned = 0, estimated_edges_pruned = NA_real_,
                estimated_fallback = FALSE, fallback_reason = "no_obs_up_to_anchor"
            )
            out <- rbindlist(list(out, out_noobs), use.names = TRUE, fill = TRUE)
        }

        # Runtime predictions
        out[, n_for_pred := ifelse(is.na(estimated_edges_pruned),
            estimated_edges_unpruned, estimated_edges_pruned)]
        out[is.na(n_for_pred), n_for_pred := 0]
        out[, predicted_seconds := predict_seconds(n_for_pred)]
        out[, predicted_hours   := predicted_seconds / 3600]
        out[, n_for_pred := NULL]

        # Sort: non-fallback first, then descending by edge count
        setorder(out, estimated_fallback, -estimated_edges_unpruned)

        keep_cols <- intersect(c(
            "Tag", "Species", "census_range", "n_censuses", "K", "max_obs",
            "max_states_per_census", "total_states",
            "estimated_edges_unpruned", "estimated_edges_pruned",
            "estimated_fallback", "fallback_reason",
            "predicted_seconds", "predicted_hours"
        ), names(out))
        cat("[estimate_dp_complexity] Done.\n"); flush.console()
        return(out[, ..keep_cols])
    }

    # ------------------------------------------------------------------
    # SLOW PATH (fast=FALSE): per-tag for loop for pruned-edge counts
    # ------------------------------------------------------------------
    setkey(data, Tag, CensusID)
    obs_summary_pre <- data[!is.na(DBH), .N, by = .(Tag, CensusID)]
    setkey(obs_summary_pre, Tag, CensusID)
    resprout_regex_pre <- "\\b(R|RP|RF|RT|QR)\\b"
    if ("ListOfTSM" %in% names(data)) {
        resprout_summary_pre <- data[
            !is.na(DBH) & !is.na(ListOfTSM) & grepl(resprout_regex_pre, ListOfTSM, perl = TRUE),
            .N, by = .(Tag, CensusID)
        ]
        setkey(resprout_summary_pre, Tag, CensusID)
    } else {
        resprout_summary_pre <- NULL
    }

    tags <- unique(data$Tag)
    n_tags <- length(tags)
    results <- vector("list", n_tags)
    report_at <- unique(round(seq(0, n_tags, length.out = 21L)))
    for (i in seq_along(tags)) {
        tg <- tags[i]
        td <- data[.(tg)]
        sp <- if (species_col %in% names(td)) unique(td[[species_col]])[1L] else NA_character_
        results[[i]] <- estimate_one_tag(td, tg, sp)
        if (i %in% report_at) {
            cat(sprintf("[estimate_dp_complexity] %d / %d tags done (%.0f%% remaining)\n",
                        i, n_tags, 100 * (1 - i / n_tags)))
            flush.console()
        }
    }
    out <- rbindlist(results, use.names = TRUE, fill = TRUE)

    out[, n_for_pred := ifelse(is.na(estimated_edges_pruned),
        estimated_edges_unpruned, estimated_edges_pruned)]
    out[, predicted_seconds := predict_seconds(n_for_pred)]
    out[, predicted_hours   := predicted_seconds / 3600]
    out[, n_for_pred := NULL]

    out[, .sort_key := ifelse(is.na(estimated_edges_pruned),
        estimated_edges_unpruned, estimated_edges_pruned)]
    setorder(out, estimated_fallback, -.sort_key)
    out[, .sort_key := NULL]
    out[]
}

#' Get detailed per-census breakdown for a specific tag
#'
#' Thin wrapper around estimate_dp_complexity that returns a list with
#' census-level detail for one tag.
#'
#' @param data RDS path, CSV path, or loaded data.table
#' @param tag Tag identifier
#' @param anchor_start Integer: anchor census (default 7)
#' @param slack_tracks Integer: slack tracks (default 1)
#' @param ... Additional arguments forwarded to estimate_dp_complexity
#' @return List with per-census detail
#' @export
get_tag_complexity_details <- function(data, tag,
                                       anchor_start = 7L,
                                       slack_tracks = 1L, ...) {
    library(data.table)

    if (is.character(data)) {
        if (grepl("\\.rds$", data, ignore.case = TRUE)) {
            data <- readRDS(data)
            setDT(data)
        } else {
            data <- fread(data)
        }
    } else {
        setDT(copy(data))
    }

    if (!(tag %in% data$Tag)) stop("Tag ", tag, " not found in data")

    result_row <- estimate_dp_complexity(
        data[Tag == tag],
        anchor_start = anchor_start,
        slack_tracks = slack_tracks, ...
    )

    # Reconstruct per-census detail for display
    tag_data <- data[Tag == tag]
    last_obs <- if (any(!is.na(tag_data$DBH))) max(tag_data$CensusID[!is.na(tag_data$DBH)]) else NA_integer_
    eff_anchor <- if (!is.na(last_obs) && last_obs < anchor_start) as.integer(last_obs) else as.integer(anchor_start)

    tmp_c <- tag_data[CensusID <= eff_anchor & !is.na(DBH), CensusID]
    if (length(tmp_c) == 0L) {
        return(list(summary = result_row, per_census = NULL))
    }

    first_obs <- as.integer(min(tmp_c))
    census_range <- seq.int(from = first_obs, to = eff_anchor)
    census_range <- census_range[vapply(census_range, function(cc) {
        nrow(tag_data[CensusID == cc & !is.na(DBH)]) > 0L || cc == eff_anchor
    }, logical(1L))]

    obs_counts <- vapply(
        census_range,
        function(cc) nrow(tag_data[CensusID == cc & !is.na(DBH)]), integer(1L)
    )

    K <- result_row$K[1L]
    count_inj <- function(K, n) {
        if (n == 0L) {
            return(1.0)
        }
        if (n > K) {
            return(0.0)
        }
        prod(seq.int(from = K, to = K - n + 1L, by = -1L))
    }
    n_states <- vapply(obs_counts, function(n) count_inj(K, n), numeric(1L))

    per_census <- data.table(
        CensusID  = census_range,
        N_Obs     = obs_counts,
        N_States  = n_states
    )

    list(
        summary    = result_row,
        per_census = per_census
    )
}

# Example usage:
if (FALSE) {
    library(data.table)
    source("dp_global/R/complexity/estimate_dp_complexity_function.R")
    d <- readRDS("bci_data/bci_multistem_xrun_debug.rds")
    comp <- estimate_dp_complexity(d, anchor_start = 7, slack_tracks = 1,
        max_states = 40000L, min_growth = -0.5, max_growth = 5,
        prune_use_bio_bounds = FALSE, recruit_max_dbh = (5 * 5) + 0.9999)
    print(head(comp, 20))
}


#' Sweep across parameter configurations and summarise predicted runtime
#'
#' Runs estimate_dp_complexity once per scenario and returns a summary table
#' showing how many tags use DP vs igraph, total predicted hours, and the
#' slowest-tag predicted time. Use this to decide how tightening pruning or
#' lowering max_states affects your batch runtime.
#'
#' @param data RDS path, CSV path, or already-loaded data.table/data.frame
#' @param scenarios A data.frame or data.table where each row is one parameter
#'   configuration.  Recognised column names (all optional — defaults used for
#'   omitted columns): \code{label}, \code{max_states}, \code{min_growth},
#'   \code{max_growth}, \code{prune_min_growth}, \code{prune_max_growth},
#'   \code{recruit_max_dbh}, \code{anchor_start}, \code{slack_tracks},
#'   \code{max_tracks}.
#' @param base_params Named list of defaults passed to estimate_dp_complexity()
#'   for any parameter not overridden by a scenario row.
#' @param return_full Logical: if TRUE, attach per-tag tables as an attribute
#'   `"full_results"` on the returned summary (default FALSE).
#' @return data.table with one row per scenario.
#' @export
sweep_dp_complexity <- function(data,
                                scenarios,
                                base_params = list(),
                                return_full = FALSE) {
    if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table required")
    library(data.table)

    # Load data once
    if (is.character(data)) {
        if (grepl("\\.rds$", data, ignore.case = TRUE)) {
            data <- readRDS(data)
            data.table::setDT(data)
        } else {
            data <- data.table::fread(data)
        }
    } else {
        data <- data.table::copy(data)
        data.table::setDT(data)
    }

    scenarios <- as.data.table(scenarios)
    if (!"label" %in% names(scenarios)) {
        scenarios[, label := paste0("scenario_", seq_len(.N))]
    }

    # Allowed overridable params with factory defaults
    defaults <- list(
        anchor_start       = 7L,
        slack_tracks       = 1L,
        max_states         = 40000L,
        max_tracks         = 9999L,
        min_growth         = -Inf,
        max_growth         = Inf,
        prune_min_growth   = NULL,
        prune_max_growth   = NULL,
        prune_use_bio_bounds = TRUE,
        recruit_max_dbh    = Inf,
        prune_recruit_max_dbh = NULL,
        prune_use_bio_recruit = TRUE,
        fast               = TRUE
    )
    # Merge user base_params over factory defaults
    for (nm in names(base_params)) defaults[[nm]] <- base_params[[nm]]

    results <- vector("list", nrow(scenarios))
    full_results <- if (isTRUE(return_full)) vector("list", nrow(scenarios)) else NULL

    for (s in seq_len(nrow(scenarios))) {
        row <- as.list(scenarios[s])
        lbl <- row$label

        # Build call args: start from defaults, then overlay scenario columns
        args <- defaults
        overridable <- setdiff(names(defaults), "fast")
        for (nm in overridable) {
            if (nm %in% names(row) && !is.na(row[[nm]])) {
                args[[nm]] <- row[[nm]]
            }
        }
        args$data <- data
        args$fast <- TRUE

        cat(sprintf("[sweep] Scenario %d/%d: %s  (max_states=%s, growth=[%.2f, %.2f])\n",
            s, nrow(scenarios), lbl,
            format(args$max_states, big.mark = ","),
            if (is.finite(args$min_growth)) args$min_growth else -Inf,
            if (is.finite(args$max_growth)) args$max_growth else Inf))
        flush.console()

        comp <- do.call(estimate_dp_complexity, args)

        n_dp     <- sum(!comp$estimated_fallback)
        n_igraph <- sum( comp$estimated_fallback)
        total_s  <- sum(comp$predicted_seconds, na.rm = TRUE)
        dp_only  <- comp[estimated_fallback == FALSE]
        slowest_tag    <- if (nrow(dp_only) > 0L) dp_only$Tag[1L] else NA
        slowest_sec    <- if (nrow(dp_only) > 0L) dp_only$predicted_seconds[1L] else NA_real_
        median_sec     <- if (nrow(dp_only) > 0L) median(dp_only$predicted_seconds, na.rm = TRUE) else NA_real_
        max_states_obs <- if (nrow(dp_only) > 0L) max(dp_only$max_states_per_census, na.rm = TRUE) else NA_real_

        results[[s]] <- data.table(
            label             = lbl,
            max_states        = as.numeric(args$max_states),
            min_growth        = args$min_growth,
            max_growth        = args$max_growth,
            prune_min         = if (!is.null(args$prune_min_growth)) args$prune_min_growth else NA_real_,
            prune_max         = if (!is.null(args$prune_max_growth)) args$prune_max_growth else NA_real_,
            recruit_max       = if (is.finite(args$recruit_max_dbh)) args$recruit_max_dbh else NA_real_,
            n_tags_total      = nrow(comp),
            n_tags_dp         = n_dp,
            n_tags_igraph     = n_igraph,
            pct_dp            = round(100 * n_dp / max(1L, nrow(comp)), 1),
            total_hours       = round(total_s / 3600, 2),
            slowest_tag       = slowest_tag,
            slowest_min       = round(slowest_sec / 60, 1),
            median_sec        = round(median_sec, 1),
            max_states_in_data = max_states_obs
        )
        if (isTRUE(return_full)) full_results[[s]] <- comp
    }
    out <- rbindlist(results, use.names = TRUE, fill = TRUE)
    if (isTRUE(return_full)) attr(out, "full_results") <- full_results
    out
}