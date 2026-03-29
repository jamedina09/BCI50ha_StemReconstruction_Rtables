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
#' @return data.table sorted descending by estimated_edges_pruned
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
                                   prune_use_bio_recruit = TRUE) {
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
        tmp_c <- tag_data[CensusID <= eff_anchor & !is.na(DBH), CensusID]
        if (length(tmp_c) == 0L) {
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
        first_obs <- as.integer(min(tmp_c))
        census_range <- seq.int(from = first_obs, to = eff_anchor)
        # Retain only censuses with at least one DBH obs (or the anchor endpoint)
        census_range <- census_range[vapply(
            census_range, function(cc) {
                nrow(tag_data[CensusID == cc & !is.na(DBH)]) > 0L || cc == eff_anchor
            },
            logical(1L)
        )]
        n_census <- length(census_range)

        obs_counts <- vapply(
            census_range,
            function(cc) nrow(tag_data[CensusID == cc & !is.na(DBH)]), integer(1L)
        )
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

        # Resprout barrier (same regex as DP)
        resprout_regex <- "\\b(R|RP|RF|RT|QR)\\b"
        n_resprout <- 0L
        if ("ListOfTSM" %in% names(tag_data)) {
            for (cc in census_range) {
                obs_tmp <- tag_data[CensusID == cc & !is.na(DBH)]
                if (nrow(obs_tmp) > 0L) {
                    n_resprout <- n_resprout + sum(
                        !is.na(obs_tmp$ListOfTSM) & grepl(resprout_regex, obs_tmp$ListOfTSM)
                    )
                }
            }
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
        #     Skipped (NA) when tag would fall back to igraph
        edges_pruned <- NA_real_
        if (!fallback && n_census >= 2L) {
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
    # 6. Run over all tags and sort
    # ------------------------------------------------------------------
    tags <- unique(data$Tag)
    results <- vector("list", length(tags))
    for (i in seq_along(tags)) {
        tg <- tags[i]
        td <- data[Tag == tg]
        sp <- if (species_col %in% names(td)) unique(td[[species_col]])[1L] else NA_character_
        results[[i]] <- estimate_one_tag(td, tg, sp)
    }
    out <- rbindlist(results, use.names = TRUE, fill = TRUE)

    # Sort: non-fallback tags first; within each group descending by pruned edges
    out[, .sort_key := ifelse(is.na(estimated_edges_pruned),
        estimated_edges_unpruned,
        estimated_edges_pruned
    )]
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
    source("dp_global/R/estimate_dp_complexity_function.R")
    d <- readRDS("bci_data/bci_multistem_xrun_debug.rds")

    # Rank all tags by estimated cost
    comp <- estimate_dp_complexity(d,
        anchor_start = 7, slack_tracks = 1,
        min_growth = -0.5, max_growth = 7.5,
        prune_use_bio_bounds = FALSE,
        recruit_max_dbh = 100
    )
    print(head(comp, 20))

    # Detailed breakdown for one tag
    det <- get_tag_complexity_details(d,
        tag = 156669, anchor_start = 7, slack_tracks = 1,
        min_growth = -0.5, max_growth = 5,
        prune_use_bio_bounds = FALSE, recruit_max_dbh = (5 * 5) + 0.9999
    )
    print(det)
}