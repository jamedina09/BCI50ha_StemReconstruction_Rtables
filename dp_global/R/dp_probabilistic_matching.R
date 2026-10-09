############################################################
# dp_probabilistic_matching.R
# Probabilistic matcher: stochastic assignment of the stems of each census
# pair, for the tags the DP does not solve
############################################################
# When the DP cannot be used (state space too large, a dead end with no
# feasible state, species or growth-form routing; see do_fallback() in
# dp_global_dp.R):
#   1. Pairwise log-likelihoods from the Bio_* parameters the DP uses
#      (compute_pairwise_log_likelihood(): Gaussian growth, survival and the
#      soft penalties; the DP's measurement-error mixture is not implemented
#      here)
#   2. Augment cost matrix with mortality/recruitment slots. Birth-death mode
#      (birth_death = TRUE, the default): every stem may continue, die or be
#      recruited in every census pair, one death cell per current stem and
#      one recruit cell per next stem, so a link is taken only when it is
#      more likely than the death of the earlier stem plus the recruitment of
#      the later one (as in the DP). Count-based slots (birth_death = FALSE):
#      only as many slots as the stem counts and the forbidden links require,
#      so with stable counts every stem is forced to continue
#   3. Draw n_samples stochastic assignments. Birth-death mode: Gumbel noise
#      on every event cell and the exact best assignment of each perturbed
#      pair (hungarian_min_rcpp). Count-based slots: Gumbel-noise greedy; a
#      pair whose greedy assignment uses a forbidden (-Inf) link is re-solved
#      exactly (enforce_feasible_assignment), so the hard limits hold; a
#      violation no assignment can avoid is chosen as the DP does (hard
#      penalty + likelihood without the hard gate). TrueStemID pins of every
#      stem are honoured while sampling: a pinned obs joins the track that
#      carries its pin and never one that carries another pin
#      (apply_pin_mask for anchor stems, apply_track_pin_mask for all); when
#      the masks forbid links a pair's slots relied on, the pair gets extra
#      death/recruit slots (pin_masked_pair), so two stems with different
#      pins are never forced into one track
#   4. Stitch per-pair assignments backward from anchor
#   5. Repair growth violations at SAMPLE level (hard-rate; the ME
#      cumulative-shrinkage layer is off when n_sigma_me = Inf); pinned
#      observations are never severed
#   6. Drop samples that break a pin to an anchor stem
#      (filter_pin_consistent_samples; all are kept when too few remain)
#   7. Compute marginal posterior probabilities (Top-K, entropy)
#   8. Export ONE coherent trajectory: the consensus (maximum-expected-
#      accuracy) sample (select_consensus_trajectory)
#   9. Pin sweep: every row with a TrueStemID (not a provisional anchor)
#      gets ReconstructedStemID = TrueStemID; the engine's value is kept in
#      ReconstructedStemID_PreSweep
#
# All Bio_* parameters are read directly from tree_data columns (no new
# estimation needed — they are already computed by dp_global_bio.R).

# ---- Main entry point ----------------------------------------------------
# match_stems_probabilistic(): reconstruct the stem identities of one tag.
#
# INPUTS
#   tree_data        data.table of one tag with CensusID, DBH, TrueStemID,
#                    ExactDate and the Bio_* columns (Bio_Mu_Growth,
#                    Bio_Sigma0_Growth, Bio_Sigma1_Growth, Bio_Max_Shrink,
#                    Bio_K_Shrink, Bio_H0_Mortality, Bio_Beta_Mortality,
#                    Bio_Recruit_Meanlog, Bio_Recruit_Sdlog,
#                    Bio_Recruit_MaxDBH_unit, Bio_Recruitment_lambda; optional
#                    Bio_Gamma_Growth, Bio_Max_Growth, Bio_K_Growth). The
#                    first value of each Bio_* column is used.
#   min_growth, max_growth   hard bounds on annual growth (cm/yr), used when
#                    the prune_* bounds are NULL
#   anchor_start     CensusID of the anchor census. When no measured row at
#                    that census has a TrueStemID, the last census with a DBH
#                    is the anchor. Rows after the anchor are not matched.
#   n_samples        number of stochastic assignments to draw
#   temperature      scale of the Gumbel noise; higher = more random
#   posterior_top_k  number of DP_PosteriorTop{k}ID / Prob columns (at least 1)
#   posterior_samples_path, posterior_samples_format
#                    passed to export_probabilistic_posteriors()
#   posterior_sample_seed    seed set before sampling (NULL: none)
#   prune_min_growth, prune_max_growth
#                    replace min_growth / max_growth when not NULL
#   prune_recruit_max_dbh    replaces Bio_Recruit_MaxDBH_unit when not NULL
#   The remaining arguments are described next to their defaults below.
#
# RETURNS  tree_data (a copy ordered by CensusID) with ReconstructedStemID on
#   its measured rows up to the anchor, ReconstructionMethod ("probabilistic",
#   "given", "provisional_dp"), the posterior columns DP_PosteriorTop{k}ID /
#   Prob, DP_PosteriorEntropy and DP_PosteriorReconstructedProb, and, when
#   pin_truestemid is TRUE and a row has a TrueStemID that is not a
#   provisional anchor ID, SweepAuditOverride and
#   ReconstructedStemID_PreSweep. ConstraintViolation and obs_row_id are added
#   when missing. The posterior samples are staged on disk, or attached as
#   attribute "DP_Posterior_Samples" when return_samples = TRUE. A tag with no
#   DBH, or with one measured census, returns before any sampling.

match_stems_probabilistic <- function(tree_data,
                                      min_growth,
                                      max_growth,
                                      anchor_start,
                                      n_samples = 200L,
                                      temperature = 1.0,
                                      posterior_top_k = 2L,
                                      posterior_samples_path = NULL,
                                      posterior_samples_format = "csv",
                                      posterior_sample_seed = NULL,
                                      prune_min_growth = NULL,
                                      prune_max_growth = NULL,
                                      prune_recruit_max_dbh = NULL,
                                      prob_lookahead_weight = 0.5, # backward conditioning weight [0,1]; 0 = independent
                                      use_bio_hard_shrink_in_prob = TRUE, # if FALSE, ignore Bio_Max_Shrink hard gate
                                      use_bio_hard_growth_in_prob = TRUE, # if FALSE, ignore Bio_Max_Growth hard gate
                                      pin_truestemid = TRUE, # pin obs with known TrueStemID to their track
                                      n_sigma_me = 3, # ME cumulative-shrinkage threshold (n * SD); lower = sever sooner
                                      dbh_round_censuses = integer(0), # CensusIDs whose small-stem DBH was rounded down (classes)
                                      dbh_round_max = 5.5, # only DBH below this (cm) was rounded
                                      dbh_round_width = 0.5, # class width (cm)
                                      return_samples = FALSE, # TRUE: attach samples as attr "DP_Posterior_Samples" instead of staging them
                                      birth_death = TRUE, # TRUE: every stem may die / be recruited in every pair (see header); FALSE: count-based slots
                                      verbose = FALSE) {
    tree_data <- tree_data[order(CensusID)]
    n_samples <- as.integer(n_samples)
    birth_death <- isTRUE(birth_death)
    posterior_top_k <- max(1L, as.integer(posterior_top_k))

    vcat <- function(...) {
        if (!isTRUE(verbose)) {
            return(invisible(NULL))
        }
        cat(..., "\n")
        flush.console()
        invisible(NULL)
    }

    tag_val <- tryCatch(
        {
            u <- unique(tree_data$Tag)
            u <- u[!is.na(u)]
            if (length(u) == 1L) u[[1L]] else NA
        },
        error = function(e) NA
    )
    prefix <- paste0("[prob_match Tag=", if (!is.na(tag_val)) tag_val else "?", "] ")

    # --- Resolve effective prune bounds (mirrors DP logic) ----------------
    eff_min_growth <- if (!is.null(prune_min_growth)) prune_min_growth else min_growth
    eff_max_growth <- if (!is.null(prune_max_growth)) prune_max_growth else max_growth

    # --- Extract bio parameters from tree_data (same as dp_global_dp.R) ---
    bio <- list(
        mu_const = unique(tree_data$Bio_Mu_Growth)[1],
        mu_gamma = if ("Bio_Gamma_Growth" %in% names(tree_data)) unique(tree_data$Bio_Gamma_Growth)[1] else 0,
        sigma0 = unique(tree_data$Bio_Sigma0_Growth)[1],
        sigma1 = unique(tree_data$Bio_Sigma1_Growth)[1],
        max_shrink = unique(tree_data$Bio_Max_Shrink)[1],
        k_shrink = unique(tree_data$Bio_K_Shrink)[1],
        max_growth_bio = if ("Bio_Max_Growth" %in% names(tree_data)) unique(tree_data$Bio_Max_Growth)[1] else Inf,
        k_growth = if ("Bio_K_Growth" %in% names(tree_data)) unique(tree_data$Bio_K_Growth)[1] else 0,
        h0 = unique(tree_data$Bio_H0_Mortality)[1],
        beta_mort = unique(tree_data$Bio_Beta_Mortality)[1],
        recruit_meanlog = unique(tree_data$Bio_Recruit_Meanlog)[1],
        recruit_sdlog = unique(tree_data$Bio_Recruit_Sdlog)[1],
        recruit_max_dbh = if (!is.null(prune_recruit_max_dbh)) prune_recruit_max_dbh else unique(tree_data$Bio_Recruit_MaxDBH_unit)[1],
        recruit_lambda = unique(tree_data$Bio_Recruitment_lambda)[1]
    )

    # --- Identify censuses and observations --------------------------------
    if (!("ReconstructionMethod" %in% names(tree_data))) {
        tree_data[, ReconstructionMethod := NA_character_]
    }
    if (!("ConstraintViolation" %in% names(tree_data))) {
        tree_data[, ConstraintViolation := NA]
    }
    if (!("obs_row_id" %in% names(tree_data))) {
        tree_data[, obs_row_id := seq_len(.N)]
    }

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
        # Single census — honor TrueStemID where present, otherwise sequential.
        # The downstream sweep would otherwise have to override seq_len(.N)
        # values on every anchored row, polluting ReconstructedStemID_PreSweep
        # and inflating SweepAuditOverride.
        tree_data[
            CensusID == anchor_start & !is.na(DBH),
            ReconstructedStemID := fifelse(
                !is.na(TrueStemID),
                as.integer(TrueStemID),
                seq_len(.N)
            )
        ]
        tree_data[is.na(ReconstructedStemID), ReconstructedStemID := NA_integer_]
        tree_data[, ReconstructionMethod := "probabilistic"]
        tree_data[
            CensusID == anchor_start & !is.na(DBH) & !is.na(TrueStemID),
            ReconstructionMethod := "given"
        ]
        return(tree_data)
    }

    # Build per-census observation data
    obs_data <- vector("list", n_census)
    for (i in seq_len(n_census)) {
        cc <- obs_census[i]
        idx <- which(tree_data$CensusID == cc & !is.na(tree_data$DBH))
        obs_data[[i]] <- list(
            census_id = cc,
            idx       = idx, # row indices in tree_data
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
        if (!is.finite(intervals[i]) || intervals[i] <= 0) intervals[i] <- 5.0 # safe fallback
    }

    # --- Anchor IDs --------------------------------------------------------
    # Build IDs for ALL anchor observations (one per row).  Where TrueStemID
    # is available, use it; where NA, assign new sequential IDs starting
    # above the max known TrueStemID so they don't collide.
    anchor_pos <- n_census # anchor is the last observed census
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
    .has_tsid <- !is.na(tree_data$TrueStemID[.anchor_rows])
    .method_vec <- ifelse(.is_provisional, "provisional_dp",
        ifelse(.has_tsid, "given", "probabilistic")
    )
    tree_data[.anchor_rows, ReconstructionMethod := .method_vec]

    # K = number of tracks (at least max obs across any census)
    max_obs <- max(vapply(obs_data, function(x) x$n, integer(1)))
    K <- max(length(anchor_ids), max_obs)

    vcat(
        prefix, "Probabilistic matching: ", n_census, " censuses, K=", K,
        ", max_obs=", max_obs, ", n_samples=", n_samples
    )

    # --- Pre-compute TrueStemID pin map for non-anchor censuses -------------
    # pin_info[[i]][j] = anchor-position index (1..n_anchor) for obs j, or NA
    pin_info <- vector("list", n_census)
    .any_pins <- FALSE
    if (isTRUE(pin_truestemid)) {
        n_anchor <- length(anchor_ids)
        for (i in seq_len(n_census)) {
            if (i == n_census) next # anchor pinned via anchor_ids directly
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
                    vcat(
                        prefix, "WARNING: duplicate TrueStemID pin at C",
                        obs_data[[i]]$census_id, " for anchor ID ", anchor_ids[tidx[.j]],
                        "; keeping first, releasing obs ", .j
                    )
                    tidx[.j] <- NA_integer_
                } else {
                    .seen <- c(.seen, tidx[.j])
                }
            }
            if (any(!is.na(tidx))) .any_pins <- TRUE
            pin_info[[i]] <- tidx
        }
    }

    # --- Pins of every observation, for the track-pin mask ------------------
    # pin_info above covers only pins that point to a stem present at the
    # anchor. pin_val also covers the pins on stems that end before the anchor
    # (Step 3a.5 / Step 3b StemID pins, trees last measured before 2010), so
    # that the samples group those observations as the pin sweep groups them
    # in the exported table. pin_val[[i]][j] is the
    # TrueStemID of obs j at census i (NA when not pinned; provisional anchor
    # IDs are not pins). While sampling backward, every next-census obs carries
    # the pin of its track (its own pin, or one inherited from a later census
    # of the same sample): a pinned obs must join the obs that carries its pin
    # and may not join an obs that carries another pin (apply_track_pin_mask).
    pin_val <- vector("list", n_census)
    .any_track_pins <- FALSE
    if (isTRUE(pin_truestemid)) {
        .is_prov_row <- tree_data$ReconstructionMethod %in% "provisional_dp"
        for (i in seq_len(n_census)) {
            .idx <- obs_data[[i]]$idx
            v <- as.integer(tree_data$TrueStemID[.idx])
            v[.is_prov_row[.idx]] <- NA_integer_
            v[!is.na(v) & duplicated(v)] <- NA_integer_ # duplicate-pin guard: keep first
            pin_val[[i]] <- v
            if (i < n_census && any(!is.na(v))) .any_track_pins <- TRUE
        }
    }

    # --- Per-pair log-likelihood matrices (backward) -----------------------
    # For each pair (c, c+1) compute the pairwise + augmented cost matrix
    pair_data <- vector("list", n_census - 1L)
    for (i in seq_len(n_census - 1L)) {
        dbh_curr <- obs_data[[i]]$dbh
        dbh_next <- obs_data[[i + 1L]]$dbh
        iv <- intervals[i]

        .round_curr <- obs_census[i] %in% dbh_round_censuses
        .round_next <- obs_census[i + 1L] %in% dbh_round_censuses
        L <- compute_pairwise_log_likelihood(dbh_curr, dbh_next, iv, bio,
            eff_min_growth, eff_max_growth,
            use_bio_hard_shrink = use_bio_hard_shrink_in_prob,
            use_bio_hard_growth = use_bio_hard_growth_in_prob,
            round_curr = .round_curr, round_next = .round_next,
            round_max = dbh_round_max, round_width = dbh_round_width
        )
        aug <- augment_cost_matrix(L, dbh_curr, dbh_next, iv, bio, birth_death = birth_death)
        .k_raised <- attr(aug, "k_raised")
        if (!is.null(.k_raised)) {
            vcat(prefix, sprintf(
                "Census pair C%d-C%d: %d death/recruit slot(s) added (K %d -> %d) so that no forbidden link is forced",
                obs_census[i], obs_census[i + 1L], .k_raised[2L] - .k_raised[1L], .k_raised[1L], .k_raised[2L]
            ))
        }
        # Same matrix without the hard gates: used only when no assignment
        # avoids every forbidden link (see enforce_feasible_assignment())
        L_free <- compute_pairwise_log_likelihood(dbh_curr, dbh_next, iv, bio,
            -Inf, Inf,
            use_bio_hard_shrink = FALSE,
            use_bio_hard_growth = FALSE,
            round_curr = .round_curr, round_next = .round_next,
            round_max = dbh_round_max, round_width = dbh_round_width
        )

        pair_data[[i]] <- list(
            log_cost = aug,
            fallback = fallback_log_cost(aug, L_free, dbh_next, iv, bio),
            n_curr   = length(dbh_curr),
            n_next   = length(dbh_next),
            # kept for pin_masked_pair(), which may rebuild the pair with more slots
            dbh_curr = dbh_curr,
            dbh_next = dbh_next,
            iv       = iv,
            L_free   = L_free
        )
    }

    # --- Draw n_samples stochastic assignments (backward from anchor) ------
    if (!is.null(posterior_sample_seed)) set.seed(as.integer(posterior_sample_seed))

    all_samples <- vector("list", n_samples)
    use_lookahead <- is.finite(prob_lookahead_weight) && prob_lookahead_weight > 0 &&
        n_census >= 3L && K >= 4L
    .n_pin_slots <- 0L # pair draws that needed extra slots for the pins (pin_masked_pair)

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
        # Pin carried by the track of each next-census obs (anchor: its own pin)
        .next_track_pin <- pin_val[[n_census]]

        # Last pair (closest to anchor): no conditioning available
        last_pair <- n_census - 1L
        .cost_last <- pair_data[[last_pair]]$log_cost
        .fb_last <- pair_data[[last_pair]]$fallback
        if (.any_pins || .any_track_pins) {
            .mp <- pin_masked_pair(.cost_last, .fb_last, pair_data[[last_pair]], bio, function(m) {
                if (.any_pins && !is.null(pin_info[[last_pair]])) {
                    m <- apply_pin_mask(
                        m, pin_info[[last_pair]], .next_obs_to_anchor_pos,
                        pair_data[[last_pair]]$n_curr, pair_data[[last_pair]]$n_next
                    )
                }
                if (.any_track_pins) {
                    m <- apply_track_pin_mask(
                        m, pin_val[[last_pair]], .next_track_pin,
                        pair_data[[last_pair]]$n_curr, pair_data[[last_pair]]$n_next
                    )
                }
                m
            })
            .cost_last <- .mp$cost
            .fb_last <- .mp$fb
            .n_pin_slots <- .n_pin_slots + .mp$grown
        }
        per_pair_assignments[[last_pair]] <- greedy_assignment_gumbel(
            .cost_last,
            temperature = temperature,
            fallback = .fb_last
        )
        if (.any_pins) {
            .next_obs_to_anchor_pos <- propagate_track_backward(
                per_pair_assignments[[last_pair]], .next_obs_to_anchor_pos,
                pair_data[[last_pair]]$n_curr, pair_data[[last_pair]]$n_next
            )
        }
        if (.any_track_pins) {
            .next_track_pin <- propagate_track_pin(
                per_pair_assignments[[last_pair]], pin_val[[last_pair]], .next_track_pin,
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

                fb_i <- pair_data[[i]]$fallback
                if (.any_pins || .any_track_pins) {
                    .mp <- pin_masked_pair(cost_i, fb_i, pair_data[[i]], bio, function(m) {
                        if (.any_pins && !is.null(pin_info[[i]])) {
                            m <- apply_pin_mask(
                                m, pin_info[[i]], .next_obs_to_anchor_pos,
                                pair_data[[i]]$n_curr, pair_data[[i]]$n_next
                            )
                        }
                        if (.any_track_pins) {
                            m <- apply_track_pin_mask(
                                m, pin_val[[i]], .next_track_pin,
                                pair_data[[i]]$n_curr, pair_data[[i]]$n_next
                            )
                        }
                        m
                    })
                    cost_i <- .mp$cost
                    fb_i <- .mp$fb
                    .n_pin_slots <- .n_pin_slots + .mp$grown
                }

                per_pair_assignments[[i]] <- greedy_assignment_gumbel(
                    cost_i,
                    temperature = temperature,
                    fallback = fb_i
                )
                if (.any_pins) {
                    .next_obs_to_anchor_pos <- propagate_track_backward(
                        per_pair_assignments[[i]], .next_obs_to_anchor_pos,
                        pair_data[[i]]$n_curr, pair_data[[i]]$n_next
                    )
                }
                if (.any_track_pins) {
                    .next_track_pin <- propagate_track_pin(
                        per_pair_assignments[[i]], pin_val[[i]], .next_track_pin,
                        pair_data[[i]]$n_curr, pair_data[[i]]$n_next
                    )
                }
            }
        }

        all_samples[[s]] <- per_pair_assignments
    }
    .n_unavoidable <- sum(vapply(all_samples, function(sa) {
        sum(vapply(sa, function(a) {
            nf <- attr(a, "n_forbidden")
            if (is.null(nf)) 0L else as.integer(nf)
        }, integer(1)))
    }, integer(1)))
    if (.n_unavoidable > 0L) {
        vcat(prefix, sprintf(
            "Unavoidable forbidden link(s) in the sampled assignments: %d across %d samples (chosen by likelihood, as the DP's hard_penalty)",
            .n_unavoidable, n_samples
        ))
    }
    if (.n_pin_slots > 0L) {
        vcat(prefix, sprintf(
            "Pins: %d census-pair draw(s) across %d samples got extra death/recruit slots so that no two stems with different pins were joined",
            .n_pin_slots, n_samples
        ))
    }

    # --- Stitch assignments backward from anchor ---------------------------
    stitched <- stitch_assignments_backward(all_samples, obs_data, anchor_ids, K)

    # --- Repair growth violations at the SAMPLE level (before marginals) ---
    # This ensures probabilities only count biologically valid paths.
    # Two layers: hard-rate bounds + ME-informed cumulative shrinkage.
    stitched <- repair_stitched_growth_violations(
        stitched, obs_data, intervals, eff_min_growth, eff_max_growth,
        me_sd1_a = 0.0062, me_sd1_b = 0.0904, n_sigma_me = n_sigma_me,
        use_bio_hard_shrink = use_bio_hard_shrink_in_prob,
        pinned = if (isTRUE(pin_truestemid)) lapply(pin_val, function(v) !is.na(v)) else NULL
    )
    .sample_breaks <- attr(stitched, "sample_level_breaks")
    .sample_me_breaks <- attr(stitched, "sample_level_me_breaks")
    if (!is.null(.sample_breaks) && .sample_breaks > 0L) {
        .msg <- paste0(
            prefix, "Sample-level repair: ", .sample_breaks,
            " growth violation(s) broken across ", n_samples, " samples",
            if (!is.null(.sample_me_breaks) && .sample_me_breaks > 0L) {
                paste0(" (", .sample_me_breaks, " from ME cumulative-shrinkage check)")
            } else {
                ""
            }
        )
        vcat(.msg)
        message(.msg) # ensure it appears on stderr / captured by log redirection
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
        max_rate = eff_max_growth
    )

    # --- Export one coherent trajectory ------------------------------------
    # ReconstructedStemID = the consensus (maximum-expected-accuracy) sample,
    # not the per-census marginal resolution above. Top-K and entropy columns
    # stay marginal; DP_PosteriorReconstructedProb = share of samples giving
    # the observation the exported ID.
    .cons <- select_consensus_trajectory(stitched, obs_data)
    .chosen <- stitched[[.cons$index]]
    for (ci in seq_len(n_census)) {
        .n <- obs_data[[ci]]$n
        if (.n == 0L) next
        .ids <- as.integer(.chosen[[ci]][seq_len(.n)])
        .p <- vapply(seq_len(.n), function(oi) {
            mean(vapply(stitched, function(s) isTRUE(as.integer(s[[ci]][oi]) == .ids[oi]), logical(1)))
        }, numeric(1))
        data.table::set(tree_data, as.integer(obs_data[[ci]]$idx), "ReconstructedStemID", .ids)
        data.table::set(tree_data, as.integer(obs_data[[ci]]$idx), "DP_PosteriorReconstructedProb", .p)
    }
    vcat(prefix, sprintf(
        "Consensus trajectory exported: sample %d of %d | predecessor agreement %.3f | partition share %.3f",
        .cons$index, length(stitched), .cons$agreement, .cons$part_freq
    ))

    # --- Diagnostic check: count growth violations left in the exported
    #     (consensus) trajectory.  The sample-level repair removes them,
    #     except on links it keeps (e.g. between two pinned observations).
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
            .msg <- paste0(
                prefix, "WARNING: ", .n_violations,
                " residual growth violation(s) in marginal trajectory ",
                "(post-marginal repair DISABLED — posteriors preserved)"
            )
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
    .has_tsid <- !is.na(tree_data$TrueStemID)
    .is_prov <- tree_data$ReconstructionMethod %in% "provisional_dp"
    .is_pinned <- .has_tsid & !.is_prov # real TrueStemID = pinned
    tree_data[, ReconstructionMethod := "probabilistic"]
    tree_data[.is_pinned & (.is_anchor | isTRUE(pin_truestemid)), ReconstructionMethod := "given"]
    tree_data[.is_prov, ReconstructionMethod := "provisional_dp"]

    # ---- Hard-invariant final sweep ----------------------------------------
    # The sampling above honours the pins as groupings: pin_info[[i]] for rows
    # whose TrueStemID is present at the anchor census (i.e. is in anchor_ids),
    # pin_val[[i]] for every pinned row.  A track that ends before the anchor
    # is labelled with a matcher ID (stitch_assignments_backward()), not with
    # its pin, and rows without DBH are not sampled at all.  This sweep
    # enforces the hard invariant
    # ReconstructedStemID == TrueStemID for ANY row with non-NA TrueStemID,
    # regardless of DBH or census position.  Mirrors the equivalent sweep
    # in dp_global_dp.R::finalize_out().
    if (isTRUE(pin_truestemid)) {
        .has_true_final <- !is.na(tree_data$TrueStemID) &
            !(tree_data$ReconstructionMethod %in% "provisional_dp")
        if (any(.has_true_final)) {
            # ---- Audit: detect engine-vs-pin disagreements -----------------
            # Mirrors the audit in dp_global_dp.R finalize_out.  Flags any
            # row where the probabilistic engine had already assigned a
            # non-NA ReconstructedStemID different from TrueStemID.
            if (!("SweepAuditOverride" %in% names(tree_data))) {
                tree_data[, SweepAuditOverride := FALSE]
            } else {
                # Per-row backfill — mirrors dp_global_dp.R finalize_out.
                tree_data[is.na(SweepAuditOverride), SweepAuditOverride := FALSE]
            }
            # Snapshot the engine's pre-sweep ReconstructedStemID.
            # Per-row backfill (see dp_global_dp.R::finalize_out for details):
            # create the column if absent, otherwise populate only NA cells.
            if (!("ReconstructedStemID_PreSweep" %in% names(tree_data))) {
                tree_data[, ReconstructedStemID_PreSweep := ReconstructedStemID]
            } else {
                tree_data[
                    is.na(ReconstructedStemID_PreSweep),
                    ReconstructedStemID_PreSweep := ReconstructedStemID
                ]
            }
            # Audit flag is computed against the PreSweep snapshot
            # (mirror dp_global_dp.R::finalize_out) so post-engine
            # renumbering doesn't falsely flag rows whose original
            # output already agreed with TrueStemID.
            .pre_snap <- tree_data$ReconstructedStemID_PreSweep
            .override <- .has_true_final &
                !is.na(.pre_snap) &
                .pre_snap != as.integer(tree_data$TrueStemID)
            if (any(.override)) {
                tree_data[.override, SweepAuditOverride := TRUE]
                .tag_id <- if ("Tag" %in% names(tree_data) && length(tree_data$Tag) > 0L) as.character(tree_data$Tag[[1L]]) else "<unknown>"
                if (isTRUE(verbose)) message(sprintf("[probabilistic] [audit] sweep overrode %d engine-assigned ReconstructedStemID value(s) for tag %s (rows flagged via SweepAuditOverride=TRUE)", sum(.override), .tag_id))
            }
            tree_data[.has_true_final, `:=`(
                ReconstructedStemID = as.integer(TrueStemID),
                ReconstructionMethod = "given"
            )]
        }
    }

    # --- Export posterior samples (same format as DP) -----------------------
    if (n_samples > 0L) {
        .samples_dt <- export_probabilistic_posteriors(
            stitched, tree_data, obs_data, obs_census,
            tag_val = tag_val,
            n_samples = length(stitched),
            posterior_samples_path = posterior_samples_path,
            posterior_samples_format = posterior_samples_format,
            verbose = verbose,
            prefix = prefix,
            vcat = vcat,
            stage = !isTRUE(return_samples)
        )
        if (isTRUE(return_samples)) attr(tree_data, "DP_Posterior_Samples") <- .samples_dt
    }

    tree_data
}

# ---- Pairwise log-likelihood matrix --------------------------------------
# Log-likelihood of every link between a stem at census c and a stem at
# census c+1: Gaussian growth likelihood of the annual growth, plus the log
# survival probability of the earlier stem, plus the soft penalties.
#
# INPUTS
#   dbh_curr, dbh_next   DBH (cm) of the stems at census c and c+1
#   interval_years       years between the two censuses
#   bio                  list built in match_stems_probabilistic(): mu_const,
#                        mu_gamma, sigma0, sigma1 (growth mean and SD), h0,
#                        beta_mort (mortality hazard), max_shrink,
#                        max_growth_bio (bio hard bounds), k_shrink, k_growth
#                        (soft penalty weights)
#   min_growth, max_growth   hard bounds on annual growth (cm/yr); a link
#                        outside them is forbidden (not applied when
#                        non-finite)
#   use_bio_hard_shrink, use_bio_hard_growth
#                        also forbid links below bio$max_shrink / above
#                        bio$max_growth_bio
#   round_curr, round_next   TRUE when census c / c+1 recorded DBH in classes,
#                        rounded down
#   round_max, round_width   only DBH below round_max (cm) was rounded, in
#                        classes of round_width (cm)
#
# RETURNS  n_curr × n_next matrix of log-likelihoods; -Inf for a forbidden
#          link or a non-finite DBH.

compute_pairwise_log_likelihood <- function(dbh_curr, dbh_next, interval_years,
                                            bio, min_growth, max_growth,
                                            use_bio_hard_shrink = TRUE,
                                            use_bio_hard_growth = TRUE,
                                            round_curr = FALSE,
                                            round_next = FALSE,
                                            round_max = 5.5,
                                            round_width = 0.5) {
    n_curr <- length(dbh_curr)
    n_next <- length(dbh_next)
    L <- matrix(-Inf, nrow = n_curr, ncol = n_next)

    mu_growth_fn <- function(d) {
        if (!is.finite(bio$mu_gamma) || bio$mu_gamma == 0 || !is.finite(d) || d <= 0) {
            return(bio$mu_const)
        }
        bio$mu_const + bio$mu_gamma * log(d)
    }

    for (i in seq_len(n_curr)) {
        d0 <- dbh_curr[i]
        if (!is.finite(d0)) next
        for (j in seq_len(n_next)) {
            d1 <- dbh_next[j]
            if (!is.finite(d1)) next

            g <- (d1 - d0) / interval_years
            # DBH rounded down to classes of round_width at flagged censuses
            # (only below round_max): the true DBH lies in [d, d + round_width),
            # as in transition_cost_rcpp.cpp: the likelihood uses the class
            # mid-points and adds round_width^2 / 12 per rounded measurement to
            # the variance. The hard gates stay on the measured DBHs.
            r0 <- isTRUE(round_curr) && d0 < round_max
            r1 <- isTRUE(round_next) && d1 < round_max
            g_mid <- g + ((if (r1) 0.5 * round_width else 0) - (if (r0) 0.5 * round_width else 0)) / interval_years
            var_round <- (r0 + r1) * round_width^2 / 12 / interval_years^2
            d0_mid <- d0 + (if (r0) 0.5 * round_width else 0)
            # Hard growth constraints — infeasible edge
            if (is.finite(min_growth) && g < min_growth) next
            if (is.finite(max_growth) && g > max_growth) next
            # Bio hard shrink gate (conditionally applied)
            if (isTRUE(use_bio_hard_shrink) && is.finite(bio$max_shrink) && g < bio$max_shrink) next
            # Bio hard growth gate (conditionally applied)
            if (isTRUE(use_bio_hard_growth) && is.finite(bio$max_growth_bio) && g > bio$max_growth_bio) next

            # Growth likelihood (Gaussian; process SD plus rounding, if any)
            sigma_d <- max(bio$sigma0 + bio$sigma1 * d0_mid, 1e-6)
            mu <- mu_growth_fn(d0_mid)
            sd_g <- if (var_round > 0) sqrt(sigma_d^2 + var_round) else sigma_d
            ll_growth <- dnorm(g_mid, mean = mu, sd = sd_g, log = TRUE)

            # Survival probability
            hazard <- bio$h0 * exp(bio$beta_mort * d0)
            p_surv <- exp(-hazard * interval_years)
            p_surv <- max(1e-12, min(1 - 1e-12, p_surv))
            ll_surv <- log(p_surv)

            # Soft penalties: shrinkage as in the C++ cost; excess growth is
            # measured above the bio hard bound max_growth_bio (the C++ cost
            # uses Bio_Max_Growth_Soft), so it only acts when the bio growth
            # gate is off
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
# K_min: optional smallest K, used by pin_masked_pair() when a sample's pin
# masks forbid survival links that the K sized here relied on.
#
# Birth-death mode (birth_death = TRUE): K = n_curr + n_next and every stem
# has its own event cells, so any combination of survivals, deaths and
# recruitments is an assignment:
#   [1:n_curr, 1:n_next]          survival link: growth + survival log-likelihood (L)
#   [i, n_next + i]               death of current stem i: log P(death)
#   [n_curr + j, j]               recruitment of next stem j: log P(recruit) + log f(size)
#   [n_curr + j, n_next + i]      empty: no event, score 0 (and no noise when sampling)
# All other cells are -Inf. A configuration scores the sum of its events, so
# a link is chosen only when it is more likely than the death of the earlier
# stem plus the recruitment of the later one, as in the DP. With one cell per
# event the Gumbel noise of greedy_assignment_gumbel() does not favour deaths
# or recruitments by their number of copies. The matrix carries attribute
# "bd" = c(n_curr, n_next), which selects the matching sampler and fallback.
#
# Count-based slots (birth_death = FALSE, the default of this function):
# K = max(n_curr, n_next), raised only so that an assignment without forbidden
# links exists; with stable stem counts no death or recruitment slot exists
# and every stem must continue. Rows above n_curr are recruit sources and
# columns above n_next death sinks; a matrix whose K was raised carries
# attribute "k_raised" = c(K before, K).
#
# INPUTS   L          n_curr × n_next survival log-likelihoods
#                     (compute_pairwise_log_likelihood())
#          dbh_curr, dbh_next, interval_years, bio   as for L; bio also gives
#                     recruit_lambda, recruit_meanlog, recruit_sdlog and
#                     recruit_max_dbh (a larger recruit is forbidden)
# RETURNS  K × K matrix of log-scores.

augment_cost_matrix <- function(L, dbh_curr, dbh_next, interval_years, bio, K_min = NULL, birth_death = FALSE) {
    n_curr <- length(dbh_curr)
    n_next <- length(dbh_next)

    if (isTRUE(birth_death)) {
        K <- n_curr + n_next
        if (!is.null(K_min) && K < K_min) K <- as.integer(K_min)
        A <- matrix(-Inf, nrow = K, ncol = K)
        if (K > 0L) {
            if (n_curr > 0L && n_next > 0L) A[seq_len(n_curr), seq_len(n_next)] <- L
            for (i in seq_len(n_curr)) {
                hazard <- bio$h0 * exp(bio$beta_mort * dbh_curr[i])
                p_death <- max(1e-12, min(1 - 1e-12, 1 - exp(-hazard * interval_years)))
                A[i, n_next + i] <- log(p_death)
            }
            p_recruit <- max(1e-12, min(1 - 1e-12, 1 - exp(-bio$recruit_lambda * interval_years)))
            for (j in seq_len(n_next)) {
                d1 <- dbh_next[j]
                if (!is.finite(d1) || d1 <= 0) next
                ll_recruit <- log(p_recruit) + dlnorm(d1, meanlog = bio$recruit_meanlog, sdlog = bio$recruit_sdlog, log = TRUE)
                if (is.finite(bio$recruit_max_dbh) && d1 > bio$recruit_max_dbh) ll_recruit <- -Inf
                A[n_curr + j, j] <- ll_recruit
            }
            if (K > n_curr && K > n_next) A[(n_curr + 1L):K, (n_next + 1L):K] <- 0
        }
        attr(A, "bd") <- c(n_curr, n_next)
        return(A)
    }

    # Adaptive K: start from max(n_curr, n_next), then ensure enough
    # death columns for rows where ALL survival entries are -Inf, and
    # enough recruit rows for cols where ALL survival entries are -Inf.
    must_die <- 0L
    if (n_curr > 0L && n_next > 0L) {
        for (i in seq_len(n_curr)) {
            if (all(L[i, ] == -Inf)) must_die <- must_die + 1L
        }
    } else if (n_curr > 0L) {
        must_die <- n_curr
    }

    must_recruit <- 0L
    if (n_curr > 0L && n_next > 0L) {
        for (j in seq_len(n_next)) {
            if (all(L[, j] == -Inf)) must_recruit <- must_recruit + 1L
        }
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

    # Enough slots for an assignment without forbidden links. With M = the
    # largest set of allowed survival links (maximum matching on finite L),
    # the n_curr - M stems left need death columns (K - n_next) and the
    # n_next - M left need recruit rows (K - n_curr): feasible iff
    # K >= n_curr + n_next - M. The count above only looks at stems with NO
    # allowed link, so on its own it misses e.g. two stems whose only allowed
    # successor is the same stem, which would force a forbidden link. At
    # K = n_curr + n_next - M every allowed assignment has exactly M
    # survivals, which keeps the maximum-survival design of this mode; pairs
    # that are feasible at the K above keep it.
    K_before <- K
    if (n_curr > 0L && n_next > 0L && K < n_curr + n_next) {
        M <- max_allowed_matching(is.finite(L))
        if (K < n_curr + n_next - M) K <- n_curr + n_next - M
    }
    if (!is.null(K_min) && K < K_min) K <- as.integer(K_min)

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
                    dlnorm(d1,
                        meanlog = bio$recruit_meanlog,
                        sdlog = bio$recruit_sdlog, log = TRUE
                    )
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
                A[ii, jj] <- -20 # small log-prob: neither recruit nor die
            }
        }
    }

    if (K > K_before) attr(A, "k_raised") <- c(K_before, K)
    A
}

# ---- Largest set of allowed survival links ---------------------------------
# Size of a maximum matching in the bipartite graph of allowed links
# (ok[i, j] = TRUE when stem i at census c may be stem j at census c+1).
# A greedy pass gives a lower bound; the exact solve (lpSolve::lp.assign)
# runs only when that bound leaves stems unmatched on the smaller side.
max_allowed_matching <- function(ok) {
    n_curr <- nrow(ok)
    n_next <- ncol(ok)
    if (!any(ok)) {
        return(0L)
    }
    used <- logical(n_next)
    m_greedy <- 0L
    for (i in seq_len(n_curr)) {
        j <- which(ok[i, ] & !used)[1L]
        if (!is.na(j)) {
            used[j] <- TRUE
            m_greedy <- m_greedy + 1L
        }
    }
    if (m_greedy == min(n_curr, n_next)) {
        return(m_greedy)
    }
    if (!requireNamespace("lpSolve", quietly = TRUE)) {
        stop("max_allowed_matching() needs the 'lpSolve' package: install.packages(\"lpSolve\")")
    }
    n <- max(n_curr, n_next)
    w <- matrix(0, n, n)
    w[seq_len(n_curr), seq_len(n_next)] <- ok * 1
    sol <- lpSolve::lp.assign(w, direction = "max")
    if (sol$status != 0L) {
        return(m_greedy)
    }
    max(m_greedy, as.integer(round(sol$objval)))
}

# ---- Cost matrix without the hard limits ----------------------------------
# The augmented matrix A with its forbidden cells filled as the DP would
# score them before adding hard_penalty: survival links with the likelihood
# without the hard growth gates (L_free = compute_pairwise_log_likelihood()
# with no gates), recruits above the size cap with their size likelihood.
# Cells with no such value (a missing DBH) stay -Inf. Used only by
# enforce_feasible_assignment(), when no assignment avoids every forbidden
# link, so that the unavoidable violation is the most likely one.
fallback_log_cost <- function(A, L_free, dbh_next, interval_years, bio) {
    n_curr <- nrow(L_free)
    n_next <- ncol(L_free)
    K <- nrow(A)
    fb <- A
    if (n_curr > 0L && n_next > 0L) {
        blk <- fb[seq_len(n_curr), seq_len(n_next), drop = FALSE]
        forb <- !is.finite(blk)
        blk[forb] <- L_free[forb]
        fb[seq_len(n_curr), seq_len(n_next)] <- blk
    }
    if (n_curr < K && n_next > 0L) {
        # same recruit likelihood as augment_cost_matrix(), without the cap
        p_recruit <- 1 - exp(-bio$recruit_lambda * interval_years)
        p_recruit <- max(1e-12, min(1 - 1e-12, p_recruit))
        ll <- rep(-Inf, n_next)
        okd <- is.finite(dbh_next) & dbh_next > 0
        ll[okd] <- log(p_recruit) + dlnorm(dbh_next[okd],
            meanlog = bio$recruit_meanlog,
            sdlog = bio$recruit_sdlog, log = TRUE
        )
        if (!is.null(attr(A, "bd"))) {
            # birth-death matrix: each next stem has one recruit cell, [n_curr + j, j]
            for (j in seq_len(n_next)) {
                if (n_curr + j <= K && !is.finite(fb[n_curr + j, j])) fb[n_curr + j, j] <- ll[j]
            }
        } else {
            rows <- (n_curr + 1L):K
            for (j in seq_len(n_next)) {
                forb <- !is.finite(fb[rows, j])
                fb[rows[forb], j] <- ll[j]
            }
        }
    }
    fb
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
#   weight     : lookahead weight (0 disables; match_stems_probabilistic()
#                passes prob_lookahead_weight, 0.5 by default). The adjustment
#                of a column is weight × (its continuity log-likelihood minus
#                the best column's), floored at -2 before weighting.
#
# Returns: modified aug_cost matrix (same dimensions)

condition_cost_matrix <- function(aug_cost, n_curr, n_next,
                                  dbh_next, next_assignment,
                                  n_next_next, dbh_further,
                                  interval_next, bio, weight) {
    if (weight <= 0 || n_next == 0L || n_next_next == 0L) {
        return(aug_cost)
    }

    # Growth mean function (same as in compute_pairwise_log_likelihood)
    mu_growth_fn <- function(d) {
        if (!is.finite(bio$mu_gamma) || bio$mu_gamma == 0 ||
            !is.finite(d) || d <= 0) {
            return(bio$mu_const)
        }
        bio$mu_const + bio$mu_gamma * log(d)
    }

    # Pass 1: compute raw continuity log-LL for each real column j
    raw_bonus <- rep(0, n_next) # 0 = neutral (no info)
    has_info <- logical(n_next)
    for (j in seq_len(n_next)) {
        k <- next_assignment[j]
        if (k > n_next_next) next # j died — no forward info

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
    if (length(info_vals) == 0L) {
        return(aug_cost)
    } # nothing to condition on

    max_bonus <- max(info_vals)
    bonus <- rep(0, n_next)
    bonus[has_info] <- pmax(raw_bonus[has_info] - max_bonus, -2)

    # Pass 2: apply weighted normalized bonus to all feasible cells
    K <- nrow(aug_cost)
    for (j in seq_len(n_next)) {
        if (bonus[j] == 0) next # no adjustment needed
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
        if (length(j_candidates) != 1L) next # target died or ambiguous — skip
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

# apply_track_pin_mask: TrueStemID pins of any stem (not only those present at
# the anchor). A pinned obs r (pin q) at the current census must join the
# next-census obs whose track carries pin q (when exactly one does), and may
# not join an obs whose track carries a different pin. Unpinned obs and obs
# whose pin no next track carries are free (they may join an unpinned track,
# e.g. the same stem under an older, renumbered StemID, or die).
#
# INPUTS
#   cost_matrix     K×K augmented log-cost matrix
#   pin_curr        integer vector length n_curr: TrueStemID of each current
#                   obs, NA when not pinned
#   next_track_pin  integer vector length n_next: pin carried by the track of
#                   each next-census obs in this sample (NA = unpinned track)
#   n_curr, n_next  number of real observations at current / next census
#
# RETURNS  modified cost_matrix
apply_track_pin_mask <- function(cost_matrix, pin_curr, next_track_pin,
                                 n_curr, n_next) {
    if (n_curr == 0L || length(pin_curr) == 0L || all(is.na(pin_curr))) {
        return(cost_matrix)
    }
    ntp <- next_track_pin[seq_len(n_next)]
    for (r in which(!is.na(pin_curr[seq_len(n_curr)]))) {
        q <- pin_curr[r]
        other <- which(!is.na(ntp) & ntp != q)
        if (length(other) > 0L) cost_matrix[r, other] <- -Inf
        same <- which(!is.na(ntp) & ntp == q)
        if (length(same) == 1L) cost_matrix[r, -same] <- -Inf
    }
    cost_matrix
}

# propagate_track_pin: after an assignment is drawn, the pin carried by the
# track of each current-census obs: its own pin, else the pin of the next obs
# it joins (NA when it dies or joins an unpinned track).
#
# RETURNS  integer vector length n_curr
propagate_track_pin <- function(assignment, pin_curr, next_track_pin,
                                n_curr, n_next) {
    out <- as.integer(pin_curr[seq_len(n_curr)])
    for (r in seq_len(n_curr)) {
        if (!is.na(out[r])) next
        col <- assignment[r]
        if (col <= n_next) out[r] <- next_track_pin[col]
    }
    out
}

# pin_masked_pair: one census pair's cost and fallback matrices with this
# sample's pin masks (apply_pin_mask, apply_track_pin_mask) applied.
# augment_cost_matrix() sizes each pair (K) for the growth limits only, so
# that the largest set of allowed survival links fits. The pin masks can then
# forbid links that this K relied on: one stem at C6 pinned to A and one stem
# at C7 pinned to B give a 1x1 matrix (13.8 -> 15.5 cm is allowed growth, so
# no death slot) whose only link the mask forbids, so that every sample would
# have to join A to B (tag 150279); the pin sweep would then separate the
# pinned rows while the unpinned rows of that track keep its label, which can
# be another stem's pin. When the masked survival links (largest allowed set
# M) need more death/recruit slots than K has, the pair is rebuilt with
# K = n_curr + n_next - M (the rule augment_cost_matrix() applies to
# growth-forbidden links) and masked again, so one pinned stem ends and the
# other starts. Pairs with enough slots are returned masked at their own size
# (so the same random numbers are drawn for them).
#
# INPUTS  cost, fb  K×K cost (lookahead included) and fallback matrices
#         pd        pair_data entry (n_curr, n_next, dbh_curr, dbh_next, iv, L_free)
#         bio       bio parameter list
#         mask      function(m): this sample's pin masks applied to a matrix
# RETURNS list(cost, fb, grown): both matrices masked; grown = TRUE when the
#         pair got extra slots
pin_masked_pair <- function(cost, fb, pd, bio, mask) {
    mc <- mask(cost)
    mf <- mask(fb)
    n_curr <- pd$n_curr
    n_next <- pd$n_next
    if (n_curr == 0L || n_next == 0L) {
        return(list(cost = mc, fb = mf, grown = FALSE))
    }
    ok <- is.finite(mc[seq_len(n_curr), seq_len(n_next), drop = FALSE])
    K_need <- n_curr + n_next - max_allowed_matching(ok)
    if (K_need <= nrow(cost)) {
        return(list(cost = mc, fb = mf, grown = FALSE))
    }
    # survival block keeps the lookahead adjustment; slots as in augment_cost_matrix()
    # (a birth-death matrix already has a death and a recruit cell per stem,
    # so it never reaches this point)
    grown <- augment_cost_matrix(cost[seq_len(n_curr), seq_len(n_next), drop = FALSE],
        pd$dbh_curr, pd$dbh_next, pd$iv, bio,
        K_min = K_need, birth_death = !is.null(attr(cost, "bd"))
    )
    list(
        cost = mask(grown),
        fb = mask(fallback_log_cost(grown, pd$L_free, pd$dbh_next, pd$iv, bio)),
        grown = TRUE
    )
}

# ---- Gumbel-noise greedy assignment --------------------------------------
# Draw one stochastic assignment from an augmented log-cost matrix using the
# Gumbel-max trick.
#   Birth-death matrix (attribute "bd", see augment_cost_matrix()): Gumbel
#   noise on every event cell (links, deaths, recruitments), none on the
#   empty cells, and the exact best assignment of the perturbed scores
#   (hungarian_min_rcpp(), lpSolve when the compiled solver is absent).
#   Matrix without that attribute (count-based slots): each row is assigned
#   to the highest-scoring available column after adding Gumbel(0,
#   temperature) noise, processed in descending order of row maxima.
#
# INPUTS
#   log_cost_matrix  K×K matrix of log-likelihoods (augmented with
#                    mortality/recruitment slots so it is square).
#   temperature      Gumbel noise scale; higher = more random.
#   fallback         optional K×K fallback_log_cost() matrix (pin-masked like
#                    log_cost_matrix), used by enforce_feasible_assignment()
#                    only when a forbidden link cannot be avoided.
#
# RETURNS
#   Integer vector of length K: assignment[row] = assigned column index.

greedy_assignment_gumbel <- function(log_cost_matrix, temperature = 1.0, fallback = NULL) {
    K <- nrow(log_cost_matrix)
    stopifnot(ncol(log_cost_matrix) == K)

    # Add Gumbel(0, temperature) noise:  -temperature * log(-log(U))
    noise <- matrix(-temperature * log(-log(runif(K * K))), nrow = K, ncol = K)

    bd <- attr(log_cost_matrix, "bd")
    if (!is.null(bd)) {
        n_curr <- bd[1L]
        n_next <- bd[2L]
        if (K > n_curr && K > n_next) noise[(n_curr + 1L):K, (n_next + 1L):K] <- 0 # empty cells: no event, no noise
        noisy <- log_cost_matrix + noise
        w <- noisy
        w[!is.finite(w)] <- -1e9 # forbidden cells: taken only when no assignment avoids them
        assignment <- if (exists("hungarian_min_rcpp", mode = "function")) {
            as.integer(hungarian_min_rcpp(-w))
        } else {
            sol <- lpSolve::lp.assign(w - min(w), direction = "max")$solution
            as.integer(apply(sol, 1L, function(r) which(r > 0.5)[1L]))
        }
        # a forbidden cell is used only when the pair has no allowed assignment:
        # choose the least bad one as the DP would (enforce_feasible_assignment)
        return(enforce_feasible_assignment(noisy, assignment, fallback = fallback, noise = noise))
    }

    noisy <- log_cost_matrix + noise

    # Greedy assignment: for each row in descending max-noisy-score order,
    # assign to the best available column
    assignment <- integer(K) # assignment[row] = col
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

    # Greedy can be forced into a forbidden (-Inf) cell when a row has no
    # allowed column left; re-solve such a pair exactly (same noisy scores).
    enforce_feasible_assignment(noisy, assignment, fallback = fallback, noise = noise)
}

# ---- Enforce the hard limits on a sampled assignment ----------------------
# A cell scored -Inf is a forbidden link: growth outside the hard bounds (the
# same bounds as the DP), a recruit above the recruit size cap, or a pinned
# row sent to the wrong track. If the greedy assignment uses one, the pair is
# re-solved exactly (lpSolve::lp.assign) on the SAME noisy scores, with -Inf
# replaced by a large penalty so that an unavoidable violation (data that no
# assignment can satisfy) is kept to a minimum, as the DP's hard_penalty does.
# As in the DP, each forbidden link costs 1e6 (fewest violations first) plus
# its log-likelihood without the hard gate and the same Gumbel noise
# (`fallback` + `noise`), so the unavoidable violation is the most likely one
# (a mild shrink before a severe one). Cells without a fallback value (a
# pinned row sent to another track, a missing DBH) cost 2e6: a pin is kept
# rather than broken. Assignments without forbidden cells are returned
# unchanged, and no random numbers are drawn, so all other samples are
# identical to the greedy result. The result carries attr "n_forbidden" when
# it still uses a forbidden link.
#
# INPUTS   noisy       K×K matrix of perturbed log-scores (may contain -Inf)
#          assignment  K-length integer vector from the greedy pass
#          fallback    optional K×K fallback_log_cost() matrix
#          noise       the K×K Gumbel noise added to form `noisy`
# RETURNS  K-length integer vector: assignment[row] = col
enforce_feasible_assignment <- function(noisy, assignment, fallback = NULL, noise = NULL) {
    K <- nrow(noisy)
    n_forbidden <- function(a) sum(!is.finite(noisy[cbind(seq_len(K), a)]))
    if (n_forbidden(assignment) == 0L) {
        return(assignment)
    }
    if (!requireNamespace("lpSolve", quietly = TRUE)) {
        stop("enforce_feasible_assignment() needs the 'lpSolve' package: install.packages(\"lpSolve\")")
    }
    score <- noisy
    bad <- !is.finite(score)
    fb <- rep(-Inf, sum(bad))
    if (!is.null(fallback)) fb <- fallback[bad] + if (is.null(noise)) 0 else noise[bad]
    score[bad] <- ifelse(is.finite(fb), -1e6 + pmax(fb, -1e5), -2e6)
    sol <- lpSolve::lp.assign(-score, direction = "min")
    if (sol$status != 0L) {
        best <- assignment
    } else {
        exact <- max.col(sol$solution > 0.5, ties.method = "first")
        best <- if (n_forbidden(exact) <= n_forbidden(assignment)) exact else assignment
    }
    n_left <- n_forbidden(best)
    if (n_left > 0L) attr(best, "n_forbidden") <- n_left
    best
}

# ---- Repair stitched samples BEFORE marginal aggregation -----------------
# Walk each sample's trajectories and break links that violate growth
# constraints.  Two layers of defense (mirroring the DP pathway):
#
#   1. Hard-rate check: annualized growth outside [min_rate, max_rate]
#      is severed immediately.
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
# it a new unique break-ID (layer 2: the observation at the start of the
# shrinkage run).  Up to max_passes iterations per sample (a
# break can shorten trajectories and expose new violations).
#
# Pinned observations (TrueStemID, `pinned`) are never severed: the database
# identity wins over the growth limits, as the pin sweep enforces in the
# exported table. If the earlier observation is pinned, the later one is
# severed instead (unless it is pinned too or sits at the anchor census); a
# link between two pinned observations is kept.
#
# INPUTS
#   stitched     list of samples; each a list of per-census ID vectors
#   obs_data     per-census observation data (uses $dbh and $n)
#   intervals    years between consecutive observed censuses
#   min_rate, max_rate   hard bounds on annual growth (cm/yr)
#   me_sd1_a, me_sd1_b   small-error SD model of layer 2
#   n_sigma_me   layer-2 threshold in SD units (Inf turns layer 2 off)
#   max_passes   passes per sample
#   use_bio_hard_shrink   FALSE skips layer 2
#   pinned       list per census of logical vectors (TRUE = pinned
#                observation), or NULL
#
# Returns the modified stitched list (same structure), with attributes
# "sample_level_breaks" (all breaks) and "sample_level_me_breaks" (layer 2).

repair_stitched_growth_violations <- function(stitched, obs_data, intervals,
                                              min_rate, max_rate,
                                              me_sd1_a = 0.0062,
                                              me_sd1_b = 0.0904,
                                              n_sigma_me = 3,
                                              max_passes = 10L,
                                              use_bio_hard_shrink = TRUE,
                                              pinned = NULL) {
    n_samples <- length(stitched)
    n_census <- length(obs_data)
    if (n_census < 2L) {
        return(stitched)
    }
    is_pinned <- function(e) !is.null(pinned) && isTRUE(pinned[[e$ci]][e$oi])

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

    total_breaks <- 0L
    total_me_breaks <- 0L

    for (s in seq_len(n_samples)) {
        for (pass in seq_len(max_passes)) {
            breaks_this_pass <- 0L

            # Build reverse map: stem_id -> list of (ci, oi, dbh)
            traj_map <- list()
            for (ci in seq_len(n_census)) {
                ids <- stitched[[s]][[ci]]
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
                cumul_shrink <- 0
                shrink_run_start <- 1L # index into entries where current shrinkage run began
                d_run_start <- entries[[1L]]$dbh # DBH at start of shrinkage run

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
                    rate <- (d_curr - d_prev) / iv

                    # --- Layer 1: hard-rate check --------------------------
                    if (rate < min_rate || rate > max_rate) {
                        sev <- entries[[r - 1L]]
                        if (is_pinned(sev)) {
                            sev <- if (!is_pinned(entries[[r]]) && entries[[r]]$ci < n_census) entries[[r]] else NULL
                        }
                        if (!is.null(sev)) {
                            break_base <- break_base + 1L
                            stitched[[s]][[sev$ci]][sev$oi] <- break_base
                            breaks_this_pass <- breaks_this_pass + 1L
                            break # re-evaluate shortened trajectory in next pass
                        }
                        # both observations pinned: the database links them; keep the link
                    }

                    # --- Layer 2: ME cumulative-shrinkage check ------------
                    # Skip when bio hard shrink gate is disabled (mirrors DP which
                    # uses soft penalties only, never a hard ME threshold).
                    if (isTRUE(use_bio_hard_shrink) && d_curr < d_prev) {
                        cumul_shrink <- cumul_shrink + (d_prev - d_curr)
                        thresh <- n_sigma_me * sqrt(me_sd(d_run_start)^2 + me_sd(d_curr)^2)
                        if (cumul_shrink > thresh && !is_pinned(entries[[shrink_run_start]])) {
                            # Sever at the start of the shrinkage run (never a pinned observation)
                            break_base <- break_base + 1L
                            stitched[[s]][[entries[[shrink_run_start]]$ci]][entries[[shrink_run_start]]$oi] <- break_base
                            breaks_this_pass <- breaks_this_pass + 1L
                            total_me_breaks <- total_me_breaks + 1L
                            break # re-evaluate shortened trajectory
                        }
                    } else {
                        # Growth step: reset cumulative shrinkage tracker
                        cumul_shrink <- 0
                        shrink_run_start <- r
                        d_run_start <- d_curr
                    }
                }
            }

            total_breaks <- total_breaks + breaks_this_pass
            if (breaks_this_pass == 0L) break # this sample converged
        }
    }

    attr(stitched, "sample_level_breaks") <- total_breaks
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
#   min_keep    safety net: all samples are kept when fewer than
#               min(min_keep, max(1, n_samples %/% 4)) are pin-consistent
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
    if (n_samples == 0L) {
        return(stitched)
    }

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
        .msg <- paste0(
            prefix, "WARNING: pin-consistent filter would keep only ",
            n_kept, "/", n_samples, " samples (min_safe=", .min_safe,
            "); keeping ALL samples (degrading to soft-pin behavior)"
        )
        vcat(.msg)
        message(.msg)
        attr(stitched, "n_pin_filtered") <- 0L
        return(stitched)
    }

    .msg <- paste0(
        prefix, "Pin-consistent filter: kept ", n_kept, "/",
        n_samples, " samples (", n_dropped, " dropped)"
    )
    vcat(.msg)
    message(.msg)

    out <- stitched[keep]
    attr(out, "n_pin_filtered") <- n_dropped
    out
}

# ---- Stitch assignments backward from anchor ----------------------------
# Turn each sample's per-pair assignments into a track ID per observation.
# Anchor observations carry anchor_ids; walking backward, an observation
# linked to a next-census observation inherits its ID, and an observation
# sent to a death column gets an ID of its own (the same in every sample).
#
# INPUTS   all_samples  list of samples; each a list (one per census pair) of
#                       assignment vectors (assignment[row] = col)
#          obs_data     per-census observation data (uses $n)
#          anchor_ids   ID of each anchor observation
#          K            not used
# RETURNS  list of samples; each a list of per-census integer ID vectors

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
        next_obs_to_track <- anchor_ids # length n_anchor
        recon_by_census[[anchor_pos]] <- next_obs_to_track

        # Walk backward
        for (i in seq.int(anchor_pos - 1L, 1L, by = -1L)) {
            assignment <- sample_assignments[[i]] # K_pair-length: assignment[row] = col
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
# Per observation, the share of samples that give it each ID, written to
# tree_data as DP_PosteriorTop{k}ID / Prob and DP_PosteriorEntropy, and one
# ID per observation resolved census by census (ReconstructedStemID with its
# share in DP_PosteriorReconstructedProb). match_stems_probabilistic() then
# replaces that resolved ID by the consensus sample
# (select_consensus_trajectory()).
#
# INPUTS   stitched         list of samples; each a list of per-census IDs
#          tree_data        the tag's data.table (modified by reference)
#          obs_data         per-census observation data ($n, $idx, $dbh)
#          obs_census       not used
#          anchor_pos       position of the anchor census in obs_data
#          posterior_top_k  number of Top-k columns to fill
#          intervals, min_rate, max_rate   when all three are given, a
#                           candidate ID is skipped if it breaks the growth
#                           bounds against the stem's nearest resolved census
# RETURNS  tree_data

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
    if (!("DP_PosteriorEntropy" %in% names(tree_data))) {
        tree_data[, DP_PosteriorEntropy := NA_real_]
    }
    if (!("DP_PosteriorReconstructedProb" %in% names(tree_data))) {
        tree_data[, DP_PosteriorReconstructedProb := NA_real_]
    }

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
            .after <- anchor_pos + .d
            if (.before >= 1L) anchor_out_order <- c(anchor_out_order, .before)
            if (.after <= n_census) anchor_out_order <- c(anchor_out_order, .after)
        }
        # Per-stem resolved track: stem_id_string -> list of (ci, dbh) entries
        .stem_tracks <- new.env(hash = TRUE, parent = emptyenv())

        # Interval between two census positions (handles gaps)
        .census_iv <- function(ci_a, ci_b) {
            lo <- min(ci_a, ci_b)
            hi <- max(ci_a, ci_b)
            if (lo == hi) {
                return(0)
            }
            idx_rng <- lo:(hi - 1L)
            idx_rng <- idx_rng[idx_rng <= length(intervals)]
            if (length(idx_rng) == 0L) {
                return(5.0)
            }
            iv <- sum(intervals[idx_rng])
            if (!is.finite(iv) || iv <= 0) 5.0 else iv
        }

        # Check growth rate against nearest resolved assignment on each side
        .growth_ok <- function(sid_key, ci_new, dbh_new) {
            trk <- .stem_tracks[[sid_key]]
            if (is.null(trk)) {
                return(TRUE)
            }
            trk_cis <- vapply(trk, function(x) x$ci, numeric(1))
            # Closest already-resolved census BEFORE ci_new
            below_idx <- which(trk_cis < ci_new)
            if (length(below_idx) > 0L) {
                j <- below_idx[which.max(trk_cis[below_idx])]
                iv <- .census_iv(trk_cis[j], ci_new)
                rate <- (dbh_new - trk[[j]]$dbh) / iv
                if (rate < min_rate || rate > max_rate) {
                    return(FALSE)
                }
            }
            # Closest already-resolved census AFTER ci_new
            above_idx <- which(trk_cis > ci_new)
            if (length(above_idx) > 0L) {
                j <- above_idx[which.min(trk_cis[above_idx])]
                iv <- .census_iv(ci_new, trk_cis[j])
                rate <- (trk[[j]]$dbh - dbh_new) / iv
                if (rate < min_rate || rate > max_rate) {
                    return(FALSE)
                }
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

        used_ids <- integer(0)
        resolved_ids <- integer(n_obs)
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
                        !.growth_ok(as.character(cand_id), ci, dbh_oi)) {
                        next
                    }
                }
                resolved_ids[oi] <- cand_id
                resolved_probs[oi] <- post$probs[j]
                used_ids <- c(used_ids, cand_id)
                assigned <- TRUE
                break
            }
            if (!assigned) {
                # All posterior alternatives are taken or growth-violated —
                # assign a new unique (break) ID
                new_id <- max(c(
                    tree_data$ReconstructedStemID, used_ids,
                    resolved_ids
                ), na.rm = TRUE) + 1L
                if (!is.finite(new_id)) new_id <- 1L
                resolved_ids[oi] <- new_id
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
            data.table::set(
                tree_data, tree_data_row, "ReconstructedStemID",
                resolved_ids[oi]
            )
            data.table::set(
                tree_data, tree_data_row, "DP_PosteriorReconstructedProb",
                resolved_probs[oi]
            )

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

# ---- Consensus trajectory from the posterior samples ---------------------
# The exported reconstruction must be ONE coherent trajectory. Resolving each
# census separately from marginal probabilities (above) mixes samples and,
# when marginals are diffuse, fragments stems (one-census stems, gaps that
# are later gap-filled and double counted).
#
# This picks the sample with maximum expected accuracy: for every observation
# its predecessor in the sample (the observation with the same track ID at the
# previous observed census, or none) is scored by how often the samples choose
# that same predecessor; the sample with the highest total wins. Ties go to
# the most frequent partition, then to the lowest sample index. When one
# partition dominates the samples it is also the MEA sample.
#
# INPUTS   stitched  list of samples; each a list of per-census ID vectors
#          obs_data  per-census observation data (uses $n)
# RETURNS  list(index = chosen sample, agreement = mean predecessor frequency
#          of the chosen sample, part_freq = frequency of its partition)
select_consensus_trajectory <- function(stitched, obs_data) {
    n_s <- length(stitched)
    n_census <- length(obs_data)
    n_obs <- vapply(obs_data, function(x) as.integer(x$n), integer(1))
    score <- numeric(n_s)
    n_scored <- 0L
    if (n_census >= 2L) {
        for (ci in 2:n_census) {
            if (n_obs[ci] == 0L) next
            # predecessor index (0 = none) of each observation in each sample
            pred <- vapply(stitched, function(s) {
                match(s[[ci]][seq_len(n_obs[ci])], s[[ci - 1L]][seq_len(n_obs[ci - 1L])], nomatch = 0L)
            }, integer(n_obs[ci]))
            if (is.null(dim(pred))) pred <- matrix(pred, nrow = 1L)
            for (oi in seq_len(n_obs[ci])) {
                freq <- tabulate(pred[oi, ] + 1L, nbins = n_obs[ci - 1L] + 1L) / n_s
                score <- score + freq[pred[oi, ] + 1L]
            }
            n_scored <- n_scored + n_obs[ci]
        }
    }
    sig <- vapply(stitched, function(s) {
        ids <- unlist(lapply(seq_len(n_census), function(ci) s[[ci]][seq_len(n_obs[ci])]))
        paste(match(ids, unique(ids)), collapse = ",")
    }, character(1))
    part_freq <- as.numeric(table(sig)[sig]) / n_s
    best <- order(-round(score, 10), -part_freq, seq_len(n_s))[1L]
    list(
        index = best,
        agreement = if (n_scored > 0L) score[best] / n_scored else 1,
        part_freq = part_freq[best]
    )
}

# ---- Repair growth violations from marginal resolution -------------------
# Not called by match_stems_probabilistic(), which exports the consensus
# sample and only counts residual violations.
#
# The per-census greedy marginal resolution can assign the same StemID to
# observations at consecutive censuses that violate growth bounds.  This
# happens because the marginals are computed independently per census.
#
# Repair: walk each StemID's trajectory.  When the growth between two
# consecutive measured rows is outside [min_rate, max_rate], the earlier
# row gets a new unique ID; up to 10 passes.  obs_data and posterior_top_k
# are not used.  Returns tree_data (modified by reference).

repair_marginal_growth_violations <- function(tree_data, obs_data, obs_census,
                                              intervals, min_rate, max_rate,
                                              posterior_top_k, vcat, prefix) {
    n_census <- length(obs_census)
    if (n_census < 2L) {
        return(tree_data)
    }

    # Build date-based interval lookup for arbitrary census pairs
    dt_dates <- tree_data[, .(MeanDate = mean(as.numeric(as.Date(ExactDate)),
        na.rm = TRUE
    )), by = CensusID]

    max_id <- suppressWarnings(max(tree_data$ReconstructedStemID, na.rm = TRUE))
    if (!is.finite(max_id)) max_id <- 0L
    total_breaks <- 0L
    max_passes <- 10L # safety limit

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
                iv <- if (length(md0) > 0L && length(md1) > 0L) {
                    (md1[1] - md0[1]) / 365.25
                } else {
                    5.0
                }
                if (!is.finite(iv) || iv <= 0) iv <- 5.0

                rate <- (dsub$DBH[r] - dsub$DBH[r - 1L]) / iv
                if (rate >= min_rate && rate <= max_rate) next

                # Violation: break the earlier census obs to a new unique ID
                row_prev <- dsub$row_idx[r - 1L]
                max_id <- max_id + 1L
                data.table::set(
                    tree_data, as.integer(row_prev),
                    "ReconstructedStemID", as.integer(max_id)
                )
                breaks_this_pass <- breaks_this_pass + 1L
                # After breaking, don't check further pairs for this stem
                # in this pass — re-evaluate in next pass
                break
            }
        }

        total_breaks <- total_breaks + breaks_this_pass
        if (breaks_this_pass == 0L) break # converged
    }

    if (total_breaks > 0L) {
        .warn_msg <- paste0(
            prefix, "WARNING: Post-marginal safety-net repair fired ",
            total_breaks, " break(s) — probabilities for these rows are stale ",
            "(greedy conflict resolution created new violations)"
        )
        vcat(.warn_msg)
        message(.warn_msg) # ensure it appears on stderr / captured by log redirection
    }
    tree_data
}

# ---- Export posterior samples (mirrors DP format) -------------------------
# Build the long samples table of a tag (Tag, Sample, CensusID,
# ReconstructedStemID, ObsRowID; one row per measured observation and sample,
# no logp column) and stage it as
# <posterior_samples_path>/posteriors/.staging/tag_<Tag>_samples_raw_<ts>.rds
# for finalize_posterior_paths(). <ts> is the global BATCH_TS when it exists,
# otherwise the current time; a NULL posterior_samples_path falls back to the
# global `out_dir`, then to tempdir(). posterior_samples_format ("rds",
# "feather" or "csv") is stored in the staging file as the format of the
# final paths file. `verbose` is not used (messages go through `vcat`).
#
# stage = FALSE: build samples_dt and return it without writing the staging
# file (used when a DP resprout-split segment falls back to this matcher; the
# split pairs the two segments' draws and stages them once per tag).
# Returns samples_dt invisibly.
export_probabilistic_posteriors <- function(stitched, tree_data, obs_data,
                                            obs_census, tag_val, n_samples,
                                            posterior_samples_path,
                                            posterior_samples_format,
                                            verbose, prefix, vcat,
                                            stage = TRUE) {
    fmt <- match.arg(posterior_samples_format, c("rds", "feather", "csv"))

    out_dir_local <- if (!is.null(posterior_samples_path)) {
        posterior_samples_path
    } else {
        get0("out_dir", ifnotfound = NULL)
    }
    if (is.null(out_dir_local) || !nzchar(out_dir_local)) out_dir_local <- tempdir()

    out_dir_post <- file.path(out_dir_local, "posteriors")
    if (isTRUE(stage) && !dir.exists(out_dir_post)) dir.create(out_dir_post, recursive = TRUE, showWarnings = FALSE)

    ts_local <- get0("BATCH_TS", ifnotfound = format(Sys.time(), "%Y%m%d_%H%M%S"))
    out_path_base <- file.path(
        out_dir_post,
        paste0(
            "tag_", ifelse(is.na(tag_val), "NA", tag_val),
            "_posterior_samples_", ts_local
        )
    )

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
    if (!isTRUE(stage)) {
        return(invisible(samples_dt))
    }

    # Stage raw samples_dt (engine ID space). The post-engine pipeline will
    # translate via renumber mapping, re-run apply_bb_invariants_to_samples
    # in renumbered space, compute path_sig / paths_summary and write the
    # final paths file. See dp_global_main.R::finalize_posterior_paths().
    staging_dir <- file.path(out_dir_post, ".staging")
    if (!dir.exists(staging_dir)) dir.create(staging_dir, recursive = TRUE, showWarnings = FALSE)
    staging_path <- file.path(staging_dir, paste0(
        "tag_", ifelse(is.na(tag_val), "NA", tag_val),
        "_samples_raw_", ts_local, ".rds"
    ))
    saveRDS(list(
        engine = "probabilistic",
        tag_val = tag_val,
        batch_ts = ts_local,
        posterior_samples_format = fmt,
        posterior_samples_path = out_dir_local,
        samples_dt = samples_dt,
        sampling_profile = list(
            posterior_samples = n_samples,
            started = Sys.time(), finished = Sys.time()
        )
    ), file = staging_path)
    vcat(prefix, sprintf(
        "Posterior sampling: staged %d samples for tag %s to %s",
        n_samples, as.character(tag_val), staging_path
    ))
    invisible(samples_dt)
}
