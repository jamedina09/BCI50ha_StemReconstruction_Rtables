############################################################
# dp_global_dp.R
# Core dynamic programming (MAP and marginal DP functions)
############################################################

# match_stems_dp_global_backward_marginals_batch()
#
# PURPOSE
#   Run backward marginal DP stem identification for one or more tags.  For
#   each tag (Tag × species group) the function tries the full DP solver; if
#   the state space exceeds `max_states` it falls back to a cheaper method
#   (probabilistic greedy matching).
#
# KEY PARAMETERS
#   tree_data            data.table — one row per (Tag, OriginalStemID, CensusID).
#   min_growth / max_growth  Hard bounds on annual growth (cm/yr).
#   anchor_start         Integer CensusID where TrueStemID is trusted.
#   max_states           Max states per census before fallback (default 50 000).
#   temperature          Marginal-DP softmax temperature; 1.0 = Bayesian.
#   use_measurement_error  Use the mixture measurement-error model in likelihoods.
#   posterior_samples    How many posterior path samples to draw (0 = none).
#   prob_n_samples       Number of stochastic samples for the probabilistic
#                        greedy matcher (used on fallback).
#   prob_species         Character vector of species to route to the probabilistic
#                        matcher instead of the default DP fallback.
#   prob_lookahead_weight  Weight [0,1] for sequential backward conditioning in the
#                        probabilistic matcher.  0 = independent per-pair sampling
#                        (original behaviour).  0.5 = default (blend future
#                        assignment info into the cost matrix).
#
# RETURNS
#   data.table with ReconstructedStemID, ReconstructionMethod, marginal columns,
#   and optional posterior path files written to posterior_samples_path.

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
                                                           # Allow DP to use a provisional anchor at the last observed DBH census when no TrueStemID exists
                                                           allow_provisional_anchor = TRUE,
                                                           use_measurement_error = FALSE,
                                                           meas_sd1_a = 0.0062,
                                                           meas_sd1_b = 0.0904,
                                                           meas_sd2 = 4.64,
                                                           meas_p_big = 0.05,
                                                           # --- growth-form based fallback ---
                                                           # vector of values in `growth_form` column that should trigger
                                                           # immediate probabilistic fallback and avoid DP entirely
                                                           fallback_growth_forms = character(0),
                                                           # --- posterior sampling options ---
                                                           posterior_samples = 0L,
                                                           posterior_samples_format = c("rds", "feather", "csv"),
                                                           posterior_samples_path = NULL,
                                                           posterior_sample_seed = NULL,
                                                           # --- pruning options (conservative hard guards) ---
                                                           prune_hard = TRUE,
                                                           prune_min_growth = NULL,
                                                           prune_max_growth = NULL,
                                                           prune_use_bio_bounds = TRUE,
                                                           prune_recruit_max_dbh = NULL,
                                                           prune_use_bio_recruit = TRUE,
                                                           # --- non-taper-corrected growth form prune bounds ---
                                                           # Growth forms whose DBH measurements are NOT taper-corrected
                                                           # (palms, strangler figs, tree ferns) exhibit real DBH growth
                                                           # plus large apparent variation from HOM changes. They receive
                                                           # wide base prune bounds (default 1.25× standard limits) and
                                                           # optional HOM-proportional widening per census pair.
                                                           non_taper_corrected_growth_forms = c("palm", "strangler_fig", "tree_fern"),
                                                           non_taper_corrected_prune_min_growth = -0.625,
                                                           non_taper_corrected_prune_max_growth = 6.25,
                                                           # HOM tolerance scale: cm of annual DBH tolerance per meter
                                                           # of HOM deviation from 1.3 m.  Set 0 to disable HOM widening.
                                                           hom_tolerance_scale = 2.0,
                                                           verbose = FALSE,
                                                           chunk_id = NULL,
                                                           allow_segment_split = TRUE,
                                                           post_segment_all_recruits = FALSE,
                                                           # --- probabilistic matcher options ---
                                                           prob_n_samples = 200L, # stochastic samples for probabilistic fallback
                                                           prob_species = character(0), # species routed to probabilistic matcher
                                                           prob_lookahead_weight = 0.5, # backward conditioning weight [0,1]
                                                           use_bio_hard_shrink_in_prob = TRUE, # use bio hard shrink gate in probabilistic
                                                           use_bio_hard_growth_in_prob = TRUE, # use bio hard growth gate in probabilistic
                                                           prob_n_sigma_me = 3,            # ME cumulative-shrinkage threshold for probabilistic matcher
                                                           # --- TrueStemID pinning at non-anchor censuses ---
                                                           pin_truestemid = TRUE) # pin obs with known TrueStemID to their track
{
    # Derive max_edges from max_states: the cross-product of two adjacent
    # census state counts can be at most max_states^2.  Using that as the
    # edge limit keeps a single user-facing knob (max_states) controlling
    # both per-census enumeration and inter-census transition budgets.
    max_edges <- as.double(max_states) * as.double(max_states)
    # Safety
    posterior_top_k <- as.integer(posterior_top_k)
    if (!is.finite(posterior_top_k) || is.na(posterior_top_k) || posterior_top_k < 1L) {
        posterior_top_k <- 1L
    }
    temperature <- suppressWarnings(as.numeric(temperature))
    if (!is.finite(temperature) || is.na(temperature) || temperature <= 0) {
        stop("temperature must be a positive finite number")
    }
    # ensure growth form list is character vector
    fallback_growth_forms <- if (is.null(fallback_growth_forms)) character(0) else as.character(fallback_growth_forms)
    # ensure prob_species list is character vector
    prob_species <- if (is.null(prob_species)) character(0) else as.character(prob_species)
    if (length(prob_species) == 1L && grepl("[,;]", prob_species)) {
        ps <- strsplit(prob_species, "[,;]")[[1L]]
        ps <- trimws(ps)
        ps <- ps[nzchar(ps)]
        prob_species <- ps
    }
    # interpret comma/semicolon-separated values in a single string element
    if (length(fallback_growth_forms) == 1L && grepl("[,;]", fallback_growth_forms)) {
        # split on commas or semicolons, trim whitespace, and drop empty
        ff <- strsplit(fallback_growth_forms, "[,;]")[[1L]]
        ff <- trimws(ff)
        ff <- ff[nzchar(ff)]
        fallback_growth_forms <- ff
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

    # Defensive early prune_stats so finalize_out (used by early returns) can always attach a value
    prune_stats <- list(
        total_examined = 0L,
        total_pruned = 0L,
        per_census = integer(0)
    )

    # Compute profiling accumulators
    transition_cost_time <- 0
    transition_cost_calls <- 0L

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
        if (!is.null(chunk_id)) paste0(" chunk=", chunk_id) else "",
        if (!is.na(tag_val)) paste0(" Tag=", tag_val) else "",
        if (!is.na(sp_val)) paste0(" species=", sp_val) else "",
        "] "
    )

    tree_data <- tree_data[order(CensusID)]

    # Helpers
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
        # Fast path: all expected columns already present with correct types.
        # After the initial setup call, subsequent calls (fallback paths, finalize_out)
        # return immediately without any data.table operations.
        .nms <- names(dt)
        .ptk <- max(1L, as.integer(posterior_top_k))
        .nids <- paste0("DP_PosteriorTop", seq_len(.ptk), "ID")
        .npbs <- paste0("DP_PosteriorTop", seq_len(.ptk), "Prob")
        if (all(c(
            .nids, .npbs, "DP_PosteriorEntropy", "DP_PosteriorReconstructedProb",
            "DP_PosteriorUnlinkedProb", "obs_row_id"
        ) %in% .nms) &&
            all(vapply(.nids, function(n) is.integer(dt[[n]]), logical(1L))) &&
            all(vapply(.npbs, function(n) is.numeric(dt[[n]]), logical(1L))) &&
            is.numeric(dt$DP_PosteriorEntropy) &&
            is.numeric(dt$DP_PosteriorReconstructedProb) &&
            is.numeric(dt$DP_PosteriorUnlinkedProb) &&
            is.integer(dt$obs_row_id)) {
            return(dt)
        }
        # Slow path: add or coerce missing/wrong-typed columns.
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
    # Defensive: add `DP_FallbackReason` column so downstream early-returns can set it
    if (!("DP_FallbackReason" %in% names(tree_data))) tree_data[, DP_FallbackReason := NA_character_]

    # Defensive initialization: ensure core tracking columns exist with correct types
    if (!("TrueStemID" %in% names(tree_data))) tree_data[, TrueStemID := as.integer(NA_integer_)]
    if (!("ReconstructedStemID" %in% names(tree_data))) tree_data[, ReconstructedStemID := as.integer(NA_integer_)]
    if (!("ReconstructionMethod" %in% names(tree_data))) tree_data[, ReconstructionMethod := NA_character_]
    if (!("ConstraintViolation" %in% names(tree_data))) tree_data[, ConstraintViolation := as.logical(rep(NA, .N))]

    # ---- MF (Missing From Field) pre-processing ----
    # Detect MF episodes and stash affected rows before the DP runs.
    # MF rows (is.na(DBH) & "MF" in ListOfTSM) represent stems confirmed alive
    # but unmeasured.  Removing them prevents the DP from interpreting the gap
    # as a spurious death-recruitment pair.
    #
    # Two detection modes (both census-level, no per-stem identifier required):
    #
    #   A. Explicit MF: rows where is.na(DBH) & ListOfTSM contains "MF".
    #      Forward propagation extends through subsequent consecutive censuses
    #      where ALL rows have is.na(DBH).
    #
    #   B. Implicit MF: censuses where ALL rows have is.na(DBH) (regardless of
    #      ListOfTSM) but the tag has at least one non-NA DBH in both an earlier
    #      and a later census.  These are sandwiched gap censuses that indicate
    #      the stem was alive but unmeasured, even without the MF code.
    mf_stash <- NULL
    has_mf_stash <- FALSE
    .n_input_rows <- nrow(tree_data) # total input rows before any MF stashing

    # Collect row indices to stash (union of explicit and implicit MF)
    episode_idx <- integer(0)

    # --- A. Explicit MF detection (ListOfTSM contains "MF", or Status == "missing") ---
    {
        is_mf_from_tsm <- if ("ListOfTSM" %in% names(tree_data)) {
            !is.na(tree_data$ListOfTSM) &
                grepl("\\bMF\\b", tree_data$ListOfTSM, perl = TRUE)
        } else {
            rep(FALSE, nrow(tree_data))
        }
        is_mf_from_status <- if ("Status" %in% names(tree_data)) {
            !is.na(tree_data$Status) & tree_data$Status == "missing"
        } else {
            rep(FALSE, nrow(tree_data))
        }
        is_mf_anchor <- is.na(tree_data$DBH) & (is_mf_from_tsm | is_mf_from_status)

        if (any(is_mf_anchor)) {
            episode_idx <- which(is_mf_anchor)

            # Census-level forward propagation: from each MF-anchor census,
            # extend through subsequent consecutive censuses where ALL rows
            # have NA DBH (stem still missing, regardless of ListOfTSM).
            all_censuses_sorted <- sort(unique(tree_data$CensusID))
            mf_anchor_censuses <- sort(unique(tree_data$CensusID[is_mf_anchor]))

            for (mf_c in mf_anchor_censuses) {
                later <- all_censuses_sorted[all_censuses_sorted > mf_c]
                for (next_c in later) {
                    rows_at_next <- which(tree_data$CensusID == next_c)
                    if (length(rows_at_next) == 0L) break
                    if (all(is.na(tree_data$DBH[rows_at_next]))) {
                        episode_idx <- c(episode_idx, rows_at_next)
                    } else {
                        break
                    }
                }
            }
        }
    }

    # --- B. Implicit MF detection (sandwiched all-NA censuses) ---
    # A census is implicitly missing if:
    #   1. ALL rows at that census have is.na(DBH)
    #   2. There exists at least one earlier census with non-NA DBH
    #   3. There exists at least one later census with non-NA DBH
    #   4. NO rows at that census carry a resprout/breakage code (R, RP, RF, RT, QR, OR)
    #      — those censuses mark a genuine biological event, not a measurement gap
    resprout_regex_mf <- "\\b(R|RP|RF|RT|QR|OR)\\b"
    all_censuses_sorted <- sort(unique(tree_data$CensusID))
    censuses_with_dbh <- sort(unique(tree_data$CensusID[!is.na(tree_data$DBH)]))

    if (length(censuses_with_dbh) >= 2L) {
        first_measured <- min(censuses_with_dbh)
        last_measured <- max(censuses_with_dbh)

        for (cc in all_censuses_sorted) {
            # Only consider censuses strictly between the first and last measured
            if (cc <= first_measured || cc >= last_measured) next
            rows_at_cc <- which(tree_data$CensusID == cc)
            if (length(rows_at_cc) == 0L) next
            if (all(is.na(tree_data$DBH[rows_at_cc]))) {
                # Skip if any row has a resprout/death code — that is a real event,
                # not a measurement gap, and must not be stashed as MF
                has_resprout_code <- (
                    "ListOfTSM" %in% names(tree_data) &&
                    any(!is.na(tree_data$ListOfTSM[rows_at_cc]) &
                        grepl(resprout_regex_mf, tree_data$ListOfTSM[rows_at_cc], perl = TRUE))
                ) | (
                    "Status" %in% names(tree_data) &&
                    any(!is.na(tree_data$Status[rows_at_cc]) &
                        tree_data$Status[rows_at_cc] == "broken below")
                )
                if (has_resprout_code) next
                # Sandwiched all-NA census with no resprout code → implicit MF
                episode_idx <- c(episode_idx, rows_at_cc)
            }
        }
    }

    # --- Stash collected MF rows (explicit + implicit) ---
    episode_idx <- sort(unique(episode_idx))

    if (length(episode_idx) > 0L) {
        mf_stash <- data.table::copy(tree_data[episode_idx])
        has_mf_stash <- TRUE
        tree_data <- tree_data[-episode_idx]
        vcat(prefix, "Stashing ", nrow(mf_stash), " missing-from-field (MF) row(s) from census(es) ", paste(sort(unique(mf_stash$CensusID)), collapse = ", "), " (will re-insert after DP)")

        # Identify censuses that became fully empty after MF removal
        mf_emptied_censuses <- integer(0)
        for (cc in unique(mf_stash$CensusID)) {
            if (nrow(tree_data[CensusID == cc & !is.na(DBH)]) == 0L) {
                mf_emptied_censuses <- c(mf_emptied_censuses, as.integer(cc))
            }
        }
        if (length(mf_emptied_censuses) > 0L) {
            vcat(
                prefix, "Censuses with no remaining observations after MF removal: ",
                paste(mf_emptied_censuses, collapse = ", ")
            )
        }
    }

    # Helper: re-insert stashed MF episode rows after the DP.
    # Attempts to infer ReconstructedStemID by finding tracks that are alive
    # in the flanking censuses but not observed at the MF census.  If exactly
    # one such candidate exists the ID is assigned; otherwise it stays NA.
    reinsert_mf_rows <- function(out, stash) {
        stash <- ensure_posterior_columns(stash)
        if (!("ReconstructionMethod" %in% names(stash))) stash[, ReconstructionMethod := NA_character_]
        if (!("ReconstructedStemID" %in% names(stash))) stash[, ReconstructedStemID := NA_integer_]
        if (!("ConstraintViolation" %in% names(stash))) stash[, ConstraintViolation := as.logical(NA)]
        stash[, ReconstructionMethod := "dp_mf_inferred"]

        # Best-effort matching: for each MF row at census C, look at tracks
        # assigned in the nearest census before C and the nearest census after C.
        # The candidate is the set of tracks present in both flanking sets BUT
        # not assigned to any observation at census C in `out`.
        # Each MF row at a given census consumes one candidate (removed from the
        # pool for subsequent MF rows at the same census) so that when there are
        # N MF rows and exactly N missing tracks the assignment is unambiguous.
        out_assigned <- out[!is.na(ReconstructedStemID)]
        dp_censuses <- sort(unique(out_assigned$CensusID))

        for (mf_c in sort(unique(stash$CensusID))) {
            mf_rows <- which(stash$CensusID == mf_c)
            if (length(mf_rows) == 0L) next

            # IDs assigned at census C in out (observed stems that were not MF)
            ids_at_c <- unique(out_assigned[CensusID == mf_c, ReconstructedStemID])

            # Nearest census with assignments before and after MF census
            before_c <- dp_censuses[dp_censuses < mf_c]
            after_c <- dp_censuses[dp_censuses > mf_c]
            ids_before <- if (length(before_c) > 0L) {
                unique(out_assigned[CensusID == max(before_c), ReconstructedStemID])
            } else {
                integer(0)
            }
            ids_after <- if (length(after_c) > 0L) {
                unique(out_assigned[CensusID == min(after_c), ReconstructedStemID])
            } else {
                integer(0)
            }

            # Candidate tracks: present in at least one flanking census but
            # not assigned at the MF census itself.
            if (length(ids_before) > 0L && length(ids_after) > 0L) {
                flanking <- intersect(ids_before, ids_after)
            } else if (length(ids_before) > 0L) {
                flanking <- ids_before
            } else {
                flanking <- ids_after
            }
            candidates <- setdiff(flanking, ids_at_c)

            # Assign one candidate per MF row (consuming from pool)
            for (ri in mf_rows) {
                if (length(candidates) == 1L) {
                    data.table::set(
                        stash, ri, "ReconstructedStemID",
                        as.integer(candidates[[1L]])
                    )
                    candidates <- candidates[-1L]
                } else if (length(candidates) > 1L) {
                    # Multiple candidates — cannot disambiguate; leave NA
                }
            }
        }

        # Ensure DP metadata columns exist on stash
        for (col in c("DP_KUsed", "DP_MaxStatesPerCensus", "DP_MaxStatesCensusID")) {
            if (!(col %in% names(stash))) stash[, (col) := NA_integer_]
        }
        if (!("DP_FallbackReason" %in% names(stash))) stash[, DP_FallbackReason := NA_character_]

        out <- data.table::rbindlist(list(out, stash), use.names = TRUE, fill = TRUE)
        if ("obs_row_id" %in% names(out)) data.table::setorder(out, obs_row_id)
        out
    }

    # Helper: mark post-anchor rows as 'given' when appropriate and normalize post rows
    propagate_post_anchor_given <- function(post) {
        if (!("ReconstructedStemID" %in% names(post))) post[, ReconstructedStemID := as.integer(NA_integer_)]
        if (!("ReconstructionMethod" %in% names(post))) post[, ReconstructionMethod := NA_character_]
        if (!("ConstraintViolation" %in% names(post))) post[, ConstraintViolation := as.logical(NA)]

        # Trust all observed TrueStemID at post-anchor censuses (live or dead)
        post[!is.na(TrueStemID), `:=`(
            ReconstructedStemID = as.integer(TrueStemID),
            ReconstructionMethod = "given"
        )]
        # Also assign broken-below / stump rows with known TrueStemID (NA-DBH R-coded rows)
        .has_status_pa <- "Status" %in% names(post)
        .has_tsm_pa    <- "ListOfTSM" %in% names(post)
        if (.has_status_pa || .has_tsm_pa) {
            .resprout_re_pa <- "\\b(R|RP|RF|RT|OR)\\b"
            .is_bb_pa    <- if (.has_status_pa) !is.na(post$Status) & post$Status == "broken below" else rep(FALSE, nrow(post))
            .is_r_tsm_pa <- if (.has_tsm_pa) !is.na(post$ListOfTSM) & grepl(.resprout_re_pa, post$ListOfTSM, perl = TRUE) else rep(FALSE, nrow(post))
            .stump_rows_pa <- which(is.na(post$ReconstructionMethod) & is.na(post$DBH) &
                                    !is.na(post$TrueStemID) & (.is_bb_pa | .is_r_tsm_pa))
            if (length(.stump_rows_pa) > 0L) {
                post[.stump_rows_pa, `:=`(
                    ReconstructedStemID = as.integer(TrueStemID),
                    ReconstructionMethod = "given"
                )]
            }
        }
        # Default remaining post-anchor rows to none_after_anchor
        post[is.na(ReconstructionMethod), ReconstructionMethod := "none_after_anchor"]
        post <- ensure_posterior_columns(post)
        post[, `:=`(
            DP_KUsed = NA_integer_,
            DP_MaxStatesPerCensus = NA_real_,
            DP_MaxStatesCensusID = NA_integer_
        )]
        post
    }

    # Helper: finalize output by appending post-anchor rows (if DP was scoped)
    finalize_out <- function(out) {
        out <- ensure_posterior_columns(out)
        # Restore original TrueStemID for provisional_dp rows (undo fabricated anchor IDs
        # that the provisional anchor mechanism wrote into tree_data during segment splits).
        .prov_rows <- which(out$ReconstructionMethod == "provisional_dp")
        if (length(.prov_rows) > 0L &&
            "obs_row_id" %in% names(out) &&
            "obs_row_id" %in% names(original_tree_data)) {
            .orig_tsid <- original_tree_data$TrueStemID[
                match(out$obs_row_id[.prov_rows], original_tree_data$obs_row_id)
            ]
            out[.prov_rows, TrueStemID := .orig_tsid]
        }
        if (isTRUE(dp_scoped_to_pre_anchor)) {
            post <- original_tree_data[CensusID > anchor_start]
            if (nrow(post) > 0L) {
                post <- propagate_post_anchor_given(post)
                # Propagate DP_FallbackReason to post rows when present on out
                if (("DP_FallbackReason" %in% names(out))) {
                    fb <- unique(na.omit(out$DP_FallbackReason))
                    if (length(fb) == 1L) {
                        post[, DP_FallbackReason := fb]
                    } else if (length(fb) > 1L) {
                        post[, DP_FallbackReason := paste(unique(fb), collapse = ";")]
                    }
                }
                out <- data.table::rbindlist(list(out, post), use.names = TRUE, fill = TRUE)
            }
        }
        # Re-insert stashed MF episode rows (if any)
        if (isTRUE(has_mf_stash) && !is.null(mf_stash) && nrow(mf_stash) > 0L) {
            out <- reinsert_mf_rows(out, mf_stash)
        }
        # ---- Hard-invariant FINAL sweep (terminal) -------------------------
        # Authoritative enforcement of the invariant
        #   ReconstructedStemID == TrueStemID  on every row with non-NA TrueStemID
        # (excluding provisional_dp rows whose TrueStemID is fabricated).
        # Runs in both DP-success and fallback paths (do_fallback) because
        # every code path ends in finalize_out.  Critically, it runs AFTER
        # the NA-R / R-boundary barrier severing inside do_fallback, which
        # would otherwise overwrite pinned IDs with synthetic ones to break
        # crossing tracks.  Synthetic IDs win for rows the DP/probabilistic
        # matcher freely chose; the supplied TrueStemID always wins where
        # it is non-NA.
        if (isTRUE(pin_truestemid) && "TrueStemID" %in% names(out)) {
            # Note: provisional_dp rows are NOT excluded here because the
            # restore block above (.prov_rows) has already reverted their
            # TrueStemID to the original input value (NA for rows that had
            # no real TrueStemID).  So `!is.na(out$TrueStemID)` correctly
            # identifies only rows where the user supplied a real ID.
            .has_true_final <- !is.na(out$TrueStemID)
            if (any(.has_true_final)) {
                # ---- Audit: detect engine-vs-pin disagreements -------------
                # If a row already had a non-NA ReconstructedStemID assigned
                # by the DP/probabilistic engine that differs from TrueStemID,
                # the sweep is silently overriding an engine decision.  Flag
                # these rows so downstream consumers can identify cases where
                # the exported posteriors describe a different ID than the
                # final ReconstructedStemID.
                if (!("SweepAuditOverride" %in% names(out))) {
                    out[, SweepAuditOverride := FALSE]
                } else {
                    # Per-row backfill: rows appended after an earlier inner sweep
                    # (segment merges, post-anchor block via rbindlist fill=TRUE)
                    # arrive with NA; set them to FALSE so only genuine overrides
                    # remain TRUE.
                    out[is.na(SweepAuditOverride), SweepAuditOverride := FALSE]
                }
                # Snapshot the engine's pre-sweep ReconstructedStemID.
                # Per-row backfill: create the column if absent, otherwise
                # populate only NA cells.  This handles the layered case
                # where the column was created by an earlier sub-problem's
                # finalize_out and later rows were appended (segment splits,
                # post-anchor block) with NA for the new column.  Each cell
                # is still snapshotted exactly once.
                if (!("ReconstructedStemID_PreSweep" %in% names(out))) {
                    out[, ReconstructedStemID_PreSweep := ReconstructedStemID]
                } else {
                    out[is.na(ReconstructedStemID_PreSweep),
                        ReconstructedStemID_PreSweep := ReconstructedStemID]
                }
                # Compute audit flag against the FIRST-EVER engine output
                # (the PreSweep snapshot), not the current ReconstructedStemID.
                # Otherwise post-engine renumbering (R-boundary track
                # severance, segment merges) can re-introduce engine-vs-pin
                # disagreements that the sweep would re-correct, falsely
                # flagging rows whose original engine output already matched
                # TrueStemID.
                .pre_snap <- out$ReconstructedStemID_PreSweep
                .override <- .has_true_final &
                    !is.na(.pre_snap) &
                    .pre_snap != as.integer(out$TrueStemID)
                if (any(.override)) {
                    out[.override, SweepAuditOverride := TRUE]
                    .tag_id <- if ("Tag" %in% names(out) && length(out$Tag) > 0L) as.character(out$Tag[[1L]]) else "<unknown>"
                    vcat(prefix, "  [audit] finalize_out sweep overrode ", sum(.override), " engine-assigned ReconstructedStemID value(s) for tag ", .tag_id, " (rows flagged via SweepAuditOverride=TRUE)")
                }
                out[.has_true_final, `:=`(
                    ReconstructedStemID  = as.integer(TrueStemID),
                    ReconstructionMethod = "given"
                )]
            }
        }
        # Attach prune stats if available (some early returns may occur before prune_stats is initialized)
        attr(out, "DP_PruneInfo") <- if (exists("prune_stats")) prune_stats else list()
        # Row-count sanity check: output must have exactly as many rows as input
        .n_out <- nrow(out)
        if (.n_out != .n_input_rows) {
            warning(sprintf(
                "DP row-count mismatch for tag %s: input=%d rows, output=%d rows (delta=%+d). Check MF stash, post-anchor reinsertion, and segment split logic.",
                unique(original_tree_data$Tag)[1L], .n_input_rows, .n_out, .n_out - .n_input_rows
            ))
        }
        out
    }

    # Fallback helper: run probabilistic matcher and return a finalized output with the
    # given fallback reason.  When K_used is provided, also sets DP metadata
    # columns on tree_data (by reference) before running the matcher.
    do_fallback <- function(reason, K_used = NULL) {
        if (!is.null(K_used)) {
            K_used <- as.integer(K_used)
            if (K_used > 0L && length(obs_counts) > 0L) {
                n_states_by_census <- vapply(obs_counts, function(n_obs) count_injective_states(K_used, n_obs), numeric(1L))
                tree_data[, `:=`(
                    DP_KUsed = K_used,
                    DP_MaxStatesPerCensus = max(n_states_by_census, na.rm = TRUE),
                    DP_MaxStatesCensusID = as.integer(census_range[which.max(n_states_by_census)])
                )]
            } else {
                tree_data[, `:=`(
                    DP_KUsed = K_used,
                    DP_MaxStatesPerCensus = 0L,
                    DP_MaxStatesCensusID = NA_integer_
                )]
            }
        }
        # Route all fallbacks to probabilistic matcher
        out <- match_stems_probabilistic(
            tree_data, min_growth, max_growth, anchor_start,
            n_samples = prob_n_samples,
            temperature = temperature,
            posterior_top_k = posterior_top_k,
            posterior_samples_path = posterior_samples_path,
            posterior_samples_format = posterior_samples_format,
            posterior_sample_seed = posterior_sample_seed,
            prune_min_growth = prune_min_growth,
            prune_max_growth = prune_max_growth,
            prune_recruit_max_dbh = prune_recruit_max_dbh,
            prob_lookahead_weight = prob_lookahead_weight,
            use_bio_hard_shrink_in_prob = use_bio_hard_shrink_in_prob,
            use_bio_hard_growth_in_prob = use_bio_hard_growth_in_prob,
            n_sigma_me = prob_n_sigma_me,
            pin_truestemid = pin_truestemid,
            verbose = verbose
        )
        if (!("DP_FallbackReason" %in% names(out))) out[, DP_FallbackReason := NA_character_]
        out[, DP_FallbackReason := reason]

        # ---- NA-R barrier post-processing for fallback matchers ----
        # The DP has built-in NA-R barrier logic (post-Viterbi), but the
        # probabilistic matcher does not. Apply it here so both paths
        # produce identical ReconstructedStemID patterns for NA-DBH R-coded
        # rows. Must run BEFORE finalize_out() so MF reinsertion sees the
        # correct flanking IDs.
        .resprout_regex_fb <- "\\b(R|RP|RF|RT|QR|OR)\\b"
        .has_tsm_fb    <- "ListOfTSM" %in% names(out)
        .has_status_fb <- "Status" %in% names(out)
        {
            .pre_censuses_fb <- sort(unique(out$CensusID[out$CensusID <= anchor_start]))
            for (.cc_fb in .pre_censuses_fb) {
                .n_live_fb <- sum(out$CensusID == .cc_fb & !is.na(out$DBH))
                if (.n_live_fb > 0L) next
                .na_rows_fb <- which(out$CensusID == .cc_fb & is.na(out$DBH))
                if (length(.na_rows_fb) == 0L) next
                .na_tsm_fb <- if (.has_tsm_fb) out$ListOfTSM[.na_rows_fb] else rep(NA_character_, length(.na_rows_fb))
                .is_r_tsm_fb    <- !is.na(.na_tsm_fb) & grepl(.resprout_regex_fb, .na_tsm_fb, perl = TRUE)
                .is_r_status_fb <- if (.has_status_fb) !is.na(out$Status[.na_rows_fb]) & out$Status[.na_rows_fb] == "broken below" else rep(FALSE, length(.na_rows_fb))
                .is_r_fb <- .is_r_tsm_fb | .is_r_status_fb
                if (!any(.is_r_fb)) next
                # This census is an NA-R barrier: 0 live + NA-DBH R-coded rows
                .cens_before_fb <- .pre_censuses_fb[.pre_censuses_fb < .cc_fb]
                if (length(.cens_before_fb) == 0L) next
                # Split tracks that cross the barrier
                .ids_before_fb <- unique(out$ReconstructedStemID[out$CensusID %in% .cens_before_fb & !is.na(out$ReconstructedStemID)])
                .ids_after_fb <- unique(out$ReconstructedStemID[out$CensusID > .cc_fb & !is.na(out$ReconstructedStemID)])
                .crossing_fb <- intersect(.ids_before_fb, .ids_after_fb)
                .cur_max_fb <- suppressWarnings(max(out$ReconstructedStemID, na.rm = TRUE))
                if (!is.finite(.cur_max_fb)) .cur_max_fb <- 0L
                if (length(.crossing_fb) > 0L) {
                    for (.old_fb in .crossing_fb) {
                        .new_fb <- as.integer(.cur_max_fb) + 1L
                        .cur_max_fb <- .new_fb
                        out[
                            CensusID %in% .cens_before_fb & ReconstructedStemID == .old_fb,
                            ReconstructedStemID := .new_fb
                        ]
                    }
                }
                # Assign dying-stem NA-R rows to pre-barrier IDs
                .barrier_rows_fb <- .na_rows_fb[.is_r_fb]
                .prev_cc_fb <- max(.cens_before_fb)
                .prev_ids_fb <- out$ReconstructedStemID[out$CensusID == .prev_cc_fb & !is.na(out$ReconstructedStemID)]
                if (length(.barrier_rows_fb) > 0L && length(.prev_ids_fb) > 0L) {
                    for (.ri_fb in seq_along(.barrier_rows_fb)) {
                        .pid_fb <- if (.ri_fb <= length(.prev_ids_fb)) .prev_ids_fb[.ri_fb] else .prev_ids_fb[length(.prev_ids_fb)]
                        data.table::set(out, .barrier_rows_fb[.ri_fb], "ReconstructedStemID", as.integer(.pid_fb))
                    }
                    out[
                        .barrier_rows_fb[is.na(out$TrueStemID[.barrier_rows_fb])],
                        ReconstructionMethod := "probabilistic"
                    ]
                }
            }

            # ---- R-boundary splitting for LIVE R-coded censuses ----
            # The NA-R barrier above handles censuses with 0 live stems.
            # This block handles censuses where R-coded stems HAVE non-NA DBH:
            # the tree broke but the field team recorded DBH > 0 on the R row.
            # The R census marks an identity boundary: the R-coded stem IS the old
            # organism (last record), and any track that continues into censuses AFTER
            # the R census must be severed (post-R rows get new synthetic IDs).
            for (.cc_fb in .pre_censuses_fb) {
                .live_rows_fb <- which(out$CensusID == .cc_fb & !is.na(out$DBH))
                if (length(.live_rows_fb) == 0L) next
                .live_tsm_fb <- if (.has_tsm_fb) out$ListOfTSM[.live_rows_fb] else rep(NA_character_, length(.live_rows_fb))
                .any_r_live_tsm    <- any(!is.na(.live_tsm_fb) & grepl(.resprout_regex_fb, .live_tsm_fb, perl = TRUE))
                .any_r_live_status <- if (.has_status_fb) any(!is.na(out$Status[.live_rows_fb]) & out$Status[.live_rows_fb] == "broken below") else FALSE
                .any_r_live <- .any_r_live_tsm | .any_r_live_status
                if (!.any_r_live) next
                # This census has live R-coded stems → sever tracks continuing after R
                .cens_before_rfb <- .pre_censuses_fb[.pre_censuses_fb < .cc_fb]
                if (length(.cens_before_rfb) == 0L) next
                .ids_before_and_r_rfb <- unique(out$ReconstructedStemID[out$CensusID <= .cc_fb & !is.na(out$ReconstructedStemID)])
                .ids_after_rfb        <- unique(out$ReconstructedStemID[out$CensusID >  .cc_fb & !is.na(out$ReconstructedStemID)])
                .crossing_rfb <- intersect(.ids_before_and_r_rfb, .ids_after_rfb)
                .cur_max_fb <- suppressWarnings(max(out$ReconstructedStemID, na.rm = TRUE))
                if (!is.finite(.cur_max_fb)) .cur_max_fb <- 0L
                if (length(.crossing_rfb) > 0L) {
                    for (.old_rfb in .crossing_rfb) {
                        .new_rfb <- as.integer(.cur_max_fb) + 1L
                        .cur_max_fb <- .new_rfb
                        out[
                            CensusID > .cc_fb & ReconstructedStemID == .old_rfb,
                            ReconstructedStemID := .new_rfb
                        ]
                    }
                }
            }
        }

        finalize_out(out)
    }

    # Preserve original dataset in case we scope DP to only pre-anchor censuses
    original_tree_data <- data.table::copy(tree_data)
    anchor_requested <- anchor_start
    dp_scoped_to_pre_anchor <- FALSE

    # Determine last census with observed DBH for this tag
    obs_census_all <- sort(unique(original_tree_data$CensusID[!is.na(original_tree_data$DBH)]))
    last_obs_census <- if (length(obs_census_all) > 0L) as.integer(max(obs_census_all)) else NA_integer_

    # Rules:
    #  - If last observed census < requested anchor: use last observed census as anchor (anchor becomes terminal)
    #  - If last observed census > requested anchor: run DP only for censuses <= anchor (preserve post-anchor rows unchanged in output)
    if (!is.na(last_obs_census) && last_obs_census < anchor_start) {
        vcat(prefix, "Anchor C", anchor_start, " is beyond last observed census C", last_obs_census, "; adjusting anchor to C", last_obs_census)
        anchor_start <- last_obs_census
    } else if (!is.na(last_obs_census) && last_obs_census > anchor_requested) {
        vcat(prefix, "Post-anchor data detected (last census=C", last_obs_census, "); running DP on censuses up to C", anchor_start, " only, preserving later rows unchanged")
        dp_scoped_to_pre_anchor <- TRUE
        # restrict tree_data to pre-anchor censuses for DP computation
        tree_data <- original_tree_data[CensusID <= anchor_start]

        # If the nominal anchor census has no living stems (e.g. all rows have NA DBH due to a
        # resprout/die event), look forward for the first post-anchor census that has living stems.
        # Prefer the first such census that also carries known TrueStemIDs (giving the DP a real
        # anchor to work from).  Without this, tree_data is scoped to a dead census and the
        # downstream anchor-id check always fails → probabilistic fallback.
        .nominal_anchor_live <- original_tree_data[CensusID == anchor_start & !is.na(DBH)]
        if (nrow(.nominal_anchor_live) == 0L) {
            .post_live <- original_tree_data[CensusID > anchor_start & !is.na(DBH)]
            if (nrow(.post_live) > 0L) {
                # First post-anchor census with TrueStemID (preferred), else first with any DBH
                .post_with_id <- .post_live[!is.na(TrueStemID)]
                .extended_anchor <- if (nrow(.post_with_id) > 0L) {
                    as.integer(min(.post_with_id$CensusID))
                } else {
                    as.integer(min(.post_live$CensusID))
                }
                vcat(
                    prefix, "Nominal anchor C", anchor_start, " has no living stems; extending DP scope to C",
                    .extended_anchor, " (first post-anchor census with living stems)"
                )
                anchor_start <- .extended_anchor
                tree_data <- original_tree_data[CensusID <= anchor_start]
                # If the extended anchor now covers all observations, the DP is no
                # longer scoped to a prefix — disable post-anchor reinsertion so
                # the same rows are not appended a second time.
                if (nrow(original_tree_data[CensusID > anchor_start]) == 0L) {
                    dp_scoped_to_pre_anchor <- FALSE
                }
            }
        }

        # If scoping removes all observations (e.g., all observations are after the requested
        # anchor), return a one-row placeholder so the Tag is represented in outputs.
        if (nrow(tree_data) == 0L) {
            vcat(prefix, "No observations before anchor C", anchor_start, "; skipping DP and returning rows as-is")
            out <- data.table::copy(original_tree_data)
            # Use helper to normalize post-anchor rows: no DP performed, so treat observed TrueStemID as given
            out <- propagate_post_anchor_given(out)
            if (isTRUE(has_mf_stash) && !is.null(mf_stash) && nrow(mf_stash) > 0L) {
                out <- reinsert_mf_rows(out, mf_stash)
            }
            attr(out, "DP_PruneInfo") <- prune_stats
            return(out)
        }
    }

    # Preserve any provided TrueStemID as hard values in output —
    # but ONLY at the anchor census, where TrueStemID is actually used
    # as a constraint by the DP solver.  Pre-anchor TrueStemID is NOT
    # used to constrain identification; labelling it "given" would be
    # misleading.  Pre-anchor rows get their method from Viterbi ("dp").
    # EXCEPTION: provisional_dp rows keep their label — TrueStemID was
    # fabricated, not from field data.
    # NOTE: use %in% not != for NA-safe comparison (NA != "x" is NA in R,
    # which data.table treats as FALSE, silently dropping rows).
    tree_data[!is.na(TrueStemID) & CensusID == anchor_start &
              !ReconstructionMethod %in% "provisional_dp", `:=`(
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

    # Ensure a defensive `prune_stats` exists so any early-return branches can
    # safely attach it to outputs without failing when pruning hasn't run.
    prune_stats <- list(
        total_examined = 0L,
        total_pruned = 0L,
        per_census = integer(0)
    )

    if (is.na(first_obs_census)) {
        vcat(prefix, "No DBH observations found up to anchor C", anchor_start, "; routing to probabilistic matcher")
        return(do_fallback("no_obs_up_to_anchor", K_used = 0L))
    }
    census_range <- seq.int(from = first_obs_census, to = anchor_start)
    # Exclude censuses fully emptied by MF removal so the DP bridges them
    if (isTRUE(has_mf_stash) && length(mf_emptied_censuses) > 0L) {
        census_range <- census_range[!census_range %in% mf_emptied_censuses]
    }
    # Remove virtual census IDs with no rows in tree_data (e.g. gaps left by an
    # upstream MF-removal pass that the sub-call's own MF detection does not see).
    census_range <- census_range[census_range %in% unique(tree_data$CensusID)]
    n_census <- length(census_range)
    obs_counts <- vapply(
        census_range,
        function(cc) nrow(tree_data[CensusID == cc & !is.na(DBH)]),
        integer(1L)
    )
    max_obs <- if (length(obs_counts) > 0L) max(obs_counts) else 0L
    vcat(prefix, "Census range: C", paste(census_range, collapse = ", C"), " (", n_census, " censuses, stems per census: ", paste(obs_counts, collapse = "/"), ", anchor=C", anchor_start, ")")

    # --- growth-form fallback check ------------------------------------
    if (length(fallback_growth_forms) > 0L && "growth_form" %in% names(tree_data)) {
        bad_idx <- which(tree_data$growth_form %in% fallback_growth_forms)
        if (length(bad_idx) > 0L) {
            vcat(prefix, "Growth form '", paste(unique(tree_data$growth_form[bad_idx]), collapse = "', '"), "' requires probabilistic matcher; skipping DP")
            return(do_fallback("growth_form_forced", K_used = as.integer(min(max_obs, max_tracks))))
        }
    }

    # --- species-based probabilistic routing check ----------------------
    if (length(prob_species) > 0L && "Species" %in% names(tree_data)) {
        sp_vals <- unique(tree_data$Species)
        sp_hit <- sp_vals[sp_vals %in% prob_species]
        if (length(sp_hit) > 0L) {
            vcat(prefix, "Species '", paste(sp_hit, collapse = "', '"), "' is in prob_species list; routing to probabilistic matcher")
            return(do_fallback("species_forced_probabilistic", K_used = as.integer(min(max_obs, max_tracks))))
        }
    }

    # Defensive initialization for prune diagnostics so early-return branches can safely set attributes
    prune_stats <- list(
        total_examined = 0L,
        total_pruned = 0L,
        per_census = integer(0)
    )

    # Initialize conservative hard-pruning diagnostics (per transition between adjacent censuses)
    prune_hard <- isTRUE(prune_hard)

    # Normalize pruning parameters
    if (!is.null(prune_min_growth)) prune_min_growth <- as.numeric(prune_min_growth)
    if (!is.null(prune_max_growth)) prune_max_growth <- as.numeric(prune_max_growth)
    prune_use_bio_bounds <- isTRUE(prune_use_bio_bounds)
    if (!is.null(prune_recruit_max_dbh)) prune_recruit_max_dbh <- as.numeric(prune_recruit_max_dbh)
    prune_use_bio_recruit <- isTRUE(prune_use_bio_recruit)
    non_taper_corrected_prune_min_growth <- as.numeric(non_taper_corrected_prune_min_growth)
    non_taper_corrected_prune_max_growth <- as.numeric(non_taper_corrected_prune_max_growth)
    hom_tolerance_scale <- as.numeric(hom_tolerance_scale)

    # --- Non-taper-corrected growth-form detection ---
    # Interpret comma/semicolon-separated string as a vector (like fallback_growth_forms).
    if (length(non_taper_corrected_growth_forms) == 1L && grepl("[,;]", non_taper_corrected_growth_forms)) {
        ff <- strsplit(non_taper_corrected_growth_forms, "[,;]")[[1L]]
        ff <- trimws(ff)
        ff <- ff[nzchar(ff)]
        non_taper_corrected_growth_forms <- ff
    }
    non_taper_corrected_growth_forms <- as.character(non_taper_corrected_growth_forms)

    is_non_taper_corrected <- FALSE
    if ("growth_form" %in% names(tree_data)) {
        gf_vals <- unique(tree_data$growth_form)
        if (length(gf_vals) > 1L) {
            warning(prefix, "Multiple growth_form values found; using most common to determine non-taper-corrected status.")
            gf_vals <- names(sort(table(tree_data$growth_form), decreasing = TRUE))[1L]
        }
        is_non_taper_corrected <- isTRUE(gf_vals %in% non_taper_corrected_growth_forms)
    }

    # --- HOM column detection (for HOM-proportional widening) ---
    use_hom_relax <- FALSE
    hom_col <- NULL
    if (is_non_taper_corrected && hom_tolerance_scale > 0) {
        hom_candidates <- intersect(tolower(names(tree_data)), c("hom"))
        if (length(hom_candidates) > 0L) {
            # Find the actual column name matching case-insensitively
            hom_col <- names(tree_data)[tolower(names(tree_data)) == "hom"][1L]
            use_hom_relax <- TRUE
            vcat(prefix, "HOM column detected ('", hom_col, "'); HOM-proportional widening enabled (scale=", hom_tolerance_scale, ").")
        } else {
            vcat(prefix, "Non-taper-corrected growth form but no HOM column found; HOM widening disabled.")
        }
    }

    prune_stats <- list(
        total_examined = 0L,
        total_pruned = 0L,
        per_census = setNames(integer(max(0, n_census - 1L)), as.character(if (n_census > 1L) census_range[-length(census_range)] else integer(0)))
    )

    # Need a fully-anchored endpoint
    # If the requested anchor census is missing entirely or all rows at that census
    # have both NA DBH and NA TrueStemID, prefer the most recent earlier census
    # that has at least one row with non-NA DBH and non-NA TrueStemID and use that
    # as the anchor instead of immediately falling back to the probabilistic matcher.
    anchor_rows_all <- tree_data[CensusID == anchor_start]
    if (nrow(anchor_rows_all) == 0L || (all(is.na(anchor_rows_all$DBH)) && all(is.na(anchor_rows_all$TrueStemID)))) {
        cand_census <- sort(unique(tree_data$CensusID[tree_data$CensusID < anchor_start & !is.na(tree_data$DBH) & !is.na(tree_data$TrueStemID)]))
        if (length(cand_census) > 0L) {
            new_anchor <- as.integer(max(cand_census))
            vcat(prefix, "Anchor C", anchor_start, " has no DBH/TrueStemID; falling back to earlier anchor at C", new_anchor)
            anchor_start <- new_anchor
            census_range <- seq.int(from = first_obs_census, to = anchor_start)
            # Exclude censuses fully emptied by MF removal
            if (isTRUE(has_mf_stash) && length(mf_emptied_censuses) > 0L) {
                census_range <- census_range[!census_range %in% mf_emptied_censuses]
            }
            census_range <- census_range[census_range %in% unique(tree_data$CensusID)]
            n_census <- length(census_range)
            vcat(prefix, "Adjusted census range: C", paste(census_range, collapse = ", C"))
            obs_counts <- vapply(
                census_range,
                function(cc) nrow(tree_data[CensusID == cc & !is.na(DBH)]),
                integer(1L)
            )
            max_obs <- if (length(obs_counts) > 0L) max(obs_counts) else 0L
        } else {
            # If configured, allow a provisional DP anchor at the last observed DBH census
            if (isTRUE(allow_provisional_anchor) && !is.na(last_obs_census) && any(!is.na(tree_data$DBH[tree_data$CensusID == last_obs_census]))) {
                provisional_anchor <- as.integer(last_obs_census)
                vcat(prefix, "No census with DBH+TrueStemID found; using provisional anchor at C", provisional_anchor)
                anchor_start <- provisional_anchor

                # Assign provisional TrueStemID/ReconstructedStemID at the anchor rows
                anchor_idx <- which(tree_data$CensusID == anchor_start & !is.na(tree_data$DBH))
                if (length(anchor_idx) > 0L) {
                    current_max <- suppressWarnings(max(tree_data$TrueStemID, na.rm = TRUE))
                    if (!is.finite(current_max)) current_max <- 0L
                    prov_ids <- as.integer(seq.int(from = current_max + 1L, length.out = length(anchor_idx)))
                    tree_data$TrueStemID[anchor_idx] <- prov_ids
                    tree_data$ReconstructedStemID[anchor_idx] <- prov_ids
                    tree_data$ReconstructionMethod[anchor_idx] <- "provisional_dp"
                    tree_data$ConstraintViolation[anchor_idx] <- FALSE
                    vcat(prefix, sprintf("Assigned %d provisional anchor ID(s) at C%d", length(anchor_idx), anchor_start))
                }

                # Recompute census_range and obs_counts now that anchor_start changed
                census_range <- seq.int(from = first_obs_census, to = anchor_start)
                # Exclude censuses fully emptied by MF removal
                if (isTRUE(has_mf_stash) && length(mf_emptied_censuses) > 0L) {
                    census_range <- census_range[!census_range %in% mf_emptied_censuses]
                }
                census_range <- census_range[census_range %in% unique(tree_data$CensusID)]
                n_census <- length(census_range)
                vcat(prefix, "Adjusted census range: C", paste(census_range, collapse = ", C"))
                obs_counts <- vapply(
                    census_range,
                    function(cc) nrow(tree_data[CensusID == cc & !is.na(DBH)]),
                    integer(1L)
                )
                max_obs <- if (length(obs_counts) > 0L) max(obs_counts) else 0L
            } else {
                vcat(prefix, "No usable anchor (no DBH or TrueStemID at any census); routing to probabilistic matcher")
                return(do_fallback("anchor_missing_truestem", K_used = as.integer(min(max_obs, max_tracks))))
            }
        }
    }

    anchor_obs <- tree_data[CensusID == anchor_start & !is.na(DBH)]
    if (nrow(anchor_obs) == 0L) {
        vcat(prefix, "Anchor C", anchor_start, " has no DBH observations; routing to probabilistic matcher")
        return(do_fallback("anchor_missing_obs", K_used = as.integer(min(max_obs, max_tracks))))
    }

    # If anchor_obs exists but some TrueStemID are missing, optionally allow a provisional DP anchor
    if (any(is.na(anchor_obs$TrueStemID))) {
        if (isTRUE(allow_provisional_anchor)) {
            vcat(prefix, "Anchor C", anchor_start, " has DBH but missing TrueStemID; assigning provisional IDs")
            anchor_idx <- which(tree_data$CensusID == anchor_start & !is.na(tree_data$DBH))
            current_max <- suppressWarnings(max(tree_data$TrueStemID, na.rm = TRUE))
            if (!is.finite(current_max)) current_max <- 0L
            prov_ids <- as.integer(seq.int(from = current_max + 1L, length.out = length(anchor_idx)))
            tree_data$TrueStemID[anchor_idx] <- prov_ids
            tree_data$ReconstructedStemID[anchor_idx] <- prov_ids
            tree_data$ReconstructionMethod[anchor_idx] <- "provisional_dp"
            tree_data$ConstraintViolation[anchor_idx] <- FALSE
            # Recompute anchor_obs and anchor_ids after provisioning
            anchor_obs <- tree_data[CensusID == anchor_start & !is.na(DBH)]
            anchor_ids <- sort(unique(anchor_obs$TrueStemID))
            anchor_ids <- anchor_ids[!is.na(anchor_ids)]
        } else {
            vcat(prefix, "No usable anchor (TrueStemID missing and provisional anchoring disabled); routing to probabilistic matcher")
            return(do_fallback("anchor_missing_truestem_prov_disabled", K_used = as.integer(min(max_obs, max_tracks))))
        }
    }
    anchor_ids <- sort(unique(anchor_obs$TrueStemID))
    anchor_ids <- anchor_ids[!is.na(anchor_ids)]

    if (length(anchor_ids) == 0L) {
        vcat(prefix, "No valid anchor IDs at anchor C", anchor_start, "; routing to probabilistic matcher")
        return(do_fallback("anchor_ids_missing", K_used = as.integer(min(max_obs, max_tracks))))
    }

    # Choose K (tracks) same logic as the MAP DP
    births_needed <- if (length(obs_counts) >= 2L) sum(pmax(0L, diff(obs_counts))) else 0L
    K_from_counts <- as.integer(if (length(obs_counts) > 0L) obs_counts[1L] + births_needed else 0L)
    K_base <- max(length(anchor_ids), max_obs, K_from_counts)

    # ---- Resprout barrier: increase K for resprout observations ----
    # Each resprout (R|RP|RF|RT|QR|OR code with non-NA DBH) forces the track into
    # phase 0 at the preceding census, effectively creating a new identity that
    # needs its own track slot.  Add one extra track per resprout observation.
    resprout_regex <- "\\b(R|RP|RF|RT|QR|OR)\\b"
    # Pre-compute per-census row indices and resprout flags ONCE here.
    # This single pass is reused in both the resprout count below and the
    # state enumeration loop after K is determined, eliminating duplicate
    # [.data.table subset calls from those two loops.
    .has_tsm <- "ListOfTSM" %in% names(tree_data)
    .obs_row_idx_pre <- vector("list", n_census)
    .is_resprout_pre <- vector("list", n_census)
    .has_na_r_barrier <- logical(n_census) # TRUE: census has 0 live stems AND NA-DBH R-coded rows
    for (.p0 in seq_len(n_census)) {
        .idx0 <- tree_data[CensusID == census_range[.p0] & !is.na(DBH), which = TRUE]
        .obs_row_idx_pre[[.p0]] <- .idx0
        if (length(.idx0) > 0L) {
            .tsm_flag0 <- if (.has_tsm) {
                .tsm0 <- tree_data$ListOfTSM[.idx0]
                !is.na(.tsm0) & grepl(resprout_regex, .tsm0, perl = TRUE)
            } else rep(FALSE, length(.idx0))
            .status_flag0 <- if ("Status" %in% names(tree_data)) {
                !is.na(tree_data$Status[.idx0]) & tree_data$Status[.idx0] == "broken below"
            } else rep(FALSE, length(.idx0))
            .is_resprout_pre[[.p0]] <- .tsm_flag0 | .status_flag0
        } else {
            .is_resprout_pre[[.p0]] <- logical(0L)
            # Detect NA-R barriers: 0 live stems AND at least one NA-DBH R-coded row
            # (these mark hard resprout boundaries where track identities must be split)
            .na_rows <- tree_data[CensusID == census_range[.p0] & is.na(DBH), which = TRUE]
            if (length(.na_rows) > 0L) {
                .is_bb_tsm0 <- if (.has_tsm) {
                    .na_tsm0 <- tree_data$ListOfTSM[.na_rows]
                    !is.na(.na_tsm0) & grepl(resprout_regex, .na_tsm0, perl = TRUE)
                } else rep(FALSE, length(.na_rows))
                .is_bb_status0 <- if ("Status" %in% names(tree_data)) {
                    !is.na(tree_data$Status[.na_rows]) & tree_data$Status[.na_rows] == "broken below"
                } else rep(FALSE, length(.na_rows))
                .has_na_r_barrier[.p0] <- any(.is_bb_tsm0 | .is_bb_status0)
            }
        }
    }
    # -----------------------------------------------------------------------
    # Resprout segment split
    # When any census p0 >= 2 (not the first in range) carries an R code and
    # lies before the anchor, split the DP into two independent sub-problems:
    #   pre-segment  : censuses 1 .. (r_boundary - 1), provisional anchor
    #   post-segment : censuses r_boundary .. anchor_start
    # Biological rationale: R means the tree resprouted — stems at the R
    # census are entirely new physical entities with no identity continuity
    # to stems at earlier censuses.  Splitting prevents pre-resprout census
    # history from contaminating post-resprout track assignments.
    # The post sub-call starts at position 1 in its range (does not satisfy
    # .p0 >= 2), so it never triggers a further split.  Additionally,
    # allow_segment_split=FALSE in recursive sub-calls prevents cascading
    # splits; downstream R codes are handled by R-recruit constraints.
    # -----------------------------------------------------------------------
    .r_boundary_pos <- NULL
    if (isTRUE(allow_segment_split)) {
        for (.chk_p in seq_len(n_census)) {
            if (.chk_p >= 2L &&
                any(.is_resprout_pre[[.chk_p]]) &&
                census_range[.chk_p] < anchor_start) {
                .r_boundary_pos <- .chk_p
                break
            }
        }
    }
    if (!is.null(.r_boundary_pos) && isTRUE(allow_segment_split)) {
        .r_boundary_census <- census_range[.r_boundary_pos]
        vcat(
            prefix, "--- RESPROUT SPLIT --- C", .r_boundary_census,
            " has resprout code; splitting into pre-resprout (C", census_range[1L],
            "-C", census_range[.r_boundary_pos - 1L], ") and post-resprout (C",
            .r_boundary_census, "-C", anchor_start, ") segments"
        )

        # Shared sub-call arguments (all params in scope, some coerced earlier)
        .sub_args <- list(
            min_growth = min_growth,
            max_growth = max_growth,
            max_tracks = max_tracks,
            max_states = max_states,
            temperature = temperature,
            posterior_top_k = posterior_top_k,
            eps_tiebreak = eps_tiebreak,
            allow_provisional_anchor = allow_provisional_anchor,
            use_measurement_error = use_measurement_error,
            meas_sd1_a = meas_sd1_a,
            meas_sd1_b = meas_sd1_b,
            meas_sd2 = meas_sd2,
            meas_p_big = meas_p_big,
            fallback_growth_forms = fallback_growth_forms,
            posterior_samples = 0L, # disable DP posteriors in sub-calls
            posterior_samples_format = posterior_samples_format,
            posterior_samples_path = posterior_samples_path,
            posterior_sample_seed = posterior_sample_seed,
            prune_hard = prune_hard,
            prune_min_growth = prune_min_growth,
            prune_max_growth = prune_max_growth,
            prune_use_bio_bounds = prune_use_bio_bounds,
            prune_recruit_max_dbh = prune_recruit_max_dbh,
            prune_use_bio_recruit = prune_use_bio_recruit,
            non_taper_corrected_growth_forms = non_taper_corrected_growth_forms,
            non_taper_corrected_prune_min_growth = non_taper_corrected_prune_min_growth,
            non_taper_corrected_prune_max_growth = non_taper_corrected_prune_max_growth,
            hom_tolerance_scale = hom_tolerance_scale,
            verbose = verbose,
            chunk_id = chunk_id,
            allow_segment_split = FALSE, # prevent cascading splits in sub-calls
            prob_n_samples = prob_n_samples,
            prob_species = prob_species,
            prob_lookahead_weight = prob_lookahead_weight,
            prob_n_sigma_me = prob_n_sigma_me,
            pin_truestemid = pin_truestemid
        )

        # When the whole-tag state space exceeds max_states, force probabilistic
        # matching on both sub-segments for consistency.  Without this, one
        # segment might use probabilistic while the other uses DP (because K
        # is smaller in the sub-segment).
        .whole_tag_max_states <- max(vapply(obs_counts, function(n_obs) {
            count_injective_states(K_base, n_obs)
        }, numeric(1)))
        if (.whole_tag_max_states > max_states) {
            vcat(
                prefix, "Whole-tag state space (",
                format(.whole_tag_max_states, big.mark = ","),
                ") exceeds max_states=", format(as.numeric(max_states), big.mark = ","),
                "; forcing probabilistic on both segments"
            )
            .sub_args$max_states <- 0L
        }

        # Post-resprout sub-call: censuses > r_boundary with original anchor
        # R census itself belongs to the pre-segment (old stem's last record).
        # Use tree_data (already scoped to <= anchor_start) so post-anchor rows
        # are only re-added once by the outer finalize_out.
        .post_data <- tree_data[CensusID > .r_boundary_census]
        out_post <- do.call(
            match_stems_dp_global_backward_marginals_batch,
            c(
                list(
                    tree_data = .post_data,
                    anchor_start = anchor_start,
                    slack_tracks = slack_tracks,
                    slack_require_anchor_recruitable = slack_require_anchor_recruitable,
                    slack_require_anchor_eps = slack_require_anchor_eps,
                    post_segment_all_recruits = TRUE
                ),
                .sub_args
            )
        )

        # Pre-resprout sub-call: censuses <= r_boundary with provisional anchor
        # Includes the R census (old stem's last record) as the terminal observation.
        .pre_data <- tree_data[CensusID <= .r_boundary_census]
        .pre_anchor <- suppressWarnings(
            max(.pre_data$CensusID[!is.na(.pre_data$DBH)], na.rm = TRUE)
        )
        if (is.finite(.pre_anchor)) {
            out_pre <- do.call(
                match_stems_dp_global_backward_marginals_batch,
                c(
                    list(
                        tree_data = .pre_data,
                        anchor_start = as.integer(.pre_anchor),
                        slack_tracks = slack_tracks,
                        # No real anchor in pre-segment: disable recruitable check
                        slack_require_anchor_recruitable = FALSE,
                        slack_require_anchor_eps = slack_require_anchor_eps
                    ),
                    .sub_args
                )
            )
            # Offset pre-segment IDs so they do not clash with post-segment IDs
            .max_post_id <- suppressWarnings(
                max(out_post$ReconstructedStemID, na.rm = TRUE)
            )
            if (!is.finite(.max_post_id)) .max_post_id <- 0L
            .offset <- as.integer(.max_post_id)
            if (.offset > 0L) {
                out_pre[
                    !is.na(ReconstructedStemID),
                    ReconstructedStemID := ReconstructedStemID + .offset
                ]
                # Offset posterior top-k IDs as well
                for (.k in seq_len(posterior_top_k)) {
                    .id_col <- paste0("DP_PosteriorTop", .k, "ID")
                    if (.id_col %in% names(out_pre)) {
                        out_pre[
                            !is.na(get(.id_col)),
                            (.id_col) := get(.id_col) + .offset
                        ]
                    }
                }
            }
            # Combine pre + post segments
            .all_cols <- union(names(out_pre), names(out_post))
            for (.c in setdiff(.all_cols, names(out_pre))) out_pre[, (.c) := NA]
            for (.c in setdiff(.all_cols, names(out_post))) out_post[, (.c) := NA]
            combined <- data.table::rbindlist(
                list(out_pre, out_post),
                use.names = TRUE, fill = TRUE
            )
        } else {
            combined <- out_post
        }
        if ("obs_row_id" %in% names(combined)) data.table::setorder(combined, obs_row_id)
        return(finalize_out(combined))
    }

    # (R-coded observations continue existing tracks — no extra track slots needed.)

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
            vcat(prefix, "Slack tracks denied: no anchor stem has DBH <= recruit max + eps (eps=", eps, ")")
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
        vcat(prefix, "K=", K, " tracks < max observed stems (", max(obs_counts), " at C", census_range[which.max(obs_counts)], "); DP infeasible, routing to probabilistic matcher")
        return(do_fallback("K_too_small"))
    }

    # ---- Hard-invariant track-set extension --------------------------------
    # Any pre-anchor row may carry a non-NA TrueStemID (e.g. via Step 2
    # terminal propagation in main_cpp_bci.R) whose value is NOT present at
    # the anchor census (typical case: a stem that died before the BCI
    # re-tagging campaign at C7).  Without including these IDs in track_ids,
    # the pin map below builds match(tsid, track_ids) -> NA and the row is
    # silently left unpinned, allowing the DP to violate the invariant
    # ReconstructedStemID == TrueStemID.  We enumerate them once here and
    # add them to the slack track set.
    .pre_anchor_tsids <- sort(unique(
        tree_data$TrueStemID[
            !is.na(tree_data$TrueStemID) &
            tree_data$CensusID < anchor_start &
            !(tree_data$TrueStemID %in% anchor_ids)
        ]
    ))
    # Make sure K is large enough to host every anchor + pre-anchor pin.
    .required_K <- length(anchor_ids) + length(.pre_anchor_tsids)
    if (.required_K > K) {
        K <- .required_K
    }

    vcat(prefix, "Using K=", K, " identity tracks (", length(anchor_ids), " anchored + ", length(.pre_anchor_tsids), " pre-anchor pinned + ", K - length(anchor_ids) - length(.pre_anchor_tsids), " extra); worst-case states/census: ", format(max(n_states_by_census, na.rm = TRUE), big.mark = ","), " at C", census_range[which.max(n_states_by_census)])

    n_extra <- K - length(anchor_ids) - length(.pre_anchor_tsids)
    current_max <- suppressWarnings(max(c(tree_data$TrueStemID, anchor_ids, .pre_anchor_tsids), na.rm = TRUE))
    if (!is.finite(current_max)) current_max <- 0
    track_ids <- c(
        anchor_ids,
        as.integer(.pre_anchor_tsids),
        if (n_extra > 0L) seq.int(from = current_max + 1L, length.out = n_extra) else integer(0)
    )

    # -----------------------------------------------------------------------
    # Constrained enumeration: backward per-interval arc consistency.
    # For each observation at census p, restrict which tracks are feasible
    # by checking whether a growth-compatible observation exists at census
    # p+1 that can also use that track (as determined at the previous
    # backward step).  Constraints compound across intervals because each
    # census inherits the tighter filtering from its successor.
    #
    # Anchor tracks vs slack tracks:
    #   Anchor tracks have a known identity at the anchor census.
    #   Slack tracks are empty at the anchor (created for recruits/deaths)
    #   and are always allowed for every observation.
    # -----------------------------------------------------------------------
    .anchor_track_set <- which(track_ids %in% anchor_ids) # tracks with anchor obs
    .slack_track_set <- which(!(track_ids %in% anchor_ids)) # tracks empty at anchor

    # Median date per census for per-interval dt computation
    .has_exact_date <- "ExactDate" %in% names(tree_data)
    .census_median_date <- setNames(rep(NA_real_, n_census), as.character(census_range))
    if (.has_exact_date) {
        for (.ci in seq_len(n_census)) {
            .dates <- tree_data[CensusID == census_range[.ci] & !is.na(DBH), ExactDate]
            if (length(.dates) > 0L) {
                .census_median_date[.ci] <- as.numeric(median(.dates, na.rm = TRUE))
            }
        }
    }

    # Conservative growth bounds for constrained enumeration.
    # Must be at least as wide as actual effective prune bounds (computed later)
    # to avoid falsely excluding feasible track assignments.
    # For non-taper-corrected species: multiply by 2 to account for possible
    # HOM-proportional widening that is applied per interval during the
    # backward pass.
    if (isTRUE(is_non_taper_corrected)) {
        .ce_min_grow <- 2.0 * non_taper_corrected_prune_min_growth
        .ce_max_grow <- 2.0 * non_taper_corrected_prune_max_growth
    } else {
        .ce_min_grow <- if (!is.null(prune_min_growth)) prune_min_growth else min_growth
        .ce_max_grow <- if (!is.null(prune_max_growth)) prune_max_growth else max_growth
    }

    # --- Pass 1 (forward): gather observation data per census ---------------
    obs_dbh <- vector("list", n_census)
    obs_row_idx <- vector("list", n_census)
    is_resprout_obs <- vector("list", n_census)
    for (p in seq_len(n_census)) {
        idx <- .obs_row_idx_pre[[p]]
        obs_row_idx[[p]] <- idx
        obs_dbh[[p]] <- tree_data$DBH[idx]
        is_resprout_obs[[p]] <- .is_resprout_pre[[p]]
    }

    # --- Pre-compute TrueStemID pin map for non-anchor censuses -------------
    # pin_tidx_at_census[[p]][j] = track index for obs j if pinned, else NA
    pin_tidx_at_census <- vector("list", n_census)
    if (isTRUE(pin_truestemid)) {
        for (p in seq_len(n_census)) {
            if (p == n_census) next  # anchor census handled separately
            idx <- obs_row_idx[[p]]
            if (length(idx) == 0L) next
            tsid <- tree_data$TrueStemID[idx]
            tidx <- match(tsid, track_ids)
            # Only pin where TrueStemID is non-NA and found in track_ids
            tidx[is.na(tsid)] <- NA_integer_
            # Duplicate-pin guard: if two obs claim the same track, keep first only
            .seen_tracks <- integer(0)
            for (.j in seq_along(tidx)) {
                if (is.na(tidx[.j])) next
                if (tidx[.j] %in% .seen_tracks) {
                    vcat(prefix, "WARNING: duplicate TrueStemID pin at C",
                         census_range[p], " for track ", track_ids[tidx[.j]],
                         "; keeping first, releasing obs ", .j)
                    tidx[.j] <- NA_integer_
                } else {
                    .seen_tracks <- c(.seen_tracks, tidx[.j])
                }
            }
            pin_tidx_at_census[[p]] <- tidx
        }
    }

    # --- Pass 2 (backward): enumerate states with per-interval constraints --
    state_mats <- vector("list", n_census)
    state_keys <- vector("list", n_census)
    .allowed_at_census <- vector("list", n_census) # per-obs allowed tracks

    for (p in seq.int(n_census, 1L, by = -1L)) {
        cc <- census_range[p]
        n_obs <- length(obs_dbh[[p]])

        # --- Determine allowed tracks for each observation ------------------
        .use_constrained <- FALSE

        if (p == n_census) {
            # Anchor census: each obs is PINNED to its specific track via
            # TrueStemID.  Setting allowed to exactly that one track (rather
            # than all K) gives the backward propagation a tight starting
            # point — the first hop checks "can obs_i grow to the exact
            # anchor DBH for track k?", and subsequent hops compound.
            if (n_obs > 0L) {
                .anchor_obs_pre <- tree_data[CensusID == anchor_start & !is.na(DBH)]
                .anchor_tidx <- match(.anchor_obs_pre$TrueStemID, track_ids)
                .allowed_at_census[[p]] <- lapply(seq_len(n_obs), function(.j) {
                    if (!is.na(.anchor_tidx[.j])) .anchor_tidx[.j] else seq_len(K)
                })
            }
        } else if (n_obs > 0L && .has_exact_date &&
            is.finite(.census_median_date[p]) &&
            is.finite(.census_median_date[p + 1L]) &&
            length(obs_dbh[[p + 1L]]) > 0L &&
            length(.anchor_track_set) > 0L) {
            # Per-interval backward propagation: for each obs i at census p,
            # track k is feasible iff there exists some obs j at census p+1
            # with k in allowed_at_census[[p+1]][[j]] AND the growth rate
            # (DBH_j - DBH_i) / dt is within bounds.
            .dt <- (.census_median_date[p + 1L] - .census_median_date[p]) / 365.25
            if (is.finite(.dt) && .dt > 0) {
                .dbh_next <- obs_dbh[[p + 1L]]
                .n_obs_next <- length(.dbh_next)
                .allowed_next <- .allowed_at_census[[p + 1L]] # already computed
                .allowed <- vector("list", n_obs)
                for (.oi in seq_len(n_obs)) {
                    .dbh_i <- obs_dbh[[p]][.oi]
                    .ok_tracks <- integer(0)
                    # Check each anchor track
                    for (.tk in .anchor_track_set) {
                        .feasible <- FALSE
                        for (.oj in seq_len(.n_obs_next)) {
                            if (!(.tk %in% .allowed_next[[.oj]])) next
                            .rate <- (.dbh_next[.oj] - .dbh_i) / .dt
                            if (.rate >= .ce_min_grow && .rate <= .ce_max_grow) {
                                .feasible <- TRUE
                                break
                            }
                        }
                        if (.feasible) .ok_tracks <- c(.ok_tracks, .tk)
                    }
                    # Slack tracks are always allowed
                    .ok_tracks <- c(.ok_tracks, .slack_track_set)
                    .allowed[[.oi]] <- sort(.ok_tracks)
                }
                .allowed_at_census[[p]] <- .allowed
                if (!all(lengths(.allowed) == K)) {
                    .use_constrained <- TRUE
                }
            } else {
                # dt not valid: all tracks allowed
                if (n_obs > 0L) {
                    .allowed_at_census[[p]] <- replicate(n_obs, seq_len(K), simplify = FALSE)
                }
            }
        } else {
            # No ExactDate or no obs at next census: all tracks allowed
            if (n_obs > 0L) {
                .allowed_at_census[[p]] <- replicate(n_obs, seq_len(K), simplify = FALSE)
            }
        }

        # --- TrueStemID pinning override (non-anchor censuses) --------------
        if (p < n_census && isTRUE(pin_truestemid) && !is.null(pin_tidx_at_census[[p]])) {
            .pins <- pin_tidx_at_census[[p]]
            .n_pinned <- sum(!is.na(.pins))
            if (.n_pinned > 0L) {
                for (.oi in which(!is.na(.pins))) {
                    .allowed_at_census[[p]][[.oi]] <- .pins[.oi]
                }
                .use_constrained <- TRUE
                vcat(prefix, "  Pinned ", .n_pinned, " obs at C", cc, " via TrueStemID")
            }
        }

        # --- Enumerate states -----------------------------------------------
        if (.use_constrained) {
            vcat(
                prefix, "  State enum C", cc, " (constrained): allowed tracks per stem = [",
                paste(lengths(.allowed_at_census[[p]]), collapse = ","), "] of K=", K
            )
            mat <- enumerate_states_constrained(K, n_obs, .allowed_at_census[[p]], max_states)
        } else {
            mat <- enumerate_states_injective(K, n_obs, max_states = max_states)
        }

        if (is.null(mat)) {
            vcat(prefix, "State space overflow at C", cc, ": ", n_obs, " stems would exceed max_states=", format(max_states, big.mark = ","), "; routing to probabilistic matcher")
            return(do_fallback("enum_exceeded"))
        }
        state_mats[[p]] <- mat
        state_keys[[p]] <- if (ncol(mat) == 0L) {
            rep("", nrow(mat))
        } else {
            do.call(paste, c(lapply(seq_len(ncol(mat)), function(j) mat[, j]), list(sep = ",")))
        }

        vcat(prefix, "  State enum C", cc, ": ", n_obs, " stem(s) -> ", format(nrow(mat), big.mark = ","), " assignment states",
            if (n_obs > 0L && any(!is.na(obs_dbh[[p]]))) paste0(" (DBH ", round(min(obs_dbh[[p]], na.rm = TRUE), 1), "-", round(max(obs_dbh[[p]], na.rm = TRUE), 1), " cm)") else "")
    }

    # Anchor state assignment (pins endpoint)
    anchor_obs_ordered <- tree_data[CensusID == anchor_start & !is.na(DBH)]
    anchor_track_idx <- match(anchor_obs_ordered$TrueStemID, track_ids)

    if (any(is.na(anchor_track_idx))) {
        vcat(prefix, "Anchor IDs not found in track list (TrueStemID mismatch); routing to probabilistic matcher")
        return(do_fallback("anchor_truestem_not_found"))
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

    # Fail-fast checks for required biological params (must be present and scalar)
    required_bio <- c(
        "Bio_Mu_Growth",
        "Bio_Sigma0_Growth",
        "Bio_Sigma1_Growth",
        "Bio_H0_Mortality",
        "Bio_Beta_Mortality",
        "Bio_Recruit_Meanlog",
        "Bio_Recruit_Sdlog",
        "Bio_Recruit_MaxDBH_unit",
        "Bio_Recruitment_lambda",
        "Bio_Max_Shrink",
        "Bio_K_Shrink"
    )
    missing_bio <- setdiff(required_bio, names(tree_data))
    bad_bio <- character(0)
    for (nm in intersect(required_bio, names(tree_data))) {
        u <- tryCatch(unique(tree_data[[nm]]), error = function(e) NA)
        if (length(u) != 1L || any(is.na(u))) bad_bio <- c(bad_bio, nm)
    }
    if (length(missing_bio) > 0L || length(bad_bio) > 0L) {
        msg_parts <- character(0)
        if (length(missing_bio) > 0L) msg_parts <- c(msg_parts, paste0("missing required bio columns: ", paste(missing_bio, collapse = ", ")))
        if (length(bad_bio) > 0L) msg_parts <- c(msg_parts, paste0("non-scalar or NA bio columns: ", paste(bad_bio, collapse = ", ")))
        stop(prefix, "Required biological parameters missing or invalid: ", paste(msg_parts, collapse = "; "))
    }

    # --- Determine effective pruning bounds (separate from biological min/max)
    user_min <- if (!is.null(prune_min_growth)) prune_min_growth else min_growth
    user_max <- if (!is.null(prune_max_growth)) prune_max_growth else max_growth
    if (isTRUE(prune_use_bio_bounds)) {
        eff_min_grow <- max(user_min, Bio_max_shrink_unit)
        eff_max_grow <- min(user_max, Bio_max_growth_unit)
    } else {
        eff_min_grow <- user_min
        eff_max_grow <- user_max
    }

    if (!is.null(prune_recruit_max_dbh)) {
        if (isTRUE(prune_use_bio_recruit) && is.finite(Bio_Recruit_MaxDBH_unit)) {
            eff_recruit_max <- min(prune_recruit_max_dbh, Bio_Recruit_MaxDBH_unit)
        } else {
            eff_recruit_max <- prune_recruit_max_dbh
        }
    } else {
        eff_recruit_max <- Bio_Recruit_MaxDBH_unit
    }

    # Non-taper-corrected override: replace effective bounds with wide values.
    # Palms, strangler figs, and tree ferns show real DBH growth AND large
    # apparent DBH variation when HOM changes between censuses. The wide base
    # bounds (default 1.25× standard limits) prevent spurious pruning.
    # HOM-proportional widening (applied later, per census pair) extends further.
    if (isTRUE(is_non_taper_corrected)) {
        eff_min_grow <- non_taper_corrected_prune_min_growth
        eff_max_grow <- non_taper_corrected_prune_max_growth
    }

    prune_stats$eff_min_growth <- as.numeric(eff_min_grow)
    prune_stats$eff_max_growth <- as.numeric(eff_max_grow)
    prune_stats$eff_recruit_max <- as.numeric(eff_recruit_max)
    prune_stats$is_non_taper_corrected <- is_non_taper_corrected
    prune_stats$use_hom_relax <- use_hom_relax

    vcat(
        prefix, "Pruning bounds: growth [", eff_min_grow, ", ", eff_max_grow, "] cm/yr",
        ", max recruit DBH=", eff_recruit_max, " cm",
        if (is_non_taper_corrected) " [non-taper-corrected species]" else "",
        if (use_hom_relax) " [HOM-proportional widening active]" else ""
    )

    # --- Precompute per-census max HOM deviation (for HOM-proportional widening) ---
    hom_max_dev_by_census <- NULL
    if (isTRUE(use_hom_relax)) {
        hom_max_dev_by_census <- setNames(numeric(n_census), as.character(census_range))
        for (p_idx in seq_len(n_census)) {
            cc_val <- census_range[p_idx]
            hom_raw <- tree_data[CensusID == cc_val, get(hom_col)]
            # NA HOM treated as 1.3 -> zero deviation
            hom_vals <- ifelse(is.na(hom_raw), 1.3, as.numeric(hom_raw))
            dev_vals <- abs(hom_vals - 1.3)
            mx <- if (length(dev_vals) > 0L) max(dev_vals, na.rm = TRUE) else 0
            if (!is.finite(mx)) mx <- 0
            hom_max_dev_by_census[as.character(cc_val)] <- mx
        }
        vcat(prefix, "HOM deviation by census: ", paste("C", names(hom_max_dev_by_census), "=", round(hom_max_dev_by_census, 3), "m", sep = "", collapse = ", "))
    }

    # Precompute track-wise DBH matrix for each census (rows = states, cols = tracks).
    # Vectorized: one matrix indexing operation per census instead of a per-state loop.
    track_dbh_by_state <- vector("list", n_census)
    for (p in seq_len(n_census)) {
        mat <- state_mats[[p]]
        n_states <- nrow(mat)
        n_obs_p <- ncol(mat)
        tdbh_mat <- matrix(NA_real_, nrow = n_states, ncol = K)
        if (n_obs_p > 0L && n_states > 0L) {
            row_idx_v <- rep(seq_len(n_states), times = n_obs_p)
            col_idx_v <- as.vector(mat) # column-major track indices
            dbh_vals <- rep(obs_dbh[[p]], each = n_states)
            tdbh_mat[cbind(row_idx_v, col_idx_v)] <- dbh_vals
        }
        track_dbh_by_state[[p]] <- tdbh_mat # n_states × K matrix
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
    vcat(prefix, "--- BACKWARD PASS --- anchor C", anchor_start, " -> C", census_range[1L], " (", (anchor_pos - 1L), " interval(s) to evaluate)")

    # Precompute census mean dates as a named numeric vector (avoids dcast overhead).
    .dt_pi <- tree_data[, .(MeanDate = mean(ExactDate, na.rm = TRUE)), by = CensusID]
    pair_interval <- setNames(as.numeric(.dt_pi$MeanDate), as.character(.dt_pi$CensusID))

    # Guard: some tags only have the anchor census (anchor_pos == 1). In that case
    # there are no earlier censuses and calling seq.int(anchor_pos - 1L, 1L, by = -1L)
    # would error with "wrong sign in 'by' argument". Skip the loop when anchor_pos == 1.
    if (anchor_pos > 1L) {
        for (p in seq.int(anchor_pos - 1L, 1L, by = -1L)) {
            # p <- 5L #position in census_range (for testing only)
            cc <- census_range[p]
            next_cc <- census_range[p + 1L]
            mat_cc <- state_mats[[p]]
            n_states_cc <- nrow(mat_cc)

            if (verbose) t_cc0 <- tic()
            vcat(prefix, "  [", (anchor_pos - p), "/", (anchor_pos - 1L), "] C", cc, " -> C", next_cc, ": ", length(obs_dbh[[p]]), " stem(s) x ", format(n_states_cc, big.mark = ","), " states vs ", length(obs_dbh[[p + 1L]]), " stem(s) x ", format(length(keys_full[[p + 1L]]), big.mark = ","), " target states")

            next_keys <- keys_full[[p + 1L]]
            n_next <- length(next_keys)
            if (n_next == 0L) {
                vcat(prefix, "  Dead end: no reachable states at C", next_cc, "; routing to probabilistic matcher")
                return(do_fallback("no_reachable_next_states"))
            }
            next_index <- seq_len(n_next)
            names(next_index) <- next_keys
            logB_next <- as.numeric(logB[[p + 1L]])
            vit_next <- as.numeric(vit_cost[[p + 1L]])

            # Fast lookup for next assignment DBHs via assignment key
            next_assign_list <- assign_full[[p + 1L]]
            # Build assignment-key -> state index for next census (since phase differs but assignment cost uses assignment).
            # Vectorized: rebuild assign list as matrix then use do.call(paste,...) instead of vapply+state_key.
            next_assign_mat <- do.call(rbind, next_assign_list)
            next_assign_key <- if (is.null(next_assign_mat) || ncol(next_assign_mat) == 0L) {
                rep("", nrow(next_assign_mat))
            } else {
                do.call(paste, c(lapply(seq_len(ncol(next_assign_mat)), function(j) next_assign_mat[, j]), list(sep = ",")))
            }
            next_assign_row_idx <- match(next_assign_key, state_keys[[p + 1L]])
            if (any(is.na(next_assign_row_idx))) {
                # Should not happen; indicates mismatch in state enumeration.
                return(do_fallback("next_assign_row_mismatch"))
            }

            # --- Cross-product guard: prevent memory crash when the product
            # of adjacent census state spaces exceeds max_edges.  The C++
            # derive_phase_prev_batch_rcpp would try to allocate n_cc × n_next
            # items, which can exceed available memory for large multi-stem tags.
            # Use >= so that the boundary case (both censuses capped at exactly
            # max_states) also triggers fallback instead of crashing.
            n_cross <- as.double(n_states_cc) * as.double(n_next)
            if (n_cross >= as.double(max_edges)) {
                vcat(
                    prefix, "  Edge overflow: ", format(n_cross, big.mark = ",", scientific = FALSE),
                    " candidate edges (", format(n_states_cc, big.mark = ","), " x ", format(n_next, big.mark = ","),
                    ") exceeds max_edges=", format(as.double(max_edges), big.mark = ",", scientific = FALSE),
                    "; routing to probabilistic matcher"
                )
                return(do_fallback("edge_count_exceeded", K_used = K))
            }

            # Pre-decode next phase vectors
            phase_tp1_by_next <- vector("list", n_next)
            for (j in seq_len(n_next)) {
                phase_tp1_by_next[[j]] <- decode_full_key(next_keys[[j]])$phase
            }

            # Interval (years) between cc and next_cc (pair_interval precomputed above the loop)
            val_next <- pair_interval[[as.character(next_cc)]]
            val_cc <- pair_interval[[as.character(cc)]]
            if (is.null(val_next) || length(val_next) == 0 || is.null(val_cc) || length(val_cc) == 0) {
                interval_val <- NA_real_
            } else {
                interval_val <- (as.numeric(val_next) - as.numeric(val_cc)) / 365.25
                if (!is.finite(interval_val) || interval_val <= 0) {
                    interval_val <- NA_real_
                }
            }
            vcat(prefix, "    dt=", sprintf("%.2f", interval_val), " yr", sep = "")

            # -----------------------------------------------------------------------
            # Batch feasibility check in C++: replaces the O(n_cc × n_next × K) R
            # inner loop.  derive_phase_prev_batch_rcpp checks phase-transition
            # constraints and hard growth-rate pruning for every (i, j) pair at once.
            # -----------------------------------------------------------------------

            # Build input matrices (rows = states/next-states, cols = tracks K)
            tdbh0_mat <- track_dbh_by_state[[p]] # already n_cc   × K
            tdbh1_mat <- track_dbh_by_state[[p + 1L]][next_assign_row_idx, , drop = FALSE] # n_next × K
            phase_tp1_mat <- do.call(rbind, phase_tp1_by_next) # n_next × K; integer

            # Resprout matrix: n_next × K (zero-row means "no resprouts")
            if (any(is_resprout_obs[[p + 1L]])) {
                resprout_flags_tp1 <- is_resprout_obs[[p + 1L]]
                resp_mat <- matrix(FALSE, nrow = n_next, ncol = K)
                for (j in seq_len(n_next)) {
                    assign_tp1 <- next_assign_list[[j]]
                    for (q in seq_along(assign_tp1)) {
                        if (resprout_flags_tp1[q]) resp_mat[j, assign_tp1[q]] <- TRUE
                    }
                }
            } else {
                resp_mat <- matrix(logical(0L), nrow = 0L, ncol = K)
            }

            # --- HOM-proportional widening of prune bounds for this census pair ---
            eff_min_grow_pair <- eff_min_grow
            eff_max_grow_pair <- eff_max_grow
            if (isTRUE(use_hom_relax) && is.finite(interval_val) && interval_val > 0) {
                dev_cc <- hom_max_dev_by_census[as.character(cc)]
                dev_next <- hom_max_dev_by_census[as.character(next_cc)]
                max_dev <- max(dev_cc, dev_next, 0, na.rm = TRUE)
                if (is.finite(max_dev) && max_dev > 0) {
                    hom_tol <- hom_tolerance_scale * max_dev / interval_val
                    eff_min_grow_pair <- eff_min_grow - hom_tol
                    eff_max_grow_pair <- eff_max_grow + hom_tol
                    vcat(
                        prefix, "    HOM widening: max_dev=", round(max_dev, 3), " m -> growth bounds widened to [",
                        round(eff_min_grow_pair, 3), ", ", round(eff_max_grow_pair, 3), "] cm/yr (+/-", round(hom_tol, 3), ")"
                    )
                }
            }

            feasible_result <- derive_phase_prev_batch_rcpp(
                tdbh0_mat = tdbh0_mat,
                tdbh1_mat = tdbh1_mat,
                phase_tp1_mat = phase_tp1_mat,
                resprout_mat = resp_mat,
                prune_hard = isTRUE(prune_hard),
                interval_val = if (is.finite(interval_val)) interval_val else NaN,
                eff_min_grow = eff_min_grow_pair,
                eff_max_grow = eff_max_grow_pair,
                eff_recruit_max = if (is.finite(eff_recruit_max)) eff_recruit_max else Inf
            )

            fe_from <- feasible_result$from_i # 1-based current assignment indices
            fe_to <- feasible_result$to_j # 1-based next full-state indices
            fe_phase <- feasible_result$phase_t # n_feasible × K integer matrix
            n_feasible <- length(fe_from)

            # -----------------------------------------------------------------------
            # Post-segment recruit continuity constraint
            # When this DP is the post-segment of a resprout split, ALL stems at the
            # first census (the R-boundary) are new organisms.  If the number of
            # stems at the first census <= the number at the next census, every stem
            # at p must be on a track that is also occupied at p+1 — a freshly
            # resprouted stem dying immediately while another independent recruit
            # appears is far less parsimonious than all recruits continuing.
            # -----------------------------------------------------------------------
            if (n_feasible > 0L && isTRUE(post_segment_all_recruits) && p == 1L) {
                .n_obs_p <- length(obs_dbh[[p]])
                .n_obs_p1 <- length(obs_dbh[[p + 1L]])
                if (.n_obs_p > 0L && .n_obs_p <= .n_obs_p1) {
                    .p_assign <- state_mats[[p]][fe_from, , drop = FALSE] # n_feasible × n_obs_p
                    .p1_tracks <- next_assign_mat[fe_to, seq_len(.n_obs_p1), drop = FALSE] # n_feasible × n_obs_p1
                    .keep_rc <- rep(TRUE, n_feasible)
                    for (.pj in seq_len(.n_obs_p)) {
                        .trk <- .p_assign[, .pj] # track of obs .pj at p
                        # Check that this track is occupied at p+1
                        .occ <- matrix(FALSE, nrow = n_feasible, ncol = .n_obs_p1)
                        for (.qj in seq_len(.n_obs_p1)) {
                            .occ[, .qj] <- .trk == .p1_tracks[, .qj]
                        }
                        .keep_rc <- .keep_rc & (rowSums(.occ) > 0L)
                    }
                    if (any(!.keep_rc)) {
                        .n_rc_pruned <- sum(!.keep_rc)
                        if (all(!.keep_rc)) {
                            # Safety: constraint would remove ALL transitions → dead end.
                            # Skip the constraint to keep the DP viable; the probabilistic
                            # fallback would lose more accuracy than a slightly relaxed
                            # recruit-continuity assumption.
                            vcat(
                                prefix, "    Recruit continuity: would remove ALL ", n_feasible,
                                " transitions at C", cc, " -> C", next_cc,
                                "; skipping constraint to avoid dead end"
                            )
                        } else {
                            vcat(
                                prefix, "    Recruit continuity: removed ", .n_rc_pruned,
                                "/", n_feasible, " transitions (all ", .n_obs_p,
                                " post-resprout recruit(s) at C", cc,
                                " must persist to C", next_cc, ")"
                            )
                            fe_from <- fe_from[.keep_rc]
                            fe_to <- fe_to[.keep_rc]
                            fe_phase <- fe_phase[.keep_rc, , drop = FALSE]
                            n_feasible <- length(fe_from)
                        }
                    }
                }
            }

            # -----------------------------------------------------------------------
            # R-recruit constraint
            # An R-coded observation at census p+1 is the OLD stem's last record —
            # the track it occupies at p+1 MUST ALREADY BE OCCUPIED at p
            # (it cannot be a new/empty track).
            # Only applies to R-coded LIVING obs (non-NA DBH).
            # Skipped when census p has 0 observations (all tracks already empty).
            # -----------------------------------------------------------------------
            if (n_feasible > 0L) {
                .r_pos_p1 <- which(is_resprout_obs[[p + 1L]]) # R-coded LIVING obs at p+1
                .n_obs_p <- length(obs_dbh[[p]])
                if (length(.r_pos_p1) > 0L && .n_obs_p > 0L) {
                    # R-coded obs at p+1 is old stem's last record: it MUST continue
                    # a track that was occupied at p (not an empty/unborn track).
                    .p_assign <- state_mats[[p]][fe_from, , drop = FALSE] # n_feasible × n_obs_p
                    .keep_r <- rep(TRUE, n_feasible)
                    for (.rj in seq_along(.r_pos_p1)) {
                        .ri_p1 <- .r_pos_p1[.rj]
                        .rt <- next_assign_mat[fe_to, .ri_p1] # track assigned to R-obs .rj at p+1
                        # pair k fails iff .rt[k] does NOT appear in .p_assign[k, ]
                        .was_occupied <- matrix(FALSE, nrow = n_feasible, ncol = .n_obs_p)
                        for (.pc in seq_len(.n_obs_p)) {
                            .was_occupied[, .pc] <- .rt == .p_assign[, .pc]
                        }
                        .keep_r <- .keep_r & (rowSums(.was_occupied) > 0L)
                    }
                    if (any(!.keep_r)) {
                        .n_r_pruned <- sum(!.keep_r)
                        vcat(
                            prefix, "    R-recruit filter: removed ", .n_r_pruned, "/", n_feasible,
                            " transitions at C", next_cc, " (broken-below: ", length(.r_pos_p1), " R-coded obs must continue occupied tracks)"
                        )
                        fe_from <- fe_from[.keep_r]
                        fe_to <- fe_to[.keep_r]
                        fe_phase <- fe_phase[.keep_r, , drop = FALSE]
                        n_feasible <- length(fe_from)
                    }
                }
            }

            # Update prune diagnostics
            if (isTRUE(prune_hard)) {
                n_examined <- n_states_cc * n_next
                n_infeasible <- n_examined - n_feasible
                prune_stats$total_examined <- prune_stats$total_examined + n_examined
                prune_stats$total_pruned <- prune_stats$total_pruned + n_infeasible
                prune_stats$per_census[[as.character(cc)]] <-
                    prune_stats$per_census[[as.character(cc)]] + n_infeasible
            }

            # Compute full-state key strings for each feasible pair using vectorized ops.
            if (n_feasible > 0L) {
                fe_assign_keys <- state_keys[[p]][fe_from] # assign key portion (precomputed)
                # Vectorized phase key encoding: avoids per-row apply + rawToChar overhead
                .phase_chars <- c("0", "1", "2")
                fe_phase_keys <- do.call(paste0, lapply(seq_len(K), function(k) .phase_chars[fe_phase[, k] + 1L]))
                fe_full_keys <- paste0(fe_assign_keys, "|", fe_phase_keys)
            }

            # Dynamic creation of current full-states (batched: no per-i loop)
            curr_keys_list <- list()
            curr_assign_list <- list()
            curr_logB <- numeric(0)
            curr_vit <- numeric(0)
            curr_ptr <- integer(0)

            # Preallocate edge arrays (at most n_feasible edges)
            from_idx <- integer(n_feasible)
            to_idx <- integer(n_feasible)
            logw <- numeric(n_feasible)
            used_edges <- 0L

            # Batched cost computation + vectorized state registration + aggregation
            if (n_feasible > 0L) {
                # --- 1. Single C++ call for ALL feasible transitions ---
                .tdbh0_all <- track_dbh_by_state[[p]][fe_from, , drop = FALSE]
                .tdbh1_all <- tdbh1_mat[fe_to, , drop = FALSE]

                if (verbose) t_tc0 <- tic()
                .all_costs <- transition_cost_paired_rcpp(
                    tdbh0_mat = .tdbh0_all,
                    tdbh1_mat = .tdbh1_all,
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
                transition_cost_calls <- transition_cost_calls + 1L
                if (verbose) transition_cost_time <- transition_cost_time + (tic() - t_tc0)

                # --- 2. Vectorized state registration (unique + match) ---
                .uniq_keys <- unique(fe_full_keys)
                .n_unique <- length(.uniq_keys)
                .edge_sidx <- match(fe_full_keys, .uniq_keys)

                # Extract assignment for each unique state from its first occurrence
                .first_occ <- match(.uniq_keys, fe_full_keys)
                .first_from <- fe_from[.first_occ]
                curr_keys_list <- as.list(.uniq_keys)
                curr_assign_list <- lapply(.first_from, function(i) as.integer(mat_cc[i, ]))

                # --- 3. Vectorized cost accumulation ---
                .logw_all <- -.all_costs / temperature
                .cand_vit_all <- .all_costs + vit_next[fe_to]
                .cand_log_all <- .logw_all + logB_next[fe_to]

                # Viterbi: min cost per state via order + first-in-group
                .ord <- order(.edge_sidx, .cand_vit_all)
                .sorted_sidx <- .edge_sidx[.ord]
                .first_in_grp <- c(TRUE, diff(.sorted_sidx) != 0L)
                .best_idx <- .ord[.first_in_grp]
                curr_vit <- .cand_vit_all[.best_idx]
                curr_ptr <- fe_to[.best_idx]

                # LogB: log-sum-exp per state
                curr_logB <- rep(-Inf, .n_unique)
                .dt_lb <- data.table::data.table(s = .edge_sidx, v = .cand_log_all)
                .lb_agg <- .dt_lb[,
                    {
                        m <- max(v)
                        list(lb = m + log(sum(exp(v - m))))
                    },
                    by = s
                ]
                curr_logB[.lb_agg$s] <- .lb_agg$lb

                # Edge storage (all feasible transitions)
                used_edges <- n_feasible
                from_idx <- .edge_sidx
                to_idx <- fe_to
                logw <- .logw_all
            }

            if (length(curr_keys_list) == 0L) {
                vcat(prefix, "  Dead end: no valid states at C", cc, " after phase/growth filtering; routing to probabilistic matcher")
                return(do_fallback("no_states_produced"))
            }

            keys_full[[p]] <- unlist(curr_keys_list, use.names = FALSE)
            assign_full[[p]] <- curr_assign_list
            logB[[p]] <- curr_logB
            vit_cost[[p]] <- curr_vit
            vit_ptr[[p]] <- curr_ptr

            if (used_edges == 0L) {
                vcat(prefix, "  Dead end: 0 feasible transitions at C", cc, "; routing to probabilistic matcher")
                return(do_fallback("no_feasible_edges"))
            }
            edges[[p]] <- data.table::data.table(
                from_idx = from_idx[seq_len(used_edges)],
                to_idx   = to_idx[seq_len(used_edges)],
                logw     = logw[seq_len(used_edges)]
            )

            vcat(prefix, "  [", (anchor_pos - p), "/", (anchor_pos - 1L), "] done: C", cc, " -> ", format(length(keys_full[[p]]), big.mark = ","), " reachable states, ", format(used_edges, big.mark = ","), " feasible edges", if (verbose) paste0(" (", sprintf("%.2fs", tic() - t_cc0), ")") else "")
        }
    }

    # -----------------
    # Decode MAP path
    # -----------------
    vcat(prefix, "--- MAP DECODING --- Extracting optimal stem identity assignments across ", n_census, " censuses ...")
    # Stable tie-breaker: among states with equal minimum cost, always pick the
    # lowest enumeration index (enumerate_states_injective is deterministic, so
    # this makes which.min deterministic across parallel runs and platforms).
    {
        .vc <- vit_cost[[1L]]
        .min_cost <- min(.vc, na.rm = TRUE)
        .tied <- which(.vc == .min_cost)
        start_idx <- .tied[1L] # lowest index wins any tie
    }
    if (length(start_idx) == 0L || !is.finite(vit_cost[[1L]][start_idx])) {
        return(do_fallback("decode_failure"))
    }
    map_idx <- integer(n_census)
    map_idx[1L] <- start_idx
    for (p in seq_len(n_census - 1L)) {
        nxt <- vit_ptr[[p]][map_idx[p]]
        if (!is.finite(nxt) || is.na(nxt) || nxt < 1L) {
            return(do_fallback("viterbi_decode_failure"))
        }
        map_idx[p + 1L] <- nxt
    }

    for (p in seq_len(n_census)) {
        cc <- census_range[p]
        obs_idx <- obs_row_idx[[p]]
        if (length(obs_idx) == 0L) next
        sv <- assign_full[[p]][[map_idx[p]]]
        if (length(sv) != length(obs_idx)) {
            return(do_fallback("assign_mismatch"))
        }
        tree_data[obs_idx, ReconstructedStemID := track_ids[sv]]
        # Mark ReconstructionMethod:
        #   - Anchor census rows with TrueStemID keep "given" (set by pre-stamp)
        #   - All other rows get "dp" (the solver determined their identity)
        obs_to_mark <- obs_idx[is.na(tree_data$TrueStemID[obs_idx]) |
                               tree_data$CensusID[obs_idx] != anchor_start]
        if (length(obs_to_mark) > 0L) {
            tree_data[obs_to_mark, ReconstructionMethod := "dp"]
        }
    }

    # Restore "given" on pinned non-anchor rows whose identity was
    # constrained by TrueStemID (real, not provisional).  The MAP decode
    # loop above blanket-stamps "dp" on all non-anchor rows but these
    # were forced by the pin constraint.
    if (isTRUE(pin_truestemid)) {
        .pinned_pre <- which(!is.na(tree_data$TrueStemID) &
                             tree_data$CensusID != anchor_start &
                             !tree_data$ReconstructionMethod %in% "provisional_dp")
        if (length(.pinned_pre) > 0L) {
            tree_data[.pinned_pre, ReconstructionMethod := "given"]
        }
        # For post-anchor rows that were never visited by the DP (CensusID > anchor_start,
        # dp_scoped_to_pre_anchor=FALSE), also copy TrueStemID into ReconstructedStemID.
        .pinned_post_unassigned <- which(!is.na(tree_data$TrueStemID) &
                                         tree_data$CensusID > anchor_start &
                                         is.na(tree_data$ReconstructedStemID))
        if (length(.pinned_post_unassigned) > 0L) {
            tree_data[.pinned_post_unassigned, ReconstructedStemID := as.integer(TrueStemID)]
        }
        # ---- Hard-invariant pre-anchor sweep -----------------------------------
        # Same fix as the post-anchor block above but for pre-anchor NA-DBH
        # rows that the DP never visited (obs_row_idx filters !is.na(DBH)).
        # A row carrying TrueStemID at C2 with DBH=NA and Status="dead"
        # would otherwise leave finalize_out with ReconstructedStemID=NA and
        # violate the invariant.  Runs BEFORE the NA-R barrier post-
        # processing below, so the existing barrier guards (which test
        # !is.na(TrueStemID)) see these rows already labelled.
        .pinned_pre_unassigned <- which(!is.na(tree_data$TrueStemID) &
                                        tree_data$CensusID < anchor_start &
                                        is.na(tree_data$ReconstructedStemID))
        if (length(.pinned_pre_unassigned) > 0L) {
            tree_data[.pinned_pre_unassigned, `:=`(
                ReconstructedStemID = as.integer(TrueStemID),
                ReconstructionMethod = "given"
            )]
        }
    }

    # -----------------------------------------------------------------
    # NA-R barrier post-processing
    # When a census has 0 live stems AND NA-DBH R-coded rows, it is a
    # hard resprout boundary: stems before the barrier cannot share the
    # same track identity as stems after it.  Re-assign pre-barrier
    # stems to new synthetic IDs, and assign the NA-R rows themselves
    # (the dying original stem) to those same pre-barrier IDs.
    # -----------------------------------------------------------------
    barrier_positions <- which(.has_na_r_barrier)
    if (length(barrier_positions) > 0L) {
        .cur_max_id <- suppressWarnings(max(tree_data$ReconstructedStemID, na.rm = TRUE))
        if (!is.finite(.cur_max_id)) .cur_max_id <- 0L
        for (.p_bar in barrier_positions) {
            .cc_bar <- census_range[.p_bar]
            .cens_before <- census_range[seq_len(.p_bar - 1L)]
            if (length(.cens_before) == 0L) next # barrier at earliest census, nothing to do
            .ids_before <- unique(tree_data[CensusID %in% .cens_before & !is.na(ReconstructedStemID), ReconstructedStemID])
            .ids_after <- unique(tree_data[CensusID > .cc_bar & !is.na(ReconstructedStemID), ReconstructedStemID])
            .crossing <- intersect(.ids_before, .ids_after)
            # Filter out anchor-validated crossings: a crossing ID is valid (should NOT be
            # severed) when it is a known post-barrier TrueStemID AND no pre-barrier row
            # with that reconstructed ID carries a conflicting (different, non-NA) TrueStemID.
            # This handles the case where the DP correctly traces a known anchor identity
            # backward through a resprout event — severing would destroy that valid connection.
            if ("TrueStemID" %in% names(tree_data) && length(.crossing) > 0L) {
                .known_post_tsids <- unique(tree_data$TrueStemID[
                    tree_data$CensusID > .cc_bar & !is.na(tree_data$TrueStemID)
                ])
                .n_cross_before <- length(.crossing)
                .crossing <- Filter(function(.cid) {
                    # Keep (sever) if this ID is not a known post-barrier anchor identity
                    if (!(.cid %in% .known_post_tsids)) return(TRUE)
                    # Keep (sever) if any pre-barrier row with this reconstructed ID
                    # has a known TrueStemID that differs (i.e. DP connected two distinct plants)
                    .pre_tsids <- tree_data$TrueStemID[
                        tree_data$CensusID %in% .cens_before &
                        tree_data$ReconstructedStemID == .cid &
                        !is.na(tree_data$TrueStemID)
                    ]
                    any(.pre_tsids != .cid)
                }, .crossing)
                .n_skip <- .n_cross_before - length(.crossing)
                if (.n_skip > 0L) {
                    vcat(prefix, "NA-R barrier at C", .cc_bar, ": ", .n_skip,
                         " anchor-validated track(s) excluded from sever",
                         " (DP correctly traced known post-barrier TrueStemID backward)")
                }
            }
            if (length(.crossing) > 0L) {
                vcat(
                    prefix, "NA-R barrier at C", .cc_bar, ": ", length(.crossing),
                    " track(s) cross the resprout boundary -> assigning new synthetic IDs to pre-barrier stems"
                )
                for (.old_id in .crossing) {
                    .new_id <- as.integer(.cur_max_id) + 1L
                    .cur_max_id <- .new_id
                    # Protect pinned rows: only reassign obs without real
                    # TrueStemID (provisional counts as non-pinned).
                    .pinned_crossing <- tree_data[
                        CensusID %in% .cens_before &
                        ReconstructedStemID == .old_id &
                        !is.na(TrueStemID) &
                        !ReconstructionMethod %in% "provisional_dp"
                    ]
                    if (nrow(.pinned_crossing) > 0L) {
                        vcat(prefix, "WARNING: pinned TrueStemID ",
                             paste(unique(.pinned_crossing$TrueStemID), collapse = ","),
                             " crosses NA-R barrier at C", .cc_bar,
                             " — respecting pin (not reassigning)")
                    }
                    tree_data[
                        CensusID %in% .cens_before &
                        ReconstructedStemID == .old_id &
                        (is.na(TrueStemID) | ReconstructionMethod %in% "provisional_dp"),
                        ReconstructedStemID := .new_id
                    ]
                    tree_data[CensusID %in% .cens_before & ReconstructedStemID == .new_id &
                        is.na(TrueStemID), ReconstructionMethod := "dp"]
                }
            }
            # Assign NA-R rows at the barrier census itself to the (now re-IDed)
            # pre-barrier identities — they represent the original stem dying.
            {
                .has_status_col <- "Status" %in% names(tree_data)
                # Compute boolean flags outside data.table[i] to avoid
                # "Object not found" errors when the column is absent.
                .n_td <- nrow(tree_data)
                .is_r_bar <- if (.has_tsm) !is.na(tree_data$ListOfTSM) & grepl(resprout_regex, tree_data$ListOfTSM, perl = TRUE) else rep(FALSE, .n_td)
                .is_bb_bar <- if (.has_status_col) !is.na(tree_data$Status) & tree_data$Status == "broken below" else rep(FALSE, .n_td)
                .barrier_na_r_rows <- which(
                    tree_data$CensusID == .cc_bar & is.na(tree_data$DBH) & (.is_r_bar | .is_bb_bar)
                )
                .prev_cc <- census_range[.p_bar - 1L]
                .prev_ids <- tree_data[
                    CensusID == .prev_cc & !is.na(ReconstructedStemID),
                    ReconstructedStemID
                ]
                if (length(.barrier_na_r_rows) > 0L && length(.prev_ids) > 0L) {
                    for (.ri in seq_along(.barrier_na_r_rows)) {
                        .pid <- if (.ri <= length(.prev_ids)) .prev_ids[.ri] else .prev_ids[length(.prev_ids)]
                        tree_data[.barrier_na_r_rows[.ri], ReconstructedStemID := .pid]
                    }
                    tree_data[
                        .barrier_na_r_rows[is.na(tree_data$TrueStemID[.barrier_na_r_rows])],
                        ReconstructionMethod := "dp"
                    ]
                    vcat(
                        prefix, "NA-R barrier at C", .cc_bar, ": linked ",
                        length(.barrier_na_r_rows), " dying-stem row(s) to pre-barrier identities"
                    )
                }
            }
        }
    }

    # -----------------
    # Forward pass for marginals
    # -----------------
    vcat(prefix, "--- FORWARD PASS --- C", census_range[1L], " -> anchor C", anchor_start, " (computing posterior marginals for uncertainty quantification)")
    # Start distribution: uniform over all reachable states at census 1.
    logalpha <- vector("list", n_census)
    logalpha[[1L]] <- rep.int(0, length(keys_full[[1L]]))

    for (p in seq_len(n_census - 1L)) {
        ed <- edges[[p]]
        if (is.null(ed) || nrow(ed) == 0L) {
            return(do_fallback("forward_edges_missing"))
        }
        la_from <- logalpha[[p]][ed$from_idx]
        vals <- la_from + ed$logw
        dt <- data.table::data.table(to_idx = ed$to_idx, v = vals)
        dt <- dt[is.finite(v)]
        if (nrow(dt) == 0L) {
            return(do_fallback("forward_no_alpha"))
        }
        la_next_dt <- dt[, .(logalpha = log_sum_exp(v)), by = to_idx]
        la_next <- rep.int(-Inf, length(keys_full[[p + 1L]]))
        la_next[la_next_dt$to_idx] <- la_next_dt$logalpha
        logalpha[[p + 1L]] <- la_next

        vcat(prefix, "  [", p, "/", (n_census - 1L), "] C", census_range[p], " -> C", census_range[p + 1L], ": ", sum(is.finite(la_next)), "/", length(la_next), " states reachable (", sprintf("%.0f%%", 100 * sum(is.finite(la_next)) / max(1L, length(la_next))), ")")
    }

    # Partition function Z = total weight of all paths ending at the fixed anchor state.
    # At anchor_start there is exactly one state.
    logZ <- logalpha[[anchor_pos]][1L]
    if (!is.finite(logZ)) {
        # Fallback: compute from backward at census 1.
        logZ <- log_sum_exp(logB[[1L]])
    }

    vcat(prefix, "Partition function: logZ = ", sprintf("%.3f", logZ), " (normalises posterior over all valid reconstructions)")

    # -----------------
    # Observation-level marginals
    # -----------------
    anchor_set <- anchor_ids
    is_anchor_track <- track_ids %in% anchor_set

    for (p in seq_len(n_census)) {
        cc <- census_range[p]
        obs_idx <- obs_row_idx[[p]]
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
        tag_local <- if (!is.na(tag_val)) tag_val else get0("WHICH_TAG", ifnotfound = NA_integer_)

        # Collect simple sampling CPU/memory metrics and attach them to the returned object
        sampling_profile <- list(posterior_samples = posterior_samples, started = Sys.time())
        t_sampling_start <- tic()

        # Build adjacency lookup for quick sampling
        t_adj_start <- tic()
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
        sampling_profile$adj_build_time <- tic() - t_adj_start

        # sampler: backwards sampling from anchor_pos to census 1
        if (!is.null(posterior_sample_seed)) set.seed(as.integer(posterior_sample_seed))
        t_gen_start <- tic()
        samples_list <- vector("list", posterior_samples)
        for (m in seq_len(posterior_samples)) {
            sampled_idx <- integer(n_census)
            sampled_idx[anchor_pos] <- 1L # anchor position has single full-state at index 1
            logp <- 0
            # sample backwards
            if (anchor_pos > 1L) {
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
            }
            # Convert sampled full-state indices to ReconstructedStemID per census
            sample_dt <- data.table::data.table(Tag = tag_local, Sample = m, CensusID = integer(0), ReconstructedStemID = integer(0), ObsRowID = integer(0))
            for (p in seq_len(n_census)) {
                cc <- census_range[p]
                assign_vec <- assign_full[[p]][[sampled_idx[p]]]
                track_ids_loc <- track_ids[assign_vec]
                # attach as multiple rows (one per observed tree)
                obs_idx <- obs_row_idx[[p]]
                if (length(obs_idx) > 0) {
                    obs_row_ids <- tree_data$obs_row_id[obs_idx]
                    sample_dt <- rbind(sample_dt, data.table::data.table(Tag = tag_local, Sample = m, CensusID = rep(cc, length(obs_idx)), ReconstructedStemID = track_ids_loc, ObsRowID = obs_row_ids))
                }
            }
            sample_dt[, logp := logp]
            samples_list[[m]] <- sample_dt
        }
        sampling_profile$sample_generation_time <- tic() - t_gen_start
        sampling_profile$samples_list_size_bytes <- as.numeric(object.size(samples_list))

        samples_dt <- data.table::rbindlist(samples_list, use.names = TRUE, fill = TRUE)
        # Ensure rows are ordered for signature construction
        samples_dt <- samples_dt[order(Sample, CensusID)]
        # Strategy A: enforce broken-below R1/R2 invariants per posterior sample
        # so MAP and posteriors agree on contract-compliant trajectories.
        # path_sig is computed AFTER this relabel.
        samples_dt <- apply_bb_invariants_to_samples(samples_dt, tree_data, verbose = FALSE)
        sampling_profile$samples_dt_size_bytes <- as.numeric(object.size(samples_dt))
        sampling_profile$after_samples_time <- tic() - t_sampling_start

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
        vcat(prefix, sprintf("Posterior sampling: drew %d samples, found %d unique reconstructions (%.0f%% path diversity)", posterior_samples, n_unique_paths, 100 * n_unique_paths / posterior_samples))
        if (uniqueN(sample_logp$logp) == 1L) {
            vcat(prefix, sprintf("Warning: all %d samples converged to same log-prob (%.3f); posterior may be concentrated on a single reconstruction", posterior_samples, sample_logp$logp[1]))
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
        # Enforce ObsRowID-based reconstructions only.
        if (!("ObsRowID" %in% names(samples_dt))) {
            stop("Posterior sampling must include ObsRowID values; regenerate with posterior_samples that preserve observation row identifiers")
        }
        recon_by_path <- samples_dt[, .(recon = paste0(ObsRowID, ":", ReconstructedStemID, collapse = ";")), by = .(path_sig, Sample)]
        # take first recon per path_sig (they are identical across samples with same path_sig)
        recon_compact <- recon_by_path[, .SD[1], by = path_sig, .SDcols = "recon"]
        paths_summary <- merge(paths_summary, recon_compact, by = "path_sig", all.x = TRUE)

        # Export only the paths summary; skip writing any summary file
        sampling_profile$export_paths <- character(0)
        sampling_profile$export_time_seconds <- 0
        if (fmt == "feather" && requireNamespace("arrow", quietly = TRUE)) {
            p2 <- paste0(out_path_base, "_paths.feather")
            t0 <- tic()
            arrow::write_feather(paths_summary, p2)
            t1 <- tic()
            sampling_profile$export_time_seconds <- sampling_profile$export_time_seconds + (t1 - t0)
            sampling_profile$export_paths <- c(sampling_profile$export_paths, p2)
            vcat(prefix, "Wrote posterior paths summary to: ", p2)
        } else if (fmt == "csv") {
            p2 <- paste0(out_path_base, "_paths.csv")
            t0 <- tic()
            data.table::fwrite(paths_summary, p2)
            t1 <- tic()
            sampling_profile$export_time_seconds <- sampling_profile$export_time_seconds + (t1 - t0)
            sampling_profile$export_paths <- c(sampling_profile$export_paths, p2)
            vcat(prefix, "Wrote posterior paths summary to: ", p2)
        } else {
            # rds variant: save only paths table
            p1 <- paste0(out_path_base, "_paths.rds")
            t0 <- tic()
            saveRDS(paths_summary, file = p1)
            t1 <- tic()
            sampling_profile$export_time_seconds <- sampling_profile$export_time_seconds + (t1 - t0)
            sampling_profile$export_paths <- c(sampling_profile$export_paths, p1)
            vcat(prefix, "Wrote posterior paths summary to: ", p1)
        }
        sampling_profile$export_total_size_bytes <- if (length(sampling_profile$export_paths) > 0) sum(file.size(sampling_profile$export_paths)) else 0L
        sampling_profile$finished <- Sys.time()

        # Attach sampling profile to the returned data.table for external inspection
        attr(tree_data, "DP_Sampling_Profile") <- sampling_profile
    }

    # If DP was scoped to pre-anchor only, merge original post-anchor rows back into the returned dataset
    if (exists("dp_scoped_to_pre_anchor", inherits = FALSE) && isTRUE(dp_scoped_to_pre_anchor)) {
        processed <- tree_data
        post_rows <- original_tree_data[CensusID > anchor_start]
        # Ensure both parts have the same set of columns (add missing ones as NA)
        all_cols <- union(names(processed), names(post_rows))
        missing_in_processed <- setdiff(all_cols, names(processed))
        missing_in_post <- setdiff(all_cols, names(post_rows))
        if (length(missing_in_processed) > 0L) {
            for (c in missing_in_processed) processed[, (c) := NA]
        }
        if (length(missing_in_post) > 0L) {
            for (c in missing_in_post) post_rows[, (c) := NA]
        }
        # Trust all observed TrueStemID at post-anchor censuses (live or dead)
        post_rows[!is.na(TrueStemID), `:=`(
            ReconstructedStemID = as.integer(TrueStemID),
            ReconstructionMethod = "given"
        )]
        # Also assign broken-below / stump rows with known TrueStemID (NA-DBH R-coded rows).
        # R-coded rows are the old stem's last record and must be assigned the same identity.
        {
            .has_status_pa2 <- "Status" %in% names(post_rows)
            .has_tsm_pa2    <- "ListOfTSM" %in% names(post_rows)
            if (.has_status_pa2 || .has_tsm_pa2) {
                .resprout_re_pa2 <- "\\b(R|RP|RF|RT|OR)\\b"
                .is_bb_pa2    <- if (.has_status_pa2) !is.na(post_rows$Status) & post_rows$Status == "broken below" else rep(FALSE, nrow(post_rows))
                .is_r_tsm_pa2 <- if (.has_tsm_pa2) !is.na(post_rows$ListOfTSM) & grepl(.resprout_re_pa2, post_rows$ListOfTSM, perl = TRUE) else rep(FALSE, nrow(post_rows))
                .stump_rows_pa2 <- which(is.na(post_rows$ReconstructionMethod) & is.na(post_rows$DBH) &
                                         !is.na(post_rows$TrueStemID) & (.is_bb_pa2 | .is_r_tsm_pa2))
                if (length(.stump_rows_pa2) > 0L) {
                    post_rows[.stump_rows_pa2, `:=`(
                        ReconstructedStemID = as.integer(TrueStemID),
                        ReconstructionMethod = "given"
                    )]
                }
            }
        }
        # Propagate the pre-anchor stem identity through the initial consecutive broken-below stump
        # chain when TrueStemID is unknown (e.g., field did not record the stump identity).
        # Only applies when the stump chain begins at the first post-anchor census and the census
        # contains no live stems (so the continued identity is unambiguous).
        {
            .has_status_chain <- "Status" %in% names(post_rows)
            .has_tsm_chain    <- "ListOfTSM" %in% names(post_rows)
            if ((.has_status_chain || .has_tsm_chain) && "CensusID" %in% names(post_rows)) {
                .resprout_re_ch <- "\\b(R|RP|RF|RT|OR)\\b"
                .is_bb_ch    <- if (.has_status_chain) !is.na(post_rows$Status) & post_rows$Status == "broken below" else rep(FALSE, nrow(post_rows))
                .is_r_tsm_ch <- if (.has_tsm_chain) !is.na(post_rows$ListOfTSM) & grepl(.resprout_re_ch, post_rows$ListOfTSM, perl = TRUE) else rep(FALSE, nrow(post_rows))
                .is_live_ch  <- !is.na(post_rows$DBH)
                # Unassigned stump rows: broken-below, NA DBH, NA TrueStemID, not yet assigned
                .is_stump_ch <- (.is_bb_ch | .is_r_tsm_ch) & !.is_live_ch & is.na(post_rows$TrueStemID) & is.na(post_rows$ReconstructionMethod)
                if (any(.is_stump_ch)) {
                    # Identify the single pre-anchor stem at the anchor census (unambiguous propagation only)
                    .pre_live_at_anchor <- processed[!is.na(DBH) & CensusID == anchor_start]
                    .pre_ids <- unique(na.omit(.pre_live_at_anchor$ReconstructedStemID))
                    if (length(.pre_ids) == 1L) {
                        .last_pre_id <- .pre_ids[1L]
                        # Walk post-anchor censuses in order; collect a consecutive stump chain
                        # that starts at the first post-anchor census and contains no live stems
                        .post_cens_ord <- sort(unique(post_rows$CensusID))
                        .chain_cens <- integer(0L)
                        for (.cc in .post_cens_ord) {
                            .rows_cc <- which(post_rows$CensusID == .cc)
                            if (any(.is_live_ch[.rows_cc])) break  # Live stem → stop propagation
                            if (any(.is_stump_ch[.rows_cc])) .chain_cens <- c(.chain_cens, .cc)
                        }
                        if (length(.chain_cens) > 0L) {
                            .chain_rows <- which(post_rows$CensusID %in% .chain_cens & .is_stump_ch)
                            post_rows[.chain_rows, `:=`(
                                ReconstructedStemID = .last_pre_id,
                                ReconstructionMethod = "given"
                            )]
                        }
                    }
                }
            }
        }
        # Default remaining post-anchor rows to 'none_after_anchor' if not already set
        if (!("ReconstructionMethod" %in% names(post_rows))) post_rows[, ReconstructionMethod := NA_character_]
        post_rows[is.na(ReconstructionMethod), ReconstructionMethod := "none_after_anchor"]
        if (!("ConstraintViolation" %in% names(post_rows))) post_rows[, ConstraintViolation := as.logical(NA)]
        post_rows <- ensure_posterior_columns(post_rows)
        post_rows[, `:=`(
            DP_KUsed = NA_integer_,
            DP_MaxStatesPerCensus = NA_real_,
            DP_MaxStatesCensusID = NA_integer_
        )]
        tree_data <- data.table::rbindlist(list(processed, post_rows), use.names = TRUE, fill = TRUE)
        # Restore original input order using obs_row_id
        if ("obs_row_id" %in% names(tree_data)) setorder(tree_data, obs_row_id)
        vcat(prefix, "Appended ", nrow(post_rows), " post-anchor row(s) to output (unchanged)")
    }

    # Final pass: propagate the anchor stem identity through any remaining unassigned broken-below
    # stump chain in tree_data. This covers the case where the anchor was auto-adjusted to the
    # last live census (dp_scoped_to_pre_anchor=FALSE) and subsequent all-stump censuses were
    # never processed by the DP (they had no DBH rows to assign).
    {
        .has_status_fp <- "Status" %in% names(tree_data)
        .has_tsm_fp    <- "ListOfTSM" %in% names(tree_data)
        if ((.has_status_fp || .has_tsm_fp) && "CensusID" %in% names(tree_data)) {
            .resprout_re_fp <- "\\b(R|RP|RF|RT|OR)\\b"
            .is_bb_fp    <- if (.has_status_fp) !is.na(tree_data$Status) & tree_data$Status == "broken below" else rep(FALSE, nrow(tree_data))
            .is_r_tsm_fp <- if (.has_tsm_fp) !is.na(tree_data$ListOfTSM) & grepl(.resprout_re_fp, tree_data$ListOfTSM, perl = TRUE) else rep(FALSE, nrow(tree_data))
            .is_live_fp  <- !is.na(tree_data$DBH)
            # Unassigned stump rows: broken-below, NA DBH, NA TrueStemID, no ReconstructionMethod yet
            .is_stump_fp <- (.is_bb_fp | .is_r_tsm_fp) & !.is_live_fp & is.na(tree_data$TrueStemID) & is.na(tree_data$ReconstructionMethod)
            if (any(.is_stump_fp)) {
                # Identify the single stem assigned at the anchor census (unambiguous propagation only)
                .anchor_live <- tree_data[!is.na(DBH) & CensusID == anchor_start]
                .anchor_ids <- unique(na.omit(.anchor_live$ReconstructedStemID))
                if (length(.anchor_ids) == 1L) {
                    .last_id_fp <- .anchor_ids[1L]
                    # Walk post-anchor censuses in order; propagate through consecutive stump chain
                    .post_cens_fp <- sort(unique(tree_data$CensusID[tree_data$CensusID > anchor_start]))
                    .chain_cens_fp <- integer(0L)
                    for (.cc in .post_cens_fp) {
                        .rows_cc <- which(tree_data$CensusID == .cc)
                        if (any(.is_live_fp[.rows_cc])) break
                        if (any(.is_stump_fp[.rows_cc])) .chain_cens_fp <- c(.chain_cens_fp, .cc)
                    }
                    if (length(.chain_cens_fp) > 0L) {
                        .chain_rows_fp <- which(tree_data$CensusID %in% .chain_cens_fp & .is_stump_fp)
                        tree_data[.chain_rows_fp, `:=`(
                            ReconstructedStemID = .last_id_fp,
                            ReconstructionMethod = "given"
                        )]
                    }
                }
            }
        }
    }

    # Attach pruning diagnostics
    attr(tree_data, "DP_PruneInfo") <- prune_stats
    if (prune_stats$total_examined > 0L) {
        vcat(prefix, "Pruning summary: removed ", format(prune_stats$total_pruned, big.mark = ","), " of ", format(prune_stats$total_examined, big.mark = ","), " candidate transitions (", sprintf("%.1f%%", 100 * prune_stats$total_pruned / prune_stats$total_examined), " pruned)")
    }

    # Attach compute profiling (transition cost timing)
    compute_profile <- list(
        transition_cost_total_seconds = as.numeric(transition_cost_time),
        transition_cost_calls = as.integer(transition_cost_calls),
        transition_cost_avg_seconds = if (transition_cost_calls > 0L) as.numeric(transition_cost_time) / as.integer(transition_cost_calls) else NA_real_
    )
    attr(tree_data, "DP_Compute_Profile") <- compute_profile

    # If a sampling profile exists (posterior sampling was requested), embed compute_profile into it for convenience
    if (exists("sampling_profile", inherits = FALSE) && is.list(sampling_profile)) {
        sampling_profile$transition_cost_total_seconds <- compute_profile$transition_cost_total_seconds
        sampling_profile$transition_cost_calls <- compute_profile$transition_cost_calls
        sampling_profile$transition_cost_avg_seconds <- compute_profile$transition_cost_avg_seconds
        attr(tree_data, "DP_Sampling_Profile") <- sampling_profile
    }

    # Re-insert stashed MF episode rows (if any)
    if (isTRUE(has_mf_stash) && !is.null(mf_stash) && nrow(mf_stash) > 0L) {
        tree_data <- reinsert_mf_rows(tree_data, mf_stash)
    }

    # Restore original TrueStemID for provisional_dp rows (undo fabricated anchor IDs).
    # Mirrors the same block in finalize_out(), which is NOT called on the normal DP path.
    {
        .prov_rows_end <- which(tree_data$ReconstructionMethod == "provisional_dp")
        if (length(.prov_rows_end) > 0L &&
            "obs_row_id" %in% names(tree_data) &&
            "obs_row_id" %in% names(original_tree_data)) {
            .orig_tsid_end <- original_tree_data$TrueStemID[
                match(tree_data$obs_row_id[.prov_rows_end], original_tree_data$obs_row_id)
            ]
            tree_data[.prov_rows_end, TrueStemID := .orig_tsid_end]
        }
    }

    vcat(prefix, "--- DONE --- ", sum(!is.na(tree_data$ReconstructedStemID)), " observations mapped to ", length(unique(tree_data$ReconstructedStemID[!is.na(tree_data$ReconstructedStemID)])), " identity track(s) in ", sprintf("%.2fs", tic() - t_start))
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
# Columns: path_sig, path_count, path_prob, recon (compact reconstruction mapping like "<ObsRowID>:<ReconstructedStemID>;..." — ObsRowID is enforced)

# How to use these for error propagation (suggestions)
# Use paths.csv directly: each row is a unique reconstruction with probability path_prob (sums to 1 across unique paths) — convenient for expectation of downstream metrics without resampling.
# Or sample reconstructions according to sample_prob in the summary file to create Monte Carlo realizations for error propagation; then expand each sample using the full long file if needed to attach per-census reconstructed IDs.
# The recon column in paths.csv is handy to quickly apply a mapping (parse "ObsRowID:ReconstructedStemID" pairs).
