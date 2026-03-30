############################################################
# dp_global_dp.R
# Core dynamic programming (MAP and marginal DP functions)
############################################################

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
                                                           # --- NEW: allow DP to use a provisional anchor at the last observed DBH census when no TrueStemID exists ---
                                                           allow_provisional_anchor = TRUE,
                                                           use_measurement_error = FALSE,
                                                           meas_sd1_a = 0.0062,
                                                           meas_sd1_b = 0.0904,
                                                           meas_sd2 = 4.64,
                                                           meas_p_big = 0.05,
                                                           # --- growth-form based fallback ---
                                                           # vector of values in `growth_form` column that should trigger
                                                           # immediate igraph fallback and avoid DP entirely
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
                                                           # NOTE: Check the TODO in the code about possibly adding a margin when using biological recruit max DBH
                                                           prune_use_bio_recruit = TRUE,
                                                           # --- palm-specific tight prune bounds ---
                                                           # When the tree is a palm (growth_form == "palm"), override eff_min_grow / eff_max_grow
                                                           # with these tight bounds (DBH is stable, not growing).
                                                           palm_prune_min_growth = -0.5,
                                                           palm_prune_max_growth =  0.5,
                                                           verbose = FALSE,
                                                           chunk_id = NULL) {
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
        .nms  <- names(dt)
        .ptk  <- max(1L, as.integer(posterior_top_k))
        .nids <- paste0("DP_PosteriorTop", seq_len(.ptk), "ID")
        .npbs <- paste0("DP_PosteriorTop", seq_len(.ptk), "Prob")
        if (all(c(.nids, .npbs, "DP_PosteriorEntropy", "DP_PosteriorReconstructedProb",
                  "DP_PosteriorUnlinkedProb", "obs_row_id") %in% .nms) &&
            all(vapply(.nids, function(n) is.integer(dt[[n]]),  logical(1L))) &&
            all(vapply(.npbs, function(n) is.numeric(dt[[n]]),  logical(1L))) &&
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

    # Collect row indices to stash (union of explicit and implicit MF)
    episode_idx <- integer(0)

    # --- A. Explicit MF detection (ListOfTSM contains "MF") ---
    if ("ListOfTSM" %in% names(tree_data)) {
        is_mf_anchor <- is.na(tree_data$DBH) &
            !is.na(tree_data$ListOfTSM) &
            grepl("\\bMF\\b", tree_data$ListOfTSM)

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
                        # Entire census is all-NA DBH → continuation of MF episode
                        episode_idx <- c(episode_idx, rows_at_next)
                    } else {
                        # Census has at least one measured DBH → episode ends
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
    all_censuses_sorted <- sort(unique(tree_data$CensusID))
    censuses_with_dbh <- sort(unique(tree_data$CensusID[!is.na(tree_data$DBH)]))

    if (length(censuses_with_dbh) >= 2L) {
        first_measured <- min(censuses_with_dbh)
        last_measured  <- max(censuses_with_dbh)

        for (cc in all_censuses_sorted) {
            # Only consider censuses strictly between the first and last measured
            if (cc <= first_measured || cc >= last_measured) next
            rows_at_cc <- which(tree_data$CensusID == cc)
            if (length(rows_at_cc) == 0L) next
            if (all(is.na(tree_data$DBH[rows_at_cc]))) {
                # Sandwiched all-NA census → implicit MF
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
        vcat(prefix, "Removed ", nrow(mf_stash), " missing-from-field (MF) row(s) before DP (will re-insert after)")

        # Identify censuses that became fully empty after MF removal
        mf_emptied_censuses <- integer(0)
        for (cc in unique(mf_stash$CensusID)) {
            if (nrow(tree_data[CensusID == cc & !is.na(DBH)]) == 0L) {
                mf_emptied_censuses <- c(mf_emptied_censuses, as.integer(cc))
            }
        }
        if (length(mf_emptied_censuses) > 0L) {
            vcat(prefix, "Censuses with no remaining observations after MF removal: ",
                 paste(mf_emptied_censuses, collapse = ", "))
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
            after_c  <- dp_censuses[dp_censuses > mf_c]
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
                    data.table::set(stash, ri, "ReconstructedStemID",
                        as.integer(candidates[[1L]]))
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
    propagate_post_anchor_given <- function(post, used_ids = NULL) {
        # used_ids: NULL => treat any TrueStemID+DBH as given (no DP performed)
        if (!("ReconstructedStemID" %in% names(post))) post[, ReconstructedStemID := as.integer(NA_integer_)]
        if (!("ReconstructionMethod" %in% names(post))) post[, ReconstructionMethod := NA_character_]
        if (!("ConstraintViolation" %in% names(post))) post[, ConstraintViolation := as.logical(NA)]

        if (is.null(used_ids)) {
            # No DP output available; treat observed TrueStemID as given
            post[!is.na(TrueStemID) & !is.na(DBH), `:=`(
                ReconstructedStemID = as.integer(TrueStemID),
                ReconstructionMethod = "given"
            )]
        } else {
            if (length(used_ids) > 0L) {
                post[!is.na(TrueStemID) & !is.na(DBH) & (TrueStemID %in% used_ids), `:=`(
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
        if (isTRUE(dp_scoped_to_pre_anchor)) {
            post <- original_tree_data[CensusID > anchor_start]
            if (nrow(post) > 0L) {
                used_ids <- unique(out$ReconstructedStemID[!is.na(out$ReconstructedStemID)])
                post <- propagate_post_anchor_given(post, used_ids)
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
        # Attach prune stats if available (some early returns may occur before prune_stats is initialized)
        attr(out, "DP_PruneInfo") <- if (exists("prune_stats")) prune_stats else list()
        out
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
        vcat(prefix, "Anchor census ", anchor_start, " is beyond last observed census ", last_obs_census, "; adjusting anchor to census ", last_obs_census)
        anchor_start <- last_obs_census
    } else if (!is.na(last_obs_census) && last_obs_census > anchor_requested) {
        vcat(prefix, "Post-anchor observations found (last census=", last_obs_census, "); running DP on censuses 1-", anchor_start, " only, preserving later rows unchanged")
        dp_scoped_to_pre_anchor <- TRUE
        # restrict tree_data to pre-anchor censuses for DP computation
        tree_data <- original_tree_data[CensusID <= anchor_start]

        # If scoping removes all observations (e.g., all observations are after the requested
        # anchor), return a one-row placeholder so the Tag is represented in outputs.
        if (nrow(tree_data) == 0L) {
            vcat(prefix, "No observations before anchor census ", anchor_start, "; skipping DP and returning rows as-is")
            out <- data.table::copy(original_tree_data)
            # Use helper to normalize post-anchor rows: no DP performed, so treat observed TrueStemID as given
            out <- propagate_post_anchor_given(out, used_ids = NULL)
            if (isTRUE(has_mf_stash) && !is.null(mf_stash) && nrow(mf_stash) > 0L) {
                out <- reinsert_mf_rows(out, mf_stash)
            }
            attr(out, "DP_PruneInfo") <- prune_stats
            return(out)
        }
    }

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

    # Ensure a defensive `prune_stats` exists so any early-return branches can
    # safely attach it to outputs without failing when pruning hasn't run.
    prune_stats <- list(
        total_examined = 0L,
        total_pruned = 0L,
        per_census = integer(0)
    )

    if (is.na(first_obs_census)) {
        vcat(prefix, "No DBH observations found up to anchor census ", anchor_start, "; falling back to igraph matcher")
        fallback_reason <- "no_obs_up_to_anchor"
        K_used <- as.integer(min(0L, max_tracks))
        tree_data[, `:=`(
            DP_KUsed = K_used,
            DP_MaxStatesPerCensus = 0L,
            DP_MaxStatesCensusID = NA_integer_
        )]
        out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, anchor_start)
        out <- ensure_posterior_columns(out)
        if (!("DP_FallbackReason" %in% names(out))) out[, DP_FallbackReason := NA_character_]
        out[, DP_FallbackReason := fallback_reason]
        attr(out, "DP_PruneInfo") <- prune_stats
        return(finalize_out(out))
    }
    census_range <- seq.int(from = first_obs_census, to = anchor_start)
    # Exclude censuses fully emptied by MF removal so the DP bridges them
    if (isTRUE(has_mf_stash) && length(mf_emptied_censuses) > 0L) {
        census_range <- census_range[!census_range %in% mf_emptied_censuses]
    }
    n_census <- length(census_range)
    vcat(prefix, "Census range: ", paste(census_range, collapse = ", "), " (first observed=", first_obs_census, ", anchor=", anchor_start, ")")
    obs_counts <- vapply(
        census_range,
        function(cc) nrow(tree_data[CensusID == cc & !is.na(DBH)]),
        integer(1L)
    )
    max_obs <- if (length(obs_counts) > 0L) max(obs_counts) else 0L

    # --- growth-form fallback check ------------------------------------
    if (length(fallback_growth_forms) > 0L && "growth_form" %in% names(tree_data)) {
        bad_idx <- which(tree_data$growth_form %in% fallback_growth_forms)
        if (length(bad_idx) > 0L) {
            vcat(prefix, "Growth form requires igraph fallback (detected: ", paste(unique(tree_data$growth_form[bad_idx]), collapse = ", "), "); skipping DP")
            fallback_reason <- "growth_form_forced"
            K_used <- as.integer(min(max_obs, max_tracks))
            n_states_by_census <- vapply(obs_counts, function(n_obs) count_injective_states(K_used, n_obs), numeric(1L))
            tree_data[, `:=`(
                DP_KUsed = K_used,
                DP_MaxStatesPerCensus = max(n_states_by_census, na.rm = TRUE),
                DP_MaxStatesCensusID = as.integer(census_range[which.max(n_states_by_census)])
            )]
            out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, anchor_start)
            out <- ensure_posterior_columns(out)
            if (!("DP_FallbackReason" %in% names(out))) out[, DP_FallbackReason := NA_character_]
            out[, DP_FallbackReason := fallback_reason]
            attr(out, "DP_PruneInfo") <- prune_stats
            return(finalize_out(out))
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
    palm_prune_min_growth <- as.numeric(palm_prune_min_growth)
    palm_prune_max_growth <- as.numeric(palm_prune_max_growth)

    # --- Palm detection (tight DBH-stability prune bounds when growth_form == "palm")
    if ("growth_form" %in% names(tree_data)) {
        gf_vals <- unique(tree_data$growth_form)
        if (length(gf_vals) > 1L) {
            warning(prefix, "Multiple growth_form values found; using most common to determine palm status.")
            gf_vals <- names(sort(table(tree_data$growth_form), decreasing = TRUE))[1L]
        }
        is_palm <- isTRUE(gf_vals == "palm")
    } else {
        is_palm <- FALSE
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
    # as the anchor instead of immediately falling back to the igraph matcher.
    anchor_rows_all <- tree_data[CensusID == anchor_start]
    if (nrow(anchor_rows_all) == 0L || (all(is.na(anchor_rows_all$DBH)) && all(is.na(anchor_rows_all$TrueStemID)))) {
        cand_census <- sort(unique(tree_data$CensusID[tree_data$CensusID < anchor_start & !is.na(tree_data$DBH) & !is.na(tree_data$TrueStemID)]))
        if (length(cand_census) > 0L) {
            new_anchor <- as.integer(max(cand_census))
            vcat(prefix, "Anchor census ", anchor_start, " has no DBH/TrueStemID; falling back to earlier anchor at census ", new_anchor)
            anchor_start <- new_anchor
            census_range <- seq.int(from = first_obs_census, to = anchor_start)
            # Exclude censuses fully emptied by MF removal
            if (isTRUE(has_mf_stash) && length(mf_emptied_censuses) > 0L) {
                census_range <- census_range[!census_range %in% mf_emptied_censuses]
            }
            n_census <- length(census_range)
            vcat(prefix, "Adjusted census range: ", paste(census_range, collapse = ", "))
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
                vcat(prefix, "No census with DBH+TrueStemID found; using provisional anchor at last observed census ", provisional_anchor)
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
                    vcat(prefix, sprintf("Assigned %d provisional anchor ID(s) at census %d", length(anchor_idx), anchor_start))
                }

                # Recompute census_range and obs_counts now that anchor_start changed
                census_range <- seq.int(from = first_obs_census, to = anchor_start)
                # Exclude censuses fully emptied by MF removal
                if (isTRUE(has_mf_stash) && length(mf_emptied_censuses) > 0L) {
                    census_range <- census_range[!census_range %in% mf_emptied_censuses]
                }
                n_census <- length(census_range)
                vcat(prefix, "Adjusted census range: ", paste(census_range, collapse = ", "))
                obs_counts <- vapply(
                    census_range,
                    function(cc) nrow(tree_data[CensusID == cc & !is.na(DBH)]),
                    integer(1L)
                )
                max_obs <- if (length(obs_counts) > 0L) max(obs_counts) else 0L
            } else {
                vcat(prefix, "No usable anchor found (no DBH or TrueStemID available); falling back to igraph matcher")
                fallback_reason <- "anchor_missing_truestem"
                K_used <- as.integer(min(max_obs, max_tracks))
                n_states_by_census <- vapply(obs_counts, function(n_obs) count_injective_states(K_used, n_obs), numeric(1L))
                tree_data[, `:=`(
                    DP_KUsed = K_used,
                    DP_MaxStatesPerCensus = max(n_states_by_census, na.rm = TRUE),
                    DP_MaxStatesCensusID = as.integer(census_range[which.max(n_states_by_census)])
                )]
                out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, anchor_start)
                out <- ensure_posterior_columns(out)
                if (!("DP_FallbackReason" %in% names(out))) out[, DP_FallbackReason := NA_character_]
                out[, DP_FallbackReason := fallback_reason]
                attr(out, "DP_PruneInfo") <- prune_stats
                return(finalize_out(out))
            }
        }
    }

    anchor_obs <- tree_data[CensusID == anchor_start & !is.na(DBH)]
    if (nrow(anchor_obs) == 0L) {
        vcat(prefix, "Anchor census ", anchor_start, " has no DBH observations; falling back to igraph matcher")
        fallback_reason <- "anchor_missing_obs"
        K_used <- as.integer(min(max_obs, max_tracks))
        n_states_by_census <- vapply(obs_counts, function(n_obs) count_injective_states(K_used, n_obs), numeric(1L))
        tree_data[, `:=`(
            DP_KUsed = K_used,
            DP_MaxStatesPerCensus = max(n_states_by_census, na.rm = TRUE),
            DP_MaxStatesCensusID = as.integer(census_range[which.max(n_states_by_census)])
        )]
        out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, anchor_start)
        out <- ensure_posterior_columns(out)
        if (!("DP_FallbackReason" %in% names(out))) out[, DP_FallbackReason := NA_character_]
        out[, DP_FallbackReason := fallback_reason]
        attr(out, "DP_PruneInfo") <- prune_stats
        return(finalize_out(out))
    }

    # If anchor_obs exists but some TrueStemID are missing, optionally allow a provisional DP anchor
    if (any(is.na(anchor_obs$TrueStemID))) {
        if (isTRUE(allow_provisional_anchor)) {
            vcat(prefix, "Anchor census ", anchor_start, " has DBH but no TrueStemID; assigning provisional IDs")
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
            vcat(prefix, "No usable anchor (TrueStemID missing and provisional anchoring disabled); falling back to igraph matcher")
            fallback_reason <- "anchor_missing_truestem_prov_disabled"
            K_used <- as.integer(min(max_obs, max_tracks))
            n_states_by_census <- vapply(obs_counts, function(n_obs) count_injective_states(K_used, n_obs), numeric(1L))
            tree_data[, `:=`(
                DP_KUsed = K_used,
                DP_MaxStatesPerCensus = max(n_states_by_census, na.rm = TRUE),
                DP_MaxStatesCensusID = as.integer(census_range[which.max(n_states_by_census)])
            )]
            out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, anchor_start)
            out <- ensure_posterior_columns(out)
            if (!("DP_FallbackReason" %in% names(out))) out[, DP_FallbackReason := NA_character_]
            out[, DP_FallbackReason := fallback_reason]
            attr(out, "DP_PruneInfo") <- prune_stats
            return(finalize_out(out))
        }
    }
    anchor_ids <- sort(unique(anchor_obs$TrueStemID))
    anchor_ids <- anchor_ids[!is.na(anchor_ids)]

    if (length(anchor_ids) == 0L) {
        vcat(prefix, "No valid anchor IDs at anchor census; falling back to igraph matcher")
        fallback_reason <- "anchor_ids_missing"
        K_used <- as.integer(min(max_obs, max_tracks))
        n_states_by_census <- vapply(obs_counts, function(n_obs) count_injective_states(K_used, n_obs), numeric(1L))
        tree_data[, `:=`(
            DP_KUsed = K_used,
            DP_MaxStatesPerCensus = max(n_states_by_census, na.rm = TRUE),
            DP_MaxStatesCensusID = as.integer(census_range[which.max(n_states_by_census)])
        )]
        out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, anchor_start)
        out <- ensure_posterior_columns(out)
        if (!("DP_FallbackReason" %in% names(out))) out[, DP_FallbackReason := NA_character_]
        out[, DP_FallbackReason := fallback_reason]
        return(finalize_out(out))
    }

    # Choose K (tracks) same logic as the MAP DP
    births_needed <- if (length(obs_counts) >= 2L) sum(pmax(0L, diff(obs_counts))) else 0L
    K_from_counts <- as.integer(if (length(obs_counts) > 0L) obs_counts[1L] + births_needed else 0L)
    K_base <- max(length(anchor_ids), max_obs, K_from_counts)

    # ---- Resprout barrier: increase K for resprout observations ----
    # Each resprout (R|RP|RF|RT|QR code with non-NA DBH) forces the track into
    # phase 0 at the preceding census, effectively creating a new identity that
    # needs its own track slot.  Add one extra track per resprout observation.
    resprout_regex <- "\\b(R|RP|RF|RT|QR)\\b"
    # Pre-compute per-census row indices and resprout flags ONCE here.
    # This single pass is reused in both the resprout count below and the
    # state enumeration loop after K is determined, eliminating duplicate
    # [.data.table subset calls from those two loops.
    .has_tsm <- "ListOfTSM" %in% names(tree_data)
    .obs_row_idx_pre <- vector("list", n_census)
    .is_resprout_pre <- vector("list", n_census)
    for (.p0 in seq_len(n_census)) {
        .idx0 <- tree_data[CensusID == census_range[.p0] & !is.na(DBH), which = TRUE]
        .obs_row_idx_pre[[.p0]] <- .idx0
        if (length(.idx0) > 0L && .has_tsm) {
            .tsm0 <- tree_data$ListOfTSM[.idx0]
            .is_resprout_pre[[.p0]] <- !is.na(.tsm0) & grepl(resprout_regex, .tsm0)
        } else {
            .is_resprout_pre[[.p0]] <- rep(FALSE, length(.idx0))
        }
    }
    n_resprout_total <- sum(vapply(.is_resprout_pre, function(x) sum(x), integer(1L)))
    if (n_resprout_total > 0L) {
        K_base <- K_base + n_resprout_total
        vcat(prefix, "Resprout barrier: adding ", n_resprout_total, " extra track(s) for resprout observations")
    }

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
            vcat(prefix, "Slack tracks not granted: no anchor DBH is small enough to be recruitable (eps=", eps, ")")
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
        vcat(prefix, "K=", K, " tracks is fewer than max observed stems (", max(obs_counts), "); falling back to igraph matcher")
        fallback_reason <- "K_too_small"
        out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, anchor_start)
        out <- ensure_posterior_columns(out)
        if (!("DP_FallbackReason" %in% names(out))) out[, DP_FallbackReason := NA_character_]
        out[, DP_FallbackReason := fallback_reason]
        attr(out, "DP_PruneInfo") <- prune_stats
        return(finalize_out(out))
    }

    vcat(prefix, "Using K=", K, " identity tracks; worst-case states per census: ", format(max(n_states_by_census, na.rm = TRUE), big.mark = ","))

    n_extra <- K - length(anchor_ids)
    current_max <- suppressWarnings(max(tree_data$TrueStemID, na.rm = TRUE))
    if (!is.finite(current_max)) current_max <- 0
    track_ids <- c(anchor_ids, if (n_extra > 0L) seq.int(from = current_max + 1L, length.out = n_extra) else integer(0))

    # Pre-enumerate assignment states (injective obs->track) for each census in census_range
    obs_dbh <- vector("list", n_census)
    obs_row_idx <- vector("list", n_census)    # precomputed row indices per census (avoids repeated [.data.table)
    is_resprout_obs <- vector("list", n_census)
    state_mats <- vector("list", n_census)
    state_keys <- vector("list", n_census)
    for (p in seq_len(n_census)) {
        cc <- census_range[p]
        idx <- .obs_row_idx_pre[[p]]    # reuse pre-computed indices (eliminates [.data.table call)
        obs_row_idx[[p]] <- idx
        obs_dbh[[p]] <- tree_data$DBH[idx]
        n_obs <- length(obs_dbh[[p]])
        is_resprout_obs[[p]] <- .is_resprout_pre[[p]]  # reuse pre-computed flags (eliminates grepl call)
        mat <- enumerate_states_injective(K, n_obs, max_states = max_states)

        if (is.null(mat)) {
            vcat(prefix, "Too many states at census ", cc, " (", n_obs, " stems, exceeds max_states=", max_states, "); falling back to igraph matcher")
            fallback_reason <- "enum_exceeded"
            out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, anchor_start)
            out <- ensure_posterior_columns(out)
            if (!("DP_FallbackReason" %in% names(out))) out[, DP_FallbackReason := NA_character_]
            out[, DP_FallbackReason := fallback_reason]
            attr(out, "DP_PruneInfo") <- prune_stats
            return(finalize_out(out))
        }
        state_mats[[p]] <- mat
        # Vectorized state key computation: avoids per-row paste() calls from apply().
        # do.call(paste, list_of_columns) produces comma-separated strings in R's C layer.
        state_keys[[p]] <- if (ncol(mat) == 0L) {
            rep("", nrow(mat))
        } else {
            do.call(paste, c(lapply(seq_len(ncol(mat)), function(j) mat[, j]), list(sep = ",")))
        }

        vcat(prefix, "Census ", cc, ": ", n_obs, " observed stem(s), ", nrow(mat), " assignment states")
    }

    # Anchor state assignment (pins endpoint)
    anchor_obs_ordered <- tree_data[CensusID == anchor_start & !is.na(DBH)]
    anchor_track_idx <- match(anchor_obs_ordered$TrueStemID, track_ids)

    if (any(is.na(anchor_track_idx))) {
        vcat(prefix, "Anchor IDs not found in track list; falling back to igraph matcher")
        fallback_reason <- "anchor_truestem_not_found"
        out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, anchor_start)
        out <- ensure_posterior_columns(out)
        if (!("DP_FallbackReason" %in% names(out))) out[, DP_FallbackReason := NA_character_]
        out[, DP_FallbackReason := fallback_reason]
        attr(out, "DP_PruneInfo") <- prune_stats
        return(finalize_out(out))
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

    # Palm override: DBH is stable in palms — collapse prune window to very tight bounds
    if (isTRUE(is_palm)) {
        eff_min_grow <- palm_prune_min_growth
        eff_max_grow <- palm_prune_max_growth
    }

    prune_stats$eff_min_growth <- as.numeric(eff_min_grow)
    prune_stats$eff_max_growth <- as.numeric(eff_max_grow)
    prune_stats$eff_recruit_max <- as.numeric(eff_recruit_max)

    vcat(prefix, "Pruning bounds: growth [", eff_min_grow, ", ", eff_max_grow, "] cm/yr, max recruit DBH=", eff_recruit_max, " cm")

    # Precompute track-wise DBH matrix for each census (rows = states, cols = tracks).
    # Vectorized: one matrix indexing operation per census instead of a per-state loop.
    track_dbh_by_state <- vector("list", n_census)
    for (p in seq_len(n_census)) {
        mat      <- state_mats[[p]]
        n_states <- nrow(mat)
        n_obs_p  <- ncol(mat)
        tdbh_mat <- matrix(NA_real_, nrow = n_states, ncol = K)
        if (n_obs_p > 0L && n_states > 0L) {
            row_idx_v <- rep(seq_len(n_states), times = n_obs_p)
            col_idx_v <- as.vector(mat)            # column-major track indices
            dbh_vals  <- rep(obs_dbh[[p]], each = n_states)
            tdbh_mat[cbind(row_idx_v, col_idx_v)] <- dbh_vals
        }
        track_dbh_by_state[[p]] <- tdbh_mat        # n_states × K matrix
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
    vcat(prefix, "Starting backward pass (anchor -> earliest census) ...")

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
            vcat(prefix, "  Backward step ", (anchor_pos - p), "/", (anchor_pos - 1L), ": census ", cc, " -> ", next_cc, " (", n_states_cc, " candidate states, ", length(keys_full[[p + 1L]]), " target states)")

        next_keys <- keys_full[[p + 1L]]
        n_next <- length(next_keys)
        if (n_next == 0L) {
            vcat(prefix, "  No reachable states at census ", next_cc, "; falling back to igraph matcher")
            fallback_reason <- "no_reachable_next_states"
            out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, anchor_start)
            out <- ensure_posterior_columns(out)
            if (!("DP_FallbackReason" %in% names(out))) out[, DP_FallbackReason := NA_character_]
            out[, DP_FallbackReason := fallback_reason]
            attr(out, "DP_PruneInfo") <- prune_stats
            return(finalize_out(out))
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
            fallback_reason <- "next_assign_row_mismatch"
            out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, anchor_start)
            out <- ensure_posterior_columns(out)
            if (!("DP_FallbackReason" %in% names(out))) out[, DP_FallbackReason := NA_character_]
            out[, DP_FallbackReason := fallback_reason]
            attr(out, "DP_PruneInfo") <- prune_stats
            return(finalize_out(out))
        }

        # Pre-decode next phase vectors
        phase_tp1_by_next <- vector("list", n_next)
        for (j in seq_len(n_next)) {
            phase_tp1_by_next[[j]] <- decode_full_key(next_keys[[j]])$phase
        }

        # Interval (years) between cc and next_cc (pair_interval precomputed above the loop)
        val_next <- pair_interval[[as.character(next_cc)]]
        val_cc   <- pair_interval[[as.character(cc)]]
        if (is.null(val_next) || length(val_next) == 0 || is.null(val_cc) || length(val_cc) == 0) {
            interval_val <- NA_real_
        } else {
            interval_val <- (as.numeric(val_next) - as.numeric(val_cc)) / 365.25
            if (!is.finite(interval_val) || interval_val <= 0) {
                interval_val <- NA_real_
            }
        }
        vcat(prefix, "  Interval: ", sprintf("%.2f", interval_val), " years (census ", cc, " -> ", next_cc, ")", sep = "")

        # -----------------------------------------------------------------------
        # Batch feasibility check in C++: replaces the O(n_cc × n_next × K) R
        # inner loop.  derive_phase_prev_batch_rcpp checks phase-transition
        # constraints and hard growth-rate pruning for every (i, j) pair at once.
        # -----------------------------------------------------------------------

        # Build input matrices (rows = states/next-states, cols = tracks K)
        tdbh0_mat     <- track_dbh_by_state[[p]]                                             # already n_cc   × K
        tdbh1_mat     <- track_dbh_by_state[[p + 1L]][next_assign_row_idx, , drop = FALSE]  # n_next × K
        phase_tp1_mat <- do.call(rbind, phase_tp1_by_next)                                  # n_next × K; integer

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

        feasible_result <- derive_phase_prev_batch_rcpp(
            tdbh0_mat      = tdbh0_mat,
            tdbh1_mat      = tdbh1_mat,
            phase_tp1_mat  = phase_tp1_mat,
            resprout_mat   = resp_mat,
            prune_hard     = isTRUE(prune_hard),
            interval_val   = if (is.finite(interval_val)) interval_val else NaN,
            eff_min_grow   = eff_min_grow,
            eff_max_grow   = eff_max_grow,
            eff_recruit_max = if (is.finite(eff_recruit_max)) eff_recruit_max else Inf
        )

        fe_from  <- feasible_result$from_i   # 1-based current assignment indices
        fe_to    <- feasible_result$to_j     # 1-based next full-state indices
        fe_phase <- feasible_result$phase_t  # n_feasible × K integer matrix
        n_feasible <- length(fe_from)

        # Update prune diagnostics
        if (isTRUE(prune_hard)) {
            n_examined   <- n_states_cc * n_next
            n_infeasible <- n_examined - n_feasible
            prune_stats$total_examined <- prune_stats$total_examined + n_examined
            prune_stats$total_pruned   <- prune_stats$total_pruned   + n_infeasible
            prune_stats$per_census[[as.character(cc)]] <-
                prune_stats$per_census[[as.character(cc)]] + n_infeasible
        }

        # Compute full-state key strings for each feasible pair using vectorized ops.
        if (n_feasible > 0L) {
            fe_assign_keys <- state_keys[[p]][fe_from]  # assign key portion (precomputed)
            # Vectorized phase key encoding: avoids per-row apply + rawToChar overhead
            .phase_chars  <- c("0", "1", "2")
            fe_phase_keys <- do.call(paste0, lapply(seq_len(K), function(k) .phase_chars[fe_phase[, k] + 1L]))
            fe_full_keys  <- paste0(fe_assign_keys, "|", fe_phase_keys)
        }

        # Dynamic creation of current full-states
        key_to_idx <- new.env(parent = emptyenv())
        curr_keys_list <- list()
        curr_assign_list <- list()
        curr_logB <- numeric(0)
        curr_vit <- numeric(0)
        curr_ptr <- integer(0)

        # Preallocate edge arrays (at most n_feasible edges)
        from_idx <- integer(n_feasible)
        to_idx   <- integer(n_feasible)
        logw     <- numeric(n_feasible)
        used_edges <- 0L

        # Process feasible pairs grouped by from_i (i-order is guaranteed by C++ loop order).
        # split() creates one group per unique from_i — only i values with ≥1 feasible j appear.
        if (n_feasible > 0L) {
            by_i <- split(seq_len(n_feasible), fe_from)

            for (i_key in names(by_i)) {
                i_val   <- as.integer(i_key)
                rows    <- by_i[[i_key]]
                tdbh0   <- track_dbh_by_state[[p]][i_val, ]
                assign0 <- mat_cc[i_val, ]
                j_vals  <- fe_to[rows]
                f_keys  <- fe_full_keys[rows]
                f_tdbh1 <- lapply(j_vals, function(jj) tdbh1_mat[jj, ])

                if (verbose) t_tc0 <- tic()
                c_trans_vec <- transition_cost_tracks_bio_batch_rcpp(
                    track_dbh_t   = tdbh0,
                    track_dbh_tp1 = f_tdbh1,
                    interval_years = interval_val,
                    # --- growth model ---
                    mu_const = Bio_Mu_Growth_unit,
                    mu_gamma = Bio_Gamma_Growth_unit,
                    sigma0   = Bio_Sigma0_unit,
                    sigma1   = Bio_Sigma1_unit,
                    max_shrink      = Bio_max_shrink_unit,
                    k_shrink        = Bio_k_shrink_unit,
                    max_growth      = Bio_max_growth_unit,
                    max_growth_soft = Bio_max_growth_soft_unit,
                    k_growth        = Bio_k_growth_unit,
                    # --- measurement error (optional) ---
                    use_measurement_error = use_measurement_error,
                    meas_sd1_a = meas_sd1_a,
                    meas_sd1_b = meas_sd1_b,
                    meas_sd2   = meas_sd2,
                    meas_p_big = meas_p_big,
                                    # --- mortality model ---
                    h0             = Bio_H0_Mortality,
                    beta           = Bio_Beta_Mortality,
                                    # --- recruitment model ---
                    recruit_meanlog = Bio_Recruit_Meanlog_unit,
                    recruit_sdlog   = Bio_Recruit_Sdlog_unit,
                    recruit_max_dbh = Bio_Recruit_MaxDBH_unit,
                    recruit_lambda  = Bio_Recruitment_lambda,
                    eps_tiebreak    = eps_tiebreak
                )
                transition_cost_calls <- transition_cost_calls + 1L
                if (verbose) transition_cost_time <- transition_cost_time + (tic() - t_tc0)

                # Step 1: sequential state registration (new keys need sequential idx)
                n_e <- length(rows)
                c_trans_num <- unlist(c_trans_vec, use.names = FALSE)
                edge_idx <- integer(n_e)
                for (e in seq_len(n_e)) {
                    curr_key <- f_keys[[e]]
                    idx <- key_to_idx[[curr_key]]
                    if (is.null(idx)) {
                        idx <- length(curr_keys_list) + 1L
                        key_to_idx[[curr_key]] <- idx
                        curr_keys_list[[idx]]  <- curr_key
                        curr_assign_list[[idx]] <- as.integer(assign0)
                        curr_logB[idx] <- -Inf
                        curr_vit[idx]  <- Inf
                        curr_ptr[idx]  <- NA_integer_
                    }
                    edge_idx[[e]] <- idx
                }
                # Step 2: vectorized cost accumulation — eliminates per-edge log_add_exp R calls
                .logw_v     <- -c_trans_num / temperature
                .cand_vit_v <- c_trans_num + vit_next[j_vals]
                .cand_log_v <- .logw_v + logB_next[j_vals]
                if (n_e == 1L) {
                    # Single-edge fast path
                    .ix <- edge_idx[1L]
                    if (!is.finite(curr_vit[.ix]) || .cand_vit_v[1L] < curr_vit[.ix]) {
                        curr_vit[.ix] <- .cand_vit_v[1L]; curr_ptr[.ix] <- j_vals[1L]
                    }
                    .lb <- curr_logB[.ix]
                    curr_logB[.ix] <- if (!is.finite(.lb)) .cand_log_v[1L] else {
                        .m <- max(.lb, .cand_log_v[1L])
                        .m + log(exp(.lb - .m) + exp(.cand_log_v[1L] - .m))
                    }
                    used_edges <- used_edges + 1L
                    from_idx[[used_edges]] <- .ix; to_idx[[used_edges]] <- j_vals[1L]
                    logw[[used_edges]] <- .logw_v[1L]
                } else if (!anyDuplicated(edge_idx)) {
                    # No duplicate indices: fully vectorized Viterbi + logB updates
                    .vit_old <- curr_vit[edge_idx]
                    .upd <- !is.finite(.vit_old) | (.cand_vit_v < .vit_old)
                    if (any(.upd)) {
                        curr_vit[edge_idx[.upd]] <- .cand_vit_v[.upd]
                        curr_ptr[edge_idx[.upd]] <- j_vals[.upd]
                    }
                    .lb_old <- curr_logB[edge_idx]
                    .fin    <- is.finite(.lb_old)
                    if (any(.fin)) {
                        .m <- pmax(.lb_old[.fin], .cand_log_v[.fin])
                        curr_logB[edge_idx[.fin]] <- .m + log(
                            exp(.lb_old[.fin] - .m) + exp(.cand_log_v[.fin] - .m))
                    }
                    if (any(!.fin)) curr_logB[edge_idx[!.fin]] <- .cand_log_v[!.fin]
                    .e_rng <- seq.int(used_edges + 1L, used_edges + n_e)
                    from_idx[.e_rng] <- edge_idx; to_idx[.e_rng] <- j_vals; logw[.e_rng] <- .logw_v
                    used_edges <- used_edges + n_e
                } else {
                    # Fallback: duplicate edge_idx values (rare — two j's yield same full-key)
                    for (e in seq_len(n_e)) {
                        .ix <- edge_idx[e]
                        if (!is.finite(curr_vit[.ix]) || .cand_vit_v[e] < curr_vit[.ix]) {
                            curr_vit[.ix] <- .cand_vit_v[e]; curr_ptr[.ix] <- j_vals[e]
                        }
                        .lb <- curr_logB[.ix]
                        curr_logB[.ix] <- if (!is.finite(.lb)) .cand_log_v[e] else {
                            .m <- max(.lb, .cand_log_v[e])
                            .m + log(exp(.lb - .m) + exp(.cand_log_v[e] - .m))
                        }
                        used_edges <- used_edges + 1L
                        from_idx[[used_edges]] <- .ix; to_idx[[used_edges]] <- j_vals[e]
                        logw[[used_edges]] <- .logw_v[e]
                    }
                }
            }
        }

        if (length(curr_keys_list) == 0L) {
            vcat(prefix, "  No valid states produced at census ", cc, "; falling back to igraph matcher")
            fallback_reason <- "no_states_produced"
            out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, anchor_start)
            out <- ensure_posterior_columns(out)
            if (!("DP_FallbackReason" %in% names(out))) out[, DP_FallbackReason := NA_character_]
            out[, DP_FallbackReason := fallback_reason]
            attr(out, "DP_PruneInfo") <- prune_stats
            return(finalize_out(out))
        }

        keys_full[[p]] <- unlist(curr_keys_list, use.names = FALSE)
        assign_full[[p]] <- curr_assign_list
        logB[[p]] <- curr_logB
        vit_cost[[p]] <- curr_vit
        vit_ptr[[p]] <- curr_ptr

        if (used_edges == 0L) {
            vcat(prefix, "  No feasible transitions at census ", cc, "; falling back to igraph matcher")
            fallback_reason <- "no_feasible_edges"
            out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, anchor_start)
            out <- ensure_posterior_columns(out)
            if (!("DP_FallbackReason" %in% names(out))) out[, DP_FallbackReason := NA_character_]
            out[, DP_FallbackReason := fallback_reason]
            attr(out, "DP_PruneInfo") <- prune_stats
            return(finalize_out(out))
        }
        edges[[p]] <- data.table::data.table(
            from_idx = from_idx[seq_len(used_edges)],
            to_idx   = to_idx[seq_len(used_edges)],
            logw     = logw[seq_len(used_edges)]
        )

        vcat(prefix, "  Backward step ", (anchor_pos - p), "/", (anchor_pos - 1L), " done: census ", cc, " has ", length(keys_full[[p]]), " reachable states, ", used_edges, " transitions", if (verbose) paste0(" (", sprintf("%.2fs", tic() - t_cc0), ")") else "")
    }
    }

    # -----------------
    # Decode MAP path
    # -----------------
    vcat(prefix, "Decoding best reconstruction path ...")
    start_idx <- which.min(vit_cost[[1L]])
    if (length(start_idx) == 0L || !is.finite(vit_cost[[1L]][start_idx])) {
        fallback_reason <- "decode_failure"
        out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, anchor_start)
        out <- ensure_posterior_columns(out)
        if (!("DP_FallbackReason" %in% names(out))) out[, DP_FallbackReason := NA_character_]
        out[, DP_FallbackReason := fallback_reason]
        attr(out, "DP_PruneInfo") <- prune_stats
        return(finalize_out(out))
    }
    map_idx <- integer(n_census)
    map_idx[1L] <- start_idx
    for (p in seq_len(n_census - 1L)) {
        nxt <- vit_ptr[[p]][map_idx[p]]
        if (!is.finite(nxt) || is.na(nxt) || nxt < 1L) {
            fallback_reason <- "viterbi_decode_failure"
            out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, anchor_start)
            out <- ensure_posterior_columns(out)
            if (!("DP_FallbackReason" %in% names(out))) out[, DP_FallbackReason := NA_character_]
            out[, DP_FallbackReason := fallback_reason]
            attr(out, "DP_PruneInfo") <- prune_stats
            return(finalize_out(out))
        }
        map_idx[p + 1L] <- nxt
    }

    for (p in seq_len(n_census)) {
        cc <- census_range[p]
        obs_idx <- obs_row_idx[[p]]
        if (length(obs_idx) == 0L) next
        sv <- assign_full[[p]][[map_idx[p]]]
        if (length(sv) != length(obs_idx)) {
            fallback_reason <- "assign_mismatch"
            out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, anchor_start)
            out <- ensure_posterior_columns(out)
            if (!("DP_FallbackReason" %in% names(out))) out[, DP_FallbackReason := NA_character_]
            out[, DP_FallbackReason := fallback_reason]
            attr(out, "DP_PruneInfo") <- prune_stats
            return(finalize_out(out))
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
    vcat(prefix, "Starting forward pass (earliest census -> anchor) for uncertainty quantification ...")
    # Start distribution: uniform over all reachable states at census 1.
    logalpha <- vector("list", n_census)
    logalpha[[1L]] <- rep.int(0, length(keys_full[[1L]]))

    for (p in seq_len(n_census - 1L)) {
        ed <- edges[[p]]
        if (is.null(ed) || nrow(ed) == 0L) {
            fallback_reason <- "forward_edges_missing"
            out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, anchor_start)
            out <- ensure_posterior_columns(out)
            if (!("DP_FallbackReason" %in% names(out))) out[, DP_FallbackReason := NA_character_]
            out[, DP_FallbackReason := fallback_reason]
            attr(out, "DP_PruneInfo") <- prune_stats
            return(finalize_out(out))
        }
        la_from <- logalpha[[p]][ed$from_idx]
        vals <- la_from + ed$logw
        dt <- data.table::data.table(to_idx = ed$to_idx, v = vals)
        dt <- dt[is.finite(v)]
        if (nrow(dt) == 0L) {
            fallback_reason <- "forward_no_alpha"
            out <- match_stems_optimal_backward(tree_data, min_growth, max_growth, anchor_start)
            out <- ensure_posterior_columns(out)
            if (!("DP_FallbackReason" %in% names(out))) out[, DP_FallbackReason := NA_character_]
            out[, DP_FallbackReason := fallback_reason]
            attr(out, "DP_PruneInfo") <- prune_stats
            return(finalize_out(out))
        }
        la_next_dt <- dt[, .(logalpha = log_sum_exp(v)), by = to_idx]
        la_next <- rep.int(-Inf, length(keys_full[[p + 1L]]))
        la_next[la_next_dt$to_idx] <- la_next_dt$logalpha
        logalpha[[p + 1L]] <- la_next

        vcat(prefix, "  Forward step ", p, "/", (n_census - 1L), ": census ", census_range[p + 1L], " — ", sum(is.finite(la_next)), "/", length(la_next), " states reachable")
    }

    # Partition function Z = total weight of all paths ending at the fixed anchor state.
    # At anchor_start there is exactly one state.
    logZ <- logalpha[[anchor_pos]][1L]
    if (!is.finite(logZ)) {
        # Fallback: compute from backward at census 1.
        logZ <- log_sum_exp(logB[[1L]])
    }

    vcat(prefix, "Log-partition function (logZ) = ", sprintf("%.3f", logZ), " (normalisation constant for posterior probabilities)")

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
        vcat(prefix, sprintf("Posterior sampling: generated %d samples (%d unique paths)", posterior_samples, n_unique_paths))
        if (uniqueN(sample_logp$logp) == 1L) {
            vcat(prefix, sprintf("Warning: all %d posterior samples have identical log-probability (%.3f); the posterior may be concentrated on a single reconstruction", posterior_samples, sample_logp$logp[1]))
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
        # Enforce ObsRowID-based reconstructions only (legacy CensusID encoding removed).
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
        post_rows <- original_tree_data[CensusID > anchor_requested]
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
        # Propagate "given" ReconstructedStemID for post-anchor rows when DP actually used those IDs
        used_ids <- unique(processed$ReconstructedStemID[!is.na(processed$ReconstructedStemID)])
        if (length(used_ids) > 0L) {
            post_rows[!is.na(TrueStemID) & !is.na(DBH) & (TrueStemID %in% used_ids), `:=`(
                ReconstructedStemID = as.integer(TrueStemID),
                ReconstructionMethod = "given"
            )]
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
        vcat(prefix, "Preserved ", nrow(post_rows), " post-anchor row(s) in output (unchanged)")
    }

    # Attach pruning diagnostics
    attr(tree_data, "DP_PruneInfo") <- prune_stats
    if (prune_stats$total_examined > 0L) {
        vcat(prefix, "Pruning summary: removed ", prune_stats$total_pruned, " of ", prune_stats$total_examined, " candidate transitions")
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

    vcat(prefix, "Finished. Total time: ", sprintf("%.2fs", tic() - t_start))
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
