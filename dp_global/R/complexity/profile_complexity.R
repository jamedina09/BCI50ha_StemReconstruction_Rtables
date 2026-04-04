## ── Complexity Profiler ──────────────────────────────────────────────────
##
## Identifies the most complex tags in a dataset: state counts,
## transition counts (edges), and estimated wall-clock time.
##
## Usage:
##   Rscript dp_global/R/complexity/profile_complexity.R                     # BCI default
##   Rscript dp_global/R/complexity/profile_complexity.R --DATA=simulated    # simulated data
##   Rscript dp_global/R/complexity/profile_complexity.R --TOP_N=30          # show top 30
##   Rscript dp_global/R/complexity/profile_complexity.R --DP_REFINE=TRUE    # refine top tags with DP-constrained enumeration
##   Rscript dp_global/R/complexity/profile_complexity.R --RUN_TAG=002216    # time a specific tag through actual DP
## ────────────────────────────────────────────────────────────────────────

library(data.table)
library(here)

# ── CLI args ─────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default) {
    pat <- paste0("^--", name, "=(.*)$")
    m   <- grep(pat, args, value = TRUE)
    if (length(m) == 0L) return(default)
    sub(pat, "\\1", m[1L])
}

DATA_SOURCE   <- tolower(get_arg("DATA", "bci"))
TOP_N         <- as.integer(get_arg("TOP_N", "20"))
DP_MAX_STATES <- as.integer(get_arg("DP_MAX_STATES", "40000"))
DP_REFINE     <- as.logical(get_arg("DP_REFINE", "FALSE"))
TOP_N_DP      <- as.integer(get_arg("TOP_N_DP", "20"))
RUN_TAG       <- get_arg("RUN_TAG", "")

# ── Load data ────────────────────────────────────────────────────────────
if (DATA_SOURCE == "simulated") {
    INPUT_FILE <- here("data_simulation", "data", "simulated_data_1.csv")
    cat("Loading simulated data:", INPUT_FILE, "\n")
    x <- fread(INPUT_FILE)
} else {
    INPUT_FILE <- here("bci_data", "bci_multistem_xrun_debug.rds")
    cat("Loading BCI data:", INPUT_FILE, "\n")
    x <- readRDS(INPUT_FILE)
    setDT(x)
}
cat(sprintf("  %s tags, %s rows\n\n",
            format(uniqueN(x$Tag), big.mark = ","),
            format(nrow(x), big.mark = ",")))

# ── Stem ID column ──────────────────────────────────────────────────────
stem_col <- intersect(c("OriginalStemID", "StemID"), names(x))[1L]
if (is.na(stem_col)) stop("No StemID/OriginalStemID column found")

species_col <- intersect(c("Species", "species"), names(x))[1L]
if (is.na(species_col)) { x[, Species := "unknown"]; species_col <- "Species" }

# ══════════════════════════════════════════════════════════════════════════
# STEP 1: Fast vectorized scan (all tags)
# ══════════════════════════════════════════════════════════════════════════
cat("=== Step 1: Fast vectorized complexity scan ===\n")

setkey(x, Tag, CensusID)

# Obs counts per (Tag, Census)
obs_sum <- x[!is.na(DBH), .N, by = .(Tag, CensusID)]
setkey(obs_sum, Tag, CensusID)

# Per-tag: effective anchor, number of unique stems
tag_meta <- obs_sum[, .(last_obs = max(CensusID)), by = Tag]
tag_meta[, eff_anchor := as.integer(pmin(last_obs, 7L))]

# Unique stems per tag
stem_counts <- x[, .(n_stem = uniqueN(get(stem_col))), by = Tag]
tag_meta <- stem_counts[tag_meta, on = "Tag"]

# Only multi-stem tags
tag_meta <- tag_meta[n_stem > 1L]
setkey(tag_meta, Tag)

# Obs up to anchor
obs_pre <- obs_sum[tag_meta[, .(Tag, eff_anchor)], on = "Tag"][CensusID <= eff_anchor]
setorder(obs_pre, Tag, CensusID)

# Per-tag census stats
tag_stats <- obs_pre[, {
    n <- N
    list(
        n_censuses    = .N,
        max_obs       = max(n),
        K_from_counts = n[1L] + sum(pmax(0L, diff(n)))
    )
}, by = Tag]

# Anchor ids (provisional anchor: if no TrueStemID, count observed stems)
data_anch <- x[!is.na(DBH)][
    tag_meta[, .(Tag, CensusID = eff_anchor)], on = .(Tag, CensusID), nomatch = NULL]
anchor_stats <- data_anch[, {
    real_ids <- unique(na.omit(TrueStemID))
    n_real   <- length(real_ids)
    n_obs    <- .N
    # Mirrors allow_provisional_anchor=TRUE: if no TrueStemID, use n_obs as anchor count
    .(n_anchor_ids = as.integer(if (n_real > 0L) n_real else n_obs))
}, by = Tag]

# Resprouts
if ("ListOfTSM" %in% names(x)) {
    resp_sum <- x[!is.na(DBH) & !is.na(ListOfTSM) &
                  grepl("\\b(R|RP|RF|RT|QR|OR)\\b", ListOfTSM, perl = TRUE),
                  .N, by = .(Tag, CensusID)]
    resp_pre <- resp_sum[tag_meta[, .(Tag, eff_anchor)], on = "Tag"][CensusID <= eff_anchor]
    resp_per_tag <- resp_pre[, .(n_resprout = sum(N)), by = Tag]
} else {
    resp_per_tag <- data.table(Tag = character(0), n_resprout = integer(0))
}

# Species
sp_tab <- x[, .(Species = as.character(get(species_col)[1L])), by = Tag]

# Assemble
out <- tag_meta[tag_stats, on = "Tag"]
out <- anchor_stats[out, on = "Tag"]
out[is.na(n_anchor_ids), n_anchor_ids := 0L]
out <- resp_per_tag[out, on = "Tag"]
out[is.na(n_resprout), n_resprout := 0L]
out <- sp_tab[out, on = "Tag"]

# Census range string
cr_str <- obs_pre[, .(census_range = paste(CensusID, collapse = ",")), by = Tag]
out <- cr_str[out, on = "Tag"]

# K computation
out[, K_base := pmax(n_anchor_ids, max_obs, K_from_counts) + n_resprout]
out[, K := as.integer(K_base + 1L)]  # +1 for slack

# State counts (vectorized via unique K,n pairs)
obs_with_K <- out[, .(Tag, K)][obs_pre, on = "Tag", nomatch = NULL]
unique_kn <- unique(obs_with_K[, .(K, N)])
unique_kn[, n_states := mapply(function(Kv, nv) {
    if (nv == 0L || Kv <= 0L) return(1.0)
    if (nv > Kv) return(0.0)
    prod(seq.int(from = as.integer(Kv), to = as.integer(Kv) - as.integer(nv) + 1L, by = -1L))
}, K, N)]
obs_with_K <- unique_kn[obs_with_K, on = .(K, N)]
setorder(obs_with_K, Tag, CensusID)

census_agg <- obs_with_K[, .(
    max_states_per_census    = max(n_states),
    total_states             = sum(n_states),
    estimated_edges_unpruned = if (.N >= 2L) sum(as.numeric(n_states[-.N]) * as.numeric(n_states[-1L])) else 0
), by = Tag]
out <- census_agg[out, on = "Tag"]

# Fallback flags
out[, estimated_fallback := (K < max_obs) |
        (max_states_per_census > as.numeric(DP_MAX_STATES))]
out[, fallback_reason := fcase(
    K < max_obs,                                      "K_too_small",
    max_states_per_census > as.numeric(DP_MAX_STATES), "enum_exceeded",
    default = NA_character_
)]

# Time estimate: linear scaling from benchmark data
# Calibrated from actual DP runs on BCI:
#   Tag 002216: 67.8M edges -> 62.6s  => ~1.08M edges/s (backward+forward)
#   Typical tags: ~30K edges -> ~0.1s  (dominated by overhead)
# Model: time = overhead + edges / throughput
predict_seconds <- function(N) {
    overhead    <- 0.03  # minimum per-tag overhead (seconds)
    throughput  <- 1.0e6 # edges per second (conservative)
    overhead + N / throughput
}
out[, predicted_seconds := predict_seconds(estimated_edges_unpruned)]

setorder(out, estimated_fallback, -max_states_per_census, -estimated_edges_unpruned)

n_total  <- nrow(out)
n_dp     <- sum(!out$estimated_fallback)
n_igraph <- sum(out$estimated_fallback)

cat(sprintf("\n=== Tag Distribution (DP_MAX_STATES = %s) ===\n",
            format(DP_MAX_STATES, big.mark = ",")))
cat(sprintf("  Total multi-stem tags:  %s\n", format(n_total, big.mark = ",")))
cat(sprintf("  Handled by DP:          %s  (%.1f%%)\n",
            format(n_dp, big.mark = ","), 100 * n_dp / n_total))
cat(sprintf("  Fallback to igraph:     %s  (%.1f%%)\n",
            format(n_igraph, big.mark = ","), 100 * n_igraph / n_total))

# Fallback breakdown
if (n_igraph > 0) {
    fb <- out[estimated_fallback == TRUE, .N, by = fallback_reason]
    cat("  Fallback reasons:\n")
    for (i in seq_len(nrow(fb)))
        cat(sprintf("    %-25s %s tags\n", fb$fallback_reason[i], format(fb$N[i], big.mark = ",")))
}

# State-space tiers
cat("\n=== State-Space Distribution (DP tags only) ===\n")
dp_tags <- out[estimated_fallback == FALSE]
breaks <- c(0, 10, 100, 1000, 5000, 10000, 40000, Inf)
labels <- c("<=10", "11-100", "101-1K", "1K-5K", "5K-10K", "10K-40K", ">40K")
dp_tags[, tier := cut(max_states_per_census, breaks = breaks, labels = labels,
                      right = TRUE, include.lowest = TRUE)]
tier_tab <- dp_tags[, .N, by = tier][order(tier)]
for (i in seq_len(nrow(tier_tab)))
    cat(sprintf("  %-12s %s tags\n", tier_tab$tier[i], format(tier_tab$N[i], big.mark = ",")))

# ══════════════════════════════════════════════════════════════════════════
# STEP 2: Top N most complex tags
# ══════════════════════════════════════════════════════════════════════════
top <- head(out[estimated_fallback == FALSE][order(-max_states_per_census, -estimated_edges_unpruned)], TOP_N)

cat(sprintf("\n=== Top %d Most Complex Tags (unconstrained) ===\n", nrow(top)))
cat(sprintf("%-10s %-10s %5s %7s %7s %15s %15s %20s %10s\n",
            "Tag", "Species", "K", "max_obs", "n_cens", "max_states", "total_states",
            "edges_unpruned", "est_time"))
cat(strrep("-", 110), "\n")
for (i in seq_len(nrow(top))) {
    r <- top[i]
    time_str <- if (r$predicted_seconds < 60) sprintf("%.1fs", r$predicted_seconds)
                else if (r$predicted_seconds < 3600) sprintf("%.1f min", r$predicted_seconds / 60)
                else sprintf("%.1f hrs", r$predicted_seconds / 3600)
    cat(sprintf("%-10s %-10s %5d %7d %7d %15s %15s %20s %10s\n",
                r$Tag, substr(r$Species, 1, 10), r$K, r$max_obs, r$n_censuses,
                format(r$max_states_per_census, big.mark = ","),
                format(r$total_states, big.mark = ","),
                format(r$estimated_edges_unpruned, big.mark = ","),
                time_str))
}

# ══════════════════════════════════════════════════════════════════════════
# STEP 3 (optional): DP-constrained refinement of top tags
# ══════════════════════════════════════════════════════════════════════════
if (isTRUE(DP_REFINE)) {
    cat(sprintf("\n=== Step 3: DP-constrained refinement (top %d tags) ===\n", min(TOP_N_DP, nrow(top))))
    source(here("dp_global", "R", "dp_global_states.R"))

    refine_tags <- head(top$Tag, min(TOP_N_DP, nrow(top)))

    for (tg in refine_tags) {
        td <- x[Tag == tg]
        r  <- out[Tag == tg]

        # Reconstruct census_range, K, anchor, track_ids
        eff_anchor <- r$eff_anchor
        K <- r$K
        obs_census <- sort(unique(td[!is.na(DBH) & CensusID <= eff_anchor, CensusID]))
        census_range <- obs_census
        if (eff_anchor > max(obs_census)) census_range <- c(census_range, eff_anchor)
        n_census <- length(census_range)

        anchor_obs <- td[CensusID == eff_anchor & !is.na(DBH)]
        anchor_ids <- sort(unique(na.omit(anchor_obs$TrueStemID)))
        if (length(anchor_ids) == 0L) { cat(sprintf("  %s: no anchor IDs, skip\n", tg)); next }

        track_ids <- c(anchor_ids,
                       if (K > length(anchor_ids))
                           seq.int(from = max(anchor_ids, 0L) + 1L, length.out = K - length(anchor_ids))
                       else integer(0))
        anchor_set <- which(track_ids %in% anchor_ids)
        slack_set  <- which(!(track_ids %in% anchor_ids))

        # Get per-census obs and DBH
        obs_dbh <- lapply(census_range, function(cc) td[CensusID == cc & !is.na(DBH), DBH])

        # Median dates for intervals
        has_date <- "ExactDate" %in% names(td)
        med_dates <- sapply(census_range, function(cc) {
            d <- td[CensusID == cc & !is.na(DBH), ExactDate]
            if (length(d) > 0L) as.numeric(median(d, na.rm = TRUE)) else NA_real_
        })

        prune_min <- -0.625; prune_max <- 6.25

        # Backward per-interval constrained enumeration (mirrors dp_global_dp.R)
        all_mats <- vector("list", n_census)
        n_states_c <- integer(n_census)
        allowed <- vector("list", n_census)
        ok <- TRUE

        for (p in seq.int(n_census, 1L, by = -1L)) {
            n_obs <- length(obs_dbh[[p]])

            if (p == n_census) {
                # Anchor census: pin observed stems to their TrueStemID tracks
                if (n_obs > 0L) {
                    anchor_tidx <- match(anchor_obs$TrueStemID, track_ids)
                    allowed[[p]] <- lapply(seq_len(n_obs), function(j) {
                        if (!is.na(anchor_tidx[j])) anchor_tidx[j] else seq_len(K)
                    })
                }
            } else if (n_obs > 0L && has_date &&
                       is.finite(med_dates[p]) && is.finite(med_dates[p + 1L]) &&
                       length(obs_dbh[[p + 1L]]) > 0L && length(anchor_set) > 0L) {
                # Non-anchor: propagate growth constraints from next census
                dt <- (med_dates[p + 1L] - med_dates[p]) / 365.25
                if (is.finite(dt) && dt > 0) {
                    dbh_next <- obs_dbh[[p + 1L]]
                    al_next  <- allowed[[p + 1L]]
                    al <- vector("list", n_obs)
                    for (oi in seq_len(n_obs)) {
                        ok_tracks <- integer(0)
                        for (tk in anchor_set) {
                            feasible <- FALSE
                            for (oj in seq_along(dbh_next)) {
                                if (!(tk %in% al_next[[oj]])) next
                                rate <- (dbh_next[oj] - obs_dbh[[p]][oi]) / dt
                                if (rate >= prune_min && rate <= prune_max) { feasible <- TRUE; break }
                            }
                            if (feasible) ok_tracks <- c(ok_tracks, tk)
                        }
                        al[[oi]] <- sort(c(ok_tracks, slack_set))
                    }
                    allowed[[p]] <- al
                } else {
                    if (n_obs > 0L) allowed[[p]] <- replicate(n_obs, seq_len(K), simplify = FALSE)
                }
            } else {
                if (n_obs > 0L) allowed[[p]] <- replicate(n_obs, seq_len(K), simplify = FALSE)
            }

            use_constrained <- !is.null(allowed[[p]]) && !all(lengths(allowed[[p]]) == K)
            if (use_constrained) {
                mat <- enumerate_states_constrained(K, n_obs, allowed[[p]], DP_MAX_STATES)
            } else {
                mat <- enumerate_states_injective(K, n_obs, max_states = DP_MAX_STATES)
            }
            if (is.null(mat)) { ok <- FALSE; break }
            all_mats[[p]] <- mat
            n_states_c[p] <- nrow(mat)
        }

        if (!ok) {
            cat(sprintf("  %-10s K=%d  max_obs=%d  EXCEEDS DP_MAX_STATES after constraining\n", tg, K, r$max_obs))
            next
        }

        # Compute constrained edges
        edges <- 0
        for (p in seq_len(n_census - 1L)) edges <- edges + as.numeric(n_states_c[p]) * as.numeric(n_states_c[p + 1L])
        max_c <- max(n_states_c)
        reduction <- round(100 * (1 - max_c / r$max_states_per_census), 1)
        est_time <- predict_seconds(edges)

        time_str <- if (est_time < 60) sprintf("%.1fs", est_time)
                    else if (est_time < 3600) sprintf("%.1f min", est_time / 60)
                    else sprintf("%.1f hrs", est_time / 3600)

        cat(sprintf("  %-10s K=%2d  max_obs=%d  states: %s -> %s  (%+.0f%%)  edges: %s  est: %s\n",
                    tg, K, r$max_obs,
                    format(r$max_states_per_census, big.mark = ","),
                    format(max_c, big.mark = ","),
                    -reduction,
                    format(edges, big.mark = ","),
                    time_str))
    }
}

# ══════════════════════════════════════════════════════════════════════════
# STEP 4: Per-census detail for hardest tag
# ══════════════════════════════════════════════════════════════════════════
hardest <- top[1L]
cat(sprintf("\n=== Hardest DP-feasible Tag: %s ===\n", hardest$Tag))
cat(sprintf("  Species:               %s\n", hardest$Species))
cat(sprintf("  K (identity tracks):   %d\n", hardest$K))
cat(sprintf("  max_obs (per census):  %d\n", hardest$max_obs))
cat(sprintf("  Censuses:              %s\n", hardest$census_range))
cat(sprintf("  Max states/census:     %s\n", format(hardest$max_states_per_census, big.mark = ",")))
cat(sprintf("  Total states:          %s\n", format(hardest$total_states, big.mark = ",")))
cat(sprintf("  Unpruned edges:        %s\n", format(hardest$estimated_edges_unpruned, big.mark = ",")))
cat(sprintf("  Estimated time:        %.1f seconds\n", hardest$predicted_seconds))

# Per-census breakdown
htd <- x[Tag == hardest$Tag]
h_census <- sort(unique(htd[!is.na(DBH) & CensusID <= hardest$eff_anchor, CensusID]))
cat("\n  Per-census breakdown:\n")
cat(sprintf("  %10s %8s %15s\n", "Census", "N_Obs", "N_States"))
for (cc in h_census) {
    n_obs <- nrow(htd[CensusID == cc & !is.na(DBH)])
    K <- hardest$K
    ns <- if (n_obs == 0L || n_obs > K) 0
          else prod(seq.int(from = K, to = K - n_obs + 1L, by = -1L))
    cat(sprintf("  %10d %8d %15s\n", cc, n_obs, format(ns, big.mark = ",")))
}

# ══════════════════════════════════════════════════════════════════════════
# STEP 5 (optional): Run actual DP on a specific tag
# ══════════════════════════════════════════════════════════════════════════
if (nzchar(RUN_TAG)) {
    cat(sprintf("\n=== Running actual DP on Tag %s ===\n", RUN_TAG))
    source(here("dp_global", "R", "dp_global_main.R"))

    td <- copy(x[Tag == RUN_TAG])
    if (nrow(td) == 0L) {
        cat(sprintf("  Tag %s not found.\n", RUN_TAG))
    } else {
        cat(sprintf("  Rows: %d\n", nrow(td)))
        t0 <- proc.time()[["elapsed"]]
        result <- tryCatch(
            match_stems_dp_global_backward_marginals_batch(
                tree_data = td,
                min_growth = -0.5, max_growth = 5,
                anchor_start = 7L, max_tracks = NULL,
                slack_tracks = 1L,
                slack_require_anchor_recruitable = TRUE,
                slack_require_anchor_eps = 1e-6,
                max_states = DP_MAX_STATES,
                temperature = 1.0, posterior_top_k = 2L,
                eps_tiebreak = 1e-6,
                allow_provisional_anchor = TRUE,
                use_measurement_error = TRUE,
                meas_sd1_a = 0.0062, meas_sd1_b = 0.0904,
                meas_sd2 = 4.64, meas_p_big = 0.05,
                fallback_growth_forms = character(0),
                posterior_samples = 0L,
                prune_hard = TRUE,
                prune_min_growth = -0.625, prune_max_growth = 6.25,
                prune_use_bio_bounds = TRUE,
                prune_recruit_max_dbh = NULL,
                prune_use_bio_recruit = TRUE,
                non_taper_corrected_growth_forms = c("palm", "strangler_fig", "tree_fern"),
                non_taper_corrected_prune_min_growth = -0.625,
                non_taper_corrected_prune_max_growth = 6.25,
                hom_tolerance_scale = 2.0,
                verbose = TRUE, chunk_id = NULL,
                allow_segment_split = TRUE,
                post_segment_all_recruits = FALSE
            ),
            error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL }
        )
        elapsed <- proc.time()[["elapsed"]] - t0
        if (!is.null(result)) {
            method <- paste(unique(result$ReconstructionMethod), collapse = ", ")
            cat(sprintf("\n  Completed in %.1f seconds (%.1f min)\n", elapsed, elapsed / 60))
            cat(sprintf("  Method: %s\n", method))
        } else {
            cat(sprintf("\n  FAILED after %.1f seconds\n", elapsed))
        }
    }
}

cat("\nDone.\n")
