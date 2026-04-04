## Timing calibration: sample tags across 5 complexity tiers and estimate
## total BCI runtime.
##
## Tiers are defined by unconstrained max_states:
##   tier1: ≤10, tier2: 10-100, tier3: 100-1K, tier4: 1K-5K, tier5: 5K-40K
##
## Usage:
##   Rscript dp_global/R/complexity/benchmark_timing_calibration.R

library(data.table)
library(here)

source(here("dp_global", "R", "dp_global_main.R"))

x <- readRDS(here("bci_data", "bci_multistem_xrun_debug.rds"))

# Compute unconstrained max_states for all multi-stem tags
tags_info <- x[, .(n_stem = uniqueN(OriginalStemID),
                   censuses = uniqueN(CensusID),
                   max_obs = max(.SD[, .N, by = CensusID]$N)), by = Tag]

tags_info <- tags_info[n_stem > 1]
tags_info[, K := n_stem + 1L]
tags_info[, max_states := factorial(K) / factorial(K - max_obs)]

# Sample across tiers
set.seed(42)
tiers <- list(
    tier1  = tags_info[max_states <= 10],
    tier2  = tags_info[max_states > 10 & max_states <= 100],
    tier3  = tags_info[max_states > 100 & max_states <= 1000],
    tier4  = tags_info[max_states > 1000 & max_states <= 5000],
    tier5  = tags_info[max_states > 5000 & max_states <= 40000]
)

sampled <- rbindlist(lapply(names(tiers), function(nm) {
    dt <- tiers[[nm]]
    n_sample <- min(3, nrow(dt))
    if (n_sample == 0) return(NULL)
    dt[sample(.N, n_sample)][, tier := nm]
}))

cat("=== Tier counts ===\n")
for (nm in names(tiers)) {
    cat(sprintf("  %s: %d tags (sampling %d)\n", nm, nrow(tiers[[nm]]),
                min(3, nrow(tiers[[nm]]))))
}

cat("\n=== Running sampled tags ===\n")
results <- list()
for (i in seq_len(nrow(sampled))) {
    tg <- sampled$Tag[i]
    td <- copy(x[Tag == tg])

    t0 <- proc.time()[["elapsed"]]
    result <- tryCatch(
        match_stems_dp_global_backward_marginals_batch(
            tree_data = td,
            min_growth = -0.5, max_growth = 5,
            anchor_start = 7L, max_tracks = NULL,
            slack_tracks = 1L,
            slack_require_anchor_recruitable = TRUE,
            slack_require_anchor_eps = 1e-6,
            max_states = 40000L,
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
            verbose = FALSE, chunk_id = NULL,
            allow_segment_split = TRUE,
            post_segment_all_recruits = FALSE
        ),
        error = function(e) { cat("ERROR on tag ", tg, ":", conditionMessage(e), "\n"); NULL }
    )
    elapsed <- proc.time()[["elapsed"]] - t0

    method <- if (!is.null(result)) paste(unique(result$ReconstructionMethod), collapse = "|") else "FAIL"
    cat(sprintf("  Tag %s (tier=%s, K=%d, max_obs=%d, max_states=%s): %.3fs [%s]\n",
                tg, sampled$tier[i], sampled$K[i], sampled$max_obs[i],
                format(sampled$max_states[i], big.mark = ","), elapsed, method))

    results[[i]] <- data.table(
        Tag = tg, tier = sampled$tier[i],
        K = sampled$K[i], max_obs = sampled$max_obs[i],
        max_states = sampled$max_states[i],
        elapsed = elapsed, method = method
    )
}

results_dt <- rbindlist(results)
cat("\n=== Timing Summary ===\n")
print(results_dt[, .(Tag, tier, K, max_obs, max_states, elapsed, method)])

# Estimate total runtime
cat("\n=== Total Runtime Estimate ===\n")
cat(sprintf("Tier sizes: %s\n", paste(sapply(tiers, nrow), collapse = ", ")))

single <- x[, .(n_stem = uniqueN(OriginalStemID)), by = Tag][n_stem == 1]
cat(sprintf("Single-stem tags (trivial): %d\n", nrow(single)))

for (nm in names(tiers)) {
    tier_results <- results_dt[tier == nm]
    if (nrow(tier_results) > 0) {
        med <- median(tier_results$elapsed)
        n <- nrow(tiers[[nm]])
        est <- n * med
        cat(sprintf("  %s: %d tags × %.3fs (median) = %.0fs = %.1f min\n",
                    nm, n, med, est, est / 60))
    }
}
