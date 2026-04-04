## Benchmark: worst-case BCI tags (K=9, max_obs=7, zero constraint reduction)
##
## These tags had integer overflow on edge counts in the constrained scan —
## truly the hardest. They typically hit resprout splits or igraph fallback.
##
## Usage:
##   Rscript dp_global/R/complexity/benchmark_worst_case_tags.R

library(data.table)
library(here)

source(here("dp_global", "R", "dp_global_main.R"))

x <- readRDS(here("bci_data", "bci_multistem_xrun_debug.rds"))

test_tags <- c("020978", "030209")

for (tg in test_tags) {
    td <- copy(x[Tag == tg])
    cat(sprintf("\n===== Running DP on Tag %s =====\n", tg))
    cat(sprintf("Rows: %d, Censuses: %s\n", nrow(td),
                paste(sort(unique(td$CensusID)), collapse = ",")))
    cat(sprintf("Obs per census: %s\n",
                paste(td[, .N, by = CensusID][order(CensusID)]$N, collapse = ",")))

    t0 <- proc.time()[["elapsed"]]
    result <- tryCatch(
        match_stems_dp_global_backward_marginals_batch(
            tree_data = td,
            min_growth = -0.5,
            max_growth = 5,
            anchor_start = 7L,
            max_tracks = NULL,
            slack_tracks = 1L,
            slack_require_anchor_recruitable = TRUE,
            slack_require_anchor_eps = 1e-6,
            max_states = 200000L,
            temperature = 1.0,
            posterior_top_k = 2L,
            eps_tiebreak = 1e-6,
            allow_provisional_anchor = TRUE,
            use_measurement_error = TRUE,
            meas_sd1_a = 0.0062,
            meas_sd1_b = 0.0904,
            meas_sd2 = 4.64,
            meas_p_big = 0.05,
            fallback_growth_forms = character(0),
            posterior_samples = 0L,
            prune_hard = TRUE,
            prune_min_growth = -0.625,
            prune_max_growth = 6.25,
            prune_use_bio_bounds = TRUE,
            prune_recruit_max_dbh = NULL,
            prune_use_bio_recruit = TRUE,
            non_taper_corrected_growth_forms = c("palm", "strangler_fig", "tree_fern"),
            non_taper_corrected_prune_min_growth = -0.625,
            non_taper_corrected_prune_max_growth = 6.25,
            hom_tolerance_scale = 2.0,
            verbose = TRUE,
            chunk_id = NULL,
            allow_segment_split = TRUE,
            post_segment_all_recruits = FALSE
        ),
        error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL }
    )
    elapsed <- proc.time()[["elapsed"]] - t0

    if (!is.null(result)) {
        method <- unique(result$ReconstructionMethod)
        fallback <- unique(result$DP_FallbackReason)
        cat(sprintf("\nTag %s: completed in %.1f seconds (%.1f min)\n", tg, elapsed, elapsed / 60))
        cat(sprintf("  Method: %s, Fallback: %s\n", paste(method, collapse = ","), paste(fallback, collapse = ",")))
    } else {
        cat(sprintf("\nTag %s: FAILED after %.1f seconds (%.1f min)\n", tg, elapsed, elapsed / 60))
    }
}
