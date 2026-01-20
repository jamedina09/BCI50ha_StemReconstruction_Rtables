# Tests for estimate_bio_pars interval handling
# Run this script with Rscript to see assertions and printed outputs
rm(list = ls())

library(data.table)
source(file.path("dp_global", "R", "dp_global_biol.R"))

run_scenario <- function(df, col_name = "Bio_IntervalYears", desc = "") {
    cat("\nScenario:", desc, "\n")
    if (!is.null(col_name)) {
        res <- estimate_bio_pars(df, interval_years = NULL, interval_col_candidates = col_name)
    } else {
        res <- estimate_bio_pars(df, interval_years = 5)
    }
    cat("Inferred interval:", res$interval$inferred_interval_years, "\n")
    cat(
        "Per-pair intervals (n =", length(res$interval$per_pair_intervals), "):",
        paste(head(res$interval$per_pair_intervals, 10), collapse = ", "), "\n"
    )
    cat("Pairs candidate count:", res$interval$pairs_candidate_count, "\n")
    cat("Pairs filled with scalar:", res$interval$pairs_filled_with_scalar_count, "\n")
    cat("Pairs dropped:", res$interval$pairs_dropped_count, "\n")
    invisible(res)
}

# Scenario A: uniform per-census column
build_tag <- function(tag, n_cens, start_dbh, inc, interval) {
    cens <- seq_len(n_cens)
    data.frame(
        Tag = tag,
        CensusID = cens,
        DBH = start_dbh + (cens - 1) * inc,
        TrueStemID = 1,
        Bio_IntervalYears = interval
    )
}

df1 <- do.call(rbind, list(
    build_tag(1, 5, 10, 1.5, 5),
    build_tag(2, 5, 20, 2.0, 5)
))

res1 <- run_scenario(as.data.table(df1), col_name = "Bio_IntervalYears", desc = "uniform per-census = 5")
stopifnot(res1$interval$inferred_interval_years == 5)
stopifnot(all(res1$interval$per_pair_intervals == 5))

# Scenario B: mixed per-census with some NA (should infer median)
df2 <- do.call(rbind, list(
    build_tag(1, 5, 10, 1.5, 5),
    build_tag(2, 5, 20, 2.0, NA),
    build_tag(3, 5, 15, 1.2, 6)
))

res2 <- run_scenario(as.data.table(df2), col_name = "Bio_IntervalYears", desc = "mixed, tag2 NA, tag1=5, tag3=6")
# inferred should be median of available (5,6) => 5.5
stopifnot(abs(res2$interval$inferred_interval_years - 5.5) < 1e-8)

# Scenario C: some per-row t1 missing but t0 available (test fallback to t0)
# We'll simulate by making Bio_IntervalYears NA for census 2 only
make_missing_t1 <- function(df) {
    dt <- as.data.table(df)
    dt[CensusID == 2, Bio_IntervalYears := NA]
    dt
}

df3 <- do.call(rbind, list(
    build_tag(1, 5, 10, 1.5, 5),
    build_tag(2, 5, 20, 2.0, 4),
    build_tag(3, 5, 15, 1.2, 6)
))
df3 <- make_missing_t1(df3)
res3 <- run_scenario(as.data.table(df3), col_name = "Bio_IntervalYears", desc = "t1 missing for census 2 rows")
# After inference, inferred should be median of available pairs (including fallback to t0)
stopifnot(is.finite(res3$interval$inferred_interval_years))

cat("\nAll interval tests passed.\n")
