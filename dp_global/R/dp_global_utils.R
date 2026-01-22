############################################################
# dp_global_utils.R
# Utility helpers for dp_global
############################################################

count_injective_states <- function(K, n_obs) {
    # PURPOSE
    # - Upper bound helper for DP state-space size.
    # - Counts the number of injective assignments of `n_obs` labeled observations
    #   into `K` labeled tracks (identity slots).
    #
    # INPUTS
    # - K: integer-ish; number of available tracks.
    # - n_obs: integer-ish; number of observed stems in a census.
    #
    # OUTPUT
    # - Numeric scalar (double): P(K, n_obs) = K*(K-1)*...*(K-n_obs+1) = K!/(K-n_obs)!
    # - Returns 0 when n_obs > K.
    # - Allows Inf on overflow.
    #
    # NOTES
    # - This grows factorially and is the main reason global DP can become expensive.
    if (!is.finite(K) || !is.finite(n_obs)) {
        return(NA_real_)
    }
    K <- as.integer(K)
    n_obs <- as.integer(n_obs)
    if (n_obs < 0L || K < 0L) {
        return(NA_real_)
    }
    if (n_obs == 0L) {
        return(1)
    }
    if (n_obs > K) {
        return(0)
    }
    prod(seq.int(from = K, to = K - n_obs + 1L, by = -1L))
}

resolve_interval_years <- function(tree_data,
                                   interval_years = NULL,
                                   interval_col_candidates = c(
                                       "Bio_IntervalYears",
                                       "IntervalYears",
                                       "interval_years",
                                       "census_interval_years",
                                       "CensusIntervalYears"
                                   )) {
    # PURPOSE
    # - Determine the census interval (years) for the DP.
    # - Supports the newer workflow where interval is stored in `tree_data` as a
    #   column (similar to Bio_* columns), so callers don't have to pass it around.
    #
    # CONTRACT
    # - The DP currently assumes a single constant interval for all adjacent census
    #   transitions within a (Tag, species) group.
    # - If a candidate interval column contains multiple distinct finite values,
    #   we error rather than silently picking one.

    if (!is.null(interval_years)) {
        interval_years <- suppressWarnings(as.numeric(interval_years))
        if (!is.finite(interval_years) || is.na(interval_years) || interval_years <= 0) {
            stop("interval_years must be a positive finite number.", call. = FALSE)
        }
        return(interval_years)
    }

    if (is.null(tree_data) || !is.data.frame(tree_data)) {
        stop(
            "interval_years not provided and tree_data is not a data.frame/data.table; ",
            "provide interval_years or add a column like 'Bio_IntervalYears'.",
            call. = FALSE
        )
    }

    for (col in interval_col_candidates) {
        if (!(col %in% names(tree_data))) next
        v <- suppressWarnings(as.numeric(tree_data[[col]]))
        v <- v[is.finite(v) & !is.na(v)]
        if (length(v) == 0L) next
        u <- unique(v)
        if (length(u) == 1L) {
            if (u[[1L]] <= 0) {
                stop("Interval column '", col, "' must be positive.", call. = FALSE)
            }
            return(u[[1L]])
        }
        stop(
            "Interval column '", col, "' has multiple distinct values (", paste(head(u, 10L), collapse = ", "),
            if (length(u) > 10L) ", ..." else "",
            "). The DP currently requires a single constant interval_years per group.",
            call. = FALSE
        )
    }

    stop(
        "interval_years not provided and no interval column found in tree_data. ",
        "Add a constant column like 'Bio_IntervalYears' (recommended) or pass interval_years explicitly.",
        call. = FALSE
    )
}

# Resolve interval years for a specific adjacent census pair (t0 -> t1).
# Preference order:
#  - explicit scalar `interval_years` argument if provided
#  - per-census value at CensusID == t1 (preferred)
#  - per-census value at CensusID == t0
#  - global constant in the interval column if present
# Errors when multiple distinct finite values are present for the pair.
resolve_interval_years_pair <- function(tree_data, t0, t1, interval_years = NULL,
                                       interval_col_candidates = c(
                                           "Bio_IntervalYears",
                                           "IntervalYears",
                                           "interval_years",
                                           "census_interval_years",
                                           "CensusIntervalYears"
                                       )) {
    if (!is.null(interval_years)) {
        interval_years <- suppressWarnings(as.numeric(interval_years))
        if (!is.finite(interval_years) || is.na(interval_years) || interval_years <= 0) {
            stop("interval_years must be a positive finite number.", call. = FALSE)
        }
        return(interval_years)
    }

    if (is.null(tree_data) || !is.data.frame(tree_data)) {
        stop(
            "interval_years not provided and tree_data is not a data.frame/data.table; ",
            "provide interval_years or add a column like 'Bio_IntervalYears'.",
            call. = FALSE
        )
    }

    for (col in interval_col_candidates) {
        if (!(col %in% names(tree_data))) next
        v <- suppressWarnings(as.numeric(tree_data[[col]]))

        # Prefer t1 values
        if (!missing(t1)) {
            v1 <- v[tree_data$CensusID == t1]
            v1 <- v1[is.finite(v1) & !is.na(v1) & v1 > 0]
            if (length(v1) > 0L) {
                u1 <- unique(v1)
                if (length(u1) == 1L) return(u1[[1L]])
                # Multiple per-census values: use the mean across rows (tolerant to small jitter)
                mean_v1 <- mean(v1, na.rm = TRUE)
                warning(
                    "Multiple interval values found for CensusID=", t1, "; using mean (", format(mean_v1, digits = 8), ")",
                    call. = FALSE
                )
                return(as.numeric(mean_v1))
            }
        }

        # Next prefer t0 values
        if (!missing(t0)) {
            v0 <- v[tree_data$CensusID == t0]
            v0 <- v0[is.finite(v0) & !is.na(v0) & v0 > 0]
            if (length(v0) > 0L) {
                u0 <- unique(v0)
                if (length(u0) == 1L) return(u0[[1L]])
                # Multiple per-census values: use the mean across rows (tolerant to small jitter)
                mean_v0 <- mean(v0, na.rm = TRUE)
                warning(
                    "Multiple interval values found for CensusID=", t0, "; using mean (", format(mean_v0, digits = 8), ")",
                    call. = FALSE
                )
                return(as.numeric(mean_v0))
            }
        }

        # Fallback: global constant across the data
        v_all <- v[is.finite(v) & !is.na(v) & v > 0]
        if (length(v_all) > 0L) {
            u_all <- unique(v_all)
            if (length(u_all) == 1L) return(u_all[[1L]])
            # Multiple global values: use mean across all rows as a tolerant fallback
            mean_all <- mean(v_all, na.rm = TRUE)
            warning(
                "Multiple interval values found globally in column '", col, "'; using mean (", format(mean_all, digits = 8), ")",
                call. = FALSE
            )
            return(as.numeric(mean_all))
        }
    }

    stop(
        "interval_years not provided and no interval information found for Census pair ", t0, "->", t1, ". Add per-census column like 'Bio_IntervalYears' or pass interval_years explicitly.",
        call. = FALSE
    )
}
