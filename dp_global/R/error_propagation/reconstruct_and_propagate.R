# Error propagation helpers for posterior path reconstructions
# Location: dp_global/R/error_propagation/reconstruct_and_propagate.R
#
# Purpose:
# - Provide utilities to load DP-produced posterior "paths" files (unique full-path
#   reconstructions and their probabilities), parse the compact recon structure,
#   apply each unique reconstruction to observed DBH data, compute per-path and
#   per-stem growth summaries, and offer Monte Carlo sampling for downstream
#   uncertainty propagation.
#
# Key functions:
# - load_posterior_paths(paths_file): read a *_paths.csv / .feather / .rds file and
#   parse the compact `recon` text into a `recon_parsed` list-column (data.table per path).
# - apply_paths_compute_growth(paths_dt, obs_dt, strict = TRUE, plot = FALSE, plot_draws = 1000):
#   apply each path to observed data and compute per-path mean/total growth;
#   optional `strict` aborts on mismatches; `plot` draws a sampled distribution
#   of mean-growth values (returns a ggplot object if ggplot2 is available).
# - sample_posterior_paths(paths_dt, n): draw Monte Carlo path indices according to path probability.
# - sample_apply_growth(paths_dt, obs_dt, n): combine sampling and per-path growth metrics
#   to produce n simulated mean-growth values useful for Monte Carlo propagation.
#
# Design notes & assumptions:
# - `paths` files contain columns: path_sig (string signature), path_prob (numeric), recon ("CensusID:StemID;...").
# - Recon `recon` strings are parsed robustly; if recon entries and observed counts
#   do not match per census the function either warns (strict=FALSE) or errors (strict=TRUE).
# - The matching between recon entries and observed rows uses a within-census
#   ordering (obs_pos) to avoid ambiguous cartesian joins.
#
# Outputs & usage:
# - Primary outputs are `paths_summary`, `per_path_growth`, and expected mean/total growth.
# - When `plot=TRUE` a histogram of sampled mean-growth (according to path_prob)
#   is rendered and a ggplot object is returned (if ggplot2 is installed).
#
# Example (interactive R):
#   library(data.table)
#   source("dp_global/R/error_propagation/reconstruct_and_propagate.R")
#   paths_dt <- load_posterior_paths("posteriors/tag_20_posterior_samples_..._paths.csv")
#   obs <- data.table::fread("data_simulation/data/simulated_data_1.csv")[Tag == 20]
#   res <- apply_paths_compute_growth(paths_dt, obs, strict = FALSE, plot = TRUE)
#   print(res$expected_mean_growth)
#   print(res$plot)  # ggplot object or NULL
#
# See README.md in this folder for more detailed guidance and command-line examples.

# Load posterior paths summary file (CSV / Feather / RDS)
#
# paths_file: Path to one of the `*_paths.csv` (or feather / rds) files produced by the DP
#            (must contain columns: `path_sig`, `path_prob`, `recon`).
# Returns a data.table with parsed `recon_parsed` list-column (each element is a
# data.table with columns `CensusID` (integer) and `ReconstructedStemID` (integer)).
# Example:
# paths_dt <- load_posterior_paths("posteriors/tag_20_posterior_samples_..._paths.csv")
load_posterior_paths <- function(paths_file) {
    stopifnot(requireNamespace("data.table", quietly = TRUE))
    ext <- tolower(tools::file_ext(paths_file))
    if (ext %in% c("csv")) {
        dt <- data.table::fread(paths_file)
    } else if (ext %in% c("feather", "feath")) {
        if (!requireNamespace("arrow", quietly = TRUE)) stop("Install 'arrow' to read feather files")
        dt <- as.data.table(arrow::read_feather(paths_file))
    } else if (ext %in% c("rds")) {
        x <- readRDS(paths_file)
        # If rds contains a list 'paths' or similar, try to find path_sig/path_prob/recon
        if (is.list(x) && !is.data.table(x)) {
            # try common layout
            if (is.data.table(x$paths)) dt <- x$paths else stop("rds file structure not recognized")
        } else if (is.data.frame(x)) {
            dt <- as.data.table(x)
        } else {
            stop("Unsupported rds structure for paths file")
        }
    } else {
        stop("Unsupported file extension: ", ext)
    }

    # Basic validation
    if (!all(c("path_sig", "path_prob", "recon") %in% names(dt))) {
        stop("paths file must contain columns: path_sig, path_prob, recon")
    }
    # Parse recon strings like "1:8;2:3;3:4" into data.tables
    parse_recon <- function(s) {
        if (is.na(s) || !nzchar(s)) {
            return(data.table::data.table(CensusID = integer(0), ReconstructedStemID = integer(0)))
        }
        # Split on semicolon, trim whitespace, and extract left/right of first ':' robustly
        parts <- strsplit(s, ";", fixed = TRUE)[[1]]
        parts <- parts[nzchar(trimws(parts))]
        if (length(parts) == 0) {
            return(data.table::data.table(CensusID = integer(0), ReconstructedStemID = integer(0)))
        }
        civ <- as.integer(vapply(parts, function(x) {
            sub(":.*$", "", trimws(x))
        }, FUN.VALUE = ""))
        riv <- as.integer(vapply(parts, function(x) {
            sub("^.*:", "", trimws(x))
        }, FUN.VALUE = ""))
        data.table::data.table(CensusID = civ, ReconstructedStemID = riv)
    }

    dt[, recon_parsed := lapply(recon, parse_recon)]
    # Ensure path_prob numeric and normalized defensively
    dt[, path_prob := as.numeric(path_prob)]
    tot <- sum(dt$path_prob, na.rm = TRUE)
    if (tot <= 0 || is.na(tot)) stop("Invalid path_prob in paths file")
    dt[, path_prob := path_prob / sum(path_prob, na.rm = TRUE)]

    data.table::setorder(dt, -path_prob)
    dt
}

# Apply each posterior path to observed data and compute growth metrics per-path
#
# paths_dt: data.table returned by `load_posterior_paths()`
# obs_dt: Observed data.table with at least columns `CensusID` and `DBH` (and optional `Tag`)
# census_col: Name of census id column in `obs_dt` (default: "CensusID")
# dbh_col: Name of DBH column in `obs_dt` (default: "DBH")
# tag_col: Optional tag column name; when present path `Tag` will be checked against it
# Returns: list with elements: paths_summary (data.table per-path aggregates) and per_path_growth (data.table per-path per-stem growths)
# Example:
# obs <- data.table::fread("data_simulation/data/simulated_data_1.csv")[Tag == 20]
# paths_dt <- load_posterior_paths("posteriors/tag_20_posterior_samples_..._paths.csv")
# res <- apply_paths_compute_growth(paths_dt, obs)
apply_paths_compute_growth <- function(paths_dt, obs_dt, census_col = "CensusID", dbh_col = "DBH", tag_col = NULL, obs_row_id_col = NULL, strict = TRUE, plot = FALSE, plot_draws = 1000) {
    stopifnot(requireNamespace("data.table", quietly = TRUE))
    # Minimal validation
    if (!("recon_parsed" %in% names(paths_dt))) stop("paths_dt must include recon_parsed list-column (use load_posterior_paths)")
    obs_dt <- data.table::as.data.table(obs_dt)
    if (!(census_col %in% names(obs_dt))) stop("obs_dt missing census column: ", census_col)
    if (!(dbh_col %in% names(obs_dt))) stop("obs_dt missing dbh column: ", dbh_col)
    # Detect optional obs_row_id column (preferred for robust matching)
    if (!is.null(obs_row_id_col) && obs_row_id_col %in% names(obs_dt)) {
        use_obs_row_id <- obs_row_id_col
    } else if ("obs_row_id" %in% names(obs_dt)) {
        use_obs_row_id <- "obs_row_id"
    } else {
        use_obs_row_id <- NULL
    }
    # strict: if TRUE, abort on any path<->census mismatches; if FALSE, summarize mismatches as warnings
    strict <- isTRUE(strict)
    plot <- isTRUE(plot)
    plot_draws <- as.integer(plot_draws)
    if (plot_draws <= 0) plot_draws <- 1000

    # For each path, join recon mapping to obs and compute stem-level growths
    per_path_list <- vector("list", nrow(paths_dt))
    per_path_stems <- vector("list", 0)
    mismatch_msgs <- character(0)

    for (i in seq_len(nrow(paths_dt))) {
        path_row <- paths_dt[i]
        # defensive try/catch per-path so one bad path doesn't abort everything
        tryCatch(
            {
                recon_dt <- path_row[["recon_parsed"]][[1]]
                if (nrow(recon_dt) == 0) {
                    per_path_list[[i]] <- data.table::data.table(path_sig = path_row$path_sig, path_prob = path_row$path_prob, mean_growth = NA_real_, total_growth = NA_real_, n_stems = 0L)
                    next
                }

                # Prepare observation positions
                obs2 <- data.table::copy(obs_dt)
                obs2[, obs_pos := seq_len(.N), by = census_col]

                # Try to detect obs_row_id-based recon entries when available
                if (!is.null(use_obs_row_id)) {
                    # Attempt to interpret recon left-values as obs_row_id where they match
                    recon_dt[, ObsRowID := NA_integer_]
                    # vectorized matching: compare as character to avoid type issues
                    match_idx <- match(as.character(recon_dt$CensusID), as.character(obs2[[use_obs_row_id]]))
                    recon_dt$ObsRowID[!is.na(match_idx)] <- as.integer(as.character(recon_dt$CensusID[!is.na(match_idx)]))
                    # Split recon entries into those using ObsRowID and those using CensusID/pos
                    recon_by_id <- recon_dt[!is.na(ObsRowID)]
                    recon_by_pos <- recon_dt[is.na(ObsRowID)]
                } else {
                    recon_by_id <- data.table::data.table()
                    recon_by_pos <- recon_dt
                }

                # Handle recon entries that match by CensusID/obs_pos (old behavior)
                if (nrow(recon_by_pos) > 0) {
                    recon_by_pos[, obs_pos := seq_len(.N), by = CensusID]

                    # Compare counts per census and collect messages
                    obs_counts <- obs2[, .N, by = census_col]
                    recon_counts <- recon_by_pos[, .N, by = "CensusID"]
                    cmp <- merge(obs_counts, recon_counts, by.x = census_col, by.y = "CensusID", all = TRUE)
                    if (nrow(cmp) > 0) {
                        for (r in seq_len(nrow(cmp))) {
                            cval <- cmp[[census_col]][r]
                            on <- cmp$N.x[r]
                            rn <- cmp$N.y[r]
                            if (is.na(on) && !is.na(rn)) {
                                mismatch_msgs <- c(mismatch_msgs, sprintf("path %s contains recon for CensusID=%s but no observations exist; entries will be ignored", path_row$path_sig, cval))
                            } else if (!is.na(on) && is.na(rn)) {
                                mismatch_msgs <- c(mismatch_msgs, sprintf("path %s missing recon entries for CensusID=%s (obs=%d); path will be matched partially", path_row$path_sig, cval, on))
                            } else if (!is.na(on) && !is.na(rn) && rn != on) {
                                mismatch_msgs <- c(mismatch_msgs, sprintf("path %s has recon entries=%d vs obs=%d at CensusID=%s; using first %d pairs", path_row$path_sig, rn, on, cval, min(on, rn)))
                            }
                        }
                    }

                    obs_counts_map <- setNames(as.integer(obs_counts$N), as.character(obs_counts[[census_col]]))
                    recon_by_pos <- recon_by_pos[as.character(CensusID) %in% names(obs_counts_map)]
                    if (nrow(recon_by_pos) > 0) {
                        recon_by_pos[, obs_count := obs_counts_map[as.character(CensusID)]]
                        recon_by_pos <- recon_by_pos[!is.na(obs_count) & obs_pos <= obs_count]
                        recon_by_pos[, obs_count := NULL]
                    } else {
                        # All recon_by_pos entries were dropped
                        recon_by_pos <- data.table::data.table()
                    }
                }

                # Now build merges
                merged_parts <- list()
                if (nrow(recon_by_id) > 0) {
                    # Join by ObsRowID directly
                    merged_id <- merge(obs2, recon_by_id[, .(ObsRowID, ReconstructedStemID)], by.x = use_obs_row_id, by.y = "ObsRowID", all = FALSE, allow.cartesian = FALSE)
                    if (nrow(merged_id) > 0) merged_parts[[length(merged_parts)+1]] <- merged_id
                    # Reporting for recon_by_id entries that had no matching obs rows
                    unmatched_ids <- recon_by_id[!as.character(CensusID) %in% as.character(obs2[[use_obs_row_id]])]
                    if (nrow(unmatched_ids) > 0) {
                        mismatch_msgs <- c(mismatch_msgs, sprintf("path %s contains recon for ObsRowIDs not present in observations; those entries will be ignored", path_row$path_sig))
                    }
                }

                if (nrow(recon_by_pos) > 0) {
                    merged_pos <- merge(obs2, recon_by_pos[, .(CensusID, obs_pos, ReconstructedStemID)], by.x = c(census_col, "obs_pos"), by.y = c("CensusID", "obs_pos"), all = FALSE, allow.cartesian = FALSE)
                    if (nrow(merged_pos) > 0) merged_parts[[length(merged_parts)+1]] <- merged_pos
                }

                if (length(merged_parts) == 0) {
                    per_path_list[[i]] <- data.table::data.table(path_sig = path_row$path_sig, path_prob = path_row$path_prob, mean_growth = NA_real_, total_growth = NA_real_, n_stems = 0L)
                    next
                }
                merged <- data.table::rbindlist(merged_parts, use.names = TRUE, fill = TRUE)

                if (nrow(merged) == 0) {
                    per_path_list[[i]] <- data.table::data.table(path_sig = path_row$path_sig, path_prob = path_row$path_prob, mean_growth = NA_real_, total_growth = NA_real_, n_stems = 0L)
                    next
                }

                # Compute growth per reconstructed stem id
                data.table::setorderv(merged, c("ReconstructedStemID", census_col))
                merged[, dbh_next := data.table::shift(get(dbh_col), type = "lead"), by = ReconstructedStemID]
                merged[, census_next := data.table::shift(get(census_col), type = "lead"), by = ReconstructedStemID]
                merged[, growth := ifelse(!is.na(dbh_next) & !is.na(get(dbh_col)), dbh_next - get(dbh_col), NA_real_)]
                merged[, census_delta := census_next - get(census_col)]
                merged_valid <- merged[!is.na(growth) & census_delta == 1]
                if (nrow(merged_valid) == 0) {
                    per_path_list[[i]] <- data.table::data.table(path_sig = path_row$path_sig, path_prob = path_row$path_prob, mean_growth = NA_real_, total_growth = NA_real_, n_stems = 0L)
                    next
                }

                stem_agg <- merged_valid[, .(stem_total_growth = sum(growth, na.rm = TRUE), stem_n = .N), by = ReconstructedStemID]
                mean_growth <- mean(stem_agg$stem_total_growth / stem_agg$stem_n, na.rm = TRUE)
                total_growth <- sum(stem_agg$stem_total_growth, na.rm = TRUE)
                per_path_list[[i]] <- data.table::data.table(path_sig = path_row$path_sig, path_prob = path_row$path_prob, mean_growth = mean_growth, total_growth = total_growth, n_stems = nrow(stem_agg))

                # attach per-stem details
                stem_agg[, path_sig := path_row$path_sig]
                stem_agg[, path_prob := path_row$path_prob]
                per_path_stems[[length(per_path_stems) + 1L]] <- stem_agg
                NULL
            },
            error = function(e) {
                warning(sprintf("Error while processing path %s: %s", paths_dt$path_sig[i], e$message))
                per_path_list[[i]] <- data.table::data.table(path_sig = paths_dt$path_sig[i], path_prob = paths_dt$path_prob[i], mean_growth = NA_real_, total_growth = NA_real_, n_stems = 0L)
                NULL
            }
        )
    }

    per_path_dt <- data.table::rbindlist(per_path_list, use.names = TRUE, fill = TRUE)
    # normalized path_prob defensive
    per_path_dt[, path_prob := as.numeric(path_prob)]
    per_path_dt[, path_prob := path_prob / sum(path_prob, na.rm = TRUE)]

    per_path_growth <- if (length(per_path_stems) > 0) data.table::rbindlist(per_path_stems, use.names = TRUE, fill = TRUE) else data.table::data.table()
    if (nrow(per_path_growth) > 0 && "path_prob" %in% names(per_path_growth)) {
        per_path_growth[, path_prob := as.numeric(path_prob)]
        per_path_growth[, path_prob := path_prob / sum(path_prob, na.rm = TRUE)]
    } else if (nrow(per_path_growth) > 0) {
        # If per_path_growth lacks path_prob (unexpected), set NA and avoid normalization
        per_path_growth[, path_prob := NA_real_]
    }

    # Aggregate mismatch messages and either warn or error depending on `strict`
    if (exists("mismatch_msgs") && length(mismatch_msgs) > 0L) {
        unique_msgs <- unique(mismatch_msgs)
        nmsgs <- length(unique_msgs)
        warn_examples <- paste(head(unique_msgs, 5), collapse = "\n")
        if (strict) {
            stop(sprintf("Encountered %d path-census mismatch messages (strict=TRUE): examples:\n%s", nmsgs, warn_examples))
        } else {
            warning(sprintf("Encountered %d path-census mismatch messages; examples:\n%s", nmsgs, warn_examples))
        }
    }

    expected_mean_growth <- sum(per_path_dt$mean_growth * per_path_dt$path_prob, na.rm = TRUE)
    expected_total_growth <- sum(per_path_dt$total_growth * per_path_dt$path_prob, na.rm = TRUE)

    # Optional plot: sample mean_growth according to path_prob to display distribution
    plot_obj <- NULL
    if (plot) {
        # require ggplot2 if available, otherwise use base hist
        draws <- rep(NA_real_, plot_draws)
        if (nrow(per_path_dt) > 0) {
            probs <- per_path_dt$path_prob
            # Some paths may have NA mean_growth (skip these in sampling)
            valid_idx <- which(!is.na(per_path_dt$mean_growth) & is.finite(per_path_dt$mean_growth) & probs > 0)
            if (length(valid_idx) > 0) {
                probs_valid <- probs[valid_idx]
                probs_valid <- probs_valid / sum(probs_valid)
                draws <- sample(per_path_dt$mean_growth[valid_idx], size = plot_draws, replace = TRUE, prob = probs_valid)
            } else {
                draws <- rep(NA_real_, plot_draws)
            }
        }
        # compute weighted sd for annotation
        wmean <- expected_mean_growth
        wsd <- NA_real_
        if (nrow(per_path_dt) > 0 && any(!is.na(per_path_dt$mean_growth))) {
            mg <- per_path_dt$mean_growth
            pp <- per_path_dt$path_prob
            valid <- !is.na(mg)
            if (sum(pp[valid]) > 0) {
                wmean2 <- sum(mg[valid] * pp[valid], na.rm = TRUE)
                wsd <- sqrt(sum(pp[valid] * (mg[valid] - wmean2)^2, na.rm = TRUE))
            }
        }

        if (requireNamespace("ggplot2", quietly = TRUE)) {
            df <- data.frame(mean_growth = draws)
            library(ggplot2)
            p <- ggplot2::ggplot(df, ggplot2::aes(x = mean_growth)) +
                ggplot2::geom_histogram(ggplot2::aes(y = ..density..), bins = 30, fill = "lightblue", color = "black") +
                ggplot2::geom_vline(xintercept = wmean, color = "red", linetype = "dashed", size = 1) +
                ggplot2::labs(title = "Posterior distribution of mean growth (sampled paths)", x = "Mean growth (cm)", y = "Density") +
                ggplot2::theme_minimal()
            if (!is.na(wsd)) {
                p <- p + ggplot2::annotate("text", x = wmean, y = Inf, label = sprintf("mean=%.3f, sd=%.3f", wmean, wsd), vjust = 2, hjust = 0.5, color = "red")
            }
            plot_obj <- p
        } else {
            hist(draws, breaks = 30, main = "Posterior distribution of mean growth (sampled paths)", xlab = "Mean growth (cm)", col = "lightblue")
            abline(v = wmean, col = "red", lty = 2)
            if (!is.na(wsd)) mtext(sprintf("mean=%.3f sd=%.3f", wmean, wsd), side = 3)
            # no plot object available
            plot_obj <- NULL
        }
    }

    list(paths_summary = per_path_dt, per_path_growth = per_path_growth, expected_mean_growth = expected_mean_growth, expected_total_growth = expected_total_growth, plot = plot_obj)
}

# Draw Monte Carlo reconstructions from the posterior paths (useful to propagate uncertainty by simulation)
#
# paths_dt: data.table from `load_posterior_paths()`
# n: number of draws (e.g., 1000)
# Returns: data.table with columns Draw (1..n), path_sig, path_prob
# Example:
# draws <- sample_posterior_paths(paths_dt, n = 1000)
sample_posterior_paths <- function(paths_dt, n = 1000) {
    stopifnot(requireNamespace("data.table", quietly = TRUE))
    if (!("path_prob" %in% names(paths_dt))) stop("paths_dt must contain path_prob")
    probs <- paths_dt$path_prob
    picks <- sample.int(nrow(paths_dt), size = n, replace = TRUE, prob = probs)
    res <- data.table::data.table(Draw = seq_len(n), path_sig = paths_dt$path_sig[picks], path_prob = paths_dt$path_prob[picks])
    res
}

# Sample posterior paths and compute growth per draw (Monte Carlo propagation)
#
# paths_dt: data.table returned by `load_posterior_paths()`
# obs_dt: observed data table
# n: number of Monte Carlo draws
# Returns: data.table with columns Draw, path_sig, path_prob, mean_growth, total_growth
sample_apply_growth <- function(paths_dt, obs_dt, n = 1000) {
    stopifnot(requireNamespace("data.table", quietly = TRUE))
    # compute per-path growth metrics (this returns per-path table and expected values)
    res <- apply_paths_compute_growth(paths_dt, obs_dt)
    per_path_dt <- res$paths_summary
    # sample draws according to path_prob
    draws <- sample_posterior_paths(paths_dt, n = n)
    draws <- merge(draws, per_path_dt[, .(path_sig, mean_growth, total_growth)], by = "path_sig", all.x = TRUE)
    # Some sampled paths may have NA growth (unmatched); keep as NA to reflect uncertainty
    draws
}

# library(data.table)
# source("./dp_global/R/error_propagation/reconstruct_and_propagate.R")
# paths_dt <- load_posterior_paths("./dp_global/output/20260122_220613_unknown_T20_DP_MB_ME_gD_sD_kgD_ksD_rcpp/posteriors/tag_20_posterior_samples_20260122_220630_paths.csv")
# obs <- data.table::fread("data_simulation/data/simulated_data_1.csv")[Tag == 20]

# res <- apply_paths_compute_growth(paths_dt, obs, strict = FALSE, plot = TRUE)
# print(res$expected_mean_growth)
# print(res$plot)  # ggplot object or NULL