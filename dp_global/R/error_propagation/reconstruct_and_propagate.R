# Minimal helpers used by `attach_paths_to_output_run.R`
# This file intentionally contains only the small set of helpers needed to run
# example attachment workflows: `parse_recon()`, `load_posterior_paths()` and
# `attach_paths_to_output()`.

stopifnot(requireNamespace("data.table", quietly = TRUE))

# Parse recon strings. Support two formats:
#  - explicit pairs: "1:8;2:3;3:4"   (CensusID:ReconstructedStemID;...)
#  - compact path: "7-4-8-..."       (one ReconstructedStemID per census, hyphen-separated)
parse_recon <- function(s) {
    if (is.na(s) || !nzchar(s)) {
        return(data.table::data.table(CensusID = integer(0), ReconstructedStemID = integer(0)))
    }
    # If string contains a ':' assume explicit "CensusID:StemID" pairs
    if (grepl(":", s, fixed = TRUE)) {
        parts <- strsplit(s, ";", fixed = TRUE)[[1]]
        parts <- parts[nzchar(trimws(parts))]
        if (length(parts) == 0) {
            return(data.table::data.table(CensusID = integer(0), ReconstructedStemID = integer(0)))
        }
        civ <- as.integer(vapply(parts, function(x) sub(":.*$", "", trimws(x)), FUN.VALUE = ""))
        riv <- as.integer(vapply(parts, function(x) sub("^.*:", "", trimws(x)), FUN.VALUE = ""))
        return(data.table::data.table(CensusID = civ, ReconstructedStemID = riv))
    }

    # Otherwise attempt compact path parsing: common separators are '-' or ',' or whitespace
    sep <- if (grepl("-", s, fixed = TRUE)) "-" else if (grepl(",", s, fixed = TRUE)) "," else " "
    parts <- strsplit(s, sep, fixed = TRUE)[[1]]
    parts <- parts[nzchar(trimws(parts))]
    if (length(parts) == 0) {
        return(data.table::data.table(CensusID = integer(0), ReconstructedStemID = integer(0)))
    }
    riv <- as.integer(trimws(parts))
    civ <- seq_along(riv)
    data.table::data.table(CensusID = civ, ReconstructedStemID = riv)
}

# Load posterior paths summary file (CSV / Feather / RDS)
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
        if (is.list(x) && !is.data.table(x)) {
            if (is.data.table(x$paths)) dt <- x$paths else stop("rds file structure not recognized")
        } else if (is.data.frame(x)) {
            dt <- as.data.table(x)
        } else {
            stop("Unsupported rds structure for paths file")
        }
    } else {
        stop("Unsupported file extension: ", ext)
    }

    if (!all(c("path_sig", "path_prob", "recon") %in% names(dt))) {
        stop("paths file must contain columns: path_sig, path_prob, recon")
    }

    dt[, recon_parsed := lapply(recon, function(x) parse_recon(x))]
    dt[, path_prob := as.numeric(path_prob)]
    tot <- sum(dt$path_prob, na.rm = TRUE)
    if (tot <= 0 || is.na(tot)) stop("Invalid path_prob in paths file")
    dt[, path_prob := path_prob / sum(path_prob, na.rm = TRUE)]
    data.table::setorder(dt, -path_prob)
    dt
}

# Attach selected posterior paths to the main reconstruction output by matching ObsRowID
attach_paths_to_output <- function(paths, out, which = c("map", "top_n", "indices", "sample"), n = 1L, indices = NULL, obs_row_id_col = "obs_row_id", require_obs_rowid_match = TRUE, write_out = FALSE, out_file = NULL) {
    stopifnot(requireNamespace("data.table", quietly = TRUE))
    which <- match.arg(which)

    # load paths if a filename provided
    if (is.character(paths) && length(paths) == 1 && file.exists(paths)) {
        paths_dt <- load_posterior_paths(paths)
    } else if (is.data.table(paths) || is.data.frame(paths)) {
        paths_dt <- data.table::as.data.table(paths)
        if (!("recon_parsed" %in% names(paths_dt))) {
            if ("recon" %in% names(paths_dt)) {
                paths_dt[, recon_parsed := lapply(recon, function(x) parse_recon(x))]
            } else {
                stop("paths input must contain recon_parsed or recon column")
            }
        }
    } else {
        stop("Unsupported 'paths' input: provide file path or data.table")
    }

    out_was_path <- FALSE
    if (is.character(out) && length(out) == 1 && file.exists(out)) {
        out_dt <- data.table::fread(out)
        out_was_path <- TRUE
    } else if (is.data.table(out) || is.data.frame(out)) {
        out_dt <- data.table::as.data.table(out)
    } else {
        stop("Unsupported 'out' input: provide file path or data.table")
    }

    if (!(obs_row_id_col %in% names(out_dt))) {
        out_dt[, (obs_row_id_col) := seq_len(.N)]
    }

    sel_idx <- integer(0)
    if (which == "map") {
        sel_idx <- which.max(paths_dt$path_prob)
    } else if (which == "top_n") {
        n <- as.integer(n)
        sel_idx <- seq_len(min(n, nrow(paths_dt)))
    } else if (which == "indices") {
        sel_idx <- as.integer(indices)
    } else if (which == "sample") {
        n <- as.integer(n)
        sel_idx <- sample.int(nrow(paths_dt), size = n, replace = FALSE, prob = paths_dt$path_prob)
    }

    k <- 0L
    for (i in sel_idx) {
        k <- k + 1L
        path_sig <- as.character(paths_dt$path_sig[i])
        recon_dt <- paths_dt$recon_parsed[[i]]
        recon_dt_copy <- data.table::copy(recon_dt)
        matches <- match(as.character(recon_dt_copy$CensusID), as.character(out_dt[[obs_row_id_col]]))
        if (any(!is.na(matches))) {
            recon_dt_copy[, ObsRowID := as.integer(as.character(CensusID))]
            map_dt <- recon_dt_copy[!is.na(ObsRowID), .(ObsRowID = ObsRowID, ReconstructedStemID = ReconstructedStemID)]
            data.table::setkeyv(out_dt, obs_row_id_col)
            data.table::setkeyv(map_dt, "ObsRowID")
            new_recon_col <- paste0("DP_ReconstructedStemID_", k)
            new_sig_col <- paste0("DP_PathSig_", k)
            out_dt[, (new_recon_col) := as.integer(NA)]
            out_dt[map_dt, (new_recon_col) := i.ReconstructedStemID]
            out_dt[, (new_sig_col) := path_sig]
        } else {
            new_recon_col <- paste0("DP_ReconstructedStemID_", k)
            new_sig_col <- paste0("DP_PathSig_", k)
            out_dt[, (new_recon_col) := as.integer(NA)]
            out_dt[, (new_sig_col) := path_sig]
            warning(sprintf("Path %s did not contain ObsRowID matches; column %s filled with NA", path_sig, new_recon_col))
        }
    }

    if (write_out && out_was_path) {
        if (is.null(out_file)) {
            out_file <- sub("\\.csv$", "", out)
            out_file <- paste0(out_file, "_with_paths.csv")
        }
        data.table::fwrite(out_dt, out_file)
        message("Wrote output with attached paths to: ", out_file)
    }

    out_dt
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
