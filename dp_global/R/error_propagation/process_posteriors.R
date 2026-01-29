# process_posteriors.R
# ---------------------
# Canonical parsing & posterior-processing helpers moved here (rename from
# `reconstruct_and_propagate.R`). This file contains:
# - `parse_recon()`
# - `load_posterior_paths()`
# - `attach_paths_to_output()`
# - `expand_draws()` / `aggregate_draws()` helpers for drawing and aggregating per-observation samples

stopifnot(requireNamespace("data.table", quietly = TRUE))

# Parse recon strings. Support explicit ObsRowID:ReconstructedStemID pairs. Legacy compact hyphen paths are supported
# for backwards compatibility but are discouraged (they lack explicit ObsRowID information).
parse_recon <- function(s) {
    if (is.na(s) || !nzchar(s)) {
        return(data.table::data.table(ObsRowID = integer(0), ReconstructedStemID = integer(0)))
    }
    # If string contains a ':' assume explicit "ObsRowID:StemID" pairs
    if (grepl(":", s, fixed = TRUE)) {
        parts <- strsplit(s, ";", fixed = TRUE)[[1]]
        parts <- parts[nzchar(trimws(parts))]
        if (length(parts) == 0) {
            return(data.table::data.table(ObsRowID = integer(0), ReconstructedStemID = integer(0)))
        }
        obsiv <- as.integer(vapply(parts, function(x) sub(":.*$", "", trimws(x)), FUN.VALUE = ""))
        riv <- as.integer(vapply(parts, function(x) sub("^.*:", "", trimws(x)), FUN.VALUE = ""))
        return(data.table::data.table(ObsRowID = obsiv, ReconstructedStemID = riv))
    }

    # Otherwise attempt compact path parsing: common separators are '-' or ',' or whitespace
    # NOTE: compact path lacks explicit ObsRowID; we return positional ObsRowID = 1..n and issue a warning.
    sep <- if (grepl("-", s, fixed = TRUE)) "-" else if (grepl(",", s, fixed = TRUE)) "," else " "
    parts <- strsplit(s, sep, fixed = TRUE)[[1]]
    parts <- parts[nzchar(trimws(parts))]
    if (length(parts) == 0) {
        return(data.table::data.table(ObsRowID = integer(0), ReconstructedStemID = integer(0)))
    }
    riv <- as.integer(trimws(parts))
    # positional indices (legacy): these are NOT real ObsRowIDs, and attachment should
    # attempt a conversion using out_dt (see attach_paths_to_output). Warn user.
    warning("Legacy compact path parsed without ObsRowID; returning positional ObsRowID = 1..n. Regenerate posterior samples to include ObsRowID for robust matching.")
    obsiv <- seq_along(riv)
    data.table::data.table(ObsRowID = obsiv, ReconstructedStemID = riv)
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
        if (is.list(x) && !data.table::is.data.table(x)) {
            if (data.table::is.data.table(x$paths)) dt <- x$paths else stop("rds file structure not recognized")
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
    } else if (data.table::is.data.table(paths) || is.data.frame(paths)) {
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
    } else if (data.table::is.data.table(out) || is.data.frame(out)) {
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

        # Primary behavior: prefer explicit ObsRowID for matching
        if ("ObsRowID" %in% names(recon_dt_copy)) {
            # Map ObsRowID to the output obs_row_id column name to make the join explicit.
            matches <- match(as.character(recon_dt_copy$ObsRowID), as.character(out_dt[[obs_row_id_col]]))
            new_recon_col <- paste0("DP_ReconstructedStemID_", k)
            new_sig_col <- paste0("DP_PathSig_", k)
            # initialize columns first (ensures they exist even if no matches)
            out_dt[, (new_recon_col) := as.integer(NA)]
            out_dt[, (new_sig_col) := path_sig]

            if (any(!is.na(matches))) {
                map_dt <- recon_dt_copy[!is.na(ObsRowID), .(ObsRowID = as.integer(as.character(ObsRowID)), ReconstructedStemID = ReconstructedStemID)]
                # rename to match the output key column for an explicit key-based join
                if (obs_row_id_col != "ObsRowID") {
                    data.table::setnames(map_dt, "ObsRowID", obs_row_id_col)
                }
                data.table::setkeyv(out_dt, obs_row_id_col)
                data.table::setkeyv(map_dt, obs_row_id_col)
                out_dt[map_dt, (new_recon_col) := i.ReconstructedStemID]
            } else {
                warning(sprintf("Path %s did not contain ObsRowID matches; column %s left as NA", path_sig, new_recon_col))
            }
        } else {
            # Enforce ObsRowID-only reconstructions. Legacy CensusID-based reconstructions
            # are no longer supported: fail loudly so users regenerate posterior samples
            # with ObsRowID encoded.
            new_recon_col <- paste0("DP_ReconstructedStemID_", k)
            new_sig_col <- paste0("DP_PathSig_", k)
            out_dt[, (new_recon_col) := as.integer(NA)]
            out_dt[, (new_sig_col) := path_sig]
            stop(sprintf("Path %s does not contain ObsRowID; legacy CensusID-based matching is removed. Regenerate posterior samples with ObsRowID encoded.", path_sig))
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

# Expand draws from per-path summaries: sample N draws using the summary table
# and convert compact 'recon' strings into long (Draw, ObsRowID, ReconstructedStemID)
expand_draws <- function(summary_dt, paths_dt, N = 1000L) {
    if (!all(c("Sample", "path_sig", "sample_prob") %in% names(summary_dt))) stop("summary_dt must contain Sample, path_sig, sample_prob")
    picks <- sample(summary_dt$Sample, size = as.integer(N), replace = TRUE, prob = summary_dt$sample_prob)
    draws <- summary_dt[J(picks), .(Draw = seq_len(N), Sample = picks, path_sig = path_sig), on = "Sample"]
    draws_paths <- merge(draws, paths_dt[, .(path_sig, recon)], by = "path_sig", all.x = TRUE)

    # Use the global `parse_recon()` (defined above) to parse recon strings
    res_list <- lapply(seq_len(nrow(draws_paths)), function(i) {
        dt <- parse_recon(draws_paths$recon[i])
        dt[, Draw := draws_paths$Draw[i]]
        dt
    })
    data.table::rbindlist(res_list, use.names = TRUE, fill = TRUE)
}

# Aggregate expanded draws into per-observation probabilities
aggregate_draws <- function(res_dt) {
    agg <- res_dt[, .(count = .N), by = .(ObsRowID, ReconstructedStemID)]
    agg[, prob := count / sum(count), by = ObsRowID]
    agg[order(ObsRowID, -prob)]
}

# Diagnostic: does the MAP joint path occur among sampled unique paths?
check_map_in_paths <- function(paths_dt, out_dt) {
    dp_rows <- out_dt[ReconstructionMethod == "dp"]
    if (nrow(dp_rows) == 0) return(list(found = FALSE, reason = "no dp rows"))

    sig <- paste0(dp_rows[order(obs_row_id), ReconstructedStemID], collapse = "-")
    which_match <- which(paths_dt$path_sig == sig)
    if (length(which_match) > 0) return(list(found = TRUE, idx = which_match[1], path_sig = sig))

    # compute best partial match for diagnostics
    main_map <- dp_rows[, .(ObsRowID = obs_row_id, ReconstructedStemID = ReconstructedStemID)]
    best <- NULL
    for (i in seq_len(nrow(paths_dt))) {
        recon_dt <- paths_dt$recon_parsed[[i]]
        common <- intersect(main_map$ObsRowID, recon_dt$ObsRowID)
        if (length(common) == 0) {
            frac <- 0
        } else {
            mm <- setNames(main_map$ReconstructedStemID, main_map$ObsRowID)
            rm <- setNames(recon_dt$ReconstructedStemID, recon_dt$ObsRowID)
            matches <- sum(vapply(common, function(k) mm[as.character(k)] == rm[as.character(k)], logical(1)))
            frac <- matches / length(common)
        }
        if (is.null(best) || frac > best$frac || (frac == best$frac && paths_dt$path_prob[i] > best$prob)) {
            best <- list(idx = i, frac = frac, prob = paths_dt$path_prob[i], path_sig = paths_dt$path_sig[i])
        }
    }
    list(found = FALSE, best = best)
}
