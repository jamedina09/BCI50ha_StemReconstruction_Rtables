########################################################################
### main_cpp_chunk.R — dp_global driver (chunked version)
###
### This file is a chunk-focused variant of `main_cpp.R`. It processes
### Tag+species groups in configurable chunks and writes per-chunk outputs
### (CSV, RDS, Feather) to avoid high memory use.
###
### Usage: same CLI as `main_cpp.R` but supports chunking options:
###  --DP_CHUNK_SIZE, --DP_CHUNK_RESUME, --DP_CHUNK_OVERWRITE, --WRITE_DP_FEATHER
########################################################################

# We largely reuse the contents of main_cpp.R but keep chunking code active
# and leave the other behavior identical.

# NOTE: This file was generated automatically as a chunked variant of
# `main_cpp.R`. Keep it in sync if you update core functionality.

source(here::here("dp_global", "scripts", "main_cpp.R"))

# The above sources the original script and defines all helpers. We now
# override the chunk-specific defaults to ensure chunked behavior.
DP_CHUNK_SIZE <- get0("DP_CHUNK_SIZE", ifnotfound = 80L)
DP_CHUNK_RESUME <- get0("DP_CHUNK_RESUME", ifnotfound = TRUE)
DP_CHUNK_OVERWRITE <- get0("DP_CHUNK_OVERWRITE", ifnotfound = FALSE)
WRITE_DP_FEATHER <- get0("WRITE_DP_FEATHER", ifnotfound = FALSE)

# Export a wrapper that runs the chunk-based branch only.
run_main_chunked <- function() {
    # Copy-paste the chunk-specific part of run_main from main_cpp.R
    ensure_dir(out_dir)

    # Create posteriors directory if requested (POSTERIOR_SAMPLES > 0)
    if (!is.null(POSTERIOR_SAMPLES) && as.integer(POSTERIOR_SAMPLES) > 0L && !is.null(POSTERIOR_SAMPLES_PATH) && nzchar(POSTERIOR_SAMPLES_PATH)) {
        ensure_dir(POSTERIOR_SAMPLES_PATH)
    }

    tryCatch({
        writeLines(as.character(Sys.time()), con = file.path(out_dir, "run_started.txt"))
    }, error = function(e) {
        message("[dp_global main_cpp_chunk.R] Warning writing run_started marker: ", conditionMessage(e))
    })

    # 5.1 Load data
    xraw <- data.table::fread(input_file)
    xraw <- ensure_species_column(xraw)
    xrun <- data.table::copy(xraw)

    # 5.2 Estimate biological parameters (reuse same workflow)
    bio_pars <- list()
    for (sp in unique(xrun$species)) {
        bio_pars[[sp]] <- estimate_bio_pars(
            xrun[species == sp],
            use_measurement_error = isTRUE(USE_MEASUREMENT_ERROR),
            max_shrink_source = MAX_SHRINK_HARD_SOURCE,
            max_shrink_fixed = MAX_SHRINK_FIXED,
            k_shrink_source = K_SHRINK_SOURCE,
            k_shrink_fixed = K_SHRINK_FIXED,
            k_growth_source = K_GROWTH_SOURCE,
            k_growth_fixed = K_GROWTH_FIXED,
            max_growth_source = MAX_GROWTH_HARD_SOURCE,
            max_growth_fixed = MAX_GROWTH_FIXED,
            shrink_hard_prob = 1e-4,
            shrink_data_quantile = 0.001,
            growth_hard_prob = 1e-4,
            growth_data_quantile = 0.999,
            growth_soft_quantile = 0.99,
            recruit_max_quantile = 0.999,
            recruit_max_source = get0("RECRUIT_MAX_SOURCE", ifnotfound = "data"),
            recruit_max_fixed = as.numeric(get0("RECRUIT_MAX_FIXED", ifnotfound = 5))
        )
    }

    # 5.3 Attach Bio columns and compute dp_max_tracks_local
    xrun <- attach_bio_columns(xrun, bio_pars)
    dp_max_tracks_local <- if (is.null(dp_max_tracks)) auto_dp_max_tracks(xrun) else as.integer(dp_max_tracks)
    dp_max_tracks_local <- as.integer(dp_max_tracks_local)

    # 5.4 Chunk processing
    if (!requireNamespace("parallel", quietly = TRUE)) stop("Package not available: parallel.")
    groups <- unique(xrun[, .(Tag, species)])
    data.table::setorder(groups, Tag, species)
    group_idx <- seq_len(nrow(groups))
    chunks <- split(group_idx, ceiling(seq_along(group_idx) / as.integer(DP_CHUNK_SIZE)))
    first_chunk <- !file.exists(DP_CSV_FILE)

    for (ci in seq_along(chunks)) {
        chunk_rds <- file.path(out_dir, sprintf("stem_reconstruction_dp_global_rcpp_chunk_%03d.rds", ci))
        if (isTRUE(DP_CHUNK_RESUME) && file.exists(chunk_rds) && !isTRUE(DP_CHUNK_OVERWRITE)) {
            message(sprintf("[dp_global main_cpp_chunk.R] Skipping chunk %d/%d — chunk RDS exists (resume enabled)", ci, length(chunks)))
            first_chunk <- !file.exists(DP_CSV_FILE)
            next
        }

        inds <- chunks[[ci]]
        groups_ci <- groups[inds]
        message(sprintf("[dp_global main_cpp_chunk.R] Chunk %d/%d — %d groups", ci, length(chunks), nrow(groups_ci)))

        res <- parallel::mclapply(seq_len(nrow(groups_ci)), function(j) {
            data.table::setDTthreads(1L)
            g <- groups_ci[j]
            dtg <- xrun[Tag == g$Tag & species == g$species]
            run_dp_one_group(dtg, dp_max_tracks = dp_max_tracks_local)
        }, mc.cores = MC_CORES)

        out_chunk <- data.table::rbindlist(res, use.names = TRUE, fill = TRUE)
        if (nrow(out_chunk) > 0L) {
            out_chunk[, DP_Chunk := ci]
            out_chunk <- maybe_add_posterior_bins(out_chunk)
            out_chunk[, out_dir := basename(out_dir)]

            if (isTRUE(WRITE_DP_CSV)) data.table::fwrite(out_chunk, file = DP_CSV_FILE, append = !first_chunk)
            if (isTRUE(WRITE_DP_FEATHER) && requireNamespace("arrow", quietly = TRUE)) arrow::write_feather(out_chunk, file.path(out_dir, sprintf("stem_reconstruction_dp_global_rcpp_chunk_%03d.feather", ci)))
            if (isTRUE(WRITE_DP_RDS)) saveRDS(out_chunk, file = chunk_rds)
        } else {
            message(sprintf("[dp_global main_cpp_chunk.R] Chunk %d returned no rows.", ci))
            if (isTRUE(WRITE_DP_RDS)) saveRDS(out_chunk, file = chunk_rds)
        }

        first_chunk <- FALSE
        rm(res, out_chunk, groups_ci)
        invisible(gc())
    }

    # Finished
    tryCatch({ writeLines(as.character(Sys.time()), con = file.path(out_dir, "run_finished.txt")) }, error = function(e) message("[dp_global main_cpp_chunk.R] Warning writing run_finished: ", conditionMessage(e)))

    invisible(list(xrun = xrun, bio_pars = bio_pars))
}

# If called directly, run the chunked driver
if (sys.nframe() == 0L) {
    message("[dp_global main_cpp_chunk.R] Starting chunked run_main_chunked()")
    run_main_chunked()
}
