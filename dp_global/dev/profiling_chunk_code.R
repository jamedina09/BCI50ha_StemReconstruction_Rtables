############################################################
### profiling_chunk_code.R — profile the chunked DP step
############################################################
# Run with:
#   Rscript --vanilla dp_global/dev/profiling_chunk_code.R

get_script_dir <- function() {
    cmd <- commandArgs(trailingOnly = FALSE)
    file_arg <- sub("^--file=", "", cmd[grep("^--file=", cmd)])
    if (length(file_arg) == 1L && nzchar(file_arg)) {
        return(dirname(normalizePath(file_arg)))
    }
    if (!is.null(sys.frames()[[1L]]$ofile)) {
        return(dirname(normalizePath(sys.frames()[[1L]]$ofile)))
    }
    getwd()
}

find_project_root <- function(start_dir) {
    d <- normalizePath(start_dir)
    for (i in 0:6) {
        cand <- d
        if (dir.exists(file.path(cand, "STEM_IDENTIFICATION_TEST"))) {
            return(cand)
        }
        d2 <- dirname(d)
        if (identical(d2, d)) break
        d <- d2
    }
    stop("Could not find project root containing STEM_IDENTIFICATION_TEST/ starting from: ", start_dir)
}

script_dir <- get_script_dir()
root_dir <- find_project_root(script_dir)
setwd(file.path(root_dir, "STEM_IDENTIFICATION_TEST"))

if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Please install the 'data.table' package to run this profiling script.")
}
library(data.table)

if (!requireNamespace("here", quietly = TRUE)) {
    stop("Please install the 'here' package to run this profiling script.")
}
library(here)

############################################################
### 1) Minimal config for chunk profiling
############################################################

input_file <- here("data_simulation", "data", "simulated_data_1.csv")

# Chunking controls (small test chunk)
DP_CHUNK_SIZE <- 2L
DP_CHUNK_START <- 1L
DP_CHUNK_END <- 1L
DP_CHUNK_RESUME <- FALSE
DP_CHUNK_OVERWRITE <- TRUE
WRITE_DP_PDF_PER_CHUNK <- TRUE
WRITE_DP_PDF <- TRUE
WRITE_DP_FEATHER <- TRUE
WRITE_DP_RDS <- TRUE
WRITE_DP_CSV <- FALSE  # avoid appending to global CSV during profiling

# Posterior sampling (enabled for profiling). Can be overridden with env var POSTERIOR_SAMPLES
POSTERIOR_SAMPLES <- suppressWarnings(as.integer(Sys.getenv("POSTERIOR_SAMPLES", "200")))
if (is.na(POSTERIOR_SAMPLES) || POSTERIOR_SAMPLES < 0) POSTERIOR_SAMPLES <- 200L
POSTERIOR_SAMPLES_FORMAT <- Sys.getenv("POSTERIOR_SAMPLES_FORMAT", "feather") # 'rds','feather','csv'
POSTERIOR_SAMPLES_PATH <- NULL
POSTERIOR_SAMPLE_SEED <- if (nzchar(Sys.getenv("POSTERIOR_SAMPLE_SEED", ""))) as.integer(Sys.getenv("POSTERIOR_SAMPLE_SEED")) else NULL


# DP controls
which_tag <- 1L
anchor_start_census <- 7L
min_annual_growth <- -0.5
max_annual_growth <- 5
DP_POSTERIOR_TEMPERATURE <- 1.0
DP_POSTERIOR_TOP_K <- 2L

# Deterministic tie-break weight
DP_EPS_TIEBREAK <- suppressWarnings(as.numeric(Sys.getenv("DP_EPS_TIEBREAK", "1e-6")))
if (!is.finite(DP_EPS_TIEBREAK) || is.na(DP_EPS_TIEBREAK) || DP_EPS_TIEBREAK < 0) {
    DP_EPS_TIEBREAK <- 1e-6
}

MANUAL_CORES <- TRUE
MANUAL_CORES_VALUE <- 1L
MC_CORES <- as.integer(MANUAL_CORES_VALUE)

dp_max_tracks <- NULL

dp_max_states <- 40000L

dp_slack_tracks <- 1L

# Profiling controls
RUN_PROFILE <- TRUE
DP_VERBOSE <- FALSE

############################################################
### 2) Source helper code (functions only) — reuse the DP helpers
############################################################
# The chunk driver (main_cpp_chunk.R) runs as a standalone script; here we source the
# core functions from the project (dp_global_main.R) so we can exercise a chunk loop
source(here("dp_global", "R", "dp_global_main.R"))

############################################################
### 3) Prepare data and biologic parameters (outside profiling)
############################################################

# Minimal helpers used by the profiling script (kept local to avoid depending on
# chunk-driver's top-level script state)
ensure_species_column <- function(x) {
    if ("species" %in% names(x)) {
        x[, species := as.character(species)]
        return(x)
    }
    candidates <- c("Species", "SP", "sp", "spcode", "sp_code", "taxon", "Taxon", "spname", "Sp")
    found <- candidates[candidates %in% names(x)]
    if (length(found) < 1L) {
        stop("No species column found. Add a 'species' column or provide a standard name like 'Species'.")
    }
    x[, species := as.character(get(found[[1L]]))]
    x
}

get_nested_numeric <- function(x, expr, fallback = NULL) {
    v <- tryCatch(eval(expr, envir = x), error = function(e) NULL)
    if (!is.null(v) && length(v) == 1L && is.finite(v)) {
        return(as.numeric(v))
    }
    fallback
}

get_growth_mu_const <- function(growth_list) {
    if (!is.null(growth_list$alpha)) {
        return(growth_list$alpha)
    }
    growth_list$mu
}

attach_bio_columns <- function(xrun, bio_pars) {
    xrun[, `:=`(
        Bio_Mu_Growth = get_growth_mu_const(bio_pars[[species]]$growth),
        Bio_Gamma_Growth = {
            g <- bio_pars[[species]]$growth
            if (!is.null(g$gamma)) g$gamma else 0
        },
        Bio_Sigma0_Growth = bio_pars[[species]]$growth$sigma0,
        Bio_Sigma1_Growth = bio_pars[[species]]$growth$sigma1,
        Bio_H0_Mortality = bio_pars[[species]]$mortality$h0,
        Bio_Beta_Mortality = bio_pars[[species]]$mortality$beta,
        Bio_Recruit_Meanlog = bio_pars[[species]]$recruitment$meanlog,
        Bio_Recruit_Sdlog = bio_pars[[species]]$recruitment$sdlog,
        Bio_Recruit_MaxDBH_unit = bio_pars[[species]]$recruitment$recruit_max_dbh,
        Bio_Recruitment_lambda = bio_pars[[species]]$recruitment$lambda,
        Bio_Max_Shrink = {
            sh <- bio_pars[[species]]$shrinkage
            get_nested_numeric(sh, quote(guardrails$hard$value), fallback = sh$max_shrink)
        },
        Bio_K_Shrink = {
            sh <- bio_pars[[species]]$shrinkage
            get_nested_numeric(sh, quote(penalties$soft$k), fallback = sh$k_shrink)
        },
        Bio_Max_Growth = {
            g <- bio_pars[[species]]$growth
            get_nested_numeric(g, quote(guardrails$hard$value), fallback = g$max_growth)
        },
        Bio_Max_Growth_Soft = {
            g <- bio_pars[[species]]$growth
            get_nested_numeric(g, quote(guardrails$soft$value), fallback = g$max_growth_soft)
        },
        Bio_K_Growth = {
            g <- bio_pars[[species]]$growth
            get_nested_numeric(g, quote(penalties$soft$k), fallback = g$k_growth)
        }
    ), by = species]
    xrun
}

message("[profiling_chunk_code] Loading data: ", input_file)
xraw <- data.table::fread(input_file)
# force single-species for profiling convenience
xraw[, species := "all"]

# Ensure species col
if (!"species" %in% names(xraw)) {
    xraw <- ensure_species_column(xraw)
}

xrun <- data.table::copy(xraw)

message("[profiling_chunk_code] Estimating bio parameters (outside profiling)")
bio_pars <- list()
for (sp in unique(xrun$species)) {
    bio_pars[[sp]] <- estimate_bio_pars(
        xrun[species == sp],
        use_measurement_error = TRUE,
        max_shrink_source = "data",
        max_shrink_fixed = -1,
        k_shrink_source = "fixed",
        k_shrink_fixed = 25,
        k_growth_source = "data",
        k_growth_fixed = 1,
        max_growth_source = "data",
        max_growth_fixed = 7.5
    )
}

xrun <- attach_bio_columns(xrun, bio_pars)

# auto_dp_max_tracks helper (copied locally)
auto_dp_max_tracks <- function(xrun) {
    max_obs_any_tag_census <- xrun[
        CensusID <= anchor_start_census & !is.na(DBH),
        .N,
        by = .(Tag, CensusID)
    ][, max(N, na.rm = TRUE)]
    if (!is.finite(max_obs_any_tag_census)) max_obs_any_tag_census <- 0L
    as.integer(max_obs_any_tag_census + 1L)
}

dp_max_tracks_local <- if (is.null(dp_max_tracks)) auto_dp_max_tracks(xrun) else as.integer(dp_max_tracks)

groups <- unique(xrun[, .(Tag, species)])
setorder(groups, Tag, species)
group_idx <- seq_len(nrow(groups))
chunks <- split(group_idx, ceiling(seq_along(group_idx) / as.integer(DP_CHUNK_SIZE)))

# Choose the first chunk for a short profiling run
ci <- if (length(chunks) >= DP_CHUNK_START) DP_CHUNK_START else 1L
if (ci > length(chunks)) stop("DP_CHUNK_START exceeds number of chunks")

groups_ci <- groups[chunks[[ci]]]
message("[profiling_chunk_code] Profiling chunk ", ci, " with ", nrow(groups_ci), " groups")

tmp_out_dir <- file.path(here("dp_global", "output"), paste0("profile_chunk_", format(Sys.time(), "%Y%m%d_%H%M%S")))
if (!dir.exists(tmp_out_dir)) dir.create(tmp_out_dir, recursive = TRUE)
message("[profiling_chunk_code] Using out_dir: ", tmp_out_dir)

# Posterior bin helper (local)
ADD_DP_POSTERIOR_BINS <- TRUE
maybe_add_posterior_bins <- function(out) {
    if (isTRUE(ADD_DP_POSTERIOR_BINS) && !is.null(out)) {
        out <- add_dp_posterior_bins(
            out,
            confident_prob = 0.95,
            unlinked_prob = 0.5,
            use_reconstructed_prob = TRUE,
            out_col = "DP_PosteriorBin"
        )
    }
    out
}

# Make a small wrapper to run the chunk (so we can profile it)
run_one_chunk <- function() {
    res <- parallel::mclapply(seq_len(nrow(groups_ci)), function(j) {
        data.table::setDTthreads(1L)
        g <- groups_ci[j]
        dtg <- xrun[Tag == g$Tag & species == g$species]
        match_stems_dp_global_backward_marginals_batch(
            tree_data = data.table::copy(dtg),
            min_growth = min_annual_growth,
            max_growth = max_annual_growth,
            anchor_start = anchor_start_census,
            max_tracks = dp_max_tracks_local,
            max_states = dp_max_states,
            slack_tracks = dp_slack_tracks,
            slack_require_anchor_recruitable = TRUE,
            temperature = DP_POSTERIOR_TEMPERATURE,
            posterior_top_k = DP_POSTERIOR_TOP_K,
            eps_tiebreak = DP_EPS_TIEBREAK,
            use_measurement_error = TRUE,
            verbose = isTRUE(DP_VERBOSE),
            posterior_samples = POSTERIOR_SAMPLES,
            prune_hard = TRUE,
            prune_min_growth = min_annual_growth,
            prune_max_growth = max_annual_growth,
            prune_use_bio_bounds = FALSE,
            prune_recruit_max_dbh = 7.5*5,
            prune_use_bio_recruit = FALSE
        )
    }, mc.cores = MC_CORES)

    out_chunk <- data.table::rbindlist(res, use.names = TRUE, fill = TRUE)
    if (nrow(out_chunk) > 0L) {
        out_chunk[, DP_Chunk := ci]
        out_chunk <- maybe_add_posterior_bins(out_chunk)

        # ensure posterior samples path if sampling requested
        if (!is.null(POSTERIOR_SAMPLES) && as.integer(POSTERIOR_SAMPLES) > 0L) {
            if (is.null(POSTERIOR_SAMPLES_PATH) || !nzchar(POSTERIOR_SAMPLES_PATH)) {
                POSTERIOR_SAMPLES_PATH <<- file.path(tmp_out_dir, "posteriors")
            }
            if (!dir.exists(POSTERIOR_SAMPLES_PATH)) dir.create(POSTERIOR_SAMPLES_PATH, recursive = TRUE)
        }

        # measure DP -> produce out_chunk time (already measured by outer profiler); now benchmark exports
        export_timings <- list()

        # Feather export
        if (isTRUE(WRITE_DP_FEATHER) && requireNamespace("arrow", quietly = TRUE)) {
            feather_path <- file.path(tmp_out_dir, sprintf("stem_reconstruction_dp_global_rcpp_chunk_%03d.feather", ci))
            t0 <- proc.time(); arrow::write_feather(out_chunk, feather_path); t1 <- proc.time()
            export_timings$feather <- as.numeric((t1 - t0)["elapsed"])
            export_timings$feather_size <- file.size(feather_path)
        } else {
            export_timings$feather <- NA_real_
            export_timings$feather_size <- NA_integer_
        }

        # RDS export
        if (isTRUE(WRITE_DP_RDS)) {
            rds_path <- file.path(tmp_out_dir, sprintf("stem_reconstruction_dp_global_rcpp_chunk_%03d.rds", ci))
            t0 <- proc.time(); saveRDS(out_chunk, file = rds_path); t1 <- proc.time()
            export_timings$rds <- as.numeric((t1 - t0)["elapsed"])
            export_timings$rds_size <- file.size(rds_path)
        } else {
            export_timings$rds <- NA_real_
            export_timings$rds_size <- NA_integer_
        }

        # CSV export (simulated combined CSV write for benchmarking)
        if (isTRUE(WRITE_DP_CSV)) {
            csv_path <- file.path(tmp_out_dir, sprintf("stem_reconstruction_dp_global_rcpp_chunk_%03d.csv", ci))
            t0 <- proc.time(); data.table::fwrite(out_chunk, csv_path); t1 <- proc.time()
            export_timings$csv <- as.numeric((t1 - t0)["elapsed"])
            export_timings$csv_size <- file.size(csv_path)
        } else {
            export_timings$csv <- NA_real_
            export_timings$csv_size <- NA_integer_
        }

        # PDF export
        if (isTRUE(WRITE_DP_PDF_PER_CHUNK) && isTRUE(WRITE_DP_PDF)) {
            pdf_path <- file.path(tmp_out_dir, sprintf("stem_reconstruction_dp_global_rcpp_chunk_%03d.pdf", ci))
            t0 <- proc.time();
            tryCatch({
                plot_tag_to_pdf(out_chunk, pdf_file = pdf_path, include_reference = TRUE)
                t1 <- proc.time()
                export_timings$pdf <- as.numeric((t1 - t0)["elapsed"])
                export_timings$pdf_size <- file.size(pdf_path)
            }, error = function(e) {
                t1 <- proc.time()
                export_timings$pdf <- as.numeric((t1 - t0)["elapsed"])
                export_timings$pdf_size <- NA_integer_
                message("[profiling_chunk_code] PDF generation failed: ", conditionMessage(e))
            })
        } else {
            export_timings$pdf <- NA_real_
            export_timings$pdf_size <- NA_integer_
        }

        # Capture posterior samples artifacts if generated by DP (scan dir)
        if (!is.null(POSTERIOR_SAMPLES) && as.integer(POSTERIOR_SAMPLES) > 0L && dir.exists(POSTERIOR_SAMPLES_PATH)) {
            p_files <- list.files(POSTERIOR_SAMPLES_PATH, full.names = TRUE)
            export_timings$posterior_files_count <- length(p_files)
            export_timings$posterior_files_total_size <- if (length(p_files) > 0) sum(file.size(p_files)) else 0L
        } else {
            export_timings$posterior_files_count <- NA_integer_
            export_timings$posterior_files_total_size <- NA_integer_
        }

        export_timings
    } else {
        NULL
    }
    invisible(out_chunk)
}

prof_file <- file.path(here("dp_global", "dev"), "dp_global_chunk_cpp.prof")
message("[profiling_chunk_code] Profiling chunk run; writing profile to: ", prof_file)

if (!isTRUE(RUN_PROFILE)) {
    message("[profiling_chunk_code] RUN_PROFILE=FALSE; performing a dry-run chunk execution")
    out <- run_one_chunk()
} else {
    Rprof(prof_file, interval = 0.01, memory.profiling = TRUE)
    out <- run_one_chunk()
    Rprof(NULL)
    message("[profiling_chunk_code] Profiling complete. Top hotspots:")
    s <- summaryRprof(prof_file)
    print(utils::head(s$by.self, 25))
    cat("\n---- by.total ----\n")
    print(utils::head(s$by.total, 25))
}

message("[profiling_chunk_code] Chunk profiling finished; artifacts written to: ", tmp_out_dir)

invisible(out)
