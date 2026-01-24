############################################################
### profiling_code.R — profile and benchmark the DP step
############################################################
# Run with:
#   Rscript --vanilla dp_global/dev/profiling_code.R

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
### 1) Minimal config for profiling
############################################################

input_file <- here("data_simulation", "data", "simulated_data_1.csv")

# /Users/medinaja/Library/CloudStorage/OneDrive-SmithsonianInstitution/STRI/STEM_IDENTIFICATION_TEST/data_simulation/data/simulated_data_1.csv

which_tag <- 1L
anchor_start_census <- 7L

# Candidate pruning bounds used by add_constraint_violation()
min_annual_growth <- -0.5
max_annual_growth <- 5

# DP controls
DP_POSTERIOR_TEMPERATURE <- 1.0
DP_POSTERIOR_TOP_K <- 2L

# Posterior sampling (env overrides)
POSTERIOR_SAMPLES <- suppressWarnings(as.integer(Sys.getenv("POSTERIOR_SAMPLES", "0")))
if (is.na(POSTERIOR_SAMPLES) || POSTERIOR_SAMPLES < 0) POSTERIOR_SAMPLES <- 0L
POSTERIOR_SAMPLES_FORMAT <- Sys.getenv("POSTERIOR_SAMPLES_FORMAT", "feather") # 'rds','feather','csv'
POSTERIOR_SAMPLES_PATH <- Sys.getenv("POSTERIOR_SAMPLES_PATH", "")
POSTERIOR_SAMPLE_SEED <- if (nzchar(Sys.getenv("POSTERIOR_SAMPLE_SEED", ""))) as.integer(Sys.getenv("POSTERIOR_SAMPLE_SEED")) else NULL

# Deterministic tie-break weight (performance knob).
# - Set to 0 to disable rank()-based tie-breaking (often much faster).
# - Override without editing:
#     DP_EPS_TIEBREAK=0 Rscript --vanilla STEM_IDENTIFICATION_TEST/dp_global/dev/profiling_code.R
DP_EPS_TIEBREAK <- suppressWarnings(as.numeric(Sys.getenv("DP_EPS_TIEBREAK", "1e-6")))
if (!is.finite(DP_EPS_TIEBREAK) || is.na(DP_EPS_TIEBREAK) || DP_EPS_TIEBREAK < 0) {
    DP_EPS_TIEBREAK <- 1e-6
}

dp_max_tracks <- NULL

dp_max_states <- 40000L

dp_slack_tracks <- 1L

# Profiling controls (minimal for CPP-only profiling)
RUN_PROFILE <- TRUE
# Turn off chatty printing during profiling (important)
DP_VERBOSE <- TRUE
# Benchmarking flags removed for a simple single-run profiler

# Export controls (env overrides accepted: 1/0, TRUE/FALSE)
env_flag <- function(name, default = FALSE) {
    v <- Sys.getenv(name, "")
    if (v == "") return(default)
    tolower(v) %in% c("1", "true", "t", "yes", "y")
}

WRITE_DP_FEATHER <- env_flag("WRITE_DP_FEATHER", TRUE)
WRITE_DP_RDS <- env_flag("WRITE_DP_RDS", TRUE)
WRITE_DP_CSV <- env_flag("WRITE_DP_CSV", TRUE)
WRITE_DP_PDF <- env_flag("WRITE_DP_PDF", TRUE)
# Control whether the DP code writes posterior summaries to disk
WRITE_POSTERIORS <- env_flag("WRITE_POSTERIORS", TRUE)


############################################################
### 2) Source DP code
############################################################

# source(here("dp_global", "R", "dp_global_biol.R"))
source(here("dp_global", "R", "dp_global_main.R"))

############################################################
### 3) Minimal data prep + Bio_* columns (done OUTSIDE profiling)
############################################################

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

auto_dp_max_tracks <- function(xrun) {
    max_obs_any_tag_census <- xrun[
        CensusID <= anchor_start_census & !is.na(DBH),
        .N,
        by = .(Tag, CensusID)
    ][, max(N, na.rm = TRUE)]
    if (!is.finite(max_obs_any_tag_census)) max_obs_any_tag_census <- 0L
    as.integer(max_obs_any_tag_census + 1L)
}

message("[profiling_code] Loading data: ", input_file)
xraw <- data.table::fread(input_file)
xraw[, species := "all"]
xraw <- ensure_species_column(xraw)
xrun <- data.table::copy(xraw)

message("[profiling_code] Estimating bio parameters (outside profiling)")
bio_pars <- list()
for (sp in unique(xrun$species)) {
    bio_pars[[sp]] <- estimate_bio_pars(
        xrun[species == sp],
        mortality_start = c(log(0.01), 0),
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

dp_max_tracks_local <- if (is.null(dp_max_tracks)) auto_dp_max_tracks(xrun) else as.integer(dp_max_tracks)

sp0 <- unique(xrun$species)
sp0 <- sp0[!is.na(sp0) & nzchar(sp0)]
sp0 <- sp0[[1L]]

dtg <- xrun[Tag == which_tag & species == sp0]
if (nrow(dtg) < 1L) stop("No rows found for which_tag=", which_tag, " and species=", sp0)

############################################################
### 4) Profile the DP call
############################################################

run_dp_cpp <- function(dt, out_dir = NULL) {
    # If caller provided an out_dir for posterior writing, map that to posterior_samples_path
    posterior_samples_path <- if (!is.null(out_dir) && nzchar(out_dir)) out_dir else if (nzchar(POSTERIOR_SAMPLES_PATH)) POSTERIOR_SAMPLES_PATH else NULL

    match_stems_dp_global_backward_marginals_batch(
        tree_data = data.table::copy(dt),
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
        posterior_samples_format = POSTERIOR_SAMPLES_FORMAT,
        posterior_samples_path = posterior_samples_path,
        posterior_sample_seed = POSTERIOR_SAMPLE_SEED,
        prune_hard = TRUE,
        prune_min_growth = min_annual_growth * 1,
        prune_max_growth = max_annual_growth * 1,
        prune_use_bio_bounds = FALSE,
        prune_recruit_max_dbh = 7.5*5,
        prune_use_bio_recruit = FALSE
    )
}

if (!isTRUE(RUN_PROFILE)) {
    message("[profiling_code] RUN_PROFILE=FALSE; skipping Rprof.")
    quit(save = "no", status = 0)
}

# Create a per-run temporary output directory for profiling artifacts
tmp_out_dir <- file.path(here("dp_global", "output"), paste0("profile_", format(Sys.time(), "%Y%m%d_%H%M%S")))
if (!dir.exists(tmp_out_dir)) dir.create(tmp_out_dir, recursive = TRUE)

prof_file <- file.path(here("dp_global", "dev"), "dp_global_cpp.prof")
message("[profiling_code] Profiling DP run (Tag=", which_tag, ", species=", sp0, ")")
message("[profiling_code] Writing profile to: ", prof_file)
message("[profiling_code] Using out_dir: ", tmp_out_dir)
message("[profiling_code] Export flags - FEATHER=", WRITE_DP_FEATHER, ", RDS=", WRITE_DP_RDS, ", CSV=", WRITE_DP_CSV, ", PDF=", WRITE_DP_PDF, ", POSTERIORS=", WRITE_POSTERIORS)

# Ensure POSTERIOR_SAMPLES_PATH maps to our tmp_out_dir when sampling requested
effective_posterior_samples <- as.integer(POSTERIOR_SAMPLES)
if (!isTRUE(WRITE_POSTERIORS)) {
    effective_posterior_samples <- 0L
}
if (!is.null(effective_posterior_samples) && effective_posterior_samples > 0L) {
    if (!nzchar(POSTERIOR_SAMPLES_PATH)) POSTERIOR_SAMPLES_PATH <- tmp_out_dir
    if (!dir.exists(file.path(POSTERIOR_SAMPLES_PATH, "posteriors"))) dir.create(file.path(POSTERIOR_SAMPLES_PATH, "posteriors"), recursive = TRUE)
}

# Temporarily override POSTERIOR_SAMPLES used by run_dp_cpp
old_POSTERIOR_SAMPLES <- POSTERIOR_SAMPLES
POSTERIOR_SAMPLES <- effective_posterior_samples
Rprof(prof_file, interval = 0.01, memory.profiling = TRUE)
out_prof <- run_dp_cpp(dtg, out_dir = tmp_out_dir)
Rprof(NULL)
POSTERIOR_SAMPLES <- old_POSTERIOR_SAMPLES


message("[profiling_code] Done. Top hotspots:\n")
s <- summaryRprof(prof_file)
print(utils::head(s$by.self, 25))
cat("\n---- by.total ----\n")
print(utils::head(s$by.total, 25))

# Export write timings and sizes for the main DP output (respect WRITE_* flags)
export_timings <- list()
if (!is.null(out_prof) && data.table::is.data.table(out_prof) && nrow(out_prof) > 0L) {
    # Feather
    if (isTRUE(WRITE_DP_FEATHER) && requireNamespace("arrow", quietly = TRUE)) {
        feather_path <- file.path(tmp_out_dir, "stem_reconstruction_dp_global_rcpp_profile.feather")
        t0 <- proc.time(); arrow::write_feather(out_prof, feather_path); t1 <- proc.time()
        export_timings$feather <- as.numeric((t1 - t0)["elapsed"])
        export_timings$feather_size <- file.size(feather_path)
    } else {
        export_timings$feather <- NA_real_
        export_timings$feather_size <- NA_integer_
    }

    # RDS
    if (isTRUE(WRITE_DP_RDS)) {
        rds_path <- file.path(tmp_out_dir, "stem_reconstruction_dp_global_rcpp_profile.rds")
        t0 <- proc.time(); saveRDS(out_prof, file = rds_path); t1 <- proc.time()
        export_timings$rds <- as.numeric((t1 - t0)["elapsed"])
        export_timings$rds_size <- file.size(rds_path)
    } else {
        export_timings$rds <- NA_real_
        export_timings$rds_size <- NA_integer_
    }

    # CSV
    if (isTRUE(WRITE_DP_CSV)) {
        csv_path <- file.path(tmp_out_dir, "stem_reconstruction_dp_global_rcpp_profile.csv")
        t0 <- proc.time(); data.table::fwrite(out_prof, csv_path); t1 <- proc.time()
        export_timings$csv <- as.numeric((t1 - t0)["elapsed"])
        export_timings$csv_size <- file.size(csv_path)
    } else {
        export_timings$csv <- NA_real_
        export_timings$csv_size <- NA_integer_
    }

    # PDF
    if (isTRUE(WRITE_DP_PDF)) {
        pdf_path <- file.path(tmp_out_dir, "stem_reconstruction_dp_global_rcpp_profile.pdf")
        t0 <- proc.time();
        tryCatch({
            plot_tag_to_pdf(out_prof, pdf_file = pdf_path, include_reference = TRUE)
            t1 <- proc.time()
            export_timings$pdf <- as.numeric((t1 - t0)["elapsed"])
            export_timings$pdf_size <- file.size(pdf_path)
        }, error = function(e) {
            t1 <- proc.time()
            export_timings$pdf <- as.numeric((t1 - t0)["elapsed"])
            export_timings$pdf_size <- NA_integer_
            message("[profiling_code] PDF generation failed: ", conditionMessage(e))
        })
    } else {
        export_timings$pdf <- NA_real_
        export_timings$pdf_size <- NA_integer_
    }
} else {
    message("[profiling_code] DP returned no rows; skipping output exports.")
}

# Posterior file stats
posterior_files_count <- NA_integer_
posterior_files_total_size <- NA_integer_
if (isTRUE(WRITE_POSTERIORS) && !is.null(POSTERIOR_SAMPLES) && as.integer(POSTERIOR_SAMPLES) > 0L && dir.exists(file.path(POSTERIOR_SAMPLES_PATH, "posteriors"))) {
    p_files <- list.files(file.path(POSTERIOR_SAMPLES_PATH, "posteriors"), full.names = TRUE)
    posterior_files_count <- length(p_files)
    posterior_files_total_size <- if (length(p_files) > 0) sum(file.size(p_files)) else 0L
} else if (!isTRUE(WRITE_POSTERIORS)) {
    posterior_files_count <- NA_integer_
    posterior_files_total_size <- NA_integer_
}

cat("\nExport timings (seconds) and sizes (bytes):\n")
print(export_timings)
cat(sprintf("posterior_files_count: %s\nposterior_files_total_size_bytes: %s\n", posterior_files_count, posterior_files_total_size))

# Inspect DP sampling & compute profiles attached to returned object (if any)
if (!is.null(out_prof)) {
    sp <- attr(out_prof, "DP_Sampling_Profile", exact = TRUE)
    cp <- attr(out_prof, "DP_Compute_Profile", exact = TRUE)
    if (!is.null(sp)) {
        cat("\nDP_Sampling_Profile:\n")
        print(sp)
    }
    if (!is.null(cp)) {
        cat("\nDP_Compute_Profile:\n")
        print(cp)
    }
}

message("[profiling_code] Profiling artifacts written to: ", tmp_out_dir)

invisible(out_prof)

# POSTERIOR_SAMPLES=20 POSTERIOR_SAMPLES_FORMAT=feather WRITE_DP_FEATHER=TRUE WRITE_DP_RDS=FALSE WRITE_DP_CSV=FALSE WRITE_DP_PDF=FALSE WRITE_POSTERIORS=FALSE Rscript --vanilla dp_global/dev/profiling_code.R