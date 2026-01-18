############################################################
### profiling_code.R — profile and benchmark the DP step
############################################################
# Run with:
#   Rscript --vanilla dp_global/dev/profiling_code.R
#
# Notes:
# - This script intentionally does NOT source main.R because main.R may rm(list=ls()).

# PROFILE_VARIANT=scalar Rscript --vanilla dp_global/dev/profiling_code.R
# PROFILE_VARIANT=batch Rscript --vanilla dp_global/dev/profiling_code.R
# PROFILE_VARIANT=cpp Rscript --vanilla dp_global/dev/profiling_code.R

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

if( !requireNamespace("here", quietly = TRUE)) {
    stop("Please install the 'here' package to run this profiling script.")
}
library(here)

############################################################
### 1) Minimal config for profiling
############################################################

input_file <- here("data_simulation", "data", "simulation_legacy_backup", "simulated_data_one_species.csv")

which_tag <- 1L
anchor_start_census <- 7L
census_interval_years <- 5

# Candidate pruning bounds used by add_constraint_violation()
min_annual_growth <- -0.5
max_annual_growth <- 5

# DP controls
RUN_DP_MARGINALS <- TRUE
DP_POSTERIOR_TEMPERATURE <- 1.0
DP_POSTERIOR_TOP_K <- 2L

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

# Benchmark / profiling controls
RUN_BENCHMARK <- TRUE
BENCHMARK_WARMUP <- TRUE
RUN_PROFILE <- TRUE

# Which variant to profile: "scalar", "batch", or "cpp"
# Override without editing: PROFILE_VARIANT=cpp Rscript --vanilla STEM_IDENTIFICATION_TEST/dp_global/dev/profiling_code.R
PROFILE_VARIANT <- Sys.getenv("PROFILE_VARIANT", "batch")
if (!(tolower(PROFILE_VARIANT) %in% c("scalar", "batch", "cpp"))) {
    stop("PROFILE_VARIANT must be 'scalar', 'batch', or 'cpp'.")
}

# Turn off chatty printing during profiling (important)
DP_VERBOSE <- FALSE

############################################################
### 2) Source DP code (no main.R)
############################################################

source(here("dp_global", "R", "dp_global_biol.R"))
source(here("dp_global", "R", "sensitivity_transition_cost_bio.R"))
source(here("dp_global", "R", "realism_calibration.R"))

if (tolower(PROFILE_VARIANT) == "cpp") {
  library(Rcpp)
  source(here("dp_global", "src", "transition_cost_rcpp.R"))
  Rcpp::sourceCpp(here("dp_global", "src", "transition_cost_rcpp.cpp"))
  
  # Redefine transition_cost_tracks_bio_batch to use C++
  original_transition_cost_tracks_bio_batch <- transition_cost_tracks_bio_batch
  transition_cost_tracks_bio_batch <- function(
    track_dbh_t,
    track_dbh_tp1,
    interval_years,
    mu_const = 0,
    mu_gamma = 0,
    sigma0 = 1,
    sigma1 = 0,
    max_shrink = -Inf,
    k_shrink = 0,
    max_growth = Inf,
    max_growth_soft = Inf,
    k_growth = 0,
    use_measurement_error = FALSE,
    meas_sd1_a = 0.0062,
    meas_sd1_b = 0.0904,
    meas_sd2 = 4.64,
    meas_p_big = 0.05,
    h0 = 0,
    beta = 0,
    recruit_meanlog = 0,
    recruit_sdlog = 1,
    recruit_max_dbh = 200,
    recruit_lambda = 0,
    eps_tiebreak = 1e-6,
    hard_penalty = 1e6
  ) {
        if (is.list(track_dbh_tp1)) {
            mat_tp1 <- do.call(rbind, track_dbh_tp1)
        } else {
            mat_tp1 <- as.matrix(track_dbh_tp1)
        }
        transition_cost_tracks_bio_batch_rcpp(
            track_dbh_t = track_dbh_t,
            mat_tp1 = mat_tp1,
            interval_years = interval_years,
            mu_const = mu_const,
            mu_gamma = mu_gamma,
            sigma0 = sigma0,
            sigma1 = sigma1,
            max_shrink = max_shrink,
            k_shrink = k_shrink,
            max_growth = max_growth,
            max_growth_soft = max_growth_soft,
            k_growth = k_growth,
            use_measurement_error = use_measurement_error,
            meas_sd1_a = meas_sd1_a,
            meas_sd1_b = meas_sd1_b,
            meas_sd2 = meas_sd2,
            meas_p_big = meas_p_big,
            h0 = h0,
            beta = beta,
            recruit_meanlog = recruit_meanlog,
            recruit_sdlog = recruit_sdlog,
            recruit_max_dbh = recruit_max_dbh,
            recruit_lambda = recruit_lambda,
            eps_tiebreak = eps_tiebreak,
            hard_penalty = hard_penalty
        )
    }
}

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
        Bio_IntervalYears = as.numeric(census_interval_years),
        Bio_Mu_Growth = get_growth_mu_const(bio_pars[[species]]$growth),
        Bio_Gamma_Growth = {
            g <- bio_pars[[species]]$growth
            if (!is.null(g$gamma)) g$gamma else 0
        },
        Bio_Sigma0_Growth = bio_pars[[species]]$growth$sigma0,
        Bio_Sigma1_Growth = bio_pars[[species]]$growth$sigma1,
        Bio_H0 = bio_pars[[species]]$mortality$h0,
        Bio_Beta = bio_pars[[species]]$mortality$beta,
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
        interval_years = census_interval_years,
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

run_dp_scalar <- function(dt) {
    match_stems_dp_global_backward_marginals(
        tree_data = data.table::copy(dt),
        min_growth = min_annual_growth,
        max_growth = max_annual_growth,
        anchor_start = anchor_start_census,
        max_tracks = dp_max_tracks_local,
        max_states = dp_max_states,
        slack_tracks = dp_slack_tracks,
        temperature = DP_POSTERIOR_TEMPERATURE,
        posterior_top_k = DP_POSTERIOR_TOP_K,
        eps_tiebreak = DP_EPS_TIEBREAK,
        use_measurement_error = TRUE,
        verbose = isTRUE(DP_VERBOSE)
    )
}

run_dp_batch <- function(dt) {
    match_stems_dp_global_backward_marginals_batch(
        tree_data = data.table::copy(dt),
        min_growth = min_annual_growth,
        max_growth = max_annual_growth,
        anchor_start = anchor_start_census,
        max_tracks = dp_max_tracks_local,
        max_states = dp_max_states,
        slack_tracks = dp_slack_tracks,
        temperature = DP_POSTERIOR_TEMPERATURE,
        posterior_top_k = DP_POSTERIOR_TOP_K,
        eps_tiebreak = DP_EPS_TIEBREAK,
        use_measurement_error = TRUE,
        verbose = isTRUE(DP_VERBOSE)
    )
}

if (isTRUE(RUN_BENCHMARK)) {
    message("[profiling_code] Benchmarking scalar vs batch (Tag=", which_tag, ", species=", sp0, ")")
    if (isTRUE(BENCHMARK_WARMUP)) {
        invisible(run_dp_scalar(dtg))
        invisible(run_dp_batch(dtg))
        gc()
    }

    gc()
    t_scalar <- system.time({
        out_scalar <- run_dp_scalar(dtg)
    })
    gc()
    t_batch <- system.time({
        out_batch <- run_dp_batch(dtg)
    })

    elapsed_scalar <- as.numeric(t_scalar[["elapsed"]])
    elapsed_batch <- as.numeric(t_batch[["elapsed"]])
    speedup <- if (is.finite(elapsed_batch) && elapsed_batch > 0) elapsed_scalar / elapsed_batch else NA_real_

    message(
        sprintf(
            "[profiling_code] Benchmark results: scalar=%.3fs, batch=%.3fs, speedup=%.2fx",
            elapsed_scalar,
            elapsed_batch,
            speedup
        )
    )

    ok <- identical(out_scalar$ReconstructedStemID, out_batch$ReconstructedStemID)
    message("[profiling_code] MAP match (ReconstructedStemID identical): ", ok)
}

if (!isTRUE(RUN_PROFILE)) {
    message("[profiling_code] RUN_PROFILE=FALSE; skipping Rprof.")
    quit(save = "no", status = 0)
}

prof_file <- file.path(
    root_dir,
    "STEM_IDENTIFICATION_TEST",
    "DP_GLOBAL",
    "dev",
    paste0("dp_global_", tolower(PROFILE_VARIANT), ".prof")
)

message("[profiling_code] Profiling DP run (Tag=", which_tag, ", species=", sp0, ", variant=", PROFILE_VARIANT, ")")
message("[profiling_code] Writing profile to: ", prof_file)

Rprof(prof_file, interval = 0.01, memory.profiling = TRUE)

out_prof <- if (isTRUE(RUN_DP_MARGINALS)) {
    if (tolower(PROFILE_VARIANT) %in% c("batch", "cpp")) {
        run_dp_batch(dtg)
    } else {
        run_dp_scalar(dtg)
    }
} else {
    match_stems_dp_global_backward(
        tree_data = data.table::copy(dtg),
        min_growth = min_annual_growth,
        max_growth = max_annual_growth,
        anchor_start = anchor_start_census,
        max_tracks = dp_max_tracks_local,
        max_states = dp_max_states,
        slack_tracks = dp_slack_tracks,
        use_measurement_error = TRUE,
        verbose = isTRUE(DP_VERBOSE)
    )
}

Rprof(NULL)

message("[profiling_code] Done. Top hotspots:\n")
s <- summaryRprof(prof_file)
print(utils::head(s$by.self, 25))
cat("\n---- by.total ----\n")
print(utils::head(s$by.total, 25))

invisible(out_prof)
