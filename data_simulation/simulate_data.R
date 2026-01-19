rm(list = ls())
################################################################################
### FOREST CENSUS DATA SIMULATION FOR STEM IDENTIFICATION TESTING
################################################################################

# PURPOSE:
# This script generates realistic forest census data for testing stem identification
# and reconstruction algorithms. It simulates multi-species tropical forest dynamics
# using biologically plausible models derived from the DP workflow in dp_global_biol.R.
#
# KEY FEATURES:
# - Multi-species forests with species-specific traits
# - Multi-stem trees with recruitment and mortality
# - Realistic growth trajectories with process variability
# - Measurement error models (diameter-dependent noise + blunders)
# - Flexible growth scaling for drought/other disturbance scenarios
# - Comprehensive output: CSV data + trajectory visualization PDFs
#
# WORKFLOW:
# 1. Configure species and simulation parameters
# 2. Generate species-specific parameter scaling
# 3. Simulate individual stem trajectories with biological processes
# 4. Apply measurement error and masking for realistic observation conditions
# 5. Export data and create diagnostic plots
#
# OUTPUT FILES:
# - simulated_data_[config].csv: Main dataset with all stem measurements
# - sp_lvl_traj_[config].pdf: Species-level growth trajectory plots
# - tg_lvl_traj_[config].pdf: Tag-level individual stem trajectory plots
#
# CONFIGURATION:
# Filenames include measurement error status, species count, and growth scaling details for easy
# identification of simulation conditions (e.g., meas_error_3sp_scaling_sp1-c3-0.5_all-c5-1.2)

# Reproducible simulation
set.seed(1234)

if (!requireNamespace("data.table", quietly = TRUE)) {
    install.packages("data.table")
}

library(data.table)

if (!requireNamespace("here", quietly = TRUE)) {
    stop("Please install the 'here' package to run this script.")
}
library(here)

################################################################################
### SIMULATION PARAMETERS
################################################################################

# This section contains all configurable parameters for the simulation.
# Parameters are organized by biological process for clarity.
# Modify these to customize your simulation scenario.

params <- list(
    # =========================================================================
    # SIMULATION STRUCTURE
    # =========================================================================
    sim = list(
        n_census = 9L, # Number of census intervals
        census_interval_years = 5 # Years between censuses
    ),

    # =========================================================================
    # SPECIES CONFIGURATION
    # =========================================================================
    # Controls the number of species, their abundance, and trait variation
    species = list(
        n_species = 3L, # Number of species to simulate
        n_trees_per_species = c(10L, 15L, 13L), # Trees per species (length must = n_species)
        max_stems = 6L, # Maximum stems per tree (sampled uniformly from 2:max_stems)
        scale_range = c(1, 1.7), # Trait scaling range (evenly distributed)
        species_names = NULL # Custom names (NULL = auto-generate sp1, sp2, ...)
    ),

    # =========================================================================
    # RECRUITMENT PROCESS
    # =========================================================================
    # Controls how new stems enter the population
    recruitment = list(
        recruit_prob = 0.3, # Probability of recruitment per available stem slot per census
        threshold_dbh = 1, # DBH threshold for stems to become observable (cm)
        meanlog = log(2), # Lognormal mean for recruit DBH (matches DP model)
        sdlog = 0.8 # Lognormal SD for recruit DBH
    ),

    # =========================================================================
    # GROWTH PROCESS
    # =========================================================================
    # Size-dependent growth model: mu(DBH) = alpha + gamma*log(DBH)
    # Variance model: sigma(DBH) = sigma0 + sigma1*DBH
    growth = list(
        alpha = 0.4, # Intercept for mean annual growth (cm/year)
        gamma = 0.2, # Slope for log(DBH) in mean growth
        sigma0 = 0.1, # Intercept for growth SD (cm/year)
        sigma1 = 0.01, # Slope for DBH in growth SD
        min_annual_growth = 0, # Minimum allowed annual growth (cm/year)
        max_annual_growth = 7.5 # Maximum allowed annual growth (cm/year)
    ),

    # =========================================================================
    # INITIAL SIZE DISTRIBUTION
    # =========================================================================
    initial_dbh = list(
        census1_meanlog = log(12), # Lognormal mean for census-1 stems
        census1_sdlog = 0.9, # Lognormal SD for census-1 stems
        recruit_meanlog = log(0.7), # Lognormal mean for recruits
        recruit_sdlog = 0.4, # Lognormal SD for recruits
        min_dbh_true = 0.1, # Minimum true DBH (cm)
        max_dbh_true = 210 # Maximum true DBH (cm)
    ),

    # =========================================================================
    # MEASUREMENT ERROR MODEL
    # =========================================================================
    # Implements the DP workflow measurement error model:
    # - Small errors: SD = a*DBH + b (diameter-dependent)
    # - Large errors (blunders): SD = constant, probability p_big
    obs = list(
        use_measurement_error = TRUE, # Enable/disable measurement error
        meas_sd1_a = 0.0062, # Slope for diameter-dependent error
        meas_sd1_b = 0.0904, # Intercept for diameter-dependent error
        meas_sd2 = 4.64, # SD for large measurement errors (blunders)
        meas_p_big = 0.05, # Probability of large measurement error
        min_dbh_obs = 1.0 # Minimum observed DBH (cm) - matches measurement threshold
    ),

    # =========================================================================
    # GROWTH SCALING EVENTS
    # =========================================================================
    # Apply growth multipliers to simulate disturbances (drought, etc.)
    # Each event specifies: species ("all" or species name), census, multiplier
    # Multiple events per species-census combination are multiplied together
    # Example: sp1 enhanced growth (1.5) in census 2 and 5, sp3 reduced growth (0.5) in census 2 and 5
    growth_scaling = list(
        events = list(
            list(species = "sp1", census = 2, multiplier = 1.7),
            list(species = "sp1", census = 5, multiplier = 1.7),
            list(species = "sp1", census = 8, multiplier = 1.7),
            list(species = "sp3", census = 2, multiplier = 0.1),
            list(species = "sp3", census = 5, multiplier = 0.1),
            list(species = "sp3", census = 8, multiplier = 0.1)
        )
    ),

    # =========================================================================
    # MORTALITY PROCESS
    # =========================================================================
    # Hazard model: hazard(DBH) = h0 * exp(beta * DBH)
    # P(death) = 1 - exp(-hazard * interval_years)
    mortality = list(
        h0 = 0.004, # Baseline hazard rate
        beta = 0.02 # DBH effect on hazard
    ),

    # =========================================================================
    # STEM ID MASKING
    # =========================================================================
    # Mimics ForestGEO protocol where IDs are only trusted in recent censuses
    mask = list(
        anchor_start_census = 7L # Censuses before this have TrueStemID masked
    ),

    # =========================================================================
    # VISUALIZATION
    # =========================================================================
    plot = list(
        make_plot = TRUE # Generate trajectory plots
    )
)

################################################################################
### SPECIES CONFIGURATION FUNCTIONS
################################################################################

#' Generate Species Configuration Table
#'
#' Creates a data.table defining species characteristics for the simulation.
#' Each species gets a trait scaling factor evenly distributed across the specified range.
#'
#' @param params List containing species configuration parameters
#' @return data.table with columns: Species (name), Scale (trait multiplier), n_trees (count)
#' @examples
#' params <- list(species = list(
#'     n_species = 2, n_trees_per_species = c(10, 15),
#'     scale_range = c(1.0, 1.5), species_names = c("Oak", "Pine")
#' ))
#' species_table <- generate_species_table(params)
generate_species_table <- function(params) {
    n_species <- params$species$n_species
    n_trees_per_species <- params$species$n_trees_per_species
    scale_range <- params$species$scale_range
    species_names <- params$species$species_names

    # Validate inputs
    if (n_species < 1) stop("n_species must be >= 1")
    if (length(scale_range) != 2) stop("scale_range must be a vector of length 2")
    if (scale_range[1] > scale_range[2]) stop("scale_range[1] must be <= scale_range[2]")

    if (length(n_trees_per_species) == 1) {
        n_trees_per_species <- rep(n_trees_per_species, n_species)
    } else if (length(n_trees_per_species) != n_species) {
        stop("Length of n_trees_per_species must equal n_species or be 1")
    }

    # Generate species names if not provided
    if (is.null(species_names)) {
        species_names <- paste0("sp", seq_len(n_species))
    } else if (length(species_names) != n_species) {
        stop("Length of species_names must equal n_species")
    }

    # Generate scales evenly distributed between scale_range[1] and scale_range[2]
    if (n_species == 1) {
        scales <- scale_range[1]
    } else {
        scales <- seq(scale_range[1], scale_range[2], length.out = n_species)
    }

    data.table(
        Species = species_names,
        Scale = scales,
        n_trees = n_trees_per_species
    )
}

species_table <- generate_species_table(params)

# Display the generated species configuration
cat("Generated species configuration:\n")
print(species_table)
cat("\n")

#' Scale Species Parameters
#'
#' Applies species-specific trait scaling to base parameters. This creates
#' interspecific variation by multiplying selected parameters by the species scale factor.
#'
#' Scaling affects:
#' - Growth parameters (alpha, gamma, sigma0, sigma1)
#' - Recruitment probability (via odds scaling to preserve [0,1] bounds)
#' - Mortality baseline rate (h0)
#'
#' @param base_params List of base parameter values
#' @param scale Numeric scaling factor for the species
#' @return Modified parameter list with species-specific scaling applied
scale_species_params <- function(base_params, scale) {
    p <- base_params

    # Helper function to scale probabilities while preserving [0,1] bounds
    # Converts to log-odds, scales, converts back
    scale_prob_odds <- function(prob, odds_multiplier) {
        if (is.na(prob) || prob <= 0) {
            return(0)
        }
        if (prob >= 1) {
            return(1)
        }
        plogis(qlogis(prob) + log(odds_multiplier))
    }

    # Apply species-specific scaling to growth parameters
    p$growth$alpha <- base_params$growth$alpha * scale
    p$growth$gamma <- base_params$growth$gamma * scale
    p$growth$sigma0 <- base_params$growth$sigma0 * scale
    p$growth$sigma1 <- base_params$growth$sigma1 * scale

    # Scale recruitment probability (preserves bounds)
    p$recruitment$recruit_prob <- scale_prob_odds(base_params$recruitment$recruit_prob, scale)

    # Scale mortality baseline rate
    p$mortality$h0 <- base_params$mortality$h0 * scale

    p
}

################################################################################
### HELPER FUNCTIONS
################################################################################

#' Truncated Lognormal Random Sampling
#'
#' Generates random samples from a lognormal distribution, truncated to specified bounds.
#' Uses rejection sampling to ensure all values fall within [min, max].
#'
#' @param n Number of samples to generate
#' @param meanlog Mean of the lognormal distribution on log scale
#' @param sdlog SD of the lognormal distribution on log scale
#' @param min Minimum allowed value
#' @param max Maximum allowed value
#' @return Numeric vector of length n with truncated lognormal samples
rtrunc_lnorm <- function(n, meanlog, sdlog, min = 0, max = Inf) {
    out <- numeric(0)
    while (length(out) < n) {
        draw <- rlnorm(n = max(100, n), meanlog = meanlog, sdlog = sdlog)
        draw <- draw[draw >= min & draw <= max]
        out <- c(out, draw)
    }
    out[seq_len(n)]
}

#' Create Growth Scaling Matrix
#'
#' Generates a matrix of growth multipliers for each species-census combination.
#' Default multiplier is 1.0 (no scaling). Growth scaling events are multiplied together,
#' allowing multiple compounding effects per species-census combination.
#'
#' @param n_census Number of census intervals
#' @param species_names Character vector of species names
#' @param growth_scaling_events List of scaling events (each with species, census, multiplier)
#' @return Matrix with species as rows, censuses as columns, values as compounded multipliers
make_growth_multipliers <- function(n_census, species_names, growth_scaling_events) {
    # Initialize matrix with default multiplier of 1.0
    multipliers <- matrix(1.0,
        nrow = length(species_names), ncol = n_census,
        dimnames = list(species_names, paste0("census_", 1:n_census))
    )

    # Apply scaling events (multiply to allow compounding effects)
    for (event in growth_scaling_events) {
        species <- event$species
        census <- event$census
        multiplier <- event$multiplier

        if (species == "all") {
            # Apply to all species for this census (multiply existing values)
            multipliers[, paste0("census_", census)] <- multipliers[, paste0("census_", census)] * multiplier
        } else {
            # Apply to specific species (multiply existing value)
            if (species %in% species_names) {
                multipliers[species, paste0("census_", census)] <- multipliers[species, paste0("census_", census)] * multiplier
            } else {
                warning(sprintf("Species '%s' not found in species list", species))
            }
        }
    }

    multipliers
}

################################################################################
### SIMULATION ENGINE
################################################################################

# Create growth scaling matrix for all species-census combinations
growth_multipliers <- make_growth_multipliers(
    n_census = params$sim$n_census,
    species_names = species_table$Species,
    growth_scaling_events = params$growth_scaling$events
)

#' Simulate Single Stem Trajectory
#'
#' Generates a complete DBH trajectory for one stem, including:
#' - Birth timing and initial size
#' - Growth with process variability and species-specific scaling
#' - Mortality based on size-dependent hazard
#' - Measurement error on observed DBH
#'
#' @param tag Tree identifier
#' @param original_stem_id Stem identifier within tree
#' @param species Species name
#' @param growth_multipliers Matrix of growth scaling factors (species × census)
#' @param params Parameter list
#' @return data.table with one row per census for this stem
simulate_one_stem <- function(tag, original_stem_id, species, growth_multipliers, params) {
    n_census <- params$sim$n_census
    interval_years <- params$sim$census_interval_years
    threshold <- params$recruitment$threshold_dbh

    # Determine birth census (recruitment timing)
    birth_census <- if (runif(1) < params$recruitment$recruit_prob) {
        sample(2:(n_census - 1L), 1L) # Recruit after census 1
    } else {
        1L # Established stem present from start
    }

    # Initialize true DBH trajectory
    dbh_true <- rep(NA_real_, n_census)
    annual_growth <- rep(NA_real_, n_census)
    death_census <- n_census

    # Set initial size based on birth timing
    if (birth_census == 1L) {
        # Established stem - larger initial size
        dbh_true[birth_census] <- rtrunc_lnorm(
            1,
            meanlog = params$initial_dbh$census1_meanlog,
            sdlog = params$initial_dbh$census1_sdlog,
            min = params$initial_dbh$min_dbh_true,
            max = params$initial_dbh$max_dbh_true
        )
    } else {
        # Recruit - smaller initial size
        dbh_true[birth_census] <- rtrunc_lnorm(
            1,
            meanlog = params$recruitment$meanlog,
            sdlog = params$recruitment$sdlog,
            min = params$initial_dbh$min_dbh_true,
            max = threshold * 0.99
        )
    }

    # Simulate forward through time until death or end
    for (t in birth_census:(n_census - 1L)) {
        if (is.na(dbh_true[t])) break

        # Calculate expected growth with species-specific scaling
        mu <- params$growth$alpha + params$growth$gamma * log(dbh_true[t])
        census_multiplier <- growth_multipliers[species, paste0("census_", t)]
        mu <- mu * census_multiplier

        # Add process variability
        sigma <- params$growth$sigma0 + params$growth$sigma1 * dbh_true[t]
        g_ann <- mu + rnorm(1, 0, sigma)
        g_ann <- pmax(params$growth$min_annual_growth, pmin(g_ann, params$growth$max_annual_growth))
        annual_growth[t] <- g_ann

        # Update DBH
        dbh_true[t + 1L] <- pmax(params$initial_dbh$min_dbh_true, dbh_true[t] + (g_ann * interval_years))

        # Check for mortality
        hazard <- params$mortality$h0 * exp(params$mortality$beta * dbh_true[t])
        p_die <- 1 - exp(-hazard * interval_years)
        if (runif(1) < p_die) {
            death_census <- t
            dbh_true[(t + 1L):n_census] <- NA_real_
            annual_growth[(t + 1L):n_census] <- NA_real_
            break
        }
    }

    # Generate observations (only after crossing size threshold)
    dbh_obs <- rep(NA_real_, n_census)
    obs_sd <- rep(NA_real_, n_census)
    obs_idx <- which(!is.na(dbh_true) & dbh_true >= threshold)
    if (length(obs_idx) > 0L) {
        obs_idx <- obs_idx[obs_idx >= birth_census & obs_idx <= death_census]
    }
    if (length(obs_idx) > 0L) {
        if (isTRUE(params$obs$use_measurement_error)) {
            # Apply measurement error model from DP workflow
            # Small errors: diameter-dependent SD
            sd1 <- pmax(params$obs$meas_sd1_a * dbh_true[obs_idx] + params$obs$meas_sd1_b, 1e-6)
            # Large errors: constant SD (blunders)
            sd2 <- params$obs$meas_sd2
            # Mixture model: small error with prob (1-p_big), large error with prob p_big
            use_large_error <- runif(length(obs_idx)) < params$obs$meas_p_big
            obs_sd[obs_idx] <- ifelse(use_large_error, sd2, sd1)
        } else {
            # Perfect observations (no measurement error)
            obs_sd[obs_idx] <- 0
        }

        # Add measurement error to true DBH
        dbh_obs[obs_idx] <- dbh_true[obs_idx] + rnorm(length(obs_idx), 0, obs_sd[obs_idx])
        dbh_obs[obs_idx] <- pmax(params$obs$min_dbh_obs, dbh_obs[obs_idx])
    }

    # Return complete stem trajectory
    data.table(
        Species = species,
        Tag = tag,
        OriginalStemID = original_stem_id,
        TrueStemID = original_stem_id, # Will be masked later for early censuses
        CensusID = seq_len(n_census),
        BirthCensus = birth_census,
        DeathCensus = death_census,
        DBH_true = dbh_true, # True DBH (for diagnostics)
        DBH = dbh_obs, # Observed DBH (with measurement error)
        AnnualGrowth = annual_growth, # Annual growth rate (cm/year)
        YearFactor = growth_multipliers[species, ], # Applied scaling factors
        ObsSD = obs_sd # Measurement error SD (for diagnostics)
    )
}

#' Simulate Single Tree
#'
#' Generates all stems for one tree by calling simulate_one_stem multiple times.
#' Each tree has a random number of stems (2 to max_stems).
#'
#' @param tag Tree identifier
#' @param species Species name
#' @param growth_multipliers Matrix of growth scaling factors
#' @param params Parameter list
#' @return data.table with trajectories for all stems of this tree
simulate_one_tree <- function(tag, species, growth_multipliers, params) {
    n_stems <- sample(2:params$species$max_stems, 1L)
    rbindlist(lapply(seq_len(n_stems), function(stem_id) {
        simulate_one_stem(
            tag = tag,
            original_stem_id = stem_id,
            species = species,
            growth_multipliers = growth_multipliers,
            params = params
        )
    }))
}

#' Simulate Single Species
#'
#' Generates all trees for one species, applying species-specific parameter scaling.
#'
#' @param species Species name
#' @param scale Species-specific trait scaling factor
#' @param n_trees Number of trees to generate for this species
#' @param tag_offset Starting tag number for this species
#' @param growth_multipliers Matrix of growth scaling factors
#' @param base_params Base parameter list (before species scaling)
#' @return data.table with all stem trajectories for this species
simulate_one_species <- function(species, scale, n_trees, tag_offset, growth_multipliers, base_params) {
    # Apply species-specific parameter scaling
    p_species <- scale_species_params(base_params = base_params, scale = scale)

    # Generate all trees for this species
    rbindlist(lapply(seq_len(n_trees), function(i) {
        tag <- as.integer(tag_offset + i)
        simulate_one_tree(tag = tag, species = species, growth_multipliers = growth_multipliers, params = p_species)
    }))
}

# Generate the complete forest dataset
tag_offset <- 0L
dt_list <- vector("list", nrow(species_table))
for (i in seq_len(nrow(species_table))) {
    sp <- species_table[i]
    dt_list[[i]] <- simulate_one_species(
        species = sp$Species,
        scale = sp$Scale,
        n_trees = as.integer(sp$n_trees),
        tag_offset = tag_offset,
        growth_multipliers = growth_multipliers,
        base_params = params
    )
    tag_offset <- tag_offset + as.integer(sp$n_trees)
}
dt <- rbindlist(dt_list)

################################################################################
### DATA PROCESSING AND OUTPUT
################################################################################

#' Create Scaling Indicator for Filenames
#'
#' Generates a descriptive string encoding measurement error status, number of species,
#' and growth scaling configuration for use in output filenames. This makes datasets
#' self-documenting and easily identifiable by simulation parameters.
#'
#' @param events List of growth scaling events
#' @param use_measurement_error Logical indicating if measurement error is enabled
#' @param n_species Number of species simulated
#' @return Character string for filename encoding (e.g., "merr_3sp_inc1_c2c5c8_p1p7_dec1_c2c5c8_p0p1")
create_scaling_indicator <- function(events, use_measurement_error = TRUE, n_species = 2L) {
    # Start with measurement error indicator
    meas_indicator <- if (isTRUE(use_measurement_error)) "merr" else "no_merr"

    # Add species count
    species_indicator <- sprintf("%dsp", n_species)

    # Add scaling information
    if (length(events) == 0) {
        scaling_part <- "no_scaling"
    } else {
        # Group events by increase/decrease
        inc_events <- events[sapply(events, function(e) e$multiplier > 1)]
        dec_events <- events[sapply(events, function(e) e$multiplier < 1)]

        parts <- c()

        if (length(inc_events) > 0) {
            inc_species <- unique(sapply(inc_events, `[[`, "species"))
            n_inc <- length(inc_species)
            inc_censuses <- sort(unique(sapply(inc_events, `[[`, "census")))
            inc_mult <- unique(sapply(inc_events, `[[`, "multiplier"))
            if (length(inc_mult) != 1) stop("Multiple multipliers for increase events")
            inc_mult_str <- gsub("\\.", "p", as.character(inc_mult))
            inc_cens_str <- paste0("c", inc_censuses, collapse = "")
            parts <- c(parts, paste0("inc", n_inc, "_", inc_cens_str, "_p", inc_mult_str))
        }

        if (length(dec_events) > 0) {
            dec_species <- unique(sapply(dec_events, `[[`, "species"))
            n_dec <- length(dec_species)
            dec_censuses <- sort(unique(sapply(dec_events, `[[`, "census")))
            dec_mult <- unique(sapply(dec_events, `[[`, "multiplier"))
            if (length(dec_mult) != 1) stop("Multiple multipliers for decrease events")
            dec_mult_str <- gsub("\\.", "p", as.character(dec_mult))
            dec_cens_str <- paste0("c", dec_censuses, collapse = "")
            parts <- c(parts, paste0("dec", n_dec, "_", dec_cens_str, "_p", dec_mult_str))
        }

        if (length(parts) == 0) {
            scaling_part <- "no_scaling"
        } else {
            scaling_part <- paste(parts, collapse = "_")
        }
    }

    # Combine all indicators
    paste(meas_indicator, species_indicator, scaling_part, sep = "_")
}

#' Build File Name for Outputs
#'
#' Generates standardized filenames for simulation outputs based on type and scaling indicator.
#'
#' @param type Character string indicating the type of file: "data", "sp_lvl_traj", "tg_lvl_traj"
#' @param scaling_indicator Character string with the scaling configuration
#' @return Character string with the filename
build_file_name <- function(type, scaling_indicator) {
    switch(type,
        "data" = sprintf("simulated_data_%s.csv", scaling_indicator),
        "sp_lvl_traj" = sprintf("sp_lvl_traj_%s.pdf", scaling_indicator),
        "tg_lvl_traj" = sprintf("tg_lvl_traj_%s.pdf", scaling_indicator),
        stop("Unknown type: ", type)
    )
}

################################################################################
### ADD A FEW EXTRA TEST CASES
################################################################################

# Simulate and add one extra tag with exactly one stem for each species
new_tag_id <- max(dt$Tag) + 1L
for (sp in c(species_table$Species)) {
    scale <- species_table[Species == sp, Scale]
    p_species <- scale_species_params(base_params = params, scale = scale)
    new_tree <- simulate_one_stem(
        tag = new_tag_id,
        original_stem_id = 1L,
        species = sp,
        growth_multipliers = growth_multipliers,
        params = p_species
    )
    dt <- rbind(dt, new_tree, fill = TRUE)
    new_tag_id <- new_tag_id + 1L
}

# Simulate and add one extra tag with exactly one stem for each species that dies before the anchor census
anchor_census <- params$mask$anchor_start_census
new_tag_id <- max(dt$Tag) + 1L
for (sp in species_table$Species) {
    scale <- species_table[Species == sp, Scale]
    p_species <- scale_species_params(base_params = params, scale = scale)
    # Simulate one stem with death before anchor census
    n_census <- params$sim$n_census
    birth_census <- 1L
    death_census <- sample(2:(anchor_census - 1L), 1L)
    # Generate DBH trajectory
    dbh_true <- rep(NA_real_, n_census)
    annual_growth <- rep(NA_real_, n_census)
    dbh_true[birth_census] <- rtrunc_lnorm(
        1,
        meanlog = params$initial_dbh$census1_meanlog,
        sdlog = params$initial_dbh$census1_sdlog,
        min = params$initial_dbh$min_dbh_true,
        max = params$initial_dbh$max_dbh_true
    )
    for (t in birth_census:(death_census - 1L)) {
        mu <- params$growth$alpha + params$growth$gamma * log(dbh_true[t])
        census_multiplier <- growth_multipliers[sp, paste0("census_", t)]
        mu <- mu * census_multiplier
        sigma <- params$growth$sigma0 + params$growth$sigma1 * dbh_true[t]
        g_ann <- mu + rnorm(1, 0, sigma)
        g_ann <- pmax(params$growth$min_annual_growth, pmin(g_ann, params$growth$max_annual_growth))
        annual_growth[t] <- g_ann
        dbh_true[t + 1L] <- pmax(params$initial_dbh$min_dbh_true, dbh_true[t] + (g_ann * params$sim$census_interval_years))
    }
    dbh_true[(death_census + 1L):n_census] <- NA_real_
    annual_growth[(death_census + 1L):n_census] <- NA_real_
    # Generate observed DBH
    dbh_obs <- rep(NA_real_, n_census)
    obs_sd <- rep(NA_real_, n_census)
    threshold <- params$recruitment$threshold_dbh
    obs_idx <- which(!is.na(dbh_true) & dbh_true >= threshold)
    if (length(obs_idx) > 0L) {
        obs_idx <- obs_idx[obs_idx >= birth_census & obs_idx <= death_census]
        if (isTRUE(params$obs$use_measurement_error)) {
            sd1 <- pmax(params$obs$meas_sd1_a * dbh_true[obs_idx] + params$obs$meas_sd1_b, 1e-6)
            sd2 <- params$obs$meas_sd2
            use_large_error <- runif(length(obs_idx)) < params$obs$meas_p_big
            obs_sd[obs_idx] <- ifelse(use_large_error, sd2, sd1)
        } else {
            obs_sd[obs_idx] <- 0
        }
        dbh_obs[obs_idx] <- dbh_true[obs_idx] + rnorm(length(obs_idx), 0, obs_sd[obs_idx])
        dbh_obs[obs_idx] <- pmax(params$obs$min_dbh_obs, dbh_obs[obs_idx])
    }
    new_tree <- data.table(
        Species = sp,
        Tag = new_tag_id,
        OriginalStemID = 1L,
        TrueStemID = 1L,
        CensusID = seq_len(n_census),
        BirthCensus = birth_census,
        DeathCensus = death_census,
        DBH_true = dbh_true,
        DBH = dbh_obs,
        AnnualGrowth = annual_growth,
        YearFactor = growth_multipliers[sp, ],
        ObsSD = obs_sd
    )
    dt <- rbind(dt, new_tree, fill = TRUE)
    new_tag_id <- new_tag_id + 1L
}

dt[CensusID < params$mask$anchor_start_census, TrueStemID := NA_integer_]
dt[is.na(DBH), TrueStemID := NA_integer_]

################################################################################
### EXPORT
################################################################################

# Display summary statistics
cat("Simulation complete. Data summary:\n")
cat(sprintf("Total stems: %d\n", length(unique(dt$Tag))))
cat(sprintf("Total observations: %d\n", nrow(dt[!is.na(DBH)])))
cat("Observations per species:\n")
dt[, .N, by = Species]

# Generate filename indicator and export main dataset
scaling_indicator <- create_scaling_indicator(params$growth_scaling$events, params$obs$use_measurement_error, params$species$n_species)
data_filename <- build_file_name("data", scaling_indicator)

# Export simulation parameters to text file
params_filename <- here("data_simulation", "data", sprintf("simulation_params_%s.txt", scaling_indicator))
capture.output(str(params), file = params_filename)

fwrite(dt, here("data_simulation", "data", data_filename))
# Apply stem ID masking to simulate ForestGEO protocol
# In early censuses, stem identities are not trusted (TrueStemID = NA)
################################################################################
### DIAGNOSTIC PLOTS
################################################################################

# Generate trajectory plots for visual inspection of simulation results
# Creates two types of PDFs:
# 1. Species-level: All stems grouped by species
# 2. Tag-level: Individual stem trajectories for each tree

if (isTRUE(params$plot$make_plot) && requireNamespace("ggplot2", quietly = TRUE)) {
    # =========================================================================
    # PDF 1: Species-level trajectories (one page per species)
    # =========================================================================
    # Shows all stem trajectories for each species on separate pages
    # Useful for checking species-specific growth patterns and scaling effects
    pdf_file_species <- here("data_simulation", "data", build_file_name("sp_lvl_traj", scaling_indicator))
    pdf(pdf_file_species, width = 8, height = 6)

    # Get unique species
    species_list <- unique(dt$Species)

    # Create one plot per species (all tags)
    for (species_name in species_list) {
        species_data <- dt[Species == species_name]

        gg <- ggplot2::ggplot(
            species_data,
            # ggplot2::aes(x = CensusID, y = DBH, group = interaction(OriginalStemID, Tag), color = factor(Tag))
            ggplot2::aes(x = CensusID, y = DBH, group = interaction(OriginalStemID, Tag))
        ) +
            ggplot2::geom_line(na.rm = TRUE, alpha = 0.7) +
            ggplot2::geom_point(size = 1.5, na.rm = TRUE, alpha = 0.7) +
            ggplot2::theme_minimal() +
            ggplot2::labs(
                title = paste("Growth Trajectories -", species_name),
                subtitle = sprintf("All tags for this species\n%s", scaling_indicator),
                x = "Census ID",
                y = "DBH (cm)"
            ) +
            ggplot2::theme(
                plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
                plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10),
                legend.position = "bottom"
            )

        print(gg)
    }

    dev.off()
    cat(sprintf("Species-level trajectories PDF saved as: %s\n", pdf_file_species))

    # =========================================================================
    # PDF 2: Tag-level trajectories (one page per tree)
    # =========================================================================
    # Shows individual stem trajectories for each tree
    # Useful for checking multi-stem dynamics and stem identification challenges
    pdf_file_tags <- here("data_simulation", "data", build_file_name("tg_lvl_traj", scaling_indicator))
    pdf(pdf_file_tags, width = 8, height = 6)

    # Create one plot per tag
    for (tag_id in unique(dt$Tag)) {
        tag_data <- dt[Tag == tag_id]
        species_name <- unique(tag_data$Species)

        gg <- ggplot2::ggplot(
            tag_data,
            ggplot2::aes(x = CensusID, y = DBH, group = interaction(OriginalStemID), color = factor(OriginalStemID))
        ) +
            ggplot2::geom_line(na.rm = TRUE, alpha = 0.8) +
            ggplot2::geom_point(size = 2, na.rm = TRUE, alpha = 0.8) +
            ggplot2::theme_minimal() +
            ggplot2::labs(
                title = sprintf("Tag %d - %s", tag_id, species_name),
                subtitle = sprintf("Individual stem trajectories\n%s", scaling_indicator),
                x = "Census ID",
                y = "DBH (cm)",
                color = "Stem ID"
            ) +
            ggplot2::theme(
                plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
                plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10),
                legend.position = "bottom"
            )

        print(gg)
    }

    dev.off()
    cat(sprintf("Tag-level trajectories PDF saved as: %s\n", pdf_file_tags))
}

# =============================================================================
# SIMULATION COMPLETE
# =============================================================================
# The script has generated:
# 1. simulated_data_[config].csv - Main dataset for stem ID testing
# 2. sp_lvl_traj_[config].pdf - Species-level trajectory plots
# 3. tg_lvl_traj_[config].pdf - Tag-level trajectory plots
#
# Filename [config] encodes: measurement error status + species count + number of species increasing/decreasing in specific censuses with multipliers

# Example: "merr_3sp_inc1_c2c5c8_p1p7_dec1_c2c5c8_p0p1"
#
# Use this data to test and validate stem identification algorithms!
