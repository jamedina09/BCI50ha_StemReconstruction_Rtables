rm(list = ls())

################################################################################
### FOREST CENSUS DATA SIMULATION FOR STEM IDENTIFICATION TESTING
################################################################################

# This script simulates forest census data for testing stem identification algorithms.
# It generates synthetic tree growth trajectories with measurement error, recruitment,
# growth, and mortality processes.

################################################################################
### SETUP
################################################################################

set.seed(1234)
library(data.table)
library(here)

################################################################################
### SIMULATION PARAMETERS
################################################################################
params <- list(
    sim = list(
        n_census = 9L,
        census_interval_years = 5
    ),
    species = list(
        n_species = 3L,
        n_trees_per_species = c(10L, 15L, 13L),
        max_stems = 6L,
        scale_range = c(1, 1.7),
        species_names = NULL
    ),
    recruitment = list(
        recruit_prob = 0.3,
        threshold_dbh = 1,
        meanlog = log(2),
        sdlog = 0.8
    ),
    growth = list(
        alpha = 0.4,
        gamma = 0.2,
        sigma0 = 0.1,
        sigma1 = 0.01,
        min_annual_growth = 0,
        max_annual_growth = 7.5
    ),
    initial_dbh = list(
        census1_meanlog = log(12),
        census1_sdlog = 0.9,
        recruit_meanlog = log(0.7),
        recruit_sdlog = 0.4,
        min_dbh_true = 0.1,
        max_dbh_true = 210
    ),
    obs = list(
        use_measurement_error = TRUE,
        meas_sd1_a = 0.0062,
        meas_sd1_b = 0.0904,
        meas_sd2 = 4.64,
        meas_p_big = 0.05,
        min_dbh_obs = 1.0
    ),
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
    mortality = list(
        h0 = 0.004,
        beta = 0.02
    ),
    mask = list(
        anchor_start_census = 7L
    ),
    plot = list(
        make_plot = TRUE
    )
)

################################################################################
### HELPER FUNCTIONS
################################################################################
generate_species_table <- function(params) {
    n_species <- params$species$n_species
    n_trees <- params$species$n_trees_per_species
    scale_range <- params$species$scale_range
    species_names <- params$species$species_names
    if (is.null(species_names)) species_names <- paste0("sp", seq_len(n_species))
    scales <- if (n_species == 1) scale_range[1] else seq(scale_range[1], scale_range[2], length.out = n_species)
    data.table(Species = species_names, Scale = scales, n_trees = n_trees)
}

make_growth_multipliers <- function(n_census, species_names, events) {
    multipliers <- matrix(1.0,
        nrow = length(species_names), ncol = n_census,
        dimnames = list(species_names, paste0("census_", 1:n_census))
    )
    for (e in events) {
        species <- e$species
        census <- e$census
        mult <- e$multiplier
        if (species == "all") {
            multipliers[, paste0("census_", census)] <- multipliers[, paste0("census_", census)] * mult
        } else if (species %in% species_names) {
            multipliers[species, paste0("census_", census)] <- multipliers[species, paste0("census_", census)] * mult
        }
    }
    multipliers
}

rtrunc_lnorm <- function(n, meanlog, sdlog, min = 0, max = Inf) {
    out <- numeric(0)
    while (length(out) < n) {
        draw <- rlnorm(max(100, n), meanlog, sdlog)
        out <- c(out, draw[draw >= min & draw <= max])
    }
    out[seq_len(n)]
}

scale_species_params <- function(base_params, scale) {
    p <- base_params
    p$growth$alpha <- p$growth$alpha * scale
    p$growth$gamma <- p$growth$gamma * scale
    p$growth$sigma0 <- p$growth$sigma0 * scale
    p$growth$sigma1 <- p$growth$sigma1 * scale
    p$recruitment$recruit_prob <- plogis(qlogis(p$recruitment$recruit_prob) + log(scale))
    p$mortality$h0 <- p$mortality$h0 * scale
    p
}

################################################################################
### SIMULATION ENGINE
################################################################################

# This section contains the core simulation logic:
# 1. Species configuration
# 2. Growth multipliers for scaling effects
# 3. Individual stem trajectory simulation
# 4. Tree-level simulation (multiple stems per tree)
# 5. Species-level simulation (multiple trees per species)
# 6. Full dataset assembly

# ============================================================================
# Species Configuration
# ============================================================================

species_table <- generate_species_table(params)
growth_multipliers <- make_growth_multipliers(params$sim$n_census, species_table$Species, params$growth_scaling$events)

# ============================================================================
# Individual Stem Trajectory Simulation
# ============================================================================

simulate_one_stem <- function(tag, original_stem_id, species, growth_multipliers, params, interval_years) {
    n_census <- params$sim$n_census
    threshold <- params$recruitment$threshold_dbh
    stopifnot(length(interval_years) == n_census - 1)
    census_dates <- c(0, cumsum(interval_years))

    birth_census <- if (runif(1) < params$recruitment$recruit_prob) sample(2:(n_census - 1L), 1L) else 1L

    dbh_true <- rep(NA_real_, n_census)
    annual_growth <- rep(NA_real_, n_census)
    death_census <- n_census

    if (birth_census == 1L) {
        dbh_true[1] <- rtrunc_lnorm(
            1, params$initial_dbh$census1_meanlog, params$initial_dbh$census1_sdlog,
            params$initial_dbh$min_dbh_true, params$initial_dbh$max_dbh_true
        )
    } else {
        dbh_true[birth_census] <- rtrunc_lnorm(
            1, params$recruitment$meanlog, params$recruitment$sdlog,
            params$initial_dbh$min_dbh_true, threshold * 0.99
        )
    }

    for (t in birth_census:(n_census - 1L)) {
        if (is.na(dbh_true[t])) break
        mu <- params$growth$alpha + params$growth$gamma * log(dbh_true[t])
        mu <- mu * growth_multipliers[species, paste0("census_", t)]
        sigma <- params$growth$sigma0 + params$growth$sigma1 * dbh_true[t]
        g_ann <- mu + rnorm(1, 0, sigma)
        g_ann <- pmax(params$growth$min_annual_growth, pmin(g_ann, params$growth$max_annual_growth))
        annual_growth[t] <- g_ann

        dbh_true[t + 1L] <- pmax(params$initial_dbh$min_dbh_true, dbh_true[t] + g_ann * interval_years[t])

        hazard <- params$mortality$h0 * exp(params$mortality$beta * dbh_true[t])
        p_die <- 1 - exp(-hazard * interval_years[t])
        if (runif(1) < p_die) {
            death_census <- t
            dbh_true[(t + 1L):n_census] <- NA_real_
            annual_growth[(t + 1L):n_census] <- NA_real_
            break
        }
    }

    data.table(
        Species = species, Tag = tag, OriginalStemID = original_stem_id, TrueStemID = original_stem_id,
        CensusID = seq_len(n_census), BirthCensus = birth_census, DeathCensus = death_census,
        DBH_true = dbh_true, AnnualGrowth = annual_growth,
        CensusInterval = c(NA_real_, interval_years), CensusDate = census_dates
    )
}

# ============================================================================
# Tree-Level Simulation (Multiple Stems per Tree)
# ============================================================================

simulate_one_tree <- function(tag, species, growth_multipliers, params) {
    n_census <- params$sim$n_census
    interval_years <- params$sim$census_interval_years + rnorm(n_census - 1, 0, 0.1)
    interval_years <- pmax(interval_years, 0.1)
    n_stems <- sample(2:params$species$max_stems, 1L)
    rbindlist(lapply(seq_len(n_stems), function(stem_id) {
        simulate_one_stem(tag, stem_id, species, growth_multipliers, params, interval_years)
    }))
}

# ============================================================================
# Species-Level Simulation (Multiple Trees per Species)
# ============================================================================

simulate_one_species <- function(species, scale, n_trees, tag_offset, growth_multipliers, base_params) {
    p_species <- scale_species_params(base_params, scale)
    rbindlist(lapply(seq_len(n_trees), function(i) {
        tag <- tag_offset + i
        simulate_one_tree(tag, species, growth_multipliers, p_species)
    }))
}

# ============================================================================
# Full Dataset Assembly
# ============================================================================

# Generate forest dataset

tag_offset <- 0L
dt_list <- vector("list", nrow(species_table))
for (i in seq_len(nrow(species_table))) {
    sp <- species_table[i]
    dt_list[[i]] <- simulate_one_species(sp$Species, sp$Scale, as.integer(sp$n_trees), tag_offset, growth_multipliers, params)
    tag_offset <- tag_offset + as.integer(sp$n_trees)
}
dt <- rbindlist(dt_list)

dt <- dt[,.(Species, Tag, OriginalStemID, TrueStemID, CensusID, DBH_true, CensusInterval, CensusDate)][order(Species, Tag, OriginalStemID, CensusID)]
dt[, DBH := DBH_true]
dt[, DBH_true := NULL]
dt[CensusID < 7, TrueStemID := NA_integer_]

dt

################################################################################
### DATA PROCESSING AND EXPORT
################################################################################

# Display summary statistics
cat("Simulation complete. Data summary:\n")
cat(sprintf("Total stems: %d\n", length(unique(dt$Tag))))
cat(sprintf("Total observations: %d\n", nrow(dt[!is.na(DBH)])))
cat("Observations per species:\n")
dt[, .N, by = Species]
dt[, ExactDate := as.Date("1980-01-01") + round(CensusDate * 365.25)]
dt[, CensusDate := NULL]
dt[, .(Species, Tag, OriginalStemID, TrueStemID, CensusID, ExactDate, DBH)][order(Species, Tag, OriginalStemID, CensusID)]

# Generate filename indicator and export main dataset
# data_filename <- build_file_name("data", "simulated_data_1")

# # Export simulation parameters to text file
# params_filename <- here("data_simulation", "data", sprintf("simulation_params_%s.txt", scaling_indicator))
# capture.output(str(params), file = params_filename)

fwrite(dt, here("data_simulation", "data", "simulated_data_1.csv"))
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
    pdf_file_species <- here("data_simulation", "data", "simulated_data_1.pdf")
    pdf(pdf_file_species, width = 8, height = 6)

    # Get unique species
    species_list <- unique(dt$Species)

    # Create one plot per species (all tags)
    for (species_name in species_list) {
        species_data <- dt[Species == species_name]

        gg <- ggplot2::ggplot(
            species_data,
            ggplot2::aes(x = CensusID, y = DBH, group = interaction(OriginalStemID, Tag))
        ) +
            ggplot2::geom_line(na.rm = TRUE, alpha = 0.7) +
            ggplot2::geom_point(size = 1.5, na.rm = TRUE, alpha = 0.7) +
            ggplot2::theme_minimal() +
            ggplot2::labs(
                title = paste("Growth Trajectories -", species_name),
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
    pdf_file_tags <- here("data_simulation", "data", "simulated_data_tag_level_trajectories_1.pdf")
    pdf(pdf_file_tags, width = 8, height = 6)

    # Create one plot per tag
    for (tag_id in unique(dt$Tag)) {
        tag_data <- dt[Tag == tag_id]
        species_name <- unique(tag_data$Species)

        gg <- ggplot2::ggplot(
            tag_data,
            ggplot2::aes(x = as.factor(CensusID), y = DBH, group = interaction(OriginalStemID), color = factor(OriginalStemID))
        ) +
            ggplot2::geom_line(na.rm = TRUE, alpha = 0.8) +
            ggplot2::geom_point(size = 2, na.rm = TRUE, alpha = 0.8) +
            ggplot2::theme_minimal() +
            ggplot2::labs(
                title = sprintf("Tag %d - %s", tag_id, species_name),
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

# The simulation has generated synthetic forest census data with the following features:
# - Multiple species with different growth rates and scaling
# - Multi-stem trees with recruitment and mortality
# - Measurement error in DBH observations
# - Variable census intervals with random noise
# - Stem identity masking for early censuses
# - Diagnostic plots for validation
