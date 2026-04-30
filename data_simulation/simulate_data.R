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

dt <- dt[, .(Species, Tag, OriginalStemID, TrueStemID, CensusID, DBH_true, CensusDate)][order(Species, Tag, OriginalStemID, CensusID)]
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

################################################################################
### INCLUDE POTENTIAL PROBLEMATIC STEMS FOR TESTING
################################################################################

setorder(dt, Species, Tag, OriginalStemID, CensusID)

## add a main stem from census 6 to 9 only
problematic_stem <- data.table(
    Species = "sp1",
    Tag = max(dt$Tag) + 1L,
    OriginalStemID = 1L,
    TrueStemID = 1L,
    CensusID = 6:9,
    ExactDate = as.Date("1980-01-01") + round(c(25, 30, 35, 40) * 365.25),
    DBH = c(30, 35, 40, 45)
)
problematic_stem <- problematic_stem[CensusID < params$mask$anchor_start_census, TrueStemID := NA_integer_]

## add a main stem with census 1 to 5 only
problematic_stem2 <- data.table(
    Species = "sp2",
    Tag = max(dt$Tag) + 2L,
    OriginalStemID = 1L,
    TrueStemID = 1L,
    CensusID = 1:5,
    ExactDate = as.Date("1980-01-01") + round(c(0, 5, 10, 15, 20) * 365.25),
    DBH = c(10, 15, 20, 25, 30)
)
problematic_stem2 <- problematic_stem2[CensusID < params$mask$anchor_start_census, TrueStemID := NA_integer_]

### add a tag with multiple OriginalStemID from census 6 to 9 only
problematic_stem3 <- data.table(
    Species = "sp3",
    Tag = max(dt$Tag) + 3L,
    OriginalStemID = rep(1:2, each = 4),
    TrueStemID = rep(1:2, each = 4),
    CensusID = rep(6:9, times = 2),
    ExactDate = as.Date("1980-01-01") + round(rep(c(25, 30, 35, 40), times = 2) * 365.25),
    DBH = c(20, 25, 30, 35, 15, 20, 25, 30)
)
problematic_stem3 <- problematic_stem3[CensusID < params$mask$anchor_start_census, TrueStemID := NA_integer_]

## add a tag with multiple stems from census 1 to 5 only
problematic_stem4 <- data.table(
    Species = "sp1",
    Tag = max(dt$Tag) + 4L,
    OriginalStemID = rep(1:2, each = 5),
    TrueStemID = rep(1:2, each = 5),
    CensusID = rep(1:5, times = 2),
    ExactDate = as.Date("1980-01-01") + round(rep(c(0, 5, 10, 15, 20), times = 2) * 365.25),
    DBH = c(12, 18, 24, 30, 36, 8, 14, 20, 26, 32)
)
problematic_stem4 <- problematic_stem4[CensusID < params$mask$anchor_start_census, TrueStemID := NA_integer_]

# Combine problematic stems with main dataset
dt_complete <- rbind(dt,
    problematic_stem,
    problematic_stem2,
    problematic_stem3,
    problematic_stem4,
    use.names = TRUE
)

dt_complete

tag_43 <- data.table(
    Species = "sp2",
    Tag = 43L,
    OriginalStemID = NA_integer_,
    TrueStemID = NA_integer_,
    CensusID = NA_integer_,
    ExactDate = as.Date("2020-03-09"),
    DBH = NA_real_
)

tag_44 <- data.table(
    Species = "sp3",
    Tag = 44L,
    OriginalStemID = 1L,
    TrueStemID = 1L,
    CensusID = 9L,
    ExactDate = as.Date("2020-03-09"),
    DBH = 4
)

tag_45 <- data.table(
    Species = "sp1",
    Tag = 45L,
    OriginalStemID = 1L,
    TrueStemID = c(1L, NA, NA),
    CensusID = 7:9,
    ExactDate = as.Date(c("2010-01-18", "2015-02-10", "2022-02-15")),
    DBH = c(1.1, NA, NA)
)

tag_46 <- data.table(
    Species = "sp1",
    Tag = 46L,
    OriginalStemID = 1L,
    TrueStemID = c(1L, 1L, NA),
    CensusID = 7:9,
    ExactDate = as.Date(c("2010-01-18", "2015-02-10", "2022-02-15")),
    DBH = c(1.1, 2.1, NA)
)

tag_47 <- data.table(
    Species = "sp1",
    Tag = 47L,
    OriginalStemID = 1L,
    TrueStemID = c(1L, 1L, 1L),
    CensusID = 7:9,
    ExactDate = as.Date(c("2010-01-18", "2015-02-10", "2022-02-15")),
    DBH = c(1.1, 2.1, 3.1)
)

#    CensusID    Tag TreeID StemTag StemID   DBH TrueStemID
#       <int> <char> <fctr>  <fctr> <char> <num>     <char>
# 1:        1 012370  10393    <NA>  10393   1.5       <NA>
# 2:        2 012370  10393    <NA> 492292    NA       <NA>
# 3:        3 012370  10393    <NA> 492292    NA       <NA>
# 4:        4 012370  10393    <NA> 492292    NA       <NA>
# 5:        5 012370  10393    <NA> 492292    NA       <NA>
# 6:        6 012370  10393    <NA> 492292    NA       <NA>
# 7:        7 012370  10393    <NA> 492292    NA       <NA>
# 8:        8 012370  10393    <NA> 492292    NA       <NA>
# 9:        9 012370  10393    <NA> 492292    NA       <NA>

tag_48 <- data.table(
    Species = "sp2",
    Tag = 48L,
    OriginalStemID = rep(NA_character_, 9),
    TrueStemID = rep(NA_character_, 9),
    CensusID = 1:9,
    ExactDate = as.Date(c(
        "1980-01-01", "1984-09-08", "1989-09-18", "1994-09-08",
        "1999-07-09", "2004-07-29", "2009-08-15", "2014-08-02", "2019-07-20"
    )),
    DBH = c(1.5, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_)
)

#    CensusID    Tag TreeID StemTag StemID   DBH TrueStemID
#       <int> <char> <fctr>  <fctr> <char> <num>     <char>
# 1:        1 007745   7620    <NA>   7620     3       7620
# 2:        1 007745   7620    <NA> 491627     3       7620
# 3:        1 007745   7620    <NA> 619926     2       7620
# 4:        1 007745   7620    <NA> 695839     2       7620
# 5:        1 007745   7620    <NA> 755314     1       7620
# 6:        1 007745   7620    <NA> 802299     1       7620
# 7:        1 007745   7620    <NA> 839676     1       7620
# 8:        1 007745   7620    <NA> 869382     1       7620
# 9:        2 007745   7620    <NA>   7620    NA       <NA>

tag_49 <- data.table(
    Species = "sp2",
    Tag = 49L,
    OriginalStemID = rep(NA_character_, 9),
    TrueStemID = c("1", "1", "1", "1", "1", "1", "1", "1", NA_character_),
    CensusID = c(1, 1, 1, 1, 1, 1, 1, 1, 2),
    ExactDate = as.Date(c(
        "1980-01-01", "1980-01-01", "1980-01-01", "1980-01-01",
        "1980-01-01", "1980-01-01", "1980-01-01", "1980-01-01", "1984-09-08"
    )),
    DBH = c(3, 3, 2, 2, 1, 1, 1, 1, NA_real_)
)

# unique(dt_complete[, .(CensusID, ExactDate)])
## get the min exactdate per census
# dt_complete[, .(MinDate = min(ExactDate, na.rm = TRUE)), by = CensusID]$MinDate

#    CensusID    Tag TreeID StemTag StemID   DBH TrueStemID
#       <int> <char> <fctr>  <fctr> <char> <num>     <char>
# 1:        1 007745   7620    <NA>   7620     3       7620
# 2:        2 007745   7620    <NA>   7620    NA       <NA>
# 3:        1 007745   7620    <NA> 491627     3       7620
# 4:        1 007745   7620    <NA> 619926     2       7620
# 5:        1 007745   7620    <NA> 695839     2       7620
# 6:        1 007745   7620    <NA> 755314     1       7620
# 7:        1 007745   7620    <NA> 802299     1       7620
# 8:        1 007745   7620    <NA> 839676     1       7620
# 9:        1 007745   7620    <NA> 869382     1       7620

tag_50 <- data.table(
    Species = "sp2",
    Tag = 50L,
    OriginalStemID = rep(NA_character_, 9),
    TrueStemID = c("1", NA_character_, "1", "1", "1", "1", "1", "1", "1"),
    CensusID = c(1, 2, 1, 1, 1, 1, 1, 1, 1),
    ExactDate = as.Date(c(
        "1980-01-01", "1984-09-08", "1980-01-01", "1980-01-01",
        "1980-01-01", "1980-01-01", "1980-01-01", "1980-01-01", "1980-01-01"
    )),
    DBH = c(3, NA_real_, 3, 2, 2, 1, 1, 1, 1)
)

tag_51 <- data.table(
    Species = "sp2",
    Tag = 51L,
    OriginalStemID = rep(NA_character_, 9),
    TrueStemID = rep(NA_character_, 9),
    CensusID = 1:9,
    ExactDate = as.Date(c(
        "1980-01-01", "1984-09-08", "1989-09-18", "1994-09-08",
        "1999-07-09", "2004-07-29", "2009-08-15", "2014-08-02", "2019-07-20"
    )),
    DBH = c(NA_real_, NA_real_, 1, 1, 1, 1, 1, 1, 1)
)

tag_52 <- data.table(
    Species = "sp2",
    Tag = 52L,
    OriginalStemID = rep(NA_character_, 9),
    TrueStemID = rep(NA_character_, 9),
    CensusID = 1:9,
    ExactDate = as.Date(c(
        "1980-01-01", "1984-09-08", "1989-09-18", "1994-09-08",
        "1999-07-09", "2004-07-29", "2009-08-15", "2014-08-02", "2019-07-20"
    )),
    DBH = c(NA_real_, NA_real_, 1, 1, 1, 1, NA_real_, NA_real_, NA_real_)
)

tag_53 <- data.table(
    Species = "sp2",
    Tag = 53L,
    OriginalStemID = rep(NA_character_, 9),
    TrueStemID = rep(NA_character_, 9),
    CensusID = 1:9,
    ExactDate = as.Date(c(
        "1980-01-01", "1984-09-08", "1989-09-18", "1994-09-08",
        "1999-07-09", "2004-07-29", "2009-08-15", "2014-08-02", "2019-07-20"
    )),
    DBH = c(1, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_)
)

#        Tag CensusID  ExactDate   DBH OriginalStemID TrueStemID
#     <char>    <int>     <Date> <num>         <char>     <char>
#  1: 140729        1 1981-12-04   1.0         122561       <NA>
#  2: 140729        2 1985-04-16   1.5         122561       <NA>
#  3: 140729        3 1990-08-03   1.3         122561       <NA>
#  4: 140729        4 1995-04-25   1.2         527958       <NA>
#  5: 140729        4 1995-04-25   1.1         642866       <NA>
#  6: 140729        5 2000-04-28   1.3         714907       <NA>
#  7: 140729        5 2000-04-28   1.2         770852       <NA>
#  8: 140729        6 2005-04-29   1.3         815104       <NA>
#  9: 140729        7 2010-04-21    NA         850376     850376
# 10: 140729        8 2015-06-16    NA         850376     850376
# 11: 140729        9 2023-06-15    NA         850376     850376

tag_54 <- data.table(
    Species = "sp1", 
    Tag = 140729L,
    OriginalStemID = c(122561, 122561, 122561, 527958, 642866, 714907, 770852, 815104, 850376, 850376, 850376),
    TrueStemID = c(NA, NA, NA, NA, NA, NA, NA, NA, 850376, 850376, 850376),
    CensusID = c(1L, 2L, 3L, 4L, 4L, 5L, 5L, 6L, 7L, 8L, 9L),
    ExactDate = as.Date(c(
        "1981-12-04", "1985-04-16", "1990-08-03", "1995-04-25", "1995-04-25", "2000-04-28", "2000-04-28", "2005-04-29", "2010-04-21", "2015-06-16", "2023-06-15"
    )),
    DBH = c(1.0, 1.5, 1.3, 1.2, 1.1, 1.3, 1.2, 1.3, NA_real_, NA_real_, NA_real_)
)

#       Tag CensusID  ExactDate   DBH OriginalStemID TrueStemID
#    <char>    <int>     <Date> <num>         <char>     <lgcl>
# 1: 176772        1 1982-04-19   2.0         153339         NA
# 2: 176772        1 1982-04-19   1.5         537473         NA
# 3: 176772        2 1985-04-21   1.5         649031         NA
# 4: 176772        3 1990-06-05    NA         649031         NA

tag_55 <- data.table(
    Species = "sp1", 
    Tag = 176772L,
    OriginalStemID = c(153339, 537473, 649031, 649031),
    TrueStemID = c(NA, NA, NA, NA),
    CensusID = c(1L, 1L, 2L, 3L),
    ExactDate = as.Date(c(
        "1982-04-19", "1982-04-19", "1985-04-21", "1990-06-05"
    )),
    DBH = c(2.0, 1.5, 1.5, NA_real_)
)

#       Tag CensusID  ExactDate   DBH OriginalStemID TrueStemID
#    <char>    <int>     <Date> <num>         <char>     <char>
# 1: 050569        1 1981-06-11   2.0          40853       <NA>
# 2: 050569        2 1985-07-03   2.5          40853       <NA>
# 3: 050569        3 1990-10-04   3.4         501696       <NA>
# 4: 050569        3 1990-10-04   2.0         626064       <NA>
# 5: 050569        4 1995-07-14   3.6         700858       <NA>
# 6: 050569        5 2000-06-25   3.6         700858       <NA>
# 7: 050569        6 2005-08-17    NA         700858       <NA>
# 8: 050569        7 2010-06-29    NA         700858     700858

tag_56 <- data.table(
    Species = "sp1", 
    Tag = 50569L,
    OriginalStemID = c(40853, 40853, 501696, 626064, 700858, 700858, 700858, 700858),
    TrueStemID = c(NA, NA, NA, NA, NA, NA, NA, 700858),
    CensusID = c(1L, 2L, 3L, 3L, 4L, 5L, 6L, 7L),
    ExactDate = as.Date(c(
        "1981-06-11", "1985-07-03", "1990-10-04", "1990-10-04", "1995-07-14", "2000-06-25", "2005-08-17", "2010-06-29"
    )),
    DBH = c(2.0, 2.5, 3.4, 2.0, 3.6, 3.6, NA_real_, NA_real_)
)

#       Tag CensusID  ExactDate   DBH OriginalStemID TrueStemID
#    <char>    <int>     <Date> <num>         <char>     <char>
# 1: 531093        4 1995-06-02   1.2         319129       <NA>
# 2: 531093        5 2000-05-31   1.8         319129       <NA>
# 3: 531093        6 2005-06-27   2.1         319129       <NA>
# 4: 531093        7 2010-06-15   2.1         319129     319129
# 5: 531093        8 2015-07-30   2.2         319129     319129
# 6: 531093        9 2023-01-25   2.8         319129     319129

tag_57 <- data.table(
    Species = "sp1", 
    Tag = 531093L,
    OriginalStemID = rep(319129, 6),
    TrueStemID = c(NA, NA, NA, 319129, 319129, 319129),
    CensusID = c(4L, 5L, 6L, 7L, 8L, 9L),
    ExactDate = as.Date(c(
        "1995-06-02", "2000-05-31", "2005-06-27", "2010-06-15", "2015-07-30", "2023-01-25"
    )),
    DBH = c(1.2, 1.8, 2.1, 2.1, 2.2, 2.8)
)

#        Tag CensusID  ExactDate   DBH OriginalStemID TrueStemID
#     <char>    <int>     <Date> <num>         <char>     <char>
#  1: 082942        1 1981-07-26   1.5          69400       <NA>
#  2: 082942        1 1981-07-26   1.0         511119       <NA>
#  3: 082942        2 1985-06-21   1.5         632151       <NA>
#  4: 082942        2 1985-06-21   1.0         705917       <NA>
#  5: 082942        3 1990-10-25   1.6         763499       <NA>
#  6: 082942        3 1990-10-25   1.1         809048       <NA>
#  7: 082942        4 1995-07-28   1.6         845302       <NA>
#  8: 082942        5 2000-07-12   1.7         845302       <NA>
#  9: 082942        6 2005-07-19    NA         845302       <NA>
# 10: 082942        7 2010-07-21    NA         845302     845302

tag_58 <- data.table(
    Species = "sp1", 
    Tag = 82942L,
    OriginalStemID = c(69400, 511119, 632151, 705917, 763499, 809048, 845302, 845302, 845302, 845302),
    TrueStemID = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, 845302),
    CensusID = c(1L, 1L, 2L, 2L, 3L, 3L, 4L, 5L, 6L, 7L),
    ExactDate = as.Date(c(
        "1981-07-26", "1981-07-26", "1985-06-21", "1985-06-21", "1990-10-25", "1990-10-25", "1995-07-28", "2000-07-12", "2005-07-19", "2010-07-21"
    )),
    DBH = c(1.5, 1.0, 1.5, 1.0, 1.6, 1.1, 1.6, 1.7, NA_real_, NA_real_)
)

#        Tag CensusID  ExactDate   DBH OriginalStemID TrueStemID
#     <char>    <int>     <Date> <num>         <char>     <char>
#  1: 111885        1 1981-11-04   4.5          96616       <NA>
#  2: 111885        2 1985-05-17   4.5         519942       <NA>
#  3: 111885        2 1985-05-17   1.5         637884       <NA>
#  4: 111885        3 1990-08-31   4.1         710724       <NA>
#  5: 111885        3 1990-08-31   1.9         767442       <NA>
#  6: 111885        4 1995-05-25   4.0         812265       <NA>
#  7: 111885        4 1995-05-25   1.8         847997       <NA>
#  8: 111885        5 2000-05-18   1.9         876277       <NA>
#  9: 111885        6 2005-06-16    NA         899199       <NA>
# 10: 111885        7 2010-05-19    NA         899199     899199
# 11: 111885        8 2015-08-11    NA         899199     899199

tag_59 <- data.table(
    Species = "sp1", 
    Tag = 111885L,
    OriginalStemID = c(96616, 519942, 637884, 710724, 767442, 812265, 847997, 876277, 899199, 899199, 899199),
    TrueStemID = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, 899199, 899199),
    CensusID = c(1L, 2L, 2L, 3L, 3L, 4L, 4L, 5L, 6L, 7L, 8L),
    ExactDate = as.Date(c(
        "1981-11-04", "1985-05-17", "1985-05-17", "1990-08-31", "1990-08-31", "1995-05-25", "1995-05-25", "2000-05-18", "2005-06-16", "2010-05-19", "2015-08-11"
    )),
    DBH = c(4.5, 4.5, 1.5, 4.1, 1.9, 4.0, 1.8, 1.9, NA_real_, NA_real_, NA_real_)
)

#        Tag CensusID  ExactDate   DBH OriginalStemID TrueStemID
#     <char>    <int>     <Date> <num>         <char>     <char>
#  1: 061350        1 1981-06-27   2.0          49762       <NA>
#  2: 061350        1 1981-06-27   1.5         504881       <NA>
#  3: 061350        2 1985-07-01   2.5         628167       <NA>
#  4: 061350        2 1985-07-01   1.5         702613       <NA>
#  5: 061350        2 1985-07-01   1.0         760769       <NA>
#  6: 061350        3 1990-11-15   3.7         806789       <NA>
#  7: 061350        3 1990-11-15   3.5         843405       <NA>
#  8: 061350        3 1990-11-15   2.7         872475       <NA>
#  9: 061350        3 1990-11-15   1.8         895994       <NA>
# 10: 061350        4 1995-07-18   2.3         915624       <NA>
# 11: 061350        4 1995-07-18   2.3         931429       <NA>
# 12: 061350        4 1995-07-18   1.9         944462       <NA>
# 13: 061350        4 1995-07-18   1.8         954604       <NA>
# 14: 061350        5 2000-07-05   3.4         962914       <NA>
# 15: 061350        5 2000-07-05   2.7         969576       <NA>
# 16: 061350        5 2000-07-05   2.2         974531       <NA>
# 17: 061350        6 2005-07-13    NA         962914       <NA>
# 18: 061350        7 2010-07-22    NA         974531     974531

tag_60 <- data.table(
    Species = "sp1", 
    Tag = 61350L,
    OriginalStemID = c(49762, 504881, 628167, 702613, 760769, 806789, 843405, 872475, 895994, 915624, 931429, 944462, 954604, 962914, 969576, 974531, 962914, 974531),
    TrueStemID = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 974531),
    CensusID = c(1L, 1L, 2L, 2L, 2L, 3L, 3L, 3L, 3L, 4L, 4L, 4L, 4L, 5L, 5L, 5L, 6L, 7L),
    ExactDate = as.Date(c(
        "1981-06-27", "1981-06-27", "1985-07-01", "1985-07-01", "1985-07-01", "1990-11-15", "1990-11-15", "1990-11-15", "1990-11-15", "1995-07-18", "1995-07-18", "1995-07-18", "1995-07-18", "2000-07-05", "2000-07-05", "2000-07-05", "2005-07-13", "2010-07-22"
    )),
    DBH = c(2.0, 1.5, 2.5, 1.5, 1.0, 3.7, 3.5, 2.7, 1.8, 2.3, 2.3, 1.9, 1.8, 3.4, 2.7, 2.2, NA_real_, NA_real_)
)

#        Tag CensusID  ExactDate   DBH OriginalStemID TrueStemID
#     <char>    <int>     <Date> <num>         <char>     <char>
#  1: 247222        1 1982-07-18   5.0         218861       <NA>
#  2: 247222        1 1982-07-18   1.0         558508       <NA>
#  3: 247222        2 1985-01-26   4.0         663129       <NA>
#  4: 247222        2 1985-01-26   1.5         731991       <NA>
#  5: 247222        3 1990-03-10   5.0         784793       <NA>
#  6: 247222        3 1990-03-10   1.5         826574       <NA>
#  7: 247222        4 1995-02-08   1.7         860033       <NA>
#  8: 247222        5 2000-01-26   1.7         886324       <NA>
#  9: 247222        5 2000-01-26   1.0         907848       <NA>
# 10: 247222        6 2005-02-02    NA         886324       <NA>
# 11: 247222        7 2010-02-03    NA         907848     907848

tag_61 <- data.table(
    Species = "sp1", 
    Tag = 247222L,
    OriginalStemID = c(218861, 558508, 663129, 731991, 784793, 826574, 860033, 886324, 907848, 886324, 907848),
    TrueStemID = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 907848),
    CensusID = c(1L, 1L, 2L, 2L, 3L, 3L, 4L, 5L, 5L, 6L, 7L),
    ExactDate = as.Date(c(
        "1982-07-18", "1982-07-18", "1985-01-26", "1985-01-26", "1990-03-10", "1990-03-10", "1995-02-08", "2000-01-26", "2000-01-26", "2005-02-02", "2010-02-03"
    )),
    DBH = c(5.0, 1.0, 4.0, 1.5, 5.0, 1.5, 1.7, 1.7, 1.0, NA_real_, NA_real_)
)

#        Tag CensusID  ExactDate   DBH OriginalStemID TrueStemID
#     <char>    <int>     <Date> <num>         <char>     <char>
#  1: 252057        1 1982-07-19   9.3         223580       <NA>
#  2: 252057        2 1985-01-23   9.5         560162       <NA>
#  3: 252057        2 1985-01-23   1.0         664302       <NA>
#  4: 252057        3 1990-03-09   9.7         732980       <NA>
#  5: 252057        3 1990-03-09   1.6         785617       <NA>
#  6: 252057        4 1995-01-31   7.9         827227       <NA>
#  7: 252057        4 1995-01-31   2.1         860569       <NA>
#  8: 252057        4 1995-01-31   2.4         886751       <NA>
#  9: 252057        4 1995-01-31   1.3         908222       <NA>
# 10: 252057        5 2000-02-03   9.8         925701       <NA>
# 11: 252057        5 2000-02-03   2.5         940180       <NA>
# 12: 252057        5 2000-02-03   2.0         951248       <NA>
# 13: 252057        6 2005-02-02   8.9         960314       <NA>
# 14: 252057        6 2005-02-02   2.0         967546       <NA>
# 15: 252057        7 2010-02-04    NA         967546     967546
# 16: 252057        8 2015-02-27    NA         967546     967546

tag_62 <- data.table(
    Species = "sp1", 
    Tag = 252057L,
    OriginalStemID = c(223580, 560162, 664302, 732980, 785617, 827227, 860569, 886751, 908222, 925701, 940180, 951248, 960314, 967546, 967546, 967546),
    TrueStemID = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 967546, 967546),
    CensusID = c(1L, 2L, 2L, 3L, 3L, 4L, 4L, 4L, 4L, 5L, 5L, 5L, 6L, 6L, 7L, 8L),
    ExactDate = as.Date(c(
        "1982-07-19", "1985-01-23", "1985-01-23", "1990-03-09", "1990-03-09", "1995-01-31", "1995-01-31", "1995-01-31", "1995-01-31", "2000-02-03", "2000-02-03", "2000-02-03", "2005-02-02", "2005-02-02", "2010-02-04", "2015-02-27"
    )),
    DBH = c(9.3, 9.5, 1.0, 9.7, 1.6, 7.9, 2.1, 2.4, 1.3, 9.8, 2.5, 2.0, 8.9, 2.0, NA_real_, NA_real_)
)

#       Tag CensusID  ExactDate   DBH OriginalStemID TrueStemID
#    <char>    <int>     <Date> <num>         <char>     <lgcl>
# 1: 002747        1 1981-10-08  20.6           2738         NA
# 2: 002747        2 1985-06-21  21.8           2738         NA
# 3: 002747        3 1990-12-18   2.7         491245         NA
# 4: 002747        3 1990-12-18   1.4         619581         NA
# 5: 002747        3 1990-12-18   1.0         695541         NA
# 6: 002747        4 1995-08-21  22.2         755073         NA
# 7: 002747        4 1995-08-21   2.7         802079         NA
# 8: 002747        5 2000-07-19    NA         755073         NA

tag_63 <- data.table(
    Species = "sp1", 
    Tag = 2747L,
    OriginalStemID = c(2738, 2738, 491245, 619581, 695541, 755073, 802079, 755073),
    TrueStemID = c(NA, NA, NA, NA, NA, NA, NA, NA),
    CensusID = c(1L, 2L, 3L, 3L, 3L, 4L, 4L, 5L),
    ExactDate = as.Date(c(
        "1981-10-08", "1985-06-21", "1990-12-18", "1990-12-18", "1990-12-18", "1995-08-21", "1995-08-21", "2000-07-19"
    )),
    DBH = c(20.6, 21.8, 2.7, 1.4, 1.0, 22.2, 2.7, NA_real_)
)

#       Tag CensusID  ExactDate   DBH OriginalStemID TrueStemID
#    <char>    <int>     <Date> <num>         <char>     <char>
# 1: 306322        2 1985-07-23   1.0         262595       <NA>
# 2: 306322        3 1990-10-25    NA         574524       <NA>
# 3: 306322        4 1995-08-09   1.1         574524       <NA>
# 4: 306322        5 2000-08-02   1.3         574524       <NA>
# 5: 306322        6 2005-09-01   1.3         574524       <NA>
# 6: 306322        7 2010-07-22   1.4         574524     574524
# 7: 306322        8 2015-10-20    NA         574524     574524
# 8: 306322        9 2023-05-02    NA         574524     574524

tag_64 <- data.table(
    Species = "sp1", 
    Tag = 306322L,
    OriginalStemID = rep(262595, 8),
    TrueStemID = c(NA, NA, NA, NA, NA, 574524, 574524, 574524),
    CensusID = c(2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L),
    ExactDate = as.Date(c(
        "1985-07-23", "1990-10-25", "1995-08-09", "2000-08-02", "2005-09-01", "2010-07-22", "2015-10-20", "2023-05-02"
    )),
    DBH = c(1.0, NA_real_, 1.1, 1.3, 1.3, 1.4, NA_real_, NA_real_)
)

#        Tag CensusID  ExactDate   DBH OriginalStemID TrueStemID
#     <char>    <int>     <Date> <num>         <char>     <char>
#  1: 071698        1 1981-07-11   2.0          59081       <NA>
#  2: 071698        2 1985-06-22   2.0          59081       <NA>
#  3: 071698        3 1990-10-18   2.8         507755       <NA>
#  4: 071698        3 1990-10-18   2.1         629999       <NA>
#  5: 071698        3 1990-10-18   1.3         704141       <NA>
#  6: 071698        4 1995-07-10   2.7         762029       <NA>
#  7: 071698        4 1995-07-10   1.3         807829       <NA>
#  8: 071698        4 1995-07-10   2.6         844274       <NA>
#  9: 071698        5 2000-06-29   3.0         873203       <NA>
# 10: 071698        6 2005-07-21    NA         873203       <NA>
# 11: 071698        7 2010-07-01    NA         873203     873203

tag_65 <- data.table(
    Species = "sp1", 
    Tag = 71698L,
    OriginalStemID = c(59081, 59081, 507755, 629999, 704141, 762029, 807829, 844274, 873203, 873203, 873203),
    TrueStemID = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 873203),
    CensusID = c(1L, 2L, 3L, 3L, 3L, 4L, 4L, 4L, 5L, 6L, 7L),
    ExactDate = as.Date(c(
        "1981-07-11", "1985-06-22", "1990-10-18", "1990-10-18", "1990-10-18", "1995-07-10", "1995-07-10", "1995-07-10", "2000-06-29", "2005-07-21", "2010-07-01"
    )),
    DBH = c(2.0, 2.0, 2.8, 2.1, 1.3, 2.7, 1.3, 2.6, 3.0, NA_real_, NA_real_)
)

#        Tag CensusID  ExactDate   DBH OriginalStemID TrueStemID
#     <char>    <int>     <Date> <num>         <char>     <char>
#  1: 229214        1 1982-07-10   1.0         202251       <NA>
#  2: 229214        2 1985-03-13   1.0         202251       <NA>
#  3: 229214        3 1990-07-18   1.2         553338       <NA>
#  4: 229214        3 1990-07-18   1.1         659595       <NA>
#  5: 229214        4 1995-03-21   1.1         729015       <NA>
#  6: 229214        4 1995-03-21   1.1         782324       <NA>
#  7: 229214        5 2000-02-25   1.2         824512       <NA>
#  8: 229214        6 2005-03-28   1.2         824512       <NA>
#  9: 229214        7 2010-03-16    NA         858289     858289
# 10: 229214        8 2015-04-01    NA         858289     858289
# 11: 229214        9 2022-06-09    NA         858289     858289

tag_66 <- data.table(
    Species = "sp1", 
    Tag = 229214L,
    OriginalStemID = c(202251, 202251, 553338, 659595, 729015, 782324, 824512, 824512, 858289, 858289, 858289),
    TrueStemID = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, 858289, 858289),
    CensusID = c(1L, 2L, 3L, 3L, 4L, 4L, 5L, 6L, 7L, 8L, 9L),
    ExactDate = as.Date(c(
        "1982-07-10", "1985-03-13", "1990-07-18", "1990-07-18", "1995-03-21", "1995-03-21", "2000-02-25", "2005-03-28", "2010-03-16", "2015-04-01", "2022-06-09"
    )),
    DBH = c(1.0, 1.0, 1.2, 1.1, 1.1, 1.1, 1.2, 1.2, NA_real_, NA_real_, NA_real_)
)

dt_complete_extra <- rbind(
    dt_complete, tag_43, tag_44, tag_45, tag_46,
    tag_47, tag_48, tag_49, tag_50,
    tag_51, tag_52, tag_53, tag_54, tag_55, tag_56, tag_57, tag_58,
    tag_59, tag_60, tag_61, tag_62, tag_63, tag_64, tag_65, tag_66
)

dt_complete_extra[, ListOfTSM := NA_character_]

#     Species    Tag OriginalStemID TrueStemID CensusID  ExactDate      DBH ListOfTSM
#      <char> <char>          <int>      <int>    <int>     <Date>    <num>    <char>
#  1:     sp1 131222         114524         NA        1 1981-12-09  3.00000      <NA>
#  2:     sp1 131222         114524         NA        2 1985-05-07       NA        MF
#  3:     sp1 131222         525593         NA        3 1990-08-01  5.10000      <NA>
#  4:     sp1 131222         641379         NA        3 1990-08-01  4.40000      <NA>
#  5:     sp1 131222         713657         NA        4 1995-05-09  7.90000      <NA>
#  6:     sp1 131222         713657         NA        5 2000-05-04 13.40000      <NA>
#  7:     sp1 131222         713657         NA        6 2005-05-27 20.40000      <NA>
#  8:     sp1 131222         713657     713657        7 2010-05-05 27.30000      <NA>
#  9:     sp1 131222         713657     713657        8 2015-06-26 33.03964      <NA>
# 10:     sp1 131222         713657     713657        9 2023-06-29 38.85152      <NA>

tag_67 <- data.table(
    Species = "sp1", 
    Tag = 131222L,
    OriginalStemID = c(114524, 114524, 525593, 641379, 713657, 713657, 713657, 713657, 713657, 713657),
    TrueStemID = c(NA, NA, NA, NA, NA, NA, NA, 713657, 713657, 713657),
    CensusID = c(1L, 2L, 3L, 3L, 4L, 5L, 6L, 7L, 8L, 9L),
    ExactDate = as.Date(c(
        "1981-12-09", "1985-05-07", "1990-08-01", "1990-08-01", "1995-05-09", "2000-05-04", "2005-05-27", "2010-05-05", "2015-06-26", "2023-06-29"
    )),
    DBH = c(3.0, NA_real_, 5.1, 4.4, 7.9, 13.4, 20.4, 27.3, 33.03964, 38.85152), 
    ListOfTSM = c(NA, "MF", NA, NA, NA, NA, NA, NA, NA, NA)
)

#     Species    Tag OriginalStemID TrueStemID CensusID  ExactDate   DBH ListOfTSM
#      <char> <char>          <int>      <int>    <int>     <Date> <num>    <char>
#  1:     sp1 166815         145077         NA        1 1982-03-28   2.0      <NA>
#  2:     sp1 166815         145077         NA        2 1985-04-25    NA        MF
#  3:     sp1 166815         145077         NA        3 1990-08-23    NA      <NA>
#  4:     sp1 166815         145077         NA        4 1995-04-05   1.5      <NA>
#  5:     sp1 166815         145077         NA        5 2000-04-05   1.6      <NA>
#  6:     sp1 166815         145077         NA        6 2005-04-26   1.6      <NA>
#  7:     sp1 166815         145077     145077        7 2010-05-06   1.7      <NA>
#  8:     sp1 166815         145077     145077        8 2015-07-01   1.7      <NA>
#  9:     sp1 166815        1038359    1038359        8 2015-07-01   1.0      <NA>
# 10:     sp1 166815         145077     145077        9 2022-09-20    NA      <NA>
# 11:     sp1 166815        1038359    1038359        9 2022-09-20    NA      <NA>

tag_68 <- data.table(
    Species = "sp1",
    Tag = 166815L,
    OriginalStemID = c(145077, 145077, 145077, 145077, 145077, 145077, 145077, 145077, 1038359, 145077, 1038359),
    TrueStemID = c(NA, NA, NA, NA, NA, NA, 145077, 145077, 1038359, 145077, 1038359),
    CensusID = c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 8L, 9L, 9L),
    ExactDate = as.Date(c(
        "1982-03-28", "1985-04-25", "1990-08-23", "1995-04-05", "2000-04-05", "2005-04-26", "2010-05-06", "2015-07-01", "2015-07-01", "2022-09-20", "2022-09-20"
    )),
    DBH = c(2.0, NA_real_, NA_real_, 1.5, 1.6, 1.6, 1.7, 1.7, 1.0, NA_real_, NA_real_),
    ListOfTSM = c(NA_character_, "MF", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
)

#    Species    Tag OriginalStemID TrueStemID CensusID  ExactDate   DBH ListOfTSM
#     <char> <char>          <int>      <int>    <int>     <Date> <num>    <char>
# 1:     sp1 246746         218395         NA        1 1982-07-16   2.0      <NA>
# 2:     sp1 246746         218395         NA        2 1985-01-22    NA        MF
# 3:     sp1 246746         218395         NA        3 1990-03-06   2.1      <NA>
# 4:     sp1 246746         558362         NA        4 1995-01-30   1.9        MF
# 5:     sp1 246746         558362         NA        5 2000-01-24   2.3      <NA>
# 6:     sp1 246746         558362         NA        6 2005-01-28   2.3      <NA>
# 7:     sp1 246746         558362     558362        7 2010-02-01   2.4      <NA>
# 8:     sp1 246746         558362     558362        8 2015-02-27   2.5      <NA>
# 9:     sp1 246746         558362     558362        9 2022-03-28   2.6      <NA>

tag_69 <- data.table(
    Species = "sp1",
    Tag = 246746L,
    OriginalStemID = c(218395, 218395, 218395, 558362, 558362, 558362, 558362, 558362, 558362),
    TrueStemID = c(NA, NA, NA, NA, NA, NA, 558362, 558362, 558362),
    CensusID = c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L),
    ExactDate = as.Date(c(
        "1982-07-16", "1985-01-22", "1990-03-06", "1995-01-30", "2000-01-24", "2005-01-28", "2010-02-01", "2015-02-27", "2022-03-28"
    )),
    DBH = c(2.0, NA_real_, 2.1, 1.9, 2.3, 2.3, 2.4, 2.5, 2.6),
    ListOfTSM = c(NA_character_, "MF", NA_character_, "R", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
)

#     Species    Tag OriginalStemID TrueStemID CensusID  ExactDate   DBH ListOfTSM
#      <char> <char>          <int>      <int>    <int>     <Date> <num>    <char>
#  1:     sp1 249573         221145         NA        1 1982-08-09   2.0      <NA>
#  2:     sp1 249573         221145         NA        2 1985-02-08    NA        MF
#  3:     sp1 249573         221145         NA        3 1990-04-18   2.2      <NA>
#  4:     sp1 249573         221145         NA        4 1995-03-28   2.4      <NA>
#  5:     sp1 249573         221145         NA        5 2000-02-24   2.5      <NA>
#  6:     sp1 249573         221145         NA        6 2005-03-04   2.7      <NA>
#  7:     sp1 249573         559342     559342        7 2010-03-17   2.7      <NA>
#  8:     sp1 249573         663693     663693        7 2010-03-17   2.4      <NA>
#  9:     sp1 249573         559342     559342        8 2015-03-19   2.7      <NA>
# 10:     sp1 249573         663693     663693        8 2015-03-19   2.5      <NA>
# 11:     sp1 249573         559342     559342        9 2022-06-03    NA       R;M
# 12:     sp1 249573         663693     663693        9 2022-06-03    NA      <NA>
# 13:     sp1 249573        1115540    1115540        9 2022-06-03   1.8      <NA>

tag_70 <- data.table(
    Species = "sp1",
    Tag = 249573L,
    OriginalStemID = c(221145, 221145, 221145, 221145, 221145, 221145, 559342, 663693, 559342, 663693, 559342, 663693, 1115540),
    TrueStemID = c(NA, NA, NA, NA, NA, NA, 559342, 663693, 559342, 663693, 559342, 663693, 1115540),
    CensusID = c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 7L, 8L, 8L, 9L, 9L, 9L),
    ExactDate = as.Date(c(
        "1982-08-09", "1985-02-08", "1990-04-18", "1995-03-28", "2000-02-24", "2005-03-04", "2010-03-17", "2010-03-17", "2015-03-19", "2015-03-19", "2022-06-03", "2022-06-03", "2022-06-03"
    )),
    DBH = c(2.0, NA_real_, 2.2, 2.4, 2.5, 2.7, 2.7, 2.4, 2.7, 2.5, NA_real_, NA_real_, 1.8),
    ListOfTSM = c(NA_character_, "MF", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, "R;M", NA_character_, NA_character_)
)

#     Species    Tag OriginalStemID TrueStemID CensusID  ExactDate   DBH ListOfTSM
#      <char> <char>          <int>      <int>    <int>     <Date> <num>    <char>
#  1:     sp1 040831          35106         NA        1 1981-05-23   1.0      <NA>
#  2:     sp1 040831          35106         NA        2 1985-07-06    NA        MF
#  3:     sp1 040831         499935         NA        3 1990-09-29   2.1      <NA>
#  4:     sp1 040831         624958         NA        3 1990-09-29   1.8      <NA>
#  5:     sp1 040831         699955         NA        4 1995-08-07   2.3      <NA>
#  6:     sp1 040831         758643         NA        4 1995-08-07   2.1      <NA>
#  7:     sp1 040831         805019         NA        5 2000-08-07   2.3      <NA>
#  8:     sp1 040831         841935         NA        5 2000-08-07   2.1      <NA>
#  9:     sp1 040831         871240         NA        6 2005-08-19   2.3      <NA>
# 10:     sp1 040831         894933         NA        6 2005-08-19   2.2      <NA>
# 11:     sp1 040831         914733     914733        7 2010-07-22   3.0      <NA>
# 12:     sp1 040831         930640     930640        7 2010-07-22   2.3      <NA>
# 13:     sp1 040831         914733     914733        8 2015-10-12   3.2      <NA>
# 14:     sp1 040831         930640     930640        8 2015-10-12    NA      <NA>
# 15:     sp1 040831        1042996    1042996        8 2015-10-12   1.3      <NA>
# 16:     sp1 040831         914733     914733        9 2023-05-02   4.2      <NA>
# 17:     sp1 040831         930640     930640        9 2023-05-02    NA      <NA>
# 18:     sp1 040831        1042996    1042996        9 2023-05-02    NA      <NA>

tag_71 <- data.table(
    Species = "sp1",
    Tag = 40831L,
    OriginalStemID = c(35106, 35106, 499935, 624958, 699955, 758643, 805019, 841935, 871240, 894933, 914733, 930640, 914733, 930640, 1042996, 914733, 930640, 1042996),
    TrueStemID = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 914733, 930640, 914733, 930640, 1042996, 914733, 930640, 1042996),
    CensusID = c(1L, 2L, 3L, 3L, 4L, 4L, 5L, 5L, 6L, 6L, 7L, 7L, 8L, 8L, 8L, 9L, 9L, 9L),
    ExactDate = as.Date(c(
        "1981-05-23", "1985-07-06", "1990-09-29", "1990-09-29", "1995-08-07", "1995-08-07", "2000-08-07", "2000-08-07", "2005-08-19", "2005-08-19", "2010-07-22", "2010-07-22", "2015-10-12", "2015-10-12", "2015-10-12", "2023-05-02", "2023-05-02", "2023-05-02"
    )),
    DBH = c(1.0, NA_real_, 2.1, 1.8, 2.3, 2.1, 2.3, 2.1, 2.3, 2.2, 3.0, 2.3, 3.2, NA_real_, NA_real_, 4.2, NA_real_, NA_real_),
    ListOfTSM = c(NA_character_, "MF", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
)

#    Species    Tag OriginalStemID TrueStemID CensusID  ExactDate   DBH ListOfTSM
#     <char> <char>          <int>      <int>    <int>     <Date> <num>    <char>
# 1:     sp1 228583         201626         NA        1 1982-07-03   2.0      <NA>
# 2:     sp1 228583         553142         NA        1 1982-07-03   1.0      <NA>
# 3:     sp1 228583         659459         NA        2 1985-03-08   2.0      <NA>
# 4:     sp1 228583         728906         NA        2 1985-03-08   1.0      <NA>
# 5:     sp1 228583         782235         NA        3 1990-07-06   1.4      <NA>
# 6:     sp1 228583         824437         NA        4 1995-03-16   1.5       L;R
# 7:     sp1 228583         824437         NA        5 2000-02-21    NA      <NA>

tag_72 <- data.table(
    Species = "sp1",
    Tag = 228583L,
    OriginalStemID = c(201626, 553142, 659459, 728906, 782235, 824437, 824437),
    TrueStemID = c(NA, NA, NA, NA, NA, NA, NA),
    CensusID = c(1L, 1L, 2L, 2L, 3L, 4L, 5L),
    ExactDate = as.Date(c(
        "1982-07-03", "1982-07-03", "1985-03-08", "1985-03-08", "1990-07-06", "1995-03-16", "2000-02-21"
    )),
    DBH = c(2.0, 1.0, 2.0, 1.0, 1.4, 1.5, NA_real_),
    ListOfTSM = c(NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, "L;R", NA_character_)
)

#     Species    Tag OriginalStemID TrueStemID CensusID  ExactDate   DBH ListOfTSM
#      <char> <char>          <int>      <int>    <int>     <Date> <num>    <char>
#  1:     sp1 031364          26232         NA        1 1981-05-04   4.5      <NA>
#  2:     sp1 031364         497233         NA        1 1981-05-04   1.5      <NA>
#  3:     sp1 031364         623343         NA        2 1985-07-12   5.1      <NA>
#  4:     sp1 031364         698622         NA        2 1985-07-12   3.0      <NA>
#  5:     sp1 031364         757565         NA        3 1990-11-16   5.2      <NA>
#  6:     sp1 031364         804157         NA        3 1990-11-16   3.7      <NA>
#  7:     sp1 031364         841234         NA        4 1995-09-11   5.3      <NA>
#  8:     sp1 031364         870659         NA        4 1995-09-11   4.7      <NA>
#  9:     sp1 031364         894441         NA        5 2000-09-12   4.8       R;M
# 10:     sp1 031364         914318         NA        5 2000-09-12   1.5      <NA>
# 11:     sp1 031364         930303         NA        6 2005-09-01   4.9      <NA>
# 12:     sp1 031364         943590         NA        6 2005-09-01   3.1      <NA>
# 13:     sp1 031364         953862     953862        7 2010-08-21   4.9      <NA>
# 14:     sp1 031364         962322     962322        7 2010-08-21   3.7      <NA>
# 15:     sp1 031364         953862     953862        8 2015-10-21   4.9      <NA>
# 16:     sp1 031364         962322     962322        8 2015-10-21   4.2      <NA>
# 17:     sp1 031364         953862     953862        9 2023-05-02    NA         R
# 18:     sp1 031364         962322     962322        9 2023-05-02   5.4      <NA>

tag_73 <- data.table(
    Species = "sp1",
    Tag = 31364L,
    OriginalStemID = c(26232, 497233, 623343, 698622, 757565, 804157, 841234, 870659, 894441, 914318, 930303, 943590, 953862, 962322, 953862, 962322, 953862, 962322),
    TrueStemID = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 953862, 962322, 953862, 962322, 953862, 962322),
    CensusID = c(1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L, 5L, 5L, 6L, 6L, 7L, 7L, 8L, 8L, 9L, 9L),
    ExactDate = as.Date(c(
        "1981-05-04", "1981-05-04", "1985-07-12", "1985-07-12", "1990-11-16", "1990-11-16", "1995-09-11", "1995-09-11", "2000-09-12", "2000-09-12", "2005-09-01", "2005-09-01", "2010-08-21", "2010-08-21", "2015-10-21", "2015-10-21", "2023-05-02", "2023-05-02"
    )),
    DBH = c(4.5, 1.5, 5.1, 3.0, 5.2, 3.7, 5.3, 4.7, 4.8, 1.5, 4.9, 3.1, 4.9, 3.7, 4.9, 4.2, NA_real_, 5.4),
    ListOfTSM = c(NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, "R;M", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, "R", NA_character_)
)

#    Species    Tag OriginalStemID TrueStemID CensusID  ExactDate   DBH ListOfTSM
#     <char> <char>          <int>      <int>    <int>     <Date> <num>    <char>
# 1:     sp1 193880         169280         NA        1 1982-05-30   1.0      <NA>
# 2:     sp1 193880         169280         NA        2 1985-03-16   1.0      <NA>
# 3:     sp1 193880         542613         NA        3 1990-07-06   1.2         R
# 4:     sp1 193880         542613         NA        4 1995-05-04   1.2      <NA>
# 5:     sp1 193880         542613         NA        5 2000-04-17   1.3      <NA>
# 6:     sp1 193880         542613         NA        6 2005-04-06   1.3      <NA>
# 7:     sp1 193880         542613     542613        7 2010-04-21   1.4      <NA>
# 8:     sp1 193880         542613     542613        8 2015-05-26    NA      <NA>
# 9:     sp1 193880         542613     542613        9 2022-08-19    NA      <NA>

tag_74 <- data.table(
    Species = "sp1",
    Tag = 193880L,
    OriginalStemID = c(169280, 169280, 542613, 542613, 542613, 542613, 542613, 542613, 542613),
    TrueStemID = c(NA, NA, NA, NA, NA, NA, 542613, 542613, 542613),
    CensusID = c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L),
    ExactDate = as.Date(c(
        "1982-05-30", "1985-03-16", "1990-07-06", "1995-05-04", "2000-04-17", "2005-04-06", "2010-04-21", "2015-05-26", "2022-08-19"
    )),
    DBH = c(1.0, 1.0, 1.2, 1.2, 1.3, 1.3, 1.4, NA_real_, NA_real_),
    ListOfTSM = c(NA_character_, NA_character_, "R", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
)

#  Species    Tag OriginalStemID TrueStemID CensusID  ExactDate   DBH ListOfTSM
#      <char> <char>          <int>      <int>    <int>     <Date> <num>    <char>
#  1:     sp1 123315         107666         NA        1 1981-11-21   2.0      <NA>
#  2:     sp1 123315         107666         NA        2 1985-05-24   3.0      <NA>
#  3:     sp1 123315         523406         NA        3 1990-09-11   3.4      <NA>
#  4:     sp1 123315         640078         NA        3 1990-09-11   1.7      <NA>
#  5:     sp1 123315         712576         NA        4 1995-06-14   3.5      <NA>
#  6:     sp1 123315         768941         NA        4 1995-06-14   1.5      <NA>
#  7:     sp1 123315         813512         NA        4 1995-06-14   1.3      <NA>
#  8:     sp1 123315         849035         NA        5 2000-06-01   3.6      <NA>
#  9:     sp1 123315         877137         NA        5 2000-06-01   1.7      <NA>
# 10:     sp1 123315         899940         NA        5 2000-06-01   1.5      <NA>
# 11:     sp1 123315         918872         NA        6 2005-06-30   3.8      <NA>
# 12:     sp1 123315         934224         NA        6 2005-06-30   1.9      <NA>
# 13:     sp1 123315         946579         NA        6 2005-06-30   1.7      <NA>
# 14:     sp1 123315         956374         NA        6 2005-06-30   1.2      <NA>
# 15:     sp1 123315         964327     964327        7 2010-05-25   1.8         R
# 16:     sp1 123315         970516     970516        7 2010-05-25   1.1      <NA>
# 17:     sp1 123315         964327     964327        8 2015-09-01    NA      <NA>
# 18:     sp1 123315         970516     970516        8 2015-09-01    NA      <NA>
# 19:     sp1 123315         964327     964327        9 2022-12-21    NA      <NA>
# 20:     sp1 123315         970516     970516        9 2022-12-21    NA      <NA>

tag_75 <- data.table(
    Species = "sp1",
    Tag = 123315L,
    OriginalStemID = c(107666, 107666, 523406, 640078, 712576, 768941, 813512, 849035, 877137, 899940, 918872, 934224, 946579, 956374, 964327, 970516, 964327, 970516, 964327, 970516),
    TrueStemID = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 964327, 970516, 964327, 970516, 964327, 970516),
    CensusID = c(1L, 2L, 3L, 3L, 4L, 4L, 4L, 5L, 5L, 5L, 6L, 6L, 6L, 6L, 7L, 7L, 8L, 8L, 9L, 9L),
    ExactDate = as.Date(c(
        "1981-11-21", "1985-05-24", "1990-09-11", "1990-09-11", "1995-06-14", "1995-06-14", "1995-06-14", "2000-06-01", "2000-06-01", "2000-06-01", "2005-06-30", "2005-06-30", "2005-06-30", "2005-06-30", "2010-05-25", "2010-05-25", "2015-09-01", "2015-09-01", "2022-12-21", "2022-12-21"
    )),
    DBH = c(2.0, 3.0, 3.4, 1.7, 3.5, 1.5, 1.3, 3.6, 1.7, 1.5, 3.8, 1.9, 1.7, 1.2, 1.8, 1.1, NA_real_, NA_real_, NA_real_, NA_real_),
    ListOfTSM = c(NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, "R", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
)

#     Species    Tag OriginalStemID TrueStemID CensusID  ExactDate   DBH ListOfTSM
#      <char> <char>          <int>      <int>    <int>     <Date> <num>    <char>
#  1:     sp1 242799         214534         NA        1 1982-06-29   9.5      <NA>
#  2:     sp1 242799         557119         NA        1 1982-06-29   2.0      <NA>
#  3:     sp1 242799         662198         NA        2 1985-01-31  11.0      <NA>
#  4:     sp1 242799         731211         NA        2 1985-01-31   3.0      <NA>
#  5:     sp1 242799         784143         NA        3 1990-03-23  14.5      <NA>
#  6:     sp1 242799         826042         NA        3 1990-03-23   5.2      <NA>
#  7:     sp1 242799         859580         NA        4 1995-01-31  14.8      <NA>
#  8:     sp1 242799         885940         NA        4 1995-01-31   5.4      <NA>
#  9:     sp1 242799         907508         NA        5 2000-02-17   1.9         R
# 10:     sp1 242799         925108         NA        6 2005-02-10    NA         R
# 11:     sp1 242799         925108     925108        7 2010-02-04    NA      <NA>
# 12:     sp1 242799         925108     925108        8 2015-03-16    NA      <NA>
# 13:     sp1 242799         925108     925108        9 2022-03-31    NA      <NA>
# 14:     sp1 242799        1115249    1115249        9 2022-03-31   1.7      <NA>

tag_76 <- data.table(
    Species = "sp1",
    Tag = 242799L,
    OriginalStemID = c(214534, 557119, 662198, 731211, 784143, 826042, 859580, 885940, 907508, 925108, 925108, 925108, 925108, 1115249),
    TrueStemID = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 925108, 925108, 925108, 1115249),
    CensusID = c(1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L, 5L, 6L, 7L, 8L, 9L, 9L),
    ExactDate = as.Date(c(
        "1982-06-29", "1982-06-29", "1985-01-31", "1985-01-31", "1990-03-23", "1990-03-23", "1995-01-31", "1995-01-31", "2000-02-17", "2005-02-10", "2010-02-04", "2015-03-16", "2022-03-31", "2022-03-31"
    )),
    DBH = c(9.5, 2.0, 11.0, 3.0, 14.5, 5.2, 14.8, 5.4, 1.9, NA_real_, NA_real_, NA_real_, NA_real_, 1.7),
    ListOfTSM = c(NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, "R", "R", NA_character_, NA_character_, NA_character_, NA_character_)
)

#    Species    Tag OriginalStemID TrueStemID CensusID  ExactDate   DBH ListOfTSM
#     <char> <char>          <int>      <int>    <int>     <Date> <num>    <char>
# 1:     sp1 159367         138454         NA        1 1982-03-23   1.5      <NA>
# 2:     sp1 159367         138454         NA        2 1985-04-20   1.0      <NA>
# 3:     sp1 159367         138454         NA        3 1990-09-20   1.6      <NA>
# 4:     sp1 159367         138454         NA        4 1995-05-30   1.9      <NA>
# 5:     sp1 159367         138454         NA        5 2000-05-17   2.0      <NA>
# 6:     sp1 159367         138454         NA        6 2005-06-08   2.0      <NA>
# 7:     sp1 159367         138454     138454        7 2010-06-03   2.3      <NA>
# 8:     sp1 159367         138454     138454        8 2015-06-30   2.3      <NA>
# 9:     sp1 159367         138454     138454        9 2022-09-28    NA         R

tag_77 <- data.table(
    Species = "sp1",
    Tag = 159367L,
    OriginalStemID = c(138454, 138454, 138454, 138454, 138454, 138454, 138454, 138454, 138454),
    TrueStemID = c(NA, NA, NA, NA, NA, NA, 138454, 138454, 138454),
    CensusID = c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L),
    ExactDate = as.Date(c(
        "1982-03-23", "1985-04-20", "1990-09-20", "1995-05-30", "2000-05-17", "2005-06-08", "2010-06-03", "2015-06-30", "2022-09-28"
    )),
    DBH = c(1.5, 1.0, 1.6, 1.9, 2.0, 2.0, 2.3, 2.3, NA_real_),
    ListOfTSM = c(NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, "R")
)

# r$> tree_167848
#    Species    Tag OriginalStemID TrueStemID CensusID  ExactDate   DBH ListOfTSM
#     <char> <char>          <int>      <int>    <int>     <Date> <num>    <char>
# 1:     sp1 167848         146104         NA        1 1982-04-08   1.0      <NA>
# 2:     sp1 167848         146104         NA        2 1985-04-27   1.0      <NA>
# 3:     sp1 167848         146104         NA        3 1990-09-08   1.5      <NA>
# 4:     sp1 167848         146104         NA        4 1995-04-12   1.5      <NA>
# 5:     sp1 167848         146104         NA        5 2000-04-13   1.5      <NA>
# 6:     sp1 167848         146104         NA        6 2005-05-17   1.6      <NA>
# 7:     sp1 167848         146104     146104        7 2010-05-18   1.6      <NA>
# 8:     sp1 167848         146104     146104        8 2015-07-15   1.7      <NA>
# 9:     sp1 167848         146104     146104        9 2022-09-27    NA         R

tag_78 <- data.table(
    Species = "sp1",
    Tag = 167848L,
    OriginalStemID = c(146104, 146104, 146104, 146104, 146104, 146104, 146104, 146104, 146104),
    TrueStemID = c(NA, NA, NA, NA, NA, NA, 146104, 146104, 146104),
    CensusID = c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L),
    ExactDate = as.Date(c(
        "1982-04-08", "1985-04-27", "1990-09-08", "1995-04-12", "2000-04-13", "2005-05-17", "2010-05-18", "2015-07-15", "2022-09-27"
    )),
    DBH = c(1.0, 1.0, 1.5, 1.5, 1.5, 1.6, 1.6, 1.7, NA_real_),
    ListOfTSM = c(NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, "R")
)

#    Species    Tag OriginalStemID TrueStemID CensusID  ExactDate   DBH ListOfTSM
#     <char> <char>          <int>      <int>    <int>     <Date> <num>    <char>
# 1:     sp1 170662         147963         NA        1 1982-02-25   1.0      <NA>
# 2:     sp1 170662         147963         NA        2 1985-04-04   1.5      <NA>
# 3:     sp1 170662         147963         NA        3 1990-06-02   2.0      <NA>
# 4:     sp1 170662         147963         NA        4 1995-03-23   2.4      <NA>
# 5:     sp1 170662         147963         NA        5 2000-03-14   2.5      <NA>
# 6:     sp1 170662         147963         NA        6 2005-03-22   2.5      <NA>
# 7:     sp1 170662         147963     147963        7 2010-03-23   2.5      <NA>
# 8:     sp1 170662         147963     147963        8 2015-05-25    NA         R
# 9:     sp1 170662         147963     147963        9 2023-04-19    NA      <NA>

tag_79 <- data.table(
    Species = "sp1",
    Tag = 170662L,
    OriginalStemID = c(147963, 147963, 147963, 147963, 147963, 147963, 147963, 147963, 147963),
    TrueStemID = c(NA, NA, NA, NA, NA, NA, 147963, 147963, 147963),
    CensusID = c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L),
    ExactDate = as.Date(c(
        "1982-02-25", "1985-04-04", "1990-06-02", "1995-03-23", "2000-03-14", "2005-03-22", "2010-03-23", "2015-05-25", "2023-04-19"
    )),
    DBH = c(1.0, 1.5, 2.0, 2.4, 2.5, 2.5, 2.5, NA_real_, NA_real_),
    ListOfTSM = c(NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, "R", NA_character_)
)

#    Species    Tag OriginalStemID TrueStemID CensusID  ExactDate   DBH ListOfTSM
#     <char> <char>          <int>      <int>    <int>     <Date> <num>    <char>
# 1:     sp1 106452          91320         NA        1 1981-09-30   1.0      <NA>
# 2:     sp1 106452          91320         NA        2 1985-05-21   1.5      <NA>
# 3:     sp1 106452          91320         NA        3 1990-08-22   1.8      <NA>
# 4:     sp1 106452          91320         NA        4 1995-05-30    NA      <NA>

tag_80 <- data.table(
    Species = "sp1",
    Tag = 106452L,
    OriginalStemID = c(91320, 91320, 91320, 91320),
    TrueStemID = c(NA, NA, NA, NA),
    CensusID = c(1L, 2L, 3L, 4L),
    ExactDate = as.Date(c(
        "1981-09-30", "1985-05-21", "1990-08-22", "1995-05-30"
    )),
    DBH = c(1.0, 1.5, 1.8, NA_real_),
    ListOfTSM = c(NA_character_, NA_character_, NA_character_, NA_character_)
)

#    Species    Tag OriginalStemID TrueStemID CensusID  ExactDate   DBH ListOfTSM
#     <char> <char>          <int>      <int>    <int>     <Date> <num>    <char>
# 1:     sp1 600156         325783         NA        5 2000-02-25   1.2      <NA>
# 2:     sp1 600156         325783         NA        6 2005-02-22    NA      <NA>
# 3:     sp1 600156         325783     325783        7 2010-02-26   1.3      <NA>
# 4:     sp1 600156         325783     325783        8 2015-03-24   1.4      <NA>
# 5:     sp1 600156         325783     325783        9 2022-03-31    NA         R

tag_81 <- data.table(
    Species = "sp1",
    Tag = 600156L,
    OriginalStemID = c(325783, 325783, 325783, 325783, 325783),
    TrueStemID = c(NA, NA, 325783, 325783, 325783),
    CensusID = c(5L, 6L, 7L, 8L, 9L),
    ExactDate = as.Date(c(
        "2000-02-25", "2005-02-22", "2010-02-26", "2015-03-24", "2022-03-31"
    )),
    DBH = c(1.2, NA_real_, 1.3, 1.4, NA_real_),
    ListOfTSM = c(NA_character_, NA_character_, NA_character_, NA_character_, "R")
)

#    Species    Tag OriginalStemID TrueStemID CensusID  ExactDate   DBH ListOfTSM
#     <char> <char>          <int>      <int>    <int>     <Date> <num>    <char>
# 1:     sp1 527220         318026         NA        4 1995-06-14   1.7      <NA>
# 2:     sp1 527220         318026         NA        5 2000-06-01    NA      <NA>
# 3:     sp1 527220         318026         NA        6 2005-06-30   3.2      <NA>
# 4:     sp1 527220         318026     318026        7 2010-05-25   4.6      <NA>
# 5:     sp1 527220         318026     318026        8 2015-09-01   6.0      <NA>
# 6:     sp1 527220         318026     318026        9 2022-12-21   8.0      <NA>

tag_82 <- data.table(
    Species = "sp1",
    Tag = 527220L,
    OriginalStemID = c(318026, 318026, 318026, 318026, 318026, 318026),
    TrueStemID = c(NA, NA, NA, 318026, 318026, 318026),
    CensusID = c(4L, 5L, 6L, 7L, 8L, 9L),
    ExactDate = as.Date(c(
        "1995-06-14", "2000-06-01", "2005-06-30", "2010-05-25", "2015-09-01", "2022-12-21"
    )),
    DBH = c(1.7, NA_real_, 3.2, 4.6, 6.0, 8.0),
    ListOfTSM = c(NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
)

tag_83 <- data.table(
    Species = "sp1",
    Tag = 83L,
    OriginalStemID = c(318026, 318026, 318026, 318026, 318026, 318026),
    TrueStemID = c(318026, 318026, 318026, 318026, 318026, 318026),
    CensusID = c(4L, 5L, 6L, 7L, 8L, 9L),
    ExactDate = as.Date(c(
        "1995-06-14", "2000-06-01", "2005-06-30", "2010-05-25", "2015-09-01", "2022-12-21"
    )),
    DBH = c(1.7, 2.3, 3.2, 4.6, 6.0, 8.0),
    ListOfTSM = c(NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
)

tag_84 <- data.table(
    Species = "sp1",
    Tag = 84L,
    OriginalStemID = c(318026, NA, 318026, 318026, 318026, 318026),
    TrueStemID = c(318026, NA, 318026, 318026, 318026, 318026),
    CensusID = c(4L, 5L, 6L, 7L, 8L, 9L),
    ExactDate = as.Date(c(
        "1995-06-14", "2000-06-01", "2005-06-30", "2010-05-25", "2015-09-01", "2022-12-21"
    )),
    DBH = c(1.7, 2.3, 3.2, 4.6, 6.0, 8.0),
    ListOfTSM = c(NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
)

tag_85 <- data.table(
    Species = "sp1",
    Tag = 85L,
    OriginalStemID = c(318026, NA, 318026, NA, 318026, 318026),
    TrueStemID = c(318026, NA, 318026, NA, 318026, 318026),
    CensusID = c(4L, 5L, 6L, 7L, 8L, 9L),
    ExactDate = as.Date(c(
        "1995-06-14", "2000-06-01", "2005-06-30", "2010-05-25", "2015-09-01", "2022-12-21"
    )),
    DBH = c(1.7, 2.3, 3.2, 4.6, 6.0, 8.0),
    ListOfTSM = c(NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
)

#     species    Tag StemID TrueStemID CensusID  ExactDate   DBH ListOfTSM
#      <char> <char>  <int>      <int>    <int>     <IDat> <num>    <char>
#  1:  alsebl 162975 141986     141986        1 1982-03-03   3.0      <NA>
#  2:  alsebl 162975 141986     141986        2 1985-04-04   3.0      <NA>
#  3:  alsebl 162975 250132     250132        2 1985-04-04   1.0      <NA>
#  4:  alsebl 162975 141986     141986        3 1990-08-24   3.8      <NA>
#  5:  alsebl 162975 250132     250132        3 1990-08-24   1.1      <NA>
#  6:  alsebl 162975 141986     141986        4 1995-06-16   3.5      <NA>
#  7:  alsebl 162975 250132     250132        4 1995-06-16   1.1      <NA>
#  8:  alsebl 162975 141986     141986        5 2000-04-05   3.7      <NA>
#  9:  alsebl 162975 250132     250132        5 2000-04-05   1.1      <NA>
# 10:  alsebl 162975 141986     141986        6 2005-05-04   3.8      <NA>
# 11:  alsebl 162975 250132     250132        6 2005-05-04   1.2         L
# 12:  alsebl 162975 141986     141986        7 2010-05-06   4.6      <NA>
# 13:  alsebl 162975 250132     250132        7 2010-05-06   1.2      <NA>
# 14:  alsebl 162975 141986     141986        8 2015-06-02   4.6         M
# 15:  alsebl 162975 250132     250132        8 2015-06-02   1.3      <NA>
# 16:  alsebl 162975 141986     141986        9 2022-09-27   4.6      <NA>
# 17:  alsebl 162975 250132     250132        9 2022-09-27   1.5      <NA>

tag_86 <- data.table(
    Species = "sp1",
    Tag = 86L,
    OriginalStemID = c(141986, 141986, 250132, 141986, 250132, 141986, 250132, 141986, 250132, 141986, 250132, 141986, 250132, 141986, 250132, 141986, 250132),
    TrueStemID = c(141986, 141986, 250132, 141986, 250132, 141986, 250132, 141986, 250132, 141986, 250132, 141986, 250132, 141986, 250132, 141986, 250132), 
    CensusID = c(1L, 2L, 2L, 3L, 3L, 4L, 4L, 5L, 5L, 6L, 6L, 7L, 7L, 8L, 8L, 9L, 9L),
    ExactDate = as.Date(c(
        "1982-03-03", "1985-04-04", "1985-04-04", "1990-08-24", "1990-08-24", "1995-06-16", "1995-06-16", "2000-04-05", "2000-04-05", "2005-05-04", "2005-05-04", "2010-05-06", "2010-05-06", "2015-06-02", "2015-06-02", "2022-09-27", "2022-09-27"
    )),
    DBH = c(3.0, 3.0, 1.0, 3.8, 1.1, 3.5, 1.1, 3.7, 1.1, 3.8, 1.2, 4.6, 1.2, 4.6, 1.3, 4.6, 1.5),
    ListOfTSM = c(NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, "L", NA_character_, NA_character_, "M", NA_character_, NA_character_, NA_character_)
)

tag_87 <- data.table(
    Species = "sp1",
    Tag = 87L,
    OriginalStemID = c(141986, 141986, 250132, 141986, NA, 141986, 250132, 141986, 250132, 141986, 250132, 141986, 250132, 141986, 250132, 141986, 250132),
    TrueStemID = c(141986, 141986, 250132, 141986, NA, 141986, 250132, 141986, 250132, 141986, 250132, 141986, 250132, 141986, 250132, 141986, 250132),
    CensusID = c(1L, 2L, 2L, 3L, 3L, 4L, 4L, 5L, 5L, 6L, 6L, 7L, 7L, 8L, 8L, 9L, 9L),
    ExactDate = as.Date(c(
        "1982-03-03", "1985-04-04", "1985-04-04", "1990-08-24", "1990-08-24", "1995-06-16", "1995-06-16", "2000-04-05", "2000-04-05", "2005-05-04", "2005-05-04", "2010-05-06", "2010-05-06", "2015-06-02", "2015-06-02", "2022-09-27", "2022-09-27"
    )),
    DBH = c(3.0, 3.0, 1.0, 3.8, 1.1, 3.5, 1.1, 3.7, 1.1, 3.8, 1.2, 4.6, 1.2, 4.6, 1.3, 4.6, 1.5),
    ListOfTSM = c(NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, "L", NA_character_, NA_character_, "M", NA_character_, NA_character_, NA_character_)
)

tag_88 <- data.table(
    Species = "sp1",
    Tag = 88L,
    OriginalStemID = c(141986, NA, 250132, 141986, NA, 141986, 250132, NA, 250132, 141986, 250132, 141986, 250132, NA, 250132, 141986, 250132),
    TrueStemID = c(141986, NA, 250132, 141986, NA, 141986, 250132, NA, 250132, 141986, 250132, 141986, 250132, NA, 250132, 141986, 250132),
    CensusID = c(1L, 2L, 2L, 3L, 3L, 4L, 4L, 5L, 5L, 6L, 6L, 7L, 7L, 8L, 8L, 9L, 9L),
    ExactDate = as.Date(c(
        "1982-03-03", "1985-04-04", "1985-04-04", "1990-08-24", "1990-08-24", "1995-06-16", "1995-06-16", "2000-04-05", "2000-04-05", "2005-05-04", "2005-05-04", "2010-05-06", "2010-05-06", "2015-06-02", "2015-06-02", "2022-09-27", "2022-09-27"
    )),
    DBH = c(3.0, 3.0, 1.0, 3.8, 1.1, 3.5, 1.1, 3.7, 1.1, 3.8, 1.2, 4.6, 1.2, 4.6, 1.3, 4.6, 1.5),
    ListOfTSM = c(NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, "L", NA_character_, NA_character_, "M", NA_character_, NA_character_, NA_character_)
)


dt_complete_extra <- dt_complete_extra[order(Tag, CensusID)]
dt_complete_extra <- rbindlist(list(
    dt_complete_extra,
    tag_67, 
    tag_68,
    tag_69,
    tag_70,
    tag_71,
    tag_72,
    tag_73,
    tag_74,
    tag_75,
    tag_76, 
    tag_77,
    tag_78,
    tag_79,
    tag_80, 
    tag_81,
    tag_82, 
    tag_83,
    tag_84,
    tag_85,
    tag_86,
    tag_87,
    tag_88
), use.names = TRUE, fill = TRUE)

############################################################
### M-CODE TEST TAGS
### Three tags designed to test the M-coded main-stem constraint.
### OriginalStemID is provided for post-hoc validation ONLY —
### the DP algorithm must not use it.
############################################################

# tag_M1 (tag 901):
# 1 stem at C1-C2, branches to 2 stems at C3.
# M on stem_M (5.2 cm) — nearly same size as stem_X (5.1 cm).
# Without M: both assignments are growth-feasible; ambiguous.
# With M: stem_M must trace back to the single C1/C2 stem.
# TrueStemID at anchor C7: stem_M = 9011, stem_X = 9012.
#
#    Tag CensusID  ExactDate  DBH OriginalStemID TrueStemID ListOfTSM
#  1: 901        1 1982-01-01  5.0          9011         NA      <NA>
#  2: 901        2 1987-01-01  5.1          9011         NA      <NA>
#  3: 901        3 1992-01-01  5.2          9011         NA         M   <- branching + M
#  4: 901        3 1992-01-01  5.1          9012         NA      <NA>   <- new branch (no M)
#  5: 901        4 1997-01-01  5.3          9011         NA      <NA>
#  6: 901        4 1997-01-01  5.2          9012         NA      <NA>
#  7: 901        5 2002-01-01  5.5          9011         NA      <NA>
#  8: 901        5 2002-01-01  5.4          9012         NA      <NA>
#  9: 901        6 2007-01-01  5.7          9011         NA      <NA>
# 10: 901        6 2007-01-01  5.6          9012         NA      <NA>
# 11: 901        7 2012-01-01  5.9          9011       9011      <NA>   <- anchor
# 12: 901        7 2012-01-01  5.8          9012       9012      <NA>   <- anchor
tag_M1 <- data.table(
    Species        = "sp1",
    Tag            = 901L,
    OriginalStemID = c(9011L, 9011L, 9011L, 9012L, 9011L, 9012L, 9011L, 9012L, 9011L, 9012L, 9011L, 9012L),
    TrueStemID     = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 9011L, 9012L),
    CensusID       = c(1L, 2L, 3L, 3L, 4L, 4L, 5L, 5L, 6L, 6L, 7L, 7L),
    ExactDate      = as.Date(c(
        "1982-01-01", "1987-01-01",
        "1992-01-01", "1992-01-01",
        "1997-01-01", "1997-01-01",
        "2002-01-01", "2002-01-01",
        "2007-01-01", "2007-01-01",
        "2012-01-01", "2012-01-01"
    )),
    DBH            = c(5.0, 5.1, 5.2, 5.1, 5.3, 5.2, 5.5, 5.4, 5.7, 5.6, 5.9, 5.8),
    ListOfTSM      = c(NA_character_, NA_character_,
                       "M", NA_character_,          # C3 branching: M on bole
                       NA_character_, NA_character_,
                       NA_character_, NA_character_,
                       NA_character_, NA_character_,
                       NA_character_, NA_character_)
)

# tag_M2 (tag 902):
# Same structure as tag_M1 but M appears at EVERY census from C3 onwards
# (legacy M annotation after branching, stable stem count).
# Only C3 is a branching event. C4-C6 stable-count M must NOT constrain.
# Expected output: identical TrueStemID assignment to tag_M1.
tag_M2 <- data.table(
    Species        = "sp1",
    Tag            = 902L,
    OriginalStemID = c(9021L, 9021L, 9021L, 9022L, 9021L, 9022L, 9021L, 9022L, 9021L, 9022L, 9021L, 9022L),
    TrueStemID     = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 9021L, 9022L),
    CensusID       = c(1L, 2L, 3L, 3L, 4L, 4L, 5L, 5L, 6L, 6L, 7L, 7L),
    ExactDate      = as.Date(c(
        "1982-01-01", "1987-01-01",
        "1992-01-01", "1992-01-01",
        "1997-01-01", "1997-01-01",
        "2002-01-01", "2002-01-01",
        "2007-01-01", "2007-01-01",
        "2012-01-01", "2012-01-01"
    )),
    DBH            = c(5.0, 5.1, 5.2, 5.1, 5.3, 5.2, 5.5, 5.4, 5.7, 5.6, 5.9, 5.8),
    ListOfTSM      = c(NA_character_, NA_character_,
                       "M", NA_character_,           # C3 branching: M on bole
                       "M", NA_character_,           # C4 stable: legacy M (must not constrain)
                       "M", NA_character_,           # C5 stable: legacy M
                       "M", NA_character_,           # C6 stable: legacy M
                       NA_character_, NA_character_)
)

# tag_M3 (tag 903):
# M is on the SMALLER stem at branching (reverse of "pick largest" heuristic).
# stem_M = 2.0 cm (M-coded bole), stem_X = 8.0 cm (no M, larger stem).
# Without M: DBH-based cost would likely mis-assign the large stem as the
# continuation of the C1/C2 single stem (5 cm). With M: stem_M is pinned as
# the continuing bole despite being smaller.
# TrueStemID at anchor C7: stem_M = 9031, stem_X = 9032.
#
#    Tag CensusID  ExactDate  DBH OriginalStemID TrueStemID ListOfTSM
#  1: 903        1 1982-01-01  5.0          9031         NA      <NA>
#  2: 903        2 1987-01-01  4.8          9031         NA      <NA>
#  3: 903        3 1992-01-01  2.0          9031         NA         M   <- branching + M (smaller bole)
#  4: 903        3 1992-01-01  8.0          9032         NA      <NA>   <- new branch (larger, no M)
#  5: 903        4 1997-01-01  2.2          9031         NA      <NA>
#  6: 903        4 1997-01-01  8.5          9032         NA      <NA>
#  7: 903        5 2002-01-01  2.5          9031         NA      <NA>
#  8: 903        5 2002-01-01  9.0          9032         NA      <NA>
#  9: 903        6 2007-01-01  2.8          9031         NA      <NA>
# 10: 903        6 2007-01-01  9.5          9032         NA      <NA>
# 11: 903        7 2012-01-01  3.1          9031       9031      <NA>   <- anchor
# 12: 903        7 2012-01-01 10.0          9032       9032      <NA>   <- anchor
tag_M3 <- data.table(
    Species        = "sp1",
    Tag            = 903L,
    OriginalStemID = c(9031L, 9031L, 9031L, 9032L, 9031L, 9032L, 9031L, 9032L, 9031L, 9032L, 9031L, 9032L),
    TrueStemID     = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 9031L, 9032L),
    CensusID       = c(1L, 2L, 3L, 3L, 4L, 4L, 5L, 5L, 6L, 6L, 7L, 7L),
    ExactDate      = as.Date(c(
        "1982-01-01", "1987-01-01",
        "1992-01-01", "1992-01-01",
        "1997-01-01", "1997-01-01",
        "2002-01-01", "2002-01-01",
        "2007-01-01", "2007-01-01",
        "2012-01-01", "2012-01-01"
    )),
    DBH            = c(5.0, 4.8, 2.0, 8.0, 2.2, 8.5, 2.5, 9.0, 2.8, 9.5, 3.1, 10.0),
    ListOfTSM      = c(NA_character_, NA_character_,
                       "M", NA_character_,          # C3 branching: M on smaller bole
                       NA_character_, NA_character_,
                       NA_character_, NA_character_,
                       NA_character_, NA_character_,
                       NA_character_, NA_character_)
)

dt_complete_extra <- rbindlist(list(
    dt_complete_extra,
    tag_M1,
    tag_M2,
    tag_M3
), use.names = TRUE, fill = TRUE)

################################################################################
### ROW-COUNT INVARIANT EDGE CASES
################################################################################
# These tags test that every input row is preserved in the output
# (no row duplication, no row loss) for every edge-case scenario.

# EC1: post-anchor only, single census C8
tag_EC1 <- data.table(
    Species = "sp1", Tag = 9901L, OriginalStemID = 1L,
    TrueStemID = 1L, CensusID = 8L,
    ExactDate = as.Date("2015-03-10"), DBH = 5.0
)

# EC2: post-anchor only, two censuses (C8-C9)
tag_EC2 <- data.table(
    Species = "sp1", Tag = 9902L, OriginalStemID = 1L,
    TrueStemID = c(1L, 1L), CensusID = c(8L, 9L),
    ExactDate = as.Date(c("2015-03-10", "2022-06-01")), DBH = c(5.0, 6.2)
)

# EC3: post-anchor only, multi-stem (2 stems at C8 and C9)
tag_EC3 <- data.table(
    Species = "sp2", Tag = 9903L,
    OriginalStemID = c(1L, 2L, 1L, 2L),
    TrueStemID = c(1L, 2L, 1L, 2L),
    CensusID = c(8L, 8L, 9L, 9L),
    ExactDate = as.Date(c("2015-03-10", "2015-03-10", "2022-06-01", "2022-06-01")),
    DBH = c(10.0, 3.0, 11.5, 4.0)
)

# EC4: single census at anchor (C7 only, 1 stem)
tag_EC4 <- data.table(
    Species = "sp1", Tag = 9904L, OriginalStemID = 1L,
    TrueStemID = 1L, CensusID = 7L,
    ExactDate = as.Date("2010-05-01"), DBH = 8.0
)

# EC5: all DBH NA but valid CensusIDs across full range
tag_EC5 <- data.table(
    Species = "sp2", Tag = 9905L,
    OriginalStemID = rep(1L, 9), TrueStemID = rep(NA_integer_, 9),
    CensusID = 1:9,
    ExactDate = as.Date(c(
        "1980-01-01", "1984-09-08", "1989-09-18", "1994-09-08",
        "1999-07-09", "2004-07-29", "2009-08-15", "2014-08-02", "2019-07-20"
    )),
    DBH = rep(NA_real_, 9)
)

# EC6: single row, all NA except Species/Tag (like tag 43)
tag_EC6 <- data.table(
    Species = "sp3", Tag = 9906L,
    OriginalStemID = NA_integer_, TrueStemID = NA_integer_,
    CensusID = NA_integer_, ExactDate = as.Date("2020-01-01"),
    DBH = NA_real_
)

# EC7: anchor + post-anchor only (C7-C9, 1 stem)
tag_EC7 <- data.table(
    Species = "sp1", Tag = 9907L, OriginalStemID = 1L,
    TrueStemID = c(1L, 1L, 1L), CensusID = c(7L, 8L, 9L),
    ExactDate = as.Date(c("2010-05-01", "2015-06-10", "2022-03-15")),
    DBH = c(12.0, 13.5, 15.0)
)

# EC8: multi-stem at anchor only (2 stems, C7 only)
tag_EC8 <- data.table(
    Species = "sp2", Tag = 9908L,
    OriginalStemID = c(1L, 2L), TrueStemID = c(1L, 2L),
    CensusID = c(7L, 7L),
    ExactDate = as.Date(c("2010-05-01", "2010-05-01")),
    DBH = c(20.0, 5.0)
)

# EC9: pre-anchor single census (C3 only, 1 stem, no TrueStemID)
tag_EC9 <- data.table(
    Species = "sp1", Tag = 9909L, OriginalStemID = 1L,
    TrueStemID = NA_integer_, CensusID = 3L,
    ExactDate = as.Date("1990-06-15"), DBH = 2.5
)

# EC10: pre-anchor + anchor + post-anchor with DBH only at post-anchor
#       (forces anchor extension like Tag 44 but with pre-anchor rows too)
tag_EC10 <- data.table(
    Species = "sp1", Tag = 9910L, OriginalStemID = 1L,
    TrueStemID = c(NA_integer_, NA_integer_, NA_integer_, 1L, 1L),
    CensusID = c(3L, 5L, 7L, 8L, 9L),
    ExactDate = as.Date(c("1990-06-15", "2000-03-10", "2010-05-01", "2015-06-10", "2022-03-15")),
    DBH = c(NA_real_, NA_real_, NA_real_, 7.0, 8.5)
)

# EC11: multi-stem, post-anchor only, with some DBH NA (mixed obs)
tag_EC11 <- data.table(
    Species = "sp2", Tag = 9911L,
    OriginalStemID = c(1L, 2L, 1L, 2L),
    TrueStemID = c(1L, 2L, 1L, 2L),
    CensusID = c(8L, 8L, 9L, 9L),
    ExactDate = as.Date(c("2015-03-10", "2015-03-10", "2022-06-01", "2022-06-01")),
    DBH = c(10.0, NA_real_, 11.5, 4.0)
)

# EC12: two rows at same census with same OriginalStemID but one has R flag
tag_EC12 <- data.table(
    Species = "sp1", Tag = 9912L,
    OriginalStemID = c(1L, 1L, 1L),
    TrueStemID = c(NA_integer_, 1L, NA_integer_),
    CensusID = c(5L, 7L, 8L),
    ExactDate = as.Date(c("2000-03-10", "2010-05-01", "2015-06-10")),
    DBH = c(3.0, 5.0, NA_real_),
    ListOfTSM = c(NA_character_, NA_character_, "R")
)

# EC13: Full span C1-C9, single stem, all DBH valid (happy-path baseline)
tag_EC13 <- data.table(
    Species = "sp1", Tag = 9913L, OriginalStemID = 1L,
    TrueStemID = rep(1L, 9), CensusID = 1:9,
    ExactDate = as.Date(c(
        "1980-01-01", "1984-09-08", "1989-09-18", "1994-09-08",
        "1999-07-09", "2004-07-29", "2009-08-15", "2014-08-02", "2019-07-20"
    )),
    DBH = c(3.0, 4.2, 5.5, 7.0, 8.8, 10.5, 12.1, 13.6, 15.0)
)

# EC14: Pre-anchor only (C1-C5), no anchor, no post-anchor
tag_EC14 <- data.table(
    Species = "sp1", Tag = 9914L, OriginalStemID = 1L,
    TrueStemID = rep(NA_integer_, 5), CensusID = 1:5,
    ExactDate = as.Date(c(
        "1980-01-01", "1984-09-08", "1989-09-18", "1994-09-08", "1999-07-09"
    )),
    DBH = c(2.0, 3.1, 4.5, 5.8, 7.0)
)

# EC15: Pre-anchor + anchor, no post-anchor, multi-stem (2 stems)
tag_EC15 <- data.table(
    Species = "sp2", Tag = 9915L,
    OriginalStemID = c(1L, 1L, 2L, 1L, 2L),
    TrueStemID     = c(NA_integer_, NA_integer_, NA_integer_, 1L, 2L),
    CensusID = c(3L, 5L, 5L, 7L, 7L),
    ExactDate = as.Date(c("1990-06-15", "2000-03-10", "2000-03-10", "2010-05-01", "2010-05-01")),
    DBH = c(4.0, 6.0, 2.0, 8.5, 3.5)
)

# EC16: Anchor has NA DBH, pre-anchor has DBH, post-anchor has DBH
#       (anchor extension forward: DP scope goes to C8 because C7 anchor is dead)
tag_EC16 <- data.table(
    Species = "sp1", Tag = 9916L, OriginalStemID = 1L,
    TrueStemID = c(NA_integer_, NA_integer_, NA_integer_, 1L, 1L),
    CensusID = c(3L, 5L, 7L, 8L, 9L),
    ExactDate = as.Date(c("1990-06-15", "2000-03-10", "2010-05-01", "2015-06-10", "2022-03-15")),
    DBH = c(4.0, 6.5, NA_real_, 9.0, 10.5)
)

# EC17: Two stems, stem 1 dies at C5 (NA after), stem 2 recruited at C5
#       Tests mortality + recruitment interplay
tag_EC17 <- data.table(
    Species = "sp1", Tag = 9917L,
    OriginalStemID = c(1L, 1L, 1L, 2L, 1L, 2L, 2L),
    TrueStemID     = c(NA_integer_, NA_integer_, NA_integer_, NA_integer_, 1L, 2L, 2L),
    CensusID = c(1L, 3L, 5L, 5L, 7L, 7L, 8L),
    ExactDate = as.Date(c("1980-01-01", "1990-06-15", "2000-03-10", "2000-03-10",
                          "2010-05-01", "2010-05-01", "2015-06-10")),
    DBH = c(10.0, 12.0, NA_real_, 1.5, NA_real_, 4.0, 5.2)
)

# EC18: Single census C1 only (earliest possible census)
tag_EC18 <- data.table(
    Species = "sp1", Tag = 9918L, OriginalStemID = 1L,
    TrueStemID = NA_integer_, CensusID = 1L,
    ExactDate = as.Date("1980-01-01"), DBH = 15.0
)

# EC19: Single census C9 only (latest possible, far post-anchor)
tag_EC19 <- data.table(
    Species = "sp2", Tag = 9919L, OriginalStemID = 1L,
    TrueStemID = 1L, CensusID = 9L,
    ExactDate = as.Date("2022-03-15"), DBH = 20.0
)

# EC20: Two censuses far apart: C1 and C9 (maximum gap, spans pre+post anchor)
tag_EC20 <- data.table(
    Species = "sp1", Tag = 9920L, OriginalStemID = 1L,
    TrueStemID = c(NA_integer_, 1L), CensusID = c(1L, 9L),
    ExactDate = as.Date(c("1980-01-01", "2022-03-15")),
    DBH = c(5.0, 50.0)
)

# EC21: All 9 censuses present, all DBH NA except anchor C7
tag_EC21 <- data.table(
    Species = "sp1", Tag = 9921L, OriginalStemID = rep(1L, 9),
    TrueStemID = c(rep(NA_integer_, 6), 1L, rep(NA_integer_, 2)),
    CensusID = 1:9,
    ExactDate = as.Date(c(
        "1980-01-01", "1984-09-08", "1989-09-18", "1994-09-08",
        "1999-07-09", "2004-07-29", "2009-08-15", "2014-08-02", "2019-07-20"
    )),
    DBH = c(NA, NA, NA, NA, NA, NA, 12.0, NA, NA)
)

# EC22: Post-anchor only, 3 stems at C8 (high K, tests state-space sizing)
tag_EC22 <- data.table(
    Species = "sp2", Tag = 9922L,
    OriginalStemID = c(1L, 2L, 3L),
    TrueStemID = c(1L, 2L, 3L),
    CensusID = c(8L, 8L, 8L),
    ExactDate = rep(as.Date("2015-03-10"), 3),
    DBH = c(10.0, 5.0, 2.5)
)

# EC23: Palm species, full span C1-C9 (different growth_form pathway)
tag_EC23 <- data.table(
    Species = "sp3", Tag = 9923L, OriginalStemID = rep(1L, 9),
    TrueStemID = c(rep(NA_integer_, 6), 1L, 1L, 1L),
    CensusID = 1:9,
    ExactDate = as.Date(c(
        "1980-01-01", "1984-09-08", "1989-09-18", "1994-09-08",
        "1999-07-09", "2004-07-29", "2009-08-15", "2014-08-02", "2019-07-20"
    )),
    DBH = c(5.0, 7.0, 9.0, 11.0, 13.0, 15.0, 17.0, 19.0, 21.0)
)

# EC24: Zero growth — identical DBH across 5 censuses (C3-C7)
tag_EC24 <- data.table(
    Species = "sp1", Tag = 9924L, OriginalStemID = rep(1L, 5),
    TrueStemID = c(rep(NA_integer_, 4), 1L),
    CensusID = 3:7,
    ExactDate = as.Date(c("1990-06-15", "1994-09-08", "1999-07-09", "2004-07-29", "2009-08-15")),
    DBH = rep(10.0, 5)
)

# EC25: Apparent shrinkage — DBH decreases between censuses (triggers negative growth checks)
tag_EC25 <- data.table(
    Species = "sp1", Tag = 9925L, OriginalStemID = rep(1L, 4),
    TrueStemID = c(NA_integer_, NA_integer_, NA_integer_, 1L),
    CensusID = c(3L, 5L, 6L, 7L),
    ExactDate = as.Date(c("1990-06-15", "2000-03-10", "2004-07-29", "2009-08-15")),
    DBH = c(15.0, 14.0, 12.5, 13.0)
)

# EC26: Post-anchor only, all DBH NA (anchor extension fails → skipped_no_data)
tag_EC26 <- data.table(
    Species = "sp1", Tag = 9926L, OriginalStemID = c(1L, 1L),
    TrueStemID = c(NA_integer_, NA_integer_),
    CensusID = c(8L, 9L),
    ExactDate = as.Date(c("2015-06-10", "2022-03-15")),
    DBH = c(NA_real_, NA_real_)
)

# EC27: 3 stems at anchor C7, 2 die post-anchor, 1 survives to C9
tag_EC27 <- data.table(
    Species = "sp2", Tag = 9927L,
    OriginalStemID = c(1L, 2L, 3L, 1L, 2L, 3L, 1L),
    TrueStemID     = c(1L, 2L, 3L, 1L, 2L, 3L, 1L),
    CensusID = c(7L, 7L, 7L, 8L, 8L, 8L, 9L),
    ExactDate = as.Date(c("2010-05-01", "2010-05-01", "2010-05-01",
                          "2015-06-10", "2015-06-10", "2015-06-10",
                          "2022-03-15")),
    DBH = c(20.0, 8.0, 3.0, 21.0, NA_real_, NA_real_, 22.5)
)

# EC28: Pre-anchor only, multi-census, multi-stem (2 stems at C1, C3, C5)
tag_EC28 <- data.table(
    Species = "sp2", Tag = 9928L,
    OriginalStemID = c(1L, 2L, 1L, 2L, 1L, 2L),
    TrueStemID     = rep(NA_integer_, 6),
    CensusID = c(1L, 1L, 3L, 3L, 5L, 5L),
    ExactDate = as.Date(c("1980-01-01", "1980-01-01", "1990-06-15", "1990-06-15",
                          "2000-03-10", "2000-03-10")),
    DBH = c(8.0, 2.0, 10.0, 3.5, 12.0, 5.0)
)

# EC29: Anchor C7 has NA DBH + single post-anchor C8 with DBH
#       (like tag 44 but with one more post-anchor row; anchor ext to C8, no further rows)
tag_EC29 <- data.table(
    Species = "sp1", Tag = 9929L, OriginalStemID = c(1L, 1L),
    TrueStemID = c(NA_integer_, 1L),
    CensusID = c(7L, 8L),
    ExactDate = as.Date(c("2010-05-01", "2015-06-10")),
    DBH = c(NA_real_, 6.0)
)

# EC30: Sparse full span — DBH only at C1, C5, C9 (large gaps, pre+post anchor)
tag_EC30 <- data.table(
    Species = "sp1", Tag = 9930L, OriginalStemID = rep(1L, 3),
    TrueStemID = c(NA_integer_, NA_integer_, 1L),
    CensusID = c(1L, 5L, 9L),
    ExactDate = as.Date(c("1980-01-01", "2000-03-10", "2022-03-15")),
    DBH = c(5.0, 10.0, 18.0)
)

# EC31: Very large DBH >150 cm (outlier, tests bounds/overflow)
tag_EC31 <- data.table(
    Species = "sp1", Tag = 9931L, OriginalStemID = rep(1L, 3),
    TrueStemID = c(NA_integer_, NA_integer_, 1L),
    CensusID = c(3L, 5L, 7L),
    ExactDate = as.Date(c("1990-06-15", "2000-03-10", "2010-05-01")),
    DBH = c(140.0, 155.0, 165.0)
)

# EC32: Multiple R-flag rows across censuses
tag_EC32 <- data.table(
    Species = "sp1", Tag = 9932L,
    OriginalStemID = c(1L, 1L, 1L, 1L),
    TrueStemID = c(NA_integer_, NA_integer_, 1L, NA_integer_),
    CensusID = c(3L, 5L, 7L, 8L),
    ExactDate = as.Date(c("1990-06-15", "2000-03-10", "2010-05-01", "2015-06-10")),
    DBH = c(3.0, NA_real_, 5.0, NA_real_),
    ListOfTSM = c(NA_character_, "R", NA_character_, "R")
)

# EC33: Tag with ListOfTSM = "B" (broken stem)
tag_EC33 <- data.table(
    Species = "sp2", Tag = 9933L,
    OriginalStemID = rep(1L, 4),
    TrueStemID = c(NA_integer_, NA_integer_, NA_integer_, 1L),
    CensusID = c(3L, 5L, 6L, 7L),
    ExactDate = as.Date(c("1990-06-15", "2000-03-10", "2004-07-29", "2009-08-15")),
    DBH = c(12.0, 15.0, 8.0, 16.0),
    ListOfTSM = c(NA_character_, NA_character_, "B", NA_character_)
)

# EC34: Two stems, identical DBH at every census (maximally confusable)
tag_EC34 <- data.table(
    Species = "sp1", Tag = 9934L,
    OriginalStemID = c(1L, 2L, 1L, 2L, 1L, 2L),
    TrueStemID     = c(NA_integer_, NA_integer_, NA_integer_, NA_integer_, 1L, 2L),
    CensusID = c(3L, 3L, 5L, 5L, 7L, 7L),
    ExactDate = as.Date(c("1990-06-15", "1990-06-15", "2000-03-10", "2000-03-10",
                          "2010-05-01", "2010-05-01")),
    DBH = c(10.0, 10.0, 12.0, 12.0, 14.0, 14.0)
)

# EC35: Pre-anchor with census gaps: C1, C3, C5, C7 (skipping even censuses)
tag_EC35 <- data.table(
    Species = "sp1", Tag = 9935L, OriginalStemID = rep(1L, 4),
    TrueStemID = c(NA_integer_, NA_integer_, NA_integer_, 1L),
    CensusID = c(1L, 3L, 5L, 7L),
    ExactDate = as.Date(c("1980-01-01", "1990-06-15", "2000-03-10", "2010-05-01")),
    DBH = c(3.0, 5.0, 7.5, 10.0)
)

# EC36: Anchor C7 NA DBH + TWO post-anchor censuses with DBH
#       (anchor extension to C8, C9 remains truly post-anchor — the exact 9902 pattern but
#        with a pre-anchor row added)
tag_EC36 <- data.table(
    Species = "sp1", Tag = 9936L, OriginalStemID = rep(1L, 4),
    TrueStemID = c(NA_integer_, NA_integer_, 1L, 1L),
    CensusID = c(5L, 7L, 8L, 9L),
    ExactDate = as.Date(c("2000-03-10", "2010-05-01", "2015-06-10", "2022-03-15")),
    DBH = c(3.0, NA_real_, 6.0, 7.5)
)

# EC37: Multi-stem anchor extension: 2 stems, anchor C7 all NA, two post-anchor censuses
tag_EC37 <- data.table(
    Species = "sp2", Tag = 9937L,
    OriginalStemID = c(1L, 2L, 1L, 2L, 1L, 2L),
    TrueStemID     = c(NA_integer_, NA_integer_, 1L, 2L, 1L, 2L),
    CensusID = c(7L, 7L, 8L, 8L, 9L, 9L),
    ExactDate = as.Date(c("2010-05-01", "2010-05-01", "2015-06-10", "2015-06-10",
                          "2022-03-15", "2022-03-15")),
    DBH = c(NA_real_, NA_real_, 10.0, 3.0, 11.5, 4.5)
)

# EC38: 3 rows at same census (C7), different OriginalStemIDs — tests high K at anchor
tag_EC38 <- data.table(
    Species = "sp2", Tag = 9938L,
    OriginalStemID = c(1L, 2L, 3L, 1L, 2L, 3L),
    TrueStemID     = c(NA_integer_, NA_integer_, NA_integer_, 1L, 2L, 3L),
    CensusID = c(5L, 5L, 5L, 7L, 7L, 7L),
    ExactDate = as.Date(c("2000-03-10", "2000-03-10", "2000-03-10",
                          "2010-05-01", "2010-05-01", "2010-05-01")),
    DBH = c(5.0, 8.0, 12.0, 7.0, 10.0, 14.0)
)

# EC39: Pre-anchor only, all DBH NA (should be skipped_no_data)
tag_EC39 <- data.table(
    Species = "sp1", Tag = 9939L, OriginalStemID = c(1L, 1L, 1L),
    TrueStemID = rep(NA_integer_, 3),
    CensusID = c(1L, 3L, 5L),
    ExactDate = as.Date(c("1980-01-01", "1990-06-15", "2000-03-10")),
    DBH = c(NA_real_, NA_real_, NA_real_)
)

# EC40: Only anchor census, DBH NA (single row, dead anchor → skip)
tag_EC40 <- data.table(
    Species = "sp1", Tag = 9940L, OriginalStemID = 1L,
    TrueStemID = NA_integer_, CensusID = 7L,
    ExactDate = as.Date("2010-05-01"), DBH = NA_real_
)

# EC41: Tag with C6 and C7 only (two consecutive censuses ending at anchor)
tag_EC41 <- data.table(
    Species = "sp1", Tag = 9941L, OriginalStemID = c(1L, 1L),
    TrueStemID = c(NA_integer_, 1L),
    CensusID = c(6L, 7L),
    ExactDate = as.Date(c("2004-07-29", "2009-08-15")),
    DBH = c(8.0, 9.5)
)

# EC42: Tag with dense post-anchor: C7, C8, C9 with 2 stems
#       (standard DP pre-anchor + multi-stem post-anchor reinsertion)
tag_EC42 <- data.table(
    Species = "sp2", Tag = 9942L,
    OriginalStemID = c(1L, 1L, 2L, 1L, 2L, 1L, 2L),
    TrueStemID     = c(NA_integer_, 1L, 2L, 1L, 2L, 1L, 2L),
    CensusID = c(5L, 7L, 7L, 8L, 8L, 9L, 9L),
    ExactDate = as.Date(c("2000-03-10", "2010-05-01", "2010-05-01",
                          "2015-06-10", "2015-06-10", "2022-03-15", "2022-03-15")),
    DBH = c(3.0, 8.0, 2.0, 9.0, 3.5, 10.0, 4.0)
)

# EC43: Anchor extension where first post-anchor census has NO TrueStemID but
#       second post-anchor census does (tests preference for TrueStemID in extension)
tag_EC43 <- data.table(
    Species = "sp1", Tag = 9943L, OriginalStemID = rep(1L, 3),
    TrueStemID = c(NA_integer_, NA_integer_, 1L),
    CensusID = c(7L, 8L, 9L),
    ExactDate = as.Date(c("2010-05-01", "2015-06-10", "2022-03-15")),
    DBH = c(NA_real_, 6.0, 7.5)
)

# EC44: Anchor extension to C9 (skip C8 which also has NA) — deepest possible extension
tag_EC44 <- data.table(
    Species = "sp1", Tag = 9944L, OriginalStemID = rep(1L, 4),
    TrueStemID = c(NA_integer_, NA_integer_, NA_integer_, 1L),
    CensusID = c(5L, 7L, 8L, 9L),
    ExactDate = as.Date(c("2000-03-10", "2010-05-01", "2015-06-10", "2022-03-15")),
    DBH = c(4.0, NA_real_, NA_real_, 8.0)
)

# EC45: Anchor + 2 post-anchor censuses, 1 stem grows then a 2nd stem recruits at C9
tag_EC45 <- data.table(
    Species = "sp1", Tag = 9945L,
    OriginalStemID = c(1L, 1L, 1L, 2L),
    TrueStemID     = c(1L, 1L, 1L, 2L),
    CensusID = c(7L, 8L, 9L, 9L),
    ExactDate = as.Date(c("2010-05-01", "2015-06-10", "2022-03-15", "2022-03-15")),
    DBH = c(10.0, 11.0, 12.0, 3.0)
)

# Add another example

#     Species  Tag  OriginalStemID.  TrueStemID CensusID   DBH ListOfTSM       Status.   ExactDate
#     <fctr> <fctr>     <fctr>    .      <int>    <int> <num>    <char>       <char> .     <IDat>
#  1: sp1    000378        391    .         NA        1   490      <NA>        alive . 1981-04-22
#  2: sp1    000378        391    .         NA        2   450      <NA>        alive . 1985-08-24
#  3: sp1    000378     491025    .         NA        3   491  M;B;cylY        alive . 1991-02-16
#  4: sp1    000378     619374    .         NA        3    40      <NA>        alive . 1991-02-16
#  5: sp1    000378     695366    .         NA        4   546  M;B;cylY        alive . 1995-09-27
#  6: sp1    000378     754928    .         NA        4    34      <NA>        alive . 1995-09-27
#  7: sp1    000378     801945    .         NA        4    11      <NA>        alive . 1995-09-27
#  8: sp1    000378     839366    .         NA        5   571    B;cylY        alive . 2000-09-27
#  9: sp1    000378     869100    .         NA        6   595  M;B;cylN        alive . 2005-10-17
# 10: sp1    000378     893101    .         NA        6   190      <NA>        alive . 2005-10-17
# 11: sp1    000378     893101    .      893101       7    NA       D;C         dead . 2010-09-14
# 12: sp1    000378     893101    .      893101       8    NA         R broken below . 2015-11-26
# 13: sp1    000378    1077577    .     1077577       8   271      <NA>        alive . 2015-11-26
# 14: sp1    000378     893101    .      893101       9    NA        OR broken below . 2023-05-30
# 15: sp1    000378    1077577    .     1077577       9   383         A        alive . 2023-05-30

tag_EC46 <- data.table(
    Species = "sp1", Tag = 000378L,
    OriginalStemID = c(391L, 391L, 491025L, 619374L, 695366L, 754928L, 801945L, 839366L, 869100L, 893101L, 893101L, 893101L, 1077577L, 893101L, 1077577L),
    TrueStemID = c(NA_integer_, NA_integer_, NA_integer_, NA_integer_, NA_integer_, NA_integer_, NA_integer_, NA_integer_, NA_integer_, NA_integer_, 893101L, 893101L, 1077577L, 893101L, 1077577L),
    CensusID = c(1L, 2L, 3L, 3L, 4L, 4L, 4L, 5L, 6L, 6L, 7L, 8L, 8L, 9L, 9L),
    ExactDate = as.Date(c(
        "1981-04-22", "1985-08-24", "1991-02-16", "1991-02-16",
        "1995-09-27", "1995-09-27", "1995-09-27", "2000-09-27", "2005-10-17", "2005-10-17", "2010-09-14", "2015-11-26", "2015-11-26", "2023-05-30", "2023-05-30"
    )),
    DBH = c(490, 450, 491, 40, 546, 34, 11, 571, 595, 190, NA, NA, 271, NA, 383),
    ListOfTSM = c(NA_character_, NA_character_, "M;B;cylY", NA_character_, "M;B;cylY", NA_character_, NA_character_, "B;cylY", "M;B;cylN", NA_character_, "D;C", "R", NA_character_, "OR", "A"),
    Status = c("alive", "alive", "alive", "alive", "alive", "alive", "alive", "alive", "alive", "alive", "dead", "broken below", "alive", "broken below", "alive")
)

#     Species    Tag  OriginalStemID TrueStemID CensusID   DBH    ListOfTSM       Status  ExactDate
#     <fctr> <fctr>     <fctr>    .      <int>    <int> <num>    <char>       <char> .     <IDat>
#  1: sp1     060145      48596       <NA>        1         25      <NA>        alive 1981-06-13
#  2: sp1     060145      48596       <NA>        2         35      <NA>        alive 1985-06-14
#  3: sp1     060145      48596       <NA>        3         35      <NA>        alive 1990-10-20
#  4: sp1     060145      48596       <NA>        4         37      <NA>        alive 1995-06-27
#  5: sp1     060145      48596       <NA>        5         37      <NA>        alive 2000-06-20
#  6: sp1     060145      48596       <NA>        6         39      <NA>        alive 2005-06-30
#  7: sp1     060145      48596      48596        7         NA       D;N         dead 2010-06-21
#  8: sp1     060145      48596       48596        8         NA         R broken below 2015-09-15
#  9: sp1     060145    1077157    1077157        8         23      <NA>        alive 2015-09-15
# 10: sp1     060145      48596     48596        9         NA        OR broken below 2023-03-21
# 11: sp1     060145    1077157    1077157        9         32      <NA>        alive 2023-03-21

tag_EC47 <- data.table(
    Species = "sp1", Tag = 060145L,
    OriginalStemID = c(48596L, 48596L, 48596L, 48596L, 48596L, 48596L, 48596L, 48596L, 1077157L, 48596L, 1077157L),
    TrueStemID = c(NA_integer_, NA_integer_, NA_integer_, NA_integer_, NA_integer_, NA_integer_, NA_integer_, NA_integer_, 1077157L, NA_integer_, 1077157L),
    CensusID = c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 8L, 9L, 9L),
    ExactDate = as.Date(c(
        "1981-06-13", "1985-06-14", "1990-10-20", "1995-06-27", "2000-06-20", "2005-06-30", "2010-06-21", "2015-09-15", "2015-09-15", "2023-03-21", "2023-03-21"
    )),
    DBH = c(25, 35, 35, 37, 37, 39, NA, NA, 23, NA, 32),
    ListOfTSM = c(NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, "D;N", "R", NA_character_, "OR", NA_character_),
    Status = c("alive", "alive", "alive", "alive", "alive", "alive", "dead", "broken below", "alive", "broken below", "alive")
)

#    Species    Tag OriginalStemID TrueStemID CensusID   DBH ListOfTSM       Status  ExactDate
#    <fctr>  <char>     <char>     <char>    <num>    <num>    <char>       <char>     <Date>
# 1: sp1     606162     328837       <NA>        5       12      <NA>        alive 2000-02-18
# 2: sp1     606162     328837       <NA>        6       13      <NA>        alive 2005-02-17
# 3: sp1     606162     328837      328837        7       NA       D;N         dead 2010-02-19
# 4: sp1     606162     328837      328837        8       NA         R broken below 2015-03-24
# 5: sp1     606162    1075846    1075846        8       19      <NA>        alive 2015-03-24
# 6: sp1     606162     328837       328837        9       NA         D         dead 2022-12-08
# 7: sp1     606162    1075846      1075846        9       NA        Ns    stem dead 2022-12-08

tag_EC48 <- data.table(
    Species = "sp1", Tag = 606162L,
    OriginalStemID = c(328837L, 328837L, 328837L, 328837L, 1075846L, 328837L, 1075846L),
    TrueStemID = c(NA_integer_, NA_integer_, NA_integer_, NA_integer_, 1075846L, NA_integer_, NA_integer_),
    CensusID = c(5L, 6L, 7L, 8L, 8L, 9L, 9L),
    ExactDate = as.Date(c(
        "2000-02-18", "2005-02-17", "2010-02-19", "2015-03-24", "2015-03-24", "2022-12-08", "2022-12-08"
    )),
    DBH = c(12, 13, NA, NA, 19, NA, NA),
    ListOfTSM = c(NA_character_, NA_character_, "D;N", "R", NA_character_, "D", "Ns"),
    Status = c("alive", "alive", "dead", "broken below", "alive", "dead", "stem dead")
)

dt_complete_extra <- rbindlist(list(
    dt_complete_extra,
    tag_EC1, tag_EC2, tag_EC3, tag_EC4, tag_EC5, tag_EC6,
    tag_EC7, tag_EC8, tag_EC9, tag_EC10, tag_EC11, tag_EC12,
    tag_EC13, tag_EC14, tag_EC15, tag_EC16, tag_EC17, tag_EC18,
    tag_EC19, tag_EC20, tag_EC21, tag_EC22, tag_EC23, tag_EC24,
    tag_EC25, tag_EC26, tag_EC27, tag_EC28, tag_EC29, tag_EC30,
    tag_EC31, tag_EC32, tag_EC33, tag_EC34, tag_EC35, tag_EC36,
    tag_EC37, tag_EC38, tag_EC39, tag_EC40, tag_EC41, tag_EC42,
    tag_EC43, tag_EC44, tag_EC45,
    tag_EC46
), use.names = TRUE, fill = TRUE)

# Include additional information about growth forms, trees, vs figs
growth_forms <- data.table(
    Species = c("sp1", "sp2", "sp3"),
    growth_form = c("tree", "tree", "palm")
)

dt_complete_extra <- merge(dt_complete_extra, growth_forms, by = "Species", all.x = TRUE)
dt_complete_extra <- dt_complete_extra[order(Tag, CensusID)]

fwrite(dt_complete_extra, here("data_simulation", "data", "simulated_data_1.csv"))
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
    species_list <- unique(dt_complete_extra$Species)

    # Create one plot per species (all tags)
    for (species_name in species_list) {
        species_data <- dt_complete_extra[Species == species_name]

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
    for (tag_id in unique(dt_complete_extra$Tag)) {
        tag_data <- dt_complete_extra[Tag == tag_id]
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
