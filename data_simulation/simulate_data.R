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

# Include additional information about growth forms, trees, vs figs
growth_forms <- data.table(
    Species = c("sp1", "sp2", "sp3"),
    GrowthForm = c("tree", "tree", "fig")
)

dt_complete_extra <- merge(dt_complete_extra, growth_forms, by = "Species", all.x = TRUE)

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
