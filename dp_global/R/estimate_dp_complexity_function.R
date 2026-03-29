############################################################
### estimate_dp_complexity.R
### Function to estimate computational complexity for DP stem reconstruction
############################################################

#' Estimate DP Computational Complexity for Stem Reconstruction
#'
#' This function analyzes a dataset of stem observations and estimates the
#' computational complexity of running the DP stem reconstruction algorithm
#' on each tag. The complexity is measured by the number of state transitions
#' that need to be evaluated during the DP algorithm.
#'
#' @param data_path Character string: path to the CSV file containing stem data
#' @param anchor_start Integer: the census ID that serves as the anchor point (default: 7)
#' @param slack_tracks Integer: number of extra tracks to allow simultaneous death+birth (default: 1)
#' @param min_growth Numeric: minimum allowed growth per year (DP model bound, default -Inf)
#' @param max_growth Numeric: maximum allowed growth per year (DP model bound, default Inf)
#' @param prune_min_growth Numeric|NULL: explicit DP prune min growth (overrides min_growth when non-NULL)
#' @param prune_max_growth Numeric|NULL: explicit DP prune max growth (overrides max_growth when non-NULL)
#' @param prune_use_bio_bounds Logical: whether to intersect with Bio bounds (default TRUE)
#' @param bio_max_shrink Numeric: biological shrink cap (from DP data, default -Inf)
#' @param bio_max_growth Numeric: biological growth cap (from DP data, default Inf)
#' @param recruit_max_dbh Numeric: biological recruitment max DBH (default Inf)
#' @param prune_recruit_max_dbh Numeric|NULL: prune recruit max DBH (overrides recruit_max_dbh)
#' @param prune_use_bio_recruit Logical: whether to use bio recruit max DBH (default TRUE)
#' @param interval_years Numeric|NULL: uniform census interval in years (default NULL => 1 or derived from ExactDate)
#' @return A data.table with complexity estimates for each tag, sorted by computational cost
#' @export
estimate_dp_complexity <- function(data_path,
                                   anchor_start = 7L,
                                   slack_tracks = 1L,
                                   min_growth = -Inf,
                                   max_growth = Inf,
                                   prune_min_growth = NULL,
                                   prune_max_growth = NULL,
                                   prune_use_bio_bounds = TRUE,
                                   bio_max_shrink = -Inf,
                                   bio_max_growth = Inf,
                                   recruit_max_dbh = Inf,
                                   prune_recruit_max_dbh = NULL,
                                   prune_use_bio_recruit = TRUE,
                                   interval_years = NULL) {

    library(data.table)

    # Load data
    if (!file.exists(data_path)) {
        stop("Data file not found: ", data_path)
    }
    data <- fread(data_path)

    # Required columns check (case-insensitive)
    required_cols <- c("Tag", "CensusID", "DBH", "TrueStemID")
    col_names <- names(data)
    missing_cols <- character(0)

    for (req in required_cols) {
        if (!(req %in% col_names)) {
            # Try case-insensitive match
            matches <- col_names[tolower(col_names) == tolower(req)]
            if (length(matches) == 1) {
                # Rename column to expected case
                setnames(data, matches, req)
            } else {
                missing_cols <- c(missing_cols, req)
            }
        }
    }

    if (length(missing_cols) > 0) {
        stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
    }

    # Handle species column (can be Species or species)
    species_col <- NULL
    if ("species" %in% col_names) {
        species_col <- "species"
    } else if ("Species" %in% col_names) {
        species_col <- "Species"
    } else {
        # Add a default species column
        data[, species := "unknown"]
        species_col <- "species"
    }

    # Effective growth/recruit bounds for pruning estimate
    eff_min_growth <- if (!is.null(prune_min_growth)) prune_min_growth else min_growth
    eff_max_growth <- if (!is.null(prune_max_growth)) prune_max_growth else max_growth
    if (isTRUE(prune_use_bio_bounds)) {
        eff_min_growth <- max(eff_min_growth, bio_max_shrink)
        eff_max_growth <- min(eff_max_growth, bio_max_growth)
    }

    eff_recruit_max <- if (!is.null(prune_recruit_max_dbh)) {
        if (isTRUE(prune_use_bio_recruit)) {
            min(recruit_max_dbh, prune_recruit_max_dbh)
        } else {
            prune_recruit_max_dbh
        }
    } else {
        recruit_max_dbh
    }

    # Determine interval years default from ExactDate if available, else 1
    if (is.null(interval_years)) {
        if ("ExactDate" %in% names(data)) {
            coords <- unique(data[, .(CensusID, ExactDate)])
            coords <- coords[order(CensusID)]
            coords[, ExactDate := as.Date(ExactDate)]
            if (nrow(coords) >= 2) {
                intv <- diff(as.numeric(coords$ExactDate)) / 365.25
                intv <- intv[is.finite(intv) & intv > 0]
                if (length(intv) > 0) {
                    interval_years <- median(intv, na.rm = TRUE)
                }
            }
        }
        if (is.null(interval_years) || !is.finite(interval_years) || interval_years <= 0) {
            interval_years <- 1.0
        }
    }

    is_transition_feasible <- function(d0, d1, interval_years, eff_min, eff_max, eff_recruit_max) {
        if (is.na(d0) && is.na(d1)) return(TRUE)
        if (!is.na(d0) && is.na(d1)) return(TRUE)
        if (is.na(d0) && !is.na(d1)) {
            if (is.finite(eff_recruit_max)) {
                return(d1 <= eff_recruit_max)
            }
            return(TRUE)
        }
        # both observed
        g <- (d1 - d0) / interval_years
        (g >= eff_min && g <= eff_max)
    }

    # Function to count injective states (P(K, n_obs) = K! / (K - n_obs)!)
    count_injective_states <- function(K, n_obs) {
        if (!is.finite(K) || !is.finite(n_obs)) {
            return(NA_real_)
        }
        K <- as.integer(K)
        n_obs <- as.integer(n_obs)
        if (n_obs < 0L || K < 0L) {
            return(NA_real_)
        }
        if (n_obs == 0L) {
            return(1)
        }
        if (n_obs > K) {
            return(0)
        }
        prod(seq.int(from = K, to = K - n_obs + 1L, by = -1L))
    }

    # Function to estimate computational complexity for a single tag
    estimate_tag_complexity <- function(tag_data, anchor_start, slack_tracks) {
        # Get observed stem counts per census up to anchor
        obs_counts <- vapply(
            seq_len(anchor_start),
            function(cc) nrow(tag_data[CensusID == cc & !is.na(DBH)]),
            integer(1L)
        )

        max_obs <- if (length(obs_counts) > 0L) max(obs_counts) else 0L

        # Determine K (number of tracks)
        births_needed <- if (length(obs_counts) >= 2L) sum(pmax(0L, diff(obs_counts))) else 0L
        K_from_counts <- as.integer(if (length(obs_counts) > 0L) obs_counts[1L] + births_needed else 0L)

        # Get anchor observations
        anchor_obs <- tag_data[CensusID == anchor_start & !is.na(DBH)]
        anchor_ids <- sort(unique(anchor_obs$TrueStemID))
        anchor_ids <- anchor_ids[!is.na(anchor_ids)]

        K_base <- max(length(anchor_ids), max_obs, K_from_counts)

        # Add slack tracks if needed
        K <- K_base
        if (slack_tracks > 0L && K_base == max_obs) {
            K <- K_base + slack_tracks
        }

        # Compute state counts per census
        n_states_by_census <- vapply(obs_counts, function(n_obs) count_injective_states(K, n_obs), numeric(1L))

        # Total states across all censuses
        total_states <- sum(n_states_by_census)

        # Maximum states in any single census
        max_states_per_census <- max(n_states_by_census, na.rm = TRUE)

        # Rough estimate of unpruned computation
        transition_computations <- sum(
            n_states_by_census[1:(length(n_states_by_census)-1)] *
            n_states_by_census[2:length(n_states_by_census)]
        )

        # Estimate pruned transition computations using growth/recruit constraints
        pruned_transition_computations <- 0
        if (length(obs_counts) >= 2L) {
            census_seq <- seq_len(anchor_start)
            for (i in seq_along(census_seq[-length(census_seq)])) {
                c0 <- census_seq[i]
                c1 <- census_seq[i + 1L]
                t0 <- tag_data[CensusID == c0]
                t1 <- tag_data[CensusID == c1]

                # using all DBH observations (including NA) in this pair
                d0 <- t0$DBH
                d1 <- t1$DBH

                if (length(d0) == 0L || length(d1) == 0L) next

                # broadcast and test all pairs in this census transition
                d0_mat <- matrix(d0, nrow = length(d0), ncol = length(d1))
                d1_mat <- matrix(d1, nrow = length(d0), ncol = length(d1), byrow = TRUE)
                feasible <- mapply(is_transition_feasible,
                                   as.vector(d0_mat),
                                   as.vector(d1_mat),
                                   MoreArgs = list(interval_years = interval_years,
                                                   eff_min = eff_min_growth,
                                                   eff_max = eff_max_growth,
                                                   eff_recruit_max = eff_recruit_max))
                pruned_transition_computations <- pruned_transition_computations + sum(feasible, na.rm = TRUE)
            }
        }

        list(
            tag = unique(tag_data$Tag),
            species = unique(tag_data[[species_col]]),
            n_censuses = anchor_start,
            max_obs_any_census = max_obs,
            K_tracks = K,
            total_states = total_states,
            max_states_per_census = max_states_per_census,
            transition_computations = transition_computations,
            transition_computations_pruned = pruned_transition_computations,
            obs_counts = obs_counts,
            n_states_by_census = n_states_by_census
        )
    }

    # Analyze all tags
    tags <- unique(data$Tag)
    results <- list()

    for (tag in tags) {
        tag_data <- data[Tag == tag]
        if (nrow(tag_data) > 0) {
            complexity <- estimate_tag_complexity(tag_data, anchor_start, slack_tracks)
            results[[as.character(tag)]] <- complexity
        }
    }

    # Convert to data.table for easy analysis
    complexity_dt <- rbindlist(lapply(results, function(x) {
        data.table(
            Tag = x$tag,
            Species = x$species,
            MaxObs = x$max_obs_any_census,
            K = x$K_tracks,
            TotalStates = x$total_states,
            MaxStatesPerCensus = x$max_states_per_census,
            TransitionComputations = x$transition_computations,
            TransitionComputationsPruned = x$transition_computations_pruned
        )
    }))

    # Sort by computational complexity (transition computations as proxy)
    complexity_dt <- complexity_dt[order(-TransitionComputations)]

    return(complexity_dt)
}

#' Get detailed complexity analysis for a specific tag
#'
#' @param data_path Character string: path to the CSV file
#' @param tag Integer: the tag number to analyze
#' @param anchor_start Integer: anchor census (default: 7)
#' @param slack_tracks Integer: slack tracks (default: 1)
#' @return List with detailed complexity information
#' @export

get_tag_complexity_details <- function(data_path, tag,
                                       anchor_start = 7L, slack_tracks = 1L) {
    library(data.table)

    data <- fread(data_path)
    tag_data <- data[Tag == tag]

    if (nrow(tag_data) == 0) {
        stop("Tag ", tag, " not found in data")
    }

    # Handle species column
    col_names <- names(data)
    species_col <- NULL
    if ("species" %in% col_names) {
        species_col <- "species"
    } else if ("Species" %in% col_names) {
        species_col <- "Species"
    } else {
        species_col <- NULL
    }

    # Reuse the complexity estimation logic
    obs_counts <- vapply(
        seq_len(anchor_start),
        function(cc) nrow(tag_data[CensusID == cc & !is.na(DBH)]),
        integer(1L)
    )

    max_obs <- max(obs_counts)
    births_needed <- if (length(obs_counts) >= 2L) sum(pmax(0L, diff(obs_counts))) else 0L
    K_from_counts <- as.integer(obs_counts[1L] + births_needed)

    anchor_obs <- tag_data[CensusID == anchor_start & !is.na(DBH)]
    anchor_ids <- sort(unique(anchor_obs$TrueStemID))
    anchor_ids <- anchor_ids[!is.na(anchor_ids)]

    K_base <- max(length(anchor_ids), max_obs, K_from_counts)
    K <- K_base
    if (slack_tracks > 0L && K_base == max_obs) {
        K <- K_base + slack_tracks
    }

    count_injective_states <- function(K, n_obs) {
        if (n_obs == 0L) return(1)
        if (n_obs > K) return(0)
        prod(seq.int(from = K, to = K - n_obs + 1L, by = -1L))
    }

    n_states_by_census <- vapply(obs_counts, function(n_obs) count_injective_states(K, n_obs), numeric(1L))

    transition_computations <- sum(
        n_states_by_census[1:(length(n_states_by_census)-1)] *
        n_states_by_census[2:length(n_states_by_census)]
    )

    list(
        tag = tag,
        species = unique(tag_data[[species_col]]),
        anchor_start = anchor_start,
        slack_tracks = slack_tracks,
        observations_per_census = data.table(CensusID = 1:anchor_start, N_Obs = obs_counts),
        max_observations = max_obs,
        K_tracks = K,
        states_per_census = data.table(CensusID = 1:anchor_start, N_States = n_states_by_census),
        max_states_per_census = max(n_states_by_census),
        total_states = sum(n_states_by_census),
        transition_computations = transition_computations
    )
}

# Example usage:
if (FALSE) {  # Set to TRUE to run example
    # Estimate complexity for all tags
    complexity <- estimate_dp_complexity("../data_simulation/data/simulation_legacy_backup/simulated_data_two_species.csv")
    print(complexity)

    # Get detailed analysis for a specific tag
    details <- get_tag_complexity_details("../data_simulation/data/simulation_legacy_backup/simulated_data_two_species.csv", tag = 11)
    print(details)
}