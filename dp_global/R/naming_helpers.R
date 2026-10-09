# naming_helpers.R
# Helper functions for building canonical output names and encoding numeric
# values for directory-safe names. Sourced by the drivers (main_cpp.R,
# main_cpp_chunk.R, main_cpp_bci.R and the BCI chunk driver).

# Encode numeric values for directory-safe names
# -0.5 -> m0p5, 7.5 -> 7p5; NULL or NA -> "NA"
encode_num <- function(x) {
    if (is.null(x) || is.na(x)) {
        return("NA")
    }
    s <- as.character(x)
    s <- gsub("-", "m", s)
    s <- gsub("\\.", "p", s)
    s
}

# Build a directory-safe output name using the run's key parameters. This
# function takes no arguments: it reads the driver's global variables
# RUN_ALL_TAGS, WHICH_TAG, DP_MODE, USE_MEASUREMENT_ERROR,
# MAX_GROWTH_HARD_SOURCE, MAX_GROWTH_FIXED, MAX_SHRINK_HARD_SOURCE,
# MAX_SHRINK_FIXED, K_GROWTH_SOURCE, K_GROWTH_FIXED, K_SHRINK_SOURCE and
# K_SHRINK_FIXED, and, when they exist, BATCH_TS, CONFIG_NAME,
# DBH_ROUND_CENSUSES_INT, PROB_BIRTH_DEATH, RECRUIT_RATE_UNIT and
# COVERAGE_RULE.
# Returns one string:
#   <timestamp>_<config>_<tags>_<DP mode>_<ME label>_<g>_<s>_<kg>_<ks>_rcpp
# e.g. 20261005_232052_unknown_allT_DP_MB_NME_R12_BD_LT_CU_g5_sm0p5_kg0_ks0_rcpp
#   tags    : allT (all tags) or T<WHICH_TAG>
#   DP mode : NO_DP (none), DP_S (map), DP_M (marginals), DP_MB
#             (marginals+bins), DP_U (other)
#   g / s   : hard growth / shrink bound: g<value> / s<value> when fixed, gD / sD
#             from data, gU / sU otherwise
#   kg / ks : soft growth / shrink penalty, coded the same way
build_out_dir_name <- function() {
    # Timestamp: use BATCH_TS if provided; else fallback to current date+time
    ts <- if (exists("BATCH_TS") && nzchar(BATCH_TS)) BATCH_TS else format(Sys.time(), "%Y%m%d_%H%M%S")

    # Config name (for output directory label)
    config_part <- if (exists("CONFIG_NAME") && !is.null(CONFIG_NAME)) {
        CONFIG_NAME
    } else {
        "unknown"
    }

    # Tag info
    tag_part <- if (isTRUE(RUN_ALL_TAGS)) {
        "allT"
    } else {
        paste0("T", as.character(WHICH_TAG))
    }

    # DP mode label
    dp_part <- switch(DP_MODE,
        "none" = "NO_DP",
        "map" = "DP_S",
        "marginals" = "DP_M",
        "marginals+bins" = "DP_MB",
        "DP_U"
    )

    # Measurement error label (+ "R<censuses>" when DBH rounded down to classes
    # at those censuses is modelled, e.g. NME_R12 for BCI 1982 and 1985)
    me_part <- if (isTRUE(USE_MEASUREMENT_ERROR)) "ME" else "NME"
    if (exists("DBH_ROUND_CENSUSES_INT") && length(DBH_ROUND_CENSUSES_INT) > 0L) {
        me_part <- paste0(me_part, "_R", paste(sort(DBH_ROUND_CENSUSES_INT), collapse = ""))
    }
    # Drivers that set these (BCI): "BD" = birth-death probabilistic matcher,
    # "LT" = recruitment rate per established tree, "CU" = species size
    # coverage by the union rule (e.g. NME_R12_BD_LT_CU)
    if (exists("PROB_BIRTH_DEATH") && isTRUE(PROB_BIRTH_DEATH)) me_part <- paste0(me_part, "_BD")
    if (exists("RECRUIT_RATE_UNIT") && identical(RECRUIT_RATE_UNIT, "tree")) {
        me_part <- paste0(me_part, "_LT")
    }
    if (exists("COVERAGE_RULE") && identical(COVERAGE_RULE, "union")) me_part <- paste0(me_part, "_CU")

    max_growth_hard_ <- switch(MAX_GROWTH_HARD_SOURCE,
        "fixed" = paste0("g", encode_num(MAX_GROWTH_FIXED)),
        "data"  = "gD",
        "gU"
    )

    max_shrink_hard_ <- switch(MAX_SHRINK_HARD_SOURCE,
        "fixed" = paste0("s", encode_num(MAX_SHRINK_FIXED)),
        "data"  = "sD",
        "sU"
    )

    soft_growth_ <- switch(K_GROWTH_SOURCE,
        "fixed" = paste0("kg", encode_num(K_GROWTH_FIXED)),
        "data"  = "kgD",
        "kgU"
    )

    soft_shrink_ <- switch(K_SHRINK_SOURCE,
        "fixed" = paste0("ks", encode_num(K_SHRINK_FIXED)),
        "data"  = "ksD",
        "ksU"
    )

    # Assemble final directory name
    dir_name <- paste(
        ts,
        config_part,
        tag_part,
        paste0(dp_part, "_", me_part),
        max_growth_hard_,
        max_shrink_hard_,
        soft_growth_,
        soft_shrink_,
        "rcpp",
        sep = "_"
    )

    return(dir_name)
}
