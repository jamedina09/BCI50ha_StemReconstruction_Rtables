# naming_helpers.R
# Helper functions for building canonical output names and encoding numeric
# values for directory-safe names used by `main_cpp.R`.

# Encode numeric values for directory-safe names
# -0.5 -> m0p5, 7.5 -> 7p5
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
# function expects the calling environment (main_cpp.R) to define the
# variables referenced (e.g., BATCH_TS, CONFIG_NAME, WHICH_TAG, DP_MODE, etc.).
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
    # "LT" = recruitment rate per established tree (e.g. NME_R12_BD_LT)
    if (exists("PROB_BIRTH_DEATH") && isTRUE(PROB_BIRTH_DEATH)) me_part <- paste0(me_part, "_BD")
    if (exists("RECRUIT_RATE_UNIT") && identical(RECRUIT_RATE_UNIT, "tree")) {
        me_part <- paste0(me_part, "_LT")
    }

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
