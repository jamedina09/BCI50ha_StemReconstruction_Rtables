############################################################
### Realism calibration helpers (DP_GLOBAL)
############################################################
#
# Goal
# - Turn DP reconstruction output into diagnostic summaries of:
#   - growth increments (DBH->DBH)
#   - shrinkage frequency and magnitude
#   - mortality frequency (DBH->NA)
#   - recruitment frequency and recruited sizes (NA->DBH)
# - Provide concrete "which parameter to tweak" suggestions.
#
# Usage
#   source("dp_global/R/dp_global_bio.R")
#   source("dp_global/R/sensitivity_transition_cost_bio.R")
#   source("dp_global/R/realism_calibration.R")
#
#   bio <- estimate_bio_pars(xraw, interval_years = 5)
#   base <- bio_pars_to_transition_args(bio)
#   rep <- realism_report_from_reconstruction(out, interval_years = 5, base_args = base)
#   rep$summary
#   rep$suggestions
#

if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required. Install it with install.packages('data.table')", call. = FALSE)
}

extract_track_timeseries <- function(out,
                                    interval_years,
                                    id_col = "ReconstructedStemID",
                                    census_col = "CensusID",
                                    dbh_col = "DBH",
                                    group_cols = c("Tag", "species")) {
    dt <- data.table::as.data.table(out)

    needed <- unique(c(group_cols, census_col, dbh_col, id_col))
    missing <- setdiff(needed, names(dt))
    if (length(missing) > 0L) {
        stop("Missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
    }

    interval_years <- as.numeric(interval_years)
    if (!is.finite(interval_years) || interval_years <= 0) {
        stop("interval_years must be positive.", call. = FALSE)
    }

    # Only observed stems define tracks.
    dt_obs <- dt[!is.na(get(dbh_col)) & !is.na(get(id_col))]
    if (nrow(dt_obs) == 0L) {
        return(data.table::data.table())
    }

    # Ensure one DBH per (group, track, census). If duplicates exist, keep the first.
    key_cols <- c(group_cols, id_col, census_col)
    data.table::setorderv(dt_obs, cols = key_cols)
    dt_obs <- dt_obs[, .SD[1L], by = key_cols]

    # Expand to a complete time series per (group, track) across its observed range.
    tracks <- dt_obs[, .(
        census_min = min(get(census_col), na.rm = TRUE),
        census_max = max(get(census_col), na.rm = TRUE)
    ), by = c(group_cols, id_col)]

    ts <- tracks[, {
        cc <- seq.int(census_min, census_max)
        data.table::data.table(tmp_census = cc)
    }, by = c(group_cols, id_col)]

    data.table::setnames(ts, "tmp_census", census_col)

    ts <- dt_obs[ts, on = key_cols]

    # Standardize column names
    ts[, `:=`(
        track_id = get(id_col),
        census = get(census_col),
        dbh = as.numeric(get(dbh_col))
    )]

    keep <- c(group_cols, "track_id", "census", "dbh")
    ts[, ..keep]
}

extract_adjacent_transitions <- function(track_ts, interval_years, group_cols = c("Tag", "species")) {
    interval_years <- as.numeric(interval_years)
    dt <- data.table::copy(data.table::as.data.table(track_ts))
    if (nrow(dt) == 0L) return(dt)

    data.table::setorderv(dt, cols = c(group_cols, "track_id", "census"))
    dt[, `:=`(
        d0 = dbh,
        d1 = data.table::shift(dbh, type = "lead"),
        c0 = census,
        c1 = data.table::shift(census, type = "lead")
    ), by = c(group_cols, "track_id")]

    # Keep adjacent censuses only
    dt <- dt[!is.na(c1) & (c1 == c0 + 1L)]

    dt[, g := (d1 - d0) / interval_years]

    dt[, case := data.table::fifelse(is.na(d0) & is.na(d1), "NA->NA",
        data.table::fifelse(is.na(d0) & !is.na(d1), "NA->DBH",
            data.table::fifelse(!is.na(d0) & is.na(d1), "DBH->NA", "DBH->DBH")
        )
    )]

    dt[, .(Tag, species, track_id, c0, c1, d0, d1, g, case)]
}

summarize_realism <- function(transitions, base_args) {
    dt <- data.table::copy(data.table::as.data.table(transitions))
    if (nrow(dt) == 0L) {
        return(list(summary = data.table::data.table(), by_group = data.table::data.table()))
    }

    # Threshold flags based on the bio model
    max_shrink <- base_args$max_shrink
    recruit_max_dbh <- base_args$recruit_max_dbh

    dt[, violates_shrink_hard := (case == "DBH->DBH" & is.finite(g) & is.finite(max_shrink) & g < max_shrink)]
    dt[, is_shrink := (case == "DBH->DBH" & is.finite(d0) & is.finite(d1) & d1 < d0)]
    dt[, violates_recruit_max := (case == "NA->DBH" & is.finite(d1) & is.finite(recruit_max_dbh) & d1 > recruit_max_dbh)]

    # Global summary
    summary <- dt[, .(
        n_steps = .N,
        frac_DBH_to_DBH = mean(case == "DBH->DBH"),
        frac_DBH_to_NA = mean(case == "DBH->NA"),
        frac_NA_to_DBH = mean(case == "NA->DBH"),
        frac_NA_to_NA = mean(case == "NA->NA"),
        growth_mean = mean(g[case == "DBH->DBH"], na.rm = TRUE),
        growth_sd = stats::sd(g[case == "DBH->DBH"], na.rm = TRUE),
        shrink_frac = mean(is_shrink, na.rm = TRUE),
        shrink_hard_violation_frac = mean(violates_shrink_hard, na.rm = TRUE),
        recruit_dbh_mean = mean(d1[case == "NA->DBH"], na.rm = TRUE),
        recruit_dbh_sd = stats::sd(d1[case == "NA->DBH"], na.rm = TRUE),
        recruit_max_violation_frac = mean(violates_recruit_max, na.rm = TRUE)
    )]

    by_group <- dt[, .(
        n_steps = .N,
        frac_DBH_to_DBH = mean(case == "DBH->DBH"),
        frac_DBH_to_NA = mean(case == "DBH->NA"),
        frac_NA_to_DBH = mean(case == "NA->DBH"),
        shrink_hard_violation_frac = mean(violates_shrink_hard),
        recruit_max_violation_frac = mean(violates_recruit_max)
    ), by = .(Tag, species)]

    list(summary = summary, by_group = by_group, transitions = dt)
}

tuning_suggestions <- function(realism, base_args) {
    # Heuristic suggestions based on the reconstruction-derived transition frequencies.
    # NOTE: This is not a formal optimizer; it provides *actionable starting points*.

    if (is.null(realism$summary) || nrow(realism$summary) == 0L) {
        return(data.table::data.table())
    }

    s <- realism$summary[1L]

    p_recruit <- 1 - exp(-base_args$recruit_lambda * base_args$interval_years)
    # interval_years might not be stored in base_args; be robust
    # fall back to 1 year if missing
    if (!is.finite(p_recruit)) {
        T <- if (!is.null(base_args$interval_years) && is.finite(base_args$interval_years)) base_args$interval_years else 1
        p_recruit <- 1 - exp(-base_args$recruit_lambda * T)
    }

    out <- list()
    add <- function(issue, param, direction, rationale) {
        out[[length(out) + 1L]] <<- data.table::data.table(
            issue = issue,
            parameter = param,
            change = direction,
            rationale = rationale
        )
    }

    # Shrinkage problems
    if (isTRUE(is.finite(s$shrink_hard_violation_frac)) && s$shrink_hard_violation_frac > 0) {
        add(
            issue = "Hard shrink violations (DBH drop too large)",
            param = "max_shrink",
            direction = "decrease (more negative) to allow larger shrink OR increase (closer to 0) to forbid shrink",
            rationale = "If violations are due to measurement error or buttress effects, allow more shrink (more negative). If they reflect ID swaps, tighten the bound (closer to 0) to prevent such matches."
        )
        add(
            issue = "Frequent shrink steps",
            param = "k_shrink",
            direction = "increase",
            rationale = "Increases the soft penalty for any shrinkage (d1<d0), pushing the DP to prefer non-shrinking matches."
        )
    } else if (isTRUE(is.finite(s$shrink_frac)) && s$shrink_frac > 0.25) {
        add(
            issue = "Many DBH decreases (shrink)",
            param = "k_shrink",
            direction = "increase",
            rationale = "Penalizes shrinkage so reconstructed trajectories look more biologically plausible."
        )
    }

    # Recruitment size problems
    if (isTRUE(is.finite(s$recruit_max_violation_frac)) && s$recruit_max_violation_frac > 0) {
        add(
            issue = "Recruits exceeding recruit_max_dbh (hard violations)",
            param = "recruit_max_dbh",
            direction = "increase (if threshold too strict) OR decrease (to forbid large recruits)",
            rationale = "If real recruits can be larger than the current bound, raise it. If the DP is using recruitment to explain large stems that should match existing tracks, lower the bound and/or reduce recruit_lambda."
        )
    }

    if (isTRUE(is.finite(s$frac_NA_to_DBH)) && s$frac_NA_to_DBH > 0.10) {
        add(
            issue = "Too many NA->DBH events (many new recruits)",
            param = "recruit_lambda",
            direction = "decrease",
            rationale = "Reduces p_recruit = 1-exp(-lambda*T), making recruitment less likely so DP prefers persistence (DBH->DBH) over breaking/creating tracks."
        )
        add(
            issue = "Recruits look too large",
            param = "recruit_meanlog / recruit_sdlog",
            direction = "decrease meanlog and/or decrease sdlog",
            rationale = "Shifts/compacts the recruit DBH distribution so NA->DBH transitions favor small recruited DBH."
        )
    }

    # Mortality frequency problems
    if (isTRUE(is.finite(s$frac_DBH_to_NA)) && s$frac_DBH_to_NA > 0.10) {
        add(
            issue = "Too many DBH->NA events (many deaths/disappearances)",
            param = "h0 / beta",
            direction = "decrease h0 and/or decrease beta (less positive)",
            rationale = "Decreases death probability p_death = 1-exp(-h0*exp(beta*d0)*T), making persistence cheaper relative to disappearance."
        )
    }

    # Growth mismatch problems (broad strokes)
    if (isTRUE(is.finite(s$growth_sd)) && isTRUE(is.finite(base_args$sigma0))) {
        # If reconstructed growth is much more variable than modeled, widen sigma.
        if (s$growth_sd > (base_args$sigma0 * 2)) {
            add(
                issue = "Growth increments are more variable than the model",
                param = "sigma0 / sigma1",
                direction = "increase",
                rationale = "Wider growth likelihood reduces over-penalization of plausible variation, which can reduce forced deaths/recruits caused by tight growth constraints."
            )
        }
    }

    if (length(out) == 0L) {
        return(data.table::data.table(
            issue = "No obvious red flags",
            parameter = NA_character_,
            change = NA_character_,
            rationale = "Reconstruction looks broadly consistent with the current parameterization; next step is to check a few tags visually and compare growth/recruit/mortality distributions to field expectations."
        ))
    }

    data.table::rbindlist(out)
}

realism_report_from_reconstruction <- function(out, interval_years, base_args) {
    # Convenience wrapper that returns a compact report.
    # Adds interval_years into base_args for internal comparisons.
    base_args2 <- base_args
    base_args2$interval_years <- as.numeric(interval_years)

    ts <- extract_track_timeseries(out, interval_years = interval_years)
    tr <- extract_adjacent_transitions(ts, interval_years = interval_years)
    rs <- summarize_realism(tr, base_args = base_args2)
    sugg <- tuning_suggestions(rs, base_args = base_args2)

    list(
        summary = rs$summary,
        by_group = rs$by_group,
        transitions = rs$transitions,
        suggestions = sugg
    )
}
