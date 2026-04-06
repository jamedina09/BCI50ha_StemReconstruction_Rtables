############################################################
### Sensitivity analysis for transition_cost_tracks_bio()
############################################################
#
# Goal
# - Understand how sensitive the transition cost is to the parameters estimated
#   by estimate_bio_pars().
# - Visualize (i) smooth sensitivity vs (ii) jump discontinuities (hard penalties)
# - Decompose the total cost into interpretable components so you can see which
#   terms dominate.
#
# Usage (typical)
#   source("dp_global/R/dp_global_bio.R")
#   source("dp_global/R/sensitivity_transition_cost_bio.R")
#
#   bio <- estimate_bio_pars(xraw, interval_years = 5)
#   base <- bio_pars_to_transition_args(bio)
#
#   sc <- make_demo_scenarios(base, interval_years = 5)
#   dt <- sweep_transition_cost(sc[["growth_ok"]]$t, sc[["growth_ok"]]$tp1,
#                              interval_years = 5, base_args = base,
#                              param = "sigma0", values = seq(0.05, 0.5, length.out = 100))
#   plot_sweep_components(dt)
#

if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required. Install it with install.packages('data.table')", call. = FALSE)
}

bio_pars_to_transition_args <- function(bio_pars) {
    # Converts output of estimate_bio_pars() into argument list for
    # transition_cost_tracks_bio_components().
    me <- bio_pars$measurement_error
    use_me <- isTRUE(bio_pars$settings$use_measurement_error)
    if (is.null(me)) {
        me <- list(sd1_a = 0.0062, sd1_b = 0.0904, sd2 = 4.64, p_big = 0.05)
    }
    g <- bio_pars$growth
    mu_const <- if (!is.null(g$alpha)) g$alpha else g$mu
    mu_gamma <- if (!is.null(g$gamma)) g$gamma else 0

    # Support both the legacy flat layout and the newer nested layout.
    max_shrink0 <- tryCatch(bio_pars$shrinkage$guardrails$hard$value, error = function(e) NULL)
    if (is.null(max_shrink0)) max_shrink0 <- bio_pars$shrinkage$max_shrink
    k_shrink0 <- tryCatch(bio_pars$shrinkage$penalties$soft$k, error = function(e) NULL)
    if (is.null(k_shrink0)) k_shrink0 <- bio_pars$shrinkage$k_shrink

    max_growth0 <- tryCatch(bio_pars$growth$guardrails$hard$value, error = function(e) NULL)
    if (is.null(max_growth0)) max_growth0 <- bio_pars$growth$max_growth
    max_growth_soft0 <- tryCatch(bio_pars$growth$guardrails$soft$value, error = function(e) NULL)
    if (is.null(max_growth_soft0)) max_growth_soft0 <- bio_pars$growth$max_growth_soft
    k_growth0 <- tryCatch(bio_pars$growth$penalties$soft$k, error = function(e) NULL)
    if (is.null(k_growth0)) k_growth0 <- bio_pars$growth$k_growth

    list(
        mu_const = mu_const,
        mu_gamma = mu_gamma,
        sigma0 = bio_pars$growth$sigma0,
        sigma1 = bio_pars$growth$sigma1,
        h0 = bio_pars$mortality$h0,
        beta = bio_pars$mortality$beta,
        recruit_meanlog = bio_pars$recruitment$meanlog,
        recruit_sdlog = bio_pars$recruitment$sdlog,
        recruit_max_dbh = bio_pars$recruitment$recruit_max_dbh,
        recruit_lambda = bio_pars$recruitment$lambda,
        max_shrink = max_shrink0,
        k_shrink = k_shrink0,
        # Extreme-growth guardrails (upper tail)
        max_growth = if (!is.null(max_growth0)) max_growth0 else Inf,
        max_growth_soft = if (!is.null(max_growth_soft0)) max_growth_soft0 else Inf,
        k_growth = if (!is.null(k_growth0)) k_growth0 else 0,
        use_measurement_error = use_me,
        meas_sd1_a = me$sd1_a,
        meas_sd1_b = me$sd1_b,
        meas_sd2 = me$sd2,
        meas_p_big = me$p_big
    )
}

default_param_grids <- function(base_args, n = 200L) {
    # Build default sweep grids for each parameter.
    #
    # Philosophy
    # - Positive scale parameters: log-spaced factors around the baseline.
    # - Location parameters: linear range around the baseline.
    # - max_shrink (typically negative): sweep from more-negative to 0.
    # - Probabilities: linear in (0,1) avoiding exact 0/1.
    #
    # You can pass your own grids into build_all_sweeps() if you want different ranges.

    n <- as.integer(n)
    if (!is.finite(n) || n < 5L) n <- 200L

    log_grid <- function(center, f_lo = 0.2, f_hi = 5, floor = 1e-8) {
        if (!is.finite(center) || is.na(center) || center <= 0) {
            return(exp(seq(log(floor), log(max(f_hi, floor * 10)), length.out = n)))
        }
        lo <- max(center * f_lo, floor)
        hi <- max(center * f_hi, lo * 1.01)
        exp(seq(log(lo), log(hi), length.out = n))
    }

    lin_grid <- function(center, span = NULL, lo = NULL, hi = NULL) {
        if (!is.finite(center) || is.na(center)) center <- 0
        if (!is.null(lo) && !is.null(hi)) {
            return(seq(lo, hi, length.out = n))
        }
        if (is.null(span)) {
            span <- if (center == 0) 1 else abs(center)
        }
        seq(center - span, center + span, length.out = n)
    }

    # Pull baselines
    b <- base_args

    # Defensive defaults for backward compatibility (in case base_args comes
    # from an older script that didn't include these fields).
    mu_gamma0 <- b$mu_gamma
    if (is.null(mu_gamma0) || !is.finite(mu_gamma0)) mu_gamma0 <- 0
    meas_sd1_a0 <- b$meas_sd1_a
    if (is.null(meas_sd1_a0) || !is.finite(meas_sd1_a0) || meas_sd1_a0 <= 0) meas_sd1_a0 <- 0.0062
    meas_sd1_b0 <- b$meas_sd1_b
    if (is.null(meas_sd1_b0) || !is.finite(meas_sd1_b0) || meas_sd1_b0 <= 0) meas_sd1_b0 <- 0.0904
    meas_sd20 <- b$meas_sd2
    if (is.null(meas_sd20) || !is.finite(meas_sd20) || meas_sd20 <= 0) meas_sd20 <- 4.64
    meas_p_big0 <- b$meas_p_big
    if (is.null(meas_p_big0) || !is.finite(meas_p_big0)) meas_p_big0 <- 0.05

    # max_shrink should generally be <= 0
    max_shrink_grid <- if (is.finite(b$max_shrink) && b$max_shrink < 0) {
        seq(min(b$max_shrink * 3, b$max_shrink - 0.05), 0, length.out = n)
    } else {
        seq(-2, 0, length.out = n)
    }

    recruit_max_grid <- if (is.finite(b$recruit_max_dbh) && b$recruit_max_dbh > 0) {
        seq(max(0.1, 0.5 * b$recruit_max_dbh), 2 * b$recruit_max_dbh, length.out = n)
    } else {
        seq(0.1, 10, length.out = n)
    }

    list(
        # growth
        mu_const = lin_grid(b$mu_const, span = if (is.finite(b$mu_const) && b$mu_const != 0) abs(b$mu_const) else 1),
        mu_gamma = lin_grid(mu_gamma0, span = if (mu_gamma0 != 0) 2 * abs(mu_gamma0) else 0.5),
        sigma0 = log_grid(b$sigma0, f_lo = 0.2, f_hi = 5),
        sigma1 = log_grid(b$sigma1, f_lo = 0.2, f_hi = 10),
        max_shrink = max_shrink_grid,
        k_shrink = log_grid(b$k_shrink, f_lo = 0.1, f_hi = 10),
        max_growth = {
            mg <- b$max_growth
            if (is.null(mg) || !is.finite(mg) || is.na(mg)) mg <- 5
            seq(0, max(mg * 2, 10), length.out = n)
        },
        max_growth_soft = {
            mg <- b$max_growth_soft
            if (is.null(mg) || !is.finite(mg) || is.na(mg)) {
                mg_h <- b$max_growth
                mg <- if (!is.null(mg_h) && is.finite(mg_h) && !is.na(mg_h)) 0.8 * mg_h else 3
            }
            seq(0, max(mg * 2, 10), length.out = n)
        },
        k_growth = log_grid({
            kg <- b$k_growth
            if (is.null(kg) || !is.finite(kg) || is.na(kg) || kg <= 0) kg <- b$k_shrink
            if (is.null(kg) || !is.finite(kg) || is.na(kg) || kg <= 0) kg <- 50
            kg
        }, f_lo = 0.1, f_hi = 10),
        # measurement error mixture (Condit-style)
        meas_sd1_a = log_grid(meas_sd1_a0, f_lo = 0.2, f_hi = 5),
        meas_sd1_b = log_grid(meas_sd1_b0, f_lo = 0.2, f_hi = 5),
        meas_sd2 = log_grid(meas_sd20, f_lo = 0.2, f_hi = 5),
        meas_p_big = {
            p0 <- meas_p_big0
            if (!is.finite(p0) || is.na(p0)) p0 <- 0.05
            p0 <- min(max(p0, 1e-4), 1 - 1e-4)
            lo <- max(1e-4, p0 * 0.1)
            hi <- min(1 - 1e-4, max(p0 * 5, lo * 1.01))
            seq(lo, hi, length.out = n)
        },
        # mortality
        h0 = log_grid(b$h0, f_lo = 0.1, f_hi = 10),
        beta = lin_grid(b$beta, span = 0.2),
        # recruitment
        recruit_meanlog = lin_grid(b$recruit_meanlog, span = 1),
        recruit_sdlog = log_grid(b$recruit_sdlog, f_lo = 0.2, f_hi = 5),
        recruit_max_dbh = recruit_max_grid,
        recruit_lambda = log_grid(b$recruit_lambda, f_lo = 0.05, f_hi = 20)
    )
}

build_all_sweeps <- function(
  scenarios,
  interval_years,
  base_args,
  grids = NULL,
  params = NULL,
  eps_tiebreak = 1e-6,
  hard_penalty = 1e6,
  abs_jump = 1e3
) {
    # Create all sensitivity dt's:
    # - for each scenario in `scenarios`
    # - for each parameter in `params`
    # - sweep across `grids[[param]]`
    #
    # Returns a list with:
    # - dts: named list of per-sweep data.tables
    # - all: one combined data.table
    # - jumps: combined jump table (where abs(diff(total)) >= abs_jump)

    if (is.null(grids)) {
        grids <- default_param_grids(base_args, n = 200L)
    }
    if (is.null(params)) {
        params <- names(grids)
    }
    params <- params[params %in% names(grids)]
    if (length(params) == 0L) {
        stop("No params to sweep.", call. = FALSE)
    }

    sc_names <- names(scenarios)
    if (is.null(sc_names) || any(sc_names == "")) {
        stop("scenarios must be a *named* list (e.g., from make_demo_scenarios()).", call. = FALSE)
    }

    dts <- list()
    jumps <- list()

    for (sc_id in sc_names) {
        sc <- scenarios[[sc_id]]
        if (!all(c("t", "tp1") %in% names(sc))) {
            next
        }

        # Store a compact text representation of the scenario transition.
        # This makes it easy to print (d0,d1,T) on every plot without carrying
        # vector columns.
        sc_d0_txt <- paste(sc$t, collapse = ",")
        sc_d1_txt <- paste(sc$tp1, collapse = ",")

        for (p in params) {
            vals <- grids[[p]]
            dt <- sweep_transition_cost(
                track_dbh_t = sc$t,
                track_dbh_tp1 = sc$tp1,
                interval_years = interval_years,
                base_args = base_args,
                param = p,
                values = vals,
                eps_tiebreak = eps_tiebreak,
                hard_penalty = hard_penalty
            )
            dt[, `:=`(
                scenario = sc_id,
                scenario_name = if (!is.null(sc$name)) as.character(sc$name) else sc_id,
                interval_years = as.numeric(interval_years),
                scenario_d0 = sc_d0_txt,
                scenario_d1 = sc_d1_txt
            )]

            key <- paste0(sc_id, "__", p)
            dts[[key]] <- dt

            j <- detect_jumps(dt, abs_jump = abs_jump)
            if (nrow(j) > 0L) {
                j[, `:=`(
                    scenario = sc_id,
                    scenario_name = if (!is.null(sc$name)) as.character(sc$name) else sc_id
                )]
                jumps[[key]] <- j
            }
        }
    }

    all_dt <- if (length(dts) > 0L) {
        data.table::rbindlist(dts, idcol = "sweep_key", fill = TRUE)
    } else {
        data.table::data.table()
    }

    jumps_dt <- if (length(jumps) > 0L) {
        data.table::rbindlist(jumps, idcol = "sweep_key", fill = TRUE)
    } else {
        data.table::data.table()
    }

    list(dts = dts, all = all_dt, jumps = jumps_dt)
}

plot_all_sweeps_to_pdf <- function(all_dts, pdf_file, facet = FALSE, y_scale = c("absolute", "delta"), subtitle = NULL) {
    # Write a multi-page PDF with one sweep per page.
    #
    # What you will see in the PDF
    # - Each page corresponds to one (scenario, parameter) sweep.
    # - X axis: the value of the parameter being swept (everything else fixed).
    # - Y axis: the transition cost returned by transition_cost_tracks_bio_components(),
    #   decomposed into components.
    #
    # Why "cost" can feel unintuitive
    # - Cost is in negative log-likelihood units (plus any hard penalties).
    # - Lower cost => higher likelihood under the model.
    # - Hard constraints create big discontinuities (e.g., +1e6).
    #
    # y_scale
    # - "absolute": plot raw costs (default).
    # - "delta": plot (cost - cost_at_baseline_parameter_value). This makes it easier
    #   to read sensitivity as "how much worse/better than baseline".
    y_scale <- match.arg(y_scale)
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required.", call. = FALSE)
    }
    if (!requireNamespace("grDevices", quietly = TRUE)) {
        stop("grDevices is required.", call. = FALSE)
    }

    if (!is.list(all_dts) || is.null(all_dts$dts)) {
        stop("Expected output of build_all_sweeps().", call. = FALSE)
    }
    dts <- all_dts$dts
    if (length(dts) == 0L) {
        stop("No sweeps to plot.", call. = FALSE)
    }

    # Slightly taller pages help keep the (3-line) subtitle and caption readable.
    grDevices::pdf(pdf_file, width = 12, height = 9, onefile = TRUE)
    on.exit(grDevices::dev.off(), add = TRUE)

    for (k in names(dts)) {
        dt <- dts[[k]]
        ttl <- paste0("Transition-cost sensitivity - ", unique(dt$scenario_name), " - ", unique(dt$param))
        p <- plot_sweep_components(dt, title = ttl, subtitle = subtitle, facet = facet, y_scale = y_scale)
        print(p)
    }
}

make_demo_scenarios <- function(base_args, interval_years) {
    # Build a few single-track scenarios to probe each case.
    #
    # NOTE: These are *illustrative*; you should also create scenarios based on
    # real observed pairs (d0,d1) from your data.

    interval_years <- as.numeric(interval_years)

    d0 <- 20
    mu_const <- base_args$mu_const
    if (is.null(mu_const) || !is.finite(mu_const)) mu_const <- 0
    mu_gamma <- base_args$mu_gamma
    if (is.null(mu_gamma) || !is.finite(mu_gamma)) mu_gamma <- 0
    mu <- mu_const + mu_gamma * log(max(d0, 1e-6))

    # A plausible growth step
    d1_ok <- d0 + mu * interval_years

    # Soft shrinkage (still above max_shrink)
    d1_shrink_soft <- d0 - 0.2

    # Hard shrinkage: force g < max_shrink
    # threshold is d1 < d0 + max_shrink * interval_years
    d1_shrink_hard <- (d0 + base_args$max_shrink * interval_years) - 1

    # Recruitment examples
    d1_rec_ok <- min(0.9 * base_args$recruit_max_dbh, base_args$recruit_max_dbh - 1e-6)
    if (!is.finite(d1_rec_ok) || d1_rec_ok <= 0) d1_rec_ok <- 1

    d1_rec_bad <- base_args$recruit_max_dbh + 0.1

    # Optional extreme-growth scenarios when guardrails are finite.
    d1_grow_soft <- NA_real_
    if (!is.null(base_args$max_growth_soft) && is.finite(base_args$max_growth_soft)) {
        d1_grow_soft <- d0 + base_args$max_growth_soft * interval_years + 0.2
    }
    d1_grow_hard <- NA_real_
    if (!is.null(base_args$max_growth) && is.finite(base_args$max_growth)) {
        d1_grow_hard <- d0 + base_args$max_growth * interval_years + 0.5
    }

    sc <- list(
        growth_ok = list(name = "DBH->DBH (near mean)", t = c(d0), tp1 = c(d1_ok)),
        shrink_soft = list(name = "DBH->DBH (soft shrink)", t = c(d0), tp1 = c(d1_shrink_soft)),
        shrink_hard = list(name = "DBH->DBH (hard shrink violation)", t = c(d0), tp1 = c(d1_shrink_hard)),
        mortality = list(name = "DBH->NA (death)", t = c(d0), tp1 = c(NA_real_)),
        recruit_ok = list(name = "NA->DBH (recruit ok)", t = c(NA_real_), tp1 = c(d1_rec_ok)),
        recruit_bad = list(name = "NA->DBH (recruit too large)", t = c(NA_real_), tp1 = c(d1_rec_bad)),
        none = list(name = "NA->NA (no recruit)", t = c(NA_real_), tp1 = c(NA_real_))
    )

    if (is.finite(d1_grow_soft)) {
        sc$grow_soft <- list(name = "DBH->DBH (soft extreme growth)", t = c(d0), tp1 = c(d1_grow_soft))
    }
    if (is.finite(d1_grow_hard)) {
        sc$grow_hard <- list(name = "DBH->DBH (hard extreme growth violation)", t = c(d0), tp1 = c(d1_grow_hard))
    }

    sc
}

summarize_components <- function(res) {
    # res is output of transition_cost_tracks_bio_components()
    dt <- res$per_track
    data.table::data.table(
        total = res$total,
        tiebreak = res$tiebreak,
        cost_recruit = sum(dt$cost_recruit),
        cost_no_recruit = sum(dt$cost_no_recruit),
        cost_mortality = sum(dt$cost_mortality),
        cost_growth_lik = sum(dt$cost_growth_lik),
        cost_shrink_soft = sum(dt$cost_shrink_soft),
        cost_growth_soft = sum(dt$cost_growth_soft),
        cost_hard = sum(dt$cost_hard),
        p_recruit = res$p_recruit
    )
}

sweep_transition_cost <- function(
  track_dbh_t,
  track_dbh_tp1,
  interval_years,
  base_args,
  param,
  values,
  eps_tiebreak = 1e-6,
  hard_penalty = 1e6
) {
    # One-at-a-time sweep for a single parameter.
    # Returns a long-ish data.table with total and component costs vs param value.

    if (!exists("transition_cost_tracks_bio_components", mode = "function")) {
        stop("transition_cost_tracks_bio_components() not found. Source dp_global_bio.R first.", call. = FALSE)
    }

    if (!(param %in% names(base_args))) {
        stop("Unknown param: ", param, ". Expected one of: ", paste(names(base_args), collapse = ", "), call. = FALSE)
    }

    values <- as.numeric(values)
    values <- values[is.finite(values)]
    if (length(values) == 0L) {
        stop("values is empty after filtering.", call. = FALSE)
    }

    out <- vector("list", length(values))
    for (i in seq_along(values)) {
        args_i <- base_args
        args_i[[param]] <- values[i]

        res <- do.call(
            transition_cost_tracks_bio_components,
            c(
                list(
                    track_dbh_t = track_dbh_t,
                    track_dbh_tp1 = track_dbh_tp1,
                    interval_years = interval_years,
                    eps_tiebreak = eps_tiebreak,
                    hard_penalty = hard_penalty
                ),
                args_i
            )
        )

        s <- summarize_components(res)
        s[, `:=`(
            param = param,
            value = values[i],
            baseline_value = base_args[[param]]
        )]
        out[[i]] <- s
    }

    data.table::rbindlist(out, fill = TRUE)
}

detect_jumps <- function(dt, abs_jump = 1e3) {
    # Identify large changes in total cost across adjacent grid points.
    # Useful for spotting where hard penalties switch on.
    dt <- data.table::copy(dt)
    data.table::setorder(dt, value)
    dt[, d_total := c(NA_real_, diff(total))]
    dt[abs(d_total) >= abs_jump]
}

param_meaning <- function(param) {
    # Human-readable meaning and (approximate) units for each parameter.
    # This is used for figure labels/captions.
    switch(
        as.character(param),
        mu_const = list(
            label = "mu_const (mean growth)",
            meaning = "Mean annual DBH increment used in the growth likelihood",
            units = "cm/year"
        ),
        mu_gamma = list(
            label = "mu_gamma (size effect on mean growth)",
            meaning = "Slope for size-dependent mean growth: mu(DBH)=mu_const + mu_gamma*log(DBH)",
            units = "(cm/year)/log(cm)"
        ),
        sigma0 = list(
            label = "sigma0 (growth SD intercept)",
            meaning = "Baseline growth standard deviation in sigma(d)=sigma0 + sigma1*d",
            units = "cm/year"
        ),
        sigma1 = list(
            label = "sigma1 (growth SD slope)",
            meaning = "DBH-dependence of growth SD in sigma(d)=sigma0 + sigma1*d",
            units = "(cm/year)/cm"
        ),
        max_shrink = list(
            label = "max_shrink (hard shrink bound)",
            meaning = "Hard lower bound on annual shrinkage; if g < max_shrink then a large penalty is added",
            units = "cm/year"
        ),
        max_growth = list(
            label = "max_growth (hard extreme-growth bound)",
            meaning = "Hard upper bound on annual growth; if g > max_growth then a large penalty is added",
            units = "cm/year"
        ),
        max_growth_soft = list(
            label = "max_growth_soft (soft extreme-growth threshold)",
            meaning = "Soft threshold for annual growth; if g exceeds this, a quadratic penalty is added for the excess",
            units = "cm/year"
        ),
        k_growth = list(
            label = "k_growth (soft extreme-growth penalty)",
            meaning = "Strength of the soft penalty applied when annual growth exceeds max_growth_soft",
            units = "penalty/(cm^2)"
        ),
        k_shrink = list(
            label = "k_shrink (soft shrink penalty)",
            meaning = "Strength of the soft penalty applied when DBH decreases (d1<d0)",
            units = "penalty/(cm^2)"
        ),
        h0 = list(
            label = "h0 (baseline mortality hazard)",
            meaning = "Baseline hazard in hazard(d)=h0*exp(beta*d)",
            units = "1/year"
        ),
        beta = list(
            label = "beta (DBH effect on mortality)",
            meaning = "Effect of DBH on mortality hazard in hazard(d)=h0*exp(beta*d)",
            units = "1/cm"
        ),
        recruit_meanlog = list(
            label = "recruit_meanlog (recruit size meanlog)",
            meaning = "Lognormal meanlog for recruited DBH (NA->DBH)",
            units = "log(cm)"
        ),
        recruit_sdlog = list(
            label = "recruit_sdlog (recruit size sdlog)",
            meaning = "Lognormal sdlog for recruited DBH (NA->DBH)",
            units = "log(cm)"
        ),
        recruit_max_dbh = list(
            label = "recruit_max_dbh (hard recruit max)",
            meaning = "Hard upper bound on recruited DBH; if d1 > recruit_max_dbh then a large penalty is added",
            units = "cm"
        ),
        recruit_lambda = list(
            label = "recruit_lambda (recruitment rate)",
            meaning = "Poisson rate per empty track per year; p_recruit = 1-exp(-lambda*T)",
            units = "1/year"
        ),
        meas_sd1_a = list(
            label = "meas_sd1_a (small-error sd1_a)",
            meaning = "Measurement error (small component): sd1(DBH)=sd1_a + sd1_b*DBH",
            units = "cm"
        ),
        meas_sd1_b = list(
            label = "meas_sd1_b (small-error sd1_b)",
            meaning = "Measurement error (small component): sd1(DBH)=sd1_a + sd1_b*DBH",
            units = "1"
        ),
        meas_sd2 = list(
            label = "meas_sd2 (big-error sd2)",
            meaning = "Measurement error (big component): constant SD for occasional large errors",
            units = "cm"
        ),
        meas_p_big = list(
            label = "meas_p_big (big-error probability)",
            meaning = "Mixture probability of the big measurement-error component",
            units = "probability"
        ),
        # fallback
        list(
            label = as.character(param),
            meaning = "Model parameter",
            units = ""
        )
    )
}

plot_sweep_components <- function(
    dt,
    title = NULL,
    subtitle = NULL,
    caption = NULL,
    facet = FALSE,
    y_scale = c("absolute", "delta")
) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required. Install it with install.packages('ggplot2')", call. = FALSE)
    }

    dt <- data.table::copy(dt)
    # Allow dt coming either from sweep_transition_cost() (no scenario columns)
    # or from build_all_sweeps() (has scenario/scenario_name).
    if (!("scenario_name" %in% names(dt))) dt[, scenario_name := NA_character_]
    if (!("interval_years" %in% names(dt))) dt[, interval_years := NA_real_]
    if (!("scenario_d0" %in% names(dt))) dt[, scenario_d0 := NA_character_]
    if (!("scenario_d1" %in% names(dt))) dt[, scenario_d1 := NA_character_]

    dt_long <- data.table::melt(
        dt,
        id.vars = c("param", "value", "baseline_value"),
        measure.vars = c(
            "total",
            "cost_hard",
            "cost_growth_lik",
            "cost_growth_soft",
            "cost_shrink_soft",
            "cost_mortality",
            "cost_recruit",
            "cost_no_recruit",
            "tiebreak"
        ),
        variable.name = "component",
        value.name = "cost"
    )

    # Human-friendly labels + ordering
    comp_levels <- c(
        "total",
        "cost_hard",
        "cost_growth_lik",
        "cost_growth_soft",
        "cost_shrink_soft",
        "cost_mortality",
        "cost_recruit",
        "cost_no_recruit",
        "tiebreak"
    )
    comp_labels <- c(
        total = "Total",
        cost_hard = "Hard penalty",
        cost_growth_lik = "Growth likelihood",
        cost_growth_soft = "Extreme growth (soft)",
        cost_shrink_soft = "Shrinkage (soft)",
        cost_mortality = "Mortality",
        cost_recruit = "Recruitment",
        cost_no_recruit = "No recruit (NA->NA)",
        tiebreak = "Tie-break"
    )
    dt_long[, component := factor(component, levels = comp_levels, labels = unname(comp_labels[comp_levels]))]

    # Visually emphasize Total
    dt_long[, is_total := (as.character(component) == "Total")]

    # ------------------------------------------------------------------
    # Make the Y axis interpretation explicit
    # ------------------------------------------------------------------
    # Users often find the absolute cost hard to interpret, because:
    # - It is not a probability; it is (approximately) -log(probability or density).
    # - It includes additive constants and (possibly) huge hard penalties.
    #
    # To make "sensitivity" clearer, we optionally plot delta-cost relative to the
    # baseline parameter value:
    #   delta_cost(x) = cost(x) - cost(baseline)
    #
    # This preserves the shape and discontinuities but makes it obvious whether
    # moving away from baseline increases/decreases cost.
    y_scale <- match.arg(y_scale)

    baseline_x <- unique(dt$baseline_value)
    baseline_x <- baseline_x[is.finite(baseline_x)]
    baseline_x <- if (length(baseline_x) > 0L) baseline_x[[1L]] else NA_real_

    # Compute per-component baseline cost at the baseline_x.
    # NOTE: baseline_x might not be exactly on the grid; we use linear interpolation
    # over the swept values.
    if (isTRUE(y_scale == "delta") && is.finite(baseline_x)) {
        base_by_comp <- dt_long[
            , {
                o <- order(value)
                bx <- baseline_x
                by <- tryCatch(
                    stats::approx(x = value[o], y = cost[o], xout = bx, rule = 2)$y,
                    error = function(e) NA_real_
                )
                list(baseline_cost = as.numeric(by))
            },
            by = component
        ]
        dt_long <- merge(dt_long, base_by_comp, by = "component", all.x = TRUE)
        dt_long[, cost_raw := cost]
        dt_long[, cost := cost_raw - baseline_cost]
    }

    p_name <- unique(dt$param)
    p_name <- as.character(p_name[[1L]])
    p_info <- param_meaning(p_name)

    truncate_vec_txt <- function(x, max_chars = 55L) {
        x <- as.character(x)
        if (is.na(x) || !nzchar(x)) return(NA_character_)
        if (nchar(x) <= max_chars) return(x)
        paste0(substr(x, 1L, max_chars - 3L), "...")
    }

    d0_txt <- truncate_vec_txt(unique(dt$scenario_d0)[1L])
    d1_txt <- truncate_vec_txt(unique(dt$scenario_d1)[1L])
    T_txt <- unique(dt$interval_years)
    T_txt <- T_txt[is.finite(T_txt)]
    T_txt <- if (length(T_txt) > 0L) as.character(T_txt[[1L]]) else NA_character_

    transition_line <- NA_character_
    if (!is.na(d0_txt) && !is.na(d1_txt)) {
        # Use 'd0'/'d1' terminology in the subtitle as requested.
        # (NA is printed as NA; values are in cm for DBH.)
        transition_line <- paste0(
            "Transition (per track): d0 = ", d0_txt, " cm  ->  d1 = ", d1_txt, " cm",
            if (!is.na(T_txt)) paste0(",  T = ", T_txt, " years") else ""
        )
    }

    if (is.null(title)) {
        title <- paste0("Sensitivity sweep: ", unique(dt$param))
    }

    if (is.null(subtitle)) {
        # scenario_name is optional; only show it if present.
        scn <- unique(dt$scenario_name)
        scn <- scn[!is.na(scn) & nzchar(scn)]
        scn_txt <- if (length(scn) > 0) paste0("Scenario: ", scn[[1L]], ". ") else ""
        # Force exactly 3 lines so the plot is self-interpreting even in PDFs.
        subtitle_line1 <- paste0(
            scn_txt,
            "Dashed vertical line = baseline estimate. Black line (Total) = sum of components."
        )
        subtitle_line2 <- if (y_scale == "delta") {
            "Delta mode: Delta cost = cost(x) - cost(baseline). Negative = better than baseline; positive = worse."
        } else {
            "Absolute mode: lower cost means a more likely transition under the model."
        }
        subtitle_line3 <- if (!is.na(transition_line)) {
            paste0(
                transition_line,
                ". If only one term applies, Total overlaps that component."
            )
        } else {
            "If only one term applies, Total overlaps that component."
        }

        subtitle <- paste0(subtitle_line1, "\n", subtitle_line2, "\n", subtitle_line3)
    }

    if (is.null(caption)) {
        line1 <- "Interpretation: one-parameter sensitivity of the transition cost."
        line2a <- "All other parameters are fixed at their baseline values (dashed line)."
        line2b <- if (y_scale == "delta") {
            "Delta mode: Delta cost = cost(x) - cost(baseline). The gray horizontal line is y=0. Negative values mean the transition is more likely than at baseline; positive values mean less likely."
        } else {
            "Absolute mode: lower cost means higher likelihood under the model (approximately -log probability/density)."
        }
        unit_txt <- if (!is.null(p_info$units) && nzchar(p_info$units)) paste0(" [", p_info$units, "]") else ""
        line3 <- paste0("X-axis (parameter meaning): ", p_info$label, unit_txt)
        line4 <- paste0("  ", p_info$meaning, ".")

        # Expand caption with an explicit "transition to what state" explanation.
        # Use more lines for readability.
        line5 <- "Transition being scored (per latent track):"
        line6 <- "  state at time t  = d0"
        line7 <- "  state at time t+1 = d1"
        line8 <- "  interval length  = T (years)"
        line9 <- if (!is.na(transition_line)) paste0("  this plot uses: ", transition_line) else ""
        line10 <- "DBH values are in cm. NA means the track is unoccupied (no observed stem assigned)."

        line11 <- "How to read the black Total curve:"
        line12 <- if (y_scale == "delta") {
            "  decrease as x increases -> larger parameter values improve this transition relative to baseline (Delta cost becomes more negative)"
        } else {
            "  smooth decrease as x increases -> larger parameter values make this transition more likely"
        }
        line13 <- if (y_scale == "delta") {
            "  increase as x increases -> larger parameter values worsen this transition relative to baseline (Delta cost becomes more positive)"
        } else {
            "  smooth increase as x increases -> larger parameter values make this transition less likely"
        }
        line14 <- if (y_scale == "delta") {
            "  sudden upward jumps/discontinuities -> a hard constraint likely activated (large penalty turns on)"
        } else {
            "  sudden jumps/discontinuities -> a hard constraint likely activated"
        }

        line15 <- "Common hard constraints (when enabled): max_shrink for DBH->DBH; recruit_max_dbh for NA->DBH."

        lines <- c(
            line1,
            line2a,
            line2b,
            line3,
            line4,
            line5,
            line6,
            line7,
            line8,
            line9,
            line10,
            line11,
            line12,
            line13,
            line14,
            line15
        )
        lines <- lines[!is.na(lines) & nzchar(lines)]
        caption <- paste(lines, collapse = "\n")
    }

    # A small, consistent palette (Total in black)
    pal <- c(
        "Total" = "#111111",
        "Hard penalty" = "#D55E00",
        "Growth likelihood" = "#0072B2",
        "Extreme growth (soft)" = "#E69F00",
        "Shrinkage (soft)" = "#CC79A7",
        "Mortality" = "#009E73",
        "Recruitment" = "#56B4E9",
        "No recruit (NA->NA)" = "#999999",
        "Tie-break" = "#F0E442"
    )

    x_lab <- paste0(p_info$label, if (!is.null(p_info$units) && nzchar(p_info$units)) paste0(" (", p_info$units, ")") else "")

    # Draw-order control: draw Total first (behind), then components.
    # This matters because in many scenarios Total == one active component.
    dt_long[, draw_order := ifelse(is_total, 0L, 1L)]
    data.table::setorder(dt_long, draw_order)

    # 2-line y label improves readability and avoids overly-long axis text.
    y_lab <- if (y_scale == "delta") {
        "Delta cost vs baseline\n(-log-likelihood units + penalties)"
    } else {
        "Cost\n(-log-likelihood units + penalties)"
    }

    p <- ggplot2::ggplot(dt_long, ggplot2::aes(x = value, y = cost, color = component)) +
        ggplot2::geom_line(ggplot2::aes(linewidth = is_total, alpha = is_total)) +
        ggplot2::scale_linewidth_manual(values = c(`FALSE` = 1.0, `TRUE` = 1.2), guide = "none") +
        ggplot2::scale_alpha_manual(values = c(`FALSE` = 1.0, `TRUE` = 0.55), guide = "none") +
        ggplot2::geom_vline(
            xintercept = unique(dt$baseline_value),
            linetype = "dashed",
            color = "black",
            linewidth = 0.6
        ) +
        {
            if (y_scale == "delta") ggplot2::geom_hline(yintercept = 0, color = "#666666", linewidth = 0.4) else NULL
        } +
        ggplot2::scale_color_manual(values = pal, drop = FALSE) +
        ggplot2::labs(
            title = title,
            subtitle = subtitle,
            x = x_lab,
            y = y_lab,
            color = "Component",
            caption = caption
        ) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
            legend.position = "bottom",
            legend.box = "vertical",
            panel.grid.minor = ggplot2::element_blank(),
            plot.title = ggplot2::element_text(face = "bold"),
            plot.caption = ggplot2::element_text(size = 9, color = "#444444")
        )

    # Optional: facet by component when users want to read each term separately.
    # This is helpful when one component dominates the scale and flattens the rest.
    if (isTRUE(facet)) {
        p <- p + ggplot2::facet_wrap(~component, scales = "free_y") +
            ggplot2::theme(legend.position = "none")
    }

    return(p)
}

transition_thresholds <- function(track_dbh_t, track_dbh_tp1, interval_years, base_args) {
    # Quick “why did it jump?” helper.
    # Returns a small table with the hard-threshold checks (and the soft growth cap).

    interval_years <- as.numeric(interval_years)

    # Shrinkage hard constraint: g < max_shrink
    # i.e. (d1-d0)/T < max_shrink  <=>  d1 < d0 + max_shrink*T
    thr_shrink_d1 <- track_dbh_t + base_args$max_shrink * interval_years

    # Extreme-growth hard constraint: g > max_growth
    # i.e. (d1-d0)/T > max_growth  <=>  d1 > d0 + max_growth*T
    thr_grow_hard_d1 <- track_dbh_t + base_args$max_growth * interval_years

    # Extreme-growth soft cap: g > max_growth_soft
    thr_grow_soft_d1 <- track_dbh_t + base_args$max_growth_soft * interval_years

    data.table::data.table(
        track = seq_along(track_dbh_t),
        d0 = track_dbh_t,
        d1 = track_dbh_tp1,
        shrink_hard_threshold_d1 = thr_shrink_d1,
        violates_shrink_hard = (!is.na(track_dbh_t) & !is.na(track_dbh_tp1) & track_dbh_tp1 < thr_shrink_d1),
        growth_soft_threshold_d1 = thr_grow_soft_d1,
        exceeds_growth_soft = (!is.na(track_dbh_t) & !is.na(track_dbh_tp1) & track_dbh_tp1 > thr_grow_soft_d1),
        growth_hard_threshold_d1 = thr_grow_hard_d1,
        violates_growth_hard = (!is.na(track_dbh_t) & !is.na(track_dbh_tp1) & track_dbh_tp1 > thr_grow_hard_d1),
        recruit_max_dbh = base_args$recruit_max_dbh,
        violates_recruit_max = (is.na(track_dbh_t) & !is.na(track_dbh_tp1) & track_dbh_tp1 > base_args$recruit_max_dbh)
    )
}

############################################################
### Optional demo (opt-in)
############################################################
# Set option and source this file to run a quick demo on the simulated CSV.
#   options(dp_global_biol.run_sensitivity_example = TRUE)
#   source("dp_global/R/sensitivity_transition_cost_bio.R")

if (isTRUE(getOption("dp_global_biol.run_sensitivity_example", FALSE))) {
    # Ensure we have the model functions
    if (!exists("estimate_bio_pars", mode = "function")) {
        stop("estimate_bio_pars() not found. Source dp_global_bio.R first.")
    }

    input_file <- "../data_simulation/data/simulation_legacy_backup/simulated_data_two_species.csv"
    xraw <- data.table::fread(input_file)
    xraw[, species := "all"]

    interval_years <- 5
    bio <- estimate_bio_pars(xraw, interval_years = interval_years, mortality_start = c(log(0.01), 0.03))
    base <- bio_pars_to_transition_args(bio)

    sc <- make_demo_scenarios(base, interval_years = interval_years)

    # Sweep a couple of key parameters and write plots
    dt1 <- sweep_transition_cost(sc$growth_ok$t, sc$growth_ok$tp1, interval_years, base, "sigma0", seq(0.02, 0.6, length.out = 120))
    dt2 <- sweep_transition_cost(sc$recruit_ok$t, sc$recruit_ok$tp1, interval_years, base, "recruit_lambda", exp(seq(log(1e-4), log(0.5), length.out = 120)))
    dt3 <- sweep_transition_cost(sc$shrink_hard$t, sc$shrink_hard$tp1, interval_years, base, "max_shrink", seq(-2, 0, length.out = 120))

    p1 <- plot_sweep_components(dt1, title = sc$growth_ok$name)
    p2 <- plot_sweep_components(dt2, title = sc$recruit_ok$name)
    p3 <- plot_sweep_components(dt3, title = sc$shrink_hard$name)

    pdf("./transition_cost_sensitivity_demo.pdf", width = 12, height = 7, onefile = TRUE)
    print(p1)
    print(p2)
    print(p3)
    grDevices::dev.off()
}

