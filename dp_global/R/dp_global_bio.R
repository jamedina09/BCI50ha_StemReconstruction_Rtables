############################################################
# dp_global_bio.R
# Biological parameter estimation and transition costs
############################################################

estimate_bio_pars <- function(
  x,
  interval_years = NULL,
  interval_col_candidates = c("Bio_IntervalYears", "IntervalYears", "interval_years", "census_interval_years", "CensusIntervalYears"),
  census_ids = NULL,
  mortality_start = c(log(0.01), 0),
  # -----------------------------------------------------------------
  # Measurement error model (DBH remeasurement)
  # -----------------------------------------------------------------
  # This is used in TWO places:
  #   (A) to correct growth-variance estimation (separating process vs measurement)
  #   (B) to set conservative guardrails (e.g., max_shrink)
  #
  # Reference (used as the basis for the measurement-error scaling):
  #   https://royalsocietypublishing.org/rstb/article/359/1443/409/20356/Error-propagation-and-scaling-for-tropical-forest
  #
  # Model assumption (mixture, per census measurement):
  #   With probability (1 - p_big): small measurement error
  #   With probability p_big: large measurement error (blunders)
  #
  # Small-error SD (cm) is diameter-dependent:
  #   SD1(D) = a * D + b
  # where D is DBH in cm.
  # (In our workflow, the fitted values are: a=0.0062, b=0.0904.)
  meas_sd1_a = 0.0062,
  meas_sd1_b = 0.0904,
  # Large-error SD (cm) and mixture weight
  meas_sd2 = 4.64,
  meas_p_big = 0.05,
  # Whether to correct growth-variance estimation for measurement error
  use_measurement_error = TRUE,
  # Quantiles used to set conservative guardrails
  shrink_hard_prob = 1e-4,
  shrink_data_quantile = 0.001,
  # Hard shrink guardrail (max_shrink)
  # - "data": estimated from observed shrink tail (with measurement-error support)
  # - "fixed": use a fixed constant bound (cm/year)
  max_shrink_source = c("data", "fixed"),
  max_shrink_fixed = -2,
  # Soft shrinkage penalty strength (k_shrink)
  # - "data": estimate from measurement-error scale (preferred) or from data variance
  # - "fixed": use a fixed constant (units: 1/cm^2)
  k_shrink_source = c("data", "fixed"),
  k_shrink_fixed = 50,
  # Soft extreme-growth penalty strength (k_growth), analogous to k_shrink
  # - "data": estimate from measurement-error scale (preferred) or from data variance
  # - "fixed": use a fixed constant (units: 1/cm^2); set 0 to disable soft penalty
  k_growth_source = c("data", "fixed"),
  k_growth_fixed = 50,
  # Extreme-growth guardrails (upper tail)
  # - growth_hard_prob is the *upper-tail* probability (e.g., 1e-4 means 99.99th percentile)
  # - growth_data_quantile is the empirical upper quantile used as a guardrail
  # - growth_soft_quantile sets a softer threshold used for a quadratic penalty
  growth_hard_prob = 1e-4,
  growth_data_quantile = 0.999,
  growth_soft_quantile = 0.99,
  # Hard growth guardrail (max_growth)
  # - "data": estimated from observed extreme-growth tail
  # - "fixed": use a fixed constant bound (cm/year)
  max_growth_source = c("data", "fixed"),
  max_growth_fixed = 7.5,
  # Recruitment max DBH (upper bound for recruits)
  recruit_max_quantile = 0.999,
  # Recruit max DBH source and fixed value (allow 'data' or 'fixed')
  recruit_max_source = c("data", "fixed"),
  recruit_max_fixed = 5
) {
    # =====================================================================
    # estimate_bio_pars()
    # =====================================================================
    # Goal
    #   Estimate a set of "biologically plausible" parameters that make the
    #   DP stem-tracking likelihood behave sensibly on a given dataset.
    #
    # What this function returns
    #   A nested list with:
    #   - growth:
    #       mu     : mean annual diameter increment (cm / year)
    #       sigma0 : baseline SD of annual increment (cm / year)
    #       sigma1 : slope for SD vs DBH (cm / year per cm DBH)
    #   - mortality:
    #       h0, beta : parameters of a hazard model for DBH -> NA transitions
    #   - recruitment:
    #       meanlog, sdlog : lognormal parameters for recruit DBH (cm)
    #       recruit_max_dbh: hard-ish upper guardrail for recruit DBH (cm)
    #       lambda         : recruit rate per available slot per year
    #   - shrinkage:
    #       k_shrink  : soft penalty weight for shrinkage (1/cm^2)
    #       max_shrink: conservative lower bound on annual growth (cm / year)
    #   - measurement_error:
    #       echo of SD1/SD2/p_big settings used
    #   - settings:
    #       echo of use_measurement_error
    #
    # How these parameters are used later
    #   They are passed into transition_cost_tracks_bio(), which scores a single
    #   transition between two adjacent censuses. The DP solver then finds the
    #   lowest-cost set of trajectories consistent with an anchor census.
    #
    # Key modeling assumptions
    #   1) Growth increments (annualized) are approximately Normal with
    #      heteroskedastic SD and a DBH-dependent mean:
    #        g = (DBH_{t+1} - DBH_t) / T
    #        g | DBH_t ~ Normal(mu(DBH_t), sigma(DBH_t)^2)
    #        mu(DBH)    = alpha + gamma * log(DBH)
    #        sigma(DBH) = sigma0 + sigma1 * DBH
    #   2) DBH measurement has nontrivial noise. If enabled, we treat measured
    #      DBH as:
    #        DBH_obs = DBH_true + epsilon
    #      where epsilon is a mixture of a "small" Normal error and a "large"
    #      Normal error (blunders). The small SD increases with DBH.
    #   3) Mortality is modeled as a hazard over the census interval:
    #        hazard(DBH) = h0 * exp(beta * DBH)
    #        P(death over interval T) = 1 - exp(-hazard(DBH) * T)
    #   4) Recruitment is modeled as NA -> DBH with:
    #        - recruit sizes ~ LogNormal(meanlog, sdlog)
    #        - recruit rate lambda per available NA slot per year
    #
    # Practical intent
    #   This is NOT meant to be a perfect ecological model. It is meant to:
    #   - make the cost function realistic enough that the DP doesn't "cheat"
    #     by preferring ID swaps / forced deaths / forced recruits.
    #   - provide reasonable defaults for sensitivity analysis.

    # ---------------------------------------------------------------------
    # Inputs / arguments (detailed)
    # ---------------------------------------------------------------------
    # x
    #   A data.frame/data.table with at least:
    #     - Tag        : group identifier (integer-like)
    #     - CensusID   : census index (integer-like, increasing)
    #     - DBH        : diameter at breast height (cm)
    #     - TrueStemID : a "ground truth" stem identity used ONLY for parameter estimation
    #     - species    : species code (optional; if missing, caller should add)
    #
    # interval_years
    #   Numeric scalar; time between adjacent censuses used for annualization.
    #   If NULL (default), the function will try to read interval information from
    #   the input data.table via one of the candidate column names
    #   ("Bio_IntervalYears", "IntervalYears", "interval_years", "census_interval_years", "CensusIntervalYears").
    #   When a per-census/row interval column is present, per-pair intervals are
    #   used to annualize increments and to compute mortality/recruitment exposure.
    #   Example: if censuses are 5 years apart, interval_years=5.
    #
    # census_ids
    #   Optional integer vector of CensusID values to use. If NULL, uses all
    #   CensusIDs present after filtering.
    #
    # mortality_start
    #   Starting values for optim() in the mortality fit:
    #     c(log(h0_start), beta_start)
    #
    # meas_sd1_a, meas_sd1_b, meas_sd2, meas_p_big
    #   Parameters of the DBH measurement error model.
    #   We use the (linear-in-DBH) small-error SD:
    #     SD1(D) = a*D + b
    #   and a "large error" SD2 (cm) with mixture probability p_big.
    #
    #   These values follow the remeasurement-based error propagation discussion
    #   in the paper linked above (Royal Society Phil. Trans. B article).
    #   We keep the model in code form only (not reproducing paper text).
    #
    # use_measurement_error
    #   If TRUE:
    #     - subtract expected measurement variance when estimating growth SD
    #     - set shrinkage guardrails informed by measurement error
    #   If FALSE:
    #     - treat all growth variability as process variability
    #     - set shrinkage guardrails purely from data quantiles
    #
    # shrink_hard_prob
    #   Lower-tail probability used for the measurement-noise-only shrink quantile.
    #   Smaller values make max_shrink more permissive (more negative).
    #
    # shrink_data_quantile
    #   Lower quantile of observed annual increments used as a data-driven guardrail.
    #   Smaller values make max_shrink more permissive (more negative).
    #
    # recruit_max_quantile
    #   Upper quantile of observed recruits used as a guardrail for recruit size.
    #   Larger values make recruit_max_dbh more permissive.

    library(data.table)
    library(MASS)

    # Interval handling
    # - If `interval_years` (scalar) is provided it will be used as the default.
    # - If `interval_years` is NULL, we will look for an interval column in the
    #   input data (one of `interval_col_candidates`). If present we will use
    #   per-pair / per-row interval values when computing annualized increments,
    #   mortality exposure and recruitment exposure. If no interval info is
    #   available, the function will error.
    interval_years_provided <- !is.null(interval_years)
    if (interval_years_provided) {
        interval_years <- as.numeric(interval_years)
        if (!is.finite(interval_years) || interval_years <= 0) {
            stop("interval_years must be positive.", call. = FALSE)
        }
    } else {
        # Defer discovery of an interval column until after input cleaning / casting so
        # we can align any candidate interval column with the wide DBH table.
        # (See later where `iw` is constructed.)
        interval_years <- NULL
    }

    # ---------------------------------------------------------------------
    # 1) Clean input
    # ---------------------------------------------------------------------
    # We only use rows where:
    # - DBH is observed and positive
    # - TrueStemID is present (this function is estimating parameters from
    #   "known tracks" based on TrueStemID)
    #
    # This is a key point: estimate_bio_pars() uses TrueStemID to assemble
    # empirical growth/mortality/recruitment signals. If TrueStemID is missing
    # or unreliable, you should not trust the resulting parameters.
    x_dt <- as.data.table(x)
    x_dt <- x_dt[!is.na(DBH) & DBH > 0 & !is.na(TrueStemID)]

    # Ensure a `species` column exists for consistent downstream handling
    if (!("species" %in% names(x_dt))) {
        x_dt[, species := NA_character_]
    } else {
        x_dt[, species := as.character(species)]
    }

    dt <- x_dt[, .(
        Tag,
        TrueStemID = as.integer(TrueStemID),
        CensusID = as.integer(CensusID),
        DBH = as.numeric(DBH),
        species = species
    )]

    if (nrow(dt) == 0) {
        stop("No usable rows after filtering.", call. = FALSE)
    }

    if (is.null(census_ids)) {
        census_ids <- sort(unique(dt$CensusID))
    }

    census_ids <- census_ids[is.finite(census_ids)]
    if (length(census_ids) < 2) {
        stop("Need at least two censuses.", call. = FALSE)
    }

    # ---------------------------------------------------------------------
    # 2) Wide format
    # ---------------------------------------------------------------------
    # Convert to one row per (Tag, TrueStemID, species) with one column per census.
    # This makes it easy to compute adjacent-census transitions.
    dw <- dcast(
        dt,
        Tag + TrueStemID + species ~ CensusID,
        value.var = "DBH"
    )

    # ---------------------------------------------------------------------
    # 3) Growth increments
    # ---------------------------------------------------------------------
    # For each adjacent census pair (t0, t1):
    #   d0 = DBH at t0 (cm)
    #   d1 = DBH at t1 (cm)
    #   g  = (d1 - d0) / T  (cm/year)  where T may be a scalar or vary per pair
    #
    # We collect:
    #   g_all           : all observed annualized increments
    #   d0_all, d1_all  : the corresponding DBHs for regression/diagnostics
    #   var_meas_g_all  : expected measurement-variance contribution to g
    #   T_all           : corresponding interval (years) used for each g
    g_all <- c()
    d0_all <- c()
    d1_all <- c()
    var_meas_g_all <- c()
    T_all <- c()
    # Diagnostic counters for interval handling (useful for testing)
    pairs_candidate_count <- 0L
    pairs_filled_with_scalar_count <- 0L
    pairs_dropped_count <- 0L

    sd1 <- function(d) {
        d <- as.numeric(d)
        pmax(meas_sd1_a * d + meas_sd1_b, 1e-6)
    }

    # Var(epsilon) under the mixture model.
    meas_var_eps <- function(d) {
        s1 <- sd1(d)
        (1 - meas_p_big) * (s1^2) + meas_p_big * (meas_sd2^2)
    }

    # If an interval column exists in the input data (one of `interval_col_candidates`),
    # cast it wide to align with `dw` so we can compute per-pair intervals. Use the
    # original input `x` (not `dt`) for this cast so we don't drop additional
    # columns during the DBH-filtering above.
    iw <- NULL
    iw_col <- NULL
    candidates_present <- intersect(interval_col_candidates, names(x))
    if (length(candidates_present) > 0L) {
        iw_col <- candidates_present[[1L]]
        # Build interval-wide table from the original `x` (not `dt`) so we
        # retain any interval columns the user provided. Filter rows to the
        # same set used to build `dw` (i.e., DBH observed & TrueStemID present).
        if (iw_col %in% names(x_dt)) {
            ix <- x_dt[, .(
                Tag,
                TrueStemID = as.integer(TrueStemID),
                species = species,
                CensusID = as.integer(CensusID),
                interval = as.numeric(get(iw_col))
            )]
            if (nrow(ix) > 0L) {
                iw <- dcast(ix, Tag + TrueStemID + species ~ CensusID, value.var = "interval")
            }
        }
    }

    for (i in seq_len(length(census_ids) - 1)) {
        t0 <- as.character(census_ids[i])
        t1 <- as.character(census_ids[i + 1])

        if (!all(c(t0, t1) %in% names(dw))) next

        ok <- !is.na(dw[[t0]]) & !is.na(dw[[t1]])
        if (!any(ok)) next

        d0 <- dw[[t0]][ok]
        d1 <- dw[[t1]][ok]

        # Determine interval(s) for this pair. Prefer per-row values from iw (t1 then t0),
        # else fall back to provided scalar `interval_years`.
        if (!is.null(iw) && (t1 %in% names(iw) || t0 %in% names(iw))) {
            n_ok <- sum(ok)
            T1 <- if (t1 %in% names(iw)) iw[[t1]][ok] else rep(NA_real_, n_ok)
            T0 <- if (t0 %in% names(iw)) iw[[t0]][ok] else rep(NA_real_, n_ok)
            # Per-row coalesce: prefer T1, fall back to T0 when T1 is missing
            Tvec <- T1
            missing_idx <- !is.finite(Tvec) | (Tvec <= 0)
            if (any(missing_idx)) {
                Tvec[missing_idx] <- T0[missing_idx]
            }
            # Fill remaining missing/invalid with scalar if provided; otherwise drop those rows
            bad <- !is.finite(Tvec) | (Tvec <= 0)
            if (any(bad)) {
                if (!is.null(interval_years) && is.finite(interval_years) && (interval_years > 0)) {
                    # Fill missing per-row intervals with scalar
                    n_fill <- sum(bad)
                    Tvec[bad] <- as.numeric(interval_years)
                    pairs_filled_with_scalar_count <- pairs_filled_with_scalar_count + n_fill
                } else {
                    # No scalar yet; drop rows with missing/invalid T for now so they
                    # don't prevent inferring a representative interval from other pairs.
                    drop_idx <- which(bad)
                    keep_idx <- which(!bad)
                    if (length(keep_idx) == 0L) {
                        # Nothing to use for this pair; skip
                        pairs_dropped_count <- pairs_dropped_count + length(drop_idx)
                        next
                    }
                    pairs_dropped_count <- pairs_dropped_count + length(drop_idx)
                    d0 <- d0[keep_idx]
                    d1 <- d1[keep_idx]
                    Tvec <- Tvec[keep_idx]
                }
            }

            g <- (d1 - d0) / Tvec
            v_meas_g <- (meas_var_eps(d0) + meas_var_eps(d1)) / (Tvec^2)
            T_all <- c(T_all, Tvec)
        } else {
            if (is.null(interval_years) || !is.finite(interval_years) || interval_years <= 0) {
                stop("No interval information found: provide 'interval_years' or add a per-census interval column (e.g., 'Bio_IntervalYears').", call. = FALSE)
            }
            T <- as.numeric(interval_years)
            g <- (d1 - d0) / T
            v_meas_g <- (meas_var_eps(d0) + meas_var_eps(d1)) / (T^2)
            T_all <- c(T_all, rep(T, length(g)))
        }

        g_all <- c(g_all, g)
        d0_all <- c(d0_all, d0)
        d1_all <- c(d1_all, d1)
        var_meas_g_all <- c(var_meas_g_all, v_meas_g)
        pairs_candidate_count <- pairs_candidate_count + length(g)
    }

    if (length(g_all) < 5) {
        stop("Not enough growth observations.", call. = FALSE)
    }

    # If interval_years was not provided as a scalar, infer a representative
    # interval from the collected per-pair intervals (median) and use that for
    # measurement-noise quantile calculations below.
    if (!interval_years_provided) {
        if (length(T_all) == 0 || !is.finite(median(T_all, na.rm = TRUE))) {
            stop("Unable to infer interval_years from data; provide 'interval_years' or add an interval column.", call. = FALSE)
        }
        interval_years <- as.numeric(median(T_all, na.rm = TRUE))
        if (!is.finite(interval_years) || interval_years <= 0) {
            stop("Inferred interval_years is not a positive finite number.", call. = FALSE)
        }
        message("[estimate_bio_pars] Inferred interval_years = ", interval_years, " from data intervals.")
    }

    # ---------------------------------------------------------------------
    # 3a) Estimate mean growth: mu(DBH) = alpha + gamma*log(DBH)
    # ---------------------------------------------------------------------
    # We fit a simple size-dependent mean model using the starting DBH (d0).
    # For robustness on small datasets, we fall back to a constant mean.
    mu_hat <- mean(g_all)

    ok_mu <- is.finite(g_all) & is.finite(d0_all) & (d0_all > 0)
    if (sum(ok_mu) >= 10L && stats::var(log(d0_all[ok_mu])) > 0) {
        fit_mu <- stats::lm(g_all[ok_mu] ~ log(d0_all[ok_mu]))
        alpha_hat <- as.numeric(stats::coef(fit_mu)[1])
        gamma_hat <- as.numeric(stats::coef(fit_mu)[2])
        if (!is.finite(alpha_hat)) alpha_hat <- mu_hat
        if (!is.finite(gamma_hat)) gamma_hat <- 0
    } else {
        alpha_hat <- mu_hat
        gamma_hat <- 0
    }

    mu_pred <- rep(alpha_hat, length(g_all))
    mu_pred[ok_mu] <- alpha_hat + gamma_hat * log(d0_all[ok_mu])

    # ---------------------------------------------------------------------
    # 3b) Estimate growth variance (sigma0, sigma1)
    # ---------------------------------------------------------------------
    # We want the *process* SD of annual increments, not the observed SD which
    # includes measurement error.
    #
    # Target model:
    #   SD_process(g | d0) = sigma0 + sigma1*d0
    #
    # Observed increments include measurement noise, so we estimate a total SD
    # and then (optionally) subtract expected measurement variance in quadrature:
    #   SD_process^2 ≈ max( SD_total^2 - Var_meas(g), 0 )
    #
    # Robust SD estimation trick
    #   For X ~ Normal(0, sd^2), E|X| = sd*sqrt(2/pi).
    #   Rearranging yields an SD proxy:
    #     sd ≈ |X| * sqrt(pi/2)
    #   We apply this to residuals (g - mu_hat) to reduce sensitivity to outliers.
    resid_abs <- abs(g_all - mu_pred)
    # For Normal(0, sd^2), E|X| = sd*sqrt(2/pi), so sd_hat ≈ |X|*sqrt(pi/2)
    sd_total_hat <- resid_abs * sqrt(pi / 2)

    if (isTRUE(use_measurement_error)) {
        sd_proc_hat <- sqrt(pmax(sd_total_hat^2 - var_meas_g_all, 1e-8))
    } else {
        sd_proc_hat <- pmax(sd_total_hat, 1e-6)
    }

    # Fit a simple linear model for SD vs DBH.
    # Notes:
    # - sigma0_hat is constrained to be positive
    # - sigma1_hat is constrained to be non-negative
    fit_sd <- lm(sd_proc_hat ~ d0_all)
    sigma0_hat <- max(coef(fit_sd)[1], 1e-4)
    sigma1_hat <- max(coef(fit_sd)[2], 0)

    # ---------------------------------------------------------------------
    # 4) Shrinkage penalty (k_shrink)
    # ---------------------------------------------------------------------
    # Shrinkage here means d1 < d0 (negative increment), which can occur due to:
    # - real biological shrinkage (limited)
    # - measurement error (common)
    # - ID swaps / bad matches (the thing we want to discourage)
    #
    # In transition_cost_tracks_bio(), shrinkage is penalized softly as:
    #   cost_shrink = k_shrink * (d0 - d1)^2
    # Units:
    #   (d0 - d1) is cm; to make cost dimensionless, k_shrink has units 1/cm^2.
    #
    # Heuristic used here
    #   If measurement error is enabled, set k_shrink so that shrinkage of about
    #   one typical measurement SD costs O(1).
    #   If measurement error is disabled, fall back to a crude estimate based on
    #   the variance of negative increments.
    # Estimate a soft shrink penalty from measurement error scale.
    # k_shrink is applied as:  k_shrink * (d0-d1)^2  (units: 1/cm^2).
    # With measurement noise, small shrinkage can be expected; we therefore scale
    # k_shrink so that shrinkage of ~1 SD costs O(1).
    k_shrink_source <- match.arg(k_shrink_source)
    k_shrink_fixed <- as.numeric(k_shrink_fixed)
    if (identical(k_shrink_source, "fixed")) {
        if (!is.finite(k_shrink_fixed) || k_shrink_fixed < 0) {
            stop("k_shrink_fixed must be a finite non-negative number when k_shrink_source='fixed'.", call. = FALSE)
        }
        k_shrink_hat_est <- NA_real_
        k_shrink_hat <- k_shrink_fixed
    } else {
        sd_meas_diff <- sqrt((meas_var_eps(d0_all) + meas_var_eps(d1_all)))
        sd_meas_diff <- sd_meas_diff[is.finite(sd_meas_diff) & sd_meas_diff > 0]
        if (isTRUE(use_measurement_error) && length(sd_meas_diff) >= 10) {
            s_typ <- stats::median(sd_meas_diff, na.rm = TRUE)
            k_shrink_hat_est <- 1 / (2 * (s_typ^2))
            k_shrink_hat_est <- min(max(k_shrink_hat_est, 1e-6), 1e6)
        } else {
            # Fallback (no measurement-error model): estimate from observed shrink magnitudes
            # on the *DBH difference* scale (cm), so k_shrink retains units 1/cm^2.
            delta_shrink <- (d0_all - d1_all)[is.finite(d0_all) & is.finite(d1_all) & (d1_all < d0_all)]
            if (length(delta_shrink) >= 5 && is.finite(stats::var(delta_shrink)) && stats::var(delta_shrink) > 0) {
                k_shrink_hat_est <- 1 / (2 * stats::var(delta_shrink))
            } else {
                k_shrink_hat_est <- 50
            }
        }
        k_shrink_hat <- k_shrink_hat_est
    }

    # ---------------------------------------------------------------------
    # 5) Mortality model
    # ---------------------------------------------------------------------
    # Mortality is inferred from TrueStemID tracks as:
    #   alive at t0 (DBH observed) and missing at t1 (DBH NA)
    #
    # We fit a simple hazard model:
    #   hazard(DBH) = h0 * exp(beta * DBH)
    #   P(death over interval T) = 1 - exp(-hazard(DBH)*T)
    #
    # mortality_start is on the unconstrained scale used by optim:
    #   par[1] = log(h0), par[2] = beta
    d0_m <- c()
    died <- c()
    T_m_vec <- c()

    for (i in seq_len(length(census_ids) - 1)) {
        t0 <- as.character(census_ids[i])
        t1 <- as.character(census_ids[i + 1])

        if (!all(c(t0, t1) %in% names(dw))) next

        at_risk <- !is.na(dw[[t0]])
        if (!any(at_risk)) next

        # Determine intervals for these at-risk rows, prefer per-row iw values (t1 then t0), else scalar
        if (!is.null(iw) && (t1 %in% names(iw) || t0 %in% names(iw))) {
            n_at <- sum(at_risk)
            T1m <- if (t1 %in% names(iw)) iw[[t1]][at_risk] else rep(NA_real_, n_at)
            T0m <- if (t0 %in% names(iw)) iw[[t0]][at_risk] else rep(NA_real_, n_at)
            Tm <- T1m
            missing_idx_m <- !is.finite(Tm) | (Tm <= 0)
            if (any(missing_idx_m)) {
                Tm[missing_idx_m] <- T0m[missing_idx_m]
            }
            bad <- !is.finite(Tm) | (Tm <= 0)
            if (any(bad)) {
                if (!is.null(interval_years) && is.finite(interval_years) && (interval_years > 0)) {
                    Tm[bad] <- as.numeric(interval_years)
                } else {
                    stop("Missing or invalid interval values for mortality rows and no scalar 'interval_years' provided.", call. = FALSE)
                }
            }
        } else {
            if (is.null(interval_years) || !is.finite(interval_years) || interval_years <= 0) {
                stop("No interval information found: provide 'interval_years' or add an interval column.", call. = FALSE)
            }
            Tm <- rep(as.numeric(interval_years), sum(at_risk))
        }

        d0_m <- c(d0_m, dw[[t0]][at_risk])
        died <- c(died, as.integer(is.na(dw[[t1]][at_risk])))
        T_m_vec <- c(T_m_vec, Tm)
        # Bookkeeping: how many mortality rows used and how many were filled/dropped
        pairs_candidate_count <- pairs_candidate_count + length(Tm)
        if (exists("T1m") && exists("T0m")) {
            pairs_filled_with_scalar_count <- pairs_filled_with_scalar_count + sum(!is.finite(T1m) & is.finite(T0m) & is.finite(Tm))
            pairs_dropped_count <- pairs_dropped_count + sum(!is.finite(T1m) & !is.finite(T0m) & !is.finite(Tm))
        }
    }

    negloglik_mort <- function(par, d0, died, Tvec) {
        h0 <- exp(par[1])
        beta <- par[2]

        hazard <- h0 * exp(beta * d0)
        p <- 1 - exp(-hazard * Tvec)
        p <- pmin(pmax(p, 1e-12), 1 - 1e-12)

        -sum(died * log(p) + (1 - died) * log(1 - p))
    }

    fit_m <- optim(
        mortality_start,
        negloglik_mort,
        d0 = d0_m,
        died = died,
        Tvec = T_m_vec,
        method = "BFGS"
    )

    h0_hat <- exp(fit_m$par[1])
    beta_hat <- fit_m$par[2]

    # ---------------------------------------------------------------------
    # 6) Recruitment model
    # ---------------------------------------------------------------------
    # We identify "recruitment events" from TrueStemID tracks as:
    #   missing at t0 (DBH NA) and observed at t1 (DBH > 0)
    #
    # We estimate:
    #   - recruit sizes via a lognormal fit
    #   - recruit_max_dbh as a high quantile guardrail
    #   - recruit rate lambda as recruits per available NA slot per year
    recruit_dbh <- c()
    n_risk <- 0
    n_rec <- 0
    total_time_at_risk <- 0

    for (i in seq_len(length(census_ids) - 1)) {
        t0 <- as.character(census_ids[i])
        t1 <- as.character(census_ids[i + 1])

        if (!all(c(t0, t1) %in% names(dw))) next

        at_risk <- is.na(dw[[t0]])
        if (!any(at_risk)) next

        d1_at_risk <- dw[[t1]][at_risk]
        # Recruits must have a positive observed size. Non-positive values will
        # break the lognormal fit and are not meaningful DBH measurements.
        rec <- at_risk
        rec[at_risk] <- !is.na(d1_at_risk) & is.finite(d1_at_risk) & (d1_at_risk > 0)

        recruit_dbh <- c(recruit_dbh, dw[[t1]][rec])

        # Count risk-years using per-row intervals if available (prefer t1 then t0 per-row), else scalar
        if (!is.null(iw) && (t1 %in% names(iw) || t0 %in% names(iw))) {
            n_at <- sum(at_risk)
            T1r <- if (t1 %in% names(iw)) iw[[t1]][at_risk] else rep(NA_real_, n_at)
            T0r <- if (t0 %in% names(iw)) iw[[t0]][at_risk] else rep(NA_real_, n_at)
            T_r <- T1r
            missing_idx_r <- !is.finite(T_r) | (T_r <= 0)
            if (any(missing_idx_r)) {
                T_r[missing_idx_r] <- T0r[missing_idx_r]
            }
            bad <- !is.finite(T_r) | (T_r <= 0)
            if (any(bad)) {
                if (!is.null(interval_years) && is.finite(interval_years) && (interval_years > 0)) {
                    T_r[bad] <- as.numeric(interval_years)
                } else {
                    stop("Missing or invalid interval values for recruitment rows and no scalar 'interval_years' provided.", call. = FALSE)
                }
            }
            total_time_at_risk <- total_time_at_risk + sum(T_r, na.rm = TRUE)
            # Also track any pair-level fills/drops as part of recruitment accounting
            if (exists("T1r") && exists("T0r")) {
                pairs_filled_with_scalar_count <- pairs_filled_with_scalar_count + sum(!is.finite(T1r) & is.finite(T0r))
                pairs_dropped_count <- pairs_dropped_count + sum(!is.finite(T1r) & !is.finite(T0r) & (!(!is.null(interval_years) && is.finite(interval_years) && (interval_years > 0))))
            }
        } else {
            if (is.null(interval_years) || !is.finite(interval_years) || interval_years <= 0) {
                stop("No interval information found: provide 'interval_years' or add an interval column.", call. = FALSE)
            }
            total_time_at_risk <- total_time_at_risk + sum(rep(as.numeric(interval_years), sum(at_risk)))
        }

        n_risk <- n_risk + sum(at_risk)
        n_rec <- n_rec + sum(rec)
    }

    recruit_dbh <- recruit_dbh[is.finite(recruit_dbh) & recruit_dbh > 0]
    # Bookkeeping for recruitment missing intervals
    # (we tracked contributions to total_time_at_risk above)

    if (length(recruit_dbh) >= 2) {
        fit_r <- fitdistr(recruit_dbh, "lognormal")
        mu_r <- fit_r$estimate["meanlog"]
        sd_r <- fit_r$estimate["sdlog"]
    } else {
        mu_r <- log(2)
        sd_r <- 0.5
    }

    # Guardrail: recruit_max_dbh can come from data or be fixed by the user.
    recruit_max_source <- match.arg(recruit_max_source)
    recruit_max_fixed <- as.numeric(recruit_max_fixed)
    if (identical(recruit_max_source, "fixed")) {
        if (!is.finite(recruit_max_fixed) || recruit_max_fixed <= 0) {
            stop("recruit_max_fixed must be a positive finite number when recruit_max_source='fixed'.", call. = FALSE)
        }
        recruit_max_dbh <- recruit_max_fixed
    } else {
        recruit_max_dbh <- if (length(recruit_dbh) > 0) {
            as.numeric(stats::quantile(recruit_dbh, recruit_max_quantile, na.rm = TRUE))
        } else {
            5
        }
    }

    # Recruitment rate (Poisson)
    lambda_hat <- if (total_time_at_risk > 0) {
        n_rec / total_time_at_risk
    } else {
        0
    }

    # ---------------------------------------------------------------------
    # 7) Shrink hard bound (max_shrink) from measurement + data
    # ---------------------------------------------------------------------
    # max_shrink is a conservative lower bound on annual growth (cm/year).
    # It is used as a guardrail to reject *absurd* shrinkage that is very
    # unlikely to be real or due to measurement error.
    #
    # Construction
    #   max_shrink_data : a small quantile of observed growth increments
    #   max_shrink_meas : a lower quantile of the measurement-noise-only
    #                    annualized difference distribution
    #   max_shrink_hat  : the minimum of these (more conservative)
    #
    # Measurement-noise-only lower quantile
    #   We compute a lower quantile of (e1 - e0)/T where e0 and e1 follow the
    #   mixture model. This yields a realistic lower tail for shrinkage driven
    #   purely by measurement error.
    #
    # Important: this is a guardrail, not a hard ecological law.
    # Measurement-informed lower bound on annual growth (mostly to prevent absurd matches).
    # We compute a conservative lower quantile of the *measurement-noise-only* annualized
    # DBH difference, using a typical diameter.
    meas_lower_quantile_g <- function(p, d_typ) {
        p <- as.numeric(p)
        if (!is.finite(p) || p <= 0 || p >= 1) {
            return(NA_real_)
        }
        d_typ <- as.numeric(d_typ)
        if (!is.finite(d_typ) || d_typ <= 0) {
            return(NA_real_)
        }

        s_small <- sd1(d_typ)
        s_big <- meas_sd2
        w_small <- 1 - meas_p_big
        w_big <- meas_p_big

        # Four-component mixture for (e1-e0)/T
        sds <- c(
            sqrt(s_small^2 + s_small^2) / interval_years,
            sqrt(s_small^2 + s_big^2) / interval_years,
            sqrt(s_big^2 + s_small^2) / interval_years,
            sqrt(s_big^2 + s_big^2) / interval_years
        )
        wts <- c(w_small * w_small, w_small * w_big, w_big * w_small, w_big * w_big)

        cdf <- function(x) sum(wts * pnorm(x, mean = 0, sd = sds))

        lo <- -10
        hi <- 0
        # Expand if needed
        if (cdf(lo) > p) {
            lo2 <- -50
            if (cdf(lo2) > p) {
                return(lo2)
            }
            lo <- lo2
        }
        if (cdf(hi) < p) {
            hi <- 10
            if (cdf(hi) < p) {
                return(hi)
            }
        }

        out <- tryCatch(
            uniroot(function(x) cdf(x) - p, lower = lo, upper = hi, tol = 1e-6)$root,
            error = function(e) NA_real_
        )
        out
    }

    d_typ <- stats::median(d0_all[is.finite(d0_all) & d0_all > 0], na.rm = TRUE)
    max_shrink_meas <- if (isTRUE(use_measurement_error)) meas_lower_quantile_g(shrink_hard_prob, d_typ) else NA_real_
    max_shrink_data <- as.numeric(stats::quantile(g_all, shrink_data_quantile, na.rm = TRUE))
    max_shrink_hat_est <- if (isTRUE(use_measurement_error) && is.finite(max_shrink_meas)) {
        min(max_shrink_data, max_shrink_meas)
    } else {
        max_shrink_data
    }

    max_shrink_source <- match.arg(max_shrink_source)
    max_shrink_fixed <- as.numeric(max_shrink_fixed)
    if (identical(max_shrink_source, "fixed")) {
        if (!is.finite(max_shrink_fixed)) {
            stop("max_shrink_fixed must be a finite number when max_shrink_source='fixed'.", call. = FALSE)
        }
        max_shrink_hat <- max_shrink_fixed
    } else {
        max_shrink_hat <- max_shrink_hat_est
    }

    # ---------------------------------------------------------------------
    # 7b) Extreme-growth guardrails (upper tail)
    # ---------------------------------------------------------------------
    # We treat very large positive increments as a likely sign of an incorrect
    # match (ID swap, mis-measurement, etc.). To avoid the DP taking such edges,
    # we set a conservative (permissive) hard upper bound, and also provide a
    # softer threshold for a quadratic penalty.
    #
    # Construction mirrors the shrinkage guardrail, but for the upper tail:
    #   max_growth_data : upper quantile of observed annualized increments
    #   max_growth_meas : upper quantile of the measurement-noise-only mixture for (e1-e0)/T
    #   max_growth_hat  : max(max_growth_data, max_growth_meas)  (more permissive)
    meas_upper_quantile_g <- function(p, d_typ) {
        p <- as.numeric(p)
        if (!is.finite(p) || p <= 0 || p >= 1) {
            return(NA_real_)
        }
        d_typ <- as.numeric(d_typ)
        if (!is.finite(d_typ) || d_typ <= 0) {
            return(NA_real_)
        }

        s_small <- sd1(d_typ)
        s_big <- meas_sd2
        w_small <- 1 - meas_p_big
        w_big <- meas_p_big

        # Four-component mixture for (e1-e0)/T
        sds <- c(
            sqrt(s_small^2 + s_small^2) / interval_years,
            sqrt(s_small^2 + s_big^2) / interval_years,
            sqrt(s_big^2 + s_small^2) / interval_years,
            sqrt(s_big^2 + s_big^2) / interval_years
        )
        wts <- c(w_small * w_small, w_small * w_big, w_big * w_small, w_big * w_big)

        cdf <- function(x) sum(wts * pnorm(x, mean = 0, sd = sds))

        lo <- 0
        hi <- 10
        if (cdf(hi) < p) {
            hi2 <- 50
            if (cdf(hi2) < p) {
                return(hi2)
            }
            hi <- hi2
        }
        if (cdf(lo) > p) {
            lo2 <- -10
            if (cdf(lo2) > p) {
                return(lo2)
            }
            lo <- lo2
        }

        out <- tryCatch(
            uniroot(function(x) cdf(x) - p, lower = lo, upper = hi, tol = 1e-6)$root,
            error = function(e) NA_real_
        )
        out
    }

    growth_hard_prob <- as.numeric(growth_hard_prob)
    growth_data_quantile <- as.numeric(growth_data_quantile)
    growth_soft_quantile <- as.numeric(growth_soft_quantile)

    max_growth_data <- as.numeric(stats::quantile(g_all, growth_data_quantile, na.rm = TRUE))
    max_growth_meas <- if (isTRUE(use_measurement_error) && is.finite(d_typ)) {
        p_hi <- 1 - growth_hard_prob
        meas_upper_quantile_g(p_hi, d_typ)
    } else {
        NA_real_
    }
    max_growth_hat_est <- if (isTRUE(use_measurement_error) && is.finite(max_growth_meas)) {
        max(max_growth_data, max_growth_meas)
    } else {
        max_growth_data
    }

    max_growth_source <- match.arg(max_growth_source)
    max_growth_fixed <- as.numeric(max_growth_fixed)
    if (identical(max_growth_source, "fixed")) {
        if (!is.finite(max_growth_fixed) || max_growth_fixed <= 0) {
            stop("max_growth_fixed must be a finite positive number when max_growth_source='fixed'.", call. = FALSE)
        }
        max_growth_hat <- max_growth_fixed
    } else {
        max_growth_hat <- max_growth_hat_est
    }

    max_growth_soft_data <- as.numeric(stats::quantile(g_all, growth_soft_quantile, na.rm = TRUE))
    max_growth_soft_hat <- if (is.finite(max_growth_hat) && is.finite(max_growth_soft_data)) {
        min(max_growth_hat, max_growth_soft_data)
    } else {
        max_growth_soft_data
    }

    # Soft extreme-growth penalty strength (units: 1/cm^2), analogous to k_shrink.
    # In transition_cost_tracks_bio(), the penalty is applied to the *excess DBH* above
    # the soft cap in cm:
    #   excess = d1 - (d0 + max_growth_soft*T)
    #   cost_growth_soft = k_growth * excess^2
    # We therefore scale k_growth so that an excess of ~1 typical measurement SD costs O(1).
    sd_meas_diff2 <- sqrt((meas_var_eps(d0_all) + meas_var_eps(d1_all)))
    sd_meas_diff2 <- sd_meas_diff2[is.finite(sd_meas_diff2) & sd_meas_diff2 > 0]
    if (isTRUE(use_measurement_error) && length(sd_meas_diff2) >= 10) {
        s_typ2 <- stats::median(sd_meas_diff2, na.rm = TRUE)
        k_growth_hat_est <- 1 / (2 * (s_typ2^2))
        k_growth_hat_est <- min(max(k_growth_hat_est, 1e-6), 1e6)
    } else {
        # Fallback (no measurement-error model): estimate from observed DBH increments (cm)
        # so k_growth retains units 1/cm^2.
        delta_pos <- (d1_all - d0_all)[is.finite(d0_all) & is.finite(d1_all) & (d1_all > d0_all)]
        if (length(delta_pos) >= 5 && is.finite(stats::var(delta_pos)) && stats::var(delta_pos) > 0) {
            k_growth_hat_est <- 1 / (2 * stats::var(delta_pos))
        } else {
            k_growth_hat_est <- 50
        }
    }

    k_growth_source <- match.arg(k_growth_source)
    k_growth_fixed <- as.numeric(k_growth_fixed)
    if (identical(k_growth_source, "fixed")) {
        if (!is.finite(k_growth_fixed) || k_growth_fixed < 0) {
            stop("k_growth_fixed must be a finite non-negative number when k_growth_source='fixed'.", call. = FALSE)
        }
        k_growth_hat <- k_growth_fixed
    } else {
        k_growth_hat <- k_growth_hat_est
    }

    # ---------------------------------------------------------------------
    # 8) Return parameters
    # ---------------------------------------------------------------------
    # The returned structure is intentionally aligned with:
    # - bio_pars_to_transition_args()
    # - transition_cost_tracks_bio()
    # - realism_calibration.R diagnostics
    list(
        growth = list(
            # Mean annual growth model: mu(DBH) = alpha + gamma*log(DBH)
            alpha = alpha_hat,
            gamma = gamma_hat,
            # Backward-compatible summary (empirical mean of g)
            mu = mu_hat,
            sigma0 = sigma0_hat,
            sigma1 = sigma1_hat,
            # Extreme-growth guardrails
            max_growth_soft = max_growth_soft_hat,
            max_growth = max_growth_hat,
            k_growth = k_growth_hat,
            k_growth_source = k_growth_source,
            k_growth_fixed = k_growth_fixed,
            k_growth_estimated_value = k_growth_hat_est,
            # Diagnostics / provenance (analogous to shrinkage)
            max_growth_data_quantile = growth_data_quantile,
            max_growth_data_value = max_growth_data,
            max_growth_meas_prob = growth_hard_prob,
            max_growth_meas_value = max_growth_meas,
            max_growth_soft_quantile = growth_soft_quantile,
            max_growth_soft_value = max_growth_soft_data,
            max_growth_source = max_growth_source,
            max_growth_fixed = max_growth_fixed,
            max_growth_estimated_value = max_growth_hat_est,
            # -----------------------------------------------------------------
            # Uniform nested layout (organizational; flat fields above remain)
            # -----------------------------------------------------------------
            guardrails = list(
                hard = list(
                    value = max_growth_hat,
                    source = max_growth_source,
                    fixed = max_growth_fixed,
                    estimated_value = max_growth_hat_est,
                    data_quantile = growth_data_quantile,
                    data_value = max_growth_data,
                    meas_prob = growth_hard_prob,
                    meas_value = max_growth_meas
                ),
                soft = list(
                    value = max_growth_soft_hat,
                    quantile = growth_soft_quantile,
                    data_value = max_growth_soft_data
                )
            ),
            penalties = list(
                soft = list(
                    k = k_growth_hat,
                    source = k_growth_source,
                    fixed = k_growth_fixed,
                    estimated_value = k_growth_hat_est
                )
            )
        ),
        mortality = list(
            h0   = h0_hat,
            beta = beta_hat
        ),
        recruitment = list(
            meanlog = mu_r,
            sdlog = sd_r,
            recruit_max_dbh = recruit_max_dbh,
            recruit_max_source = recruit_max_source,
            recruit_max_fixed = recruit_max_fixed,
            lambda = lambda_hat
        ),
        shrinkage = list(
            k_shrink = k_shrink_hat,
            max_shrink = max_shrink_hat,
            k_shrink_source = k_shrink_source,
            k_shrink_fixed = k_shrink_fixed,
            k_shrink_estimated_value = k_shrink_hat_est,
            max_shrink_source = max_shrink_source,
            max_shrink_fixed = max_shrink_fixed,
            max_shrink_estimated_value = max_shrink_hat_est,
            max_shrink_data_quantile = shrink_data_quantile,
            max_shrink_data_value = max_shrink_data,
            max_shrink_meas_prob = shrink_hard_prob,
            max_shrink_meas_value = max_shrink_meas,
            # -----------------------------------------------------------------
            # Uniform nested layout (organizational; flat fields above remain)
            # -----------------------------------------------------------------
            guardrails = list(
                hard = list(
                    value = max_shrink_hat,
                    source = max_shrink_source,
                    fixed = max_shrink_fixed,
                    estimated_value = max_shrink_hat_est,
                    data_quantile = shrink_data_quantile,
                    data_value = max_shrink_data,
                    meas_prob = shrink_hard_prob,
                    meas_value = max_shrink_meas
                ),
                # For shrinkage, the soft penalty begins as soon as DBH decreases (d1 < d0),
                # i.e., the soft "threshold" is effectively 0 cm shrink.
                soft = list(
                    value = 0
                )
            ),
            penalties = list(
                soft = list(
                    k = k_shrink_hat,
                    source = k_shrink_source,
                    fixed = k_shrink_fixed,
                    estimated_value = k_shrink_hat_est
                )
            )
        ),
        measurement_error = list(
            sd1_a = meas_sd1_a,
            sd1_b = meas_sd1_b,
            sd2 = meas_sd2,
            p_big = meas_p_big
        ),
        settings = list(
            use_measurement_error = use_measurement_error
        ),
        interval = list(
            inferred_interval_years = interval_years,
            per_pair_intervals = T_all,
            pairs_candidate_count = as.integer(pairs_candidate_count),
            pairs_filled_with_scalar_count = as.integer(pairs_filled_with_scalar_count),
            pairs_dropped_count = as.integer(pairs_dropped_count)
        )
    )
}

estimate_bio_pars_test_v2 <- function(
  x,
  mortality_start = c(log(0.01), 0),
  # -----------------------------------------------------------------
  # Measurement error model (DBH remeasurement)
  # -----------------------------------------------------------------
  # This is used in TWO places:
  #   (A) to correct growth-variance estimation (separating process vs measurement)
  #   (B) to set conservative guardrails (e.g., max_shrink)
  #
  # Reference (used as the basis for the measurement-error scaling):
  #   https://royalsocietypublishing.org/rstb/article/359/1443/409/20356/Error-propagation-and-scaling-for-tropical-forest
  #
  # Model assumption (mixture, per census measurement):
  #   With probability (1 - p_big): small measurement error
  #   With probability p_big: large measurement error (blunders)
  #
  # Small-error SD (cm) is diameter-dependent:
  #   SD1(D) = a * D + b
  # where D is DBH in cm.
  # (In our workflow, the fitted values are: a=0.0062, b=0.0904.)
  meas_sd1_a = 0.0062,
  meas_sd1_b = 0.0904,
  # Large-error SD (cm) and mixture weight
  meas_sd2 = 4.64,
  meas_p_big = 0.05,
  # Whether to correct growth-variance estimation for measurement error
  use_measurement_error = TRUE,
  # Recruit max DBH source and fixed value (allow 'data' or 'fixed')
  recruit_max_source = c("data", "fixed"),
  recruit_max_fixed = 5,
  # Quantiles used to set conservative guardrails
  shrink_hard_prob = 1e-4,
  shrink_data_quantile = 0.001,
  # Hard shrink guardrail (max_shrink)
  # - "data": estimated from observed shrink tail (with measurement-error support)
  # - "fixed": use a fixed constant bound (cm/year)
  max_shrink_source = c("data", "fixed"),
  max_shrink_fixed = -2,
  # Extreme-growth guardrails (upper tail)
  # - growth_hard_prob is the *upper-tail* probability (e.g., 1e-4 means 99.99th percentile)
  # - growth_data_quantile is the empirical upper quantile used as a guardrail
  # - growth_soft_quantile sets a softer threshold used for a quadratic penalty
  growth_hard_prob = 1e-4,
  growth_data_quantile = 0.999,
  growth_soft_quantile = 0.99,
  # Hard growth guardrail (max_growth)
  # - "data": estimated from observed extreme-growth tail
  # - "fixed": use a fixed constant bound (cm/year)
  max_growth_source = c("data", "fixed"),
  max_growth_fixed = 7.5,
  # Soft extreme-growth penalty strength (k_growth), analogous to k_shrink
  # - "data": estimate from measurement-error scale (preferred) or from data variance
  # - "fixed": use a fixed constant (units: 1/cm^2); set 0 to disable soft penalty
  k_growth_source = c("data", "fixed"),
  k_growth_fixed = 50,
  # Soft shrinkage penalty strength (k_shrink)
  # - "data": estimate from measurement-error scale (preferred) or from data variance
  # - "fixed": use a fixed constant (units: 1/cm^2)
  k_shrink_source = c("data", "fixed"),
  k_shrink_fixed = 50,
  # Recruitment max DBH (upper bound for recruits)
  recruit_max_quantile = 0.999 # ,
) {
    # =====================================================================
    # estimate_bio_pars()
    # =====================================================================
    # Estimate biologically plausible parameters (growth, mortality, recruitment)
    # from tracked stems (requires TrueStemID). Returns a nested list with:
    #  - growth: alpha, gamma, sigma0, sigma1, guardrails, penalties
    #  - mortality: h0, beta
    #  - recruitment: lognormal size params and rate lambda
    #  - shrinkage and measurement_error settings
    #
    # Brief models:
    #  - mu(DBH) = alpha + gamma * log(DBH)
    #  - sigma(DBH) = sigma0 + sigma1 * DBH
    #  - measurement error: small/large Normal mixture (SD1(D)=a*D+b)
    #  - mortality hazard: h0 * exp(beta * DBH)
    # Notes:
    #  - This function supplies reasonable defaults and guardrails for DP.
    #  - It is not intended as a full ecological analysis.

    # ---------------------------------------------------------------------
    # Inputs / arguments (detailed)
    # ---------------------------------------------------------------------
    # x
    #   A data.frame/data.table with at least:
    #     - Tag        : group identifier (integer-like)
    #     - CensusID   : census index (integer-like, increasing)
    #     - DBH        : diameter at breast height (cm)
    #     - TrueStemID : a "ground truth" stem identity used ONLY for parameter estimation
    #     - species    : species code (optional; if missing, caller should add)
    #
    # interval_years
    #   Numeric scalar; time between adjacent censuses used for annualization.
    #   If NULL (default), the function will try to read interval information from
    #   the input data.table via one of the candidate column names
    #   ("Bio_IntervalYears", "IntervalYears", "interval_years", "census_interval_years", "CensusIntervalYears").
    #   When a per-census/row interval column is present, per-pair intervals are
    #   used to annualize increments and to compute mortality/recruitment exposure.
    #   Example: if censuses are 5 years apart, interval_years=5.
    #
    # census_ids
    #   Optional integer vector of CensusID values to use. If NULL, uses all
    #   CensusIDs present after filtering.
    #
    # mortality_start
    #   Starting values for optim() in the mortality fit:
    #     c(log(h0_start), beta_start)
    #
    # meas_sd1_a, meas_sd1_b, meas_sd2, meas_p_big
    #   Parameters of the DBH measurement error model.
    #   We use the (linear-in-DBH) small-error SD:
    #     SD1(D) = a*D + b
    #   and a "large error" SD2 (cm) with mixture probability p_big.
    #
    #   These values follow the remeasurement-based error propagation discussion
    #   in the paper linked above (Royal Society Phil. Trans. B article).
    #   We keep the model in code form only (not reproducing paper text).
    #
    # use_measurement_error
    #   If TRUE:
    #     - subtract expected measurement variance when estimating growth SD
    #     - set shrinkage guardrails informed by measurement error
    #   If FALSE:
    #     - treat all growth variability as process variability
    #     - set shrinkage guardrails purely from data quantiles
    #
    # shrink_hard_prob
    #   Lower-tail probability used for the measurement-noise-only shrink quantile.
    #   Smaller values make max_shrink more permissive (more negative).
    #
    # shrink_data_quantile
    #   Lower quantile of observed annual increments used as a data-driven guardrail.
    #   Smaller values make max_shrink more permissive (more negative).
    #
    # recruit_max_quantile
    #   Upper quantile of observed recruits used as a guardrail for recruit size.
    #   Larger values make recruit_max_dbh more permissive.

    ### SETUP & LIBRARIES
    library(data.table)
    library(MASS)

    # =========================================================================
    # 1. INPUT CLEANING & VALIDATION
    # =========================================================================
    # We only use rows where:
    # - DBH is observed and positive
    # - TrueStemID is present (this function is estimating parameters from
    #   "known tracks" based on TrueStemID)
    #
    # This is a key point: estimate_bio_pars() uses TrueStemID to assemble
    # empirical growth/mortality/recruitment signals. If TrueStemID is missing
    # or unreliable, you should not trust the resulting parameters.
    x_dt <- as.data.table(x)
    x_dt <- x_dt[!is.na(DBH) & DBH > 0 & !is.na(TrueStemID)]
    # Ensure a `species` column exists for consistent downstream handling
    if (!("species" %in% names(x_dt))) {
        x_dt[, species := NA_character_]
    } else {
        x_dt[, species := as.character(species)]
    }

    x_dt[, Tag := as.character(Tag)]
    x_dt[, OriginalStemID := as.character(OriginalStemID)]
    x_dt[, TrueStemID := as.character(TrueStemID)]
    x_dt[, CensusID := as.integer(CensusID)]
    x_dt[, DBH := as.numeric(DBH)]
    x_dt[, ExactDate := as.IDate(ExactDate)]

    x_dt[, Tag_TrueStemID := paste(Tag, TrueStemID, sep = "_")]

    # =========================================================================
    # 2. DATE PROCESSING & COMPLETE-GRID CREATION
    # =========================================================================
    # mean date per Tag x Census
    mean_date_tag_census <- x_dt[
        , .(MeanExactDate = as.IDate(mean(as.numeric(ExactDate)))),
        by = .(Tag, CensusID)
    ]

    # mean date per Census
    mean_date_census <- x_dt[
        , .(MeanExactDate_Census = as.IDate(mean(as.numeric(ExactDate)))),
        by = CensusID
    ]

    ## ---- Complete grid (faster than expand.grid) --------------------------
    full_anchor_data <- CJ(
        Tag_TrueStemID = unique(x_dt$Tag_TrueStemID),
        CensusID       = unique(x_dt[CensusID >= anchor_start_census, CensusID])
    )

    # split IDs
    full_anchor_data[, c("Tag", "TrueStemID") := tstrsplit(Tag_TrueStemID, "_")]
    full_anchor_data[, Tag_TrueStemID := NULL]

    ## ---- Join dates (keyed joins = faster) ---------------------------------
    full_anchor_data <- mean_date_tag_census[
        full_anchor_data,
        on = .(Tag, CensusID)
    ]
    full_anchor_data <- mean_date_census[
        full_anchor_data,
        on = .(CensusID)
    ]
    # fill missing dates
    full_anchor_data[
        is.na(MeanExactDate),
        MeanExactDate := MeanExactDate_Census
    ][, MeanExactDate_Census := NULL]

    ## ---- Join anchor values -----------------------------------------------
    dt <- x_dt[, .(
        Tag,
        TrueStemID,
        CensusID,
        DBH,
        species,
        ExactDate
    )]

    if (nrow(dt) == 0) {
        stop("No usable rows after filtering.", call. = FALSE)
    }

    anchor_data_complete <- dt[
        full_anchor_data,
        on = .(Tag, TrueStemID, CensusID)
    ]

    anchor_data_complete <- anchor_data_complete[, ExactDate := fifelse(
        is.na(ExactDate),
        MeanExactDate,
        ExactDate
    )][, MeanExactDate := NULL]

    ## ---- Add mnemonic (single join) ---------------------------------------
    tag_mnemonic <- unique(dt[, .(Tag, species)])

    anchor_data_complete <- tag_mnemonic[
        anchor_data_complete,
        on = .(Tag)
    ]
    anchor_data_complete[, i.species := NULL]

    census_ids <- sort(unique(dt$CensusID))

    census_ids <- census_ids[is.finite(census_ids)]
    if (length(census_ids) < 2) {
        stop("Need at least two censuses.", call. = FALSE)
    }

    # =========================================================================
    # 3. WIDE-FORMAT CONVERSION
    # =========================================================================
    # Convert to one row per (Tag, TrueStemID, species) with one column per census.
    # This makes it easy to compute adjacent-census transitions.
    dw <- dcast(
        anchor_data_complete,
        Tag + TrueStemID + species ~ CensusID,
        value.var = "DBH"
    )

    iw <- dcast(
        anchor_data_complete,
        Tag + TrueStemID + species ~ CensusID,
        value.var = "ExactDate"
    )

    # =========================================================================
    # 4. GROWTH INCREMENT COLLECTION & DIAGNOSTICS
    # =========================================================================
    # For each adjacent census pair (t0, t1):
    #   d0 = DBH at t0 (cm)
    #   d1 = DBH at t1 (cm)
    #   g  = (d1 - d0) / T  (cm/year)  where T may be a scalar or vary per pair
    #
    # We collect:
    #   g_all           : all observed annualized increments
    #   d0_all, d1_all  : the corresponding DBHs for regression/diagnostics
    #   var_meas_g_all  : expected measurement-variance contribution to g
    #   T_all           : corresponding interval (years) used for each g
    g_all <- c()
    d0_all <- c()
    d1_all <- c()
    var_meas_g_all <- c()
    T_all <- c()
    # Diagnostic counters for interval handling (useful for testing)
    pairs_candidate_count <- 0L
    pairs_filled_with_scalar_count <- 0L
    pairs_dropped_count <- 0L

    sd1 <- function(d) {
        d <- as.numeric(d)
        pmax(meas_sd1_a * d + meas_sd1_b, 1e-6)
    }

    # Var(epsilon) under the mixture model.
    meas_var_eps <- function(d) {
        s1 <- sd1(d)
        (1 - meas_p_big) * (s1^2) + meas_p_big * (meas_sd2^2)
    }

    for (i in seq_len(length(census_ids) - 1)) {
        t0 <- as.character(census_ids[i])
        t1 <- as.character(census_ids[i + 1])

        if (!all(c(t0, t1) %in% names(dw))) next

        ok <- !is.na(dw[[t0]]) & !is.na(dw[[t1]])
        if (!any(ok)) next

        d0 <- dw[[t0]][ok]
        d1 <- dw[[t1]][ok]

        n_ok <- sum(ok)
        T1 <- if (t1 %in% names(iw)) iw[[t1]][ok] else rep(NA_real_, n_ok)
        T0 <- if (t0 %in% names(iw)) iw[[t0]][ok] else rep(NA_real_, n_ok)
        # Per-row coalesce: prefer T1, fall back to T0 when T1 is missing
        Tvec <- as.numeric(T1 - T0) / 365.25
        g <- (d1 - d0) / Tvec
        v_meas_g <- (meas_var_eps(d0) + meas_var_eps(d1)) / (Tvec^2)
        T_all <- c(T_all, Tvec)

        g_all <- c(g_all, g)
        d0_all <- c(d0_all, d0)
        d1_all <- c(d1_all, d1)
        var_meas_g_all <- c(var_meas_g_all, v_meas_g)
        pairs_candidate_count <- pairs_candidate_count + length(g)
    }

    if (length(g_all) < 5) {
        stop("Not enough growth observations.", call. = FALSE)
    }

    # =========================================================================
    # 5. GROWTH MODEL ESTIMATION — MEAN
    # =========================================================================
    # We fit a simple size-dependent mean model using the starting DBH (d0).
    # For robustness on small datasets, we fall back to a constant mean.
    mu_hat <- mean(g_all)

    ok_mu <- is.finite(g_all) & is.finite(d0_all) & (d0_all > 0)
    if (sum(ok_mu) >= 10L && stats::var(log(d0_all[ok_mu])) > 0) {
        fit_mu <- stats::lm(g_all[ok_mu] ~ log(d0_all[ok_mu]))
        alpha_hat <- as.numeric(stats::coef(fit_mu)[1])
        gamma_hat <- as.numeric(stats::coef(fit_mu)[2])
        if (!is.finite(alpha_hat)) alpha_hat <- mu_hat
        if (!is.finite(gamma_hat)) gamma_hat <- 0
    } else {
        alpha_hat <- mu_hat
        gamma_hat <- 0
    }

    mu_pred <- rep(alpha_hat, length(g_all))
    mu_pred[ok_mu] <- alpha_hat + gamma_hat * log(d0_all[ok_mu])

    # =========================================================================
    # 6. GROWTH MODEL ESTIMATION — VARIANCE
    # =========================================================================
    # We want the *process* SD of annual increments, not the observed SD which
    # includes measurement error.
    #
    # Target model:
    #   SD_process(g | d0) = sigma0 + sigma1*d0
    #
    # Observed increments include measurement noise, so we estimate a total SD
    # and then (optionally) subtract expected measurement variance in quadrature:
    #   SD_process^2 ≈ max( SD_total^2 - Var_meas(g), 0 )
    #
    # Robust SD estimation trick
    #   For X ~ Normal(0, sd^2), E|X| = sd*sqrt(2/pi).
    #   Rearranging yields an SD proxy:
    #     sd ≈ |X| * sqrt(pi/2)
    #   We apply this to residuals (g - mu_hat) to reduce sensitivity to outliers.
    resid_abs <- abs(g_all - mu_pred)
    # For Normal(0, sd^2), E|X| = sd*sqrt(2/pi), so sd_hat ≈ |X|*sqrt(pi/2)
    sd_total_hat <- resid_abs * sqrt(pi / 2)

    if (isTRUE(use_measurement_error)) {
        sd_proc_hat <- sqrt(pmax(sd_total_hat^2 - var_meas_g_all, 1e-8))
    } else {
        sd_proc_hat <- pmax(sd_total_hat, 1e-6)
    }

    # Fit a simple linear model for SD vs DBH.
    # Notes:
    # - sigma0_hat is constrained to be positive
    # - sigma1_hat is constrained to be non-negative
    fit_sd <- lm(sd_proc_hat ~ d0_all)
    sigma0_hat <- max(coef(fit_sd)[1], 1e-4)
    sigma1_hat <- max(coef(fit_sd)[2], 0)

    # =========================================================================
    # 7. SHRINKAGE PENALTY ESTIMATION (k_shrink)
    # =========================================================================
    # Shrinkage here means d1 < d0 (negative increment), which can occur due to:
    # - real biological shrinkage (limited)
    # - measurement error (common)
    # - ID swaps / bad matches (the thing we want to discourage)
    #
    # In transition_cost_tracks_bio(), shrinkage is penalized softly as:
    #   cost_shrink = k_shrink * (d0 - d1)^2
    # Units:
    #   (d0 - d1) is cm; to make cost dimensionless, k_shrink has units 1/cm^2.
    #
    # Heuristic used here
    #   If measurement error is enabled, set k_shrink so that shrinkage of about
    #   one typical measurement SD costs O(1).
    #   If measurement error is disabled, fall back to a crude estimate based on
    #   the variance of negative increments.
    # Estimate a soft shrink penalty from measurement error scale.
    # k_shrink is applied as:  k_shrink * (d0-d1)^2  (units: 1/cm^2).
    # With measurement noise, small shrinkage can be expected; we therefore scale
    # k_shrink so that shrinkage of ~1 SD costs O(1).
    k_shrink_source <- match.arg(k_shrink_source)
    k_shrink_fixed <- as.numeric(k_shrink_fixed)
    if (identical(k_shrink_source, "fixed")) {
        if (!is.finite(k_shrink_fixed) || k_shrink_fixed < 0) {
            stop("k_shrink_fixed must be a finite non-negative number when k_shrink_source='fixed'.", call. = FALSE)
        }
        k_shrink_hat_est <- NA_real_
        k_shrink_hat <- k_shrink_fixed
    } else {
        sd_meas_diff <- sqrt((meas_var_eps(d0_all) + meas_var_eps(d1_all)))
        sd_meas_diff <- sd_meas_diff[is.finite(sd_meas_diff) & sd_meas_diff > 0]
        if (isTRUE(use_measurement_error) && length(sd_meas_diff) >= 10) {
            s_typ <- stats::median(sd_meas_diff, na.rm = TRUE)
            k_shrink_hat_est <- 1 / (2 * (s_typ^2))
            k_shrink_hat_est <- min(max(k_shrink_hat_est, 1e-6), 1e6)
        } else {
            # Fallback (no measurement-error model): estimate from observed shrink magnitudes
            # on the *DBH difference* scale (cm), so k_shrink retains units 1/cm^2.
            delta_shrink <- (d0_all - d1_all)[is.finite(d0_all) & is.finite(d1_all) & (d1_all < d0_all)]
            if (length(delta_shrink) >= 5 && is.finite(stats::var(delta_shrink)) && stats::var(delta_shrink) > 0) {
                k_shrink_hat_est <- 1 / (2 * stats::var(delta_shrink))
            } else {
                k_shrink_hat_est <- 50
            }
        }
        k_shrink_hat <- k_shrink_hat_est
    }

    # =========================================================================
    # 8. MORTALITY MODEL ESTIMATION
    # =========================================================================
    # Mortality is inferred from TrueStemID tracks as:
    #   alive at t0 (DBH observed) and missing at t1 (DBH NA)
    #
    # We fit a simple hazard model:
    #   hazard(DBH) = h0 * exp(beta * DBH)
    #   P(death over interval T) = 1 - exp(-hazard(DBH)*T)
    #
    # mortality_start is on the unconstrained scale used by optim:
    #   par[1] = log(h0), par[2] = beta
    d0_m <- c()
    died <- c()
    T_m_vec <- c()

    for (i in seq_len(length(census_ids) - 1)) {
        t0 <- as.character(census_ids[i])
        t1 <- as.character(census_ids[i + 1])

        if (!all(c(t0, t1) %in% names(dw))) next

        at_risk <- !is.na(dw[[t0]])
        if (!any(at_risk)) next

        n_at <- sum(at_risk)
        T1m <- if (t1 %in% names(iw)) iw[[t1]][at_risk] else rep(NA_real_, n_at)
        T0m <- if (t0 %in% names(iw)) iw[[t0]][at_risk] else rep(NA_real_, n_at)
        Tm <- as.numeric(T1m - T0m) / 365.25

        d0_m <- c(d0_m, dw[[t0]][at_risk])
        died <- c(died, as.integer(is.na(dw[[t1]][at_risk])))
        T_m_vec <- c(T_m_vec, Tm)
        # Bookkeeping: how many mortality rows used and how many were filled/dropped
        pairs_candidate_count <- pairs_candidate_count + length(Tm)
        if (exists("T1m") && exists("T0m")) {
            pairs_filled_with_scalar_count <- pairs_filled_with_scalar_count + sum(!is.finite(T1m) & is.finite(T0m) & is.finite(Tm))
            pairs_dropped_count <- pairs_dropped_count + sum(!is.finite(T1m) & !is.finite(T0m) & !is.finite(Tm))
        }
    }

    negloglik_mort <- function(par, d0, died, Tvec) {
        h0 <- exp(par[1])
        beta <- par[2]

        hazard <- h0 * exp(beta * d0)
        p <- 1 - exp(-hazard * Tvec)
        p <- pmin(pmax(p, 1e-12), 1 - 1e-12)

        -sum(died * log(p) + (1 - died) * log(1 - p))
    }

    fit_m <- optim(
        mortality_start,
        negloglik_mort,
        d0 = d0_m,
        died = died,
        Tvec = T_m_vec,
        method = "BFGS"
    )

    h0_hat <- exp(fit_m$par[1])
    beta_hat <- fit_m$par[2]

    # =========================================================================
    # 9. RECRUITMENT MODEL ESTIMATION
    # =========================================================================
    # We identify "recruitment events" from TrueStemID tracks as:
    #   missing at t0 (DBH NA) and observed at t1 (DBH > 0)
    #
    # We estimate:
    #   - recruit sizes via a lognormal fit
    #   - recruit_max_dbh as a high quantile guardrail
    #   - recruit rate lambda as recruits per available NA slot per year
    recruit_dbh <- c()
    n_risk <- 0
    n_rec <- 0
    total_time_at_risk <- 0

    for (i in seq_len(length(census_ids) - 1)) {
        t0 <- as.character(census_ids[i])
        t1 <- as.character(census_ids[i + 1])

        if (!all(c(t0, t1) %in% names(dw))) next

        at_risk <- is.na(dw[[t0]])
        if (!any(at_risk)) next

        d1_at_risk <- dw[[t1]][at_risk]
        # Recruits must have a positive observed size. Non-positive values will
        # break the lognormal fit and are not meaningful DBH measurements.
        rec <- at_risk
        rec[at_risk] <- !is.na(d1_at_risk) & is.finite(d1_at_risk) & (d1_at_risk > 0)

        recruit_dbh <- c(recruit_dbh, dw[[t1]][rec])

        n_at <- sum(at_risk)
        T1r <- if (t1 %in% names(iw)) iw[[t1]][at_risk] else rep(NA_real_, n_at)
        T0r <- if (t0 %in% names(iw)) iw[[t0]][at_risk] else rep(NA_real_, n_at)
        T_r <- as.numeric(T1r - T0r) / 365.25

        total_time_at_risk <- total_time_at_risk + sum(T_r, na.rm = TRUE)
        # Also track any pair-level fills/drops as part of recruitment accounting
        if (exists("T1r") && exists("T0r")) {
            pairs_filled_with_scalar_count <- pairs_filled_with_scalar_count + sum(!is.finite(T1r) & is.finite(T0r))
            pairs_dropped_count <- pairs_dropped_count + sum(!is.finite(T1r) & !is.finite(T0r))
        }

        n_risk <- n_risk + sum(at_risk)
        n_rec <- n_rec + sum(rec)
    }
    recruit_dbh <- recruit_dbh[is.finite(recruit_dbh) & recruit_dbh > 0]
    # Bookkeeping for recruitment missing intervals
    # (we tracked contributions to total_time_at_risk above)

    if (length(recruit_dbh) >= 2) {
        fit_r <- fitdistr(recruit_dbh, "lognormal")
        mu_r <- fit_r$estimate["meanlog"]
        sd_r <- fit_r$estimate["sdlog"]
    } else {
        mu_r <- log(2)
        sd_r <- 0.5
    }

    # Guardrail: recruit_max_dbh can come from data or be fixed by the user.
    recruit_max_source <- match.arg(recruit_max_source)
    recruit_max_fixed <- as.numeric(recruit_max_fixed)
    if (identical(recruit_max_source, "fixed")) {
        if (!is.finite(recruit_max_fixed) || recruit_max_fixed <= 0) {
            stop("recruit_max_fixed must be a positive finite number when recruit_max_source='fixed'.", call. = FALSE)
        }
        # NOTE: recruit_max_dbh is on the DBH (cm) scale
        recruit_max_dbh <- recruit_max_fixed
    } else {
        recruit_max_dbh <- if (length(recruit_dbh) > 0) {
            as.numeric(stats::quantile(recruit_dbh, recruit_max_quantile, na.rm = TRUE))
        } else {
            5
        }
    }

    # Recruitment rate (Poisson)
    lambda_hat <- if (total_time_at_risk > 0) {
        n_rec / total_time_at_risk
    } else {
        0
    }

    # =========================================================================
    # 10. GUARDRAILS: SHRINKAGE (max_shrink)
    # =========================================================================
    # max_shrink is a conservative lower bound on annual growth (cm/year).
    # It is used as a guardrail to reject *absurd* shrinkage that is very
    # unlikely to be real or due to measurement error.
    #
    # Construction
    #   max_shrink_data : a small quantile of observed growth increments
    #   max_shrink_meas : a lower quantile of the measurement-noise-only
    #                    annualized difference distribution
    #   max_shrink_hat  : the minimum of these (more conservative)
    #
    # Measurement-noise-only lower quantile
    #   We compute a lower quantile of (e1 - e0)/T where e0 and e1 follow the
    #   mixture model. This yields a realistic lower tail for shrinkage driven
    #   purely by measurement error.
    #
    # Important: this is a guardrail, not a hard ecological law.
    # Measurement-informed lower bound on annual growth (mostly to prevent absurd matches).
    # We compute a conservative lower quantile of the *measurement-noise-only* annualized
    # DBH difference, using a typical diameter.

    long_time_interval <- c()

    for (i in seq_len(length(census_ids) - 1)) {
        t0 <- as.character(census_ids[i])
        t1 <- as.character(census_ids[i + 1])
        if (!all(c(t0, t1) %in% names(dw))) next
        T1mt <- iw[[t1]]
        T0mt <- iw[[t0]]
        long_time_interval <- as.numeric(T1mt - T0mt) / 365.25
    }

    time_interval_median <- median(long_time_interval, na.rm = TRUE)

    meas_lower_quantile_g <- function(p, d_typ, median_time_interval = time_interval_median) {
        p <- as.numeric(p)
        if (!is.finite(p) || p <= 0 || p >= 1) {
            return(NA_real_)
        }
        d_typ <- as.numeric(d_typ)
        if (!is.finite(d_typ) || d_typ <= 0) {
            return(NA_real_)
        }

        s_small <- sd1(d_typ)
        s_big <- meas_sd2
        w_small <- 1 - meas_p_big
        w_big <- meas_p_big

        # Four-component mixture for (e1-e0)/T
        sds <- c(
            sqrt(s_small^2 + s_small^2) / median_time_interval,
            sqrt(s_small^2 + s_big^2) / median_time_interval,
            sqrt(s_big^2 + s_small^2) / median_time_interval,
            sqrt(s_big^2 + s_big^2) / median_time_interval
        )
        wts <- c(w_small * w_small, w_small * w_big, w_big * w_small, w_big * w_big)

        cdf <- function(x) sum(wts * pnorm(x, mean = 0, sd = sds))

        lo <- -10
        hi <- 0
        # Expand if needed
        if (cdf(lo) > p) {
            lo2 <- -50
            if (cdf(lo2) > p) {
                return(lo2)
            }
            lo <- lo2
        }
        if (cdf(hi) < p) {
            hi <- 10
            if (cdf(hi) < p) {
                return(hi)
            }
        }

        out <- tryCatch(
            uniroot(function(x) cdf(x) - p, lower = lo, upper = hi, tol = 1e-6)$root,
            error = function(e) NA_real_
        )
        out
    }

    d_typ <- stats::median(d0_all[is.finite(d0_all) & d0_all > 0], na.rm = TRUE)
    max_shrink_meas <- if (isTRUE(use_measurement_error)) meas_lower_quantile_g(shrink_hard_prob, d_typ) else NA_real_
    max_shrink_data <- as.numeric(stats::quantile(g_all, shrink_data_quantile, na.rm = TRUE))
    max_shrink_hat_est <- if (isTRUE(use_measurement_error) && is.finite(max_shrink_meas)) {
        min(max_shrink_data, max_shrink_meas)
    } else {
        max_shrink_data
    }

    max_shrink_source <- match.arg(max_shrink_source)
    max_shrink_fixed <- as.numeric(max_shrink_fixed)
    if (identical(max_shrink_source, "fixed")) {
        if (!is.finite(max_shrink_fixed)) {
            stop("max_shrink_fixed must be a finite number when max_shrink_source='fixed'.", call. = FALSE)
        }
        max_shrink_hat <- max_shrink_fixed
    } else {
        max_shrink_hat <- max_shrink_hat_est
    }

    # -------------------------------------------------------------------------
    # GUARDRAILS: EXTREME GROWTH (upper tail)
    # -------------------------------------------------------------------------
    # We treat very large positive increments as a likely sign of an incorrect
    # match (ID swap, mis-measurement, etc.). To avoid the DP taking such edges,
    # we set a conservative (permissive) hard upper bound, and also provide a
    # softer threshold for a quadratic penalty.
    #
    # Construction mirrors the shrinkage guardrail, but for the upper tail:
    #   max_growth_data : upper quantile of observed annualized increments
    #   max_growth_meas : upper quantile of the measurement-noise-only mixture for (e1-e0)/T
    #   max_growth_hat  : max(max_growth_data, max_growth_meas)  (more permissive)
    meas_upper_quantile_g <- function(p, d_typ, median_time_interval = time_interval_median) {
        p <- as.numeric(p)
        if (!is.finite(p) || p <= 0 || p >= 1) {
            return(NA_real_)
        }
        d_typ <- as.numeric(d_typ)
        if (!is.finite(d_typ) || d_typ <= 0) {
            return(NA_real_)
        }

        s_small <- sd1(d_typ)
        s_big <- meas_sd2
        w_small <- 1 - meas_p_big
        w_big <- meas_p_big

        # Four-component mixture for (e1-e0)/T
        sds <- c(
            sqrt(s_small^2 + s_small^2) / median_time_interval,
            sqrt(s_small^2 + s_big^2) / median_time_interval,
            sqrt(s_big^2 + s_small^2) / median_time_interval,
            sqrt(s_big^2 + s_big^2) / median_time_interval
        )
        wts <- c(w_small * w_small, w_small * w_big, w_big * w_small, w_big * w_big)

        cdf <- function(x) sum(wts * pnorm(x, mean = 0, sd = sds))

        lo <- 0
        hi <- 10
        if (cdf(hi) < p) {
            hi2 <- 50
            if (cdf(hi2) < p) {
                return(hi2)
            }
            hi <- hi2
        }
        if (cdf(lo) > p) {
            lo2 <- -10
            if (cdf(lo2) > p) {
                return(lo2)
            }
            lo <- lo2
        }

        out <- tryCatch(
            uniroot(function(x) cdf(x) - p, lower = lo, upper = hi, tol = 1e-6)$root,
            error = function(e) NA_real_
        )
        out
    }

    growth_hard_prob <- as.numeric(growth_hard_prob)
    growth_data_quantile <- as.numeric(growth_data_quantile)
    growth_soft_quantile <- as.numeric(growth_soft_quantile)

    max_growth_data <- as.numeric(stats::quantile(g_all, growth_data_quantile, na.rm = TRUE))
    max_growth_meas <- if (isTRUE(use_measurement_error) && is.finite(d_typ)) {
        p_hi <- 1 - growth_hard_prob
        meas_upper_quantile_g(p_hi, d_typ)
    } else {
        NA_real_
    }
    max_growth_hat_est <- if (isTRUE(use_measurement_error) && is.finite(max_growth_meas)) {
        max(max_growth_data, max_growth_meas)
    } else {
        max_growth_data
    }

    max_growth_source <- match.arg(max_growth_source)
    max_growth_fixed <- as.numeric(max_growth_fixed)
    if (identical(max_growth_source, "fixed")) {
        if (!is.finite(max_growth_fixed) || max_growth_fixed <= 0) {
            stop("max_growth_fixed must be a finite positive number when max_growth_source='fixed'.", call. = FALSE)
        }
        max_growth_hat <- max_growth_fixed
    } else {
        max_growth_hat <- max_growth_hat_est
    }

    max_growth_soft_data <- as.numeric(stats::quantile(g_all, growth_soft_quantile, na.rm = TRUE))
    max_growth_soft_hat <- if (is.finite(max_growth_hat) && is.finite(max_growth_soft_data)) {
        min(max_growth_hat, max_growth_soft_data)
    } else {
        max_growth_soft_data
    }

    # Soft extreme-growth penalty strength (units: 1/cm^2), analogous to k_shrink.
    # In transition_cost_tracks_bio(), the penalty is applied to the *excess DBH* above
    # the soft cap in cm:
    #   excess = d1 - (d0 + max_growth_soft*T)
    #   cost_growth_soft = k_growth * excess^2
    # We therefore scale k_growth so that an excess of ~1 typical measurement SD costs O(1).
    sd_meas_diff2 <- sqrt((meas_var_eps(d0_all) + meas_var_eps(d1_all)))
    sd_meas_diff2 <- sd_meas_diff2[is.finite(sd_meas_diff2) & sd_meas_diff2 > 0]
    if (isTRUE(use_measurement_error) && length(sd_meas_diff2) >= 10) {
        s_typ2 <- stats::median(sd_meas_diff2, na.rm = TRUE)
        k_growth_hat_est <- 1 / (2 * (s_typ2^2))
        k_growth_hat_est <- min(max(k_growth_hat_est, 1e-6), 1e6)
    } else {
        # Fallback (no measurement-error model): estimate from observed DBH increments (cm)
        # so k_growth retains units 1/cm^2.
        delta_pos <- (d1_all - d0_all)[is.finite(d0_all) & is.finite(d1_all) & (d1_all > d0_all)]
        if (length(delta_pos) >= 5 && is.finite(stats::var(delta_pos)) && stats::var(delta_pos) > 0) {
            k_growth_hat_est <- 1 / (2 * stats::var(delta_pos))
        } else {
            k_growth_hat_est <- 50
        }
    }

    k_growth_source <- match.arg(k_growth_source)
    k_growth_fixed <- as.numeric(k_growth_fixed)
    if (identical(k_growth_source, "fixed")) {
        if (!is.finite(k_growth_fixed) || k_growth_fixed < 0) {
            stop("k_growth_fixed must be a finite non-negative number when k_growth_source='fixed'.", call. = FALSE)
        }
        k_growth_hat <- k_growth_fixed
    } else {
        k_growth_hat <- k_growth_hat_est
    }

    # =========================================================================
    # 11. RETURN PARAMETERS (structured list)
    # =========================================================================
    # The returned structure is intentionally aligned with:
    # - bio_pars_to_transition_args()
    # - transition_cost_tracks_bio()
    # - realism_calibration.R diagnostics
    list(
        growth = list(
            # Mean annual growth model: mu(DBH) = alpha + gamma*log(DBH)
            alpha = alpha_hat,
            gamma = gamma_hat,
            # Backward-compatible summary (empirical mean of g)
            mu = mu_hat,
            sigma0 = sigma0_hat,
            sigma1 = sigma1_hat,
            # Extreme-growth guardrails
            max_growth_soft = max_growth_soft_hat,
            max_growth = max_growth_hat,
            k_growth = k_growth_hat,
            k_growth_source = k_growth_source,
            k_growth_fixed = k_growth_fixed,
            k_growth_estimated_value = k_growth_hat_est,
            # Diagnostics / provenance (analogous to shrinkage)
            max_growth_data_quantile = growth_data_quantile,
            max_growth_data_value = max_growth_data,
            max_growth_meas_prob = growth_hard_prob,
            max_growth_meas_value = max_growth_meas,
            max_growth_soft_quantile = growth_soft_quantile,
            max_growth_soft_value = max_growth_soft_data,
            max_growth_source = max_growth_source,
            max_growth_fixed = max_growth_fixed,
            max_growth_estimated_value = max_growth_hat_est,
            # -----------------------------------------------------------------
            # Uniform nested layout (organizational; flat fields above remain)
            # -----------------------------------------------------------------
            guardrails = list(
                hard = list(
                    value = max_growth_hat,
                    source = max_growth_source,
                    fixed = max_growth_fixed,
                    estimated_value = max_growth_hat_est,
                    data_quantile = growth_data_quantile,
                    data_value = max_growth_data,
                    meas_prob = growth_hard_prob,
                    meas_value = max_growth_meas
                ),
                soft = list(
                    value = max_growth_soft_hat,
                    quantile = growth_soft_quantile,
                    data_value = max_growth_soft_data
                )
            ),
            penalties = list(
                soft = list(
                    k = k_growth_hat,
                    source = k_growth_source,
                    fixed = k_growth_fixed,
                    estimated_value = k_growth_hat_est
                )
            )
        ),
        mortality = list(
            h0   = h0_hat,
            beta = beta_hat
        ),
        recruitment = list(
            meanlog = mu_r,
            sdlog = sd_r,
            recruit_max_dbh = recruit_max_dbh,
            recruit_max_source = recruit_max_source,
            recruit_max_fixed = recruit_max_fixed,
            lambda = lambda_hat
        ),
        shrinkage = list(
            k_shrink = k_shrink_hat,
            max_shrink = max_shrink_hat,
            k_shrink_source = k_shrink_source,
            k_shrink_fixed = k_shrink_fixed,
            k_shrink_estimated_value = k_shrink_hat_est,
            max_shrink_source = max_shrink_source,
            max_shrink_fixed = max_shrink_fixed,
            max_shrink_estimated_value = max_shrink_hat_est,
            max_shrink_data_quantile = shrink_data_quantile,
            max_shrink_data_value = max_shrink_data,
            max_shrink_meas_prob = shrink_hard_prob,
            max_shrink_meas_value = max_shrink_meas,
            # -----------------------------------------------------------------
            # Uniform nested layout (organizational; flat fields above remain)
            # -----------------------------------------------------------------
            guardrails = list(
                hard = list(
                    value = max_shrink_hat,
                    source = max_shrink_source,
                    fixed = max_shrink_fixed,
                    estimated_value = max_shrink_hat_est,
                    data_quantile = shrink_data_quantile,
                    data_value = max_shrink_data,
                    meas_prob = shrink_hard_prob,
                    meas_value = max_shrink_meas
                ),
                # For shrinkage, the soft penalty begins as soon as DBH decreases (d1 < d0),
                # i.e., the soft "threshold" is effectively 0 cm shrink.
                soft = list(
                    value = 0
                )
            ),
            penalties = list(
                soft = list(
                    k = k_shrink_hat,
                    source = k_shrink_source,
                    fixed = k_shrink_fixed,
                    estimated_value = k_shrink_hat_est
                )
            )
        ),
        measurement_error = list(
            sd1_a = meas_sd1_a,
            sd1_b = meas_sd1_b,
            sd2 = meas_sd2,
            p_big = meas_p_big
        ),
        settings = list(
            use_measurement_error = use_measurement_error
        ),
        interval = list(
            # inferred_interval_years = interval_years,
            per_pair_intervals = T_all,
            pairs_candidate_count = as.integer(pairs_candidate_count),
            pairs_filled_with_scalar_count = as.integer(pairs_filled_with_scalar_count),
            pairs_dropped_count = as.integer(pairs_dropped_count)
        )
    )
}

transition_cost_tracks_bio_components <- function(
  track_dbh_t,
  track_dbh_tp1,
  interval_years,
  # -----------------------
  # GROWTH MODEL PARAMETERS
  # -----------------------
  mu_const,
  mu_gamma = 0,
  sigma0,
  sigma1,
  max_shrink,
  k_shrink,
  max_growth = Inf,
  max_growth_soft = Inf,
  k_growth = 0,
  # -------------------------
  # MORTALITY MODEL PARAMETERS
  # -------------------------
  h0,
  beta,
  # ----------------------------
  # RECRUITMENT MODEL PARAMETERS
  # ----------------------------
  recruit_meanlog,
  recruit_sdlog,
  recruit_max_dbh,
  recruit_lambda,
  # -----------------
  # MEASUREMENT ERROR (optional)
  # -----------------
  use_measurement_error = FALSE,
  meas_sd1_a = 0.0062,
  meas_sd1_b = 0.0904,
  meas_sd2 = 4.64,
  meas_p_big = 0.05,
  # -----------------
  # DETERMINISTIC TIE-BREAK
  # -----------------
  eps_tiebreak = 1e-6,
  hard_penalty = 1e6
) {
    # PURPOSE
    # - Diagnostics companion to transition_cost_tracks_bio().
    # - Returns a per-track breakdown of the *same* terms used in the scalar cost,
    #   plus the tie-break contribution.
    #
    # OUTPUT
    # - list(per_track=..., tiebreak=..., total=..., p_recruit=...)
    #   where per_track is a data.table with case labels and component costs.

    if (length(track_dbh_tp1) != length(track_dbh_t)) {
        stop("track_dbh_t and track_dbh_tp1 must have the same length.", call. = FALSE)
    }

    K <- length(track_dbh_t)
    interval_years <- as.numeric(interval_years)
    if (!is.finite(interval_years) || interval_years <= 0) {
        stop("interval_years must be positive.", call. = FALSE)
    }

    # Recruitment probability over interval
    p_recruit <- 1 - exp(-recruit_lambda * interval_years)
    p_recruit <- pmin(pmax(p_recruit, 1e-12), 1 - 1e-12)

    # Pre-allocate breakdown table
    if (!requireNamespace("data.table", quietly = TRUE)) {
        stop("Package 'data.table' is required for transition_cost_tracks_bio_components().")
    }

    dt <- data.table::data.table(
        track = seq_len(K),
        d0 = as.numeric(track_dbh_t),
        d1 = as.numeric(track_dbh_tp1),
        case = NA_character_,
        # Components (all non-negative; some may be 0)
        cost_recruit = 0,
        cost_no_recruit = 0,
        cost_mortality = 0,
        cost_growth_lik = 0,
        cost_shrink_soft = 0,
        cost_growth_soft = 0,
        cost_hard = 0,
        stringsAsFactors = FALSE
    )

    mu_growth <- function(d) {
        if (!is.finite(mu_gamma) || mu_gamma == 0 || !is.finite(d) || d <= 0) {
            return(mu_const)
        }
        mu_const + mu_gamma * log(d)
    }

    log_sum_exp <- function(x) {
        m <- max(x)
        if (!is.finite(m)) {
            return(m)
        }
        m + log(sum(exp(x - m)))
    }

    meas_sd1 <- function(d) pmax(meas_sd1_a * d + meas_sd1_b, 1e-6)

    for (k in seq_len(K)) {
        d0 <- dt$d0[k]
        d1 <- dt$d1[k]

        # CASE 1: NA -> NA
        if (is.na(d0) && is.na(d1)) {
            dt$case[k] <- "NA->NA"
            dt$cost_no_recruit[k] <- -log(1 - p_recruit)
            next
        }

        # CASE 2: NA -> DBH (recruitment)
        if (is.na(d0) && !is.na(d1)) {
            dt$case[k] <- "NA->DBH"
            if (d1 > recruit_max_dbh) {
                dt$cost_hard[k] <- hard_penalty
            } else {
                dt$cost_recruit[k] <- -log(p_recruit) - dlnorm(d1, recruit_meanlog, recruit_sdlog, log = TRUE)
            }
            next
        }

        # CASE 3: DBH -> NA (mortality)
        if (!is.na(d0) && is.na(d1)) {
            dt$case[k] <- "DBH->NA"
            hazard <- h0 * exp(beta * d0)
            p_death <- 1 - exp(-hazard * interval_years)
            p_death <- pmin(pmax(p_death, 1e-12), 1 - 1e-12)
            dt$cost_mortality[k] <- -log(p_death)
            next
        }

        # CASE 4: DBH -> DBH (growth)
        dt$case[k] <- "DBH->DBH"
        g <- (d1 - d0) / interval_years

        # Hard constraint on shrinkage
        if (is.finite(max_shrink) && (g < max_shrink)) {
            dt$cost_hard[k] <- hard_penalty
            next
        }

        # Hard constraint on extreme positive growth
        if (is.finite(max_growth) && (g > max_growth)) {
            dt$cost_hard[k] <- hard_penalty
            next
        }

        sigma_d <- sigma0 + sigma1 * d0
        sigma_d <- pmax(sigma_d, 1e-6)
        mu <- mu_growth(d0)

        if (isTRUE(use_measurement_error)) {
            s_small0 <- meas_sd1(d0)
            s_small1 <- meas_sd1(d1)
            s_big <- meas_sd2
            w_small <- 1 - meas_p_big
            w_big <- meas_p_big

            sd_meas_mix <- c(
                sqrt(s_small0^2 + s_small1^2) / interval_years,
                sqrt(s_small0^2 + s_big^2) / interval_years,
                sqrt(s_big^2 + s_small1^2) / interval_years,
                sqrt(s_big^2 + s_big^2) / interval_years
            )
            wt_meas_mix <- c(w_small * w_small, w_small * w_big, w_big * w_small, w_big * w_big)
            sd_tot <- sqrt(sigma_d^2 + sd_meas_mix^2)

            ll <- log(wt_meas_mix) + stats::dnorm(g, mean = mu, sd = sd_tot, log = TRUE)
            dt$cost_growth_lik[k] <- -log_sum_exp(ll)
        } else {
            dt$cost_growth_lik[k] <-
                (g - mu)^2 / (2 * sigma_d^2) +
                log(sigma_d) +
                0.5 * log(2 * pi)
        }

        if (d1 < d0) {
            dt$cost_shrink_soft[k] <- k_shrink * (d0 - d1)^2
        }

        if (is.finite(max_growth_soft) && is.finite(k_growth) && k_growth > 0) {
            d1_soft_cap <- d0 + max_growth_soft * interval_years
            if (is.finite(d1_soft_cap) && d1 > d1_soft_cap) {
                dt$cost_growth_soft[k] <- k_growth * (d1 - d1_soft_cap)^2
            }
        }
    }

    tie_cost <- 0
    if (eps_tiebreak > 0) {
        r0 <- rank(track_dbh_t, ties.method = "first")
        r1 <- rank(track_dbh_tp1, ties.method = "first")
        both_obs <- !is.na(track_dbh_t) & !is.na(track_dbh_tp1)
        if (any(both_obs)) {
            tie_cost <- eps_tiebreak * sum(abs(r0[both_obs] - r1[both_obs]))
        }
    }

    dt[, total_track := cost_recruit + cost_no_recruit + cost_mortality + cost_growth_lik + cost_shrink_soft + cost_growth_soft + cost_hard]

    list(
        per_track = dt,
        tiebreak = tie_cost,
        total = sum(dt$total_track) + tie_cost,
        p_recruit = p_recruit
    )
}
