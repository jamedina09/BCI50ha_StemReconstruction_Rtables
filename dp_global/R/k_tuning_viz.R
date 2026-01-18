############################################################
### k_tuning_viz.R — visualize k_shrink/k_growth effects
############################################################
#
# Goal
# - Help interpret the *soft* penalty coefficients k_shrink and k_growth.
# - For a given (d0,d1,T) and a fitted bio parameter set, compare:
#     (A) JOIN: DBH->DBH (stay on same track)
#     (B) SPLIT: DBH->NA + NA->DBH (death + recruit)
# - Sweep k values and find the crossing point where JOIN becomes as costly as SPLIT.
#
# Notes
# - In the transition cost, soft penalties are applied on *DBH differences* (cm), not cm/year:
#     shrink: cost += k_shrink * (d0-d1)^2 only if d1 < d0
#     growth: cost += k_growth * (d1 - (d0 + max_growth_soft*T))^2 only if exceed cap
# - SPLIT does not depend on k_shrink/k_growth (it depends on mortality+recruitment terms).

k_threshold_from_costs <- function(join_base, split_cost, delta_cm) {
    join_base <- as.numeric(join_base)
    split_cost <- as.numeric(split_cost)
    delta_cm <- as.numeric(delta_cm)

    if (!is.finite(join_base) || !is.finite(split_cost) || !is.finite(delta_cm) || delta_cm <= 0) {
        return(NA_real_)
    }
    # solve: join_base + k * delta^2 == split_cost
    max(0, (split_cost - join_base) / (delta_cm^2))
}

k_sweep_join_vs_split <- function(
  scenarios,
  interval_years,
  bio,
  temperature = 1,
  k_grid = exp(seq(log(1e-3), log(1e3), length.out = 200)),
  which_k = c("auto", "shrink", "growth"),
  prune_min_annual_growth = NULL,
  prune_max_annual_growth = NULL,
  subtitle = NULL
) {
    if (!requireNamespace("data.table", quietly = TRUE)) {
        stop("Package 'data.table' is required for k_sweep_join_vs_split().")
    }
    which_k <- match.arg(which_k)

    if (is.null(scenarios) || nrow(as.data.frame(scenarios)) == 0L) {
        stop("scenarios must be a data.frame/data.table with columns d0 and d1 (and optional label).", call. = FALSE)
    }

    sc <- data.table::as.data.table(scenarios)
    if (!("d0" %in% names(sc)) || !("d1" %in% names(sc))) {
        stop("scenarios must contain columns 'd0' and 'd1'.", call. = FALSE)
    }
    if (!("label" %in% names(sc))) {
        sc[, label := paste0("d0=", d0, ", d1=", d1)]
    }

    interval_years <- as.numeric(interval_years)
    if (!is.finite(interval_years) || interval_years <= 0) {
        stop("interval_years must be positive.", call. = FALSE)
    }

    temperature <- as.numeric(temperature)
    if (!is.finite(temperature) || temperature <= 0) {
        stop("temperature must be positive.", call. = FALSE)
    }

    if (!exists("transition_cost_tracks_bio_components", mode = "function")) {
        stop("transition_cost_tracks_bio_components() not found. Did you source dp_global_biol.R?", call. = FALSE)
    }
    if (!exists("bio_pars_to_transition_args", mode = "function")) {
        stop("bio_pars_to_transition_args() not found. Did you source R/sensitivity_transition_cost_bio.R?", call. = FALSE)
    }

    base <- bio_pars_to_transition_args(bio)

    # Helper: compute join/split base costs with *soft penalties disabled*.
    base0 <- base
    base0$k_shrink <- 0
    base0$k_growth <- 0
    base0$eps_tiebreak <- 0

    # Candidate pruning bounds are DP-enumerator constraints (not part of cost).
    prune_min <- if (!is.null(prune_min_annual_growth)) as.numeric(prune_min_annual_growth) else NA_real_
    prune_max <- if (!is.null(prune_max_annual_growth)) as.numeric(prune_max_annual_growth) else NA_real_

    out_list <- vector("list", nrow(sc))

    for (i in seq_len(nrow(sc))) {
        d0 <- as.numeric(sc$d0[i])
        d1 <- as.numeric(sc$d1[i])
        label <- as.character(sc$label[i])

        join <- do.call(
            transition_cost_tracks_bio_components,
            c(list(track_dbh_t = c(d0), track_dbh_tp1 = c(d1), interval_years = interval_years), base0)
        )
        split_death <- do.call(
            transition_cost_tracks_bio_components,
            c(list(track_dbh_t = c(d0), track_dbh_tp1 = c(NA_real_), interval_years = interval_years), base0)
        )
        split_recruit <- do.call(
            transition_cost_tracks_bio_components,
            c(list(track_dbh_t = c(NA_real_), track_dbh_tp1 = c(d1), interval_years = interval_years), base0)
        )

        join_base <- as.numeric(join$total)
        split_cost <- as.numeric(split_death$total) + as.numeric(split_recruit$total)

        annual_growth <- if (is.finite(d0) && is.finite(d1)) (d1 - d0) / interval_years else NA_real_
        join_pruned <- is.finite(annual_growth) && (
            (is.finite(prune_min) && annual_growth < prune_min) ||
                (is.finite(prune_max) && annual_growth > prune_max)
        )

        # Which soft penalty is relevant for this (d0,d1)?
        delta_shrink <- if (is.finite(d0) && is.finite(d1)) max(0, d0 - d1) else NA_real_

        d1_soft_cap <- if (is.finite(base$max_growth_soft)) d0 + base$max_growth_soft * interval_years else Inf
        delta_growth <- if (is.finite(d0) && is.finite(d1) && is.finite(d1_soft_cap)) max(0, d1 - d1_soft_cap) else NA_real_

        mode <- which_k
        if (identical(mode, "auto")) {
            # Prefer shrink if there is shrink; else prefer growth if exceed soft cap.
            if (is.finite(delta_shrink) && delta_shrink > 0) {
                mode <- "shrink"
            } else if (is.finite(delta_growth) && delta_growth > 0) {
                mode <- "growth"
            } else {
                mode <- "shrink" # arbitrary; delta will be 0 => k has no effect
            }
        }

        delta_cm <- if (identical(mode, "shrink")) delta_shrink else delta_growth

        k_star <- k_threshold_from_costs(join_base = join_base, split_cost = split_cost, delta_cm = delta_cm)

        dt <- data.table::data.table(
            scenario = label,
            d0 = d0,
            d1 = d1,
            interval_years = interval_years,
            annual_growth = annual_growth,
            join_pruned = isTRUE(join_pruned),
            mode = mode,
            delta_cm = delta_cm,
            k = as.numeric(k_grid)
        )

        # JOIN cost under swept k.
        join_cost <- join_base
        if (is.finite(delta_cm) && delta_cm > 0) {
            join_cost <- join_cost + dt$k * (delta_cm^2)
        }

        dt[, `:=`(
            join_base = join_base,
            split_cost = split_cost,
            join_cost = join_cost,
            delta_join_minus_split = join_cost - split_cost,
            join_preferred = join_cost < split_cost,
            k_cross = k_star,
            # NOTE: exp(-x) can underflow to 0 for large x, which breaks log-scale plots.
            weight_ratio_join_over_split = exp(-(join_cost - split_cost) / temperature),
            log10_weight_ratio_join_over_split = -((join_cost - split_cost) / temperature) / log(10),
            temperature = temperature,
            recruit_max_dbh = as.numeric(base$recruit_max_dbh),
            max_growth_soft = as.numeric(base$max_growth_soft),
            max_shrink = as.numeric(base$max_shrink)
        )]

        out_list[[i]] <- dt
    }

    data.table::rbindlist(out_list, use.names = TRUE, fill = TRUE)
}

plot_k_sweep_join_vs_split <- function(
  dt,
  show_weight_ratio = TRUE,
  k_max = 1000,
  out_path = out_path, 
  subtitle = NULL
) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required for plot_k_sweep_join_vs_split(). Install ggplot2 or use the returned data.table.")
    }
    if (!requireNamespace("data.table", quietly = TRUE)) {
        stop("Package 'data.table' is required.")
    }

    d <- data.table::as.data.table(dt)
    d <- d[k <= k_max]

    # Cost panel
    d_cost <- data.table::copy(d)
    d_cost[, `:=`(
        join_line = join_cost,
        split_line = split_cost,
        scenario_label = paste0(
            scenario, "\n",
            "d0=", signif(d0, 4), ", d1=", signif(d1, 4),
            ", T=", interval_years, " yr\n",
            "growth/shrink=", sprintf("%.3f", annual_growth), " cm/yr"
        )
    )]

    p_cost <- ggplot2::ggplot(d_cost, ggplot2::aes(x = k)) +
        ggplot2::geom_line(ggplot2::aes(y = join_line, color = "join (DBH->DBH)"), linewidth = 0.8) +
        ggplot2::geom_line(ggplot2::aes(y = split_line, color = "split (DBH->NA + NA->DBH)"), linewidth = 0.8, linetype = "dashed") +
        ggplot2::scale_x_log10(
            breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000, 1e5),
            labels = scales::label_number(scale_cut = scales::cut_short_scale())
        ) +
        ggplot2::facet_wrap(~scenario_label, scales = "free_y") +
        ggplot2::labs(
            x = "k",
            y = "transition cost (negative log-likelihood)",
            color = NULL,
            title = "Soft-penalty sweep: when does JOIN lose to SPLIT?",
            subtitle = subtitle,
            caption = paste(
                "Solid colored curve: JOIN cost (DBH->DBH) as k increases.",
                "Dashed colored line: SPLIT cost (DBH->NA + NA->DBH), does not depend on k.",
                "Vertical gray line: k_cross (JOIN cost == SPLIT cost for this scenario)."
            )
        ) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
            legend.position = "bottom",
            plot.caption = ggplot2::element_text(hjust = 0),
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 11),
            axis.title.x = ggplot2::element_text(size = 12),
            axis.title.y = ggplot2::element_text(size = 12)
        )

    # Add vertical crosspoint line if finite
    k_cross_dt <- unique(d_cost[, .(scenario, k_cross, mode, delta_cm, join_pruned, scenario_label)])
    k_cross_dt <- k_cross_dt[is.finite(k_cross) & k_cross > 0]

    # for each facet in p_cost, add vertical line at k_cross
    if (nrow(k_cross_dt) > 0) {
        p_cost <- p_cost +
            ggplot2::geom_vline(
                data = k_cross_dt,
                mapping = ggplot2::aes(xintercept = k_cross, group = scenario_label), # ← key fix
                linewidth = 0.6,
                alpha = 0.6
            ) +
            ggplot2::geom_text(
                data = k_cross_dt,
                mapping = ggplot2::aes(
                    x = k_cross,
                    y = Inf,
                    label = paste0("k_cross=", signif(k_cross, 3)),
                    group = scenario_label # ← ensures correct panel
                ),
                angle = 90,
                vjust = 1.1,
                hjust = 1,
                size = 3.2,
                color = "gray30",
                inherit.aes = FALSE
            )
    }

    if (!isTRUE(show_weight_ratio)) {
        return(list(cost = p_cost))
    }

    # Weight ratio panel (marginal-DP interpretation).
    # Plot log10(weight ratio) to avoid underflow-to-zero issues.
    p_w <- ggplot2::ggplot(d, ggplot2::aes(x = k, y = log10_weight_ratio_join_over_split)) +
        ggplot2::geom_hline(yintercept = 0, linewidth = 0.6, linetype = "dotted") +
        ggplot2::geom_line(linewidth = 0.8) +
        ggplot2::scale_x_log10(
            breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000, 1e5),
            labels = scales::label_number(scale_cut = scales::cut_short_scale())
        ) +
        ggplot2::facet_wrap(~scenario) +
        ggplot2::labs(
            x = "k",
            y = "log10(join/split weight ratio)",
            title = "Marginal-DP interpretation: log10 relative weight of JOIN vs SPLIT",
            subtitle = subtitle,
            caption = paste(
                "Curve: log10( w_join / w_split ), where w ~ exp(-cost / temperature).",
                "Dotted line at 0 means equal weight (JOIN and SPLIT equally likely under the posterior). \n",
                ">0 favors JOIN; <0 favors SPLIT."
            )
        ) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
            plot.caption = ggplot2::element_text(hjust = 0),
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 11),
            axis.title.x = ggplot2::element_text(size = 12),
            axis.title.y = ggplot2::element_text(size = 12)
        )

    # Add k_cross vertical line to the weight panel too (for alignment).
    if (nrow(k_cross_dt) > 0) {
        p_w <- p_w +
            ggplot2::geom_vline(
                data = k_cross_dt,
                ggplot2::aes(xintercept = k_cross),
                linewidth = 0.5,
                alpha = 0.4
            )
    }

    list(cost = p_cost, weight_ratio = p_w)

    # export p_cost and p_w per page in a pdf
    if (!missing(out_path) && is.character(out_path) && nchar(out_path) > 0) {
        pdf(out_path, width = 14, height = 10)
        print(p_cost)
        print(p_w)
        dev.off()
    }
}

k_sweep_crosspoints <- function(dt) {
    if (!requireNamespace("data.table", quietly = TRUE)) {
        stop("Package 'data.table' is required.")
    }
    d <- data.table::as.data.table(dt)
    unique(d[, .(scenario, d0, d1, interval_years, mode, delta_cm, k_cross, join_pruned, join_base, split_cost)])
}
