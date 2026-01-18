# Utilities to inspect and visualize biological mortality parameters
#
# Functions:
# - mortality_params_df(bio_pars)
# - plot_mortality_params_bars(bio_pars, species = NULL, show = TRUE)
# - plot_mortality_hazard_curve(bio_pars, species = NULL, dbh_range = c(1,100), n = 200, log_y = FALSE, show = TRUE)
#
# Usage examples:
# bio <- ... # list of species -> list(mortality = list(h0 = ..., beta = ...), ...)
# print(mortality_params_df(bio))
# plot_mortality_params_bars(bio)
# plot_mortality_hazard_curve(bio, species = c("sp1", "sp3"))

mortality_params_df <- function(bio_pars) {
    if (!is.list(bio_pars) || length(bio_pars) == 0L) {
        stop("bio_pars must be a non-empty list of species parameter lists")
    }
    sp <- names(bio_pars)
    if (is.null(sp)) sp <- seq_along(bio_pars)

    df <- data.frame(
        species = character(0),
        h0 = numeric(0),
        beta = numeric(0),
        stringsAsFactors = FALSE
    )

    for (i in seq_along(bio_pars)) {
        b <- bio_pars[[i]]
        nm <- if (!is.null(names(bio_pars))) names(bio_pars)[[i]] else as.character(i)
        h0 <- NA_real_
        beta <- NA_real_
        if (!is.null(b) && is.list(b) && !is.null(b$mortality)) {
            h0 <- as.numeric(b$mortality$h0)
            beta <- as.numeric(b$mortality$beta)
        }
        df <- rbind(df, data.frame(species = nm, h0 = h0, beta = beta, stringsAsFactors = FALSE))
    }
    df
}

# Bar plot for h0 and beta across species
plot_mortality_params_bars <- function(bio_pars, species = NULL, show = TRUE) {
    df <- mortality_params_df(bio_pars)
    if (!is.null(species)) {
        df <- df[df$species %in% species, , drop = FALSE]
        if (nrow(df) == 0L) stop("No requested species found in bio_pars")
    }

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        warning("ggplot2 not available; falling back to base plotting. Install ggplot2 for nicer plots.")
        par(mfrow = c(1, 2))
        barplot(df$h0, names.arg = df$species, main = "h0 (mortality baseline)", ylab = "h0")
        barplot(df$beta, names.arg = df$species, main = "beta (mortality DBH coeff)", ylab = "beta")
        invisible(df)
    } else {
        library(ggplot2)
        library(reshape2)
        dfm <- reshape2::melt(df, id.vars = "species", measure.vars = c("h0", "beta"))
        p <- ggplot(dfm, aes(x = species, y = value, fill = variable)) +
            geom_col(position = position_dodge(width = 0.75)) +
            facet_wrap(~variable, scales = "free_y", nrow = 1) +
            theme_minimal() +
            xlab("Species") +
            ylab("Value") +
            ggtitle("Mortality parameters by species") +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))

        if (isTRUE(show)) print(p)
        invisible(p)
    }
}

# Plot hazard curve (hazard(DBH) = h0 * exp(beta * DBH)).
# Optionally show the interval mortality probability p = 1 - exp(-hazard * interval_years).
plot_mortality_hazard_curve <- function(bio_pars, species = NULL, dbh_range = c(1, 100), n = 200, interval_years = NULL, log_y = FALSE, show = TRUE) {
    df <- mortality_params_df(bio_pars)
    if (!is.null(species)) df <- df[df$species %in% species, , drop = FALSE]
    if (nrow(df) == 0L) stop("No species selected or species not found in bio_pars")

    dbh <- seq(dbh_range[1], dbh_range[2], length.out = n)
    out <- do.call(rbind, lapply(seq_len(nrow(df)), function(i) {
        sp <- df$species[i]
        h0 <- df$h0[i]
        beta <- df$beta[i]
        if (is.na(h0) || is.na(beta)) {
            return(NULL)
        }
        hazard <- h0 * exp(beta * dbh)
        if (!is.null(interval_years)) {
            prob <- 1 - exp(-hazard * interval_years)
            data.frame(species = sp, dbh = dbh, hazard = hazard, prob = prob, stringsAsFactors = FALSE)
        } else {
            data.frame(species = sp, dbh = dbh, hazard = hazard, stringsAsFactors = FALSE)
        }
    }))

    if (nrow(out) == 0L) stop("No valid mortality parameters (h0/beta) found for selected species")

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        warning("ggplot2 not available; falling back to base plotting. Install ggplot2 for nicer plots.")
        species_list <- unique(out$species)
        colors <- grDevices::rainbow(length(species_list))

        if (is.null(interval_years)) {
            plot(NULL, xlim = range(out$dbh), ylim = range(out$hazard), xlab = "DBH", ylab = "Hazard (annual)", main = "Mortality hazard vs DBH")
            for (i in seq_along(species_list)) {
                sp <- species_list[i]
                dsub <- out[out$species == sp, ]
                lines(dsub$dbh, dsub$hazard, col = colors[i], lwd = 2)
            }
            legend("topright", legend = species_list, col = colors, lwd = 2)
        } else {
            par(mfrow = c(2, 1), mar = c(4, 4, 2, 2))
            # hazard
            plot(NULL, xlim = range(out$dbh), ylim = range(out$hazard), xlab = "DBH", ylab = "Hazard (annual)", main = "Mortality hazard vs DBH")
            for (i in seq_along(species_list)) {
                sp <- species_list[i]
                dsub <- out[out$species == sp, ]
                lines(dsub$dbh, dsub$hazard, col = colors[i], lwd = 2)
            }
            legend("topright", legend = species_list, col = colors, lwd = 2)
            # interval probability
            plot(NULL, xlim = range(out$dbh), ylim = range(out$prob), xlab = "DBH", ylab = paste0("Probability over ", interval_years, " yr"), main = paste0(interval_years, "-year death probability"))
            for (i in seq_along(species_list)) {
                sp <- species_list[i]
                dsub <- out[out$species == sp, ]
                lines(dsub$dbh, dsub$prob, col = colors[i], lwd = 2)
            }
            legend("topright", legend = species_list, col = colors, lwd = 2)
        }

        invisible(out)
    } else {
        library(ggplot2)
        if (is.null(interval_years)) {
            p <- ggplot(out, aes(x = dbh, y = hazard, color = species)) +
                geom_line(linewidth = 1) +
                theme_minimal() +
                xlab("DBH") +
                ylab("Mortality hazard (annual)") +
                ggtitle("Mortality hazard vs DBH")

            if (isTRUE(log_y)) p <- p + scale_y_log10()

            if (isTRUE(show)) print(p)
            invisible(p)
        } else {
            library(reshape2)
            outm <- reshape2::melt(out, id.vars = c("species", "dbh"), measure.vars = c("hazard", "prob"), variable.name = "metric", value.name = "value")
            outm$metric_label <- ifelse(outm$metric == "hazard", "Hazard (annual)", paste0("Probability over ", interval_years, " yr"))
            p <- ggplot(outm, aes(x = dbh, y = value, color = species)) +
                geom_line(linewidth = 1) +
                facet_wrap(~metric_label, scales = "free_y", ncol = 1) +
                theme_minimal() +
                xlab("DBH") +
                ylab(NULL) +
                ggtitle(paste0("Mortality hazard and ", interval_years, "-year mortality probability")) +
                theme(strip.text = element_text(size = 11), legend.position = "right")

            if (isTRUE(log_y)) p <- p + scale_y_log10()

            if (isTRUE(show)) print(p)
            invisible(p)
        }
    }
}

# ---- Growth parameter helpers ----------------------------------------------

growth_params_df <- function(bio_pars) {
    if (!is.list(bio_pars) || length(bio_pars) == 0L) stop("bio_pars must be a non-empty list")
    sp <- names(bio_pars)
    if (is.null(sp)) sp <- seq_along(bio_pars)

    df <- data.frame(species = character(0), mu_const = numeric(0), gamma = numeric(0), sigma0 = numeric(0), sigma1 = numeric(0), max_growth = numeric(0), max_growth_soft = numeric(0), k_growth = numeric(0), stringsAsFactors = FALSE)

    for (i in seq_along(bio_pars)) {
        b <- bio_pars[[i]]
        nm <- if (!is.null(names(bio_pars))) names(bio_pars)[[i]] else as.character(i)
        mu_const <- NA_real_
        gamma <- NA_real_
        sigma0 <- NA_real_
        sigma1 <- NA_real_
        max_growth <- NA_real_
        max_growth_soft <- NA_real_
        k_growth <- NA_real_
        if (!is.null(b) && is.list(b) && !is.null(b$growth)) {
            g <- b$growth
            mu_const <- if (!is.null(g$alpha)) as.numeric(g$alpha) else as.numeric(g$mu)
            gamma <- if (!is.null(g$gamma)) as.numeric(g$gamma) else 0
            sigma0 <- if (!is.null(g$sigma0)) as.numeric(g$sigma0) else NA_real_
            sigma1 <- if (!is.null(g$sigma1)) as.numeric(g$sigma1) else NA_real_
            max_growth <- if (!is.null(g$max_growth)) as.numeric(g$max_growth) else NA_real_
            max_growth_soft <- if (!is.null(g$max_growth_soft)) as.numeric(g$max_growth_soft) else NA_real_
            k_growth <- if (!is.null(g$k_growth)) as.numeric(g$k_growth) else NA_real_
        }
        df <- rbind(df, data.frame(species = nm, mu_const = mu_const, gamma = gamma, sigma0 = sigma0, sigma1 = sigma1, max_growth = max_growth, max_growth_soft = max_growth_soft, k_growth = k_growth, stringsAsFactors = FALSE))
    }
    df
}

plot_growth_params_bars <- function(bio_pars, species = NULL, show = TRUE) {
    df <- growth_params_df(bio_pars)
    if (!is.null(species)) {
        df <- df[df$species %in% species, , drop = FALSE]
        if (nrow(df) == 0L) stop("No requested species found in bio_pars")
    }

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        warning("ggplot2 not available; falling back to base plotting. Install ggplot2 for nicer plots.")
        par(mfrow = c(2, 3))
        barplot(df$mu_const, names.arg = df$species, main = "mu (growth mean)", ylab = "mu")
        barplot(df$gamma, names.arg = df$species, main = "gamma (growth slope)", ylab = "gamma")
        barplot(df$sigma0, names.arg = df$species, main = "sigma0", ylab = "sigma0")
        barplot(df$sigma1, names.arg = df$species, main = "sigma1", ylab = "sigma1")
        barplot(df$max_growth, names.arg = df$species, main = "max_growth", ylab = "max_growth")
        barplot(df$k_growth, names.arg = df$species, main = "k_growth", ylab = "k_growth")
        invisible(df)
    } else {
        library(ggplot2)
        library(reshape2)
        dfm <- reshape2::melt(df, id.vars = "species", measure.vars = c("mu_const", "gamma", "sigma0", "sigma1", "max_growth", "k_growth"))
        p <- ggplot(dfm, aes(x = species, y = value, fill = variable)) +
            geom_col(position = position_dodge(width = 0.75)) +
            facet_wrap(~variable, scales = "free_y", nrow = 1) +
            theme_minimal() +
            xlab("Species") +
            ylab("Value") +
            ggtitle("Growth parameters by species") +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))

        if (isTRUE(show)) print(p)
        invisible(p)
    }
}

# Plot growth mean curve and growth SD: mu(DBH) = alpha + gamma * log(DBH), sigma(DBH) = sigma0 + sigma1 * DBH
plot_growth_mean_curve <- function(bio_pars, species = NULL, dbh_range = c(1, 100), n = 200, show = TRUE) {
    df <- growth_params_df(bio_pars)
    if (!is.null(species)) df <- df[df$species %in% species, , drop = FALSE]
    if (nrow(df) == 0L) stop("No species selected or species not found in bio_pars")

    dbh <- seq(dbh_range[1], dbh_range[2], length.out = n)
    out_mu <- do.call(rbind, lapply(seq_len(nrow(df)), function(i) {
        sp <- df$species[i]
        mu_const <- df$mu_const[i]
        gamma <- df$gamma[i]
        if (is.na(mu_const)) {
            return(NULL)
        }
        mu <- if (!is.na(gamma) && gamma != 0) mu_const + gamma * log(dbh) else rep(mu_const, length(dbh))
        data.frame(species = sp, dbh = dbh, metric = "mu", value = mu, stringsAsFactors = FALSE)
    }))
    out_sigma <- do.call(rbind, lapply(seq_len(nrow(df)), function(i) {
        sp <- df$species[i]
        s0 <- df$sigma0[i]
        s1 <- df$sigma1[i]
        if (is.na(s0) && is.na(s1)) {
            return(NULL)
        }
        sigma <- (if (is.na(s0)) 0 else s0) + (if (is.na(s1)) 0 else s1) * dbh
        data.frame(species = sp, dbh = dbh, metric = "sigma", value = sigma, stringsAsFactors = FALSE)
    }))

    out <- rbind(out_mu, out_sigma)
    if (nrow(out) == 0L) stop("No valid growth parameters found for selected species")

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        warning("ggplot2 not available; falling back to base plotting. Install ggplot2 for nicer plots.")
        species_list <- unique(out$species)
        colors <- grDevices::rainbow(length(species_list))
        par(mfrow = c(2, 1), mar = c(4, 4, 2, 2))
        # mu
        plot(NULL, xlim = range(out$dbh), ylim = range(out$value[out$metric == "mu"]), xlab = "DBH", ylab = "Mean growth (cm)", main = "Growth mean vs DBH")
        for (i in seq_along(species_list)) {
            sp <- species_list[i]
            dsub <- out[out$species == sp & out$metric == "mu", ]
            lines(dsub$dbh, dsub$value, col = colors[i], lwd = 2)
        }
        legend("topright", legend = species_list, col = colors, lwd = 2)
        # sigma
        plot(NULL, xlim = range(out$dbh), ylim = range(out$value[out$metric == "sigma"]), xlab = "DBH", ylab = "SD of growth (cm)", main = "Growth SD vs DBH")
        for (i in seq_along(species_list)) {
            sp <- species_list[i]
            dsub <- out[out$species == sp & out$metric == "sigma", ]
            lines(dsub$dbh, dsub$value, col = colors[i], lwd = 2)
        }
        legend("topright", legend = species_list, col = colors, lwd = 2)
        invisible(out)
    } else {
        library(ggplot2)
        p <- ggplot(out, aes(x = dbh, y = value, color = species)) +
            geom_line(linewidth = 1) +
            facet_wrap(~metric, scales = "free_y", ncol = 1, labeller = as_labeller(c(mu = "Mean growth (mu)", sigma = "Growth SD (sigma)"))) +
            theme_minimal() +
            xlab("DBH") +
            ylab(NULL) +
            ggtitle("Growth mean and SD vs DBH") +
            theme(strip.text = element_text(size = 11), legend.position = "right")

        if (isTRUE(show)) print(p)
        invisible(p)
    }
}

# ---- Recruitment parameter helpers -----------------------------------------

recruitment_params_df <- function(bio_pars) {
    if (!is.list(bio_pars) || length(bio_pars) == 0L) stop("bio_pars must be a non-empty list")
    sp <- names(bio_pars)
    if (is.null(sp)) sp <- seq_along(bio_pars)

    df <- data.frame(species = character(0), meanlog = numeric(0), sdlog = numeric(0), recruit_max_dbh = numeric(0), lambda = numeric(0), stringsAsFactors = FALSE)

    for (i in seq_along(bio_pars)) {
        b <- bio_pars[[i]]
        nm <- if (!is.null(names(bio_pars))) names(bio_pars)[[i]] else as.character(i)
        meanlog <- NA_real_
        sdlog <- NA_real_
        recruit_max_dbh <- NA_real_
        lambda <- NA_real_
        if (!is.null(b) && is.list(b) && !is.null(b$recruitment)) {
            r <- b$recruitment
            meanlog <- if (!is.null(r$meanlog)) as.numeric(r$meanlog) else NA_real_
            sdlog <- if (!is.null(r$sdlog)) as.numeric(r$sdlog) else NA_real_
            recruit_max_dbh <- if (!is.null(r$recruit_max_dbh)) as.numeric(r$recruit_max_dbh) else NA_real_
            lambda <- if (!is.null(r$lambda)) as.numeric(r$lambda) else NA_real_
        }
        df <- rbind(df, data.frame(species = nm, meanlog = meanlog, sdlog = sdlog, recruit_max_dbh = recruit_max_dbh, lambda = lambda, stringsAsFactors = FALSE))
    }
    df
}

plot_recruitment_params_bars <- function(bio_pars, species = NULL, show = TRUE) {
    df <- recruitment_params_df(bio_pars)
    if (!is.null(species)) {
        df <- df[df$species %in% species, , drop = FALSE]
        if (nrow(df) == 0L) stop("No requested species found in bio_pars")
    }

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        warning("ggplot2 not available; falling back to base plotting. Install ggplot2 for nicer plots.")
        par(mfrow = c(1, 3))
        barplot(df$meanlog, names.arg = df$species, main = "meanlog", ylab = "meanlog")
        barplot(df$sdlog, names.arg = df$species, main = "sdlog", ylab = "sdlog")
        barplot(df$recruit_max_dbh, names.arg = df$species, main = "recruit_max_dbh", ylab = "recruit_max_dbh")
        invisible(df)
    } else {
        library(ggplot2)
        library(reshape2)
        dfm <- reshape2::melt(df, id.vars = "species", measure.vars = c("meanlog", "sdlog", "recruit_max_dbh"))
        p <- ggplot(dfm, aes(x = species, y = value, fill = variable)) +
            geom_col(position = position_dodge(width = 0.75)) +
            facet_wrap(~variable, scales = "free_y", nrow = 1) +
            theme_minimal() +
            xlab("Species") +
            ylab("Value") +
            ggtitle("Recruitment parameters by species") +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))

        if (isTRUE(show)) print(p)
        invisible(p)
    }
}

# Plot recruit DBH PDF: lognormal with meanlog/sdlog
plot_recruitment_pdf_curve <- function(bio_pars, species = NULL, dbh_range = c(0.1, 50), n = 200, show = TRUE) {
    df <- recruitment_params_df(bio_pars)
    if (!is.null(species)) df <- df[df$species %in% species, , drop = FALSE]
    if (nrow(df) == 0L) stop("No species selected or species not found in bio_pars")

    dbh <- seq(dbh_range[1], dbh_range[2], length.out = n)
    out <- do.call(rbind, lapply(seq_len(nrow(df)), function(i) {
        sp <- df$species[i]
        meanlog <- df$meanlog[i]
        sdlog <- df$sdlog[i]
        recruit_max <- df$recruit_max_dbh[i]
        lambda <- df$lambda[i]
        if (is.na(meanlog) || is.na(sdlog)) {
            return(NULL)
        }
        dens <- dlnorm(dbh, meanlog, sdlog)
        data.frame(species = sp, dbh = dbh, density = dens, recruit_max_dbh = recruit_max, lambda = lambda, stringsAsFactors = FALSE)
    }))

    if (nrow(out) == 0L) stop("No valid recruitment parameters found for selected species")

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        warning("ggplot2 not available; falling back to base plotting. Install ggplot2 for nicer plots.")
        species_list <- unique(out$species)
        colors <- grDevices::rainbow(length(species_list))
        plot(NULL, xlim = range(out$dbh), ylim = range(out$density), xlab = "DBH", ylab = "Density", main = "Recruit DBH density (lognormal)")
        for (i in seq_along(species_list)) {
            sp <- species_list[i]
            dsub <- out[out$species == sp, ]
            lines(dsub$dbh, dsub$density, col = colors[i], lwd = 2)
            rm <- unique(dsub$recruit_max_dbh)
            if (!is.na(rm)) abline(v = rm, col = colors[i], lty = 2)
        }
        legend("topright", legend = species_list, col = colors, lwd = 2)
        invisible(out)
    } else {
        library(ggplot2)
        p <- ggplot(out, aes(x = dbh, y = density, color = species)) +
            geom_line(linewidth = 1) +
            theme_minimal() +
            xlab("DBH") +
            ylab("Density") +
            ggtitle("Recruit DBH density (lognormal)") +
            geom_vline(aes(xintercept = recruit_max_dbh, color = species), linetype = "dashed") +
            labs(subtitle = "dashed line = recruit_max_dbh; lambda (recruit rate) printed in legend")

        # Add lambda to legend by modifying species labels
        labs_df <- unique(out[, c("species", "lambda")])
        labs_df$label <- paste0(labs_df$species, " (lambda=", signif(labs_df$lambda, 3), ")")
        if (nrow(labs_df) > 0) {
            p <- p + scale_color_discrete(labels = labs_df$label)
        }

        if (isTRUE(show)) print(p)
        invisible(p)
    }
}

# Plot hazard and corresponding interval mortality probability together
plot_mortality_with_interval_prob <- function(bio_pars, species = NULL, dbh_range = c(1, 100), n = 200, interval_years = 5, log_y = FALSE, show = TRUE) {
    df <- mortality_params_df(bio_pars)
    if (!is.null(species)) df <- df[df$species %in% species, , drop = FALSE]
    if (nrow(df) == 0L) stop("No species selected or species not found in bio_pars")

    dbh <- seq(dbh_range[1], dbh_range[2], length.out = n)
    out <- do.call(rbind, lapply(seq_len(nrow(df)), function(i) {
        sp <- df$species[i]
        h0 <- df$h0[i]
        beta <- df$beta[i]
        if (is.na(h0) || is.na(beta)) {
            return(NULL)
        }
        hazard <- h0 * exp(beta * dbh)
        prob <- 1 - exp(-hazard * interval_years)
        data.frame(species = sp, dbh = dbh, hazard = hazard, prob = prob, stringsAsFactors = FALSE)
    }))

    if (nrow(out) == 0L) stop("No valid mortality parameters (h0/beta) found for selected species")

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        warning("ggplot2 not available; falling back to base plotting. Install ggplot2 for nicer plots.")
        species_list <- unique(out$species)
        # Two-panel base plot: top = hazard, bottom = probability
        par(mfrow = c(2, 1), mar = c(4, 4, 2, 2))
        plot(NULL, xlim = range(out$dbh), ylim = range(out$hazard), xlab = "DBH", ylab = "Hazard (annual)", main = "Mortality hazard")
        colors <- grDevices::rainbow(length(species_list))
        for (i in seq_along(species_list)) {
            sp <- species_list[i]
            dsub <- out[out$species == sp, ]
            lines(dsub$dbh, dsub$hazard, col = colors[i], lwd = 2)
        }
        legend("topright", legend = species_list, col = colors, lwd = 2)

        plot(NULL, xlim = range(out$dbh), ylim = range(out$prob), xlab = "DBH", ylab = paste0("Probability over ", interval_years, " yr"), main = paste0(interval_years, "-year death probability"))
        for (i in seq_along(species_list)) {
            sp <- species_list[i]
            dsub <- out[out$species == sp, ]
            lines(dsub$dbh, dsub$prob, col = colors[i], lwd = 2)
        }
        legend("topright", legend = species_list, col = colors, lwd = 2)

        invisible(out)
    } else {
        library(ggplot2)
        library(reshape2)
        outm <- reshape2::melt(out, id.vars = c("species", "dbh"), measure.vars = c("hazard", "prob"), variable.name = "metric", value.name = "value")
        # Label metrics nicely and facet into two rows
        outm$metric_label <- ifelse(outm$metric == "hazard", "Hazard (annual)", paste0("Probability over ", interval_years, " yr"))

        p <- ggplot(outm, aes(x = dbh, y = value, color = species)) +
            geom_line(linewidth = 1) +
            facet_wrap(~metric_label, scales = "free_y", ncol = 1) +
            theme_minimal() +
            xlab("DBH") +
            ylab(NULL) +
            ggtitle(paste0("Mortality hazard and ", interval_years, "-year mortality probability")) +
            theme(strip.text = element_text(size = 11), legend.position = "right")

        if (isTRUE(log_y)) p <- p + scale_y_log10()

        if (isTRUE(show)) print(p)
        invisible(p)
    }
}

# Export a multi-page PDF report with diagnostics for each species
# - Each plot call produces one page in the PDF (base or ggplot)
# - By default we include a small set of summary pages (parameter bars) and,
#   then, for each species: mortality (hazard + interval prob), growth (mu + sigma), recruitment (PDF).
if (!requireNamespace("here", quietly = TRUE)) {
    stop("Please install the 'here' package to use check_functions.r");
}

export_bio_pars_report <- function(
  bio_pars,
  species = NULL,
  out_file = here("dp_global", "output", "bio_pars_report.pdf"),
  dbh_ranges = list(mortality = c(1, 100), growth = c(1, 100), recruit = c(0.1, 50)),
  interval_years = 5,
  include_param_summary = TRUE,
  include_individual_pages = FALSE,
  width = 9,
  height = 5,
  open = FALSE
) {
    if (!is.list(bio_pars) || length(bio_pars) == 0L) stop("bio_pars must be a non-empty list")

    if (is.null(species)) {
        species <- names(bio_pars)
        if (is.null(species) || any(nzchar(species) == FALSE)) {
            species <- seq_along(bio_pars)
        }
    }
    species <- as.character(species)

    out_dir <- dirname(out_file)
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    pdf(out_file, width = width, height = height)
    on.exit(
        {
            try(dev.off(), silent = TRUE)
        },
        add = TRUE
    )

    # Summary pages (one page per summary figure)
    if (isTRUE(include_param_summary)) {
        try(plot_mortality_params_bars(bio_pars, species = NULL, show = TRUE), silent = TRUE)
        try(plot_growth_params_bars(bio_pars, species = NULL, show = TRUE), silent = TRUE)
        try(plot_recruitment_params_bars(bio_pars, species = NULL, show = TRUE), silent = TRUE)
    }

    # Combined curves: one page per type showing all requested species together
    try(
        {
            # Mortality combined: hazard + interval probability (if interval provided)
            if (!is.null(interval_years)) {
                plot_mortality_with_interval_prob(bio_pars, species = species, dbh_range = dbh_ranges$mortality, interval_years = interval_years, show = TRUE)
            } else {
                plot_mortality_hazard_curve(bio_pars, species = species, dbh_range = dbh_ranges$mortality, show = TRUE)
            }
        },
        silent = TRUE
    )

    try(
        {
            # Growth combined: mu and sigma panels
            plot_growth_mean_curve(bio_pars, species = species, dbh_range = dbh_ranges$growth, show = TRUE)
        },
        silent = TRUE
    )

    try(
        {
            # Recruitment combined: lognormal densities with recruit_max_dbh annotated
            plot_recruitment_pdf_curve(bio_pars, species = species, dbh_range = dbh_ranges$recruit, show = TRUE)
        },
        silent = TRUE
    )

    # Optional: include individual per-species pages (one species per trio of plots)
    if (isTRUE(include_individual_pages)) {
        for (sp in species) {
            try(plot_mortality_with_interval_prob(bio_pars, species = sp, dbh_range = dbh_ranges$mortality, interval_years = interval_years, show = TRUE), silent = TRUE)
            try(plot_growth_mean_curve(bio_pars, species = sp, dbh_range = dbh_ranges$growth, show = TRUE), silent = TRUE)
            try(plot_recruitment_pdf_curve(bio_pars, species = sp, dbh_range = dbh_ranges$recruit, show = TRUE), silent = TRUE)
        }
    }

    # Ensure device is closed before optionally opening file
    try(dev.off(), silent = TRUE)

    if (isTRUE(open)) {
        # macOS 'open' - best-effort
        try(system2("open", args = shQuote(normalizePath(out_file))), silent = TRUE)
    }

    invisible(normalizePath(out_file))
}

# End of file
