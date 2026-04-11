############################################################
### basal_area_uncertainty.R — Individual-level basal area
###   with posterior uncertainty from stem identification
###
### Computes basal area (BA) and BA change at the individual (tag)
### level by tracking stems through posterior path samples from
### the DP / probabilistic reconstruction.
###
### For each tag (individual tree):
###   1. Per-stem BA = pi/4 * DBH^2 (cm^2), using reconstructed IDs
###   2. Tag-level total BA = sum of all stem BAs per census
###   3. BA change between consecutive censuses, decomposed into:
###      - Growth: BA change in surviving stems
###      - Loss:   BA removed by stem mortality
###      - Gain:   BA added by stem recruitment
###   4. Uncertainty quantified across posterior path samples
###
### Key insight: Tag-level total BA per census is INVARIANT to
### stem identity assignment (same DBH values sum identically).
### The decomposition into growth/loss/gain IS identity-dependent
### and is where posterior uncertainty manifests.
###
### Usage (from project root):
###   Rscript dp_global/scripts/basal_area_uncertainty.R \
###     --RUN_DIR=dp_global/output/<run_dir>
###
### Outputs (written to RUN_DIR/):
###   basal_area_tag_census.csv  — per-tag x census BA and stem count
###   basal_area_tag_change.csv  — per-tag BA change decomposition
###                                with posterior uncertainty
###   basal_area_figures.pdf     — multi-page diagnostic figures
############################################################

library(data.table)

# ---- 0. CLI and file discovery ------------------------------------------

cli_args <- commandArgs(trailingOnly = TRUE)
cli_run_dir <- NULL
for (a in cli_args) {
    m <- regmatches(a, regexec("^--RUN_DIR=(.+)$", a, ignore.case = TRUE))[[1]]
    if (length(m) == 2L) cli_run_dir <- m[2]
}
if (!is.null(cli_run_dir)) RUN_DIR <- cli_run_dir

if (!exists("RUN_DIR") || !nzchar(RUN_DIR)) {
    stop("RUN_DIR must be set (--RUN_DIR=<path> or assign before sourcing).")
}

if (!grepl("^/", RUN_DIR)) {
    root <- tryCatch(here::here(), error = function(e) getwd())
    RUN_DIR <- file.path(root, RUN_DIR)
}

recon_file <- file.path(RUN_DIR, "stem_reconstruction_dp_global_rcpp.csv")
post_dir <- file.path(RUN_DIR, "posteriors")

stopifnot(file.exists(recon_file))
has_posteriors <- dir.exists(post_dir) &&
    length(list.files(post_dir, pattern = "_paths\\.csv$")) > 0L

cat("[BA] Run directory :", RUN_DIR, "\n")
cat("[BA] Posteriors    :", if (has_posteriors) "found" else "none", "\n")

# ---- 1. Load reconstruction --------------------------------------------

rec <- fread(recon_file, na.strings = c("", "NA"))
cat("[BA] Loaded", nrow(rec), "rows,", uniqueN(rec$Tag), "tags\n")

census_dates <- rec[!is.na(ExactDate),
    .(MeanDate = mean(as.numeric(as.Date(ExactDate)), na.rm = TRUE)),
    by = CensusID
]
setkey(census_dates, CensusID)

ba_cm2 <- function(dbh) pi / 4 * dbh^2

# ---- 2. Map posterior files to tags -------------------------------------

parse_recon <- function(recon_str) {
    pairs <- strsplit(recon_str, ";", fixed = TRUE)[[1]]
    parts <- strsplit(pairs, ":", fixed = TRUE)
    data.table(
        obs_row_id = as.integer(vapply(parts, `[`, character(1), 1L)),
        StemID     = as.integer(vapply(parts, `[`, character(1), 2L))
    )
}

all_tags <- sort(unique(rec$Tag))
post_tag_map <- list()

if (has_posteriors) {
    post_files <- list.files(post_dir,
        pattern = "_paths\\.csv$",
        full.names = TRUE
    )
    for (pf in post_files) {
        tag_match <- regmatches(
            basename(pf),
            regexec("^tag_([^_]+)_posterior_samples", basename(pf))
        )[[1]]
        if (length(tag_match) < 2L) next
        tag_str <- tag_match[2]
        if (tag_str == "NA") {
            # Single-tag run: infer tag from reconstruction data
            if (length(all_tags) == 1L) {
                post_tag_map[[as.character(all_tags[1])]] <- pf
            }
        } else {
            post_tag_map[[tag_str]] <- pf
        }
    }
}

cat("[BA] Posterior files mapped:", length(post_tag_map), "tags\n")

# ---- 3. MAP tag-level BA per census -------------------------------------

map_stem_ba <- rec[!is.na(DBH) & !is.na(ReconstructedStemID),
    .(BA = ba_cm2(DBH)),
    by = .(Tag, CensusID, ReconstructedStemID)
]

tag_census <- map_stem_ba[, .(
    TotalBA_cm2 = sum(BA),
    NumStems    = .N
), by = .(Tag, CensusID)]

setorder(tag_census, Tag, CensusID)
tag_census <- merge(tag_census, census_dates, by = "CensusID", all.x = TRUE)
tag_census[, Year := 1970 + MeanDate / 365.25]

cat("[BA] Tag x census rows:", nrow(tag_census), "\n")

# ---- 4. BA change decomposition (MAP) ----------------------------------

decompose_ba_change <- function(stem_dt, census_pairs) {
    results <- vector("list", nrow(census_pairs))
    for (i in seq_len(nrow(census_pairs))) {
        cf <- census_pairs$c_from[i]
        ct <- census_pairs$c_to[i]

        sf <- stem_dt[CensusID == cf, .(StemID, BA_from = BA)]
        st <- stem_dt[CensusID == ct, .(StemID, BA_to = BA)]
        merged <- merge(sf, st, by = "StemID", all = TRUE)
        merged[, status := fifelse(
            !is.na(BA_from) & !is.na(BA_to), "survivor",
            fifelse(!is.na(BA_from), "death", "recruit")
        )]

        results[[i]] <- data.table(
            CensusID_from = cf,
            CensusID_to = ct,
            Growth_BA = sum(merged[status == "survivor", BA_to - BA_from],
                na.rm = TRUE
            ),
            Loss_BA = -sum(merged[status == "death", BA_from],
                na.rm = TRUE
            ),
            Gain_BA = sum(merged[status == "recruit", BA_to],
                na.rm = TRUE
            ),
            DeltaBA_total = sum(st$BA_to, na.rm = TRUE) -
                sum(sf$BA_from, na.rm = TRUE),
            NumSurvivors = sum(merged$status == "survivor"),
            NumDeaths = sum(merged$status == "death"),
            NumRecruits = sum(merged$status == "recruit")
        )
    }
    rbindlist(results)
}

map_change_list <- list()
for (tg in all_tags) {
    tag_rec <- rec[Tag == tg & !is.na(DBH) & !is.na(ReconstructedStemID)]
    censuses <- sort(unique(tag_rec$CensusID))
    if (length(censuses) < 2L) next

    census_pairs <- data.table(
        c_from = censuses[-length(censuses)],
        c_to   = censuses[-1]
    )

    stem_dt <- tag_rec[, .(
        StemID = ReconstructedStemID, CensusID,
        BA = ba_cm2(DBH)
    )]
    decomp <- decompose_ba_change(stem_dt, census_pairs)
    decomp[, Tag := tg]
    map_change_list[[length(map_change_list) + 1L]] <- decomp
}

map_change <- if (length(map_change_list) > 0L) {
    rbindlist(map_change_list, use.names = TRUE)
} else {
    data.table()
}

if (nrow(map_change) > 0L) {
    map_change <- merge(map_change,
        census_dates[, .(CensusID_from = CensusID, Date_from = MeanDate)],
        by = "CensusID_from", all.x = TRUE
    )
    map_change <- merge(map_change,
        census_dates[, .(CensusID_to = CensusID, Date_to = MeanDate)],
        by = "CensusID_to", all.x = TRUE
    )
    map_change[, Interval_yr := (Date_to - Date_from) / 365.25]
    map_change[is.na(Interval_yr) | Interval_yr <= 0, Interval_yr := 5.0]
}

cat(
    "[BA] MAP decomposition:", nrow(map_change), "intervals across",
    uniqueN(map_change$Tag), "tags\n"
)

# ---- 5. Posterior uncertainty in decomposition --------------------------

weighted_quantile <- function(vals, w, prob) {
    ord <- order(vals)
    cw <- cumsum(w[ord])
    vals[ord][which(cw >= prob)[1]]
}

post_change_list <- list()

for (tg_str in names(post_tag_map)) {
    tg <- tryCatch(as.integer(tg_str), warning = function(w) tg_str)
    pf <- post_tag_map[[tg_str]]

    post <- fread(pf)
    if (nrow(post) == 0L) next

    tag_rec <- rec[Tag == tg]
    if (nrow(tag_rec) == 0L) next

    obs_lookup <- tag_rec[, .(obs_row_id, CensusID, DBH)]
    censuses <- sort(unique(tag_rec$CensusID))
    if (length(censuses) < 2L) next

    census_pairs <- data.table(
        c_from = censuses[-length(censuses)],
        c_to   = censuses[-1]
    )

    n_paths <- nrow(post)
    cat(sprintf(
        "[BA] Tag %s: %d posterior paths, %d censuses\n",
        tg_str, n_paths, length(censuses)
    ))

    # Parse all paths into one table
    all_paths <- rbindlist(lapply(seq_len(n_paths), function(i) {
        dt <- parse_recon(post$recon[i])
        dt[, `:=`(path_idx = i, path_prob = post$path_prob[i])]
        dt
    }))
    all_paths <- merge(all_paths, obs_lookup, by = "obs_row_id", all.x = TRUE)
    all_paths[, BA := ba_cm2(DBH)]
    all_paths <- all_paths[!is.na(BA)]

    # Decompose per path (vectorised over census pairs)
    path_decomp_list <- vector("list", n_paths)
    for (ip in seq_len(n_paths)) {
        pdata <- all_paths[path_idx == ip]
        stem_dt <- pdata[, .(StemID, CensusID, BA)]
        decomp <- decompose_ba_change(stem_dt, census_pairs)
        decomp[, `:=`(path_idx = ip, path_prob = pdata$path_prob[1])]
        path_decomp_list[[ip]] <- decomp
    }

    path_decomp <- rbindlist(path_decomp_list, use.names = TRUE)
    path_decomp[, Tag := tg]

    # Weighted summary across paths
    post_summary <- path_decomp[,
        {
            w <- path_prob / sum(path_prob)
            list(
                Growth_mean = sum(w * Growth_BA),
                Growth_sd = sqrt(sum(w * (Growth_BA - sum(w * Growth_BA))^2)),
                Growth_q025 = weighted_quantile(Growth_BA, w, 0.025),
                Growth_q975 = weighted_quantile(Growth_BA, w, 0.975),
                Loss_mean = sum(w * Loss_BA),
                Loss_sd = sqrt(sum(w * (Loss_BA - sum(w * Loss_BA))^2)),
                Loss_q025 = weighted_quantile(Loss_BA, w, 0.025),
                Loss_q975 = weighted_quantile(Loss_BA, w, 0.975),
                Gain_mean = sum(w * Gain_BA),
                Gain_sd = sqrt(sum(w * (Gain_BA - sum(w * Gain_BA))^2)),
                Gain_q025 = weighted_quantile(Gain_BA, w, 0.025),
                Gain_q975 = weighted_quantile(Gain_BA, w, 0.975),
                DeltaBA_check = sum(w * DeltaBA_total),
                NumSurvivors_mean = sum(w * NumSurvivors),
                NumDeaths_mean = sum(w * NumDeaths),
                NumRecruits_mean = sum(w * NumRecruits),
                NumPaths = .N
            )
        },
        by = .(Tag, CensusID_from, CensusID_to)
    ]

    post_change_list[[length(post_change_list) + 1L]] <- post_summary
}

# ---- 6. Merge MAP + posterior and write CSVs ----------------------------

tag_change <- copy(map_change)
if (nrow(tag_change) > 0L) {
    tag_change[, c("Date_from", "Date_to") := NULL]
}

if (length(post_change_list) > 0L) {
    post_change <- rbindlist(post_change_list, use.names = TRUE, fill = TRUE)
    tag_change <- merge(tag_change, post_change,
        by = c("Tag", "CensusID_from", "CensusID_to"), all.x = TRUE
    )
}

setorder(tag_change, Tag, CensusID_from)

tag_census_out <- file.path(RUN_DIR, "basal_area_tag_census.csv")
fwrite(
    tag_census[, .(Tag, CensusID, Year, TotalBA_cm2, NumStems)],
    tag_census_out
)
cat("[BA] Wrote:", tag_census_out, "\n")

tag_change_out <- file.path(RUN_DIR, "basal_area_tag_change.csv")
fwrite(tag_change, tag_change_out)
cat("[BA] Wrote:", tag_change_out, "\n")

# ---- 7. Figures ---------------------------------------------------------

pdf_out <- file.path(RUN_DIR, "basal_area_figures.pdf")
pdf(pdf_out, width = 11, height = 8.5)

# Colour palette
COL_G <- "#2ca02c" # growth  — green
COL_L <- "#d62728" # loss    — red
COL_R <- "#1f77b4" # recruit — blue
COL_T <- "grey30" # total
COL_CI <- adjustcolor("grey50", alpha.f = 0.25)

# --- 7a. Summary page: all tags BA trajectory ----------------------------

if (length(all_tags) <= 30L) {
    # Small dataset: one subplot per tag
    nc <- ceiling(sqrt(length(all_tags)))
    nr <- ceiling(length(all_tags) / nc)
    par(
        mfrow = c(nr, nc), mar = c(3, 3.5, 2, 0.5), mgp = c(2, 0.6, 0),
        cex.main = 0.9, cex.lab = 0.8, cex.axis = 0.7
    )
    for (tg in all_tags) {
        tc <- tag_census[Tag == tg]
        if (nrow(tc) == 0L) next
        plot(tc$CensusID, tc$TotalBA_cm2,
            type = "b", pch = 19,
            col = COL_T, lwd = 1.5,
            xlab = "Census", ylab = expression("BA (cm"^2 * ")"),
            main = paste0("Tag ", tg), xaxt = "n"
        )
        axis(1, at = tc$CensusID, labels = tc$CensusID)
        grid(col = "grey90")
    }
} else {
    # Many tags: distribution summary
    par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1), mgp = c(2.5, 0.7, 0))

    # Distribution of total BA (latest census per tag)
    latest <- tag_census[, .SD[which.max(CensusID)], by = Tag]
    hist(latest$TotalBA_cm2,
        breaks = 30,
        col = adjustcolor("steelblue", 0.6), border = "white",
        main = "Total BA per Individual (Latest Census)",
        xlab = expression("BA (cm"^2 * ")"), ylab = "Number of Tags"
    )

    # Distribution of number of stems
    hist(latest$NumStems,
        breaks = max(latest$NumStems),
        col = adjustcolor("darkorange", 0.6), border = "white",
        main = "Stems per Individual (Latest Census)",
        xlab = "Number of Stems", ylab = "Number of Tags"
    )

    # Total BA over censuses (all tags superimposed)
    cids <- sort(unique(tag_census$CensusID))
    agg <- tag_census[, .(
        BA_mean = mean(TotalBA_cm2),
        BA_med = median(TotalBA_cm2),
        BA_q25 = quantile(TotalBA_cm2, 0.25),
        BA_q75 = quantile(TotalBA_cm2, 0.75)
    ),
    by = CensusID
    ]
    plot(agg$CensusID, agg$BA_med,
        type = "b", pch = 19,
        ylim = range(c(agg$BA_q25, agg$BA_q75)),
        xlab = "Census", ylab = expression("BA (cm"^2 * ")"),
        main = "Median BA per Individual Over Censuses", xaxt = "n"
    )
    axis(1, at = agg$CensusID)
    polygon(c(agg$CensusID, rev(agg$CensusID)),
        c(agg$BA_q25, rev(agg$BA_q75)),
        col = COL_CI, border = NA
    )
    lines(agg$CensusID, agg$BA_med, type = "b", pch = 19, lwd = 2)
    legend("topleft", "IQR", fill = COL_CI, border = NA, cex = 0.8)

    # Number of tags per census
    ntags <- tag_census[, .N, by = CensusID]
    barplot(ntags$N,
        names.arg = ntags$CensusID,
        col = "steelblue", border = NA,
        main = "Tags with BA Data per Census",
        xlab = "Census", ylab = "Number of Tags"
    )
}

# --- 7b. Per-tag detail pages (tags with posteriors only) ----------------

tags_with_post <- as.integer(names(post_tag_map))
tags_with_post <- tags_with_post[!is.na(tags_with_post)]

for (tg in sort(tags_with_post)) {
    tc <- tag_census[Tag == tg]
    ch <- tag_change[Tag == tg]
    if (nrow(tc) == 0L) next

    par(
        mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1), mgp = c(2.5, 0.7, 0),
        cex.main = 1.0, cex.lab = 0.9, cex.axis = 0.8
    )

    # Panel A: Total BA trajectory
    plot(tc$CensusID, tc$TotalBA_cm2,
        type = "b", pch = 19,
        col = COL_T, lwd = 2,
        xlab = "Census", ylab = expression("Total BA (cm"^2 * ")"),
        main = paste0("Tag ", tg, " -- Individual BA"),
        xaxt = "n"
    )
    axis(1, at = tc$CensusID)
    grid(col = "grey90")

    # Panel B: Number of stems
    plot(tc$CensusID, tc$NumStems,
        type = "s", lwd = 2, col = "steelblue",
        xlab = "Census", ylab = "Stems",
        main = paste0("Tag ", tg, " -- Stem Count"),
        xaxt = "n", ylim = c(0, max(tc$NumStems) * 1.2)
    )
    axis(1, at = tc$CensusID)
    points(tc$CensusID, tc$NumStems, pch = 19, col = "steelblue", cex = 1.2)
    grid(col = "grey90")

    if (nrow(ch) > 0L) {
        has_post <- "Growth_mean" %in% names(ch) && any(!is.na(ch$Growth_mean))

        growth_vals <- if (has_post) ch$Growth_mean else ch$Growth_BA
        loss_vals <- if (has_post) ch$Loss_mean else ch$Loss_BA
        gain_vals <- if (has_post) ch$Gain_mean else ch$Gain_BA

        # Panel C: BA decomposition bars
        n_int <- nrow(ch)
        intervals <- paste0("C", ch$CensusID_from, "->", ch$CensusID_to)
        xpos <- seq_len(n_int)
        bw <- 0.25

        all_vals <- c(growth_vals, loss_vals, gain_vals, ch$DeltaBA_total)
        if (has_post) {
            all_vals <- c(
                all_vals,
                ch$Growth_q025, ch$Growth_q975,
                ch$Loss_q025, ch$Loss_q975,
                ch$Gain_q025, ch$Gain_q975
            )
        }
        all_vals <- all_vals[is.finite(all_vals)]
        yr <- range(all_vals, na.rm = TRUE)
        yr <- yr + c(-1, 1) * diff(yr) * 0.15

        plot(NULL,
            xlim = c(0.3, n_int + 0.7), ylim = yr,
            xlab = "", ylab = expression(Delta * "BA (cm"^2 * ")"),
            main = paste0("Tag ", tg, " -- BA Change Decomposition"),
            xaxt = "n"
        )
        axis(1, at = xpos, labels = intervals, cex.axis = 0.75)
        abline(h = 0, col = "grey60", lty = 2)
        grid(col = "grey93", nx = NA, ny = NULL)

        rect(xpos - 1.5 * bw, 0, xpos - 0.5 * bw, growth_vals,
            col = COL_G, border = NA
        )
        rect(xpos - 0.5 * bw, 0, xpos + 0.5 * bw, loss_vals,
            col = COL_L, border = NA
        )
        rect(xpos + 0.5 * bw, 0, xpos + 1.5 * bw, gain_vals,
            col = COL_R, border = NA
        )

        if (has_post) {
            arrows(xpos - bw, ch$Growth_q025, xpos - bw, ch$Growth_q975,
                angle = 90, code = 3, length = 0.04, col = COL_G, lwd = 1.5
            )
            arrows(xpos, ch$Loss_q025, xpos, ch$Loss_q975,
                angle = 90, code = 3, length = 0.04, col = COL_L, lwd = 1.5
            )
            arrows(xpos + bw, ch$Gain_q025, xpos + bw, ch$Gain_q975,
                angle = 90, code = 3, length = 0.04, col = COL_R, lwd = 1.5
            )
        }

        points(xpos, ch$DeltaBA_total, pch = 18, col = COL_T, cex = 1.5)
        if (n_int > 1L) {
            lines(xpos, ch$DeltaBA_total, col = COL_T, lwd = 1.5, lty = 2)
        }

        legend("topright",
            legend = c(
                "Survivor growth", "Mortality loss",
                "Recruitment gain", expression(Delta * "BA total")
            ),
            fill = c(COL_G, COL_L, COL_R, NA),
            border = NA,
            pch = c(NA, NA, NA, 18),
            lty = c(NA, NA, NA, 2),
            col = c(NA, NA, NA, COL_T),
            cex = 0.75, bg = "white"
        )

        # Panel D: Stem demographics
        n_surv <- if (has_post) ch$NumSurvivors_mean else ch$NumSurvivors
        n_dead <- if (has_post) ch$NumDeaths_mean else ch$NumDeaths
        n_recr <- if (has_post) ch$NumRecruits_mean else ch$NumRecruits

        cmat <- rbind(Survivors = n_surv, Deaths = n_dead, Recruits = n_recr)
        colnames(cmat) <- intervals

        barplot(cmat,
            beside = TRUE,
            col = c(COL_G, COL_L, COL_R), border = NA,
            main = paste0("Tag ", tg, " -- Stem Demographics"),
            xlab = "", ylab = "Number of Stems",
            legend.text = TRUE,
            args.legend = list(cex = 0.75, bg = "white")
        )
    } else {
        plot.new()
        text(0.5, 0.5, "Single census", cex = 1.2)
        plot.new()
        text(0.5, 0.5, "Single census", cex = 1.2)
    }
}

# --- 7c. Uncertainty summary page (across all tags with posteriors) ------

if (length(post_change_list) > 0L) {
    post_all <- rbindlist(post_change_list, use.names = TRUE, fill = TRUE)

    par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1), mgp = c(2.5, 0.7, 0))

    nonzero_g <- post_all$Growth_sd[post_all$Growth_sd > 1e-8]
    nonzero_l <- post_all$Loss_sd[post_all$Loss_sd > 1e-8]
    nonzero_r <- post_all$Gain_sd[post_all$Gain_sd > 1e-8]

    if (length(nonzero_g) > 0L) {
        hist(nonzero_g,
            breaks = 25,
            col = adjustcolor(COL_G, 0.5), border = "white",
            main = "Growth Component -- Posterior SD",
            xlab = expression("SD (cm"^2 * ")"), ylab = "Intervals"
        )
    } else {
        plot.new()
        text(0.5, 0.5, "No growth uncertainty", cex = 1.1)
    }

    if (length(nonzero_l) > 0L) {
        hist(nonzero_l,
            breaks = 25,
            col = adjustcolor(COL_L, 0.5), border = "white",
            main = "Loss Component -- Posterior SD",
            xlab = expression("SD (cm"^2 * ")"), ylab = "Intervals"
        )
    } else {
        plot.new()
        text(0.5, 0.5, "No loss uncertainty", cex = 1.1)
    }

    if (length(nonzero_r) > 0L) {
        hist(nonzero_r,
            breaks = 25,
            col = adjustcolor(COL_R, 0.5), border = "white",
            main = "Recruitment Component -- Posterior SD",
            xlab = expression("SD (cm"^2 * ")"), ylab = "Intervals"
        )
    } else {
        plot.new()
        text(0.5, 0.5, "No recruitment uncertainty", cex = 1.1)
    }

    # CI widths for growth
    ci_w <- abs(post_all$Growth_q975 - post_all$Growth_q025)
    ci_w <- ci_w[is.finite(ci_w) & ci_w > 1e-8]
    if (length(ci_w) > 0L) {
        hist(ci_w,
            breaks = 25,
            col = adjustcolor("grey50", 0.5), border = "white",
            main = "Growth Component -- 95% CI Width",
            xlab = expression("CI Width (cm"^2 * ")"), ylab = "Intervals"
        )
    } else {
        plot.new()
        text(0.5, 0.5, "No CI width variation", cex = 1.1)
    }
}

invisible(dev.off())
cat("[BA] Wrote figures:", pdf_out, "\n")

# ---- 8. Summary statistics ----------------------------------------------

cat("\n=========== BASAL AREA SUMMARY ===========\n\n")
cat(sprintf("Tags total: %d\n", length(all_tags)))
cat(sprintf("Tags with posteriors: %d\n", length(post_tag_map)))
cat(sprintf("Census intervals: %d\n", nrow(tag_change)))

if (nrow(tag_change) > 0L && "Growth_sd" %in% names(tag_change)) {
    uncertain <- tag_change[!is.na(Growth_sd) & Growth_sd > 1e-8]

    cat(sprintf(
        "\nIntervals with decomposition uncertainty: %d / %d (%.1f%%)\n",
        nrow(uncertain), nrow(tag_change),
        100 * nrow(uncertain) / max(nrow(tag_change), 1)
    ))

    if (nrow(uncertain) > 0L) {
        cat("\n  Growth component:\n")
        cat(sprintf("    Mean SD  : %.3f cm^2\n", mean(uncertain$Growth_sd)))
        cat(sprintf("    Max  SD  : %.3f cm^2\n", max(uncertain$Growth_sd)))
        cat(sprintf(
            "    Mean 95%% CI width: %.3f cm^2\n",
            mean(abs(uncertain$Growth_q975 - uncertain$Growth_q025),
                na.rm = TRUE
            )
        ))

        cat("\n  Loss component:\n")
        cat(sprintf(
            "    Mean SD  : %.3f cm^2\n",
            mean(uncertain$Loss_sd, na.rm = TRUE)
        ))
        cat(sprintf(
            "    Max  SD  : %.3f cm^2\n",
            max(uncertain$Loss_sd, na.rm = TRUE)
        ))

        cat("\n  Recruitment component:\n")
        cat(sprintf(
            "    Mean SD  : %.3f cm^2\n",
            mean(uncertain$Gain_sd, na.rm = TRUE)
        ))
        cat(sprintf(
            "    Max  SD  : %.3f cm^2\n",
            max(uncertain$Gain_sd, na.rm = TRUE)
        ))

        cat("\n  Per-tag detail:\n")
        for (tg in sort(unique(uncertain$Tag))) {
            td <- uncertain[Tag == tg]
            cat(sprintf(
                "    Tag %s: %d intervals, Growth SD [%.3f, %.3f],",
                tg, nrow(td), min(td$Growth_sd), max(td$Growth_sd)
            ))
            cat(sprintf(
                " Loss SD [%.3f, %.3f],",
                min(td$Loss_sd, na.rm = TRUE), max(td$Loss_sd, na.rm = TRUE)
            ))
            cat(sprintf(
                " Gain SD [%.3f, %.3f]\n",
                min(td$Gain_sd, na.rm = TRUE), max(td$Gain_sd, na.rm = TRUE)
            ))
        }
    } else {
        cat("\nNo decomposition uncertainty detected (all paths agree).\n")
    }
}

cat("\nDone.\n")

# Rscript dp_global/scripts/basal_area_uncertainty.R --RUN_DIR="dp_global/output/20260411_133614_BCI_tag171486_T171486_DP_MB_ME_g5_sm0p5_kg0_ks0_rcpp"