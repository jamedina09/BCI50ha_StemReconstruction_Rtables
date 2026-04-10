############################################################
### basal_area_uncertainty.R — Posterior uncertainty in basal area
############################################################
# Quantify uncertainty in basal area estimates using posterior
# path samples from the DP / probabilistic reconstruction.
#
# For each tag × census, the total basal area (sum of pi/4 * DBH^2
# across living stems) is invariant to identity assignment.  The
# uncertainty that identity reconstruction introduces is in:
#
#   (a) Per-stem BA trajectories (which DBH values belong to which stem)
#   (b) BA growth rates derived from those trajectories
#   (c) Demographic rates (mortality / recruitment) for BA accounting
#
# This script computes (a) and (b) across posterior path samples and
# summarises the resulting distributions.
#
# Usage (from project root):
#   Rscript dp_global/scripts/basal_area_uncertainty.R \
#     --RUN_DIR=dp_global/output/<run_dir>
#
# Or source interactively and set RUN_DIR before sourcing.
#
# Outputs (written to RUN_DIR/):
#   basal_area_uncertainty_tag.csv   — per-tag × census BA summary
#   basal_area_uncertainty_stem.csv  — per-stem × census BA posterior
#   basal_area_uncertainty_growth.csv — per-stem BA growth rate posterior
############################################################

library(data.table)

# ---- 0) Locate run directory -------------------------------------------

# Accept CLI --RUN_DIR=... or use the variable if already set
cli_args <- commandArgs(trailingOnly = TRUE)
cli_run_dir <- NULL
for (a in cli_args) {
    m <- regmatches(a, regexec("^--RUN_DIR=(.+)$", a, ignore.case = TRUE))[[1]]
    if (length(m) == 2L) cli_run_dir <- m[2]
}
if (!is.null(cli_run_dir)) RUN_DIR <- cli_run_dir

if (!exists("RUN_DIR") || !nzchar(RUN_DIR))
    stop("RUN_DIR must be set (--RUN_DIR=<path> or assign before sourcing).")

# Resolve relative paths against workspace root
if (!grepl("^/", RUN_DIR)) {
    root <- tryCatch(here::here(), error = function(e) getwd())
    RUN_DIR <- file.path(root, RUN_DIR)
}

recon_file <- file.path(RUN_DIR, "stem_reconstruction_dp_global_rcpp.csv")
post_dir   <- file.path(RUN_DIR, "posteriors")

stopifnot(
    file.exists(recon_file),
    dir.exists(post_dir)
)

cat("[BA uncertainty] Run directory:", RUN_DIR, "\n")

# ---- 1) Load reconstruction -------------------------------------------

rec <- fread(recon_file)
cat("[BA uncertainty] Loaded", nrow(rec), "rows,",
    length(unique(rec$Tag)), "tags\n")

# Compute per-census mean dates for interval calculation
census_dates <- rec[!is.na(ExactDate),
    .(MeanDate = mean(as.numeric(as.Date(ExactDate)), na.rm = TRUE)),
    by = CensusID]
setkey(census_dates, CensusID)

# ---- 2) Parse posterior paths ------------------------------------------

post_files <- list.files(post_dir, pattern = "_paths\\.csv$", full.names = TRUE)
cat("[BA uncertainty] Found", length(post_files), "posterior files\n")

# Parse a recon string "1:3;4:5;..." into a data.table(obs_row_id, stemid)
parse_recon <- function(recon_str) {
    pairs <- strsplit(recon_str, ";", fixed = TRUE)[[1]]
    parts <- strsplit(pairs, ":", fixed = TRUE)
    data.table(
        obs_row_id = as.integer(vapply(parts, `[`, character(1), 1L)),
        StemID     = as.integer(vapply(parts, `[`, character(1), 2L))
    )
}

# ---- 3) Basal area helper (cm^2) ----------------------------------------
#   BA = pi/4 * DBH^2   (DBH in cm => BA in cm^2)
ba_cm2 <- function(dbh) pi / 4 * dbh^2

# ---- 4) Process each tag ------------------------------------------------

tag_results  <- list()
stem_results <- list()
growth_results <- list()

tags_with_posteriors <- integer(0)

for (pf in post_files) {
    # Extract tag from filename: tag_<id>_posterior_samples__paths.csv
    tag_match <- regmatches(basename(pf),
        regexec("^tag_([^_]+)_posterior_samples__paths\\.csv$", basename(pf)))[[1]]
    if (length(tag_match) < 2L) next
    tag_id <- as.integer(tag_match[2])

    post <- fread(pf)
    if (nrow(post) == 0L) next

    # Reconstruction rows for this tag
    tag_rec <- rec[Tag == tag_id]
    if (nrow(tag_rec) == 0L) next

    tags_with_posteriors <- c(tags_with_posteriors, tag_id)

    # Lookup table: obs_row_id -> (CensusID, DBH)
    obs_lookup <- tag_rec[, .(obs_row_id, CensusID, DBH, OriginalStemID)]
    setkey(obs_lookup, obs_row_id)

    censuses <- sort(unique(tag_rec$CensusID))
    n_paths  <- nrow(post)

    # ---- 4a) Tag-level BA per census per posterior path ------------------
    # Note: tag-level BA is constant across paths (same DBH values summed)
    # but per-STEM BA changes.  We compute both for verification.

    # Pre-compute MAP (main reconstruction) per-stem BA trajectories
    map_obs <- tag_rec[!is.na(DBH) & !is.na(ReconstructedStemID),
        .(obs_row_id, CensusID, DBH, ReconstructedStemID)]

    # Per-posterior-path stem-level BA
    path_stem_ba_list <- vector("list", n_paths)

    for (ip in seq_len(n_paths)) {
        path_dt <- parse_recon(post$recon[ip])
        # Merge with obs_lookup to get CensusID and DBH
        path_dt <- merge(path_dt, obs_lookup, by = "obs_row_id", all.x = TRUE)
        path_dt <- path_dt[!is.na(DBH)]
        if (nrow(path_dt) == 0L) next

        path_dt[, BA := ba_cm2(DBH)]
        path_dt[, path_prob := post$path_prob[ip]]
        path_dt[, path_idx := ip]

        path_stem_ba_list[[ip]] <- path_dt[, .(
            Tag = tag_id,
            CensusID,
            StemID,
            DBH,
            BA,
            path_idx,
            path_prob
        )]
    }

    psb <- rbindlist(path_stem_ba_list, use.names = TRUE)
    if (nrow(psb) == 0L) next

    # ---- 4b) Per-stem BA posterior summary --------------------------------
    # Weighted mean and quantiles of BA per (StemID, CensusID)
    stem_summary <- psb[, {
        w <- path_prob / sum(path_prob)
        ba_vals <- BA
        wm <- sum(w * ba_vals)
        wv <- sum(w * (ba_vals - wm)^2)
        # Weighted quantiles: order by BA, cumsum weights
        ord <- order(ba_vals)
        cw  <- cumsum(w[ord])
        q025 <- ba_vals[ord][which(cw >= 0.025)[1]]
        q975 <- ba_vals[ord][which(cw >= 0.975)[1]]
        q50  <- ba_vals[ord][which(cw >= 0.50)[1]]
        list(
            BA_mean = wm,
            BA_sd   = sqrt(wv),
            BA_median = q50,
            BA_q025 = q025,
            BA_q975 = q975,
            n_unique_dbh = uniqueN(ba_vals),
            n_paths = .N
        )
    }, by = .(Tag, StemID, CensusID)]

    # Also attach the MAP BA for reference
    map_stem <- map_obs[, .(Tag = tag_id, StemID = ReconstructedStemID,
        CensusID, BA_map = ba_cm2(DBH))]
    stem_summary <- merge(stem_summary, map_stem,
        by = c("Tag", "StemID", "CensusID"), all.x = TRUE)

    stem_results[[length(stem_results) + 1L]] <- stem_summary

    # ---- 4c) Tag-level BA per census summary -----------------------------
    tag_ba_per_path <- psb[, .(BA_total = sum(BA)), by = .(Tag, CensusID, path_idx, path_prob)]

    tag_summary <- tag_ba_per_path[, {
        w <- path_prob / sum(path_prob)
        ba_vals <- BA_total
        wm <- sum(w * ba_vals)
        wv <- sum(w * (ba_vals - wm)^2)
        ord <- order(ba_vals)
        cw <- cumsum(w[ord])
        list(
            BA_total_mean = wm,
            BA_total_sd   = sqrt(wv),
            BA_total_q025 = ba_vals[ord][which(cw >= 0.025)[1]],
            BA_total_q975 = ba_vals[ord][which(cw >= 0.975)[1]],
            n_paths = .N
        )
    }, by = .(Tag, CensusID)]

    # Attach MAP total BA
    map_tag_ba <- map_obs[, .(BA_total_map = sum(ba_cm2(DBH))), by = .(CensusID)]
    map_tag_ba[, Tag := tag_id]
    tag_summary <- merge(tag_summary, map_tag_ba,
        by = c("Tag", "CensusID"), all.x = TRUE)

    tag_results[[length(tag_results) + 1L]] <- tag_summary

    # ---- 4d) Per-stem BA growth rates ------------------------------------
    # For each posterior path, compute annualised BA change per stem
    # (dBA/dt between consecutive censuses for the same StemID)
    setorder(psb, path_idx, StemID, CensusID)

    growth_entries <- list()
    for (pidx in unique(psb$path_idx)) {
        pdata <- psb[path_idx == pidx]
        pw <- pdata$path_prob[1]
        for (sid in unique(pdata$StemID)) {
            sdata <- pdata[StemID == sid]
            if (nrow(sdata) < 2L) next
            for (r in 2:nrow(sdata)) {
                c0 <- sdata$CensusID[r - 1]
                c1 <- sdata$CensusID[r]
                md0 <- census_dates[J(c0), MeanDate]
                md1 <- census_dates[J(c1), MeanDate]
                if (length(md0) == 0 || length(md1) == 0) next
                iv <- (md1[1] - md0[1]) / 365.25
                if (!is.finite(iv) || iv <= 0) iv <- 5.0
                ba0 <- sdata$BA[r - 1]
                ba1 <- sdata$BA[r]
                dba_rate <- (ba1 - ba0) / iv
                growth_entries[[length(growth_entries) + 1L]] <- data.table(
                    Tag = tag_id,
                    StemID = sid,
                    CensusID_from = c0,
                    CensusID_to   = c1,
                    BA_from = ba0,
                    BA_to   = ba1,
                    interval_yr = iv,
                    dBA_rate = dba_rate,
                    path_idx = pidx,
                    path_prob = pw
                )
            }
        }
    }

    if (length(growth_entries) > 0L) {
        ge <- rbindlist(growth_entries, use.names = TRUE)

        # Summarise growth rate posterior per (StemID, CensusID_from, CensusID_to)
        growth_summary <- ge[, {
            w <- path_prob / sum(path_prob)
            vals <- dBA_rate
            wm <- sum(w * vals)
            wv <- sum(w * (vals - wm)^2)
            ord <- order(vals)
            cw <- cumsum(w[ord])
            list(
                dBA_rate_mean = wm,
                dBA_rate_sd   = sqrt(wv),
                dBA_rate_q025 = vals[ord][which(cw >= 0.025)[1]],
                dBA_rate_q975 = vals[ord][which(cw >= 0.975)[1]],
                n_paths = .N
            )
        }, by = .(Tag, StemID, CensusID_from, CensusID_to)]

        growth_results[[length(growth_results) + 1L]] <- growth_summary
    }
}

# ---- 5) Combine and write outputs -------------------------------------

cat("[BA uncertainty] Processed", length(tags_with_posteriors),
    "tags with posterior samples\n")

# Tag-level BA
if (length(tag_results) > 0L) {
    tag_dt <- rbindlist(tag_results, use.names = TRUE, fill = TRUE)
    setorder(tag_dt, Tag, CensusID)
    tag_out <- file.path(RUN_DIR, "basal_area_uncertainty_tag.csv")
    fwrite(tag_dt, tag_out)
    cat("[BA uncertainty] Wrote tag-level BA summary:", tag_out, "\n")
    cat("  ", nrow(tag_dt), "rows\n")
} else {
    cat("[BA uncertainty] WARNING: no tag-level results produced.\n")
}

# Stem-level BA
if (length(stem_results) > 0L) {
    stem_dt <- rbindlist(stem_results, use.names = TRUE, fill = TRUE)
    setorder(stem_dt, Tag, StemID, CensusID)
    stem_out <- file.path(RUN_DIR, "basal_area_uncertainty_stem.csv")
    fwrite(stem_dt, stem_out)
    cat("[BA uncertainty] Wrote stem-level BA summary:", stem_out, "\n")
    cat("  ", nrow(stem_dt), "rows\n")
} else {
    cat("[BA uncertainty] WARNING: no stem-level results produced.\n")
}

# BA growth rates
if (length(growth_results) > 0L) {
    growth_dt <- rbindlist(growth_results, use.names = TRUE, fill = TRUE)
    setorder(growth_dt, Tag, StemID, CensusID_from, CensusID_to)
    growth_out <- file.path(RUN_DIR, "basal_area_uncertainty_growth.csv")
    fwrite(growth_dt, growth_out)
    cat("[BA uncertainty] Wrote BA growth rate summary:", growth_out, "\n")
    cat("  ", nrow(growth_dt), "rows\n")
} else {
    cat("[BA uncertainty] WARNING: no growth rate results produced.\n")
}

# ---- 6) Print summary statistics ----------------------------------------

if (length(stem_results) > 0L) {
    cat("\n============ BASAL AREA UNCERTAINTY SUMMARY ============\n\n")

    # How many stems have non-trivial uncertainty?
    uncertain <- stem_dt[n_unique_dbh > 1]
    cat(sprintf("Stems × censuses with identity uncertainty: %d / %d (%.1f%%)\n",
        nrow(uncertain), nrow(stem_dt),
        100 * nrow(uncertain) / max(nrow(stem_dt), 1)))

    if (nrow(uncertain) > 0L) {
        cat(sprintf("  Mean BA SD (uncertain only): %.2f cm^2\n",
            mean(uncertain$BA_sd)))
        cat(sprintf("  Max BA SD: %.2f cm^2 (Tag %d, Stem %d, Census %d)\n",
            max(uncertain$BA_sd),
            uncertain$Tag[which.max(uncertain$BA_sd)],
            uncertain$StemID[which.max(uncertain$BA_sd)],
            uncertain$CensusID[which.max(uncertain$BA_sd)]))

        # 95% CI width as % of mean BA
        uncertain[, ci_width_pct := 100 * (BA_q975 - BA_q025) / pmax(BA_mean, 1e-6)]
        cat(sprintf("  Mean 95%% CI width (as %% of mean BA): %.1f%%\n",
            mean(uncertain$ci_width_pct, na.rm = TRUE)))
        cat(sprintf("  Max 95%% CI width: %.1f%%\n",
            max(uncertain$ci_width_pct, na.rm = TRUE)))
    }

    # Tags where tag-level BA is constant (verification)
    if (length(tag_results) > 0L) {
        tag_const <- tag_dt[BA_total_sd < 1e-8]
        cat(sprintf("\nTag × census pairs with zero BA uncertainty: %d / %d\n",
            nrow(tag_const), nrow(tag_dt)))
        cat("  (Tag-level BA is invariant to identity assignment when all\n")
        cat("   observations are present in all posterior paths.)\n")
    }
}

if (length(growth_results) > 0L) {
    cat("\n--- BA Growth Rate Uncertainty ---\n")
    uncertain_gr <- growth_dt[dBA_rate_sd > 1e-8]
    cat(sprintf("Growth intervals with uncertainty: %d / %d (%.1f%%)\n",
        nrow(uncertain_gr), nrow(growth_dt),
        100 * nrow(uncertain_gr) / max(nrow(growth_dt), 1)))
    if (nrow(uncertain_gr) > 0L) {
        cat(sprintf("  Mean dBA/dt SD (uncertain only): %.2f cm^2/yr\n",
            mean(uncertain_gr$dBA_rate_sd)))
    }
}

cat("\nDone.\n")
