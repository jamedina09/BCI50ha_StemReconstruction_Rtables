############################################################
### basal_area_uncertainty.R — Individual-level basal area
###   with posterior uncertainty from stem identification
###
### Computes basal area (BA) and BA change at the individual (tag)
### level by tracking stems through posterior path samples from
### the DP / probabilistic reconstruction.
###
### For each tag (individual tree):
###   1. Per-stem BA = pi/4 * (DBH/100)^2 (m^2), using reconstructed IDs
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

rm(list = ls())

library(data.table)

suppressPackageStartupMessages({
    library(ggplot2)
    library(patchwork)
    library(scales)
})

bci_stem_nums <- as.character(1:9)

census_list <- lapply(bci_stem_nums, function(num) {
  filepath <- paste0("./BCI_stem_reconstruction/DATA/RTABLES/bci.stem", num, ".Rdata")
  if (!file.exists(filepath)) stop("Missing census file: ", filepath)
  load(filepath)
  get(paste0("bci.stem", num))
})
names(census_list) <- paste0("bci.stem", bci_stem_nums)

rec <- data.table::rbindlist(census_list, fill = TRUE, idcol = "censusID")
rm(census_list, bci_stem_nums)

rec[, date_quad_census := median(ExactDate, na.rm = TRUE), .(quadrat, CensusID)]
rec[, date_plot_census := median(ExactDate, na.rm = TRUE), .(CensusID)]

rec[, ExactDate := fifelse(is.na(ExactDate), date_quad_census, ExactDate)]
rec[is.na(ExactDate), ExactDate := date_plot_census]

rec[, `:=`(date_quad_census = NULL, date_plot_census = NULL)]

post_file <- file.path("./BCI_stem_reconstruction/DATA/POSTERIORS/posterior_sampled_paths.rds")

# ---- 1. Load reconstruction --------------------------------------------

# census_dates <- rec[!is.na(ExactDate),
#     .(MeanDate = mean(as.numeric(as.Date(ExactDate)), na.rm = TRUE)),
#     by = CensusID
# ]
# setkey(census_dates, CensusID)

ba_m2 <- function(dbh) pi / 4 * (dbh / 100)^2

# ---- 2. Map posterior files to tags -------------------------------------

parse_recon <- function(recon_str) {
    pairs <- strsplit(recon_str, ";", fixed = TRUE)[[1]]
    parts <- strsplit(pairs, ":", fixed = TRUE)
    data.table(
        obs_row_id = as.integer(vapply(parts, `[`, character(1), 1L)),
        StemID     = as.integer(vapply(parts, `[`, character(1), 2L))
    )
}

set.seed(123)
all_trees <- sample(unique(rec[!is.na(obs_row_id)]$treeID), 10)
# post_tag_map <- list()

# split all trees in rec by treeid in post_tag_map
# post_tag_map <- split(rec, rec$treeID) # placeholder: map each treeID to itself for now


# has_posteriors <- file.exists(post_file)

# if (has_posteriors) {
#     post_files <- list.files(post_dir,
#         pattern = "_paths\\.csv$",
#         full.names = TRUE
#     )

# post_file <- data.table(readRDS(post_file))

#     for (pf in post_files) {
#         tag_match <- regmatches(
#             basename(pf),
#             regexec("^tag_([^_]+)_posterior_samples", basename(pf))
#         )[[1]]
#         if (length(tag_match) < 2L) next
#         tag_str <- tag_match[2]
#         if (tag_str == "NA") {
#             # Single-tag run: infer tag from reconstruction data
#             if (length(all_tags) == 1L) {
#                 post_tag_map[[as.character(all_tags[1])]] <- pf
#             }
#         } else {
#             post_tag_map[[tag_str]] <- pf
#         }
#     }
# }

# cat("[BA] Posterior files mapped:", length(post_tag_map), "tags\n")

# ---- 3. MAP tag-level BA per census -------------------------------------

map_stem_ba <- rec[!is.na(dbh) & !is.na(stemID),
    .(BA = ba_m2(dbh)),
    by = .(treeID, CensusID, stemID)
]

tree_census <- map_stem_ba[, .(
    TotalBA_m2 = sum(BA),
    NumStems   = .N
), by = .(treeID, CensusID)]

setorder(tree_census, treeID, CensusID)
# tree_census <- merge(tree_census, census_dates, by = "CensusID", all.x = TRUE)
# tree_census[, Year := 1970 + MeanDate / 365.25]

cat("[BA] Tag x census rows:", nrow(tree_census), "\n")

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

for (tg in all_trees) {
    # tg <- 100L
    tag_rec <- rec[treeID == tg & !is.na(dbh) & !is.na(stemID)]
    censuses <- sort(unique(tag_rec$CensusID))
    dates <- tag_rec[, .(Date = median(ExactDate, na.rm = TRUE)), by = CensusID]
    if (length(censuses) < 2L) next
census_pairs <- data.table(
    c_from = dates$CensusID[-length(dates$CensusID)],
    c_to   = dates$CensusID[-1],
    Date_from = dates$Date[-length(dates$Date)],
    Date_to   = dates$Date[-1]
)
    stem_dt <- tag_rec[, .(
        StemID = stemID, CensusID,
        BA = ba_m2(dbh)
    )]
    decomp <- decompose_ba_change(stem_dt, census_pairs)
    decomp[, treeID := tg]
# merge dates
decomp <- merge(decomp, dates[, .(CensusID, Date)], by.x = "CensusID_from", by.y = "CensusID", all.x = TRUE)
decomp <- merge(decomp, dates[, .(CensusID, Date)], by.x = "CensusID_to", by.y = "CensusID", all.x = TRUE, suffixes = c("_from", "_to"))
    map_change_list[[length(map_change_list) + 1L]] <- decomp
}

map_change <- if (length(map_change_list) > 0L) {
    rbindlist(map_change_list, use.names = TRUE)
} else {
    data.table()
}

if (nrow(map_change) > 0L) {
    map_change[, Interval_yr := as.numeric(difftime(Date_to, Date_from, units = "days")) / 365.25]
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
all_path_decomp_list <- list()
all_paths_by_tag  <- list()  # tag -> data.table(path_idx, path_prob, obs_row_id, StemID, CensusID, DBH, BA)

post_full <- data.table(readRDS(post_file))

# names(post_tag_map)[1:10]

# for (tg_str in names(post_tag_map)[1:10]) {
for (tg_str in all_trees) {
    # tg_tr <- 43329
    cat("[BA] Processing tag:", tg_str, "\n")
    tg <- tg_str
    # pf <- post_file

    tag_rec <- rec[treeID == tg]
    if (nrow(tag_rec) == 0L) next

    obs_lookup <- tag_rec[, .(obs_row_id, CensusID, dbh)]
    censuses <- sort(unique(tag_rec$CensusID))
    if (length(censuses) < 2L) next

    census_pairs <- data.table(
        c_from = censuses[-length(censuses)],
        c_to   = censuses[-1]
    )

post <- post_full[Tag == unique(tag_rec$tag)]
    n_paths <- nrow(post)

    # Parse all paths into one table
    all_paths <- rbindlist(lapply(seq_len(n_paths), function(i) {
        dt <- parse_recon(post$recon[i])
        dt[, `:=`(path_idx = i, path_prob = post$path_prob[i])]
        dt
    }))
    all_paths <- merge(all_paths, obs_lookup, by = "obs_row_id", all.x = TRUE)
    all_paths[, BA := ba_m2(dbh)]
    all_paths <- all_paths[!is.na(dbh)]

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
    all_path_decomp_list[[length(all_path_decomp_list) + 1L]] <- path_decomp

    # Retain per-path stem assignments for trajectory figures
    ap <- copy(all_paths)
    ap[, treeID := tg]
    all_paths_by_tag[[as.character(tg)]] <- ap

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
    tag_census[, .(Tag, CensusID, Year, TotalBA_m2, NumStems)],
    tag_census_out
)
cat("[BA] Wrote:", tag_census_out, "\n")

tag_change_out <- file.path(RUN_DIR, "basal_area_tag_change.csv")
fwrite(tag_change, tag_change_out)
cat("[BA] Wrote:", tag_change_out, "\n")

