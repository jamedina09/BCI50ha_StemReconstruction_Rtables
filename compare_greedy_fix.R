#!/usr/bin/env Rscript
# Compare pre- vs post-growth-aware-greedy-fix outputs
library(data.table)

pre_dir  <- "baseline_pre_greedy_fix"
# Find latest post-change output dirs
all_dirs <- list.dirs("dp_global/output", recursive = FALSE)
all_dirs <- all_dirs[grepl("^20260411_123", basename(all_dirs))]  # post-change runs

# --- Helper: compare two CSV files ----------------------------------------
compare_csv <- function(pre_file, post_file, label) {
    if (!file.exists(pre_file))  { cat(label, ": PRE file missing\n"); return(invisible()) }
    if (!file.exists(post_file)) { cat(label, ": POST file missing\n"); return(invisible()) }
    a <- fread(pre_file,  na.strings = "")
    b <- fread(post_file, na.strings = "")

    cat("\n====", label, "====\n")
    cat("  Rows: pre=", nrow(a), " post=", nrow(b), "\n")

    # Compare ReconstructedStemID
    if ("ReconstructedStemID" %in% names(a) && "ReconstructedStemID" %in% names(b)) {
        id_a <- a$ReconstructedStemID
        id_b <- b$ReconstructedStemID
        diffs <- which(id_a != id_b | xor(is.na(id_a), is.na(id_b)))
        # filter out rows where BOTH are NA (post-anchor / NA-DBH)
        diffs <- diffs[!(is.na(id_a[diffs]) & is.na(id_b[diffs]))]
        cat("  ReconstructedStemID diffs:", length(diffs), "\n")
    }

    # Compare ReconstructionMethod
    if ("ReconstructionMethod" %in% names(a) && "ReconstructionMethod" %in% names(b)) {
        m_a <- a$ReconstructionMethod
        m_b <- b$ReconstructionMethod
        mdiffs <- which(m_a != m_b | xor(is.na(m_a), is.na(m_b)))
        mdiffs <- mdiffs[!(is.na(m_a[mdiffs]) & is.na(m_b[mdiffs]))]
        cat("  ReconstructionMethod diffs:", length(mdiffs))
        if (length(mdiffs) > 0 && length(mdiffs) <= 20) {
            changes <- paste0(m_a[mdiffs], "->", m_b[mdiffs])
            cat("  [", paste(unique(changes), collapse = ", "), "]")
        }
        cat("\n")
    }

    # Compare posteriors (Top1 ID, Top1 Prob, Entropy)
    post_diffs <- 0L
    for (col in c("DP_PosteriorTop1ID", "DP_PosteriorTop1Prob",
                   "DP_PosteriorTop2ID", "DP_PosteriorTop2Prob",
                   "DP_PosteriorEntropy")) {
        if (col %in% names(a) && col %in% names(b)) {
            va <- a[[col]]; vb <- b[[col]]
            d <- which(abs(va - vb) > 1e-6 | xor(is.na(va), is.na(vb)))
            d <- d[!(is.na(va[d]) & is.na(vb[d]))]
            post_diffs <- post_diffs + length(d)
        }
    }
    cat("  Posterior column diffs:", post_diffs, "\n")

    # Check for growth violations in POST output
    if ("ReconstructedStemID" %in% names(b) && "CensusID" %in% names(b) && "DBH" %in% names(b)) {
        prob_rows <- which(b$ReconstructionMethod == "probabilistic")
        if (length(prob_rows) > 0) {
            bp <- b[prob_rows]
            sids <- unique(bp$ReconstructedStemID[!is.na(bp$ReconstructedStemID)])
            n_viol <- 0L
            for (sid in sids) {
                rows <- bp[ReconstructedStemID == sid]
                if (nrow(rows) < 2) next
                setorder(rows, CensusID)
                for (r in 2:nrow(rows)) {
                    d1 <- rows$DBH[r-1]; d2 <- rows$DBH[r]
                    if (is.na(d1) || is.na(d2)) next
                    c1 <- rows$CensusID[r-1]; c2 <- rows$CensusID[r]
                    date1 <- as.Date(rows$ExactDate[r-1]); date2 <- as.Date(rows$ExactDate[r])
                    iv <- as.numeric(date2 - date1) / 365.25
                    if (iv <= 0) iv <- 5
                    rate <- (d2 - d1) / iv
                    if (rate < -0.5 || rate > 5.0) n_viol <- n_viol + 1L
                }
            }
            cat("  POST hard-rate violations (probabilistic rows):", n_viol, "\n")
        }
    }

    # Also check PRE for comparison
    if ("ReconstructedStemID" %in% names(a) && "CensusID" %in% names(a) && "DBH" %in% names(a)) {
        prob_rows <- which(a$ReconstructionMethod == "probabilistic")
        if (length(prob_rows) > 0) {
            ap <- a[prob_rows]
            sids <- unique(ap$ReconstructedStemID[!is.na(ap$ReconstructedStemID)])
            n_viol <- 0L
            for (sid in sids) {
                rows <- ap[ReconstructedStemID == sid]
                if (nrow(rows) < 2) next
                setorder(rows, CensusID)
                for (r in 2:nrow(rows)) {
                    d1 <- rows$DBH[r-1]; d2 <- rows$DBH[r]
                    if (is.na(d1) || is.na(d2)) next
                    date1 <- as.Date(rows$ExactDate[r-1]); date2 <- as.Date(rows$ExactDate[r])
                    iv <- as.numeric(date2 - date1) / 365.25
                    if (iv <= 0) iv <- 5
                    rate <- (d2 - d1) / iv
                    if (rate < -0.5 || rate > 5.0) n_viol <- n_viol + 1L
                }
            }
            cat("  PRE  hard-rate violations (probabilistic rows):", n_viol, "\n")
        }
    }
}

# --- Simulated data --------------------------------------------------------
sim_post_dir <- all_dirs[grepl("unknown_allT", all_dirs)]
if (length(sim_post_dir) > 0) {
    sim_post_dir <- sort(sim_post_dir, decreasing = TRUE)[1]
    compare_csv(
        file.path(pre_dir, "simulated_data.csv"),
        file.path(sim_post_dir, "stem_reconstruction_dp_global_rcpp.csv"),
        "Simulated Data"
    )
}

# --- BCI tags ---------------------------------------------------------------
tags <- c("115427", "119453", "123375", "115203", "242799", "246746",
          "277120", "190932", "171486", "220311")
for (tag in tags) {
    pre_file <- file.path(pre_dir, paste0("bci_tag_", tag, ".csv"))
    post_dir <- all_dirs[grepl(paste0("tag", tag), all_dirs)]
    if (length(post_dir) > 0) {
        post_dir <- sort(post_dir, decreasing = TRUE)[1]
        post_file <- file.path(post_dir, "stem_reconstruction_dp_global_rcpp.csv")
        compare_csv(pre_file, post_file, paste0("BCI tag ", tag))
    } else {
        cat("\n==== BCI tag", tag, "==== POST dir not found\n")
    }
}
