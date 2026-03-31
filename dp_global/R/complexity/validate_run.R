library(data.table)

# Find the latest all-tags run
dirs <- list.dirs("dp_global/output", recursive = FALSE)
all_tag_dirs <- dirs[grepl("allT", dirs)]
new_dir <- tail(sort(all_tag_dirs), 1L)
cat("Using:", new_dir, "\n")

new <- fread(file.path(new_dir, "stem_reconstruction_dp_global_rcpp.csv"))

# Non-M tags accuracy (tags < 900)
non_m <- new[Tag < 900 & !is.na(DBH) & !is.na(OriginalStemID) & !is.na(ReconstructedStemID)]
acc_non_m <- mean(non_m$OriginalStemID == non_m$ReconstructedStemID)
cat(sprintf("Non-M accuracy : %.1f%% (%d/%d rows)\n",
    100 * acc_non_m,
    sum(non_m$OriginalStemID == non_m$ReconstructedStemID),
    nrow(non_m)))

# M test tags (901, 902, 903)
m_tags <- new[Tag %in% c(901L, 902L, 903L) & !is.na(OriginalStemID) & !is.na(ReconstructedStemID)]
if (nrow(m_tags) > 0L) {
    acc_m <- mean(m_tags$OriginalStemID == m_tags$ReconstructedStemID)
    cat(sprintf("M-tag accuracy : %.1f%% (%d/%d rows)\n",
        100 * acc_m,
        sum(m_tags$OriginalStemID == m_tags$ReconstructedStemID),
        nrow(m_tags)))
    print(m_tags[order(Tag, CensusID),
                 .(Tag, CensusID, OriginalStemID, ReconstructedStemID, DBH, ListOfTSM)])
} else {
    cat("M-tags not found in output\n")
}

# Overall
all_rows <- new[!is.na(OriginalStemID) & !is.na(ReconstructedStemID)]
acc_all <- mean(all_rows$OriginalStemID == all_rows$ReconstructedStemID)
cat(sprintf("Overall        : %.1f%% (%d/%d rows)\n",
    100 * acc_all,
    sum(all_rows$OriginalStemID == all_rows$ReconstructedStemID),
    nrow(all_rows)))

# Show any wrong M-tag rows
if (nrow(m_tags) > 0L) {
    wrong <- m_tags[OriginalStemID != ReconstructedStemID]
    if (nrow(wrong) > 0L) {
        cat("\nWRONG M-tag rows:\n")
        print(wrong[order(Tag, CensusID),
                    .(Tag, CensusID, OriginalStemID, ReconstructedStemID, DBH, ListOfTSM)])
    } else {
        cat("\nAll M-tag rows correct.\n")
    }
}
