library(data.table); library(here)
source(here("dp_global","R","dp_global_utils.R"))
source(here("dp_global","R","dp_global_states.R"))
source(here("dp_global","R","dp_global_matchers.R"))
source(here("dp_global","R","dp_global_bio.R"))
source(here("dp_global","R","dp_global_dp.R"))

d <- as.data.table(readRDS(here("bci_data","bci_multistem_xrun_debug.rds")))
tag_data <- d[Tag == "115427"]

# Trace what the DP sees for census_range
anchor_start <- 7L
resprout_regex <- "\\b(R|RP|RF|RT|QR)\\b"
obs_censuses <- sort(unique(tag_data[!is.na(DBH), CensusID]))
cat("Observed censuses (non-NA DBH):", obs_censuses, "\n")
all_census <- sort(unique(tag_data$CensusID))
cat("All censuses:", all_census, "\n")

r4 <- tag_data[CensusID == 4 & is.na(DBH)]
cat("C4 rows:", nrow(r4), "\n")
cat("C4 ListOfTSM:", paste(r4$ListOfTSM), "\n")
cat("C4 match perl=TRUE:", grepl(resprout_regex, r4$ListOfTSM, perl=TRUE), "\n")

# Now simulate the N_census and barrier detection loop as done in the DP
# The DP uses: census_range = first_obs .. anchor_start
first_obs_census <- min(tag_data[!is.na(DBH), CensusID])
last_obs_census  <- max(tag_data[!is.na(DBH), CensusID])
eff_anchor <- if (last_obs_census < anchor_start) last_obs_census else anchor_start
first_obs <- min(tag_data[CensusID <= eff_anchor & !is.na(DBH), CensusID])
census_range <- sort(unique(c(
    tag_data[CensusID >= first_obs & CensusID <= eff_anchor, CensusID],
    eff_anchor
)))
census_range <- census_range[census_range <= eff_anchor]
n_census <- length(census_range)
anchor_pos <- which(census_range == eff_anchor)
cat("census_range:", census_range, "\n")
cat("anchor_pos:", anchor_pos, "\n")

.has_tsm <- "ListOfTSM" %in% names(tag_data)
.has_na_r_barrier <- logical(n_census)
for (.p0 in seq_len(n_census)) {
    .idx0 <- tag_data[CensusID == census_range[.p0] & !is.na(DBH), which = TRUE]
    if (!.has_tsm || length(.idx0) > 0L) next   # only enter else branch when idx0 is empty
    # .idx0 is empty
    .na_rows <- tag_data[CensusID == census_range[.p0] & is.na(DBH), which = TRUE]
    if (length(.na_rows) > 0L) {
        .na_tsm <- tag_data$ListOfTSM[.na_rows]
        match_result <- !is.na(.na_tsm) & grepl(resprout_regex, .na_tsm, perl = TRUE)
        cat(sprintf("  p=%d census=%d: na_rows=%d, tsm=%s, match=%s\n",
            .p0, census_range[.p0], length(.na_rows),
            paste(.na_tsm, collapse=","), paste(match_result, collapse=",")))
        .has_na_r_barrier[.p0] <- any(match_result)
    }
}
cat("has_na_r_barrier:", .has_na_r_barrier, "\n")
barrier_positions <- which(.has_na_r_barrier & seq_len(n_census) < anchor_pos)
cat("barrier_positions:", barrier_positions, "\n")
