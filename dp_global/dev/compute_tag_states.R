library(data.table)
source("dp_global/R/dp_global_utils.R")
xs <- fread("data_simulation/data/simulated_for_benchmark.csv")
anchor_start <- 7L
out <- list()

for (tg in unique(xs$Tag)) {
    td <- xs[Tag == tg & CensusID <= anchor_start]
    oc_all <- sapply(seq_len(anchor_start), function(cc) nrow(td[CensusID == cc & !is.na(DBH)]))
    oc <- oc_all[oc_all > 0]
    if (length(oc) == 0) next
    births_needed <- if (length(oc) >= 2) sum(pmax(0L, diff(as.integer(oc)))) else 0L
    K_from_counts <- as.integer(oc[1]) + as.integer(births_needed)
    K_base <- max(as.integer(max(oc)), as.integer(K_from_counts))
    K <- K_base
    n_states <- vapply(oc, function(n) count_injective_states(K, n), numeric(1))
    out[[tg]] <- data.table(Tag = tg, K = K, MaxStates = max(n_states, na.rm = TRUE), TransitionComputations = sum(head(n_states, -1) * tail(n_states, -1)), ObsCounts = paste(oc, collapse = ","))
}

res <- rbindlist(out)

predict_time <- function(N, unit = c("sec", "min", "hour")) {
    unit <- match.arg(unit)
    logN <- log10(N)
    logT <- -1.45817 +
        (-0.07313 * logN) +
        (0.14164 * logN^2)
    secs <- 10^logT
    out <- switch(unit,
        sec = secs,
        min = secs / 60,
        hour = secs / 3600
    )
    out
}

res[, PredictedSec := round(predict_time(TransitionComputations, unit = "hour"), 4)]

print(res)
write.csv(res, file = file.path("..", "..", "dp_global", "dev", "tag_state_summary.csv"), row.names = FALSE)
message("Wrote tag_state_summary.csv")
