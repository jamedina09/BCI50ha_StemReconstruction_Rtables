d <- read.csv("/Users/medinaja/Library/CloudStorage/OneDrive-SmithsonianInstitution/STRI/STEM_TABLES_FORESTGEO/3_BCI/2_STEM_IDENTIFICATION/example.csv", stringsAsFactors = FALSE)

compute <- function(d, anchor, slack_tracks = 1L, max_tracks = 30L, max_states = 50000, window = NULL) {
  d <- d[order(d$CensusID), ]
  first_obs_census <- if (any(d$CensusID <= anchor & !is.na(d$DBH))) min(unique(d$CensusID[d$CensusID <= anchor & !is.na(d$DBH)])) else NA_integer_
  if (is.na(first_obs_census)) return(list(reason = "no_obs"))
  census_range <- seq(from = first_obs_census, to = anchor)
  if (!is.null(window)) {
    low <- max(first_obs_census, anchor - window + 1L)
    census_range <- seq(from = low, to = anchor)
  }

  obs_counts <- sapply(census_range, function(cc) sum(d$CensusID == cc & !is.na(d$DBH)))
  max_obs <- if (length(obs_counts) > 0) max(obs_counts) else 0
  births_needed <- if (length(obs_counts) >= 2) sum(pmax(0, diff(obs_counts))) else 0
  K_from_counts <- if (length(obs_counts) > 0) obs_counts[1L] + births_needed else 0
  anchor_obs <- d[d$CensusID == anchor & !is.na(d$DBH), ]
  anchor_ids <- unique(anchor_obs$TrueStemID)
  K_base <- max(length(anchor_ids), max_obs, K_from_counts)
  K_target <- K_base + as.integer(slack_tracks)
  K <- min(as.integer(K_target), as.integer(max_tracks))

  Pvals <- sapply(obs_counts, function(n_obs) {
    if (n_obs > 0 && n_obs > K) return(NA_real_)
    if (n_obs == 0) return(1)
    nstates <- 1
    for (j in 0:(n_obs - 1)) nstates <- nstates * (K - j)
    nstates
  })

  exceeds <- which(sapply(seq_along(Pvals), function(i) {
    n <- obs_counts[i]
    if (n > K) return(TRUE)
    if (is.na(Pvals[i])) return(TRUE)
    if (is.finite(Pvals[i]) && Pvals[i] > max_states) return(TRUE)
    FALSE
  }))

  list(anchor = anchor,
       params = list(slack = slack_tracks, max_tracks = max_tracks, max_states = max_states, window = window),
       first_obs_census = first_obs_census,
       census_range = census_range,
       obs_counts = obs_counts,
       K = K,
       Pvals = Pvals,
       first_exceed = if (length(exceeds) > 0) census_range[exceeds[1]] else NA_integer_,
       all_ok = length(exceeds) == 0)
}

anchors <- 7:9
scenarios <- list(
  list(slack = 1, max_states = 50000, window = NULL),
  list(slack = 0, max_states = 50000, window = NULL),
  list(slack = 1, max_states = 1e7, window = NULL),
  list(slack = 0, max_states = 1e7, window = NULL),
  list(slack = 1, max_states = 50000, window = 3),
  list(slack = 0, max_states = 50000, window = 3),
  list(slack = 1, max_states = 50000, window = 2)
)

cat("Enumeration diagnostics for example.csv\n")
for (a in anchors) {
  cat(sprintf("\n=== Anchor %d ===\n", a))
  for (s in scenarios) {
    res <- compute(d, a, s$slack, 30L, s$max_states, s$window)
    cat(sprintf("params: slack=%d max_states=%g window=%s\n", s$slack, s$max_states, if (is.null(s$window)) "all" else as.character(s$window)))
    cat("census_range:", paste(res$census_range, collapse = ","), "\n")
    cat("obs_counts:", paste(res$obs_counts, collapse = ","), "\n")
    cat("K:", res$K, "\n")
    cat("Pvals per census (or NA if n_obs>K):\n")
    for (i in seq_along(res$census_range)) {
      cc <- res$census_range[i]
      pv <- res$Pvals[i]
      cat(sprintf("  Census %d: n_obs=%d P=%s\n", cc, res$obs_counts[i], if (is.na(pv)) "NA" else format(pv, digits = 6, scientific = TRUE)))
    }
    if (!is.na(res$first_exceed)) cat("-> First census exceeding max_states or impossible:", res$first_exceed, "\n") else cat("-> All censuses OK for enumeration under these params\n")
    cat("---\n")
  }
}
