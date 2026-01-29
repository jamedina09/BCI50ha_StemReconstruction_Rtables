# Quick check: call DP on small toy case where anchor census has DBH but missing TrueStemID
r_files <- list.files('dp_global/R', full.names = TRUE)
r_files <- r_files[endsWith(r_files, '.R')]
for (f in r_files) source(f)

library(data.table)
cat('Function exists?', exists('match_stems_dp_global_backward_marginals_batch', inherits = TRUE), '\n')
cat('Search path:\n')
print(search())
dt <- data.table(
  Tag = rep(999L, 5),
  CensusID = c(3L,4L,5L,6L,7L),
  DBH = c(NA_real_, NA_real_, 1.0, 1.5, 2.0),
  TrueStemID = as.integer(NA),
  ReconstructedStemID = as.integer(NA),
  ReconstructionMethod = NA_character_
)

cat('Running DP with allow_provisional_anchor=TRUE\n')
res <- match_stems_dp_global_backward_marginals_batch(copy(dt), anchor_start = 7L, min_growth=-10, max_growth=30, allow_provisional_anchor = TRUE, verbose=TRUE)
cat('ReconstructionMethod:', paste(unique(res$ReconstructionMethod), collapse=', '), '\n')
print(res)

cat('\nRunning DP with allow_provisional_anchor=FALSE\n')
res2 <- match_stems_dp_global_backward_marginals_batch(copy(dt), anchor_start = 7L, min_growth=-10, max_growth=30, allow_provisional_anchor = FALSE, verbose=TRUE)
cat('ReconstructionMethod:', paste(unique(res2$ReconstructionMethod), collapse=', '), '\n')
print(res2)
