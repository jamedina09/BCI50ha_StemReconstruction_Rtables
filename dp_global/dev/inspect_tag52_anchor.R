# Inspect anchor assignment for Tag=52
r_files <- list.files('dp_global/R', full.names = TRUE)
r_files <- r_files[endsWith(r_files, '.R')]
for (f in r_files) source(f)

library(data.table)

xraw <- data.table::fread(here('data_simulation','data','simulated_data_1.csv'))
# normalize species to the forced single-species label used by the runner
xraw[, species := 'all']

dtg <- xraw[Tag == 52]
cat('Rows for Tag 52:', nrow(dtg), '\n')
cat('Unique CensusID:', sort(unique(dtg$CensusID)), '\n')

# Prepare minimal copy of tag data (don't call attach_bio_columns in this standalone inspect)
xrun <- dtg

anchor_start <- 7L
obs_census_all <- sort(unique(xrun$CensusID[!is.na(xrun$DBH)]))
last_obs_census <- if (length(obs_census_all) > 0L) as.integer(max(obs_census_all)) else NA_integer_
cat('last_obs_census =', last_obs_census, '\n')
if (!is.na(last_obs_census) && last_obs_census < anchor_start) anchor_start <- last_obs_census
cat('Using anchor_start =', anchor_start, '\n')

# Find anchor rows all and anchor_idx
anchor_rows_all <- xrun[CensusID == anchor_start]
cat('nrow(anchor_rows_all)=', nrow(anchor_rows_all), '\n')
cat('DBH NA count =', sum(is.na(anchor_rows_all$DBH)), 'TrueStemID NA count =', sum(is.na(anchor_rows_all$TrueStemID)), '\n')

anchor_idx <- which(xrun$CensusID == anchor_start & !is.na(xrun$DBH))
cat('anchor_idx length =', length(anchor_idx), '\n')
cat('anchor_idx values =', paste(anchor_idx, collapse=','), '\n')

current_max <- suppressWarnings(max(xrun$TrueStemID, na.rm = TRUE))
cat('current_max TrueStemID =', current_max, '\n')
if (!is.finite(current_max)) current_max <- 0L
prov_ids <- as.integer(seq.int(from = current_max + 1L, length.out = length(anchor_idx)))
cat('prov_ids len =', length(prov_ids), 'prov_ids =', paste(prov_ids, collapse=','), '\n')

cat('ConstraintViolation length =', if ('ConstraintViolation' %in% names(xrun)) length(xrun$ConstraintViolation) else 'missing', '\n')

# Try assigning (simulate what DP does)
tryCatch({
  xrun$TrueStemID[anchor_idx] <- prov_ids
  xrun$ReconstructedStemID[anchor_idx] <- prov_ids
  xrun$ReconstructionMethod[anchor_idx] <- 'provisional_dp'
  xrun$ConstraintViolation[anchor_idx] <- FALSE
  cat('Assignment succeeded\n')
  cat('ConstraintViolation slice:', paste(xrun$ConstraintViolation[anchor_idx], collapse=','), '\n')
}, error = function(e) {
  cat('Assignment error:', conditionMessage(e), '\n')
  # Show types and lengths for debugging
  sapply(c('TrueStemID','ReconstructedStemID','ReconstructionMethod','ConstraintViolation'), function(col) {
    if (col %in% names(xrun)) {
      vc <- xrun[[col]]
      cat(sprintf('%s: len=%d type=%s NAcount=%d\n', col, length(vc), typeof(vc), sum(is.na(vc))))
    } else {
      cat(sprintf('%s: <missing>\n', col))
    }
  })
})
