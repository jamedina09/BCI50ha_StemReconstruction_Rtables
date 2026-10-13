# ========================================================================
# SCRIPT: compare_viewfulltable_bci_old_and_oct23_2026.R
# PURPOSE: Compare the previous BCI `ViewFullTable` with the one received on
#          2026-10-13 and report which columns have changed.
#
# Table of Contents:
#   0. Setup
#   1. Configuration (paths of the two versions, output folder)
#   2. Load input files (everything as character, so values are compared
#      exactly as written in the files)
#   3. Structure: columns and rows added / removed (rows keyed by DBHID)
#   4. Column-by-column comparison on the rows shared by both versions
#   5. Detail: every changed cell (DBHID, column, old value, new value)
#   6. Rows that exist only in one version
#   7. Final check: compare DBH distributions using density (plots + KS test)
#   8. Save outputs to DATA/CHECKS
#
# Main results (objects to look at):
#   - changed_columns : names of the columns whose values changed
#   - comparison      : one row per column with the number of changed cells
#   - changes         : one row per changed cell
#   - rows_only_old / rows_only_new : records present in only one version
#   - report_densities_*  : DBH density comparison (L1/L2 distances, KS test)
#
# Notes:
#  - `DBHID` is the measurement key of ViewFullTable; rows are matched on it,
#    so the comparison does not depend on row order.
#  - NA vs NA counts as equal; NA vs value counts as a change.
# ========================================================================

# ---- 0. Setup ----
rm(list = ls())

library(data.table)
library(here)

setDTthreads(16) # tune this value for your machine

# ---- 1. Configuration ----
site <- "bci"
key_col <- "DBHID"

file_old <- here::here(
  "BCI_stem_reconstruction", "DATA", "RAW", "ViewFiles_bci_allcensuses",
  paste0("ViewFullTable_", site, ".csv")
)
file_new <- here::here(
  "BCI_stem_reconstruction", "DATA", "RAW", "oct_13_2026",
  paste0("ViewFullTable_", site, ".csv")
)

CHECKS_folder <- here::here("BCI_stem_reconstruction", "DATA", "CHECKS")
dir.create(CHECKS_folder, showWarnings = FALSE, recursive = TRUE)
out_prefix <- file.path(CHECKS_folder, "compare_viewfulltable_old_vs_oct13_2026")
log_file <- paste0(out_prefix, "_log.txt")
cat("Comparison run:", format(Sys.time()), "\n", file = log_file)

# ---- Helper: check ----
# Loud check (the script is run line by line): prints ❌, raises an immediate
# warning, writes to the log, and only then stops.
check <- function(ok, msg) {
  if (isTRUE(all(ok))) {
    cat("✓", msg, "\n")
    return(invisible(TRUE))
  }
  full_msg <- paste0("CHECK FAILED: ", msg)
  cat("❌", full_msg, "\n")
  warning(full_msg, call. = FALSE, immediate. = TRUE)
  cat(full_msg, "\n", file = log_file, append = TRUE)
  stop(full_msg, call. = FALSE)
}

# ---- 2. Load input files ----
check(file.exists(file_old), paste("old file exists:", file_old))
check(file.exists(file_new), paste("new file exists:", file_new))

read_view <- function(file) {
  fread(
    file = file,
    sep = "\t",
    quote = "",
    colClasses = "character",
    na.strings = c("NA", "NULL", "")
  )
}

old <- read_view(file_old)
new <- read_view(file_new)

check(key_col %in% names(old) && key_col %in% names(new), paste(key_col, "present in both versions"))
check(!anyNA(old[[key_col]]) && !anyNA(new[[key_col]]), paste(key_col, "has no missing values"))
check(!anyDuplicated(old[[key_col]]) && !anyDuplicated(new[[key_col]]), paste(key_col, "is unique in both versions"))

# ---- 3. Structure: columns and rows ----
cols_only_old <- setdiff(names(old), names(new))
cols_only_new <- setdiff(names(new), names(old))
cols_shared <- setdiff(intersect(names(old), names(new)), key_col)

ids_only_old <- setdiff(old[[key_col]], new[[key_col]])
ids_only_new <- setdiff(new[[key_col]], old[[key_col]])
ids_shared <- intersect(old[[key_col]], new[[key_col]])

structure_summary <- data.table(
  item = c(
    "rows old", "rows new", "rows shared", "rows only in old", "rows only in new",
    "columns only in old", "columns only in new"
  ),
  value = c(
    nrow(old), nrow(new), length(ids_shared), length(ids_only_old), length(ids_only_new),
    paste(cols_only_old, collapse = ", "), paste(cols_only_new, collapse = ", ")
  )
)
structure_summary

# ---- 4. Column-by-column comparison (shared rows) ----
# Align both versions on the key so row i is the same measurement in both.
old_s <- old[match(ids_shared, old[[key_col]])]
new_s <- new[match(ids_shared, new[[key_col]])]
check(identical(old_s[[key_col]], new_s[[key_col]]), "shared rows aligned on the key")

# TRUE where a cell differs (NA vs NA is equal, NA vs value is a change)
cell_differs <- function(x, y) {
  (is.na(x) != is.na(y)) | (!is.na(x) & !is.na(y) & x != y)
}

# TRUE where both values are numbers that are equal (e.g. "1.30" vs "1.3")
same_number <- function(x, y) {
  xn <- suppressWarnings(as.numeric(x))
  yn <- suppressWarnings(as.numeric(y))
  !is.na(xn) & !is.na(yn) & xn == yn
}

comparison <- rbindlist(lapply(cols_shared, function(col) {
  x <- old_s[[col]]
  y <- new_s[[col]]
  d <- cell_differs(x, y)
  data.table(
    column = col,
    n_changed = sum(d),
    pct_changed = round(100 * mean(d), 4),
    n_value_to_NA = sum(d & !is.na(x) & is.na(y)),
    n_NA_to_value = sum(d & is.na(x) & !is.na(y)),
    n_value_to_value = sum(d & !is.na(x) & !is.na(y)),
    n_format_only = sum(d & same_number(x, y)), # same number, written differently
    n_rows_changed_ids = uniqueN(old_s[[key_col]][d])
  )
}))
comparison[, result := fifelse(n_changed == 0, "Identical", "Different")]
setorder(comparison, -n_changed)

# >>> RESULT: columns that have changed <<<
changed_columns <- comparison[result == "Different", column]
changed_columns
comparison[result == "Different"]

# ---- 5. Detail: every changed cell ----
id_cols <- intersect(c("TreeID", "StemID", "Tag", "StemTag", "CensusID", "PlotCensusNumber"), names(new_s))

changes <- rbindlist(lapply(changed_columns, function(col) {
  d <- which(cell_differs(old_s[[col]], new_s[[col]]))
  cbind(
    new_s[d, c(key_col, id_cols), with = FALSE],
    data.table(column = col, old_value = old_s[[col]][d], new_value = new_s[[col]][d])
  )
}))
changes

# Most frequent old -> new changes per column (quick view of what happened)
changes_top <- changes[, .N, by = .(column, old_value, new_value)][order(column, -N)]
changes_top[, head(.SD, 10), by = column]

# Changed cells per column and census
if (nrow(changes) > 0) dcast(changes, column ~ PlotCensusNumber, value.var = "DBHID", fun.aggregate = length)

# ---- 6. Rows that exist only in one version ----
rows_only_old <- old[get(key_col) %in% ids_only_old]
rows_only_new <- new[get(key_col) %in% ids_only_new]

rows_only_old[, .N, keyby = .(CensusID, PlotCensusNumber, Status)]
rows_only_new[, .N, keyby = .(CensusID, PlotCensusNumber, Status)]

# ---- 7. Final check: compare DBH distributions using density ----
# DBH (and HOM) are part of the column-by-column comparison in Section 4.
# Here the DBH distributions of both versions are compared visually and with
# L1/L2 distances plus a Kolmogorov-Smirnov test.
dbh_cols <- intersect(c("DBH", "HOM"), cols_shared)
comparison[column %in% dbh_cols]

# DBH as numeric (the tables were loaded as character)
dbh_old_shared <- as.numeric(old_s$DBH) # old version, rows shared by both
dbh_new_shared <- as.numeric(new_s$DBH) # new version, rows shared by both
dbh_old_all <- as.numeric(old$DBH) # old version, all rows
dbh_new_all <- as.numeric(new$DBH) # new version, all rows (includes rows only in new)
dbh_only_new <- as.numeric(rows_only_new$DBH) # rows that exist only in the new version

# Measured DBH per version (rows with a DBH value) and largest cell difference
data.table(
  set = c("old (shared rows)", "new (shared rows)", "old (all rows)", "new (all rows)", "rows only in new"),
  n_rows = c(length(dbh_old_shared), length(dbh_new_shared), length(dbh_old_all), length(dbh_new_all), length(dbh_only_new)),
  n_with_DBH = c(sum(!is.na(dbh_old_shared)), sum(!is.na(dbh_new_shared)), sum(!is.na(dbh_old_all)), sum(!is.na(dbh_new_all)), sum(!is.na(dbh_only_new))),
  median_DBH = c(median(dbh_old_shared, na.rm = TRUE), median(dbh_new_shared, na.rm = TRUE), median(dbh_old_all, na.rm = TRUE), median(dbh_new_all, na.rm = TRUE), median(dbh_only_new, na.rm = TRUE)),
  max_DBH = c(max(dbh_old_shared, na.rm = TRUE), max(dbh_new_shared, na.rm = TRUE), max(dbh_old_all, na.rm = TRUE), max(dbh_new_all, na.rm = TRUE), max(dbh_only_new, na.rm = TRUE))
)
max(abs(dbh_old_shared - dbh_new_shared), na.rm = TRUE) # 0 = DBH identical on shared rows

# ---- Helper: compare_densities ----
# Purpose: Compute and plot density estimates for two numeric vectors and
#          report L1/L2 differences plus a Kolmogorov-Smirnov test.
# Args:
#   x1, x2 - numeric vectors (NA will be removed)
#   n_grid - grid size for density interpolation (default 2048)
#   plot_title - main title for the plots
# Returns: a list with L1/L2 distances, KS test object, density grids and
#          difference curve and summaries.
compare_densities <- function(x1, x2, n_grid = 2048, plot_title = "Density Comparison") {
  # Remove NAs
  x1 <- na.omit(x1)
  x2 <- na.omit(x2)
  # Compute densities
  d1 <- density(x1, n = n_grid)
  d2 <- density(x2, n = n_grid)
  # Align densities on same grid
  x_grid <- seq(min(c(d1$x, d2$x)), max(c(d1$x, d2$x)), length.out = n_grid)
  d1y <- approx(d1$x, d1$y, xout = x_grid, rule = 2)$y
  d2y <- approx(d2$x, d2$y, xout = x_grid, rule = 2)$y
  # Compute difference curves
  diff_curve <- d1y - d2y
  # Compute L1 and L2 differences
  library(pracma)
  L1_diff <- trapz(x_grid, abs(diff_curve))
  L2_diff <- sqrt(trapz(x_grid, diff_curve^2))
  # Kolmogorov-Smirnov test
  ks <- ks.test(x1, x2)
  # Plot densities and difference
  par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))
  plot(x_grid, d1y,
    type = "l", col = "blue", lwd = 2,
    main = paste(plot_title, "- Densities"), xlab = "Value", ylab = "Density"
  )
  lines(x_grid, d2y, col = "red", lwd = 2)
  legend("topright", legend = c("x1", "x2"), col = c("blue", "red"), lwd = 2)
  plot(x_grid, diff_curve,
    type = "l", col = "darkgreen", lwd = 2,
    main = "Difference Curve (x1 - x2)", xlab = "Value", ylab = "Difference"
  )
  abline(h = 0, col = "gray", lty = 2)
  # Return full report as a list
  report <- list(
    L1_difference = L1_diff,
    L2_difference = L2_diff,
    KS_test = ks,
    density_grid = x_grid,
    density_x1 = d1y,
    density_x2 = d2y,
    difference_curve = diff_curve,
    summary_x1 = summary(x1),
    summary_x2 = summary(x2)
  )
  return(report)
}

# ---- Helper: save_density_plot ----
# Purpose: Redraw the two panels of `compare_densities()` from its report and
#          save them as a PNG (no recomputation of the densities).
save_density_plot <- function(report, file, plot_title = "Density Comparison") {
  png(file, width = 1600, height = 1600, res = 200)
  par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))
  plot(report$density_grid, report$density_x1,
    type = "l", col = "blue", lwd = 2,
    main = paste(plot_title, "- Densities"), xlab = "Value", ylab = "Density"
  )
  lines(report$density_grid, report$density_x2, col = "red", lwd = 2)
  legend("topright", legend = c("x1", "x2"), col = c("blue", "red"), lwd = 2)
  plot(report$density_grid, report$difference_curve,
    type = "l", col = "darkgreen", lwd = 2,
    main = "Difference Curve (x1 - x2)", xlab = "Value", ylab = "Difference"
  )
  abline(h = 0, col = "gray", lty = 2)
  dev.off()
  invisible(file)
}

# Run density comparison and print summary statistics
# (x1 = old version, x2 = new version; each plot is also saved to DATA/CHECKS)

# a) Shared rows: old vs new DBH of the same measurements
report_densities_shared <- compare_densities(dbh_old_shared, dbh_new_shared, plot_title = "DBH old vs new (shared rows)")
save_density_plot(report_densities_shared, paste0(out_prefix, "_dbh_density_shared_rows.png"), plot_title = "DBH old vs new (shared rows)")
report_densities_shared$L1_difference
report_densities_shared$L2_difference
report_densities_shared$KS_test

# b) All rows: full old table vs full new table (new includes the added rows)
report_densities_all <- compare_densities(dbh_old_all, dbh_new_all, plot_title = "DBH old vs new (all rows)")
save_density_plot(report_densities_all, paste0(out_prefix, "_dbh_density_all_rows.png"), plot_title = "DBH old vs new (all rows)")
report_densities_all$L1_difference
report_densities_all$L2_difference
report_densities_all$KS_test

# c) Rows only in new vs the old table (are the added records comparable?)
report_densities_only_new <- compare_densities(dbh_old_all, dbh_only_new, plot_title = "DBH old vs rows only in new")
save_density_plot(report_densities_only_new, paste0(out_prefix, "_dbh_density_rows_only_new.png"), plot_title = "DBH old vs rows only in new")
report_densities_only_new$L1_difference
report_densities_only_new$L2_difference
report_densities_only_new$KS_test

# d) Same three comparisons on log10(DBH): DBH is strongly right-skewed, so on
#    the raw scale the curves collapse near zero; the log scale shows the shape.
report_log_shared <- compare_densities(log10(dbh_old_shared), log10(dbh_new_shared), plot_title = "log10 DBH old vs new (shared rows)")
save_density_plot(report_log_shared, paste0(out_prefix, "_dbh_density_log10_shared_rows.png"), plot_title = "log10 DBH old vs new (shared rows)")

report_log_all <- compare_densities(log10(dbh_old_all), log10(dbh_new_all), plot_title = "log10 DBH old vs new (all rows)")
save_density_plot(report_log_all, paste0(out_prefix, "_dbh_density_log10_all_rows.png"), plot_title = "log10 DBH old vs new (all rows)")

report_log_only_new <- compare_densities(log10(dbh_old_all), log10(dbh_only_new), plot_title = "log10 DBH old vs rows only in new")
save_density_plot(report_log_only_new, paste0(out_prefix, "_dbh_density_log10_rows_only_new.png"), plot_title = "log10 DBH old vs rows only in new")

# ---- 8. Save outputs ----
fwrite(structure_summary, paste0(out_prefix, "_structure.csv"))
fwrite(comparison, paste0(out_prefix, "_columns.csv"))
fwrite(changes, paste0(out_prefix, "_changed_cells.csv"))
fwrite(changes_top, paste0(out_prefix, "_changed_values_top.csv"))
fwrite(rows_only_old, paste0(out_prefix, "_rows_only_old.csv"))
fwrite(rows_only_new, paste0(out_prefix, "_rows_only_new.csv"))

cat(
  "Columns changed:", if (length(changed_columns)) paste(changed_columns, collapse = ", ") else "none", "\n",
  "Rows only in old:", length(ids_only_old), "| rows only in new:", length(ids_only_new), "\n",
  file = log_file, append = TRUE
)
cat("✓ Outputs written to", CHECKS_folder, "\n")
