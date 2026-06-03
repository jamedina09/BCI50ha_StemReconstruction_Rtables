# =============================================================================
# 1_prepare_posteriors_BCI.R
#
# Purpose: Consolidate posterior Feather outputs from a BCI DP run into one
#          RDS file with an explicit Tag column.
# =============================================================================

# =============================================================================
# SETUP
# =============================================================================

rm(list = ls())

# Load required packages (will error if not available; install before running)
library(arrow) # read_feather()
library(data.table) # fast data manipulation, rbindlist()

# =============================================================================
# CONFIGURATION
# =============================================================================

workspace_root <- getwd()
home_dir <- "/Users/medinaja/outputs_bci_stem_identification"

# Choose the run folder inside the output directory and locate its posterior files.
run_code <- list.files(home_dir)
post_dir <- normalizePath(file.path(home_dir, run_code, "posteriors"), winslash = "/", mustWork = FALSE)

# Output directory where consolidated posteriors will be saved
output_dir_posteriors <- file.path(workspace_root, "BCI_stem_reconstruction", "DATA", "POSTERIORS")

# Ensure run-specific output directory exists
if (!dir.exists(output_dir_posteriors)) {
    dir.create(output_dir_posteriors, recursive = TRUE, showWarnings = TRUE)
}

# =============================================================================
# 1. HELPER FUNCTION
# =============================================================================

# read_and_bind_feathers() ------------------------------------------------
# Read posterior Feather files, extract Tag from the filename, and combine
# all files into a single keyed data.table.
read_and_bind_feathers <- function(file_paths) {
    if (!is.character(file_paths) || length(file_paths) == 0) {
        stop("file_paths must be a non-empty character vector")
    }
    tags <- sub(".*(\\d{6}).*", "\\1", basename(file_paths))
    if (anyNA(nchar(tags)) || any(nchar(tags) != 6)) {
        warning("Failed to extract valid tags from some filenames; check filename format.")
    }
    dt_list <- lapply(seq_along(file_paths), function(i) {
        dt <- as.data.table(arrow::read_feather(file_paths[i]))
        dt[, Tag := tags[i]]
        dt
    })
    result <- rbindlist(dt_list, use.names = TRUE, fill = TRUE)
    setcolorder(result, c("Tag", setdiff(names(result), "Tag")))
    setkey(result, "Tag")
    result
}

# =============================================================================
# 2. MAIN EXECUTION
# =============================================================================

# Locate all feather files matching the posterior filename pattern
post_files <- list.files(
    post_dir,
    pattern = "_paths\\.feather$",
    full.names = TRUE
)

# Validate that at least one file was found
if (length(post_files) == 0) {
    stop("No posterior feather files found in directory: ", post_dir)
}

cat("Found", length(post_files), "posterior files\n")

# Read and bind all posterior files into a single data.table
dt_posteriors <- read_and_bind_feathers(post_files)

# link treeid to tag
tag_treeid_map <- as.data.table(readRDS(file.path(
    workspace_root, "BCI_stem_reconstruction", "DATA", "PROCESSED", "ViewFullTable_single_vs_multiple_stem_tags.rds"
)))
tag_treeid_map <- unique(tag_treeid_map[, .(TreeID, Tag)])
tag_treeid_map[, Tag := as.character(Tag)]
tag_treeid_map[, TreeID := as.character(TreeID)]

tag_post <- data.table(Tag = unique(dt_posteriors$Tag))

weird_tags <- unique(setdiff(tag_post$Tag, tag_treeid_map$Tag))

dt_posteriors <- merge(
    dt_posteriors,
    tag_treeid_map,
    by = "Tag",
    all.x = TRUE
)

col_order <- c("Tag", "TreeID", "path_sig", "path_count", "path_prob", "recon")

dt_posteriors <- dt_posteriors[, ..col_order][!is.na(TreeID)]

# rename Tag to tag and TreeID to treeID for consistency with the rest of the codebase
setnames(dt_posteriors, old = c("Tag", "TreeID"), new = c("tag", "treeID"))

# Report the consolidated table size and preview the top rows.
dt_size_mb <- object.size(dt_posteriors) / (1024^2)
cat("Size of consolidated data.table:", round(dt_size_mb, 2), "MB\n")
head(dt_posteriors)

# Deduplicate if duplicate rows were introduced during file binding.
dt_posteriors <- unique(dt_posteriors)

inspectdf::inspect_na(dt_posteriors)

cat(
    "Consolidated table dimensions:", nrow(dt_posteriors), "rows ×",
    ncol(dt_posteriors), "columns\n"
)

# =============================================================================
# 3. OUTPUT
# =============================================================================

output_file <- file.path(output_dir_posteriors, "posterior_sampled_paths.rds")
saveRDS(dt_posteriors, file = output_file)

cat("Posterior samples saved to:", output_file, "\n")
cat("Done.\n")
