## 0_PREPARE_POSTERIORS.R
## =============================================================================
## Purpose:
##   Aggregate posterior samples from individual tree stem identification runs
##   into a single consolidated RDS file. Each feather file contains MCMC chain
##   output (stem-identity paths) for one tree tag; this script reads all such
##   files, adds the tag as an identifier column, combines them into a single
##   data.table, and saves to disk for downstream analysis.
##
##   Input:  Directory of feather files with suffix '_paths.feather', each named
##           as 'tag_XXXXXX_*_paths.feather' where XXXXXX is a 6-digit tag ID.
##   Output: Single RDS file containing all posterior samples with Tag column.
##
## Processing pipeline:
##   1. Validate package availability (arrow, data.table).
##   2. Extract tag IDs from feather filenames using regex.
##   3. Read each feather file as data.table and append tag column.
##   4. Combine all data.tables using fast rbindlist.
##   5. Remove duplicate rows if present.
##   6. Save consolidated table to RDS.
## =============================================================================

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

home_dir <- "/Users/medinaja/outputs_bci_stem_identification"

# Select the run subfolder by index from the list of directories in `home_dir`.
# Run list.files(home_dir) interactively to inspect available runs and confirm
# the correct index BEFORE executing the rest of this script.
run_code <- list.files(home_dir)[2]

# Directory containing posterior feather files from the stem identification run
# (expand tilde to user home directory for portability)
post_dir <- path.expand(
    file.path(
        home_dir,
        run_code,
        "posteriors"
    )
)

# Output directory where consolidated posteriors will be saved
output_dir_posteriors <- "./BCI_stem_reconstruction/DATA/POSTERIORS"

# Ensure run-specific output directory exists
if (!dir.exists(output_dir_posteriors)) {
    dir.create(output_dir_posteriors, recursive = TRUE, showWarnings = TRUE)
}

# =============================================================================
# 1. HELPER FUNCTION: Read and Bind Feather Files
# =============================================================================

# read_and_bind_feathers() ------------------------------------------------
# Reads a collection of Apache Arrow feather files, extracts tree tag IDs from
# filenames, and combines them into a single data.table with the tag as an
# explicit column. Designed for posterior samples from the stem-identification
# pipeline, where each file contains MCMC chain output for one tree.
#
# Parameters:
#   file_paths : character vector of full paths to feather files;
#                filenames should match pattern 'tag_XXXXXX_*_paths.feather'.
# Returns:
#   data.table with columns: Tag (character), plus all columns from the feather
#   files (typically including MCMC chain index, stem-identity sequences, etc.).
#   Rows with duplicate values across all columns are NOT removed by this
#   function; deduplication can be applied separately if needed.
#   The table is keyed on Tag for improved join/filter performance.
read_and_bind_feathers <- function(file_paths) {
    # --- Input validation ---
    if (!is.character(file_paths) || length(file_paths) == 0) {
        stop("file_paths must be a non-empty character vector")
    }
    # ---- Extract tag IDs from filenames ---
    # Filename pattern: 'tag_XXXXXX_...', extract the 6 digits as the tag
    tags <- sub(pattern = ".*(\\d{6}).*", replacement = "\\1", x = basename(file_paths))
    # Verify extraction succeeded (all tags should match the pattern)
    if (anyNA(nchar(tags)) || any(nchar(tags) != 6)) {
        warning("Failed to extract valid tags from some filenames; check filename format.")
    }
    # ---- Read feather files and add tag column ---
    dt_list <- lapply(seq_along(file_paths), function(i) {
        dt <- as.data.table(arrow::read_feather(file_paths[i]))
        dt[, Tag := tags[i]]
        return(dt)
    })
    # ---- Combine all data.tables using fast rbindlist ---
    result <- rbindlist(dt_list, use.names = TRUE, fill = TRUE)
    # ---- Reorder columns with Tag first for readability and convenience ---
    setcolorder(result, c("Tag", setdiff(names(result), "Tag")))
    # ---- Set key on Tag for improved query and join performance ---
    setkey(result, "Tag")
    return(result)
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

## check size in mb for dt_posteriors
dt_size_mb <- object.size(dt_posteriors) / (1024^2)
cat("Size of consolidated data.table:", round(dt_size_mb, 2), "MB\n")

head(dt_posteriors)

# Remove duplicate rows if any posteriors have been processed multiple times
dt_posteriors <- unique(dt_posteriors)

cat(
    "Consolidated table dimensions:", nrow(dt_posteriors), "rows ×",
    ncol(dt_posteriors), "columns\n"
)

# =============================================================================
# 3. OUTPUT
# =============================================================================

# Save consolidated posteriors to RDS for fast loading and downstream analysis
output_file <- file.path(output_dir_posteriors, "posterior_sampled_paths.rds")
saveRDS(dt_posteriors, file = output_file)

cat("Posterior samples saved to:", output_file, "\n")
cat("Done.\n")
