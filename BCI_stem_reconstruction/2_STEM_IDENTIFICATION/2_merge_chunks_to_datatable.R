## 2_merge_chunks_to_datatable.R
## =============================================================================
## Purpose:
##   Combine individual Feather chunks produced during the stem-identification
##   DP (dynamic-programming) run into a single merged dataset and save it as
##   both Parquet and RDS for downstream analysis.
##
##   Processing pipeline:
##     1. Read Feather chunk files in batches to limit peak memory usage.
##     2. Write each batch as an intermediate Parquet part.
##     3. Merge all Parquet parts into one final Parquet + RDS output.
##     4. Validate the merged multi-stem result against the raw input table.
##     5. Reattach single-stem records and save the complete dataset.
##
## Inputs:
##   - Feather chunks       : <home_dir>/<run_code>/
##   - Raw processed table  : DATA/PROCESSED/6_ViewFullTable_taper_corrected_growth_forms.rds
##
## Outputs:
##   - 2_STEM_IDENTIFICATION/output/<run_code>/merged_output.parquet
##   - 2_STEM_IDENTIFICATION/output/<run_code>/merged_output.rds
##   - DATA/PROCESSED/7_complete_dataset_with_reconstructed_stemids.rds
## =============================================================================

# =============================================================================
# 0. SETUP
# =============================================================================

# Clear workspace to ensure no leftover objects affect the run
rm(list = ls())

# Packages
library(arrow) # read/write Feather/Parquet and manage Arrow datasets
library(here) # construct file paths relative to the project root
library(data.table) # fast in-memory data manipulation

# =============================================================================
# 1. CONFIGURATION
# =============================================================================

# Root directory on the local machine where the Feather chunk files are stored.
# This folder is outside the project tree; update the path if you relocate the
# output from a previous DP run.
home_dir <- "/Users/medinaja/outputs_bci_stem_identification"

# Select the run subfolder by index from the list of directories in `home_dir`.
# Run list.files(home_dir) interactively to inspect available runs and confirm
# the correct index BEFORE executing the rest of this script.
run_code <- list.files(home_dir)[1]

# Derived paths -----------------------------------------------------------------
chunks_path <- file.path(home_dir, run_code) # Feather input directory
outputs_path <- here("2_STEM_IDENTIFICATION", "output", run_code) # merged output directory

# Create the output directory if it does not already exist
if (!dir.exists(outputs_path)) {
    dir.create(outputs_path, recursive = TRUE)
}

# Final output file names
out_file <- here(outputs_path, "merged_output.parquet")
out_file_rds <- here(outputs_path, "merged_output.rds")

# =============================================================================
# 2. BATCH-PROCESS FEATHER CHUNKS INTO PARQUET PARTS
# =============================================================================

# Collect the full list of Feather files produced by the DP run
feathers <- list.files(chunks_path, "\\.feather$", full.names = TRUE)
cat("Total Feather files found:", length(feathers), "\n")

# Extract the numeric sequence number (1–4 digits) embedded in each file name.
# This serves two purposes:
#   (a) Validates that every file matches the expected naming convention.
#   (b) Provides the sort key so chunks are processed in the correct order,
#       regardless of how list.files() orders them on the filesystem.
feathers_seq <- as.integer(gsub(".*_(\\d{1,4})\\.feather$", "\\1", feathers))
if (any(is.na(feathers_seq))) {
    stop("Error: Some Feather files do not match the expected naming pattern.")
}

# Sort files by their embedded sequence number to guarantee correct ordering
# before the batched Parquet conversion and subsequent merge.
feathers <- feathers[order(feathers_seq)]
feathers_seq <- feathers_seq[order(feathers_seq)]

# Split the sorted file list into fixed-size batches to cap peak memory consumption
group_size <- 500 # files per batch
n_groups <- ceiling(length(feathers) / group_size)

# Temporary directory for intermediate Parquet parts (removed after final merge)
temp_dir <- here("2_STEM_IDENTIFICATION", "output", "temp_parts")
dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)

# Process each batch sequentially:
#   1. Compute the index range for the current batch.
#   2. Slice the sorted Feather list to get this batch's files.
#   3. Open as a lazy Arrow dataset — no data loaded into RAM until written.
#   4. Write the batch to a zero-padded, numbered Parquet part file.
# ZSTD level-3 compression offers a good balance between file size and speed.
for (i in seq_len(n_groups)) {
    cat("Processing group", i, "of", n_groups, "\n")
    # Compute the index range for this batch
    start_idx <- (i - 1) * group_size + 1
    end_idx <- min(i * group_size, length(feathers))
    group_files <- feathers[start_idx:end_idx]
    # open_dataset streams lazily; actual reads happen only when writing to disk
    ds_group <- open_dataset(group_files, format = "feather", unify_schemas = FALSE)
    # Zero-padded part name ensures correct lexicographic ordering when re-read
    temp_file <- file.path(temp_dir, sprintf("part_%03d.parquet", i))
    write_parquet(ds_group, temp_file,
        compression = "zstd", compression_level = 3, version = "2.6"
    )
}

# =============================================================================
# 3. MERGE PARQUET PARTS INTO FINAL DATASET
# =============================================================================

# Open all intermediate parts as one unified Arrow dataset
temp_files <- list.files(temp_dir, "\\.parquet$", full.names = TRUE)
ds_final <- open_dataset(temp_files, format = "parquet", unify_schemas = FALSE)

# Persist the merged dataset as Parquet (primary output)
write_parquet(ds_final, out_file,
    compression = "zstd", compression_level = 3, version = "2.6"
)

# Also save as RDS so downstream R scripts can load the result without Arrow
saveRDS(data.table::as.data.table(ds_final), file = out_file_rds)

cat("Merging complete!\n")
cat("Number of records in merged dataset:", nrow(data.table::as.data.table(ds_final)), "\n")

# Materialise the Arrow dataset into a data.table for all in-memory operations
# that follow (column inspection, row- binding, validation).  After this point
# ds_final is a regular data.table and Arrow is no longer needed.
ds_final <- data.table::as.data.table(ds_final)

# Remove the intermediate Parquet parts — they are no longer needed now that
# the single merged Parquet file has been written to disk.
unlink(temp_dir, recursive = TRUE)

chk <- ds_final[, .(TreeID, Tag, StemID, TrueStemID, ReconstructedStemID, StemTag, CensusID, DBH, obs_row_id)][!is.na(StemTag)]

unique(chk[StemID != ReconstructedStemID]$CensusID)

# =============================================================================
# 4. VALIDATION HELPER
# =============================================================================

# compare_columns() performs a column-by-column equality check between two
# data.tables.  Strategy:
#   1. Find the set of column names present in both tables.
#   2. Coerce each column pair to factors with a unified level set, which
#      correctly handles mixed types (e.g. character vs. numeric) and
#      explicit NA values (addNA() adds NA as a proper factor level).
#   3. Compare integer codes of the aligned factors; if every code matches,
#      the column is "Similar", otherwise "Different".
# Columns present in only one table are silently ignored.
# Returns: data.table with columns 'column' and 'result' (one row per shared column).
compare_columns <- function(dt1, dt2) {
    cols <- intersect(names(dt1), names(dt2))
    rbindlist(lapply(cols, function(col) {
        v1 <- addNA(as.factor(dt1[[col]]))
        v2 <- addNA(as.factor(dt2[[col]]))
        all_lvls <- union(levels(v1), levels(v2))
        f1 <- factor(v1, levels = all_lvls)
        f2 <- factor(v2, levels = all_lvls)
        data.table(
            column = col,
            result = if (identical(levels(f1), levels(f2)) &&
                identical(as.integer(f1), as.integer(f2))) {
                "Similar"
            } else {
                "Different"
            }
        )
    }))
}

# =============================================================================
# 5. LOAD RAW INPUT AND PREPARE MULTI-STEM SUBSET FOR VALIDATION
# =============================================================================

# The raw table produced in step 6 of the data-preparation pipeline; it
# contains both single-stem and multi-stem records and is the reference for
# all validation steps below.
xraw <- as.data.table(readRDS(here(
    "DATA", "PROCESSED", "6_ViewFullTable_taper_corrected_growth_forms.rds"
)))
xraw[, growth_form := Lifeform]
xraw[, Lifeform := NULL]

# Sanity check: the DP algorithm was run exclusively on multi-stem records
# (single_stem_tags == FALSE). Therefore ds_final must contain exactly as many
# rows as there are multi-stem records in xraw.  A mismatch here means some
# chunks were missing, duplicated, or the wrong run_code was selected above.
cat("Row count breakdown by stem type in raw data:\n")
print(table(xraw$single_stem_tags))
cat("Row count in merged DP output:", nrow(ds_final), "\n")

# Extract Tag values for multi-stem records; used interactively to verify
# which tags the DP algorithm processed and cross-check against ds_final.
tags_multistem_raw <- xraw[single_stem_tags == FALSE]$Tag
cat("Number of multi-stem tag-census records in raw data:", length(tags_multistem_raw), "\n")

# Isolate multi-stem input records and sort by RowID for aligned comparison
input_multi_stem_data <- xraw[single_stem_tags == FALSE]
setorder(input_multi_stem_data, RowID)

# The DP algorithm stores the species mnemonic in a column named 'species',
# which duplicates the original 'Mnemonic' column from the raw input.
# Verify they are identical for every row before dropping the redundant copy.
# NOTE: isTRUE(unique(...)) returns TRUE only when the comparison vector
# contains exactly one distinct value and that value is TRUE (all rows agree).
if (isTRUE(unique(ds_final$Mnemonic == ds_final$species))) {
    cat("OK: 'Mnemonic' and 'species' columns are identical; dropping 'species'.\n")
    ds_final[, species := NULL]
} else {
    stop("PROBLEM: 'Mnemonic' and 'species' differ — investigate before proceeding.")
}

# Sort the merged output to match the ordering of the raw multi-stem subset
setorder(ds_final, RowID)

length(unique(ds_final$species))

# =============================================================================
# 6. FIRST VALIDATION PASS: DP OUTPUT VS. RAW MULTI-STEM DATA
# =============================================================================

# Expected columns reported as "Different" after the DP run:
#   - DBH       : the DP algorithm used a taper-corrected value internally;
#                 the original measurement is preserved in DBH_mm_original_backup
#                 and will be restored in Section 7.
#   - ExactDate : missing measurement dates were imputed during the DP run.
# Any column other than these two appearing as "Different" warrants investigation
# before proceeding.

cat("\n--- Validation pass 1: DP output vs. raw multi-stem input ---\n")
comparison1 <- compare_columns(input_multi_stem_data, ds_final)
print(comparison1[result == "Different"])

# =============================================================================
# 7. RESTORE ORIGINAL DBH AND RE-VALIDATE
# =============================================================================

# Replace the algorithm's internal (taper-corrected) DBH with the original
# field measurement stored in DBH_mm_original_backup, then discard the backup
# column.  The chained [:=] calls modify ds_final in place in two steps:
#   step 1: overwrite DBH with DBH_mm_original_backup values.
#   step 2: remove the now-redundant DBH_mm_original_backup column.

if ("DBH_mm_original_backup" %in% names(ds_final)) {
    ds_final[, DBH := DBH_mm_original_backup][, DBH_mm_original_backup := NULL]
    cat("Original DBH values restored in ds_final.\n")
}

cat("\n--- Validation pass 2: after restoring original DBH ---\n")
comparison2 <- compare_columns(input_multi_stem_data, ds_final)
print(comparison2[result == "Different"])
# Expected: only ExactDate remains "Different" (imputed NAs during DP run).

# =============================================================================
# 8. ASSEMBLE COMPLETE DATASET (SINGLE-STEM + MULTI-STEM)
# =============================================================================

ds_final
names(ds_final)

# Column names from the original raw table
original_names <- names(xraw)

head(ds_final, 2)

# Columns added by the DP stem-identification algorithm that do not exist in
# the original raw table.  Only the four active columns are retained in the
# final output; posterior-probability columns are excluded to reduce file size.
# To include them, uncomment the relevant lines and ensure the rbind schemas
# in the single-stem and multi-stem blocks are updated accordingly.
new_names <- c(
    "TrueStemID", # definitive stem identifier assigned by the DP run
    "ReconstructedStemID", # stem ID reconstructed across historical censuses
    "SweepAuditOverride", # boolean flag for manual override of DP assignment based on sweep-audit evidence
    "ReconstructedStemID_PreSweep", # stem ID reconstructed by the DP before applying sweep-audit overrides
    "SweepRollbackToPreSweep", # boolean flag indicating whether the final ReconstructedStemID was rolled back to the pre-sweep version due to an override
    "ReconstructionMethod", # algorithm branch taken (e.g., "DP", "fallback")
    "DP_FallbackReason", # reason a fallback was triggered; NA when DP succeeded
    "obs_row_id", # sequential row number within each tag, assigned by the DP during processing
    # "DP_PosteriorTop1ID" ,    # top posterior assignment: stem ID
    # "DP_PosteriorTop1Prob" # ,  # top posterior assignment: posterior probability
    # "DP_PosteriorTop2ID",    # second-best posterior assignment: stem ID
    # "DP_PosteriorTop2Prob"   # second-best posterior assignment: probability
    "DP_PosteriorReconstructedProb"
)

# Select only the original columns plus the new DP-specific columns from the
# merged output (drops any intermediate / diagnostic columns from the DP run)
original_new_names <- c(original_names, new_names)
output_multistem_data <- ds_final[, ..original_new_names]
names(ds_final)

# Add the four new DP columns to xraw with NA values so that both the
# single-stem and multi-stem subsets share an identical column schema before
# row-binding.  This mutates xraw in place.
xraw[, (new_names) := NA]

# Isolate single-stem records from the raw table
input_single_stem_data <- xraw[single_stem_tags == TRUE]

## Verify that the single-stem subset is genuinely single-stem:
## every Tag should have exactly one unique StemID across all censuses.
## A UniqueStemID > 1 for any tag indicates a mis-classification during
## the single_stem_tags assignment step and must be investigated.
chk_correctness_single <- copy(input_single_stem_data)
chk_correctness_single[!is.na(CensusID), UniqueCensusIDs := uniqueN(CensusID), by = Tag]
chk_correctness_single[!is.na(StemID), UniqueStemID := uniqueN(StemID), by = Tag]
cat("\nDistribution of unique StemIDs per tag in the single-stem subset:\n")
print(table(chk_correctness_single$UniqueStemID))
# Expected: all tags have UniqueStemID == 1 (or NA for tags with no recorded StemID)

# For single-stem trees there is only one stem, so StemID from the original
# data serves directly as the ReconstructedStemID — no reconstruction needed.
input_single_stem_data[, ReconstructedStemID := StemID]
# Ensure ReconstructionMethod is character before assigning the label
# (it may have been read in as a different type if xraw had a factor column).
input_single_stem_data[, ReconstructionMethod := as.character(ReconstructionMethod)]
input_single_stem_data[, ReconstructionMethod := "single_stem_tag_no_reconstructed"]

input_single_stem_data[is.na(ReconstructedStemID)]

# Verify that single-stem records with NA StemID also lack a DBH value.
# If any NA-StemID record has a non-NA DBH, it would represent an ambiguous
# placement that needs investigation before the final merge.
unique(input_single_stem_data[is.na(StemID)]$DBH)
# Expected: a single NA — confirms all NA-StemID records also have no DBH measurement.

# Combine single-stem (from raw) and multi-stem (from DP output) into one table
complete_dataset <- rbind(input_single_stem_data, output_multistem_data)
complete_dataset[, ReconstructedStemID := as.character(ReconstructedStemID)]
setorder(complete_dataset, RowID) # restore original row ordering
setorder(xraw, RowID)

if (nrow(complete_dataset) == nrow(xraw)) {
    cat("\nRow count check passed: complete dataset has the same number of rows as the raw table.\n")
} else {
    stop("\nRow count mismatch: complete dataset has ", nrow(complete_dataset), " rows, but raw table has ", nrow(xraw), " rows.")
}

inc <- unique(output_multistem_data[is.na(DBH), ReconstructedStemID])
output_multistem_data[ReconstructedStemID %in% inc[5], .(Tag, StemID, TrueStemID, ReconstructedStemID, StemTag, CensusID, DBH)]
rm(inc)

# =============================================================================
# 9. FINAL VALIDATION: COMPLETE DATASET VS. RAW TABLE (SHARED COLUMNS)
# =============================================================================

# The shared columns of complete_dataset should match xraw exactly:
#   - Single-stem rows were taken directly from xraw (after NA-filling new cols).
#   - Multi-stem rows had their original columns restored in Section 7.

cat("\n--- Final validation: complete dataset vs. raw table (shared columns) ---\n")
comparison_final <- compare_columns(complete_dataset, xraw)
print(comparison_final[result == "Different"])

# those different values are the new columns I addedd and the exactdate

# =============================================================================
# 10. CHECK RECONSTRUCTEDSTEMID ASSIGNMENTS
# =============================================================================
table(complete_dataset$ReconstructionMethod, useNA = "ifany")

check_recstemid <- complete_dataset[ReconstructionMethod != "single_stem_tag_no_reconstructed" & !is.na(DBH), .(CensusID, Tag, StemTag, StemID, DBH, TrueStemID, ReconstructedStemID, SweepAuditOverride, ReconstructedStemID_PreSweep, ReconstructionMethod, DP_FallbackReason, Status)]

check_recstemid[, c("StemID_01", "TrueStemID_01", "ReconstructedStemID_01", "ReconstructedStemID_PreSweep_01") := list(
    ifelse(is.na(StemID), 0, 1),
    ifelse(is.na(TrueStemID), 0, 1),
    ifelse(is.na(ReconstructedStemID), 0, 1),
    ifelse(is.na(ReconstructedStemID_PreSweep), 0, 1)
)]

# NOTE: All reconstructed stemid have one value
check_recstemid[ReconstructedStemID_01 != StemID_01]

# =============================================================================
# 11. SAVE FINAL COMPLETE DATASET
# =============================================================================

# Print a concise summary of the final dataset before writing to disk.
cat("\nComplete dataset summary before saving:\n")
cat("  Total rows                :", nrow(complete_dataset), "\n")
cat("  Single-stem rows          :", nrow(complete_dataset[single_stem_tags == TRUE]), "\n")
cat("  Multi-stem rows           :", nrow(complete_dataset[single_stem_tags == FALSE]), "\n")
cat("  Total columns             :", ncol(complete_dataset), "\n")
cat("  ReconstructionMethod dist.:\n")
print(table(complete_dataset$ReconstructionMethod, useNA = "ifany"))

table(complete_dataset$ReconstructionMethod, useNA = "ifany")
# NOTE: There are 47 rows with ReconstructionMethod == "skipped_no_data"
# which means that DP did not assign a ReconstructedStemID to these rows because they had no data
# Fill in ReconstructedStemID for these rows using StemID

chk_skipped_no_data_tags <- unique(complete_dataset[ReconstructionMethod == "skipped_no_data"]$Tag)
complete_dataset[Tag %in% chk_skipped_no_data_tags, .(Tag, CensusID, ReconstructedStemID, StemID, DBH, ReconstructionMethod)][order(Tag, CensusID)]
# NOTE: all of them in census 9 and all have NA in ReconstructedStemID but not
# in StemID, so we can fill them with StemID without introducing any errors

complete_dataset[Tag %in% chk_skipped_no_data_tags & is.na(ReconstructedStemID), ReconstructedStemID := StemID]

complete_dataset[is.na(obs_row_id) & single_stem_tags == FALSE]
# all this are skipped tags because they didnt have any data (no DBH, no
# censusid, etc) and were not processed by the DP algorithm. This is expected.

# We refilled the ReconstructedStemID for the skipped_no_data rows with the original StemID, so now there should be no NA values in ReconstructedStemID for the multi-stem subset.

name_out_directory <- list.files(home_dir)[3]

# Directory containing posterior feather files from the stem identification run
# (expand tilde to user home directory for portability)
post_dir <- path.expand(
    file.path(
        "./DATA/PROCESSED",
        name_out_directory
    )
)

# Ensure output directory exists
if (!dir.exists(post_dir)) {
    dir.create(post_dir, recursive = TRUE, showWarnings = TRUE)
}

saveRDS(
    complete_dataset,
    here(post_dir, "7_complete_dataset_with_reconstructed_stemids.rds")
)

cat(
    "\nDone. Final dataset saved to:\n  ",
    here("DATA", "PROCESSED", "7_complete_dataset_with_reconstructed_stemids.rds"), "\n"
)
