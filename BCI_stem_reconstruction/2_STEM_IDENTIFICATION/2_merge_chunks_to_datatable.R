# =============================================================================
# 2_merge_chunks_to_datatable.R
#
# Purpose: Merge stem-identification DP Feather chunks into final datasets,
#          validate the merged multi-stem result, and reassemble the complete
#          dataset with reconstructed stem IDs.
# =============================================================================

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

# Root directory where the Feather chunk outputs from the DP run are stored.
# Update this path if the prior run folder has moved.
home_dir <- "/Users/medinaja/outputs_bci_stem_identification"

# Choose the run subfolder by its directory name.
run_code <- list.files(home_dir)

# Derived paths -----------------------------------------------------------------
chunks_path <- file.path(home_dir, run_code) # Feather input directory
outputs_path <- here("BCI_stem_reconstruction", "DATA", run_code) # merged output directory

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

# Collect Feather chunk files and verify the expected numeric suffix.
feathers <- list.files(chunks_path, "\\.feather$", full.names = TRUE)
cat("Total Feather files found:", length(feathers), "\n")

# Extract the embedded sequence number from each filename for correct ordering.
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
temp_dir <- here("BCI_stem_reconstruction", "DATA", "temp_parts")
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

# Write the merged result to Parquet and also save an RDS copy for R users.
write_parquet(ds_final, out_file,
    compression = "zstd", compression_level = 3, version = "2.6"
)
saveRDS(data.table::as.data.table(ds_final), file = out_file_rds)

cat("Merging complete!\n")
cat("Number of records in merged dataset:", nrow(data.table::as.data.table(ds_final)), "\n")

# Convert the Arrow dataset to a regular data.table for the remaining in-memory checks.
ds_final <- data.table::as.data.table(ds_final)

# Remove the intermediate Parquet parts — they are no longer needed now that
# the single merged Parquet file has been written to disk.
unlink(temp_dir, recursive = TRUE)

# =============================================================================
# 4. VALIDATION HELPER
# =============================================================================

# compare_columns() performs column-by-column equality checks between two
# data.tables using unified factor levels to handle mixed types and NAs
# correctly. Columns present in only one table are silently ignored.
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
    "BCI_stem_reconstruction", "DATA", "PROCESSED", "ViewFullTable_taper_corrected_growth_forms.rds"
)))
xraw[, growth_form := Lifeform]
xraw[, Lifeform := NULL]

# Validate that the DP output row count matches the multi-stem records in raw input.
cat("Row count breakdown by stem type in raw data:\n")
print(table(xraw$single_stem_tags))
cat("Row count in merged DP output:", nrow(ds_final), "\n")

# Extract Tag values for multi-stem records; used interactively to verify
# which tags the DP algorithm processed and cross-check against ds_final.
tags_multistem_raw <- xraw[single_stem_tags == FALSE]$Tag
cat("Number of multi-stem tag-census records in raw data:", length(tags_multistem_raw), "\n")
rm(tags_multistem_raw)

# Isolate multi-stem input records and sort by RowID for aligned comparison
input_multi_stem_data <- xraw[single_stem_tags == FALSE]
setorder(input_multi_stem_data, RowID)

# Verify the DP output's duplicate species column matches the original Mnemonic.
if (isTRUE(unique(ds_final$Mnemonic == ds_final$Species))) {
    cat("OK: 'Mnemonic' and 'species' columns are identical; dropping 'species'.\n")
    ds_final[, Species := NULL]
} else {
    stop("PROBLEM: 'Mnemonic' and 'species' differ — investigate before proceeding.")
}

# Sort the merged output to match the ordering of the raw multi-stem subset
setorder(ds_final, RowID)

# =============================================================================
# 6. FIRST VALIDATION PASS: DP OUTPUT VS. RAW MULTI-STEM DATA
# =============================================================================

# Known differences after DP output:
#   - DBH: internal taper-corrected values were used; original DBH is preserved in backup.
#   - ExactDate: some missing dates were imputed.
#   - growth_form: categories were simplified for internal DP parameter estimation.
# Any other differences should be investigated.

cat("\n--- Validation pass 1: DP output vs. raw multi-stem input ---\n")
comparison1 <- compare_columns(input_multi_stem_data, ds_final)
print(comparison1[result == "Different"])

# check growth forms
unique(merge(unique(input_multi_stem_data[, .(Tag, growth_form)]),
    unique(ds_final[, .(Tag, growth_form)]),
    by = "Tag", all = TRUE
)[, .(growth_form.x, growth_form.y)])

# =============================================================================
# 7. RESTORE ORIGINAL DBH AND RE-VALIDATE
# =============================================================================

# Restore original DBH values from DBH_mm_original_backup and drop the backup column.

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

# Keep the original raw table columns plus the DP-specific output columns.
original_names <- names(xraw)

# Retain only the active DP output columns; omit posterior-probability diagnostics.
new_names <- c(
    "TrueStemID", # definitive stem identifier assigned by the DP run
    "ReconstructedStemID", # stem ID reconstructed across historical censuses
    "SweepAuditOverride", # boolean flag for manual override of DP assignment based on sweep-audit evidence
    "ReconstructedStemID_PreSweep", # stem ID reconstructed by the DP before applying sweep-audit overrides
    "SweepRollbackToPreSweep", # boolean flag indicating whether the final ReconstructedStemID was rolled back to the pre-sweep version due to an override
    "ReconstructionMethod", # algorithm branch taken (e.g., "DP", "fallback")
    "DP_FallbackReason", # reason a fallback was triggered; NA when DP succeeded
    "obs_row_id", # sequential row number within each tag, assigned by the DP during processing
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
input_single_stem_data_raw <- xraw[single_stem_tags == TRUE]

# Verify that the single-stem subset contains one StemID per Tag.
chk_correctness_single <- copy(input_single_stem_data_raw)
chk_correctness_single[!is.na(CensusID), UniqueCensusIDs := uniqueN(CensusID), by = Tag]
chk_correctness_single[!is.na(StemID), UniqueStemID := uniqueN(StemID), by = Tag]
cat("\nDistribution of unique StemIDs per tag in the single-stem subset:\n")
print(table(chk_correctness_single$UniqueStemID, useNA = "ifany"))
# Expected: all tags have UniqueStemID == 1 (or NA for tags with no recorded StemID)

unique(chk_correctness_single[is.na(StemID)]$Tag)

# Simplified renumbering for single-stem tags: renames ReconstructedStemID from 1 to N per tag.
# Emits a warning if a tag has more than one unique stem.
renumber_single_stem_ids <- function(dt, verbose = TRUE) {
    if (is.null(dt) || nrow(dt) == 0L) {
        return(dt)
    }
    needed <- c("Tag", "ReconstructedStemID")
    if (!all(needed %in% names(dt))) {
        if (isTRUE(verbose)) message("[renumber_single_stem_ids] Required columns missing; skipping.")
        return(dt)
    }
    dt <- copy(dt)
    # Cast column to integer so assignments don't get coerced to logical
    dt[, ReconstructedStemID := as.integer(ReconstructedStemID)]
    # Populate ReconstructedStemID from StemID where it is NA
    if ("StemID" %in% names(dt)) {
        dt[
            is.na(ReconstructedStemID) & !is.na(StemID),
            ReconstructedStemID := as.integer(factor(StemID))
        ]
    }
    # Renumber sequentially within each Tag
    dt[!is.na(ReconstructedStemID),
        ReconstructedStemID := {
            stem_ids <- unique(ReconstructedStemID)
            n_stems <- length(stem_ids)
            if (n_stems > 1 && isTRUE(verbose)) {
                warning(sprintf(
                    "[renumber_single_stem_ids] Tag %s has %d unique non-NA stem IDs; expected 1. Renumbering anyway.",
                    .BY$Tag, n_stems
                ))
            }
            id_map <- setNames(seq_len(n_stems), as.character(stem_ids))
            as.integer(id_map[as.character(ReconstructedStemID)])
        },
        by = Tag
    ]
    dt
}

input_single_stem_data <- renumber_single_stem_ids(input_single_stem_data_raw)
input_single_stem_data[, .(Tag, CensusID, StemID, ReconstructedStemID)]

input_single_stem_data[is.na(ReconstructedStemID) & !is.na(StemID), .(Tag, CensusID, StemID, ReconstructedStemID)]
input_single_stem_data[is.na(ReconstructedStemID), .(Tag, CensusID, StemID, ReconstructedStemID)]
input_single_stem_data[Tag == "528106", .(Tag, CensusID, StemID, ReconstructedStemID)]

# Ensure ReconstructionMethod is character before assigning the label
# (it may have been read in as a different type if xraw had a factor column).
input_single_stem_data[, ReconstructionMethod := as.character(ReconstructionMethod)]
input_single_stem_data[, ReconstructionMethod := "single_stem_tag_no_reconstructed"]

# Verify that single-stem records with NA StemID also lack a DBH value
# (NA StemID with non-NA DBH would indicate an ambiguous placement).
unique(input_single_stem_data[is.na(StemID)]$DBH)

# Check reconstructed IDs assigned to single-stem tags.
input_single_stem_data[ReconstructedStemID > 1]
table(output_multistem_data$ReconstructedStemID, useNA = "ifany")

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

# =============================================================================
# 9. FINAL VALIDATION: COMPLETE DATASET VS. RAW TABLE (SHARED COLUMNS)
# =============================================================================

# The shared columns of complete_dataset should match xraw exactly:
#   - Single-stem rows were taken directly from xraw (after NA-filling new cols).
#   - Multi-stem rows had their original columns restored in Section 7.

cat("\n--- Final validation: complete dataset vs. raw table (shared columns) ---\n")
comparison_final <- compare_columns(complete_dataset, xraw)
print(comparison_final[result == "Different"])

# Differences should be limited to the new DP-specific columns and imputed ExactDate values.

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

# Check that rows with reconstructed IDs also have a stem ID indicator.
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
# For skipped_no_data rows, use StemID to reconstruct missing ReconstructedStemID values.
# These rows had no DP assignment but do have valid StemID values.

chk_skipped_no_data_tags <- unique(complete_dataset[ReconstructionMethod == "skipped_no_data"]$Tag)
complete_dataset[Tag %in% chk_skipped_no_data_tags, .(Tag, CensusID, ReconstructedStemID, StemID, DBH, ReconstructionMethod)][order(Tag, CensusID)]

complete_dataset[
    Tag %in% chk_skipped_no_data_tags & is.na(ReconstructedStemID)
]

complete_dataset[is.na(ReconstructedStemID) & !is.na(StemID), .(Tag, CensusID, ReconstructedStemID, StemID, DBH, ReconstructionMethod)]

complete_dataset_final_raw <- complete_dataset[!ReconstructionMethod %in% "skipped_no_data"]

skipped_no_data <- complete_dataset[ReconstructionMethod %in% "skipped_no_data"]
skipped_no_data <- renumber_single_stem_ids(skipped_no_data)

complete_dataset_final <- rbind(complete_dataset_final_raw, skipped_no_data)
setorder(complete_dataset_final, RowID)

# Sanity: skipped tags (no DBH/CensusID data) are expected to have NA obs_row_id.
complete_dataset_final[is.na(obs_row_id) & single_stem_tags == FALSE]

# Review rows that still lack a ReconstructionMethod.
table(complete_dataset_final$ReconstructionMethod, useNA = "ifany")

complete_dataset_final[is.na(ReconstructionMethod)]
# These rows are missing DBH and/or ReconstructionMethod; they require follow-up review.

chk_skipped_full_na <- unique(complete_dataset_final[is.na(ReconstructionMethod)]$Tag)
print(complete_dataset_final[Tag %in% chk_skipped_full_na, .(Tag, CensusID, ReconstructedStemID, StemID, DBH, ReconstructionMethod, ListOfTSM)][order(Tag, CensusID)], nrows = 200)
# These rows correspond to dead tags appearing later as R/broken-below.
# They are expected to be corrected during subsequent RTable reconstruction.

# Output directory for the final complete dataset.
post_dir <- path.expand(
    file.path(
        "./BCI_stem_reconstruction/DATA/PROCESSED"
    )
)

# Ensure output directory exists
if (!dir.exists(post_dir)) {
    dir.create(post_dir, recursive = TRUE, showWarnings = TRUE)
}

saveRDS(
    complete_dataset_final,
    here(post_dir, "complete_dataset_final_with_reconstructed_stemids.rds")
)
