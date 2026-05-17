## 3_compare_dryad_and_me.R
## =============================================================================
## Purpose:
##   Evaluate the quality of stem-identity reconstruction by comparing the
##   dataset produced in step 7 of the pipeline (reconstructed StemIDs) against
##   the Dryad/Condit reference measurements for BCI.
##
##   The comparison is based on cumulative per-stem DBH trajectories:
##   if a tag's reconstructed stems rank in the same order and sum to the same
##   total DBH as in the Dryad data, the reconstruction is considered correct.
##
##   Processing pipeline:
##     1. Load the reconstructed dataset and the Dryad reference table.
##     2. Summarize total DBH per tag × stem ID for both datasets (censuses 1–8).
##     3. Compare trajectory summaries to classify each tag as "Match" / "Mismatch".
##     4. Annotate tags with their reconstruction method (DP, probabilistic, both).
##     5. Split mismatched tags by method and draw a random 100-tag sample each.
##     6. Render diagnostic PDFs with side-by-side growth curves and data tables.
##
## Inputs:
##   - Reconstructed dataset : DATA/PROCESSED/7_complete_dataset_with_reconstructed_stemids.rds
##   - Dryad reference       : DATA/DRYAD_CONDIT/bci_dryad_condit.rds
##
## Outputs:
##   - 2_STEM_IDENTIFICATION/differences_using_dp.pdf
##   - 2_STEM_IDENTIFICATION/differences_using_probabilistic.pdf
##   - 2_STEM_IDENTIFICATION/differences_using_both_methods.pdf
## =============================================================================

# =============================================================================
# 0. SETUP
# =============================================================================

# Clear workspace to ensure no leftover objects affect the run
rm(list = ls())

# Packages
library(data.table) # fast in-memory data manipulation
library(here) # construct file paths relative to the project root
library(plotrix) # addtable2plot(): embed data tables inside base-R plot panels

# =============================================================================
# 1. HELPER FUNCTIONS
# =============================================================================

# load_dryad() -----------------------------------------------------------------
# Reads the BCI Dryad/Condit RDS file, converts DBH from mm (as stored) to cm,
# and harmonises column names to match the reconstructed dataset:
#   - DBH  : numeric, diameter at breast height in centimetres
#   - Tag  : character, matching the 'Tag' column in the reconstructed data
#   - tag  : dropped after renaming to avoid duplicate columns
#
# Parameters:
#   dryad_path : path to the Dryad RDS file.  Falls back to a sibling-directory
#                path if the default is not found (useful when running from a
#                different working directory).
# Returns: data.table with the Dryad measurements ready for comparison.
load_dryad <- function(dryad_path = here("BCI_stem_reconstruction/DATA/RAW/bci_dryad_condit.rds")) {
    if (!file.exists(dryad_path)) {
        candidate <- here("..", "DATA", "DRYAD_CONDIT", "bci_dryad_condit.rds")
        if (!file.exists(candidate)) {
            stop("Dryad file not found: ", dryad_path, " or ", candidate)
        }
        dryad_path <- candidate
    }
    dryad <- as.data.table(readRDS(dryad_path))
    dryad[, DBH := as.numeric(dbh) / 10]
    dryad[, Tag := as.character(tag)]
    dryad[, tag := NULL]
    dryad
}

# plot_stem_growth() -----------------------------------------------------------
# Draws a base-R scatter/line plot of DBH over census time for all stems
# belonging to one tree tag.  Each unique stem ID is drawn in a distinct colour
# so that crossing or merging trajectories stand out visually.
#
# Parameters:
#   tag         : character, the Tag value to filter on.
#   data        : data.table containing at minimum the columns Tag, CensusID,
#                 the column named by group_by, and the column named by dbh_var.
#   group_by    : column name (character) that identifies individual stems;
#                 defaults to "ReconstructedStemID".
#   dbh_var     : column name (character) for the DBH values to plot;
#                 defaults to "DBH".
#   show_legend : logical; whether to draw the legend (default TRUE).
#   who         : label prefix for the plot title (e.g. "Dryad" or "Reconstructed").
# Returns: invisibly NULL; called for its side-effect (drawing a plot).
plot_stem_growth <- function(tag,
                             data,
                             group_by = "ReconstructedStemID",
                             dbh_var = "DBH",
                             show_legend = TRUE,
                             who = "unknown") {
    sub <- data[Tag == tag][order(get(group_by), CensusID)]
    if (nrow(sub) == 0) {
        return(invisible(NULL))
    }
    plot(sub[[dbh_var]] ~ as.integer(sub$CensusID),
        type = "n",
        xlab = "CensusID",
        ylab = paste(dbh_var, "(cm)"),
        main = paste(who, "-", group_by, "-", "Tag:", tag),
        xlim = c(1, 9)
    )
    unique_groups <- unique(sub[[group_by]])
    cols <- setNames(rainbow(length(unique_groups)), unique_groups)
    sub[,
        {
            col <- cols[as.character(get(group_by))]
            points(get(dbh_var) ~ as.integer(CensusID), col = col, pch = 19)
            lines(get(dbh_var) ~ as.integer(CensusID), col = col, lwd = 1.5)
        },
        by = group_by
    ]
    if (show_legend) {
        legend("bottomright",
            legend = names(cols),
            col = cols,
            pch = 19,
            lwd = 1.5,
            title = group_by,
            bty = "n"
        )
    }
    invisible(NULL)
}

# summarize_trajectories() -----------------------------------------------------
# Computes the cumulative (total) DBH per tag × stem ID for both the
# reconstructed dataset and the Dryad reference, restricting to censuses
# 1 through max_census to keep the comparison within the period covered by
# the Dryad data.
#
# Parameters:
#   indat      : data.table, the reconstructed dataset (from step 7).
#   dryad      : data.table, the Dryad reference data (from load_dryad()).
#   max_census : integer, upper census bound for the comparison (default 8,
#                the last census included in the Condit Dryad release).
# Returns: a named list with two data.tables:
#   $indat_summary : columns Tag, ReconstructedStemID, TotalDBH
#   $dryad_summary : columns Tag, stemID, TotalDBH_Condit
summarize_trajectories <- function(indat, dryad, max_census = 8) {
    indat_compare <- indat[CensusID <= max_census]
    indat_summary <- indat_compare[!is.na(DBH), .(TotalDBH = sum(DBH, na.rm = TRUE)), by = .(Tag, ReconstructedStemID)]
    dryad_summary <- dryad[!is.na(DBH), .(TotalDBH_Condit = sum(DBH, na.rm = TRUE)), by = .(Tag, stemID)]
    list(indat_summary = indat_summary, dryad_summary = dryad_summary)
}

# compare_trajectories() -------------------------------------------------------
# Determines whether the reconstructed stem trajectories for each tag match the
# Dryad reference trajectories.  The comparison strategy:
#   1. Sort both summaries by ascending TotalDBH within each Tag (rank-based
#      comparison, so stems are aligned by size rather than by ID label).
#   2. Assign a within-tag sequence number (seq) to each ranked stem.
#   3. Join on (Tag, seq) and verify that:
#        (a) the number of stems matches between datasets (N_dryad == N_indat), and
#        (b) every aligned TotalDBH pair is exactly equal.
#   4. Use fcase() to consolidate the two checks into a single three-valued
#      flag: TRUE (match), FALSE (mismatch), or NA (tag absent from indat).
#
# Parameters:
#   indat_summary : data.table, as returned by summarize_trajectories()$indat_summary.
#   dryad_summary : data.table, as returned by summarize_trajectories()$dryad_summary.
# Returns: data.table with columns Tag and all_equal_trajectories (logical or NA).
compare_trajectories <- function(indat_summary, dryad_summary) {
    dryad_sorted <- copy(dryad_summary)
    indat_sorted <- copy(indat_summary)
    setorder(dryad_sorted, Tag, TotalDBH_Condit)
    setorder(indat_sorted, Tag, TotalDBH)
    dryad_sorted[, TotalDBH_Condit := as.numeric(as.character(TotalDBH_Condit))]
    indat_sorted[, TotalDBH := as.numeric(as.character(TotalDBH))]
    dryad_sorted[, seq := seq_len(.N), by = Tag]
    indat_sorted[, seq := seq_len(.N), by = Tag]
    dryad_counts <- dryad_sorted[, .(N_dryad = .N), by = Tag]
    indat_counts <- indat_sorted[, .(N_indat = .N), by = Tag]
    counts_merged <- merge(dryad_counts, indat_counts, by = "Tag", all.x = TRUE)
    joined <- indat_sorted[dryad_sorted,
        on = c("Tag", "seq"), nomatch = NA,
        .(Tag = i.Tag, TotalDBH_Condit = i.TotalDBH_Condit, TotalDBH = TotalDBH)
    ]
    value_check <- joined[, .(all_equal_val = !anyNA(TotalDBH) && all(TotalDBH_Condit == TotalDBH)), by = Tag]
    results <- merge(counts_merged, value_check, by = "Tag", all.x = TRUE)
    results[, all_equal_trajectories := fcase(is.na(N_indat), NA, N_dryad != N_indat, FALSE, default = all_equal_val)]
    results[, .(Tag, all_equal_trajectories)][]
}

# render_comparison_pdf() ------------------------------------------------------
# Generates a multi-page PDF with four panels per page, one page per tag:
#   Panel 1 (top-left)  : Dryad DBH-over-census growth curves, coloured by stemID.
#   Panel 2 (top-right) : Reconstructed DBH-over-census growth curves, coloured
#                         by ReconstructedStemID.
#   Panel 3 (bottom-left)  : Data table of Dryad records for that tag.
#   Panel 4 (bottom-right) : Data table of reconstructed records for that tag.
# Tags with no Dryad records (after filtering out P/V status rows) are skipped.
#
# Parameters:
#   tags      : character vector of Tag values to render.
#   file_name : output PDF path (character).
#   indat     : data.table, the reconstructed dataset.
#   dryad     : data.table, the Dryad reference data.
#   n_rows    : reserved parameter (currently unused); retained for future use.
# Returns: invisibly NULL; called for its side-effect (writing the PDF to disk).
render_comparison_pdf <- function(tags, file_name, indat, dryad, n_rows = NULL) {
    if (length(tags) == 0) {
        message("Skipping empty tag list for ", file_name)
        return(invisible(NULL))
    }
    pdf(file_name, width = 16, height = 10)
    for (tag in tags) {
        tab_dryad <- dryad[
            Tag == tag & status != "P" & status != "V" & (!is.na(DBH) | !is.na(codes)),
            .(CensusID, sp, Tag, StemTag, stemID, DBH, codes, status, DFstatus)
        ]
        if (nrow(tab_dryad) == 0) next

        par(mfrow = c(2, 2))
        plot_stem_growth(tag, data = dryad, group_by = "stemID", dbh_var = "DBH", who = paste("Dryad", file_name))
        # indat[, DBH_orig_cm := DBH_mm_original_backup / 10]
        plot_stem_growth(tag, data = indat, group_by = "ReconstructedStemID", dbh_var = "DBH", who = paste("Reconstructed", file_name))

        plot.new()
        addtable2plot(0.01, 0.5, tab_dryad, bty = "o", display.rownames = FALSE, cex = 0.8)
        plot.new()
        tab_recon <- indat[Tag == tag, .(
            CensusID, Mnemonic, Tag, StemTag, StemID, ReconstructedStemID, DBH, # DBH_orig_c
            ListOfTSM, Status, ReconstructionMethod
        )]
        addtable2plot(0.01, 0.5, tab_recon, bty = "o", display.rownames = FALSE, cex = 0.8)
    }
    dev.off()
}

# =============================================================================
# 2. MAIN EXECUTION
# =============================================================================
# Step 1: Load input datasets ---------------------------------------------------

# Reconstructed dataset produced by 2_merge_chunks_to_datatable.R.
# Pipeline stores DBH in mm; convert to cm to match Dryad reference scale.
indat <- as.data.table(readRDS("./BCI_stem_reconstruction/DATA/PROCESSED/complete_dataset_with_reconstructed_stemids.rds"))

indat[, DBH := as.numeric(as.character(DBH)) / 10] # mm → cm
quantile(indat$DBH, probs = seq(0, 1, 0.25), na.rm = TRUE) # confirm cm scale

# Dryad/Condit reference data, loaded and pre-processed by load_dryad().
dryad <- load_dryad()

# Restrict both datasets to censuses 1–8: the period covered by the Dryad release.
# All comparisons (single-stem classification, trajectory matching) are done
# within this window, matching the Quarto report logic exactly.
indat_c8 <- indat[CensusID <= 8]
dryad_c8 <- dryad[CensusID <= 8]

# Classify single vs multi-stem tags using uniqueN(stemID) within census 1–8.
# A tag is single-stem if it has exactly one unique (non-NA) stem ID across
# all its records in the comparison window — same logic as the Quarto report.
id_single_stem_tags_indat_c8 <- indat_c8[
    , .(one_stemid = uniqueN(StemID[!is.na(StemID)]) == 1),
    by = Tag
]
id_single_stem_tags_dryad_c8 <- dryad_c8[
    , .(one_stemid = uniqueN(stemID[!is.na(stemID)]) == 1),
    by = Tag
]

cat("Reconstructed — Single-stem tags:", uniqueN(id_single_stem_tags_indat_c8[one_stemid == TRUE]$Tag), "\n")
cat("Reconstructed — Multi-stem tags :", uniqueN(id_single_stem_tags_indat_c8[one_stemid == FALSE]$Tag), "\n")
cat("Dryad         — Single-stem tags:", uniqueN(id_single_stem_tags_dryad_c8[one_stemid == TRUE]$Tag), "\n")
cat("Dryad         — Multi-stem tags :", uniqueN(id_single_stem_tags_dryad_c8[one_stemid == FALSE]$Tag), "\n")

# Identify tags where single vs multi-stem classification disagrees
dryad_single_stem_tags <- unique(id_single_stem_tags_dryad_c8[one_stemid == TRUE]$Tag)
indat_single_stem_tags <- unique(id_single_stem_tags_indat_c8[one_stemid == TRUE]$Tag)
cat("Single-stem in Dryad but multi-stem in Reconstructed:", length(setdiff(dryad_single_stem_tags, indat_single_stem_tags)), "\n")
cat("Single-stem in Reconstructed but multi-stem in Dryad:", length(setdiff(indat_single_stem_tags, dryad_single_stem_tags)), "\n")

single_in_dryad_not_in_indat <- setdiff(dryad_single_stem_tags, indat_single_stem_tags)
# sample_single_in_dryad_not_in_indat <- sample(single_in_dryad_not_in_indat, min(200, length(single_in_dryad_not_in_indat)))

chk1 <- indat_c8[Tag %in% single_in_dryad_not_in_indat[1], .(CensusID, Tag, StemTag, StemID, DBH, ListOfTSM, Status, ReconstructionMethod)]
chk2 <- dryad_c8[Tag %in% single_in_dryad_not_in_indat[1], .(CensusID, Tag, StemTag, stemID, DBH, codes, status, DFstatus)]

setorder(chk1, Tag, CensusID)
setorder(chk2, Tag, CensusID)

# Step 2: Summarize cumulative DBH per tag × stem ---------------------------
# Totals are computed only for censuses 1–8, the period covered by the Dryad
# release.  Records with NA DBH (dead or missing stems) are excluded from
# the sum so that trajectory shapes are based on measured values only.
sums <- summarize_trajectories(indat, dryad)
indat_summary <- sums$indat_summary
dryad_summary <- sums$dryad_summary
# Sort both summaries by Tag for consistent display in interactive inspection.
setorder(indat_summary, Tag)
setorder(dryad_summary, Tag)

# Step 3: Compare reconstructed trajectories against Dryad reference -----------
# compare_trajectories() returns one row per tag with all_equal_trajectories:
#   TRUE  — stem count and cumulative DBH values match exactly.
#   FALSE — stem counts differ or at least one DBH rank mismatches.
#   NA    — tag is present in Dryad but absent from the reconstructed data
comparison_dt <- compare_trajectories(sums$indat_summary, sums$dryad_summary)

# Step 4: Annotate comparison results with stem morphology --------------------
# Merge single_stem_tags from the reconstructed data so we can examine whether
# reconstruction accuracy differs between single-stem and multi-stem trees.
# We expect single-stem trees to have a 100% match rate since no reconstruction
# was needed — their ReconstructedStemID was copied directly from StemID.
morpho_tag <- unique(indat[, .(Tag, single_stem_tags)])
comparison_dt <- merge(comparison_dt, morpho_tag, by = "Tag", all.x = TRUE)

# Step 5: Summarise match rates by stem morphology ----------------------------
# Recode the boolean flags to readable labels for tabulation.
compare_matches <- comparison_dt[
    , .(
        single_stem = ifelse(single_stem_tags, "Single-stem", "Multi-stem"),
        all_equal_trajectories = ifelse(
            is.na(all_equal_trajectories), NA,
            ifelse(all_equal_trajectories, "Match", "Mismatch")
        )
    )
]

# Number of unique tags in the Dryad repository (reference total)
cat("Total unique tags in Dryad reference:", length(unique(dryad$Tag)), "\n")
# Expected: ~423617 (multi-census BCI tags as of Condit et al. Dryad release)

# Match-rate breakdown: single-stem tags match by construction (no
# reconstruction needed); multi-stem tags split between Match and Mismatch —
# the mismatched subset feeds the diagnostic PDFs below.
round(table(compare_matches$single_stem, compare_matches$all_equal_trajectories) /
    nrow(comparison_dt) * 100, 4)

# Collect the Tags whose reconstructed trajectories do not match Dryad.
differences_tags <- comparison_dt[all_equal_trajectories == FALSE]$Tag
cat("Number of mismatched tags:", length(differences_tags), "\n")

# Step 6: Classify mismatched tags by reconstruction method --------------------
# For each tag, determine which reconstruction method(s) were used across its
# records: only_dp, only_probabilistic, both_probabilistic_and_dp, or NA
# (NA is unexpected for multi-stem tags and warrants investigation).
indat[, method_type := {
    has_prob <- any(grepl("probabilistic", ReconstructionMethod, fixed = TRUE))
    has_dp <- any(grepl("dp", ReconstructionMethod, fixed = TRUE))
    if (has_prob & has_dp) {
        "both_probabilistic_and_dp"
    } else if (has_prob) {
        "only_probabilistic"
    } else if (has_dp) {
        "only_dp"
    } else {
        NA_character_
    }
}, by = Tag]

# Sanity check: print the distinct method_type values actually present.
cat("\nDistinct method_type values assigned:\n")
print(unique(indat$method_type))

# Isolate mismatched tag sets by method, restricting to one type per set so
# that each PDF targets a homogeneous group of reconstruction failures.
differences_tags_dp <- unique(indat[Tag %in% differences_tags & method_type == "only_dp"]$Tag)
differences_tags_prob <- unique(indat[Tag %in% differences_tags & method_type == "only_probabilistic"]$Tag)
differences_tags_dp_prob <- unique(indat[Tag %in% differences_tags & method_type == "both_probabilistic_and_dp"]$Tag)

cat(
    "Mismatched tags by method  —  DP only:", length(differences_tags_dp),
    "| Probabilistic only:", length(differences_tags_prob),
    "| Both:", length(differences_tags_dp_prob), "\n"
)

# Step 7: Draw random samples and generate diagnostic PDFs --------------------
# Cap each sample at 100 tags so the resulting PDFs remain tractable in size
# and review time.  set.seed() is omitted intentionally so that re-runs produce
# different samples, broadening the coverage of manual inspections over time.
sample_differences_tags_prob <- sample(differences_tags_prob, min(300, length(differences_tags_prob)))
sample_differences_tags_dp <- sample(differences_tags_dp, min(300, length(differences_tags_dp)))
sample_differences_tags_dp_prob <- sample(differences_tags_dp_prob, min(300, length(differences_tags_dp_prob)))

## create a data frame in long format with two columns, one the sample difference type and the tags in the other
explore_tags_long <- rbind(
    data.table(method = "sample_differences_tags_dp", Tag = sample_differences_tags_dp),
    data.table(method = "sample_differences_tags_prob", Tag = sample_differences_tags_prob),
    data.table(method = "sample_differences_tags_dp_prob", Tag = sample_differences_tags_dp_prob)
)

# render_comparison_pdf() is called with match() to preserve the original ordering
# from differences_tags_* rather than the random sampling order, making PDFs
# consistently ordered across re-runs for the same sampled tags.
# indat_c8 is used (census 1–8 only) to match what the Quarto report passes
# to its render function.
render_comparison_pdf(
    differences_tags_dp[match(sample_differences_tags_dp, differences_tags_dp)],
    "./BCI_stem_reconstruction/2_STEM_IDENTIFICATION/comparissons/differences_using_dp.pdf", indat_c8, dryad_c8
)
render_comparison_pdf(
    differences_tags_prob[match(sample_differences_tags_prob, differences_tags_prob)],
    "./BCI_stem_reconstruction/2_STEM_IDENTIFICATION/comparissons/differences_using_probabilistic.pdf", indat_c8, dryad_c8
)
render_comparison_pdf(
    differences_tags_dp_prob[match(sample_differences_tags_dp_prob, differences_tags_dp_prob)],
    "./BCI_stem_reconstruction/2_STEM_IDENTIFICATION/comparissons/differences_using_both_methods.pdf", indat_c8, dryad_c8
)

cat("\nDone. Diagnostic PDFs written to: 2_STEM_IDENTIFICATION/\n")

# Step 8: Compare overall tag-level differences per year --------------------
indat_tags_per_census <- indat_c8[!is.na(DBH), .(N_tags = uniqueN(Tag)), by = CensusID]
dryad_tags_per_census <- dryad_c8[!is.na(DBH), .(N_tags = uniqueN(Tag)), by = CensusID]

indat_tags_per_census[, CensusID := as.character(CensusID)]
dryad_tags_per_census[, CensusID := as.character(CensusID)]

# merge the two summaries and compute the difference in tag counts per census
tags_comparison <- merge(
    indat_tags_per_census, dryad_tags_per_census,
    by = "CensusID", all = TRUE, suffixes = c("_indat", "_dryad")
)
tags_comparison[, Tag_count_difference := N_tags_indat - N_tags_dryad]

tags_comparison
