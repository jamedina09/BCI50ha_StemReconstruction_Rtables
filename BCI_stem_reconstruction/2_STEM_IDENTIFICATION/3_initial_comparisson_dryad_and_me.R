## 3_initial_comparisson_dryad_and_me.R
## =============================================================================
## Purpose: Compare reconstructed BCI stem trajectories against the Dryad/Condit
##          reference. Identify tags whose reconstructed DBH trajectories
##          match or mismatch the reference and generate diagnostic PDFs.
##
## Sections:
##   0. Setup
##   1. Helper functions
##   2. Load datasets
##   3. Summarize tag × stem DBH trajectories
##   4. Compare reconstructed and reference trajectories
##   5. Classify mismatches by reconstruction method
##   6. Generate diagnostic plots
##   7. Compare tag counts per census
##
## Inputs:
##   - BCI_stem_reconstruction/DATA/PROCESSED/complete_dataset_final_with_reconstructed_stemids.rds
##   - BCI_stem_reconstruction/DATA/RAW/bci_dryad_condit.rds
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
# Load the Dryad/Condit reference file and convert DBH from mm to cm.
# If the default path is missing, fall back to a sibling directory location.
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
# Plot DBH over census time for a single tag, colouring by stem identifier.
# Supports either Dryad stem IDs or reconstructed stem IDs.
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
# Sum DBH per tag × stem for the reconstructed and Dryad datasets.
# Only censuses 1–8 are included, matching the Dryad release coverage.
summarize_trajectories <- function(indat, dryad, max_census = 8) {
    indat_compare <- indat[CensusID <= max_census]
    indat_summary <- indat_compare[!is.na(DBH), .(TotalDBH = sum(DBH, na.rm = TRUE)), by = .(Tag, ReconstructedStemID)]
    dryad_summary <- dryad[!is.na(DBH), .(TotalDBH_Condit = sum(DBH, na.rm = TRUE)), by = .(Tag, stemID)]
    list(indat_summary = indat_summary, dryad_summary = dryad_summary)
}

# compare_trajectories() -------------------------------------------------------
# Compare reconstructed and Dryad trajectories by ranking stems within each tag
# and checking whether the total DBH values align.
# Result: one row per tag with a match/mismatch/NA flag.
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
# Create a PDF showing Dryad vs reconstructed DBH trajectories and summary tables
# for a set of tags.
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
# Step 1: Load datasets -------------------------------------------------------

# Load reconstructed output and convert DBH from mm to cm.
indat <- as.data.table(readRDS("./BCI_stem_reconstruction/DATA/PROCESSED/complete_dataset_final_with_reconstructed_stemids.rds"))
indat[, DBH := as.numeric(DBH) / 10]

# Load the Dryad reference data and prepare it for comparison.
dryad <- load_dryad()

# Quick unit sanity check for a sample tag.
indat[Tag == "-05599", .(Tag, DBH)]
dryad[Tag == "-05599", .(Tag, DBH)]

# Restrict both datasets to censuses 1–8, matching the Dryad release window.
indat_c8 <- indat[CensusID <= 8]
dryad_c8 <- dryad[CensusID <= 8]

# Determine single- vs multi-stem tags in each dataset using unique StemIDs.
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
# Compute total DBH per stem within the comparison window, excluding missing DBH.
sums <- summarize_trajectories(indat, dryad)
indat_summary <- sums$indat_summary
dryad_summary <- sums$dryad_summary
# Sort both summaries by Tag for consistent display in interactive inspection.
setorder(indat_summary, Tag)
setorder(dryad_summary, Tag)

# Step 3: Compare reconstructed trajectories against Dryad reference -----------
# Classify each tag as Match, Mismatch, or NA based on stem count and total DBH.
comparison_dt <- compare_trajectories(sums$indat_summary, sums$dryad_summary)

# Step 4: Annotate comparison results with stem morphology --------------------
# Attach single- vs multi-stem labels from the reconstructed data.
morpho_tag <- unique(indat[, .(Tag, single_stem_tags)])
comparison_dt <- merge(comparison_dt, morpho_tag, by = "Tag", all.x = TRUE)

# Step 5: Summarise match rates by stem morphology ----------------------------
# Convert comparison flags into readable match/mismatch labels.
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

#               Match Mismatch
#   Multi-stem  18.95    11.73
#   Single-stem 69.32     0.00

# Collect the Tags whose reconstructed trajectories do not match Dryad.
differences_tags <- comparison_dt[all_equal_trajectories == FALSE]$Tag
cat("Number of mismatched tags:", length(differences_tags), "\n")

# Step 6: Classify mismatched tags by reconstruction method --------------------
# Label each tag by whether it used DP, probabilistic, or both methods.
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
# Sample mismatched tags by method and generate diagnostic reports.
set.seed(123)
sample_differences_tags_prob <- sample(differences_tags_prob, min(300, length(differences_tags_prob)))
set.seed(123)
sample_differences_tags_dp <- sample(differences_tags_dp, min(300, length(differences_tags_dp)))
set.seed(123)
sample_differences_tags_dp_prob <- sample(differences_tags_dp_prob, min(300, length(differences_tags_dp_prob)))

## create a data frame in long format with two columns, one the sample difference type and the tags in the other
explore_tags_long <- rbind(
    data.table(method = "sample_differences_tags_dp", Tag = sample_differences_tags_dp),
    data.table(method = "sample_differences_tags_prob", Tag = sample_differences_tags_prob),
    data.table(method = "sample_differences_tags_dp_prob", Tag = sample_differences_tags_dp_prob)
)

# Preserve the original tag ordering when rendering the sampled PDFs.
# Use the census 1–8 subset to match the Quarto report logic.
comparissons_path <- "./BCI_stem_reconstruction/2_STEM_IDENTIFICATION/comparissons"
render_comparison_pdf(
    differences_tags_dp[match(sample_differences_tags_dp, differences_tags_dp)],
    file.path(comparissons_path, "differences_using_dp.pdf"), indat_c8, dryad_c8
)
render_comparison_pdf(
    differences_tags_prob[match(sample_differences_tags_prob, differences_tags_prob)],
    file.path(comparissons_path, "differences_using_probabilistic.pdf"), indat_c8, dryad_c8
)
render_comparison_pdf(
    differences_tags_dp_prob[match(sample_differences_tags_dp_prob, differences_tags_dp_prob)],
    file.path(comparissons_path, "differences_using_both_methods.pdf"), indat_c8, dryad_c8
)

cat("\nDone. Diagnostic PDFs written to: 2_STEM_IDENTIFICATION/\n")

# Step 8: Compare overall tag counts by census -------------------------------
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
