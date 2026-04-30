############################################################
### main_cpp_bci.R — BCI debug driver: one tag
###
### Loads bci_multistem_xrun_debug.rds and runs the full DP
### pipeline for a single tag.
###
### Run from project root:
###   Rscript dp_global/scripts/main_cpp_bci.R
###   Rscript dp_global/scripts/main_cpp_bci.R --WHICH_TAG=12345
###   Rscript dp_global/scripts/main_cpp_bci.R --WHICH_TAG=187064 --DP_VERBOSE=FALSE
############################################################

if (sys.nframe() == 0L) rm(list = ls())

library(data.table)
library(here)

# ----------------------------------------------------------------
# 1. Source main_cpp.R
#    This defines all helper functions (ensure_dir, log_msg,
#    ensure_species_column, attach_bio_columns, run_dp_one_group,
#    maybe_add_posterior_bins, etc.) and sources the infrastructure
#    R modules.  run_main() is NOT called because sys.nframe() != 0
#    at this point.
#    parse_args() also runs here, capturing any CLI overrides.
# ----------------------------------------------------------------
source(here("dp_global", "scripts", "main_cpp.R"))

# ----------------------------------------------------------------
# 2. BCI-specific defaults
#    These override the simulated-data defaults set by main_cpp.R.
#    CLI args (already captured in `overrides`) are re-applied
#    below so they can still override these BCI defaults.
# ----------------------------------------------------------------
# INPUT_FILE                          <- here("bci_data", "bci_multistem_xrun_debug.rds")
# INPUT_FILE                          <- here("bci_data", "wrong_tags.rds")
INPUT_FILE <- here("bci_data", "multistem_tags.rds")
WHICH_TAG <- "115203"
FORCE_ONE_SPECIES_PARAMETERS <- FALSE # use real species from BCI data
MAX_GROWTH_FIXED <- 5.0
MAX_SHRINK_FIXED <- -0.5
RECRUIT_MAX_FIXED <- (MAX_GROWTH_FIXED * 5) + 0.9999
DP_MAX_STATES <- 1039L
DP_SLACK_TRACKS <- 1L
DP_SLACK_REQUIRE_ANCHOR_RECRUITABLE <- TRUE
DP_SLACK_REQUIRE_ANCHOR_EPS <- 1e-6
ANCHOR_START_CENSUS <- 7L
DP_MODE <- "marginals+bins"
DP_VERBOSE <- TRUE
DP_POSTERIOR_TOP_K <- 2L
DP_MAX_TRACKS <- NULL # auto-detect
DP_FALLBACK_GROWTH_FORMS <- character(0)
PRUNE_BOUND_FACTOR <- 5
NON_TAPER_CORRECTED_GROWTH_FORMS <- c("palm", "strangler_fig", "tree_fern")
NON_TAPER_CORRECTED_PRUNE_MIN_GROWTH <- PRUNE_BOUND_FACTOR * MAX_SHRINK_FIXED
NON_TAPER_CORRECTED_PRUNE_MAX_GROWTH <- PRUNE_BOUND_FACTOR * MAX_GROWTH_FIXED
HOM_TOLERANCE_SCALE <- 2.0
USE_MEASUREMENT_ERROR <- FALSE
ALLOW_PROVISIONAL_DP_ANCHOR <- TRUE
POSTERIOR_SAMPLES <- 0L # disable posteriors; set >0 to enable
POSTERIOR_SAMPLES_FORMAT <- "csv"
POSTERIOR_SAMPLES_PATH <- NULL
POSTERIOR_SAMPLE_SEED <- NULL
RUN_ALL_TAGS <- FALSE
WRITE_DP_CSV <- TRUE
WRITE_DP_RDS <- TRUE
WRITE_DP_PDF <- TRUE
DP_PDF_INCLUDE_REFERENCE <- TRUE

# Bio parameter estimation sources (same as main_cpp_chunk.R)
MAX_GROWTH_HARD_SOURCE <- "fixed"
MAX_SHRINK_HARD_SOURCE <- "fixed"
K_SHRINK_SOURCE <- "fixed"
K_SHRINK_FIXED <- 0
K_GROWTH_SOURCE <- "fixed"
K_GROWTH_FIXED <- 0
RECRUIT_MAX_SOURCE <- "fixed"

# Species that should bypass DP and go directly to the probabilistic greedy
# matcher.  Empty vector means all species go through DP normally.
PROB_SPECIES <- character(0)
PROB_LOOKAHEAD_WEIGHT <- 0.5
PIN_TRUESTEMID <- TRUE

# ----------------------------------------------------------------
# 3. Re-apply CLI overrides on top of BCI defaults
#    (overrides were parsed by main_cpp.R's parse_args();
#     re-applying them here lets users still override BCI defaults)
# ----------------------------------------------------------------
for (.nm in names(overrides)) {
    if (.nm %in% c("help", "h")) next
    .norm <- normalize_key(.nm)
    .match <- find_matching_var(.norm)
    if (!is.null(.match)) {
        assign(.match, overrides[[.nm]], envir = globalenv())
        message("[main_cpp_bci.R] CLI override: ", .match, " = ", as.character(overrides[[.nm]]))
    }
}

# Recompute derived values after overrides
# WHICH_TAG         <- as.integer(WHICH_TAG)
RECRUIT_MAX_FIXED <- (MAX_GROWTH_FIXED * 5) + 0.9999
if (!is.null(overrides$RECRUIT_MAX_FIXED)) RECRUIT_MAX_FIXED <- as.numeric(overrides$RECRUIT_MAX_FIXED)

ADD_DP_POSTERIOR_BINS <- DP_MODE == "marginals+bins"
PLOT_PDF_ONE_TAG_ONLY <- TRUE

# ----------------------------------------------------------------
# 4. Compute output directory
# ----------------------------------------------------------------
CONFIG_NAME <- paste0("BCI_tag", WHICH_TAG)
out_dir <- file.path(base_out_dir, build_out_dir_name())
message("[main_cpp_bci.R] out_dir: ", out_dir)

# ----------------------------------------------------------------
# 5. Load BCI data (RDS)
# ----------------------------------------------------------------
message("[main_cpp_bci.R] Loading: ", INPUT_FILE)
if (!file.exists(INPUT_FILE)) stop("INPUT_FILE not found: ", INPUT_FILE)
xraw <- readRDS(INPUT_FILE)
setDT(xraw)
message("[main_cpp_bci.R] Loaded ", nrow(xraw), " rows, ", uniqueN(xraw$Tag), " tags.")

# add rownumber
xraw[, rownum := .I]

# ---- TrueStemID reconstruction ----
#
# TrueStemID is the ground-truth stem identity used as a hard anchor by the DP.
# When TrueStemID is non-NA for a row, the DP MUST assign that exact stem ID and
# cannot deviate. When it is NA, the DP resolves it freely using biological costs.
#
# -----------------------------------------------------------------------
# STEP 1 — Certain identity anchors (no ambiguity possible)
# -----------------------------------------------------------------------
#
#   (a) StemTag: any row where the field crew physically tagged the stem.
#       The OriginalStemID is unambiguous at any census, including pre-C7.
#
#   (b) CensusID >= 7 (year 2010 onward): from C7 the BCI database assigned
#       OriginalStemIDs via a systematic re-tagging campaign. Every stem
#       present at C7+ has a trustworthy, reliable database ID.
#
#   All other rows are left as NA — the DP resolves them.
#
# -----------------------------------------------------------------------
# STEP 2 — Conservative terminal propagation (post-last-DBH zone only)
# -----------------------------------------------------------------------
#
#   Within each (Tag, OriginalStemID) group, once a stem has made its last
#   live measurement (last non-NA DBH), all subsequent rows are in the
#   terminal phase: they can only record death, resprout, or missing status.
#   In that zone the OriginalStemID is unambiguous — the database does not
#   reassign IDs for simple death / carry-forward records.
#
#   2a. Identify the boundary: the last census with a non-NA DBH per
#       (Tag, OriginalStemID). Rows strictly after this are the terminal
#       phase. Stems that never recorded a DBH get NA and are excluded
#       from all propagation by the guards in 2b and 2c.
#
#   2b. DIRECT ANCHOR terminal-event rows to their own OriginalStemID.
#       Any post-last-DBH row carrying an explicit death, broken-below, or
#       R-family resprout code is safe to anchor. No prior Step-1 anchor is
#       required — a death/resprout record for a given OriginalStemID is
#       unambiguously about that biological individual.
#       Handles:
#         • Pure pre-C7 stems with no StemTag (Case 1): e.g. last DBH at C1,
#           Status="dead" at C2 → anchored here; 2c fills any later gaps.
#         • Spans-C7 stems with gaps before the C7 anchor (Case 2): e.g.
#           Status="dead" at C5 anchored here; C6 gap filled by 2c.
#
#   2c. BIDIRECTIONAL FILL of remaining post-last-DBH NA gaps.
#       After 2b, rows with no explicit terminal status (e.g. "missing"
#       carry-forward rows, or gaps between a dead row and a later C7+
#       anchor) may still be NA. LOCF carries anchors forward; NOCB carries
#       a later C7+ anchor backward. The CensusID > last_dbh_census filter
#       strictly limits the operation to the terminal phase: pre-last-DBH
#       rows are never modified.
#
#   Pre-last-DBH rows are deliberately left as NA. The DP must be free to
#   resolve ambiguous early-census identity assignments.

# ----- Step 1: certain anchors -----
xraw[, TrueStemID := NA_integer_]
xraw[!is.na(StemTag),                    TrueStemID := OriginalStemID]  # (a) physical tag
xraw[is.na(TrueStemID) & CensusID >= 7L, TrueStemID := OriginalStemID]  # (b) C7+ re-tagging

xraw[Tag == "004808", .(Tag, CensusID, OriginalStemID, TrueStemID, DBH, StemTag, ListOfTSM, Status)]
xraw[Tag == "006160", .(Tag, CensusID, OriginalStemID, TrueStemID, DBH, StemTag, ListOfTSM, Status)]
xraw[Tag == "264355", .(Tag, CensusID, OriginalStemID, TrueStemID, DBH, StemTag, ListOfTSM, Status)]

# ----- Step 2: terminal propagation -----

# Word-boundary regex for R-family resprout codes in ListOfTSM.
# \\b = word boundary (perl=TRUE). Must match the pattern in dp_global_dp.R.
resprout_regex <- "\\b(R|RP|RF|RT|QR|OR)\\b"

# 2a. Last census with a non-NA DBH per (Tag, OriginalStemID).
#     CensusID > .last_dbh_census defines the terminal phase for that stem.
#     setorder ensures rows are sorted before the LOCF/NOCB fill in 2c.
setorder(xraw, Tag, OriginalStemID, CensusID)
.last_dbh <- xraw[
    !is.na(DBH),
    .(last_dbh_census = max(CensusID)),
    by = .(Tag, OriginalStemID)
]
xraw[.last_dbh, on = .(Tag, OriginalStemID), .last_dbh_census := i.last_dbh_census]

xraw[Tag == "004808", .(Tag, CensusID, OriginalStemID, TrueStemID, DBH, StemTag, ListOfTSM, Status)]
xraw[Tag == "006160", .(Tag, CensusID, OriginalStemID, TrueStemID, DBH, StemTag, ListOfTSM, Status)]
xraw[Tag == "264355", .(Tag, CensusID, OriginalStemID, TrueStemID, DBH, StemTag, ListOfTSM, Status)]

# 2b. Anchor terminal-event rows directly to their own OriginalStemID.
#     All five conditions must hold:
#       • is.na(TrueStemID)            — not already anchored by Step 1
#       • !is.na(.last_dbh_census)     — stem has a DBH history (was ever measured)
#       • CensusID > .last_dbh_census  — strictly in the terminal phase
#       • terminal event present:
#           – R-family resprout code in ListOfTSM, OR
#           – Status is "dead", "stem dead", or "broken below"
xraw[
    is.na(TrueStemID) &
        !is.na(.last_dbh_census) &
        CensusID > .last_dbh_census &
        (
            (!is.na(ListOfTSM) & grepl(resprout_regex, ListOfTSM, perl = TRUE)) |
            (!is.na(Status)    & Status %in% c("broken below", "dead", "stem dead"))
        ),
    TrueStemID := OriginalStemID
]

xraw[Tag == "004808", .(Tag, CensusID, OriginalStemID, TrueStemID, DBH, StemTag, ListOfTSM, Status)]
xraw[Tag == "006160", .(Tag, CensusID, OriginalStemID, TrueStemID, DBH, StemTag, ListOfTSM, Status)]
xraw[Tag == "264355", .(Tag, CensusID, OriginalStemID, TrueStemID, DBH, StemTag, ListOfTSM, Status)]

# 2c. Bidirectional fill of remaining post-last-DBH NA gaps within each group.
#     LOCF: carry a 2b anchor (or Step-1 anchor) forward to later terminal rows.
#     NOCB: carry a later C7+ anchor backward to earlier terminal rows that had
#           no explicit terminal status code.
#     The row filter (CensusID > .last_dbh_census) guarantees pre-last-DBH rows
#     are never touched; nafill operates only on the terminal-phase subset.
xraw[
    !is.na(.last_dbh_census) & CensusID > .last_dbh_census,
    TrueStemID := nafill(nafill(TrueStemID, type = "locf"), type = "nocb"),
    by = .(Tag, OriginalStemID)
]

# 2d. Remove temporary column.
xraw[, .last_dbh_census := NULL]

setorder(xraw, rownum)

# xraw[Tag == "119453", .(Tag, CensusID, OriginalStemID, TrueStemID, DBH, StemTag, ListOfTSM, Status)]
# xraw[Tag == "004808", .(Tag, CensusID, OriginalStemID, TrueStemID, DBH, StemTag, ListOfTSM, Status)]
# xraw[Tag == "006160", .(Tag, CensusID, OriginalStemID, TrueStemID, DBH, StemTag, ListOfTSM, Status)]
# xraw[Tag == "264355", .(Tag, CensusID, OriginalStemID, TrueStemID, DBH, StemTag, ListOfTSM, Status)]
# xraw[Tag == "246746", .(Tag, CensusID, OriginalStemID, TrueStemID, DBH, StemTag, ListOfTSM, Status)]
# xraw[Tag == "567547", .(Tag, CensusID, OriginalStemID, TrueStemID, DBH, StemTag, ListOfTSM, Status)]

# xraw[Tag %in% unique(xraw[is.na(TrueStemID) & CensusID >= 7L, .(Tag, CensusID, OriginalStemID, TrueStemID, DBH, StemTag, ListOfTSM, Status)]$Tag)]

# ----------------------------------------------------------------
# 6. Ensure species column
# ----------------------------------------------------------------
xraw <- ensure_species_column(xraw)
xrun <- copy(xraw)

# ----------------------------------------------------------------
# 7. Validate tag
# ----------------------------------------------------------------
# if (!(WHICH_TAG %in% xrun$Tag)) {
#     stop("Tag ", WHICH_TAG, " not found in data. ",
#          "First few tags: ", paste(head(sort(unique(xrun$Tag))), collapse = ", "), " ...")
# }
tag_sp <- unique(xrun[Tag == WHICH_TAG, species])
message("[main_cpp_bci.R] Tag ", WHICH_TAG, " — species: ", paste(tag_sp, collapse = ", "))
message("[main_cpp_bci.R] Tag ", WHICH_TAG, " — rows: ", nrow(xrun[Tag == WHICH_TAG]))

# ----------------------------------------------------------------
# 8. Attach bio columns and prepare for DP
#    Subset to the tag's species only — bio_pars only covers that
#    species, so attach_bio_columns would fail on other species.
# ----------------------------------------------------------------
xrun_tag <- xrun[species %in% tag_sp]

dp_max_tracks_local <- if (is.null(DP_MAX_TRACKS)) {
    auto_dp_max_tracks(xrun_tag)
} else {
    as.integer(DP_MAX_TRACKS)
}
message("[main_cpp_bci.R] dp_max_tracks: ", dp_max_tracks_local)

# ----------------------------------------------------------------
# 9. Run DP for WHICH_TAG
# ----------------------------------------------------------------
ensure_dir(out_dir)
writeLines(as.character(Sys.time()), file.path(out_dir, "run_started.txt"))
message("[main_cpp_bci.R] Running DP for Tag ", WHICH_TAG, "...")

dtg <- xrun_tag[Tag == WHICH_TAG]
out <- run_dp_one_group(dtg, dp_max_tracks = dp_max_tracks_local)

out <- maybe_add_posterior_bins(out)
if (!is.null(out)) {
    out[, run_out_dir := basename(out_dir)]
    message("[main_cpp_bci.R] DP done. ", nrow(out), " rows returned.")
} else {
    message("[main_cpp_bci.R] DP returned NULL (tag may have been skipped).")
}

# ----------------------------------------------------------------
# 10. Write outputs
# ----------------------------------------------------------------
DP_BASE <- "stem_reconstruction_dp_global_rcpp"
make_path <- function(ext) file.path(out_dir, paste0(DP_BASE, ".", ext))

if (isTRUE(WRITE_DP_CSV) && !is.null(out)) {
    p <- make_path("csv")
    fwrite(out, p)
    message("[main_cpp_bci.R] CSV: ", p)
}
if (isTRUE(WRITE_DP_RDS) && !is.null(out)) {
    p <- make_path("rds")
    saveRDS(out, p)
    message("[main_cpp_bci.R] RDS: ", p)
}
if (isTRUE(WRITE_DP_PDF) && !is.null(out) && nrow(out) > 0L) {
    p <- make_path("pdf")
    plot_tag_to_pdf(
        out,
        pdf_file          = p,
        include_reference = DP_PDF_INCLUDE_REFERENCE,
        tag               = WHICH_TAG
    )
    message("[main_cpp_bci.R] PDF: ", p)
}

# ----------------------------------------------------------------
# 11. Finished
# ----------------------------------------------------------------
writeLines(as.character(Sys.time()), file.path(out_dir, "run_finished.txt"))
message("[main_cpp_bci.R] Done. Output dir: ", out_dir)

# Rscript dp_global/scripts/main_cpp_bci.R --WHICH_TAG=119453
# Rscript dp_global/scripts/main_cpp_bci.R --WHICH_TAG=115427
# Rscript dp_global/scripts/main_cpp_bci.R --WHICH_TAG=123375
# Rscript dp_global/scripts/main_cpp_bci.R --WHICH_TAG=115203
# Rscript dp_global/scripts/main_cpp_bci.R --WHICH_TAG=242799
# Rscript dp_global/scripts/main_cpp_bci.R --WHICH_TAG=246746
# Rscript dp_global/scripts/main_cpp_bci.R --WHICH_TAG=277120
# Rscript dp_global/scripts/main_cpp_bci.R --WHICH_TAG=190932
# Rscript dp_global/scripts/main_cpp_bci.R --WHICH_TAG=171486
# Rscript dp_global/scripts/main_cpp_bci.R --WHICH_TAG=220311
# Rscript dp_global/scripts/main_cpp_bci.R --WHICH_TAG=204785
# Rscript dp_global/scripts/main_cpp_bci.R --WHICH_TAG=242114
# Rscript dp_global/scripts/main_cpp_bci.R --WHICH_TAG=001080
# Rscript dp_global/scripts/main_cpp_bci.R --WHICH_TAG=005558