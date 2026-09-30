# =============================================================================
# 1_prepare_posteriors_BCI.R
#
# Purpose: Consolidate posterior Feather outputs from a BCI DP run into one
#          RDS file with an explicit Tag column.
#
# Posterior columns (one row per unique reconstruction path of a tag):
#   path_sig   : signature of the path (sequence of ReconstructedStemIDs)
#   path_count : number of the DP's posterior samples that produced this path.
#                The samples are drawn by backward sampling from the DP
#                posterior, so path_count / sum(path_count) is the posterior
#                probability estimate to use for Monte Carlo sampling.
#   path_prob  : kept for reference only. The engine re-weights each sample by
#                exp(logp) before summing, so path_prob is proportional to
#                count x p (roughly p^2) and must NOT be used for sampling.
#   recon      : "ObsRowID:ReconstructedStemID;..." identity of every
#                observation in the path
# =============================================================================

# =============================================================================
# SETUP
# =============================================================================

rm(list = ls())

# Load required packages (will error if not available; install before running)
library(arrow) # read_feather()
library(data.table) # fast data manipulation, rbindlist()

# ── Hard check that stays visible in interactive (line-by-line) runs ────────
# On failure: prints a ❌ line (count + examples), raises an immediate
# warning, and only then stops. On success it prints a ✓ line.
bio_check <- function(ok, msg, examples = NULL, n_bad = NULL) {
    if (isTRUE(all(ok))) {
        cat("✓", msg, "\n")
        return(invisible(TRUE))
    }
    n_txt <- if (!is.null(n_bad)) sprintf(" [%d case(s)]", n_bad) else ""
    ex_txt <- if (length(examples) > 0L) {
        paste0(" | examples: ", paste(head(unique(examples), 10), collapse = ", "))
    } else {
        ""
    }
    full_msg <- paste0("CHECK FAILED: ", msg, n_txt, ex_txt)
    cat("❌", full_msg, "\n")
    warning(full_msg, call. = FALSE, immediate. = TRUE)
    stop(full_msg, call. = FALSE)
}

# =============================================================================
# CONFIGURATION
# =============================================================================

workspace_root <- getwd()
home_dir <- "/Users/medinaja/outputs_bci_stem_identification"

# Choose the run folder inside the output directory and locate its posterior
# files. Only directories are candidates (the folder may also hold a .zip
# archive of the same run).
run_code <- basename(list.dirs(home_dir, recursive = FALSE, full.names = TRUE))
bio_check(
    length(run_code) == 1L,
    sprintf("Exactly one DP run folder in %s", home_dir),
    examples = run_code,
    n_bad = length(run_code)
)
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

# tag_from_filename() ------------------------------------------------------
# Posterior files are named "tag_<Tag>_posterior_samples_<ts>_paths.feather".
# The tag is everything between "tag_" and "_posterior_samples", whatever its
# length.
tag_from_filename <- function(file_paths) {
    sub("^tag_(.+?)_posterior_samples.*$", "\\1", basename(file_paths))
}

# read_and_bind_feathers() ------------------------------------------------
# Read posterior Feather files, extract Tag from the filename, and combine
# all files into a single keyed data.table.
read_and_bind_feathers <- function(file_paths) {
    if (!is.character(file_paths) || length(file_paths) == 0) {
        stop("file_paths must be a non-empty character vector")
    }
    tags <- tag_from_filename(file_paths)
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

# One posterior file per tag: a tag parsed twice would merge two trees.
file_tags <- tag_from_filename(post_files)
bio_check(
    !anyDuplicated(file_tags) && !any(file_tags == basename(post_files)),
    "Every posterior file yields a distinct, correctly parsed tag",
    examples = file_tags[duplicated(file_tags) | file_tags == basename(post_files)],
    n_bad = sum(duplicated(file_tags) | file_tags == basename(post_files))
)
cat("Tag length distribution:\n")
print(table(nchar(file_tags)))

# Read and bind all posterior files into a single data.table
dt_posteriors <- read_and_bind_feathers(post_files)

# link treeid to tag
tag_treeid_map <- as.data.table(readRDS(file.path(
    workspace_root, "BCI_stem_reconstruction", "DATA", "PROCESSED", "ViewFullTable_single_vs_multiple_stem_tags.rds"
)))
tag_treeid_map <- unique(tag_treeid_map[, .(TreeID, Tag)])
tag_treeid_map[, Tag := as.character(Tag)]
tag_treeid_map[, TreeID := as.character(TreeID)]

# Every posterior tag must belong to a known tree; otherwise its paths would be
# silently dropped below.
weird_tags <- unique(setdiff(unique(dt_posteriors$Tag), tag_treeid_map$Tag))
bio_check(
    length(weird_tags) == 0L,
    "Every posterior tag exists in the Tag <-> TreeID map",
    examples = weird_tags,
    n_bad = length(weird_tags)
)

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

# One tree, one posterior: all paths of a tree must describe the SAME set of
# observations (they differ only in how the observations are linked). Two
# observation sets in one tree mean paths from two different trees were mixed.
obs_sets <- dt_posteriors[, .(obs_set = vapply(
    strsplit(recon, ";", fixed = TRUE),
    function(p) paste(sort(as.integer(sub(":.*", "", p))), collapse = ","),
    character(1)
)), by = .(treeID, path_sig)][, .(n_sets = uniqueN(obs_set)), by = treeID]
bio_check(
    obs_sets[, all(n_sets == 1L)],
    "All paths of each tree cover the same set of observations",
    examples = obs_sets[n_sets > 1L, treeID],
    n_bad = obs_sets[n_sets > 1L, .N]
)
rm(obs_sets)

# Posterior samples per tree (the DP draws 200 per tag; fewer means some
# samples failed). Reported, not enforced.
cat("Posterior samples per tree (sum of path_count):\n")
print(dt_posteriors[, .(n_samples = sum(path_count)), by = treeID][, .N, by = n_samples][order(-N)])

inspectdf::inspect_na(dt_posteriors)

cat(
    "Consolidated table dimensions:", nrow(dt_posteriors), "rows ×",
    ncol(dt_posteriors), "columns |", uniqueN(dt_posteriors$treeID), "trees\n"
)

# =============================================================================
# 3. OUTPUT
# =============================================================================

output_file <- file.path(output_dir_posteriors, "posterior_sampled_paths.rds")
saveRDS(dt_posteriors, file = output_file)

cat("Posterior samples saved to:", output_file, "\n")
cat("Done.\n")
