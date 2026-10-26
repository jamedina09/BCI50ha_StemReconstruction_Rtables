# =============================================================================
# build-data-package.R — copies the data files of the Ecology Data Paper into
# repository/data/ from the stage-3 outputs, unchanged, and checks them:
#   - the census tables and the species table (R data files, .Rdata) from
#     DATA/RTABLES, and the posterior samples of the stem reconstruction (.rds)
#     from DATA/POSTERIORS;
#   - every copy is identical to its source (MD5);
#   - file names follow Ecology Data Papers 2.0 (letters, digits, "-", "_" and
#     "." only; at most 35 characters);
#   - each file holds the documented object, columns and column types;
#   - the keys between files hold (stems, species, posterior samples).
# It also writes checks/data-file-inventory.csv (size, MD5, rows, columns of
# every file) and checks/field-codes-to-define.csv (the field codes found in
# the data, for the data manager to define).
#
# Run from the project root:
#   Rscript data_paper_and_repo_publication/ecology_data_paper/build/build-data-package.R
# Needs data.table and digest.
# =============================================================================

suppressMessages({
    library(data.table)
    library(digest)
})

PKG <- "data_paper_and_repo_publication/ecology_data_paper"
SRC <- "BCI_stem_reconstruction/DATA"
DATA_OUT <- file.path(PKG, "repository", "data")
CHECK_DIR <- file.path(PKG, "checks")
LOG <- file.path(CHECK_DIR, "data-checks.txt")

if (!dir.exists(file.path(SRC, "RTABLES"))) stop("Run from the project root (", SRC, "/RTABLES not found)")
dir.create(DATA_OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(CHECK_DIR, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("build-data-package.R checks, %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), file = LOG)

# A failed check prints ❌, raises a warning, writes to the log and stops.
check <- function(ok, msg, examples = NULL) {
    ok <- isTRUE(ok)
    line <- sprintf("%s %s %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), if (ok) "✓" else "❌", msg)
    if (!ok && length(examples)) line <- paste0(line, " | e.g. ", paste(head(examples, 5), collapse = ", "))
    cat(line, "\n")
    cat(line, "\n", file = LOG, append = TRUE)
    if (!ok) {
        cat("❌ CHECK FAILED:", msg, "\n")
        warning(msg, call. = FALSE)
        stop("CHECK FAILED: ", msg, call. = FALSE)
    }
    invisible(TRUE)
}

# ---- Sources ------------------------------------------------------------------
# The outputs of stage 3, written after the stage-2 table, with the posterior
# samples of the engine run of that table.
sources <- c(
    file.path(SRC, "RTABLES", c(paste0("bci.stem", 1:9, ".Rdata"), "bci.spptable.rdata")),
    file.path(SRC, "POSTERIORS", "posterior_sampled_paths.rds")
)
files <- basename(sources)
stage2_table <- file.path(SRC, "PROCESSED", "complete_dataset_final_with_reconstructed_stemids.rds")
check(all(file.exists(sources)), "The stage-3 outputs exist (run stage 3 first)", files[!file.exists(sources)])
check(file.exists(stage2_table) && all(file.mtime(sources) > file.mtime(stage2_table)),
      "The stage-3 outputs are newer than the stage-2 table (rerun stage 3 after the merge)")

# ---- File names --------------------------------------------------------------------
check(all(grepl("^[A-Za-z0-9._-]+$", files)), "File names use only letters, digits, '-', '_' and '.'", files[!grepl("^[A-Za-z0-9._-]+$", files)])
check(all(nchar(files) <= 35), "File names have at most 35 characters", files[nchar(files) > 35])
check(!anyDuplicated(tolower(files)), "File names are unique (ignoring case)")

# ---- Copy, unchanged ------------------------------------------------------------------
old <- list.files(DATA_OUT, pattern = "[.](Rdata|rdata|rds)$", full.names = TRUE)
unlink(old[!basename(old) %in% files])
ok_copy <- file.copy(sources, file.path(DATA_OUT, files), overwrite = TRUE, copy.date = TRUE)
check(all(ok_copy), "Every data file was copied", files[!ok_copy])
md5_src <- vapply(sources, digest, character(1), file = TRUE, algo = "md5", USE.NAMES = FALSE)
md5_out <- vapply(file.path(DATA_OUT, files), digest, character(1), file = TRUE, algo = "md5", USE.NAMES = FALSE)
check(identical(md5_src, md5_out), "Every data file is identical to its stage-3 source (MD5)", files[md5_src != md5_out])

# Google Drive and Finder add hidden files ("Icon\r", .DS_Store) to synced
# folders; they are not part of the package and must not be uploaded.
present <- list.files(DATA_OUT, all.files = TRUE, no.. = TRUE)
system_files <- grepl("^Icon\r?$", present) | present %in% ".DS_Store"
if (any(system_files)) {
    cat("System files in repository/data (leave them out of the upload):",
        paste(gsub("\r", "\\\\r", present[system_files]), collapse = ", "), "\n", file = LOG, append = TRUE)
}
check(setequal(present[!system_files], files), "The data folder holds only the files of the package (system files aside)",
      setdiff(present[!system_files], files))

# ---- Objects, columns and types ---------------------------------------------------------
load_one <- function(path) {
    e <- new.env(parent = emptyenv())
    nm <- load(path, envir = e)
    list(name = nm, value = get(nm[1], envir = e))
}
STEM_TYPES <- c(CensusID = "numeric", treeID = "character", stemID = "integer", tag = "character",
                StemTag = "character", sp = "character", quadrat = "character", gx = "numeric",
                gy = "numeric", dbh = "numeric", hom = "numeric", ExactDate = "Date", date = "numeric",
                codes = "character", DFstatus = "character", Rstatus = "character", StemPaths = "integer")
types_of <- function(df) vapply(df, function(v) class(v)[1], character(1))

census <- vector("list", 9)
for (k in 1:9) {
    o <- load_one(file.path(DATA_OUT, sprintf("bci.stem%d.Rdata", k)))
    check(identical(o$name, sprintf("bci.stem%d", k)) && is.data.frame(o$value),
          sprintf("bci.stem%d.Rdata holds one data.frame named bci.stem%d", k, k), o$name)
    check(identical(types_of(o$value), STEM_TYPES), sprintf("bci.stem%d has the 17 documented columns, in order, with their types", k),
          names(STEM_TYPES)[types_of(o$value)[names(STEM_TYPES)] != STEM_TYPES])
    check(all(o$value$CensusID == k), sprintf("bci.stem%d holds census %d only", k, k))
    census[[k]] <- as.data.table(o$value)
}
o <- load_one(file.path(DATA_OUT, "bci.spptable.rdata"))
check(identical(o$name, "bci.spptable") && is.data.frame(o$value) && all(vapply(o$value, is.character, logical(1))),
      "bci.spptable.rdata holds one data.frame named bci.spptable, all columns text")
SPP_COLS <- c("sp", "Order", "Family", "Genus", "SpeciesName", "Latin", "InfraspecificRank", "InfraspecificEpithet",
              "Authority", "Synonyms", "Lifeform_RFoster", "Lifeform_RPerez_SAguilar", "CommonName", "Herbarium",
              "NameCondit2019", "NameChange", "Notes")
check(identical(names(o$value), SPP_COLS), "bci.spptable has the 17 documented columns in order")
spp <- as.data.table(o$value)
post <- readRDS(file.path(DATA_OUT, "posterior_sampled_paths.rds"))
check(is.data.frame(post) && identical(names(post), c("tag", "treeID", "path_sig", "path_count", "path_prob", "recon")) &&
          identical(unname(types_of(post)), c("character", "character", "character", "integer", "numeric", "character")),
      "posterior_sampled_paths.rds holds a data frame with the 6 documented columns and types")
run_code <- sub("^run_code: ", "", grep("^run_code: ", readLines(file.path(SRC, "PROCESSED", "engine_run.txt")), value = TRUE))
check(identical(attr(post, "engine_run"), run_code),
      "The posterior samples come from the engine run of the stage-2 table (attribute engine_run = engine_run.txt)",
      c(attr(post, "engine_run"), run_code))
cat("Engine run:", run_code, "\n", file = LOG, append = TRUE)
post <- as.data.table(post)

# ---- Keys between files -------------------------------------------------------------
key1 <- census[[1]][, .(treeID, stemID)]
check(!anyDuplicated(key1), "treeID + stemID identify one row of each census table")
for (k in 2:9) check(identical(census[[k]][, .(treeID, stemID)], key1), sprintf("bci.stem%d lists the same stems in the same order as bci.stem1", k))
stem_info <- rbindlist(lapply(census, function(d) d[, .(treeID, tag, sp)]))
check(stem_info[, uniqueN(tag), by = treeID][, all(V1 == 1L)], "Each tree has one tag in every table")
check(stem_info[, uniqueN(sp), by = treeID][, all(V1 == 1L)], "Each tree has one species code in every table")
check(all(unique(stem_info$sp) %in% spp$sp), "Every species code of the census tables is in the species table", setdiff(unique(stem_info$sp), spp$sp))
check(all(spp$sp %in% unique(stem_info$sp)), "Every code of the species table occurs in the census tables", setdiff(spp$sp, unique(stem_info$sp)))
check(!anyDuplicated(spp$sp), "Species codes are unique in the species table")
tree_sp <- unique(stem_info[, .(treeID, tag, sp)])
rm(stem_info)

check(all(post$treeID %in% tree_sp$treeID), "Every tree of the posterior samples is in the census tables")
check(post[, sum(path_count), by = treeID][, all(V1 == 200L)], "The path counts of each tree sum to its 200 samples")
pairs <- post[, .(pair = unlist(strsplit(recon, ";", fixed = TRUE))), by = .(treeID, path = seq_len(nrow(post)))]
pairs[, c("StemPaths", "label") := tstrsplit(pair, ":", fixed = TRUE)]
pairs[, StemPaths := as.integer(StemPaths)]
n_sig <- lengths(strsplit(post$path_sig, "-", fixed = TRUE))
n_rec <- pairs[, .N, by = path][order(path), N]
check(identical(as.integer(n_sig), n_rec), "path_sig has one label per record listed in recon")
meas <- unique(rbindlist(lapply(census, function(d) d[!is.na(StemPaths) & !is.na(dbh), .(treeID, StemPaths)])))
miss <- unique(pairs[, .(treeID, StemPaths)])[!meas, on = .(treeID, StemPaths)]
check(nrow(miss) == 0L, "Every record of the posterior samples is a measured record of its tree in the census tables (treeID + StemPaths)", miss$treeID)
keyp <- unique(rbindlist(lapply(census, function(d) d[!is.na(StemPaths), .(treeID, StemPaths, CensusID)])))
check(!anyDuplicated(keyp[, .(treeID, StemPaths)]), "StemPaths is unique within its tree across the nine tables")
n_pairs <- nrow(pairs)
rm(pairs, meas, miss, keyp)

# ---- Inventory ------------------------------------------------------------------------
dims <- c(lapply(census, dim), list(dim(spp), dim(post)))
inv <- data.table(
    file = files,
    object = c(paste0("bci.stem", 1:9), "bci.spptable", "(readRDS)"),
    bytes = file.size(file.path(DATA_OUT, files)),
    md5 = md5_out,
    rows = vapply(dims, `[`, numeric(1), 1),
    columns = vapply(dims, `[`, numeric(1), 2)
)
inv[, size_MB := round(bytes / 1e6, 2)]
fwrite(inv, file.path(CHECK_DIR, "data-file-inventory.csv"))
print(inv)
cat(sprintf("Total: %.1f MB in %d files; posterior samples hold %s record-label pairs\n",
            sum(inv$bytes) / 1e6, nrow(inv), format(n_pairs, big.mark = ",")), file = LOG, append = TRUE)

# ---- Field codes for the data manager ----------------------------------------------------
# Every single code found in `codes` (split at ";"), with its records by census
# and by field status. The definitions are to be supplied by the data manager.
codes <- rbindlist(lapply(census, function(d) {
    cd <- d[!is.na(codes), .(CensusID, DFstatus, codes)]
    cd[, .(code = unlist(strsplit(codes, ";", fixed = TRUE))), by = .(row = seq_len(nrow(cd)), CensusID, DFstatus)][, row := NULL]
}))
by_census <- dcast(codes[, .N, by = .(code, CensusID)], code ~ paste0("census_", CensusID), value.var = "N", fill = 0L)
by_status <- codes[, .(records = .N, DFstatus = paste(sort(unique(DFstatus)), collapse = "; ")), by = code]
fc <- merge(by_status, by_census, by = "code")[order(-records)]
fc[, definition := NA_character_]
fwrite(fc, file.path(CHECK_DIR, "field-codes-to-define.csv"), na = "")
cat(sprintf("Field codes: %d distinct codes written to checks/field-codes-to-define.csv\n", nrow(fc)), file = LOG, append = TRUE)

check(TRUE, sprintf("Data package written: %d files, %.1f MB", nrow(inv), sum(inv$bytes) / 1e6))
