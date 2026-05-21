suppressMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
b <- fread(args[1])
n <- fread(args[2])
cat("rows base/new:", nrow(b), "/", nrow(n), "\n")
key <- c("Tag", "OriginalStemID", "CensusID")
m <- merge(
  b[, .(Tag, OriginalStemID, CensusID, RecB = ReconstructedStemID, MetB = ReconstructionMethod)],
  n[, .(Tag, OriginalStemID, CensusID, RecN = ReconstructedStemID, MetN = ReconstructionMethod)],
  by = key, all = TRUE
)
cat("merged rows:", nrow(m), "\n")
ne <- function(a, b) {
  (is.na(a) != is.na(b)) | (!is.na(a) & !is.na(b) & a != b)
}
rd <- ne(m$RecB, m$RecN)
md <- ne(m$MetB, m$MetN)
cat("Recon diffs:", sum(rd), "\n")
cat("Method diffs:", sum(md), "\n")
cat("Rows missing in baseline:", sum(is.na(m$RecB) & !is.na(m$RecN)), "\n")
cat("Rows missing in new:", sum(!is.na(m$RecB) & is.na(m$RecN)), "\n")
cat("\nMethod table (baseline):\n"); print(b[, .N, by = ReconstructionMethod][order(-N)])
cat("\nMethod table (new):\n"); print(n[, .N, by = ReconstructionMethod][order(-N)])
if (sum(md) > 0) {
  cat("\nMethod transitions (baseline -> new):\n")
  print(m[md, .N, by = .(MetB, MetN)][order(-N)])
}
if (sum(rd) > 0) {
  cat("\nFirst 10 Recon diffs:\n")
  print(head(m[rd], 10))
}
