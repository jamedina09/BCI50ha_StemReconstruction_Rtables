############################################################
# dp_global_diag.R
# Diagnostics, posterior bins and plotting helpers
############################################################

add_constraint_violation <- function(x, id_col = "ReconstructedStemID", min_growth, max_growth, pair_interval) {
    # PURPOSE
    # - Post-hoc diagnostic: flag potentially implausible links along each reconstructed
    #   track when the *implied* per-year growth between adjacent censuses falls outside
    #   [min_growth, max_growth].
    #
    # INPUTS
    # - x: data.table with at least id_col, DBH, CensusID.
    # - id_col: which ID column defines a track (defaults to ReconstructedStemID).
    # - min_growth/max_growth: allowable annual growth bounds (cm/year).
    # - pair_interval: years between consecutive censuses (assumed constant here).
    #
    # OUTPUT
    # - `x` with/updated `ConstraintViolation` logical column (TRUE for flagged rows).
    #
    # NOTES
    # - Only evaluates consecutive censuses (CensusID increases by 1).
    # - Flags the earlier observation in the violating pair (so the “bad link” is easy
    #   to see in time series plots).
    if (!("ConstraintViolation" %in% names(x))) {
        x[, ConstraintViolation := NA]
    }
    if (!all(c(id_col, "DBH", "CensusID") %in% names(x))) {
        return(x)
    }

    data.table::setorder(x, CensusID)
    ids <- unique(x[[id_col]])
    ids <- ids[!is.na(ids)]
    if (length(ids) == 0L) {
        return(x)
    }
    for (sid in ids) {
        # sid <- ids[1] # for testing
        ii <- which(x[[id_col]] == sid & !is.na(x$DBH))
        if (length(ii) < 2L) next
        ii <- ii[order(x$CensusID[ii])]
        for (k in seq_len(length(ii) - 1L)) {
            # k <- 1L  # for testing
            i0 <- ii[k]
            i1 <- ii[k + 1L]
            if (x$CensusID[i1] != x$CensusID[i0] + 1L) next
            pair_T <- (as.numeric(pair_interval[[as.character(x$CensusID[i1])]]) - as.numeric(pair_interval[[as.character(x$CensusID[i0])]])) / 365.25
            # pair_T <- get_pair_interval(x$CensusID[i0], x$CensusID[i1])
            g <- (x$DBH[i1] - x$DBH[i0]) / pair_T
            cond <- isTRUE((g < min_growth) | (g > max_growth))
            if (cond || isTRUE(x$ConstraintViolation[i0])) {
                x$ConstraintViolation[i0] <- TRUE
            }
        }
    }
    x
}

add_dp_posterior_bins <- function(
  x,
  confident_prob = 0.95,
  unlinked_prob = 0.50,
  use_reconstructed_prob = TRUE,
  out_col = "DP_PosteriorBin"
) {
    # PURPOSE
    # - Convenience label per observation based on the DP marginal posterior.
    #
    # BINS
    # - "unlinked-likely": posterior probability of being unlinked >= unlinked_prob
    # - "confident": posterior probability of the chosen reconstructed ID >= confident_prob
    #   (or Top1Prob if use_reconstructed_prob=FALSE)
    # - "ambiguous": everything else (posterior spread across multiple IDs)
    #
    # NOTES
    # - Rows with missing posterior columns (e.g., anchor / given IDs) get NA.
    # - This is a summary; it does not change `ReconstructedStemID`.
    if (!data.table::is.data.table(x)) {
        x <- data.table::as.data.table(x)
    }

    if (!is.character(out_col) || length(out_col) != 1L || !nzchar(out_col)) {
        stop("out_col must be a single non-empty column name")
    }

    # Minimal required posterior columns
    if (!("DP_PosteriorUnlinkedProb" %in% names(x))) {
        x[, (out_col) := NA_character_]
        x[, (out_col) := factor(get(out_col), levels = c("confident", "ambiguous", "unlinked-likely"))]
        return(x)
    }

    score_col <- if (isTRUE(use_reconstructed_prob) && ("DP_PosteriorReconstructedProb" %in% names(x))) {
        "DP_PosteriorReconstructedProb"
    } else if ("DP_PosteriorTop1Prob" %in% names(x)) {
        "DP_PosteriorTop1Prob"
    } else {
        NA_character_
    }

    if (!is.character(score_col) || is.na(score_col)) {
        x[, (out_col) := NA_character_]
        x[, (out_col) := factor(get(out_col), levels = c("confident", "ambiguous", "unlinked-likely"))]
        return(x)
    }

    # Classify
    x[, (out_col) := {
        p_unlinked <- DP_PosteriorUnlinkedProb
        p_score <- get(score_col)

        # Only label when we actually have posterior numbers
        has_post <- is.finite(p_unlinked) | is.finite(p_score)

        out <- rep(NA_character_, .N)
        out[has_post] <- "ambiguous"
        out[has_post & is.finite(p_unlinked) & (p_unlinked >= unlinked_prob)] <- "unlinked-likely"
        out[has_post & is.finite(p_score) & (p_score >= confident_prob) &
            !(is.finite(p_unlinked) & (p_unlinked >= unlinked_prob))] <- "confident"
        out
    }]

    x[, (out_col) := factor(get(out_col), levels = c("confident", "ambiguous", "unlinked-likely"))]
    x
}

# ----------------------------
# Plotting: one PDF page per Tag
# ----------------------------
# This helper is NOT part of the matching algorithm.
# It is only a visualization tool to help you see whether reconstructed IDs
# create smooth DBH trajectories across censuses.
#
# For each Tag, we optionally create two side-by-side plots:
# 1) DBH trajectories grouped by `ReferenceStemID` (only if present and `include_reference=TRUE`)
# 2) DBH trajectories grouped by `ReconstructedStemID`
#
# Each line is a "stem ID group" and points are DBH measurements by census.
# If DP posterior bins are present (DP_PosteriorBin), point shapes show
# confident/ambiguous/unlinked-likely, and constraint violations are overlaid
# as an "X" symbol.

plot_tag_to_pdf <- function(out, pdf_file, include_reference = FALSE, tag = NULL) {
    ## ---- Package checks (no attach unless needed) ----
    # We use requireNamespace() so we don't attach these packages globally.
    pkgs <- c("ggplot2", "cowplot")
    invisible(lapply(pkgs, function(p) {
        if (!requireNamespace(p, quietly = TRUE)) {
            stop("Package not installed: ", p, call. = FALSE)
        }
    }))

    ## ---- Required columns ----
    # Plot 2 (reconstructed IDs) always requires these.
    required_cols <- c(
        "Tag", "CensusID", "DBH", "species",
        "ReconstructedStemID", "ReconstructionMethod", "ConstraintViolation"
    )
    # Plot 1 (reference/original IDs) is optional.
    if (isTRUE(include_reference)) {
        required_cols <- c(required_cols, "OriginalStemID")
    }

    # Find which required columns are absent.
    missing_cols <- setdiff(required_cols, names(out))
    if (length(missing_cols) > 0) {
        stop("Missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
    }

    # Get out_dir for subtitle if available
    out_dir_value <- if ("out_dir" %in% names(out)) unique(out$out_dir) else NULL

    ## ---- Optional filter: one tag (or subset of tags) ----
    if (!is.null(tag)) {
        tag <- sort(unique(as.integer(tag)))
        out <- out[Tag %in% tag]
        if (nrow(out) == 0L) {
            stop("No rows found for requested tag(s): ", paste(tag, collapse = ", "), call. = FALSE)
        }
    }

    ## ---- Output filename ----
    # If `tag` is specified, we automatically write to a different filename than
    # the "all tags" PDF so you don't overwrite the full report.
    resolve_pdf_file <- function(path, tag) {
        if (!is.character(path) || length(path) != 1L || !nzchar(path)) {
            stop("pdf_file must be a single non-empty path", call. = FALSE)
        }

        if (dir.exists(path)) {
            base <- if (is.null(tag)) {
                "stem_reconstruction_all_tags.pdf"
            } else if (length(tag) == 1L) {
                paste0("stem_reconstruction_tag_", tag, ".pdf")
            } else {
                paste0("stem_reconstruction_tags_", paste(tag, collapse = "_"), ".pdf")
            }
            return(file.path(path, base))
        }

        ext <- tolower(tools::file_ext(path))
        if (ext != "pdf") {
            stop("pdf_file must be a .pdf file or an existing directory", call. = FALSE)
        }
        dirn <- dirname(path)
        base0 <- tools::file_path_sans_ext(basename(path))

        if (is.null(tag)) {
            return(path)
        }
        if (length(tag) == 1L) {
            return(file.path(dirn, paste0(base0, "_tag_", tag, ".pdf")))
        }
        file.path(dirn, paste0(base0, "_tags_", paste(tag, collapse = "_"), ".pdf"))
    }

    pdf_file_final <- resolve_pdf_file(pdf_file, tag)

    ## ---- Helper: common theme & scales ----
    # A common theme makes the two panels comparable and consistent.
    common_theme <- ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
            plot.margin = ggplot2::margin(10, 10, 10, 10),
            # legend to bottom and horizontal
            # legend.position = "bottom"#,
            # legend.direction = "horizontal"
        )

    # Shape scales
    # - If DP posterior bins exist, shapes show the bin.
    # - Otherwise, shapes show constraint violations.
    bin_shape_scale <- ggplot2::scale_shape_manual(
        values = c(
            `confident` = 16,
            `ambiguous` = 17,
            `unlinked-likely` = 1,
            `NA` = 4
        ),
        drop = FALSE
    )

    violation_shape_scale <- ggplot2::scale_shape_manual(
        values = c(`FALSE` = 16, `TRUE` = 17),
        labels = c(`FALSE` = "ok", `TRUE` = "violation")
    )

    ## ---- Open PDF ----
    # onefile=TRUE means: one PDF file with multiple pages.
    grDevices::pdf(pdf_file_final, width = 13, height = 8, onefile = TRUE)

    ## ---- Loop over Tag ----
    # Each iteration creates ONE page in the PDF.
    keys <- data.table::data.table(Tag = sort(unique(out$Tag)))

    for (ii in seq_len(nrow(keys))) {
        # ii <- 1
        tag <- keys$Tag[[ii]]
        species_code <- unique(out[Tag == tag]$species)
        # Filter to the current Tag (and species, if applicable).
        tag_data <- data.table::copy(out[Tag == tag])

        tag_data[, ConstraintViolationFlag := (!is.na(ConstraintViolation) & as.logical(ConstraintViolation))]

        has_bins <- ("DP_PosteriorBin" %in% names(tag_data)) && any(!is.na(tag_data$DP_PosteriorBin))

        p1 <- NULL
        if (isTRUE(include_reference)) {
            # Plot 1: Original/reference IDs (if present in `out`)
            if (isTRUE(has_bins)) {
                p1 <- ggplot2::ggplot(
                    tag_data,
                    ggplot2::aes(
                        x = as.factor(CensusID),
                        y = DBH,
                        group = OriginalStemID,
                        color = factor(OriginalStemID),
                        shape = DP_PosteriorBin
                    )
                ) +
                    ggplot2::geom_line(na.rm = TRUE) +
                    ggplot2::geom_point(size = 3, na.rm = TRUE) +
                    ggplot2::geom_point(
                        data = tag_data[ConstraintViolationFlag == TRUE],
                        ggplot2::aes(x = as.factor(CensusID), y = DBH, group = OriginalStemID),
                        inherit.aes = FALSE,
                        shape = 4,
                        color = "black",
                        size = 3.5,
                        stroke = 1,
                        na.rm = TRUE
                    ) +
                    bin_shape_scale +
                    ggplot2::labs(
                        title = paste("Original Stem IDs - Tag", tag, "(species:", species_code, ")"),
                        subtitle = out_dir_value,
                        color = "OriginalStemID",
                        shape = "Posterior bin"
                    ) +
                    common_theme
            } else {
                p1 <- ggplot2::ggplot(
                    tag_data,
                    ggplot2::aes(
                        x = as.factor(CensusID),
                        y = DBH,
                        group = OriginalStemID,
                        color = factor(OriginalStemID),
                        shape = ConstraintViolationFlag
                    )
                ) +
                    ggplot2::geom_line(na.rm = TRUE) +
                    ggplot2::geom_point(size = 3, na.rm = TRUE) +
                    violation_shape_scale +
                    ggplot2::labs(
                        title = paste("Original Stem IDs - Tag", tag, "(species:", species_code, ")"),
                        subtitle = out_dir_value,
                        color = "OriginalStemID",
                        shape = "ConstraintViolation"
                    ) +
                    common_theme
            }
        }

        ## --- Plot 2: ReconstructedStemID ---
        # Group by ReconstructedStemID: this shows your algorithm's result.
        # Ideally, DBH trajectories should look smooth/biologically plausible.
        if (isTRUE(has_bins)) {
            p2 <- ggplot2::ggplot(
                tag_data,
                ggplot2::aes(
                    x = as.factor(CensusID),
                    y = DBH,
                    group = ReconstructedStemID,
                    color = factor(ReconstructedStemID),
                    shape = DP_PosteriorBin
                )
            ) +
                ggplot2::geom_line(na.rm = TRUE) +
                ggplot2::geom_point(size = 3, na.rm = TRUE) +
                ggplot2::geom_point(
                    data = tag_data[ConstraintViolationFlag == TRUE],
                    ggplot2::aes(x = as.factor(CensusID), y = DBH, group = ReconstructedStemID),
                    inherit.aes = FALSE,
                    shape = 4,
                    color = "black",
                    size = 3.5,
                    stroke = 1,
                    na.rm = TRUE
                ) +
                bin_shape_scale +
                ggplot2::labs(
                    title = paste("Reconstructed Stem IDs - Tag", tag, "(species:", species_code, ")"),
                    subtitle = out_dir_value,
                    color = "ReconstructedStemID",
                    shape = "Posterior bin"
                ) +
                common_theme
        } else {
            p2 <- ggplot2::ggplot(
                tag_data,
                ggplot2::aes(
                    x = as.factor(CensusID),
                    y = DBH,
                    group = ReconstructedStemID,
                    color = factor(ReconstructedStemID),
                    shape = ConstraintViolationFlag
                )
            ) +
                ggplot2::geom_line(na.rm = TRUE) +
                ggplot2::geom_point(size = 3, na.rm = TRUE) +
                violation_shape_scale +
                ggplot2::labs(
                    title = paste("Reconstructed Stem IDs - Tag", tag, "(species:", species_code, ")"),
                    subtitle = out_dir_value,
                    color = "ReconstructedStemID",
                    shape = "ConstraintViolation"
                ) +
                common_theme
        }

        # Combine panels
        if (isTRUE(include_reference)) {
            combined_plot <- cowplot::plot_grid(
                p1, p2,
                ncol = 2,
                align = "h",
                axis = "tb",
                rel_widths = c(1, 1)
            )
            print(combined_plot)
        } else {
            print(p2)
        }
    }

    # Always close the PDF device.
    grDevices::dev.off()
    invisible(pdf_file_final)
}
