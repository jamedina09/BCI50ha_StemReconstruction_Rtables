#' Estimate DBH at 1.3 m using a taper model
#'
#' This function computes the expected DBH at 1.3 m (DBHc) from a diameter
#' measured at some other height using published taper parameterizations.
#'
#' Equation obtained from:
#' Cushman, K. C., S. Bunyavejchewin, D. Cárdenas, et al. “Variation in Trunk
#' Taper of Buttressed Trees within and among Five Lowland Tropical Forests.”
#' Biotropica 53, no. 5 (2021): 1442–53. Scopus.
#' https://doi.org/10.1111/btp.12994.
#'
#' @param dbh_mm numeric vector of diameters at measurement height in millimetres
#' @param hom numeric vector of heights of measurement in metres; NA values will be set to common_hom
#' @param wsg numeric vector of wood specific gravity in g cm^-3, or NULL to use the model without WSG
#' @param common_hom numeric minimum height (metres) used for calculation (default 1.3)
#' @return numeric vector of estimated DBH at 1.3 m (millimetres)
#' @export
taper <- function(dbh_mm, hom, wsg = NULL, common_hom = 1.3) {
    # Defensive checks
    if (length(dbh_mm) != length(hom) || (!is.null(wsg) && length(wsg) != length(dbh_mm))) {
        stop("'dbh_mm', 'hom' and 'wsg' (if provided) must have the same length")
    }

    # Copy inputs so we don't modify caller's vectors
    dbh_mm <- as.numeric(dbh_mm)
    hom <- as.numeric(hom)
    if (!is.null(wsg)) wsg <- as.numeric(wsg)

    # Replace NA heights with common_hom (do not modify valid measured heights)
    hom_na <- is.na(hom)
    hom[hom_na] <- common_hom

    # Convert dbh from mm to cm for the model
    dbh_cm <- dbh_mm / 10

    # Compute taper parameter b depending on whether wsg is provided
    # Protect against log(0) or negative inputs by coercing non-positive values to NA
    dbh_cm[dbh_cm <= 0] <- NA_real_
    hom_for_log <- hom
    hom_for_log[hom_for_log <= 0] <- NA_real_
    if (!is.null(wsg)) {
        wsg[wsg <= 0] <- NA_real_
        b <- 0.151 - 0.025 * log(dbh_cm) - 0.02 * log(hom_for_log) - 0.021 * log(wsg)
    } else {
        b <- 0.156 - 0.023 * log(dbh_cm) - 0.021 * log(hom_for_log)
    }

    # Where b or inputs are NA, propagate NA
    invalid <- is.na(b) | is.na(dbh_cm)

    # Estimate DBHc in cm: dbh_cm / exp(-b * (hom - common_hom))
    denom <- exp(-b * (hom - common_hom))
    dbhc_cm <- dbh_cm / denom

    # Convert back to mm and set invalid values to NA
    dbhc_mm <- dbhc_cm * 10
    dbhc_mm[invalid] <- NA_real_

    return(dbhc_mm)
}

#' Apply taper correction to a data.frame/data.table
#'
#' This function applies `taper()` to a dataset with customizable column names.
#'
#' @param df data.frame or data.table containing the measurements
#' @param dbh_col name of column with measured diameter (default "dbh")
#' @param hom_col name of column with height of measurement (default "hom")
#' @param wsg_col name of column with wood specific gravity (default NULL). If NULL the model without WSG is used.
#' @param output_col name of output column to store corrected DBH (default "dbhc")
#' @param taper_correction logical, whether to apply taper correction (default TRUE)
#' @param common_hom numeric value used when hom is NA (default 1.3 m)
#' @param convert_units logical, try to auto-convert hom units when suspicious (default TRUE)
#' @param verbose logical, show messages (default TRUE)
#' @param overwrite logical, if TRUE will overwrite an existing column with name output_col (default TRUE)
#' @return data.frame or data.table of same class as `df` with new column `output_col`
#' @export
apply_taper_correction <- function(
  df,
  dbh_col = "dbh",
  hom_col = "hom",
  wsg_col = NULL,
  output_col = "dbhc",
  taper_correction = TRUE,
  common_hom = 1.3,
  convert_units = TRUE,
  verbose = TRUE,
  overwrite = TRUE
) {
    # Basic input validation
    if (!is.data.frame(df)) stop("df must be a data.frame or data.table")
    cols <- names(df)
    if (!(dbh_col %in% cols)) stop(sprintf("Diameter column '%s' not found in df", dbh_col))
    if (!(hom_col %in% cols)) stop(sprintf("Height column '%s' not found in df", hom_col))
    if (!is.null(wsg_col) && !(wsg_col %in% cols)) stop(sprintf("WSG column '%s' not found in df", wsg_col))

    # Copy to data.table for fast in-place operations (preserve class for return)
    is_dt <- inherits(df, "data.table")
    DT <- data.table::as.data.table(df)

    # Warn about overwriting
    if (output_col %in% names(DT) && !overwrite) stop(sprintf("Column '%s' already exists. Set overwrite = TRUE to replace it.", output_col))

    # Coerce numeric columns and check for nonsense values
    DT[, (dbh_col) := as.numeric(get(dbh_col))]
    DT[, (hom_col) := as.numeric(get(hom_col))]
    if (!is.null(wsg_col)) DT[, (wsg_col) := as.numeric(get(wsg_col))]

    # Minimal checks for negative or zero dbh
    if (any(DT[[dbh_col]] <= 0, na.rm = TRUE)) {
        warning(sprintf("Some values in '%s' are <= 0; these will return NA after correction.", dbh_col))
    }

    # Prepare hom column for calculation
    DT[, hom_m := get(hom_col)]

    if (convert_units) {
        # Heuristic: if many hom > 25 assume hom is in cm -> convert to m
        if (sum(DT$hom_m > 25, na.rm = TRUE) > 0) {
            if (verbose) message("Converting hom values > 25 to meters (assumed cm -> m)")
            DT[hom_m > 25 & !is.na(hom_m), hom_m := hom_m / 100]
        }
        # Heuristic: if hom between 2 and 25 but larger than dbh_cm + 5 assume decimeters -> convert
        dbh_cm <- DT[[dbh_col]] / 10
        idx_dm <- which(!is.na(DT$hom_m) & DT$hom_m < 25 & DT$hom_m > dbh_cm + 5)
        if (length(idx_dm) > 0) {
            if (verbose) message(sprintf("Converting %d hom values from decimeters to meters (dividing by 10)", length(idx_dm)))
            DT[idx_dm, hom_m := hom_m / 10]
        }
    }

    # Replace NA hom values with common_hom (do not modify measured heights)
    DT[is.na(hom_m), hom_m := common_hom]

    # Initialize output column with original dbh
    DT[, (output_col) := get(dbh_col)]

    if (taper_correction) {
        # Compute corrected DBH using taper model; use wsg if available and requested
        if (!is.null(wsg_col)) {
            DT[, (output_col) := taper(get(dbh_col), hom_m, get(wsg_col), common_hom = common_hom)]
        } else {
            DT[, (output_col) := taper(get(dbh_col), hom_m, NULL, common_hom = common_hom)]
        }
        if (verbose) message("Taper correction applied and stored in column: ", output_col)
    } else {
        if (verbose) message("Taper correction not requested; original dbh copied to ", output_col)
    }

    # Clean temporary column
    DT[, hom_m := NULL]

    # Return object in the same class as input
    if (is_dt) {
        return(DT)
    } else {
        return(as.data.frame(DT))
    }
}

# End of taper_correction.R
