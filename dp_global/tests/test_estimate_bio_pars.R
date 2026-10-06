# Tests for estimate_bio_pars() (dp_global/R/dp_global_bio.R): the growth-SD
# fit when the SD falls with size, and the recruitment rate units.
# Run from the project root:
#   testthat::test_file("dp_global/tests/test_estimate_bio_pars.R")
suppressMessages({
  library(testthat)
  library(data.table)
})
root <- c(".", "..", "../..")
root <- normalizePath(root[file.exists(file.path(root, "dp_global", "R", "dp_global_bio.R"))][1])
source(file.path(root, "dp_global", "R", "dp_global_bio.R"))

census_date <- function(c) as.IDate("2010-01-01") + round((c - 7L) * 5 * 365.25)
fit <- function(x, ...) suppressWarnings(estimate_bio_pars(x, anchor_start_census = 7L, use_measurement_error = FALSE,
  max_shrink_source = "fixed", max_shrink_fixed = -0.5, max_growth_source = "fixed", max_growth_fixed = 5,
  k_shrink_source = "fixed", k_shrink_fixed = 0, k_growth_source = "fixed", k_growth_fixed = 0,
  recruit_max_source = "fixed", recruit_max_fixed = 26, ...))
# stems measured at censuses 7, 8, 9 with growth noise sd_fun(DBH)
stems <- function(n, sd_fun, seed) {
  set.seed(seed)
  rbindlist(lapply(seq_len(n), function(i) {
    d <- runif(1, 2, 14); out <- list()
    for (c in 7:9) { out[[length(out) + 1]] <- data.table(Tag = sprintf("T%04d", i), TrueStemID = i, CensusID = c, DBH = d, ExactDate = census_date(c), species = "x")
      d <- max(1, d + 5 * rnorm(1, 0.05, sd_fun(d))) }
    rbindlist(out)
  }))
}
# the estimator's own SD proxy, recomputed outside it (pairs within the stems)
proxy <- function(x) {
  w <- dcast(x, Tag + TrueStemID ~ CensusID, value.var = "DBH")
  dt <- dcast(x, Tag + TrueStemID ~ CensusID, value.var = "ExactDate")
  g <- c((w[["8"]] - w[["7"]]) / (as.numeric(dt[["8"]] - dt[["7"]]) / 365.25), (w[["9"]] - w[["8"]]) / (as.numeric(dt[["9"]] - dt[["8"]]) / 365.25))
  d0 <- c(w[["7"]], w[["8"]])
  r <- residuals(lm(g ~ log(d0)))
  list(sd = abs(r) * sqrt(pi / 2), d0 = d0)
}

test_that("SD falling with size: refitted as a constant (the mean of the SD proxies), not the intercept", {
  x <- stems(400, function(d) 0.3 / d, seed = 1)
  b <- fit(x)
  p <- proxy(x)
  raw <- coef(lm(p$sd ~ p$d0))
  expect_lt(raw[2], 0)                                    # the data have a negative slope
  expect_equal(b$growth$sigma1, 0)
  expect_equal(b$growth$sigma0, mean(p$sd), tolerance = 1e-8)
  expect_lt(b$growth$sigma0, unname(raw[1]))               # smaller than the intercept the old code kept
})

test_that("SD rising with size: the linear fit is kept unchanged", {
  x <- stems(400, function(d) 0.02 + 0.01 * d, seed = 2)
  b <- fit(x)
  p <- proxy(x)
  raw <- coef(lm(p$sd ~ p$d0))
  expect_gt(raw[2], 0)
  expect_equal(b$growth$sigma0, max(unname(raw[1]), 0.01), tolerance = 1e-8)
  expect_equal(b$growth$sigma1, unname(raw[2]), tolerance = 1e-8)
})

test_that("recruitment rate: new stems per established tree per year, and the legacy per-slot rate", {
  x <- stems(10, function(d) 0.05, seed = 3)               # 10 established trees, one stem each, censuses 7-9
  new_in_tree <- rbindlist(lapply(1:4, function(i) data.table(Tag = sprintf("T%04d", i), TrueStemID = 100L + i, CensusID = 8:9, DBH = c(1.5, 1.6), ExactDate = census_date(8:9), species = "x")))
  new_trees <- rbindlist(lapply(1:3, function(i) data.table(Tag = sprintf("N%04d", i), TrueStemID = 200L + i, CensusID = 8:9, DBH = c(1.4, 1.5), ExactDate = census_date(8:9), species = "x")))
  x <- rbind(x, new_in_tree, new_trees)
  T78 <- as.numeric(census_date(8L) - census_date(7L)) / 365.25
  T89 <- as.numeric(census_date(9L) - census_date(8L)) / 365.25
  b <- fit(x)
  r <- b$recruitment
  expect_equal(r$lambda_unit, "tree")
  expect_equal(r$n_new_stems_established_trees, 4L)          # the 3 stems of new trees are not new stems of established trees
  expect_equal(r$established_tree_years, 10 * T78 + 13 * T89, tolerance = 1e-9)
  expect_equal(r$lambda, (4 + 0.5) / (10 * T78 + 13 * T89), tolerance = 1e-9)
  expect_equal(r$lambda_slot, 7 / (7 * T78), tolerance = 1e-9) # legacy: 7 recruits over 7 empty slots
  b2 <- fit(x, recruit_rate_unit = "slot")
  expect_equal(b2$recruitment$lambda, r$lambda_slot)
  expect_error(fit(x, recruit_rate_unit = "plot"))
})
