# Tests for the birth-death mode of the probabilistic matcher
# (augment_cost_matrix(), fallback_log_cost(), greedy_assignment_gumbel(),
# match_stems_probabilistic(birth_death = ...)) and the exact assignment solver
# hungarian_min_rcpp() (dp_global/src/transition_cost_rcpp.cpp).
# Run from the project root:
#   testthat::test_file("dp_global/tests/test_probabilistic_birth_death.R")
suppressMessages({
  library(testthat)
  library(data.table)
})
# dp_global_main.R finds its modules from the working directory (project root);
# testthat runs this file from dp_global/tests
root <- c(".", "..", "../..")
root <- normalizePath(root[file.exists(file.path(root, "dp_global", "R", "dp_global_main.R"))][1])
.owd <- setwd(root)
invisible(capture.output(suppressMessages(source(file.path(root, "dp_global", "R", "dp_global_main.R")))))
setwd(.owd)

# Bio parameters of a palm-like species (Oenocarpus-like: slow growth, recruits ~7 cm)
bio <- list(mu_const = 0.567, mu_gamma = -0.231, sigma0 = 0.095, sigma1 = 0, max_shrink = -0.5, k_shrink = 0,
  max_growth_bio = 5, k_growth = 0, h0 = 0.0442, beta_mort = -0.0127, recruit_meanlog = 1.95, recruit_sdlog = 0.24,
  recruit_max_dbh = 26, recruit_lambda = 0.166)

test_that("hungarian_min_rcpp finds the optimal assignment", {
  brute <- function(w) {
    n <- nrow(w)
    perms <- function(v) if (length(v) <= 1) list(v) else do.call(c, lapply(seq_along(v), function(i) lapply(perms(v[-i]), function(p) c(v[i], p))))
    min(sapply(perms(seq_len(n)), function(p) sum(w[cbind(seq_len(n), p)])))
  }
  set.seed(1)
  for (k in 1:150) {
    n <- sample(1:6, 1)
    w <- matrix(rnorm(n * n), n)
    a <- hungarian_min_rcpp(w)
    expect_equal(sort(a), seq_len(n))
    expect_equal(sum(w[cbind(seq_len(n), a)]), brute(w), tolerance = 1e-9)
  }
  expect_error(hungarian_min_rcpp(matrix(c(0, Inf, 1, 2), 2)), "finite")
})

test_that("birth-death matrix: one death cell per current stem, one recruit cell per next stem", {
  d0 <- c(9.5, 8.0); d1 <- c(7.6, 8.1, 8.0)
  L <- compute_pairwise_log_likelihood(d0, d1, 5, bio, -0.5, 5)
  A <- augment_cost_matrix(L, d0, d1, 5, bio, birth_death = TRUE)
  expect_equal(dim(A), c(5L, 5L))
  expect_equal(attr(A, "bd"), c(2L, 3L))
  expect_equal(A[1:2, 1:3], L)
  death <- A[1:2, 4:5]
  expect_true(all(is.finite(diag(death))) && all(!is.finite(death[row(death) != col(death)])))
  rec <- A[3:5, 1:3]
  expect_true(all(is.finite(diag(rec))) && all(!is.finite(rec[row(rec) != col(rec)])))
  expect_true(all(A[3:5, 4:5] == 0))
  # legacy mode is unchanged: K = max(n_curr, n_next), no "bd" attribute
  Al <- augment_cost_matrix(L, d0, d1, 5, bio)
  expect_null(attr(Al, "bd"))
  expect_equal(nrow(Al), 3L)
})

test_that("the birth-death sampler never uses a forbidden cell when an allowed assignment exists", {
  d0 <- c(10, 5); d1 <- c(40)   # 10 -> 40 and 5 -> 40 exceed 5 cm/yr; 40 cm exceeds the recruit cap
  L <- compute_pairwise_log_likelihood(d0, d1, 5, bio, -0.5, 5)
  A <- augment_cost_matrix(L, d0, d1, 5, bio, birth_death = TRUE)
  fb <- fallback_log_cost(A, compute_pairwise_log_likelihood(d0, d1, 5, bio, -Inf, Inf, use_bio_hard_shrink = FALSE, use_bio_hard_growth = FALSE), d1, 5, bio)
  set.seed(3)
  a <- greedy_assignment_gumbel(A, fallback = fb)
  # no assignment avoids every forbidden cell here: the least bad one is taken and flagged
  expect_false(is.null(attr(a, "n_forbidden")))
  d1 <- c(10.2)
  L <- compute_pairwise_log_likelihood(d0, d1, 5, bio, -0.5, 5)
  A <- augment_cost_matrix(L, d0, d1, 5, bio, birth_death = TRUE)
  for (s in 1:50) {
    a <- greedy_assignment_gumbel(A)
    expect_true(all(is.finite(A[cbind(seq_len(nrow(A)), a)])))
  }
})

# A palm clump: 3 stems in 2000, 3 in 2005 and 2010 (anchor, pinned). Between
# 2000 and 2005 the 9.5 cm stem dies and a new 7.6 cm stem appears; the other
# two barely change (palms do not grow in diameter).
clump <- function() {
  rows <- list(
    list(1L, 9.5, NA), list(1L, 8.0, NA), list(1L, 9.0, NA),
    list(2L, 7.6, NA), list(2L, 8.02, NA), list(2L, 9.03, NA),
    list(3L, 7.7, 101L), list(3L, 8.03, 102L), list(3L, 9.05, 103L)
  )
  x <- rbindlist(lapply(rows, function(r) data.table(CensusID = r[[1]], DBH = r[[2]], TrueStemID = r[[3]])))
  x[, `:=`(Tag = "P1", ExactDate = as.Date("2000-06-01") + round((CensusID - 1L) * 5 * 365.25),
    Bio_Mu_Growth = bio$mu_const, Bio_Gamma_Growth = bio$mu_gamma, Bio_Sigma0_Growth = bio$sigma0, Bio_Sigma1_Growth = bio$sigma1,
    Bio_Max_Shrink = bio$max_shrink, Bio_K_Shrink = 0, Bio_Max_Growth = 5, Bio_K_Growth = 0, Bio_H0_Mortality = bio$h0,
    Bio_Beta_Mortality = bio$beta_mort, Bio_Recruit_Meanlog = bio$recruit_meanlog, Bio_Recruit_Sdlog = bio$recruit_sdlog,
    Bio_Recruit_MaxDBH_unit = 26, Bio_Recruitment_lambda = bio$recruit_lambda)]
  x
}
run_clump <- function(bd) {
  out <- match_stems_probabilistic(clump(), min_growth = -0.5, max_growth = 5, anchor_start = 3L, n_samples = 100L,
    posterior_sample_seed = 7L, prob_lookahead_weight = 1, n_sigma_me = Inf, birth_death = bd, return_samples = TRUE)
  id <- function(c, d) out[CensusID == c & abs(DBH - d) < 1e-9, ReconstructedStemID]
  list(dying = id(1L, 9.5), recruit = id(2L, 7.6), stable = c(id(1L, 8.0) == id(2L, 8.02), id(1L, 9.0) == id(2L, 9.03)))
}

test_that("birth-death mode: a dying palm is not continued by the new stem; stable stems stay linked", {
  r <- run_clump(TRUE)
  expect_false(isTRUE(r$dying == r$recruit))
  expect_true(all(r$stable))
})

test_that("legacy mode: with stable stem counts the dying stem is forced to continue", {
  r <- run_clump(FALSE)
  expect_equal(r$dying, r$recruit)
  expect_true(all(r$stable))
})
