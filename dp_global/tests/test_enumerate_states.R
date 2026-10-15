# Tests for the state enumeration of the DP (dp_global/R/dp_global_states.R):
# enumerate_states_constrained() and enumerate_states_injective(). Above
# max_states both must return NULL, so that the DP sends the tag to the
# probabilistic matcher (enum_exceeded) instead of solving it on a truncated
# set of states.
# Run from the project root:
#   testthat::test_file("dp_global/tests/test_enumerate_states.R")
suppressMessages(library(testthat))
root <- c(".", "..", "../..")
root <- normalizePath(root[file.exists(file.path(root, "dp_global", "R", "dp_global_states.R"))][1])
source(file.path(root, "dp_global", "R", "dp_global_states.R"))

# All injective assignments of n_obs observations to K tracks that respect
# allowed_tracks, one string per assignment (brute force)
brute_force <- function(K, n_obs, allowed_tracks) {
  grid <- as.matrix(expand.grid(allowed_tracks))
  grid <- grid[apply(grid, 1, function(r) !anyDuplicated(r)), , drop = FALSE]
  sort(apply(grid, 1, paste, collapse = ","))
}
as_keys <- function(mat) sort(apply(mat, 1, paste, collapse = ","))

test_that("below max_states every feasible assignment is returned once", {
  allowed <- list(c(1, 2, 3), c(2, 3), c(1, 3, 4))
  mat <- enumerate_states_constrained(4L, 3L, allowed, max_states = 100L)
  expect_equal(as_keys(mat), brute_force(4L, 3L, allowed))
  expect_false(anyDuplicated(as_keys(mat)) > 0)
})

test_that("exactly max_states feasible assignments are all returned; one more gives NULL", {
  allowed <- list(c(1, 2, 3), c(2, 3), c(1, 3, 4))
  n_feasible <- length(brute_force(4L, 3L, allowed))
  mat <- enumerate_states_constrained(4L, 3L, allowed, max_states = n_feasible)
  expect_equal(nrow(mat), n_feasible)
  expect_equal(as_keys(mat), brute_force(4L, 3L, allowed))
  # the defect fixed here: the first max_states assignments used to be
  # returned as if they were all of them
  expect_null(enumerate_states_constrained(4L, 3L, allowed, max_states = n_feasible - 1L))
})

test_that("a large constrained census above max_states gives NULL, not a truncated set", {
  K <- 9L
  n_obs <- 6L
  allowed <- c(list(1L), replicate(n_obs - 1L, seq_len(K), simplify = FALSE)) # one pinned observation
  # 8 * 7 * 6 * 5 * 4 = 6,720 feasible assignments
  expect_equal(nrow(enumerate_states_constrained(K, n_obs, allowed, max_states = 6720L)), 6720L)
  expect_null(enumerate_states_constrained(K, n_obs, allowed, max_states = 6719L))
  expect_null(enumerate_states_constrained(K, n_obs, allowed, max_states = 1000L))
})

test_that("edge cases: no observation, no feasible assignment, max_states = 0", {
  expect_equal(dim(enumerate_states_constrained(3L, 0L, list(), max_states = 10L)), c(1L, 0L))
  expect_null(enumerate_states_constrained(2L, 3L, list(1:2, 1:2, 1:2), max_states = 10L)) # n_obs > K
  expect_null(enumerate_states_constrained(3L, 2L, list(integer(0), 1:3), max_states = 10L)) # an observation without a track
  expect_null(enumerate_states_constrained(3L, 2L, list(2L, 2L), max_states = 10L)) # two observations, one track
  expect_null(enumerate_states_constrained(3L, 2L, list(1:3, 1:3), max_states = 0L))
  expect_error(enumerate_states_constrained(3L, 2L, list(1:3), max_states = 10L), "allowed_tracks must have length")
})

test_that("with every track allowed, the constrained and the injective enumeration agree, also above max_states", {
  K <- 5L
  n_obs <- 3L
  all_tracks <- replicate(n_obs, seq_len(K), simplify = FALSE)
  inj <- enumerate_states_injective(K, n_obs, max_states = 1000L) # P(5, 3) = 60
  con <- enumerate_states_constrained(K, n_obs, all_tracks, max_states = 1000L)
  expect_equal(as_keys(con), as_keys(inj))
  expect_null(enumerate_states_injective(K, n_obs, max_states = 59L))
  expect_null(enumerate_states_constrained(K, n_obs, all_tracks, max_states = 59L))
})
