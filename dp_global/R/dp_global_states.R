############################################################
# dp_global_states.R
# State enumeration & small state helpers
############################################################

enumerate_states_injective <- function(K, n_obs, max_states) {
    # PURPOSE
    # - Enumerate the DP "state space" for a single census: all injective mappings
    #   of n_obs labeled observations into K labeled tracks.
    #
    # INPUTS
    # - K: integer; number of tracks (identity slots).
    # - n_obs: integer; number of observed stems in this census.
    # - max_states: hard cap on number of states to enumerate.
    #
    # OUTPUT
    # - Matrix with one row per state and n_obs columns.
    #   Each row is a length-n_obs integer vector of track indices; e.g. c(2,5,1)
    #   means obs1->track2, obs2->track5, obs3->track1.
    # - Returns NULL when:
    #   - n_obs > K (impossible to assign injectively), or
    #   - the estimated number of states P(K, n_obs) exceeds max_states.
    #
    # COMPLEXITY
    # - Time and memory scale as O(P(K, n_obs) * n_obs).
    # - This can explode quickly; callers must provide a small max_states and fall
    #   back when enumeration is refused.
    if (n_obs == 0L) {
        return(matrix(integer(0), nrow = 1L, ncol = 0L))
    }

    # If we have more observations than tracks, it is impossible to assign
    # each observation to a unique track.
    if (n_obs > K) {
        return(NULL)
    }

    # Before we try to build the list, estimate how big it will be.
    # Number of injective mappings: P(K, n_obs) = K*(K-1)*...*(K-n_obs+1)
    n_states_est <- prod(seq.int(from = K, by = -1L, length.out = n_obs))

    # Safety guard: if it's too big, we refuse and make the caller fall back.
    if (!is.finite(n_states_est) || n_states_est > max_states) {
        return(NULL)
    }
    n_states <- as.integer(n_states_est)

    # Pre-allocate the output matrix and fill it using a shared row buffer + counter.
    # This avoids the O(n_states^2) list-copying that recursive c(out, build(...)) incurs.
    mat <- matrix(0L, nrow = n_states, ncol = n_obs)
    row_ctr <- 1L
    current_row <- integer(n_obs)

    fill <- function(col, avail) {
        if (col > n_obs) {
            mat[row_ctr, ] <<- current_row
            row_ctr <<- row_ctr + 1L
            return()
        }
        for (t in avail) {
            current_row[col] <<- t
            fill(col + 1L, avail[avail != t])
        }
    }
    fill(1L, seq_len(K))
    mat
}

enumerate_states_constrained <- function(K, n_obs, allowed_tracks, max_states) {
    # PURPOSE
    # - Like enumerate_states_injective(), but restricts each observation to a
    #   subset of tracks.  This eliminates provably infeasible assignments
    #   (e.g., growth bounds violation over the cumulative span to the anchor)
    #   BEFORE enumerating, dramatically reducing the state count when many
    #   tracks are biologically impossible for a given observation.
    #
    # INPUTS
    # - K: integer; number of tracks.
    # - n_obs: integer; number of observed stems.
    # - allowed_tracks: list of length n_obs.  allowed_tracks[[i]] is an integer
    #   vector of track indices that observation i may be assigned to.
    # - max_states: hard cap on enumerated states.
    #
    # OUTPUT
    # - Matrix with one row per feasible injective assignment, n_obs columns.
    # - Returns NULL when no feasible assignment exists or max_states exceeded.

    if (n_obs == 0L) {
        return(matrix(integer(0), nrow = 1L, ncol = 0L))
    }
    if (n_obs > K) return(NULL)
    if (length(allowed_tracks) != n_obs) {
        stop("allowed_tracks must have length n_obs (", n_obs, "), got ", length(allowed_tracks))
    }
    # Fast check: any obs with 0 allowed tracks means no solution
    if (any(lengths(allowed_tracks) == 0L)) return(NULL)

    # Pre-allocate output matrix and fill via recursive backtracking
    mat <- matrix(0L, nrow = max_states, ncol = n_obs)
    row_ctr <- 1L
    current_row <- integer(n_obs)
    avail <- rep(TRUE, K)  # tracks not yet used by earlier obs

    fill <- function(col) {
        if (col > n_obs) {
            if (row_ctr > max_states) return()
            mat[row_ctr, ] <<- current_row
            row_ctr <<- row_ctr + 1L
            return()
        }
        for (t in allowed_tracks[[col]]) {
            if (!avail[t]) next  # injectivity: track already used
            if (row_ctr > max_states) return()
            current_row[col] <<- t
            avail[t] <<- FALSE
            fill(col + 1L)
            avail[t] <<- TRUE
        }
    }
    fill(1L)

    n_actual <- row_ctr - 1L
    if (n_actual == 0L) return(NULL)
    if (n_actual > max_states) return(NULL)
    mat[seq_len(n_actual), , drop = FALSE]
}

state_key <- function(state_vec) {
    # PURPOSE
    # - Canonical string key for a census "state" (injective observation->track map).
    #
    # INPUT
    # - state_vec: integer vector of length n_obs; entries are track indices.
    #
    # OUTPUT
    # - Single string like "2,5,1" used to index DP tables (named vectors/lists).
    if (length(state_vec) == 0L) {
        return("")
    }
    paste(state_vec, collapse = ",")
}

state_to_track_dbh <- function(state_vec, obs_dbh, K) {
    # PURPOSE
    # - Convert a census "state" into a track-indexed DBH vector so transition costs
    #   can be computed track-by-track.
    #
    # INPUTS
    # - state_vec: integer vector of length n_obs mapping obs index -> track index.
    # - obs_dbh: numeric vector of length n_obs; DBH values for each observation.
    # - K: total number of tracks.
    #
    # OUTPUT
    # - Numeric vector of length K where unused tracks are NA and used tracks receive
    #   the DBH of the observation assigned to that track.
    out <- rep(NA_real_, K)
    if (length(state_vec) == 0L) {
        return(out)
    }
    out[state_vec] <- obs_dbh
    out
}
