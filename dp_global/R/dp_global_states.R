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
    n_states <- 1
    for (i in 0:(n_obs - 1L)) n_states <- n_states * (K - i)

    # Safety guard: if it's too big, we refuse and make the caller fall back.
    if (is.finite(n_states) && n_states > max_states) {
        return(NULL)
    }
    tracks <- seq_len(K)

    build <- function(prefix, remaining) {
        # Recursive constructor.
        #
        # prefix    = the track choices we have already made for obs 1..length(prefix)
        # remaining = which tracks are still available to use
        #
        # Base case:
        # - If we've chosen n_obs tracks, we have completed one full state.
        if (length(prefix) == n_obs) {
            return(list(prefix))
        }

        # Otherwise, we need to pick a track for the next observation.
        # We'll try each available track and recurse.
        out <- vector("list", 0L)
        for (t in remaining) {
            # Choose track t for the next observation.
            # Then remove t from remaining (cannot reuse it).
            out <- c(out, build(c(prefix, t), remaining[remaining != t]))
        }
        out
    }
    combos <- build(integer(0), tracks)
    do.call(rbind, lapply(combos, function(v) matrix(v, nrow = 1L)))
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
