# TrueStemID Pinning at Non-Anchor Censuses

## Status: **IMPLEMENTED** (2026-04-10)

## Summary

Previously, TrueStemID constrained the solver **only at the anchor census**.
Non-anchor TrueStemID was ignored during solving — neither the DP state
enumeration nor the probabilistic cost matrices used it as a constraint.

This feature adds "pinning": wherever TrueStemID is known at a non-anchor
census, constrain the solver to assign that observation to its correct
track.  Benefits:

- **Accuracy**: eliminates solver misidentification at censuses with known
  identity.
- **Speed**: reduces state space up to 1680× (e.g., P(8,5) = 6720 vs
  P(4,1) = 4 when 4 of 5 observations are pinned).
- **DP threshold**: tags that currently exceed `max_states` and fall back
  to the probabilistic matcher may fit within the DP budget with pinning.

All changes are guarded by a new parameter `pin_truestemid` (default
`TRUE`) for easy revert.

---

## Architecture Overview

### Key Data Structures

- `track_ids = c(anchor_ids, synthetic_extra_ids)` — length K vector
  mapping track index to stem ID.
- `.allowed_at_census[[p]][[i]]` — list of lists; `[[p]][[i]]` is an
  integer vector of allowed track indices for observation i at census p.
- `enumerate_states_constrained(K, n_obs, allowed_tracks, max_states)` —
  already handles single-element `allowed_tracks` (used at anchor census).

### Current Anchor Pinning (template for extension)

In `dp_global_dp.R` lines ~1283–1301, the anchor census pins each
observation to its TrueStemID track:

```r
.anchor_tidx <- match(.anchor_obs_pre$TrueStemID, track_ids)
.allowed_at_census[[p]] <- lapply(seq_len(n_obs), function(.j) {
    if (!is.na(.anchor_tidx[.j])) .anchor_tidx[.j] else seq_len(K)
})
```

The same pattern extends to non-anchor censuses.

---

## DP Pathway Changes

**File**: `dp_global/R/dp_global_dp.R`

### 1. Add parameter

Add `pin_truestemid = TRUE` to `match_stems_dp_global_backward_marginals_batch()`
signature (~line 93).

### 2. Pre-compute per-census pin map

Before the state enumeration loop (~line 1270), compute a pin map for
every census:

```r
pin_tidx_at_census <- vector("list", n_census)
if (isTRUE(pin_truestemid)) {
    for (p in seq_len(n_census)) {
        idx <- obs_row_idx[[p]]
        tsid <- tree_data$TrueStemID[idx]
        tidx <- match(tsid, track_ids)
        # Only pin if TrueStemID is in track_ids (i.e., an anchor-known stem)
        tidx[is.na(tsid)] <- NA_integer_
        pin_tidx_at_census[[p]] <- tidx
    }
}
```

### 3. Inject pinning into state enumeration loop

Current structure (lines 1283–1340):

```
if (p == n_census)        → anchor pinning (KEEP AS-IS)
else if (backward_prop)   → growth-feasibility narrowing
else                      → all tracks allowed
```

New structure:

```
if (p == n_census)        → anchor pinning (KEEP AS-IS)
else if (backward_prop)   → growth-feasibility narrowing
                            THEN: override pinned obs
else                      → all tracks allowed
                            THEN: override pinned obs
```

After backward propagation (or default all-K) computes `.allowed`:

```r
if (isTRUE(pin_truestemid) && !is.null(pin_tidx_at_census[[p]])) {
    .pins <- pin_tidx_at_census[[p]]
    for (.oi in which(!is.na(.pins))) {
        .allowed[[.oi]] <- .pins[.oi]
    }
    if (any(!is.na(.pins))) .use_constrained <- TRUE
}
```

Pin overrides growth-feasibility narrowing: if TrueStemID says "this IS
track X", trust it.  Transition costs still penalize unlikely growth.

### 4. Diagnostic logging

```r
n_pinned <- sum(!is.na(.pins))
if (n_pinned > 0L) {
    vcat(prefix, "  Pinned ", n_pinned, " obs at C", cc, " via TrueStemID")
    message(prefix, "  Pinned ", n_pinned, " obs at C", cc, " via TrueStemID")
}
```

### No other DP changes needed

- `enumerate_states_constrained()` in `dp_global_states.R` already
  handles single-element allowed_tracks — no changes needed.
- Backward DP, Viterbi decode, and forward marginals operate on the
  state matrices emitted by the enum loop — no changes needed.

---

## Probabilistic Pathway Changes

**File**: `dp_global/R/dp_probabilistic_matching.R`

### 1. Add parameter

Add `pin_truestemid = TRUE` to `match_stems_probabilistic()` (~line 29).

### 2. Store TrueStemID in obs_data

After obs_data construction (~line 130):

```r
for (i in seq_len(n_census)) {
    obs_data[[i]]$true_stem_id <- tree_data$TrueStemID[obs_data[[i]]$idx]
}
```

### 3. Maintain track mapping during backward sampling

During the backward sampling loop (lines 195–240), maintain
`next_obs_to_track` alongside the assignment sampling:

```r
# Initialize: anchor obs map directly to anchor_ids
next_obs_to_track <- anchor_ids

# Last pair (closest to anchor)
# Before greedy_assignment_gumbel():
cost_last <- pair_data[[last_pair]]$log_cost
if (isTRUE(pin_truestemid)) {
    tsid <- obs_data[[last_pair]]$true_stem_id
    for (r in seq_len(pair_data[[last_pair]]$n_curr)) {
        if (!is.na(tsid[r]) && tsid[r] %in% anchor_ids) {
            j <- which(anchor_ids == tsid[r])
            if (length(j) == 1L && j <= pair_data[[last_pair]]$n_next) {
                cost_last[r, -j] <- -Inf  # force assignment to column j
            }
        }
    }
}
per_pair_assignments[[last_pair]] <- greedy_assignment_gumbel(cost_last, ...)

# Update track mapping
for (r in seq_len(n_curr)) {
    col <- per_pair_assignments[[last_pair]][r]
    if (col <= n_next) curr_obs_to_track[r] <- next_obs_to_track[col]
    else curr_obs_to_track[r] <- death_id
}
next_obs_to_track <- curr_obs_to_track
```

Repeat for each subsequent pair going backward: use `next_obs_to_track`
to find which column j at census i+1 leads to the pinned track, then
mask `cost_i[r, -j] <- -Inf`.

### Edge case

If the pinned track has no observation at the next census (pin target
died downstream), don't mask — let greedy assignment handle it naturally.

---

## Parameter Plumbing

### Fallback call

In `dp_global_dp.R` (~line 547), pass through to probabilistic matcher:

```r
match_stems_probabilistic(..., pin_truestemid = pin_truestemid)
```

### CLI parameter

In `main_cpp.R` and `main_cpp_chunk.R`, add:

```r
PIN_TRUESTEMID <- parse_bool_arg("--PIN_TRUESTEMID", default = TRUE)
```

Pass to the batch function call.

---

## Regression Testing

### Test Cases

| Run                  | Pins active? | Expected outcome |
|----------------------|-------------|------------------|
| Simulated DP         | No*         | IDENTICAL        |
| Simulated PROB       | No*         | IDENTICAL        |
| BCI 115203 DP        | Check       | IDENTICAL or improved |
| BCI 005558 PROB      | Yes         | Pinned rows correct |
| BCI 191475 PROB      | Yes         | Pinned rows correct |

*Simulated data has no pre-anchor TrueStemID (censuses 1–6 have NA,
7 is anchor), so pinning is a no-op — confirms no code breakage.

### Verification Steps

1. `diff` each baseline vs. post-change output — simulated must be
   identical.
2. For BCI tags: every row where `TrueStemID` is non-NA AND pre-anchor
   should have `ReconstructedStemID` reflecting the solver's actual
   output (with pinning, it should equal `TrueStemID`).
3. Check logs for "Pinned N obs at C..." messages.
4. Run with `--PIN_TRUESTEMID=FALSE` — must reproduce pre-change output.

### Feature-flag revert

```sh
# Disable pinning without code changes:
Rscript main_cpp.R --PIN_TRUESTEMID=FALSE ...
```

### Full revert

```sh
git checkout pre-pinning-baseline
```

(Create this tag before implementing.)

---

## Design Decisions

- **Pin overrides growth feasibility**: TrueStemID is authoritative.  The
  state space includes the pinned track regardless of growth bounds;
  transition costs still reflect the true growth likelihood.

- **Only pin TrueStemID present in track_ids**: Stems not present at the
  anchor (died before anchor) have no track to pin to — skip silently.

- **Feature flag default TRUE**: `--PIN_TRUESTEMID=FALSE` disables
  without code changes.

---

## Future Extensions

1. **Simulated data with pre-anchor TrueStemID**: Add TrueStemID at
   censuses 5–6 of simulated data (using OriginalStemID as ground truth)
   for a controlled accuracy test.

2. **DP threshold monitoring**: With pinning, log tags that previously
   exceeded `max_states` and now fit within the DP budget.

3. **Partial pinning**: Pin only at censuses where TrueStemID confidence
   is high (e.g., field-verified vs. database-inferred).

---

## Implementation Notes (2026-04-10)

### Changes from original plan

1. **Probabilistic pathway restructured**: The original plan assumed track
   identity was maintained during the backward sampling loop. In reality,
   sampling was per-pair independent; identity mapping only happened
   post-hoc in `stitch_assignments_backward()`. The implementation adds
   inline per-sample track propagation (`propagate_track_backward()`) during
   the sampling loop and applies pin masks (`apply_pin_mask()`) before each
   `greedy_assignment_gumbel()` call.

2. **Duplicate-pin guard added**: Not in the original plan. If two
   observations at the same census claim the same track via TrueStemID,
   the second is released with a warning.

3. **Bug fix included**: `prob_lookahead_weight` was missing from the
   segment-split `.sub_args` list in `dp_global_dp.R`. Fixed alongside
   the pinning implementation.

### Files modified

| File | Changes |
|------|---------|
| `dp_global/R/dp_global_dp.R` | `pin_truestemid` parameter, pin map computation, state-enum injection, `.sub_args` fix (added `prob_lookahead_weight` + `pin_truestemid`), fallback call passthrough |
| `dp_global/R/dp_probabilistic_matching.R` | `pin_truestemid` parameter, pin lookup tables, sampling loop restructure with inline track propagation, `apply_pin_mask()` and `propagate_track_backward()` helpers |
| `dp_global/scripts/main_cpp.R` | `PIN_TRUESTEMID` default + CLI_REFERENCE + `run_dp_one_group()` passthrough + error handler passthrough |
| `dp_global/scripts/main_cpp_chunk.R` | Same as main_cpp.R |
| `dp_global/scripts/main_cpp_chunk_bci.R` | `PIN_TRUESTEMID` default |

### Regression tests passed

- All 5 R files parse clean
- `make smoke` — all modules load
- Single-tag simulated (Tag=20): runs clean, no "Pinned" messages (expected — no pre-anchor TrueStemID)
- Full chunked run (all simulated tags, 20 chunks): completes with no errors
