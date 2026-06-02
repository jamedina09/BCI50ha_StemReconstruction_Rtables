# BCI Stem Identification: Chunked Run Instructions

This guide explains how to launch and manage chunked stem-identification runs for BCI data using the provided R scripts. All commands should be run from the **project root** (the directory containing both `dp_global/` and `BCI_stem_reconstruction/`).

---

## Launch a New Full Run

To start a new chunked run, use the following command:

```sh
Rscript BCI_stem_reconstruction/2_STEM_IDENTIFICATION/1_main_cpp_chunk_bci.R \
  --DP_MAX_STATES=10000 \
  --PROB_SPECIES="oenoma,bactma,ficuob,ficupo,ficuc2,ficubu,ficuc1,ficuci,ficupe" \
  --DP_FALLBACK_GROWTH_FORMS="strangler" \
  --POSTERIOR_SAMPLE_SEED=42 \
  --MANUAL_CORES=TRUE \
  --MANUAL_CORES_VALUE=18 \
  --DP_CHUNK_SIZE=18 \
  --USE_MEASUREMENT_ERROR=FALSE \
  --BASE_OUT_DIR=/Users/medinaja/outputs_bci_stem_identification
```

**Flag explanations:**

- `DP_MAX_STATES`: Cap on DP state-space size per tag (higher = more accuracy, longer runtime).
- `PROB_SPECIES`: Comma-separated species codes to force probabilistic fallback (e.g., for strangler figs).
- `DP_FALLBACK_GROWTH_FORMS`: Growth forms that trigger fallback to the probabilistic matcher.
- `POSTERIOR_SAMPLE_SEED`: Integer seed for reproducible sampling.
- `MANUAL_CORES`/`MANUAL_CORES_VALUE`: Enable and set the number of parallel workers.
- `DP_CHUNK_SIZE`: Number of tags processed per parallel chunk (match to core count for efficiency).
- `USE_MEASUREMENT_ERROR`: Whether to use per-census DBH measurement error (set FALSE for speed).
- `BASE_OUT_DIR`: Root directory for all output (a timestamped subdirectory is created for each run).

After all chunks finish, run `2_merge_chunks_to_datatable.R` to merge per-chunk Feather files into the final dataset.

---

## Resume a Partial Run

If a run is interrupted, you can resume it. Chunks with a `_done.txt` marker are skipped automatically.

```sh
Rscript BCI_stem_reconstruction/2_STEM_IDENTIFICATION/1_main_cpp_chunk_bci.R \
  --OUT_DIR_OVERRIDE=/Users/medinaja/outputs_bci_stem_identification/<timestamp_run_code> \
  --DP_CHUNK_RESUME=TRUE \
  --DP_MAX_STATES=10000 \
  --PROB_SPECIES="oenoma,bactma,ficuob,ficupo,ficuc2,ficubu,ficuc1,ficuci,ficupe" \
  --DP_FALLBACK_GROWTH_FORMS="strangler" \
  --POSTERIOR_SAMPLE_SEED=42 \
  --MANUAL_CORES=TRUE \
  --MANUAL_CORES_VALUE=18 \
  --DP_CHUNK_SIZE=18 \
  --USE_MEASUREMENT_ERROR=FALSE \
  [--DP_CHUNK_START=<start>] [--DP_CHUNK_END=<end>]
```

- Replace `<timestamp_run_code>` with the actual directory name from your interrupted run.
- Optionally specify `DP_CHUNK_START` and/or `DP_CHUNK_END` to process a specific range of chunks.

---

## Utility: Check Maximum DP States for a Tag

To estimate whether a tag's stem count will fit within your `DP_MAX_STATES` budget, run this in R:

```r
n <- 5          # max observed stems in any census for the tag
K <- n + 2      # K = n + slack_tracks + 1 (default slack = 1)
states <- prod(K:(K - n + 1))   # permutations P(K, n)
cat("Max DP states for n =", n, ":", states, "\n")
# Example: n=5  →  P(7,5) = 2,520  (well within 10,000)
#          n=6  →  P(8,6) = 20,160 (exceeds 10,000 — DP falls back)
```

---

## Output Structure

- Each run creates a timestamped output directory inside `BASE_OUT_DIR`.
- Per-chunk Feather files are written to this directory.
- After all chunks finish, merge them using `2_merge_chunks_to_datatable.R`.

---

For more details, see comments in the driver scripts and the main project README files.
