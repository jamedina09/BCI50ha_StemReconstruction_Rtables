# STEM_IDENTIFICATION_TEST

## Overview

Biologically informed dynamic-programming (DP) solver that reconstructs stem identities across forest censuses backward in time from a known anchor census. Given multi-stem tree measurements and a late-census anchor with trusted `TrueStemID`, the algorithm assigns each earlier observation to a latent identity track by minimising negative log-likelihood costs that encode growth, mortality, and recruitment biology. Uncertainty is quantified via forward-backward marginals and optional posterior sampling.

When the exact DP state space is too large (combinatorial explosion from many stems per tree), the solver automatically falls back to a **probabilistic greedy matching** module that uses the same biological cost model and hard pruning bounds with Gumbel-noise stochastic sampling and **sequential backward conditioning** to produce approximate reconstructions with per-observation posterior probabilities. The probabilistic matcher applies a two-layer **sample-level growth repair** before computing marginals: (1) hard-rate bound enforcement and (2) **measurement-error-informed cumulative-shrinkage detection** that severs runs of small decreases exceeding the ME threshold — mirroring the DP's global cost accumulation that naturally penalises consecutive shrinkage, but adapted for the per-pair greedy context. At the anchor census, known `TrueStemID` identities are always preserved as hard constraints; pre-anchor observations with `TrueStemID` are solved by the DP and can receive different identity assignments based on biological likelihood. The posterior path samples from both methods are used downstream for **basal area uncertainty quantification**.

## How the Two Algorithms Work (Plain-Language Summary)

### The Problem

In long-term forest census plots, individual trees can have multiple stems. Each stem is measured every few years, but **stem identity labels are only reliable at one late census** (the "anchor"). For all earlier censuses, we need to figure out which measurement belongs to which stem — a problem that gets harder when stems die, new ones appear, or measurement errors create confusing size sequences.

### Exact DP (Dynamic Programming) Solver

Think of the DP solver like solving a jigsaw puzzle backward from a finished picture. At the anchor census we know exactly which measurement belongs to which stem. Working backward one census at a time, the DP **evaluates every possible way** to connect earlier measurements to the known stem identities. Each candidate connection is scored using biology: how fast trees actually grow, how likely a stem is to die or a new one to appear, and how noisy the measurement tools are. The algorithm picks the single assignment with the **best overall score** across the entire history — not just one census at a time, but jointly optimised over all censuses at once.

Because it examines every possibility, the DP always finds the mathematically optimal answer. The downside is that the number of possibilities explodes factorially with the number of stems: a tree with 3 stems per census has a manageable puzzle, but a tree with 7+ stems per census has billions of combinations — too many to enumerate.

### Probabilistic Greedy Matcher (Fallback)

When there are too many stems for the DP to try every combination, the probabilistic matcher takes over. Instead of exhaustive enumeration, it uses a **sampling strategy**: it generates hundreds of plausible random assignments (using Gumbel-noise perturbations of the same biological scores), stitches each sample backward from the anchor to build full identity histories, removes any links that violate growth constraints, and then **takes a vote** across all surviving samples. The identity that wins the most votes for each observation becomes the final assignment, and the vote share becomes the posterior probability.

To compensate for its greedy per-pair nature (which lacks the DP's global cost accumulation), the probabilistic matcher applies two extra safeguards: a hard-rate bound check and a measurement-error-informed cumulative-shrinkage detector that severs runs of small decreases that collectively exceed what measurement noise could plausibly explain. Both the bio hard gates (`Bio_Max_Shrink`, `Bio_Max_Growth`) in pairwise edge construction and the ME cumulative-shrinkage check in trajectory repair can be selectively disabled via the `USE_BIO_HARD_SHRINK_IN_PROB` and `USE_BIO_HARD_GROWTH_IN_PROB` flags — useful when a known large-shrinkage or large-growth event (e.g., storm damage) would otherwise cause the matcher to artificially split a continuous stem.

### When Each Algorithm Runs

The choice is made **per tag** (i.e., per tree) and is controlled by the `DP_MAX_STATES` parameter:

| Scenario | Algorithm | Why |
|----------|-----------|-----|
| Tree has ≤ 6 stems per census (default settings) | **Exact DP** | State space fits within `DP_MAX_STATES` (40,000) — all combinations are enumerable |
| Tree has 7+ stems in any census | **Probabilistic matcher** | Factorial explosion exceeds the state budget — DP would run out of memory or time |
| Certain species or growth forms (configured via `PROB_SPECIES` / `FALLBACK_GROWTH_FORMS`) | **Probabilistic matcher** | User routes specific groups to the fallback regardless of stem count |
| DP solver hits a runtime error (e.g., memory) | **Probabilistic matcher** | Automatic error-recovery fallback |

### Can Both Run on the Same Tag?

For a simple tag with no resprout events, **only one** algorithm produces the reconstruction. However, when a tag contains an **R event** (a resprout or breakage code such as R, RP, RF, RT, QR, or OR), the tag is split into a pre-resprout segment and a post-resprout segment, each solved independently. Because the two segments have different stem counts and therefore different state-space sizes, **each segment chooses its algorithm independently** — so it is possible for the post-resprout segment to be solved by the exact DP while the pre-resprout segment falls back to the probabilistic matcher (or vice versa). Both algorithms' outputs are then combined into a single reconstruction for the tag.

Across a full dataset, most tags will use one algorithm throughout; only tags with R events can have mixed-method output rows. The two algorithms use **identical biological models** and hard pruning bounds, so their outputs are directly comparable and mergeable regardless of which segments used which method. Both produce the same output schema including posterior probabilities and optional posterior path samples for downstream uncertainty quantification.

## Directory Layout

```
├── dp_global/                 # Core algorithm, drivers, and C++ acceleration
│   ├── R/                     # R modules (sourced by dp_global_main.R)
│   │   ├── dp_global_main.R           # Module loader (sources all R modules in order)
│   │   ├── dp_global_bio.R            # Biological parameter estimation
│   │   ├── dp_global_states.R         # State enumeration & track-DBH helpers
│   │   ├── dp_probabilistic_matching.R # Probabilistic greedy matching fallback
│   │   ├── dp_global_dp.R            # Core DP solver (backward/forward pass, marginals)
│   │   ├── dp_global_utils.R         # Shared utilities
│   │   ├── dp_global_diag.R          # Diagnostics & PDF plotting
│   │   ├── naming_helpers.R          # Output directory naming
│   │   ├── complexity/               # DP complexity estimator
│   │   └── dpglobal_bundle/          # Portable deployment bundle builder
│   ├── scripts/               # CLI driver scripts (main_cpp*.R) + basal area uncertainty
│   └── src/                   # C++ transition cost (Rcpp)
├── data_simulation/           # Simulated forest-census data generator
│   └── data/                  # Generated test datasets (CSV)
├── bci_data/                  # BCI census data (not tracked by git)
└── Makefile                   # Convenience targets (smoke test)
```

**Not tracked by git** (see `.gitignore`): `dp_global/output/`, `dp_global/ForestGEO_codes/`, `dp_global/examples/`, `bci_data/`, `*.rds`, `*.pdf`, `*.log`.

## Prerequisites

R ≥ 4.0 with packages: `data.table`, `Rcpp`, `here`.
Optional: `ggplot2`, `cowplot` (plotting), `arrow` (feather output), `withr` (bundle sourcing).

## Quickstart

```bash
# Verify modules load
make smoke

# Single-tag run on simulated data
Rscript dp_global/scripts/main_cpp.R --WHICH_TAG=20

# Full run (all tags, chunked output with resume support)
Rscript dp_global/scripts/main_cpp_chunk.R --RUN_ALL_TAGS=TRUE --DP_CHUNK_SIZE=7

# BCI data (single tag)
Rscript dp_global/scripts/main_cpp_bci.R --WHICH_TAG=123375
```

Key CLI parameters for controlling solver behavior:

| Flag | Default | Purpose |
|------|---------|---------|
| `--DP_MAX_STATES` | `40000` | Max injective states per census before probabilistic fallback |
| `--PROB_N_SAMPLES` | `200` | Number of Gumbel-noise samples for probabilistic matching |
| `--PROB_LOOKAHEAD_WEIGHT` | `0.5` | Sequential backward conditioning weight (0 = disabled) |
| `--POSTERIOR_SAMPLES` | `200` | Number of posterior path samples (0 to disable) |
| `--USE_BIO_HARD_SHRINK_IN_PROB` | `TRUE` | Apply `Bio_Max_Shrink` hard gate and ME cumulative-shrinkage check in probabilistic matcher. Set `FALSE` to allow shrinkage beyond the bio bound (soft penalty only; useful for confirmed large-shrinkage events) |
| `--USE_BIO_HARD_GROWTH_IN_PROB` | `TRUE` | Apply `Bio_Max_Growth` hard gate in probabilistic matcher. Set `FALSE` to allow growth beyond the bio bound (soft penalty only) |

### Understanding `DP_MAX_STATES`

`DP_MAX_STATES` controls the maximum number of assignment states the DP solver will enumerate at any single census before falling back to the probabilistic greedy matcher. It also controls the inter-census transition budget: the cross-product of states between any two adjacent censuses must not exceed `DP_MAX_STATES²`.

#### How states are counted

At each census, the DP enumerates all injective (one-to-one) assignments of $n$ observed stems to $K$ identity tracks. The number of such assignments is the falling factorial:

$$P(K, n) = K \times (K-1) \times \cdots \times (K - n + 1) = \frac{K!}{(K-n)!}$$

where $K$ is the number of tracks (determined by the anchor stem count, births needed, and slack). Typically $K = n + 1$ (with `slack_tracks = 1` and no births) or $K = n + 2$ (with 1 birth track added).

This grows **factorially**, so even modest increases in stem count cause explosive growth in the state space.

#### Two fallback triggers

1. **Per-census enumeration (`enum_exceeded`):** If $P(K, n)$ exceeds `DP_MAX_STATES` at any single census, the solver cannot enumerate states and falls back.
2. **Inter-census transitions (`edge_count_exceeded`):** If $P(K, n_1) \times P(K, n_2)$ exceeds `DP_MAX_STATES`² for any pair of adjacent censuses, the transition matrix is too large and the solver falls back.

In practice, the per-census limit is reached first because the state counts grow so rapidly.

#### Fallback thresholds by `DP_MAX_STATES` value

The tables below show when fallback occurs for different `DP_MAX_STATES` values. "Stems observed" is the number of stems with non-NA DBH in a single census. $K = n + 1$ assumes `slack_tracks = 1` with no birth tracks needed.

**With $K = n + 1$ (minimum realistic tracks):**

| Stems ($n$) | Tracks ($K$) | States $P(K,n)$ | 1,000 | 20,000 | 40,000 |
|:-:|:-:|--:|:-:|:-:|:-:|
| 2 | 3 | 6 | DP | DP | DP |
| 3 | 4 | 24 | DP | DP | DP |
| 4 | 5 | 120 | DP | DP | DP |
| 5 | 6 | 720 | DP | DP | DP |
| 6 | 7 | 5,040 | fallback | DP | DP |
| 7 | 8 | 40,320 | fallback | fallback | fallback |
| 8 | 9 | 362,880 | fallback | fallback | fallback |

**With $K = n + 2$ (when 1 birth track is needed):**

| Stems ($n$) | Tracks ($K$) | States $P(K,n)$ | 1,000 | 20,000 | 40,000 |
|:-:|:-:|--:|:-:|:-:|:-:|
| 2 | 4 | 12 | DP | DP | DP |
| 3 | 5 | 60 | DP | DP | DP |
| 4 | 6 | 360 | DP | DP | DP |
| 5 | 7 | 2,520 | DP | DP | DP |
| 6 | 8 | 20,160 | fallback | DP | DP |
| 7 | 9 | 181,440 | fallback | fallback | fallback |
| 8 | 10 | 1,814,400 | fallback | fallback | fallback |

**Summary — maximum stems per census handled by exact DP:**

| `DP_MAX_STATES` | `max_edges` (= `DP_MAX_STATES`²) | Max stems ($K = n+1$) | Max stems ($K = n+2$) |
|--:|--:|:-:|:-:|
| 1,000 | 1,000,000 | 5 | 5 |
| 20,000 | 400,000,000 | 6 | 6 |
| 40,000 | 1,600,000,000 | 6 | 6 |

**Key insight:** With the default `DP_MAX_STATES = 40,000`, the DP handles tags with up to **6 observed stems per census** exactly. Tags with **7 or more stems** in any census are routed to the probabilistic greedy matcher. Increasing `DP_MAX_STATES` to 50,000 would not help — the next factorial step (40,320 for 7 stems with $K=8$) requires `DP_MAX_STATES ≥ 40,321` AND the inter-census product must fit, which it does since $40{,}320^2 = 1.6 \times 10^9 < 40{,}321^2$.

#### How to choose a value

```r
# 1. Find the most complex tags in your dataset
library(data.table)
dt <- fread("your_data.csv")
obs_per_census <- dt[!is.na(DBH), .N, by = .(Tag, CensusID)]
max_obs <- obs_per_census[, .(max_n = max(N)), by = Tag][order(-max_n)]
head(max_obs, 10)  # top 10 most complex tags

# 2. Compute states for a specific stem count
n <- 6   # max observed stems in any census
K <- 8   # n + 2 (slack + 1 birth)
states <- prod(K:(K - n + 1))  # P(8, 6) = 20,160
cat("States:", states, "\n")

# 3. Set DP_MAX_STATES above that to guarantee exact DP
# Rscript dp_global/scripts/main_cpp_chunk.R --DP_MAX_STATES=25000
```

**Trade-off:** Higher values → exact DP for more tags (slower, more memory). Lower values → more tags use the probabilistic fallback (faster, approximate but uses the same biological model and pruning bounds).

See `dp_global/README.md` for the full algorithm description and `dp_global/scripts/README.md` for the CLI flag reference.

## Key Documentation

| Document | Contents |
|----------|----------|
| `dp_global/README.md` | Algorithm details, cost model, data requirements, parameter estimation, fallback mechanisms |
| `dp_global/scripts/README.md` | CLI flags, chunking, resume, example invocations, basal area uncertainty |
| `dp_global/src/README.md` | C++ acceleration API and validation |
| `data_simulation/README.md` | Simulation parameters, biological models, output format |

## Reconstruction Methods

Each observation in the output receives a `ReconstructionMethod` label indicating how its `ReconstructedStemID` was determined:

| Method | Description |
|--------|-------------|
| `given` | Identity known from input `TrueStemID` (anchor, pre-anchor pin, post-anchor row, or hard-invariant sweep) |
| `dp` | Assigned by the exact DP solver |
| `probabilistic` | Assigned by the probabilistic greedy matching fallback |
| `provisional_dp` | Provisional anchor assigned by DP |
| `dp_mf_inferred` | Missing-from-field census identity inferred from flanking DP assignments |
| `carried_terminal` | Orphan terminal-event row (`Status` ∈ {`dead`, `stem dead`, `broken below`}, `DBH = NA`, engine returned `NA`) backfilled by post-engine LOCF from the most recent prior `ReconstructedStemID` in the same `(Tag, OriginalStemID)` group. Applied uniformly across all driver scripts via `apply_carried_terminal_backfill()` in `dp_global/R/dp_global_main.R`. |
| `none_after_anchor` | Post-anchor row without assignment |
| `skipped_no_data` | Tag had no usable data for reconstruction |

### Hard-invariant sweep and the `SweepAuditOverride` audit column

When `PIN_TRUESTEMID = TRUE` (default), every output row with a non-NA `TrueStemID` is forced to `ReconstructedStemID = TrueStemID` and `ReconstructionMethod = "given"` by an idempotent sweep that runs at three sites: inside `finalize_out()` (DP path), inside `match_stems_probabilistic()` (probabilistic fallback), and at the script level in `run_dp_one_group()`. This sweep guarantees the invariant even on rows the DP never visits (NA-DBH terminal rows anchored by the pre-DP propagation in `main_cpp_bci.R` Steps 2/3, MF re-insertion edge cases, and probabilistic-fallback leaks).

In the rare case where the DP or probabilistic engine had already assigned a non-NA `ReconstructedStemID` that disagrees with `TrueStemID`, the sweep silently overrides it. Those rows are flagged in the **`SweepAuditOverride`** logical column. Downstream consumers that propagate posterior uncertainty should treat any row where `SweepAuditOverride == TRUE` as observed (P=1, entropy=0) — the values in the `DP_PosteriorTop*` / `DP_PosteriorReconstructedProb` columns describe the *engine's* (overridden) choice, not the final `ReconstructedStemID`. For all other `given` rows the posterior collapses naturally to the pinned ID and the posterior columns are consistent with the final assignment.

The engine's pre-sweep choice is preserved alongside the final value in the **`ReconstructedStemID_PreSweep`** column. It equals `ReconstructedStemID` everywhere `SweepAuditOverride == FALSE` and carries the original engine-assigned ID where `SweepAuditOverride == TRUE`. The column is populated once by the first sweep layer that fires (engine `finalize_out` → probabilistic matcher → script-level backstop) and preserved unchanged by later sweeps, so it always reflects the pre-sweep state regardless of which code path produced the row.

## Conventions

- Driver scripts: `dp_global/scripts/` — `main_cpp.R` (single-tag/small), `main_cpp_chunk.R` (large chunked), `main_cpp_bci.R` (BCI-specific with `withr` bundle sourcing), `basal_area_uncertainty.R` (posterior BA uncertainty).
- Module internals: `dp_global/R/` — sourced in order by `dp_global_main.R`.
- Output directories: auto-created under `dp_global/output/` (not tracked by git).

## Uncertainty Estimation Workflow

After running the reconstruction pipeline, posterior path samples can be used to quantify uncertainty in derived quantities such as basal area:

```bash
# 1. Run reconstruction pipeline
Rscript dp_global/scripts/main_cpp_chunk.R --POSTERIOR_SAMPLES=250

# 2. Quantify basal area uncertainty from posterior paths
Rscript dp_global/scripts/basal_area_uncertainty.R \
  --RUN_DIR=dp_global/output/<run_dir>
```

The BA uncertainty script reads the main reconstruction and the `posteriors/` directory, then produces two summary CSVs and a multi-page PDF of diagnostic figures:

- `basal_area_tag_census.csv` — per-tag per-census total BA (m²) and stem count.
- `basal_area_tag_change.csv` — per-tag BA change between consecutive censuses, decomposed into **growth** (survivors), **loss** (mortality), and **gain** (recruitment), with posterior means, SDs, and 95% credible intervals. All BA values in m².
- `basal_area_figures.pdf` — per-tag panels (BA trajectory, stem count, decomposition bars with uncertainty whiskers, stem demographics), uncertainty histograms, and posterior density plots (kernel densities of Growth/Loss/Gain/DeltaBA pooled across all intervals and per census interval, with weighted-mean vertical lines).

Tag-level total BA per census is **invariant** to identity assignment — the same DBH values sum identically regardless of which stem they are assigned to. The decomposition into growth, loss, and gain components **is** identity-dependent and is where posterior uncertainty manifests. See `dp_global/scripts/README.md` for details.
