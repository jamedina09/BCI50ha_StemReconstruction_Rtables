# Forest Census Data Simulation for Stem Identification Testing

This directory contains the `simulate_data.R` script that generates realistic forest census data for testing stem identification and reconstruction algorithms.

## Overview

Simulates biologically plausible multi-species tropical forest census data (growth, mortality, recruitment, measurement error, and stem ID masking) compatible with the DP workflow. The final dataset contains 136 tags:

- **42 simulated trees** (Tags 1–42): 38 randomly generated trees across 3 species + 4 problematic-stem variants for testing confusable growth trajectories
- **46 hardcoded diagnostic tags** (Tags 43–88 and BCI-derived tag numbers): Real-world-inspired multi-stem patterns extracted from BCI data for regression testing
- **3 M-code test tags** (Tags 901–903): Validate M-coded main-stem constraints in the DP solver
- **45 row-count invariant edge case tags** (Tags 9901–9945): Systematic coverage of every combination of census span, stem count, DBH availability, and special flags

Key features include:
- **Variable census intervals**: ~5 years between censuses with random noise
- **Date-stamped censuses**: Each census has an exact date starting from 1980-01-01
- **Organized code structure**: Script divided into clear sections for easy maintenance
- **Multi-species forests** with species-specific trait variation
- **Multi-stem trees** with recruitment and mortality processes
- **Realistic growth trajectories** with size-dependent growth and process variability
- **Measurement error models** matching the DP workflow
- **Flexible growth scaling** for simulating disturbances
- **Stem ID masking** to simulate ForestGEO protocols
- **Comprehensive diagnostics** with trajectory visualization

## Quickstart

Run the simulation:

```bash
Rscript simulate_data.R
```

Generated outputs are placed in `data/`.

## Output Files

The simulation generates several files in the `data/` subdirectory:

- `simulated_data_1.csv`: Main dataset with all stem measurements
- `simulated_data_1.pdf`: Species-level growth trajectory plots
- `simulated_data_tag_level_trajectories_1.pdf`: Tag-level individual stem trajectory plots

### Main Dataset: `simulated_data_1.csv`
Columns:
- `Species`: Species name
- `Tag`: Tree identifier
- `OriginalStemID`: True stem ID within tree
- `TrueStemID`: Observed stem ID (masked in early censuses)
- `CensusID`: Census number (1-9)
- `ExactDate`: Date of census (YYYY-MM-DD format)
- `DBH`: Observed DBH (cm, with measurement error)
- `CensusInterval`: Years between this census and the previous one (variable, ~5 years ± noise)

## Biological Models

### Growth Model
Size-dependent annual growth follows:
```
μ(DBH) = α + γ × log(DBH)
```
Where:
- `α = 0.4` cm/year (intercept)
- `γ = 0.2` (slope for log(DBH))

Growth variability:
```
σ(DBH) = σ₀ + σ₁ × DBH
```
Where:
- `σ₀ = 0.1` cm/year (intercept)
- `σ₁ = 0.01` (slope for DBH)

### Mortality Model
Size-dependent hazard rate:
```
hazard(DBH) = h₀ × exp(β × DBH)
P(death) = 1 - exp(-hazard × interval_years)
```
Where:
- `h₀ = 0.004` (baseline hazard)
- `β = 0.02` (DBH effect)

### Recruitment Model
New stems enter with lognormal DBH distribution:
```
DBH ~ Lognormal(μ = log(2), σ = 0.8)
```
Recruitment probability per available stem slot: 0.3 per census

### Measurement Error Model

References: Chave, J., Condit, R., Aguilar, S., Hernandez, A., Lao, S., & Perez, R. (2004). Error propagation and scaling for tropical forest biomass estimates. *Philosophical Transactions of the Royal Society of London. Series B: Biological Sciences*, 359(1443), 409-420.

Mixture model matching DP workflow:
- **Small errors** (95% probability): `SD = 0.0062 × DBH + 0.0904` cm
- **Large errors/blunders** (5% probability): `SD = 4.64` cm

## Simulation Process

The `simulate_data.R` script is organized into the following sections:

1. **Setup**: Load required libraries and set random seed
2. **Simulation Parameters**: Define all biological and methodological parameters
3. **Helper Functions**: Utility functions for parameter scaling and truncated distributions
4. **Simulation Engine**:
   - **Species Configuration**: Generate species with trait scaling factors
   - **Individual Stem Trajectory Simulation**: Simulate growth, mortality, and recruitment for each stem
   - **Tree-Level Simulation**: Generate multi-stem trees with variable numbers of stems
   - **Species-Level Simulation**: Create populations of trees for each species
   - **Full Dataset Assembly**: Combine all simulated data into final dataset
5. **Hardcoded Diagnostic Tags**: Append BCI-derived multi-stem patterns, M-code test tags, and row-count invariant edge cases
6. **Data Processing and Export**: Clean data, convert dates, and write output files
7. **Diagnostic Plots**: Generate trajectory visualizations for validation

## Parameters

### Simulation Structure
- `n_census`: Number of census intervals (9)
- `census_interval_years`: Base years between censuses (5)

### Species Configuration
- `n_species`: Number of species (3)
- `n_trees_per_species`: Trees per species [10, 15, 13] (38 total trees)
- `max_stems`: Maximum stems per tree (6)
- `scale_range`: Trait scaling range [1.0, 1.7]
- Species get evenly distributed scaling factors

### Recruitment Process
- `recruit_prob`: Probability of recruitment per available stem slot (0.3)
- `threshold_dbh`: DBH threshold for stems to become observable (1 cm)
- `meanlog`: Lognormal mean for recruit DBH (log(2) ≈ 0.69)
- `sdlog`: Lognormal SD for recruit DBH (0.8)

### Growth Process
- `alpha`: Growth model intercept (0.4 cm/year)
- `gamma`: Growth model slope for log(DBH) (0.2)
- `sigma0`: Growth variability intercept (0.1 cm/year)
- `sigma1`: Growth variability slope for DBH (0.01)
- `min_annual_growth`: Minimum allowed annual growth (0 cm/year)
- `max_annual_growth`: Maximum allowed annual growth (7.5 cm/year)

### Initial Size Distribution
- `census1_meanlog`: Lognormal mean for established stems (log(12) ≈ 2.48)
- `census1_sdlog`: Lognormal SD for established stems (0.9)
- `recruit_meanlog`: Lognormal mean for recruits (log(0.7) ≈ -0.36)
- `recruit_sdlog`: Lognormal SD for recruits (0.4)
- `min_dbh_true`: Biological minimum DBH (0.1 cm)
- `max_dbh_true`: Biological maximum DBH (210 cm)

### Measurement Error Model
- `use_measurement_error`: Enable/disable measurement error (TRUE)
- `meas_sd1_a`: Diameter-dependent error slope (0.0062)
- `meas_sd1_b`: Diameter-dependent error intercept (0.0904)
- `meas_sd2`: Large error SD (4.64 cm)
- `meas_p_big`: Probability of large errors (0.05)
- `min_dbh_obs`: Minimum observed DBH after error (1.0 cm)

### Mortality Process
- `h0`: Baseline hazard rate (0.004)
- `beta`: DBH effect on hazard (0.02)

### Growth Scaling Events
Simulates disturbances by applying multipliers. Multiple events per species-census combination are multiplied together for compounding effects:
```r
events = list(
    list(species = "sp1", census = 2, multiplier = 1.7),  # sp1 enhanced growth in census 2
    list(species = "sp1", census = 5, multiplier = 1.7),  # sp1 enhanced growth in census 5
    list(species = "sp1", census = 8, multiplier = 1.7),  # sp1 enhanced growth in census 8
    list(species = "sp3", census = 2, multiplier = 0.1),  # sp3 suppressed growth in census 2
    list(species = "sp3", census = 5, multiplier = 0.1),  # sp3 suppressed growth in census 5
    list(species = "sp3", census = 8, multiplier = 0.1)   # sp3 suppressed growth in census 8
)
```

### Stem ID Masking
- `anchor_start_census`: Censuses from this point have trusted IDs (7)
- Earlier censuses have `TrueStemID = NA`

### M-Code Test Tags (901–903)

Three additional tags (901, 902, 903) are appended to the simulated dataset to validate the M-coded main-stem constraint in the DP solver:

- **Tag 901**: 1 stem at C1–C2, branches to 2 stems at C3. Stem with M code (5.2 cm) is nearly the same size as the other stem (5.1 cm). Without M-pinning, the assignment is ambiguous; with M-pinning, the M-coded stem must trace back to the single pre-branching stem.
- **Tag 902**: M code present only in the anchor census (legacy scenario). Tests that M codes at the anchor do not interfere with the DP solver.
- **Tag 903**: M code applied to the smaller of two stems (2.0 cm vs 8.0 cm). Tests that M-pinning overrides growth-based preference (which would favor matching the larger stem).

All three tags include a `ListOfTSM` column carrying the `M` code. The `OriginalStemID` column provides ground truth for validation but is not used by the DP algorithm.

### Row-Count Invariant Edge Cases (Tags 9901–9945)

45 additional tags (EC1–EC45, Tags 9901–9945) are appended to exercise every combination of census span, stem count, DBH availability, and special flags. These tags verify that the DP workflow preserves a strict row-count invariant: every input row must appear exactly once in the output (no duplication, no loss).

The edge cases are grouped by scenario type:

**Single-census / minimal tags:**
- **EC1** (9901): Post-anchor only, single census C8
- **EC4** (9904): Single census at anchor C7
- **EC6** (9906): Single row, all fields NA except Species/Tag
- **EC18** (9918): Single census C1 only (earliest)
- **EC19** (9919): Single census C9 only (latest)
- **EC40** (9940): Anchor census only, DBH NA (dead anchor → skip)

**Post-anchor only tags:**
- **EC2** (9902): Two censuses C8–C9, single stem
- **EC3** (9903): Multi-stem (2 stems) at C8–C9
- **EC7** (9907): Anchor + post-anchor C7–C9
- **EC11** (9911): Multi-stem, post-anchor only, mixed DBH/NA
- **EC22** (9922): 3 stems at C8 (high K, tests state-space sizing)
- **EC26** (9926): Post-anchor only, all DBH NA → `skipped_no_data`

**Pre-anchor only tags:**
- **EC9** (9909): Pre-anchor single census C3
- **EC14** (9914): Pre-anchor C1–C5 only, no anchor
- **EC28** (9928): Pre-anchor multi-census multi-stem (2 stems at C1, C3, C5)
- **EC39** (9939): Pre-anchor only, all DBH NA → `skipped_no_data`

**Full-span and cross-boundary tags:**
- **EC13** (9913): Full span C1–C9, single stem, all DBH valid (happy-path baseline)
- **EC15** (9915): Pre-anchor + anchor, no post-anchor, multi-stem
- **EC20** (9920): C1 and C9 only (maximum gap spanning pre+post anchor)
- **EC21** (9921): All 9 censuses, all DBH NA except anchor C7
- **EC23** (9923): Palm species, full span (different `growth_form` pathway)
- **EC30** (9930): Sparse full span — DBH only at C1, C5, C9

**Anchor extension tags:**
- **EC10** (9910): Pre-anchor + anchor + post-anchor, DBH only at post-anchor (forces anchor extension)
- **EC16** (9916): Anchor has NA DBH, pre- and post-anchor have DBH
- **EC29** (9929): Anchor C7 NA DBH + single post-anchor C8 with DBH
- **EC36** (9936): Anchor C7 NA + two post-anchor censuses
- **EC37** (9937): Multi-stem anchor extension (2 stems, anchor all NA)
- **EC43** (9943): Anchor extension where first post-anchor has no TrueStemID but second does
- **EC44** (9944): Deepest anchor extension — skip C8 (also NA), extend to C9

**Mortality and recruitment tags:**
- **EC17** (9917): Stem 1 dies at C5, stem 2 recruited at C5
- **EC27** (9927): 3 stems at anchor, 2 die post-anchor
- **EC45** (9945): 1 stem grows, 2nd stem recruits at C9

**Growth edge cases:**
- **EC24** (9924): Zero growth — identical DBH across 5 censuses
- **EC25** (9925): Apparent shrinkage — DBH decreases between censuses
- **EC31** (9931): Very large DBH >150 cm (outlier/bounds test)
- **EC34** (9934): Two stems, identical DBH at every census (maximally confusable)

**Special flags:**
- **EC5** (9905): All DBH NA but valid CensusIDs across full range
- **EC8** (9908): Multi-stem at anchor only (2 stems, C7 only)
- **EC12** (9912): R-flag (resprout) at post-anchor
- **EC32** (9932): Multiple R-flag rows across censuses
- **EC33** (9933): ListOfTSM = "B" (broken stem)

**Census pattern tags:**
- **EC35** (9935): Pre-anchor with gaps (C1, C3, C5, C7 — skipping even censuses)
- **EC38** (9938): 3 stems at same census C5, then 3 at anchor C7
- **EC41** (9941): C6 + C7 only (two consecutive censuses ending at anchor)
- **EC42** (9942): Dense post-anchor with 2 stems at C7, C8, C9

### Visualization
- `make_plot`: Generate trajectory plots (TRUE)

## Output Files

### Main Dataset: `simulated_data_1.csv`
This CSV contains the longitudinal stem-level measurements simulated by `simulate_data.R`. Each row is a single (Tag, OriginalStemID, CensusID) observation. Below are the columns, types, units, and NA semantics to help you load and use the dataset:

| Column | Type | Units / format | Description & NA semantics |
|--------|------|----------------|----------------------------|
| `Species` | character | N/A | Species code/name for the tree. |
| `Tag` | integer / character | N/A | Unique tree identifier (grouping for stems). Use together with `OriginalStemID` to refer to a physical stem. |
| `OriginalStemID` | integer | N/A | The true stem index within the tree (1,2,3,...). Identifies a specific stem across censuses. |
| `TrueStemID` | integer or NA | N/A | Observed/masked stem ID used to emulate field-recorded IDs. For censuses before `anchor_start_census` this column is intentionally set to `NA` to simulate unreliable early IDs; values become available (non-NA) at and after the anchor census. |
| `CensusID` | integer | N/A | Census ordinal (1 = first census, 2 = second, ...). Used for ordering and grouping observations. |
| `CensusInterval` | numeric or NA | years | Interval (years) between this census and the previous census. `NA` for the first census or if interval is unknown. Use `interval_years = NULL` in downstream DP functions to prefer per-row intervals when present. |
| `DBH` | numeric or NA | cm | Observed DBH value after applying measurement error. Missing (`NA`) means unobserved (stem not recorded in that census). |
| `ExactDate` | character (YYYY-MM-DD) | date | Exact calendar date of the census for this row. Useful to compute intervals directly if `CensusInterval` is missing or to verify time continuity. |

Notes & usage tips:
- Rows with `DBH = NA` indicate the stem was not observed at that census (true absence or below detection threshold); downstream routines should handle missing DBH appropriately.
- `TrueStemID` is the field-style ID used to emulate realistic matching challenges. The simulation sets `TrueStemID = NA` before the configured `anchor_start_census` to simulate masked early IDs.
- `CensusInterval` is provided per-row but you can compute intervals from `ExactDate` if preferred (e.g., for irregular sampling). When running the DP/marginal inference, set `interval_years = NULL` to allow functions to detect per-pair intervals automatically from the data (see `dp_global/R/dp_global_utils.R::resolve_interval_years`).

Quick R example (validate file & anchor census):

```r
library(data.table)
dt <- fread("data/simulated_data_1.csv")
# basic summary
str(dt)
# identify anchor census(es) (where TrueStemID becomes available)
sort(unique(dt[!is.na(TrueStemID), CensusID]))
# check first few rows for a Tag
dt[Tag == 1L, .SD[1:12]]
```

This detailed schema should help you connect simulation outputs to the DP workflow and tests.

### Diagnostic Plots
- `simulated_data_1.pdf`: Species-level growth trajectories (one page per species)
- `simulated_data_tag_level_trajectories_1.pdf`: Tag-level individual stem trajectories (one page per tree)

### File Naming
Output files use simple sequential numbering (e.g., `simulated_data_1.csv`, `simulated_data_2.csv`) for easy reference. The simulation parameters are embedded in the script and can be inspected directly.

## Usage

### Running the Simulation
```bash
cd /path/to/data_simulation
Rscript simulate_data.R
```

The script generates synthetic forest census data with:
- Variable census intervals (~5 years)
- Date-stamped censuses starting from 1980-01-01
- Multi-species dynamics with realistic growth, mortality, and recruitment
- Measurement error and stem identification challenges

### Modifying Parameters
Edit the `params` list in `simulate_data.R` to customize:
- Number of species and trees per species
- Biological parameters (growth rates, mortality, recruitment)
- Census timing and intervals
- Measurement error settings
- Growth scaling events for disturbance simulation

### Dependencies
- R packages: `data.table`, `ggplot2`, `here`
- Compatible with DP workflow in `../dp_global/R/dp_global_bio.R`

## Applications

This simulated data is designed for:
- **Stem identification algorithm testing**: Validate reconstruction methods
- **Measurement error impact assessment**: Study effects of observation noise
- **Growth scaling scenario analysis**: Test disturbance response algorithms
- **Multi-stem dynamics**: Evaluate handling of trees with multiple stems

## Biological Realism

- The simulation incorporates key tropical forest dynamics and ensures comprehensive edge case coverage for algorithm testing:
  - **Size-dependent processes**: Growth and mortality scale with tree size
  - **Stochastic variability**: Process noise in growth and recruitment
  - **Species differences**: Trait variation creates diverse responses
  - **Measurement challenges**: Realistic observation errors and ID uncertainties
  - **Successional dynamics**: Recruitment balances mortality over time
  - **Single-stem edge cases**: Every species includes both a persistent and an early-dead single-stem tree for robust testing
  - **Row-count invariant edge cases**: 45 additional tags (EC1–EC45) systematically cover every combination of census span, anchor position, DBH availability, stem count, and special flags (R-code, B-code, M-code)

## File Organization

```
data_simulation/
├── simulate_data.R          # Main simulation script (organized in sections)
├── README.md                # This documentation
└── data/                    # Output directory
    ├── simulated_data_1.csv                            # Main dataset (136 tags)
    ├── report_run_simulated_data_1.csv                 # DP run report
    ├── stem_reconstruction_dp_global_rcpp.csv          # DP reconstruction output
    ├── simulated_data_1.pdf                            # Species-level trajectory plots (not tracked by git)
    └── simulated_data_tag_level_trajectories_1.pdf     # Tag-level trajectory plots (not tracked by git)
```

## Notes

- Random seed (1234) ensures reproducible results
- Census intervals are 5 years (fixed base interval)
- Census dates start from 1980-01-01
- All parameters match DP workflow conventions
- Measurement error model validated against field data
- Growth scaling enables controlled disturbance experiments
- Stem ID masking creates realistic identification challenges
- Code is organized into clear sections for maintainability