# Forest Census Data Simulation for Stem Identification Testing

This directory contains the `simulate_data.R` script, which builds the test dataset for the stem identification engines (`dp_global/`), and `sample_data_BCI/`, which holds the BCI multiple-stem tags table read by `dp_global/scripts/main_cpp_bci.R` and analyses of record patterns in the BCI data.

## Overview

`simulate_data.R` simulates multi-species forest census data (growth, mortality, recruitment and stem ID masking) in the column layout the DP workflow reads, and appends hand-written tags. The final dataset contains 142 tags:

- **42 simulated-style trees** (Tags 1–42): 38 randomly generated trees across 3 species + 4 hand-written "problematic" tags (39–42: tags present only before or only from census 6 on, with one or two stems)
- **52 hardcoded diagnostic tags** (17 tags in 43–88 and 35 BCI tag numbers such as 378, 2747 or 606162): multi-stem patterns copied from BCI data for regression testing, several with `ListOfTSM` codes and `Status` values (dead, broken below, stem dead)
- **3 M-code test tags** (Tags 901–903): tags whose main stem carries the `M` code
- **45 row-count invariant edge case tags** (Tags 9901–9945): combinations of census span, stem count, DBH availability, and special flags

Key features include:
- **Variable census intervals**: about 5 years between censuses, with random noise per tree
- **Date-stamped censuses**: each tree's first census is dated 1980-01-01; later dates follow from its intervals
- **Multi-species forests** with species-specific trait scaling
- **Multi-stem trees** with recruitment and mortality processes
- **Size-dependent growth** with process variability
- **Growth scaling events** for simulating disturbances
- **Stem ID masking**: `TrueStemID` is `NA` before census 7
- **Trajectory plots** per species and per tag

The measurement-error settings in the parameter list (`params$obs`) are **not applied**: the exported `DBH` is the simulated true DBH.

## Quickstart

Run the simulation from the project root (the script builds its paths with `here()`):

```bash
Rscript data_simulation/simulate_data.R
```

Generated outputs are placed in `data_simulation/data/` (the folder must exist). The random seed is fixed, so every run writes the same dataset.

## Output Files

- `simulated_data_1.csv`: Main dataset with all stem records
- `simulated_data_1.pdf`: Species-level growth trajectory plots (one page per species)
- `simulated_data_tag_level_trajectories_1.pdf`: Tag-level stem trajectory plots (one page per tag)

The PDFs are written when `ggplot2` is installed and `params$plot$make_plot` is `TRUE`.

### Main Dataset: `simulated_data_1.csv`

Each row is one record of a stem in a census.

| Column | Type | Units / format | Description & NA semantics |
|--------|------|----------------|----------------------------|
| `Species` | character | N/A | Species code (`sp1`, `sp2`, `sp3`). |
| `Tag` | integer | N/A | Tree identifier (grouping for stems). |
| `OriginalStemID` | integer | N/A | Stem identifier within the tree: the simulated identity for Tags 1–42, a database-style StemID in the hardcoded tags (`NA` on some of their rows). Ground truth for validation; the engines do not use it to link stems. |
| `TrueStemID` | integer or NA | N/A | Stem identity the engines may trust. `NA` before census 7 for the simulated trees; in the hardcoded tags it is set row by row. |
| `CensusID` | integer | N/A | Census number (1–9; `NA` on one edge-case row). |
| `DBH` | numeric or NA | cm | Diameter. `NA` means the stem has no measurement in that census (not yet recruited, dead, or an unmeasured record in the hardcoded tags). |
| `ExactDate` | date (YYYY-MM-DD) | date | Date of the record; census intervals are computed from it. |
| `ListOfTSM` | character | codes | Field codes on the hardcoded tags (e.g. `R`, `OR`, `M`, `B`, `MF`); empty elsewhere. |
| `Status` | character | N/A | Field status on some hardcoded tags (`alive`, `dead`, `stem dead`, `broken below`); empty elsewhere. |
| `growth_form` | character | N/A | `tree` (`sp1`, `sp2`) or `palm` (`sp3`). |

Quick R example (validate file & anchor census):

```r
library(data.table)
dt <- fread("data_simulation/data/simulated_data_1.csv")
str(dt)
# censuses in which TrueStemID is available
sort(unique(dt[!is.na(TrueStemID), CensusID]))
# records of one tag
dt[Tag == 1L]
```

## Biological Models

The values below are the base parameters; each species multiplies the growth coefficients, the baseline hazard and the odds of late recruitment by its scaling factor.

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

The annual growth of an interval is drawn from Normal(μ, σ), multiplied (its mean) by the growth-scaling factor of the species and census, and clipped to `[min_annual_growth, max_annual_growth]`.

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
Each stem is either present at the first census or born later: with probability `recruit_prob` (0.3) it is born in a census drawn at random from 2 to 8. A stem born later starts with

```
DBH ~ Lognormal(μ = log(2), σ = 0.8), truncated to [0.1, 0.99] cm
```

A stem present at the first census starts with `DBH ~ Lognormal(log(12), 0.9)`, truncated to [0.1, 210] cm.

### Measurement Error

No measurement error is added. The parameter list holds the settings of the mixture model used in the DP workflow (small errors with `SD = 0.0062 × DBH + 0.0904` cm, large errors with `SD = 4.64` cm and probability 0.05), but the script does not read them.

Reference for that model: Chave, J., Condit, R., Aguilar, S., Hernandez, A., Lao, S., & Perez, R. (2004). Error propagation and scaling for tropical forest biomass estimates. *Philosophical Transactions of the Royal Society of London. Series B: Biological Sciences*, 359(1443), 409-420.

## Simulation Process

The `simulate_data.R` script is organized into the following sections:

1. **Setup**: Load required libraries and set the random seed (1234)
2. **Simulation Parameters**: The `params` list
3. **Helper Functions**: Species table, growth multipliers, truncated lognormal, per-species parameter scaling
4. **Simulation Engine**:
   - **Species Configuration**: Species with trait scaling factors
   - **Individual Stem Trajectory Simulation**: Birth, growth and mortality of each stem
   - **Tree-Level Simulation**: Trees with 2 to `max_stems` stems sharing the tree's census intervals
   - **Species-Level Simulation**: The trees of each species
   - **Full Dataset Assembly**: All simulated trees, with `TrueStemID` masked before census 7
5. **Data Processing**: Dates from the census intervals
6. **Hardcoded Tags**: The four "problematic" tags, the BCI-derived multi-stem patterns, the M-code tags and the row-count edge cases; growth forms; export of the CSV
7. **Diagnostic Plots**: Trajectory plots per species and per tag

## Parameters

### Simulation Structure
- `n_census`: Number of censuses (9)
- `census_interval_years`: Base years between censuses (5); each tree draws its own intervals as 5 + Normal(0, 0.1) years

### Species Configuration
- `n_species`: Number of species (3)
- `n_trees_per_species`: Trees per species [10, 15, 13] (38 total trees)
- `max_stems`: Maximum stems per tree (6); every simulated tree has 2 to 6 stems
- `scale_range`: Trait scaling range [1.0, 1.7]
- Species get evenly distributed scaling factors

### Recruitment Process
- `recruit_prob`: Probability that a stem is born after the first census (0.3)
- `threshold_dbh`: Upper limit of a recruit's first DBH (1 cm; the first DBH is drawn below 0.99 × this value)
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
- `census1_meanlog`: Lognormal mean for stems present at the first census (log(12) ≈ 2.48)
- `census1_sdlog`: Lognormal SD for stems present at the first census (0.9)
- `recruit_meanlog`, `recruit_sdlog`: defined (log(0.7), 0.4) but not read by the script
- `min_dbh_true`: Biological minimum DBH (0.1 cm)
- `max_dbh_true`: Biological maximum DBH (210 cm)

### Measurement Error Model (defined, not applied)
- `use_measurement_error` (TRUE), `meas_sd1_a` (0.0062), `meas_sd1_b` (0.0904), `meas_sd2` (4.64 cm), `meas_p_big` (0.05), `min_dbh_obs` (1.0 cm): none of these is read by the script

### Mortality Process
- `h0`: Baseline hazard rate (0.004)
- `beta`: DBH effect on hazard (0.02)

### Growth Scaling Events
Simulates disturbances by applying multipliers to the mean growth. Multiple events per species-census combination are multiplied together for compounding effects:
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
- `anchor_start_census`: Censuses from this point have trusted IDs (7). The value is read for the four "problematic" tags; the simulated trees are masked before census 7 directly
- Earlier censuses have `TrueStemID = NA`

### Visualization
- `make_plot`: Generate trajectory plots (TRUE)

## Hardcoded Tags

### M-Code Test Tags (901–903)

Three tags (901, 902, 903) carry the `M` code (multiple stems) in `ListOfTSM`. The `dp_global` engines do not read the `M` code; the tags keep such records in the test data.

- **Tag 901**: 1 stem at C1–C2, branches to 2 stems at C3. The stem with the M code (5.2 cm) is nearly the same size as the other stem (5.1 cm), so growth alone does not tell which one continues the single earlier stem.
- **Tag 902**: Same structure as Tag 901, with the M code kept on the main stem in the censuses after the branching.
- **Tag 903**: M code on the smaller of the two stems (2.0 cm vs 8.0 cm) at the branching census.

The `OriginalStemID` column provides ground truth for validation but is not used by the DP algorithm.

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
- **EC12** (9912): Unmeasured row with the R code (resprout) at post-anchor C8
- **EC32** (9932): Multiple R-flag rows across censuses
- **EC33** (9933): ListOfTSM = "B" (broken stem)

**Census pattern tags:**
- **EC35** (9935): Pre-anchor with gaps (C1, C3, C5, C7 — skipping even censuses)
- **EC38** (9938): 3 stems at same census C5, then 3 at anchor C7
- **EC41** (9941): C6 + C7 only (two consecutive censuses ending at anchor)
- **EC42** (9942): Dense post-anchor with 2 stems at C7, C8, C9

## Usage

### Modifying Parameters
Edit the `params` list in `simulate_data.R` to customize:
- Number of species and trees per species
- Biological parameters (growth rates, mortality, recruitment)
- Census timing and intervals
- Growth scaling events for disturbance simulation

### Dependencies
- R packages: `data.table`, `here`; `ggplot2` for the plots
- The dataset is the default input of `dp_global/scripts/main_cpp.R` and `main_cpp_chunk.R`

## Applications

This simulated data is designed for:
- **Stem identification algorithm testing**: Validate reconstruction methods
- **Growth scaling scenario analysis**: Test how the engines handle disturbed growth
- **Multi-stem dynamics**: Evaluate handling of trees with multiple stems
- **Edge cases**: Check that every input row comes back exactly once (row-count invariant) and that field codes are handled

## File Organization

```
data_simulation/
├── simulate_data.R          # Main simulation script (organized in sections)
├── README.md                # This documentation
├── data/                    # Output directory
│   ├── simulated_data_1.csv                            # Main dataset (142 tags)
│   ├── simulated_data_1.pdf                            # Species-level trajectory plots
│   └── simulated_data_tag_level_trajectories_1.pdf     # Tag-level trajectory plots
└── sample_data_BCI/
    ├── multistem_tags.rds   # BCI multiple-stem tags with species parameters (Bio_* columns);
    │                        # input of dp_global/scripts/main_cpp_bci.R
    └── general_data/        # Not under version control
        ├── _status_pattern_lib.R, broken_below_pattern.R   # Helpers of the pattern reports
        ├── broken_below_pattern.qmd, dead_pattern.qmd, missing_pattern.qmd
        │                    # How broken-below, dead and missing records relate to StemIDs in the BCI data
        └── ForestGEO_codes/ # forestgeo_codes_reference.qmd, improvement_proposals.md
```

## Notes

- Random seed (1234) ensures reproducible results
- Census intervals are 5 years plus a small random deviation per tree
- Each tree's census dates start from 1980-01-01
- Column names and units (DBH in cm) are those the DP workflow reads
- No measurement error is added to the DBH
- Growth scaling enables controlled disturbance experiments
- Stem ID masking creates the identification problem the engines solve
