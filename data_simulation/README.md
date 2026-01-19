# Forest Census Data Simulation for Stem Identification Testing

This directory contains the `simulate_data.R` script that generates realistic forest census data for testing stem identification and reconstruction algorithms.

## Overview

The simulation creates biologically plausible multi-species tropical forest dynamics using models derived from the DP (Dynamic Programming) workflow in `dp_global_biol.R`. It produces datasets that mimic ForestGEO census protocols with realistic growth, mortality, recruitment, measurement error, and stem identification challenges.

## Key Features

- **Multi-species forests** with species-specific trait variation
- **Multi-stem trees** with recruitment and mortality processes
- **Guaranteed single-stem trees**: For every species, the simulation always includes:
  - One tag (tree) with exactly one stem that survives the entire census period
  - One tag (tree) with exactly one stem that dies before the anchor census
- **Realistic growth trajectories** with size-dependent growth and process variability
- **Measurement error models** matching the DP workflow (diameter-dependent noise + blunders)
- **Flexible growth scaling** for simulating disturbances (drought, etc.)
- **Stem ID masking** to simulate ForestGEO protocols
- **Comprehensive diagnostics** with trajectory visualization

## Output Files


The simulation generates several files with self-documenting filenames. The main dataset always includes:

- All simulated multi-stem trees per species (with random number of stems)
- For each species, one single-stem tree that survives all censuses
- For each species, one single-stem tree that dies before the anchor census

- `simulated_data_[config].csv`: Main dataset with all stem measurements
- `sp_lvl_traj_[config].pdf`: Species-level growth trajectory plots
- `tg_lvl_traj_[config].pdf`: Tag-level individual stem trajectory plots
- `simulation_params_[config].txt`: Text file containing the simulation parameters used

### Filename Configuration

The `[config]` part encodes the simulation settings:

- **Measurement error**: `merr` or `no_merr`
- **Species count**: `3sp` (for 3 species)
- **Growth scaling**: Describes disturbances applied
  - `inc{N}_{censuses}_p{M}`: N species with increased growth (multiplier M) in specified censuses
  - `dec{N}_{censuses}_p{M}`: N species with decreased growth (multiplier M) in specified censuses
  - Censuses are listed as `c2c5c8` (census 2, 5, 8)
  - Multipliers use `p` for decimal point (e.g., `p1p7` = 1.7)

Example: `merr_3sp_inc1_c2c5c8_p1p7_dec1_c2c5c8_p0p1`

This indicates:
- Measurement error enabled
- 3 species total
- 1 species with 1.7x growth increase in censuses 2, 5, 8
- 1 species with 0.1x growth decrease in censuses 2, 5, 8

## Biological Models

### Growth Model
Size-dependent annual growth follows:
```
μ(DBH) = α + γ × log(DBH)
```
Where:
- `α = 0.45` cm/year (intercept)
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
- `h₀ = 0.005` (baseline hazard)
- `β = 0.02` (DBH effect)

### Recruitment Model
New stems enter with lognormal DBH distribution:
```
DBH ~ Lognormal(μ = log(2), σ = 0.8)
```
Recruitment probability per available stem slot: 0.4 per census

### Measurement Error Model

References: Chave, J., Condit, R., Aguilar, S., Hernandez, A., Lao, S., & Perez, R. (2004). Error propagation and scaling for tropical forest biomass estimates. *Philosophical Transactions of the Royal Society of London. Series B: Biological Sciences*, 359(1443), 409-420.

Mixture model matching DP workflow:
- **Small errors** (95% probability): `SD = 0.0062 × DBH + 0.0904` cm
- **Large errors/blunders** (5% probability): `SD = 4.64` cm

## Simulation Process

1. **Species Configuration**: Generate species with trait scaling factors
2. **Parameter Scaling**: Apply species-specific scaling to growth, mortality, and recruitment
3. **Stem Simulation**: For each stem:
   - Determine birth timing (established vs. recruit)
   - Simulate growth trajectory with process variability
   - Apply species-specific growth scaling events
   - Model mortality based on size-dependent hazard
4. **Observation Generation**: Only stems ≥1 cm DBH are measured
5. **Measurement Error**: Add realistic measurement noise to observed DBH
6. **Stem ID Masking**: Simulate ForestGEO protocol (IDs unreliable in early censuses)

## Parameters

### Simulation Structure
- `n_census`: Number of census intervals (9)
- `census_interval_years`: Years between censuses (5). If you need per-pair/per-census intervals, add a per-row column (for example `Bio_IntervalYears`) with the interval in years for that measurement pair — `estimate_bio_pars()` will detect and use per-row intervals when `interval_years = NULL`.

### Species Configuration
- `n_species`: Number of species (2)
- `n_trees_per_species`: Trees per species [20, 15] (35 total trees)
- `max_stems`: Maximum stems per tree (2-7 sampled uniformly)
- `scale_range`: Trait scaling range [1.0, 1.45]
- Species get evenly distributed scaling factors

### Recruitment Process
- `recruit_prob`: Probability of recruitment per available stem slot (0.4)
- `threshold_dbh`: DBH threshold for stems to become observable (1 cm)
- `meanlog`: Lognormal mean for recruit DBH (log(2) ≈ 0.69)
- `sdlog`: Lognormal SD for recruit DBH (0.8)

### Growth Process
- `alpha`: Growth model intercept (0.45 cm/year)
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
- `h0`: Baseline hazard rate (0.005)
- `beta`: DBH effect on hazard (0.02)

### Growth Scaling Events
Simulates disturbances by applying multipliers. Multiple events per species-census combination are multiplied together for compounding effects:
```r
events = list(
    list(species = "sp1", census = 3L, multiplier = 0.5),  # sp1 drought in census 3
    list(species = "sp1", census = 3L, multiplier = 1.2),  # sp1 partial recovery in census 3 (final: 0.5 × 1.2 = 0.6)
    list(species = "all", census = 5L, multiplier = 1.2)   # All species enhanced growth in census 5
)
```

### Stem ID Masking
- `anchor_start_census`: Censuses from this point have trusted IDs (7)
- Earlier censuses have `TrueStemID = NA`

### Visualization
- `make_plot`: Generate trajectory plots (TRUE)

### Visualization
- `make_plot`: Generate trajectory plots (TRUE)

## Output Files

### Main Dataset: `simulated_data_[config].csv`
Columns:
- `Species`: Species name
- `Tag`: Tree identifier
- `OriginalStemID`: True stem ID within tree
- `TrueStemID`: Observed stem ID (masked in early censuses)
- `CensusID`: Census number (1-9)
- `BirthCensus`: Census when stem was born/recruited
- `DeathCensus`: Census when stem died (if applicable)
- `DBH_true`: True DBH (cm, for diagnostics)
- `DBH`: Observed DBH (cm, with measurement error)
- `Bio_IntervalYears`: (optional) Interval in years used for this observation/pair. When present, `estimate_bio_pars()` can use this per-row interval (set `interval_years = NULL`) to compute annualized increments and mortality probabilities.
- `AnnualGrowth`: Annual growth rate (cm/year)
- `YearFactor`: Applied growth scaling multiplier
- `ObsSD`: Measurement error standard deviation

**Special rows:**
- For each species, there will always be:
  - One tag with a single stem that never dies (DeathCensus = n_census)
  - One tag with a single stem that dies before the anchor census (DeathCensus < anchor_start_census)
- `BirthCensus`: Census when stem was born/recruited
- `DeathCensus`: Census when stem died (if applicable)
- `DBH_true`: True DBH (cm, for diagnostics)
- `DBH`: Observed DBH (cm, with measurement error)
- `AnnualGrowth`: Annual growth rate (cm/year)
- `YearFactor`: Applied growth scaling multiplier
- `ObsSD`: Measurement error standard deviation

### Diagnostic Plots
- `sp_lvl_traj_[config].pdf`: Species-level growth trajectories
- `tg_lvl_traj_[config].pdf`: Tag-level individual stem trajectories

### Filename Encoding
Filenames include configuration details for self-documentation:
```
meas_error_3sp_scaling_sp1-c3-0.5_all-c5-1.2
```
- `meas_error`: Measurement error enabled
- `3sp`: 3 species simulated (general format: `{n}sp`)
- `scaling_sp1-c3-0.5`: sp1 reduced to 50% growth in census 3
- `all-c5-1.2`: All species increased to 120% growth in census 5

**Complete naming convention:**
- **Prefix**: `meas_error` (if measurement error enabled) or `no_error`
- **Species count**: `{n}sp` where n = number of species (e.g., `2sp`, `3sp`, `5sp`)
- **Scaling events**: `scaling_{species}-{census}-{multiplier}_...`
  - `sp1-c3-0.5`: species sp1, census 3, 0.5x growth multiplier
  - `all-c5-1.2`: all species, census 5, 1.2x growth multiplier
  - **Multiple events**: All events are listed (e.g., `sp1-c3-0.5_sp1-c6-0.5` for multiple events on sp1)
- **Examples**:
  - `meas_error_2sp_scaling_sp1-c3-0.5_all-c5-1.2.csv`
  - `no_error_1sp_scaling_all-c4-0.8.pdf`
  - `meas_error_3sp_scaling_sp1-c3-0.5_sp1-c6-0.5_sp3-c3-1.2_sp3-c6-1.2.pdf`

## Usage

### Running the Simulation
```bash
cd /path/to/data_simulation
Rscript simulate_data.R
```

### Modifying Parameters
Edit the `params` list in `simulate_data.R` to customize:
- Species composition and traits
- Growth scaling events
- Measurement error settings
- Mortality/recruitment rates

### Dependencies
- R packages: `data.table`, `ggplot2`, `here`
- Compatible with DP workflow in `../dp_global/R/dp_global_biol.R`

## Applications

This simulated data is designed for:
- **Stem identification algorithm testing**: Validate reconstruction methods
- **Measurement error impact assessment**: Study effects of observation noise
- **Growth scaling scenario analysis**: Test disturbance response algorithms
- **Multi-stem dynamics**: Evaluate handling of trees with multiple stems

## Biological Realism

- The simulation incorporates key tropical forest dynamics and ensures edge cases for algorithm testing:
  - **Size-dependent processes**: Growth and mortality scale with tree size
  - **Stochastic variability**: Process noise in growth and recruitment
  - **Species differences**: Trait variation creates diverse responses
  - **Measurement challenges**: Realistic observation errors and ID uncertainties
  - **Successional dynamics**: Recruitment balances mortality over time
  - **Single-stem edge cases**: Every species includes both a persistent and an early-dead single-stem tree for robust testing

## File Organization

```
data_simulation/
├── simulate_data.R          # Main simulation script
├── README.md               # This documentation
├── data/                   # Output directory
│   ├── simulated_data_*.csv    # Main datasets
│   ├── sp_lvl_traj_*.pdf       # Species-level plots
│   └── tg_lvl_traj_*.pdf       # Tag-level plots
└── [other files]
```

## Notes

- Random seed (123) ensures reproducible results
- All parameters match DP workflow conventions
- Measurement error model validated against field data
- Growth scaling enables controlled disturbance experiments
- Stem ID masking creates realistic identification challenges


## Building This Documentation

**Fully offline HTML (recommended):**
```bash
pandoc README.md --standalone --mathml \
  --embed-resources -o README.html
```

**MathJax HTML (requires internet):**
```bash
pandoc README.md --standalone --mathjax \
  --embed-resources -o README.html
```