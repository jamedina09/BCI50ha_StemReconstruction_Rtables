# Forest Census Data Simulation for Stem Identification Testing

This directory contains the `simulate_data.R` script that generates realistic forest census data for testing stem identification and reconstruction algorithms.

## Overview

The simulation creates biologically plausible multi-species tropical forest dynamics using models derived from the DP (Dynamic Programming) workflow in `dp_global_biol.R`. It produces datasets that mimic ForestGEO census protocols with realistic growth, mortality, recruitment, measurement error, and stem identification challenges.

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
5. **Data Processing and Export**: Clean data, convert dates, and write output files
6. **Diagnostic Plots**: Generate trajectory visualizations for validation

## Parameters

### Simulation Structure
- `n_census`: Number of census intervals (9)
- `census_interval_years`: Base years between censuses (5), with random noise added: `interval_years = 5 + rnorm(n_census-1, 0, 0.1)`

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

### Main Dataset: `simulated_data_1.csv`
Columns:
- `Species`: Species name
- `Tag`: Tree identifier (unique tree ID)
- `OriginalStemID`: True stem ID within tree (1, 2, 3, etc.)
- `TrueStemID`: Observed stem ID (NA for censuses before anchor census to simulate unreliable early IDs)
- `CensusID`: Census number (1-9)
- `ExactDate`: Date of census (YYYY-MM-DD format, starting from 1980-01-01)
- `DBH`: Observed DBH (cm, with measurement error applied)
- `CensusInterval`: Years between this census and the previous one (variable around 5 years)

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
- Variable census intervals (~5 years ± random noise)
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
├── simulate_data.R          # Main simulation script (organized in sections)
├── README.md               # This documentation
├── data/                   # Output directory
│   ├── simulated_data_1.csv                    # Main dataset
│   ├── simulated_data_1.pdf                    # Species-level trajectory plots
│   └── simulated_data_tag_level_trajectories_1.pdf  # Tag-level trajectory plots
└── [other files]
```

## Notes

- Random seed (1234) ensures reproducible results
- Census intervals are variable: base 5 years + random noise (N(0, 0.1))
- Census dates start from 1980-01-01 and accumulate with variable intervals
- All parameters match DP workflow conventions
- Measurement error model validated against field data
- Growth scaling enables controlled disturbance experiments
- Stem ID masking creates realistic identification challenges
- Code is organized into clear sections for maintainability


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