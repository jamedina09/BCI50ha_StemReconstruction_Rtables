############################################################
### test_complexity_estimator.R
### Demonstrate usage of the DP complexity estimation functions
############################################################
rm(list = ls())

# Source the functions
library(here)

source(here("dp_global", "R", "estimate_dp_complexity_function.R"))

# Filenames
list.files(here("data_simulation", "data"), pattern = "^simulated.*\\.csv$")

filename <- "simulated_data_merr_3sp_inc1_c2c5c8_p1p7_dec1_c2c5c8_p0p1.csv"
data_path <- here("data_simulation", "data", filename)

cat("Estimating DP computational complexity for all tags...\n\n")

# Get complexity estimates for all tags
complexity <- estimate_dp_complexity(data_path)

cat("Tags ranked by estimated computational complexity:\n")
print(complexity)

cat("\nTop 5 slowest tags:\n")
print(head(complexity, 5))

# cat("\nDetailed analysis for tag (the slowest):\n")
# details_tag <- get_tag_complexity_details(data_path, tag = 15)
# print(details_tag$observations_per_census)
# cat("\nStates per census for the tag:\n")
# print(details_tag$states_per_census)

# cat("\nTag summary:\n")
# cat("  Max observations per census:", details_tag$max_observations, "\n")
# cat("  Number of tracks (K):", details_tag$K_tracks, "\n")
# cat("  Max states per census:", details_tag$max_states_per_census, "\n")
# cat("  Total transition computations:", format(details_tag$transition_computations, big.mark = ","), "\n")
################################################################################
### EXPORT RESULTS
################################################################################

# Export complexity results table
output_filename <- paste0("report_run_", filename)
output_path <- here("data_simulation", "data", output_filename)

# Ensure data.table is loaded for fwrite
if (!requireNamespace("data.table", quietly = TRUE)) {
    install.packages("data.table")
}
library(data.table)

# Export the complexity table
fwrite(complexity, output_path)

cat(sprintf("\nComplexity results exported to: %s\n", output_path))
