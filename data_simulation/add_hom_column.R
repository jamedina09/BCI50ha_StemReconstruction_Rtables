#!/usr/bin/env Rscript
# Add HOM column to simulated data for testing HOM-proportional widening.
library(data.table)
set.seed(42)

dt <- fread(here::here("data_simulation", "data", "simulated_data_1.csv"))

# Add HOM column:
# - trees: mostly 1.3 (standard), some NA
# - palms: mix of NA (20%), 1.3 (30%), and varying HOM from 0.5-3.0 (50%)
dt[, hom := NA_real_]

# Trees: mostly 1.3, some NA
dt[growth_form == "tree", hom := ifelse(runif(.N) < 0.3, NA_real_, 1.3)]

# Palms: mix of NA, standard, and varying HOM
palm_n <- dt[growth_form == "palm", .N]
palm_choice <- sample(c("na", "std", "vary"), palm_n, replace = TRUE, prob = c(0.2, 0.3, 0.5))
dt[growth_form == "palm", hom := fifelse(
    palm_choice == "na", NA_real_,
    fifelse(palm_choice == "std", 1.3,
            round(runif(.N, 0.5, 3.0), 2))
)]

cat("HOM summary for palms:\n")
print(dt[growth_form == "palm", .(
    n = .N,
    n_na = sum(is.na(hom)),
    mean_hom = round(mean(hom, na.rm = TRUE), 3),
    max_dev = round(max(abs(hom - 1.3), na.rm = TRUE), 3)
), by = Tag])

fwrite(dt, here::here("data_simulation", "data", "simulated_data_1_hom.csv"))
cat("\nWrote simulated_data_1_hom.csv with", nrow(dt), "rows\n")
