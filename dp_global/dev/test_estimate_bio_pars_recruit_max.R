# Test recruit_max_source/fixed behavior
library(data.table)
library(here)
# Source the function under test (assume repo root is working dir when invoked)
source("./dp_global/R/dp_global_biol.R")

# Build a dataset with >=5 growth pairs and one recruit
# Six stems with observed DBH in both censuses (growth pairs)
df <- data.table(
  Tag = c(1,1,2,2,3,3,4,4,5,5,6,6,7,7),
  TrueStemID = c(1,1,2,2,3,3,4,4,5,5,6,6,7,7),
  CensusID = c(1,2,1,2,1,2,1,2,1,2,1,2,1,2),
  DBH = c(10,12,8,9,15,16,7,7.5,20,21,4,5,NA,6),
  species = rep("sp", 14)
)

# Basic smoke: should run with an explicit interval
res_fixed <- estimate_bio_pars(df, interval_years = 5, recruit_max_source = "fixed", recruit_max_fixed = 3)
if (!is.finite(res_fixed$recruitment$recruit_max_dbh) || as.numeric(res_fixed$recruitment$recruit_max_dbh) != 3) {
  stop("Fixed recruit_max_dbh did not take effect: ", res_fixed$recruitment$recruit_max_dbh)
}

res_data <- estimate_bio_pars(df, interval_years = 5, recruit_max_source = "data")
if (!is.finite(res_data$recruitment$recruit_max_dbh) || res_data$recruitment$recruit_max_dbh <= 0) {
  stop("Data recruit_max_dbh seems invalid: ", res_data$recruitment$recruit_max_dbh)
}

cat("test_estimate_bio_pars_recruit_max: OK\n")
