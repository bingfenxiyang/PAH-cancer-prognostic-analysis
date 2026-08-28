# Script: 01_data_preprocessing.R
# Purpose: Validate and match patient-level expression and survival tables.
# Input/output: See data/README.md. Random seed: not applicable.
# Dependency: R/pipeline_functions.R.
source(here::here("R", "pipeline_functions.R"))
message("Data are validated within each cohort script to prevent use of stale workspace objects.")

