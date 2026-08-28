# Script: 03_LGG_Figure7.R | Figure 7A-J | LGG prognostic signature
# Inputs/outputs: config/cohorts.csv and results/Figure7_LGG. Seed: 123.
# Packages: here, survival, glmnet, randomForestSRC. Depends on run_cohort_analysis.R.
source(here::here("scripts", "run_cohort_analysis.R")); run_configured_cohort("LGG", 123L)
