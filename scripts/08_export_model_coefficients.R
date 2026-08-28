# Script: 08_export_model_coefficients.R
# Purpose: Combine per-cohort coefficients and risk formulas after cohort scripts finish.
files <- list.files(here::here("results"), pattern = "_model_coefficients\\.csv$", recursive = TRUE, full.names = TRUE)
if (!length(files)) stop("No model coefficient files found; run cohort scripts first")
tables <- lapply(files, function(f) transform(read.csv(f), source_file = basename(f)))
write.csv(do.call(rbind, tables), here::here("results", "all_model_coefficients.csv"), row.names = FALSE)

