# Internal runner used by the five figure-specific scripts.
run_cohort_analysis <- function(cancer, figure_id, training_expression, training_clinical,
                                validation_expression, validation_clinical, output_dir,
                                candidate_genes, seed = 123L, univariable_p_threshold = 0.05) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  source(here::here("R", "pipeline_functions.R"))
  train <- prepare_survival_data(training_expression, training_clinical, candidate_genes)

  # Panel A: Univariable Cox regression. Raw-P filtering reproduces the supplied draft;
  # BH-adjusted P values are exported for transparency and sensitivity assessment.
  univariable <- run_univariable_cox(train, intersect(candidate_genes, names(train)))
  write.csv(univariable, file.path(output_dir, paste0(cancer, "_univariable_cox.csv")), row.names = FALSE)
  screened <- univariable$gene[univariable$raw_p_value < univariable_p_threshold]
  if (length(screened) < 2L) stop("Fewer than two genes passed the prespecified univariable threshold")

  # Panels B-C: 10-fold LASSO-Cox using lambda.min.
  lasso <- run_lasso_cox(train, screened, seed)
  saveRDS(lasso$cv, file.path(output_dir, paste0(cancer, "_lasso_cv.rds")))
  write.csv(lasso$coefficients, file.path(output_dir, paste0(cancer, "_lasso_nonzero_coefficients.csv")), row.names = FALSE)

  # Panel D: Random survival forest (ntree=1000, nodesize=14, log-rank splitting).
  rsf <- run_rsf(train, screened, seed)
  saveRDS(rsf$fit, file.path(output_dir, paste0(cancer, "_rsf_model.rds")))
  write.csv(rsf$importance, file.path(output_dir, paste0(cancer, "_rsf_variable_importance.csv")), row.names = FALSE)

  # Panel E and signature: intersection followed by multivariable Cox regression.
  selected <- intersect(lasso$selected_genes, rsf$selected_genes)
  write.csv(data.frame(gene = selected), file.path(output_dir, paste0(cancer, "_selected_genes.csv")), row.names = FALSE)
  signature <- fit_signature(train, selected)
  write.csv(signature$coefficients, file.path(output_dir, paste0(cancer, "_model_coefficients.csv")), row.names = FALSE)
  write_model_formula(signature$coefficients, file.path(output_dir, paste0(cancer, "_risk_formula.txt")))
  write.csv(signature$scored, file.path(output_dir, paste0(cancer, "_training_scores.csv")), row.names = FALSE)
  saveRDS(signature$model, file.path(output_dir, paste0(cancer, "_multivariable_cox_model.rds")))

  # Validation: apply frozen training genes, coefficients, and training cutoff; do not refit.
  if (file.exists(validation_expression) && file.exists(validation_clinical)) {
    valid <- prepare_survival_data(validation_expression, validation_clinical, selected)
    valid_scored <- score_external_cohort(valid, signature)
    write.csv(valid_scored, file.path(output_dir, paste0(cancer, "_validation_scores.csv")), row.names = FALSE)
  } else {
    warning("Validation files are unavailable; external validation was not run")
  }

  metadata <- list(cancer = cancer, figure = figure_id, seed = seed, nfolds = 10,
                   lambda = "lambda.min", rsf_ntree = 1000, rsf_nodesize = 14,
                   rsf_splitrule = "logrank", univariable_p_threshold = univariable_p_threshold,
                   training_n = nrow(train), selected_genes = selected)
  dput(metadata, file = file.path(output_dir, paste0(cancer, "_analysis_metadata.R")))
  invisible(metadata)
}

run_configured_cohort <- function(cancer_code, seed = 123L) {
  cfg <- read.csv(here::here("config", "cohorts.csv"), stringsAsFactors = FALSE)
  cohort_config <- cfg[cfg$cancer == cancer_code, , drop = FALSE]
  if (nrow(cohort_config) != 1L) stop("Expected exactly one configuration row for ", cancer_code)
  resolve <- function(path) here::here(strsplit(path, "/", fixed = TRUE)[[1]])
  candidate_genes <- scan(here::here("config", "candidate_genes.txt"), what = "character", quiet = TRUE)
  run_cohort_analysis(cancer = cohort_config$cancer, figure_id = cohort_config$figure,
    training_expression = resolve(cohort_config$training_expression), training_clinical = resolve(cohort_config$training_clinical),
    validation_expression = resolve(cohort_config$validation_expression), validation_clinical = resolve(cohort_config$validation_clinical),
    output_dir = resolve(cohort_config$output_dir), candidate_genes = candidate_genes, seed = seed)
}
