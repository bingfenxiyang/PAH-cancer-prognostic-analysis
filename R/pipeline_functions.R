assert_columns <- function(data, required, label) {
  missing <- setdiff(required, names(data))
  if (length(missing)) stop(label, " is missing columns: ", paste(missing, collapse = ", "))
}

read_patient_table <- function(path, label) {
  if (!file.exists(path)) stop("Missing ", label, ": ", path)
  x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  assert_columns(x, "sample_id", label)
  x$sample_id <- gsub("\\.", "-", trimws(x$sample_id))
  if (anyNA(x$sample_id) || anyDuplicated(x$sample_id)) stop(label, " contains missing or duplicated sample_id values")
  x
}

prepare_survival_data <- function(expression_path, clinical_path, candidate_genes) {
  expression <- read_patient_table(expression_path, "expression table")
  clinical <- read_patient_table(clinical_path, "clinical table")
  assert_columns(clinical, c("OS_time", "OS_status"), "clinical table")
  clinical$OS_time <- suppressWarnings(as.numeric(clinical$OS_time))
  clinical$OS_status <- suppressWarnings(as.numeric(clinical$OS_status))
  if (any(!clinical$OS_status %in% c(0, 1), na.rm = TRUE)) stop("OS_status must be coded 0/1")
  common <- intersect(expression$sample_id, clinical$sample_id)
  if (length(common) < 10L) stop("Fewer than 10 matched patients; verify sample identifiers")
  expression <- expression[match(common, expression$sample_id), , drop = FALSE]
  clinical <- clinical[match(common, clinical$sample_id), , drop = FALSE]
  if (!identical(expression$sample_id, clinical$sample_id)) stop("Sample-order matching failed")
  available <- intersect(candidate_genes, names(expression))
  missing <- setdiff(candidate_genes, available)
  if (length(missing)) warning("Candidate genes absent from expression data: ", paste(missing, collapse = ", "))
  if (!length(available)) stop("No candidate genes are available")
  out <- cbind(clinical[, c("sample_id", "OS_time", "OS_status")], expression[, available, drop = FALSE])
  out <- out[complete.cases(out[, c("OS_time", "OS_status")]) & out$OS_time > 0, , drop = FALSE]
  attr(out, "missing_candidate_genes") <- missing
  out
}

run_univariable_cox <- function(data, genes) {
  rows <- lapply(genes, function(gene) {
    model <- survival::coxph(stats::reformulate(gene, response = "survival::Surv(OS_time, OS_status)"), data = data)
    s <- summary(model)
    data.frame(gene = gene, coefficient = unname(stats::coef(model)), hazard_ratio = s$conf.int[1, "exp(coef)"],
               ci_lower_95 = s$conf.int[1, "lower .95"], ci_upper_95 = s$conf.int[1, "upper .95"],
               standard_error = s$coefficients[1, "se(coef)"], z_statistic = s$coefficients[1, "z"],
               raw_p_value = s$coefficients[1, "Pr(>|z|)"], stringsAsFactors = FALSE)
  })
  result <- do.call(rbind, rows)
  result$fdr_bh <- stats::p.adjust(result$raw_p_value, method = "BH")
  result
}

run_lasso_cox <- function(data, genes, seed = 123L) {
  set.seed(seed)
  x <- as.matrix(data[, genes, drop = FALSE])
  y <- survival::Surv(data$OS_time, data$OS_status)
  cv <- glmnet::cv.glmnet(x, y, family = "cox", nfolds = 10, type.measure = "deviance")
  fit <- glmnet::glmnet(x, y, family = "cox")
  beta <- as.matrix(stats::coef(cv, s = "lambda.min"))[, 1]
  selected <- names(beta)[beta != 0]
  list(cv = cv, fit = fit, selected_genes = selected,
       coefficients = data.frame(gene = selected, coefficient = unname(beta[selected])), seed = seed)
}

run_rsf <- function(data, genes, seed = 123L) {
  set.seed(seed)
  form <- stats::reformulate(genes, response = "survival::Surv(OS_time, OS_status)")
  fit <- randomForestSRC::rfsrc(form, data = data[, c("OS_time", "OS_status", genes)], ntree = 1000,
                               nodesize = 14, splitrule = "logrank", importance = TRUE, forest = TRUE)
  selection <- randomForestSRC::var.select(fit)
  top <- selection$topvars
  importance <- data.frame(gene = names(fit$importance), variable_importance = unname(fit$importance))
  list(fit = fit, selected_genes = top, importance = importance, seed = seed)
}

fit_signature <- function(data, genes) {
  if (!length(genes)) stop("LASSO-RSF intersection is empty")
  model <- survival::coxph(stats::reformulate(genes, response = "survival::Surv(OS_time, OS_status)"), data = data, x = TRUE)
  s <- summary(model)
  coefficients <- data.frame(gene = names(stats::coef(model)), coefficient = unname(stats::coef(model)),
    hazard_ratio = s$conf.int[, "exp(coef)"], standard_error = s$coefficients[, "se(coef)"],
    z_statistic = s$coefficients[, "z"], raw_p_value = s$coefficients[, "Pr(>|z|)"], row.names = NULL)
  linear_predictor <- as.numeric(stats::predict(model, newdata = data, type = "lp"))
  cutoff <- stats::median(linear_predictor, na.rm = TRUE)
  list(model = model, coefficients = coefficients, training_cutoff = cutoff,
       scored = transform(data, linear_predictor = linear_predictor, relative_risk = exp(linear_predictor),
                          risk_group = factor(ifelse(linear_predictor > cutoff, "High", "Low"), levels = c("High", "Low"))))
}

score_external_cohort <- function(data, signature) {
  required <- signature$coefficients$gene
  assert_columns(data, required, "validation cohort")
  lp <- as.numeric(stats::predict(signature$model, newdata = data, type = "lp", reference = "sample"))
  transform(data, linear_predictor = lp, relative_risk = exp(lp),
            risk_group = factor(ifelse(lp > signature$training_cutoff, "High", "Low"), levels = c("High", "Low")))
}

write_model_formula <- function(coefficients, path) {
  terms <- sprintf("(%0.8f x %s)", coefficients$coefficient, coefficients$gene)
  writeLines(c(paste("LP =", paste(terms, collapse = " + ")), "Relative risk = exp(LP)"), path)
}

