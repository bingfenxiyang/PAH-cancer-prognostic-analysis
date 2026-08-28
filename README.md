# Reproducible pan-cancer prognostic analysis of PAH-associated candidate genes

## Study overview

This repository curates the supplied draft code for univariable Cox screening, LASSO-Cox selection, random survival forest (RSF) selection, multivariable Cox signature construction, Kaplan-Meier analysis, time-dependent ROC analysis, risk visualization, and risk-group expression comparisons. It is a reproducibility scaffold, not evidence of full computational reproduction: the original patient-level inputs, cohort-specific selected genes, final coefficients, and a runnable R environment were not supplied.

## Repository structure

- `original_code/`: English provenance record and checksum for the private source draft retained offline by the authors.
- `R/`: reusable validation and modeling functions.
- `scripts/`: numbered analysis scripts and figure export code.
- `config/`: cohort-to-figure map and candidate genes.
- `data/`: input schema and provenance fields; raw data are excluded.
- `results/`: expected cohort-specific outputs.
- `environment/`: version-capture placeholders to replace after execution.
- `reports/`: code audit, reproduction status, and publication statements.

## System requirements and installation

Use an R release and package versions captured from the actual reproduction environment. Run `Rscript scripts/00_install_packages.R`, then replace `environment/sessionInfo.txt` with `capture.output(sessionInfo())` and record package versions. No versions have been invented in this repository.

## Input format

See `data/README.md`. Expression and clinical tables are matched explicitly by unique `sample_id`; row order is never assumed. Overall survival must be supplied in days as `OS_time`, with `OS_status` coded 0 (censored/alive) or 1 (death).

## Reproduction workflow

1. Complete the public-data provenance table in `data/README.md`.
2. Place processed expression and clinical tables at the paths in `config/cohorts.csv`.
3. Verify that training and validation expression units and preprocessing are compatible.
4. Run `scripts/03_LGG_Figure7.R` through `scripts/07_LAML_SupplementaryFigure3.R`.
5. Source `scripts/09_generate_prognostic_panels.R` and call `plot_cohort_panels()` for each completed cohort.
6. Run `scripts/08_export_model_coefficients.R`.
7. Capture `sessionInfo()`, package versions, checksums, console logs, and any warnings.
8. Compare all regenerated sample counts, selected genes, coefficients, hazard ratios, P values, AUCs, and plots with the submitted manuscript.

## Figure-to-script mapping

| Manuscript result | Cancer type | Script | Main outputs |
|---|---|---|---|
| Figure 7 | LGG | `03_LGG_Figure7.R` | Cox, LASSO, RSF, KM, ROC, risk and expression plots |
| Figure 8 | LUAD | `04_LUAD_Figure8.R` | Corresponding LUAD analyses |
| Figure 9 | LIHC | `05_LIHC_Figure9.R` | Corresponding LIHC analyses |
| Supplementary Figure 2 | PAAD | `06_PAAD_SupplementaryFigure2.R` | Corresponding PAAD analyses |
| Supplementary Figure 3 | LAML | `07_LAML_SupplementaryFigure3.R` | Corresponding LAML analyses |

Within Figure 7, panel A is univariable Cox regression; B-C are LASSO-Cox cross-validation and coefficient paths; D is RSF error/importance; E is the LASSO-RSF intersection; F is training-cohort KM analysis; G is time-dependent ROC; H is validation-cohort KM analysis; I is the risk/status/heat-map composite; and J shows signature-gene expression by risk group. This mapping is supported by the supplied image and draft code but must be checked against the final legends for Figures 8-9 and supplementary figures.

## Model development and validation

The training workflow uses 10-fold LASSO-Cox cross-validation and `lambda.min`. RSF uses `randomForestSRC`, `ntree=1000`, `nodesize=14`, `splitrule="logrank"`, and variable importance. The final multivariable Cox signature uses the intersection of LASSO- and RSF-selected genes. The linear predictor is `LP = sum(beta_i * expression_i)`; relative risk is `exp(LP)`. Training groups are split at the training median LP. External cohorts must use the frozen training genes, coefficients, and training cutoff. Re-fitting a Cox model or deriving a new median in validation data is model redevelopment, not direct external validation.

## Statistical methods and known limitations

Univariable screening retains the supplied raw-P threshold of 0.05 for fidelity, while also exporting Benjamini-Hochberg values. This two-stage, outcome-guided procedure has optimism and multiple-testing risk. The same data are used for screening, feature selection, coefficient estimation, and apparent performance unless an independent validation cohort is supplied. Ten-fold CV tunes lambda but does not by itself provide an unbiased estimate of the complete pipeline's performance. No calibration, C-index with uncertainty, proportional-hazards diagnostics, bootstrap optimism correction, competing-risk assessment, or missing-expression strategy was present in the supplied draft.

## Expected outputs

Each cohort produces complete univariable Cox results, LASSO coefficients, RSF importance, selected genes, final Cox coefficients, an explicit LP formula, scored patients, analysis metadata, and optional external-validation scores. Figure exports use editable SVG/PDF and 600-dpi TIFF where implemented.

## Random seeds

The seed is 123 for LASSO fold assignment and RSF. The seed is recorded in each cohort metadata file. Reproducibility can still depend on R/package versions and parallel settings.

## Data and code availability

Complete the statements in `reports/publication_statements.md` only after the data provenance is verified and this repository has been uploaded. The software repository should remain publicly accessible for at least five years after publication, subject to the final journal policy and institutional retention requirements.

## Citation and contact

Complete `CITATION.cff` and add the corresponding author's institutional email before public release. Do not publish private contact information without author approval.
