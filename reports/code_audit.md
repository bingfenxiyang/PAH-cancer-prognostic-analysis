# Code audit report

## Overall assessment

**Not yet reproducible; major scientific and implementation issues require resolution before submission.** The original file is a working notebook-like draft with Chinese prose embedded as invalid R syntax, repeated blocks, undefined objects, positional column selection, and no cohort-specific data provenance.

## P0 — may change results or scientific interpretation

1. **Expression/clinical matching is not implemented.** The draft refers to merged data but does not show an ID-based join; row-order mismatch could invalidate every model. The curated workflow performs one-to-one matching by unique `sample_id` and stops on duplicates.
2. **Validation refits the Cox model and recalculates the median.** This is not direct external validation and can inflate apparent validation performance. The curated workflow applies frozen training coefficients and the training cutoff.
3. **Risk-score definition is conflated.** `predict(..., type="risk")` returns `exp(LP)`, while the manuscript formula is `LP=sum(beta*expression)`. Both are now exported separately.
4. **Invalid/undefined columns and objects.** Examples include `cli`, `data$OS` after transposition, `colnames(dt)[？]`, fixed positional columns, and an undefined palette in one plotting path. The original file cannot run end to end.
5. **Potential leakage/optimism.** Univariable screening, LASSO/RSF selection, coefficient fitting, risk cutoff selection, and apparent performance use the same training observations. Ten-fold CV tunes LASSO only and does not validate the complete pipeline.

## P1 — important statistical/reporting risks

1. Thirty-seven univariable Cox tests use raw P<0.05 without an explicit multiplicity strategy. Raw and BH-adjusted P values are now both exported; changing the selection rule requires author approval because it may change results.
2. No proportional-hazards diagnostics, calibration, C-index uncertainty, bootstrap optimism correction, or complete-case sensitivity analysis is reported.
3. The violin plots run multiple gene-wise Wilcoxon tests without correction. The revised script exports raw and BH-adjusted P values.
4. The code converts survival time to years for ROC but uses days elsewhere. The revised ROC script keeps days and evaluates 365, 1095, and 1825 days.
5. `var.select()` top variables depend on package version and method defaults that were not recorded.

## P2 — code quality and figure reproducibility

- Repeated package loading, duplicated plotting code, debugging prose, Unicode symbols outside comments, calls to `dev.off()` without a guaranteed open device, absolute/implicit paths, and ad hoc output names were removed from executable scripts.
- English names, relative paths, fixed seeds, deterministic output directories, model metadata, coefficients, and formulas were added.
- The original file is retained unchanged under `original_code/`.

## Author input needed before full reproduction

- Exact public database/project identifiers, versions, access dates, preprocessing, normalization, and sample-selection rules for all cohorts.
- Five training and five validation expression/clinical inputs, or the exact subset actually used.
- Confirmation of whether the displayed validation curves were produced by model refitting or frozen-model application.
- Final cohort-specific selected genes, coefficients, risk cutoffs, sample counts, event counts, AUCs, and complete figure legends.
- Whether raw-P screening must remain the primary analysis or whether an FDR/sensitivity analysis can be added.

