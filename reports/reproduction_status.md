# Reproduction status

## Verified from supplied materials

- Figure 7 panel roles are consistent with the visible image and major code blocks.
- The draft specifies 10-fold LASSO-Cox cross-validation with `lambda.min`.
- RSF parameters are `ntree=1000`, `nodesize=14`, and `splitrule="logrank"`.
- The candidate list contains 37 genes.
- Overall-survival status is described as 0=censored/alive and 1=death.

## Not verified

- Patient-level sample matching, exclusions, missing values, duplicate samples, and expression units.
- Cohort-specific selected genes and coefficients for LGG, LUAD, LIHC, PAAD, and LAML.
- Exact correspondence of Figures 8-9 and supplementary panels, because those figures were not supplied.
- Reproduction of any numerical result or plot, because processed inputs and R were unavailable.
- Package versions and operating-system details.

## Current reproduction conclusion

The supplied materials are insufficient to regenerate all target figures from raw or processed input. The curated repository provides deterministic guards and a documented workflow, but full reproduction requires the missing cohort files, final cohort definitions, original preprocessing code, and an R environment.

