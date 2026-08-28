# Script: 00_install_packages.R
# Purpose: Install the packages required by all analysis scripts.
# Input: None. Output: Installed R packages. Random seed: not applicable.
cran <- c("here", "survival", "survminer", "glmnet", "randomForestSRC", "timeROC", "ggplot2", "ggvenn", "pheatmap", "patchwork", "ggpubr", "reshape2", "svglite", "ragg")
missing <- cran[!vapply(cran, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")

