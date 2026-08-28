# Script: 09_generate_prognostic_panels.R
# Purpose: Generate publication-ready KM, time-dependent ROC, risk-distribution,
# survival-status, heat-map and violin panels from frozen cohort outputs.
# Corresponding panels: Figure 7F-J and analogous panels in Figures 8-9 and Supplementary Figures 2-3.
# Inputs: *_training_scores.csv, *_model_coefficients.csv. Outputs: SVG/PDF/TIFF panels.
# Packages: survival, survminer, timeROC, ggplot2, pheatmap, patchwork, ggpubr, reshape2, ragg, svglite.
# Random seed: not applicable. Dependency: run the relevant cohort script first.

ggplot2::theme_set(ggplot2::theme_classic(base_size = 7, base_family = "Arial") +
  ggplot2::theme(axis.line = ggplot2::element_line(linewidth = 0.35),
                 axis.ticks = ggplot2::element_line(linewidth = 0.35),
                 panel.grid = ggplot2::element_blank()))

save_plot <- function(plot, stem, width_mm = 89, height_mm = 75, dpi = 600) {
  width <- width_mm / 25.4; height <- height_mm / 25.4
  svglite::svglite(paste0(stem, ".svg"), width = width, height = height); print(plot); grDevices::dev.off()
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = width, height = height, family = "sans"); print(plot); grDevices::dev.off()
  ragg::agg_tiff(paste0(stem, ".tiff"), width = width, height = height, units = "in", res = dpi); print(plot); grDevices::dev.off()
}

plot_cohort_panels <- function(score_file, coefficient_file, output_dir, cancer) {
  scored_data <- read.csv(score_file, check.names = FALSE)
  required <- c("OS_time", "OS_status", "linear_predictor", "risk_group")
  if (!all(required %in% names(scored_data))) stop("Training score file lacks required columns")
  scored_data$risk_group <- factor(scored_data$risk_group, levels = c("High", "Low"))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  pal <- c(High = "#D55E00", Low = "#0072B2")

  # Panel F: Kaplan-Meier curves and risk table; log-rank P value from survdiff.
  km <- survival::survfit(survival::Surv(OS_time, OS_status) ~ risk_group, data = scored_data)
  km_plot <- survminer::ggsurvplot(km, data = scored_data, pval = TRUE, risk.table = TRUE,
    palette = unname(pal), legend.title = "Risk group", legend.labs = c("High", "Low"),
    xlab = "Follow-up time (days)", ylab = "Overall survival probability")
  save_plot(km_plot$plot / km_plot$table + patchwork::plot_layout(heights = c(3, 1)),
            file.path(output_dir, paste0(cancer, "_KM_training")), 120, 100)

  # Panel G: 1-, 3-, and 5-year time-dependent ROC; OS_time remains in days.
  roc <- timeROC::timeROC(T = scored_data$OS_time, delta = scored_data$OS_status, marker = scored_data$linear_predictor,
    cause = 1, weighting = "marginal", times = c(365, 1095, 1825), iid = TRUE)
  tiff_path <- file.path(output_dir, paste0(cancer, "_timeROC.tiff"))
  ragg::agg_tiff(tiff_path, width = 89/25.4, height = 89/25.4, units = "in", res = 600)
  graphics::par(pty = "s", mar = c(4, 4, 1, 1)); plot(roc, time = 365, col = pal["High"], lwd = 1.5)
  plot(roc, time = 1095, col = pal["Low"], add = TRUE, lwd = 1.5); plot(roc, time = 1825, col = "#E69F00", add = TRUE, lwd = 1.5)
  graphics::abline(0, 1, lty = 2); graphics::legend("bottomright", sprintf("%d year AUC = %.2f", c(1,3,5), roc$AUC), col = c(pal, "#E69F00"), lty = 1, bty = "n"); grDevices::dev.off()
  write.csv(data.frame(year = c(1,3,5), auc = roc$AUC), file.path(output_dir, paste0(cancer, "_timeROC_AUC.csv")), row.names = FALSE)

  # Panel I: Risk score, survival status, and row-scaled signature expression.
  scored_data <- scored_data[order(scored_data$linear_predictor), ]; scored_data$patient_index <- seq_len(nrow(scored_data)); scored_data$status <- factor(ifelse(scored_data$OS_status == 1, "Death", "Alive"), levels = c("Death", "Alive"))
  p1 <- ggplot2::ggplot(scored_data, ggplot2::aes(patient_index, linear_predictor, colour = risk_group)) + ggplot2::geom_point(size=.6) + ggplot2::scale_colour_manual(values=pal) + ggplot2::theme_classic(base_size=7) + ggplot2::labs(x=NULL,y="Linear predictor",colour="Risk group")
  p2 <- ggplot2::ggplot(scored_data, ggplot2::aes(patient_index, OS_time, colour = status)) + ggplot2::geom_point(size=.6) + ggplot2::scale_colour_manual(values=c(Death="#D55E00",Alive="#0072B2")) + ggplot2::theme_classic(base_size=7) + ggplot2::labs(x="Patient index",y="Overall survival (days)",colour="Status")
  genes <- read.csv(coefficient_file)$gene; expression_matrix <- as.matrix(scored_data[, genes, drop=FALSE]); scaled_expression <- t(scale(expression_matrix));
  annotation <- data.frame(Risk_group=scored_data$risk_group); rownames(annotation) <- scored_data$sample_id
  pheatmap::pheatmap(scaled_expression, cluster_rows=FALSE, cluster_cols=FALSE, show_colnames=FALSE,
    annotation_col=annotation, annotation_colors=list(Risk_group=pal), filename=file.path(output_dir,paste0(cancer,"_risk_heatmap.tiff")), width=7.2, height=2.4)
  save_plot(p1 / p2, file.path(output_dir,paste0(cancer,"_risk_and_status")), 120, 90)

  # Panel J: Two-sided Wilcoxon rank-sum tests; BH values are exported across signature genes.
  long <- reshape2::melt(scored_data[, c("risk_group", genes), drop=FALSE], id.vars="risk_group", variable.name="gene", value.name="expression")
  tests <- do.call(rbind, lapply(split(long, long$gene), function(d) data.frame(gene=unique(d$gene), raw_p_value=stats::wilcox.test(expression~risk_group,d,exact=FALSE)$p.value)))
  tests$fdr_bh <- stats::p.adjust(tests$raw_p_value, method="BH"); write.csv(tests,file.path(output_dir,paste0(cancer,"_risk_group_expression_tests.csv")),row.names=FALSE)
  violin <- ggplot2::ggplot(long,ggplot2::aes(risk_group,expression,fill=risk_group))+ggplot2::geom_violin(trim=FALSE)+ggplot2::geom_boxplot(width=.16,fill="white",outlier.shape=NA)+ggplot2::facet_wrap(~gene,scales="free_y")+ggplot2::scale_fill_manual(values=pal)+ggplot2::theme_classic(base_size=7)+ggplot2::labs(x=NULL,y="Expression",fill="Risk group")
  save_plot(violin,file.path(output_dir,paste0(cancer,"_signature_violin")),183,120)
  invisible(list(km=km,roc=roc,tests=tests))
}
