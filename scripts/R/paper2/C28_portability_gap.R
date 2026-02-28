cat("C28 PORTABILITY GAP STARTED\n")
cat("WD:", getwd(), "\n")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# Inputs
loco_lr_path   <- "results/paper2/table_loco_logreg_summary.csv"
loco_rf_path   <- "results/paper2/table_loco_rf_summary.csv"
within_lr_path <- "results/paper2/table_within_cohort_logreg_summary.csv"
within_rf_path <- "results/paper2/table_within_cohort_rf_summary.csv"

stopifnot(file.exists(loco_lr_path), file.exists(loco_rf_path),
          file.exists(within_lr_path), file.exists(within_rf_path))

loco_lr   <- fread(loco_lr_path)
loco_rf   <- fread(loco_rf_path)
within_lr <- fread(within_lr_path)
within_rf <- fread(within_rf_path)

# Standardize cohort column name for joining
setnames(loco_lr, "heldout_cohort", "cohort")
setnames(loco_rf, "heldout_cohort", "cohort")

# Add model_family label
loco_lr[, model_family := "LogReg (ridge)"]
loco_rf[, model_family := "RF (ranger)"]
within_lr[, model_family := "LogReg (ridge)"]
within_rf[, model_family := "RF (ranger)"]

# Keep only needed columns
loco_lr_k <- loco_lr[, .(model_family, cohort,
                         auroc_loco = auroc, auprc_loco = auprc, f1_loco = f1,
                         balacc_loco = balacc, sens_loco = sensitivity, spec_loco = specificity,
                         n_test)]
loco_rf_k <- loco_rf[, .(model_family, cohort,
                         auroc_loco = auroc, auprc_loco = auprc, f1_loco = f1,
                         balacc_loco = balacc, sens_loco = sensitivity, spec_loco = specificity,
                         n_test)]

within_lr_k <- within_lr[, .(model_family, cohort,
                             auroc_within = auroc_mean, auroc_within_sd = auroc_sd,
                             auprc_within = auprc_mean, auprc_within_sd = auprc_sd,
                             f1_within = f1_mean, f1_within_sd = f1_sd,
                             balacc_within = balacc_mean, balacc_within_sd = balacc_sd,
                             n)]
within_rf_k <- within_rf[, .(model_family, cohort,
                             auroc_within = auroc_mean, auroc_within_sd = auroc_sd,
                             auprc_within = auprc_mean, auprc_within_sd = auprc_sd,
                             f1_within = f1_mean, f1_within_sd = f1_sd,
                             balacc_within = balacc_mean, balacc_within_sd = balacc_sd,
                             n)]

loco_all   <- rbindlist(list(loco_lr_k, loco_rf_k), use.names=TRUE, fill=TRUE)
within_all <- rbindlist(list(within_lr_k, within_rf_k), use.names=TRUE, fill=TRUE)

gap <- merge(within_all, loco_all, by=c("model_family","cohort"), all=FALSE)

# Portability gaps (within minus LOCO)
gap[, `:=`(
  d_auroc  = auroc_within - auroc_loco,
  d_auprc  = auprc_within - auprc_loco,
  d_f1     = f1_within - f1_loco,
  d_balacc = balacc_within - balacc_loco
)]

dir.create("results/paper2", recursive=TRUE, showWarnings=FALSE)
out_tab <- "results/paper2/table_portability_gap.csv"
fwrite(gap, out_tab)
cat("Saved:", out_tab, "\n")

# Figure: AUROC within vs LOCO, faceted by model family
long <- rbindlist(list(
  gap[, .(model_family, cohort, setting="Within-cohort", metric="AUROC", value=auroc_within, sd=auroc_within_sd)],
  gap[, .(model_family, cohort, setting="LOCO",         metric="AUROC", value=auroc_loco,   sd=NA_real_)]
), use.names=TRUE, fill=TRUE)

p <- ggplot(long, aes(x=cohort, y=value, fill=setting)) +
  geom_col(position="dodge") +
  facet_wrap(~model_family, nrow=1) +
  ylim(0, 1) +
  labs(x=NULL, y="AUROC") +
  theme_bw(base_size=12)

fig_path <- "results/paper2/fig_auroc_within_vs_loco.png"
ggsave(fig_path, p, width=9.5, height=3.3, dpi=300)
cat("Saved:", fig_path, "\n")

# Figure: portability gap (ΔAUROC)
p2 <- ggplot(gap, aes(x=cohort, y=d_auroc, fill=model_family)) +
  geom_col(position="dodge") +
  labs(x=NULL, y="ΔAUROC (Within − LOCO)") +
  theme_bw(base_size=12)

fig2_path <- "results/paper2/fig_delta_auroc.png"
ggsave(fig2_path, p2, width=7.0, height=3.6, dpi=300)
cat("Saved:", fig2_path, "\n")

cat("C28 PORTABILITY GAP DONE\n")