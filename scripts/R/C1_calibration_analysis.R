#!/usr/bin/env Rscript

# ============================================================
# C1_calibration_analysis.R
# Computes calibration metrics under OOF and LOCO
# ============================================================

library(data.table)
library(dplyr)
library(ggplot2)

cat("Loading metadata with predictions...\n")

meta <- fread("data/metadata.csv")
meta <- as.data.frame(meta)

if(!all(c("PD","prob_oof","prob_loco") %in% colnames(meta))) {
  stop("metadata.csv must contain PD, prob_oof, prob_loco columns.")
}

y_true <- meta$PD

# ------------------------------
# BRIER SCORE
# ------------------------------

brier_score <- function(p, y) {
  mean((p - y)^2)
}

brier_oof  <- brier_score(meta$prob_oof,  y_true)
brier_loco <- brier_score(meta$prob_loco, y_true)

cat("\nBrier Scores:\n")
cat("OOF  :", brier_oof,  "\n")
cat("LOCO :", brier_loco, "\n")

write.csv(data.frame(
  Brier_OOF = brier_oof,
  Brier_LOCO = brier_loco
), "results/brier_scores.csv", row.names = FALSE)

# ------------------------------
# CALIBRATION FUNCTION
# ------------------------------

compute_calibration <- function(prob, y, label) {

  df <- data.frame(prob = prob, y = y)

  df <- df %>%
    mutate(bin = ntile(prob, 10)) %>%
    group_by(bin) %>%
    summarise(
      mean_pred = mean(prob),
      observed  = mean(y),
      .groups = "drop"
    )

  df$model <- label
  return(df)
}

cal_oof  <- compute_calibration(meta$prob_oof,  y_true, "OOF")
cal_loco <- compute_calibration(meta$prob_loco, y_true, "LOCO")

calibration_all <- rbind(cal_oof, cal_loco)

# ------------------------------
# PLOT
# ------------------------------

p <- ggplot(calibration_all,
            aes(x = mean_pred, y = observed, color = model)) +
  geom_point(size = 3) +
  geom_line() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(
    x = "Mean Predicted Probability",
    y = "Observed Event Rate",
    title = "Calibration Curves: OOF vs LOCO"
  ) +
  theme_minimal()

ggsave("results/calibration_plot.png", p, width = 6, height = 5)

write.csv(calibration_all,
          "results/calibration_table.csv",
          row.names = FALSE)

cat("\nCalibration analysis complete.\n")
