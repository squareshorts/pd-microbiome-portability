#!/usr/bin/env Rscript

# ============================================================
# C1_calibration_loco_ipd.R
# LOCO calibration analysis (elastic net)
# ============================================================

library(data.table)
library(compositions)
library(glmnet)
library(dplyr)
library(ggplot2)

counts_path <- "data/processed/ipd_sex/ipd_genus_counts_union.csv"
meta_path   <- "data/processed/ipd_sex/ipd_metadata.csv"

counts <- fread(counts_path)
meta   <- fread(meta_path)

counts <- as.data.frame(counts)
meta   <- as.data.frame(meta)

# ------------------------------------------------------------
# Remove non-numeric columns
# ------------------------------------------------------------

numeric_cols <- sapply(counts, is.numeric)
counts <- counts[, numeric_cols]
counts <- as.matrix(counts)

if(!all(c("PD","Cohort") %in% colnames(meta))) {
  stop("Metadata must contain PD and Cohort columns.")
}

meta$cohort <- factor(meta$Cohort)
y <- meta$PD

# ------------------------------------------------------------
# CLR transform
# ------------------------------------------------------------

X <- clr(counts + 1)

cohorts <- unique(meta$cohort)
prob_loco <- rep(NA, length(y))

cat("Running LOCO elastic net models...\n")

for(c in cohorts) {

  train_idx <- which(meta$cohort != c)
  test_idx  <- which(meta$cohort == c)

  X_train <- X[train_idx, ]
  y_train <- y[train_idx]
  X_test  <- X[test_idx, ]

  cv_fit <- cv.glmnet(
    X_train, y_train,
    family = "binomial",
    alpha = 0.5,
    nfolds = 5
  )

  best_lambda <- cv_fit$lambda.1se

  model <- glmnet(
    X_train, y_train,
    family = "binomial",
    alpha = 0.5,
    lambda = best_lambda
  )

  prob_loco[test_idx] <- predict(
    model, X_test,
    type = "response"
  )[,1]
}

# ------------------------------------------------------------
# Brier Score
# ------------------------------------------------------------

brier <- mean((prob_loco - y)^2)

cat("LOCO Brier score:", brier, "\n")

write.csv(data.frame(Brier_LOCO = brier),
          "results/paper2_portability/brier_loco.csv",
          row.names = FALSE)

# ------------------------------------------------------------
# Calibration Curve
# ------------------------------------------------------------

cal_df <- data.frame(prob = prob_loco, y = y)

cal_df <- cal_df %>%
  mutate(bin = ntile(prob, 10)) %>%
  group_by(bin) %>%
  summarise(
    mean_pred = mean(prob),
    observed  = mean(y),
    .groups = "drop"
  )

p <- ggplot(cal_df,
            aes(x = mean_pred, y = observed)) +
  geom_point(size = 3) +
  geom_line() +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed") +
  labs(
    x = "Mean Predicted Probability",
    y = "Observed Event Rate",
    title = "LOCO Calibration Curve"
  ) +
  theme_minimal()

ggsave("results/paper2_portability/calibration_loco.png",
       p, width = 6, height = 5)

write.csv(cal_df,
          "results/paper2_portability/calibration_loco_table.csv",
          row.names = FALSE)

cat("Calibration analysis complete.\n")
