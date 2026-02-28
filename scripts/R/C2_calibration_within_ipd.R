#!/usr/bin/env Rscript

library(data.table)
library(compositions)
library(glmnet)

counts_path <- "data/processed/ipd_sex/ipd_genus_counts_union.csv"
meta_path   <- "data/processed/ipd_sex/ipd_metadata.csv"

counts <- fread(counts_path)
meta   <- fread(meta_path)

counts <- as.data.frame(counts)
meta   <- as.data.frame(meta)

numeric_cols <- sapply(counts, is.numeric)
counts <- as.matrix(counts[, numeric_cols])

meta$cohort <- factor(meta$Cohort)
y <- meta$PD

X <- clr(counts + 1)

prob_within <- rep(NA, length(y))

for(c in unique(meta$cohort)) {

  idx <- which(meta$cohort == c)

  X_sub <- X[idx, ]
  y_sub <- y[idx]

  cv_fit <- cv.glmnet(
    X_sub, y_sub,
    family = "binomial",
    alpha = 0.5,
    nfolds = 5
  )

  best_lambda <- cv_fit$lambda.1se

  model <- glmnet(
    X_sub, y_sub,
    family = "binomial",
    alpha = 0.5,
    lambda = best_lambda
  )

  prob_within[idx] <- predict(
    model, X_sub,
    type = "response"
  )[,1]
}

brier_within <- mean((prob_within - y)^2)

cat("Within-cohort Brier score:", brier_within, "\n")

write.csv(data.frame(Brier_Within = brier_within),
          "results/paper2_portability/brier_within.csv",
          row.names = FALSE)
