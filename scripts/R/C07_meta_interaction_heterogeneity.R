library(data.table)

cohort_res <- fread("results/sex_paper/cohortwise_interactions.csv")

# Fixed-effect meta + heterogeneity stats
meta_res <- cohort_res[, {
  w <- 1/(se^2)
  beta_fe <- sum(w*beta)/sum(w)
  se_fe <- sqrt(1/sum(w))
  z <- beta_fe/se_fe
  p_fe <- 2*pnorm(-abs(z))

  # Cochran's Q and I^2
  Q <- sum(w * (beta - beta_fe)^2)
  df <- .N - 1
  p_Q <- pchisq(Q, df=df, lower.tail=FALSE)
  I2 <- if (Q > 0) max(0, (Q - df)/Q) else 0

  .(beta_FE=beta_fe, se_FE=se_fe, z=z, p_FE=p_fe,
    Q=Q, df=df, p_Q=p_Q, I2=I2, k=.N)
}, by=.(Genus)]

meta_res[, q_FE := p.adjust(p_FE, method="BH")]
setorder(meta_res, p_FE)

dir.create("results/sex_paper", showWarnings=FALSE, recursive=TRUE)
fwrite(meta_res, "results/sex_paper/meta_interaction_FE_heterogeneity.csv")

# Forest-plot-ready table: cohort betas + SE + meta + heterogeneity
forest <- merge(cohort_res, meta_res, by="Genus", all.x=TRUE)
setorder(forest, p_FE, Genus, Cohort)
fwrite(forest, "results/sex_paper/forest_table_interactions.csv")

# Print top 20 with heterogeneity
cat("Top 20 meta interaction hits with heterogeneity:\n")
print(meta_res[1:20, .(Genus, beta_FE, se_FE, p_FE, q_FE, Q, df, p_Q, I2, k)])

# Also print the subset with low heterogeneity among the most significant
cat("\nTop hits with I2 <= 0.25 (low heterogeneity):\n")
print(meta_res[I2 <= 0.25][1:20, .(Genus, beta_FE, se_FE, p_FE, q_FE, I2, p_Q, k)])
