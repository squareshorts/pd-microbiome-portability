library(data.table)

meta <- fread("data/processed/ipd_sex/ipd_metadata.csv")
X <- fread("data/processed/ipd_sex/ipd_genus_counts_union.csv")
setkey(meta, SampleID); setkey(X, SampleID)
dt <- X[meta]

genus_cols <- setdiff(names(dt), c("SampleID","PD","Sex","Cohort"))

# Prevalence filter to avoid near-zero taxa
prev <- colMeans(as.matrix(dt[, ..genus_cols]) > 0)
keep <- genus_cols[prev >= 0.10]  # 10% overall prevalence
cat("Genera kept (>=10% prevalence overall):", length(keep), "\n")

# Function: CLR within cohort (using cohort-specific row means)
clr_transform <- function(mat) {
  mat_pc <- mat + 1
  log_mat <- log(mat_pc)
  log_mat - rowMeans(log_mat)
}

cohorts <- sort(unique(dt$Cohort))

# Cohortwise fits
cohort_res <- rbindlist(lapply(cohorts, function(coh) {
  d <- dt[Cohort == coh]
  mat <- as.matrix(d[, ..keep]); storage.mode(mat) <- "numeric"
  clr <- clr_transform(mat)

  df <- data.table(PD=factor(d$PD), Sex=factor(d$Sex))
  df[, Sex := relevel(Sex, ref="Female")]
  df[, PD := relevel(PD, ref="0")]

  rbindlist(lapply(seq_along(keep), function(j){
    g <- keep[j]
    y <- clr[, j]
    fit <- lm(y ~ PD*Sex, data=df)
    coefs <- summary(fit)$coefficients
    if (!("PD1:SexMale" %in% rownames(coefs))) return(NULL)
    data.table(
      Cohort=coh,
      Genus=g,
      beta=coefs["PD1:SexMale",1],
      se=coefs["PD1:SexMale",2],
      p=coefs["PD1:SexMale",4]
    )
  }))
}))

# Fixed-effect meta-analysis on interaction betas
# beta_FE = sum(w*beta)/sum(w), w=1/se^2; se_FE=sqrt(1/sum(w))
meta_res <- cohort_res[, {
  w <- 1/(se^2)
  beta_fe <- sum(w*beta)/sum(w)
  se_fe <- sqrt(1/sum(w))
  z <- beta_fe/se_fe
  p_fe <- 2*pnorm(-abs(z))
  .(beta_FE=beta_fe, se_FE=se_fe, z=z, p_FE=p_fe,
    k=.N)
}, by=.(Genus)]

meta_res[, q_FE := p.adjust(p_FE, method="BH")]
setorder(meta_res, p_FE)

dir.create("results/sex_paper", showWarnings=FALSE, recursive=TRUE)
fwrite(cohort_res, "results/sex_paper/cohortwise_interactions.csv")
fwrite(meta_res, "results/sex_paper/meta_interaction_FE.csv")

cat("Top 20 meta-analyzed interaction hits:\n")
print(meta_res[1:20])
