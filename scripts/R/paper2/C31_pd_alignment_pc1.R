cat("C31 PD ALIGNMENT WITH PC1 STARTED\n")
cat("WD:", getwd(), "\n")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# Inputs
x_path    <- "data/processed/ipd_paper2/ipd_genus_shared2_clr.csv"
meta_path <- "data/processed/ipd_paper2/ipd_metadata_paper2.csv"
stopifnot(file.exists(x_path), file.exists(meta_path))

Xdt  <- fread(x_path)
meta <- fread(meta_path)
setnames(Xdt, 1, "SampleID")

stopifnot(all(c("SampleID","PD","Cohort") %in% names(meta)))

# Align rows strictly
setkey(meta, SampleID)
setkey(Xdt, SampleID)
dt <- Xdt[meta[, .(SampleID)], on="SampleID", nomatch=0]
stopifnot(nrow(dt) == nrow(meta))

feat_cols <- setdiff(names(dt), "SampleID")
X <- as.matrix(dt[, ..feat_cols])

# PCA
pca <- prcomp(X, center=TRUE, scale.=FALSE)
var_exp <- (pca$sdev^2) / sum(pca$sdev^2)

scores <- data.table(
  SampleID = dt$SampleID,
  PC1 = as.numeric(pca$x[,1]),
  PC2 = as.numeric(pca$x[,2]),
  PC3 = as.numeric(pca$x[,3])
)[meta, on="SampleID"]

cat(sprintf("Variance explained: PC1=%.4f PC2=%.4f PC3=%.4f\n", var_exp[1], var_exp[2], var_exp[3]))

# Helper: Cohen's d for two groups (PD=1 vs PD=0)
cohen_d <- function(x, g01) {
  x0 <- x[g01==0]; x1 <- x[g01==1]
  n0 <- length(x0); n1 <- length(x1)
  m0 <- mean(x0); m1 <- mean(x1)
  s0 <- sd(x0); s1 <- sd(x1)
  sp <- sqrt(((n0-1)*s0^2 + (n1-1)*s1^2) / (n0+n1-2))
  d <- (m1 - m0) / sp
  list(d=d, m0=m0, m1=m1, s0=s0, s1=s1, n0=n0, n1=n1)
}

# 1) Pooled model: PD ~ PC1 + Cohort (tests PD alignment with PC1 controlling cohort)
# Use logistic regression with cohort fixed effects.
fit1 <- glm(PD ~ PC1 + Cohort, data=scores, family=binomial())
s1 <- summary(fit1)$coefficients
pc1_beta <- s1["PC1","Estimate"]
pc1_se   <- s1["PC1","Std. Error"]
pc1_p    <- s1["PC1","Pr(>|z|)"]
pc1_or   <- exp(pc1_beta)

# Also show Cohort-only baseline for context
fit0 <- glm(PD ~ Cohort, data=scores, family=binomial())
ll0 <- as.numeric(logLik(fit0))
ll1 <- as.numeric(logLik(fit1))
lr_stat <- 2*(ll1-ll0)
lr_p <- pchisq(lr_stat, df=1, lower.tail=FALSE)

cat("\nPooled logistic model: PD ~ PC1 + Cohort\n")
cat(sprintf("PC1 beta=%.4f (SE=%.4f) OR=%.3f p=%.3g\n", pc1_beta, pc1_se, pc1_or, pc1_p))
cat(sprintf("LR test vs Cohort-only (adds PC1): LR=%.3f p=%.3g\n", lr_stat, lr_p))

# 2) Within-cohort PD separation along PCs (effect sizes + t-test)
cohorts <- sort(unique(scores$Cohort))
within_rows <- list()
for (coh in cohorts) {
  dcoh <- scores[Cohort==coh]
  if (length(unique(dcoh$PD)) < 2) next

  # PC1 / PC2 / PC3
  for (pc in c("PC1","PC2","PC3")) {
    x <- dcoh[[pc]]
    g <- dcoh$PD
    dd <- cohen_d(x, g)
    tt <- t.test(x ~ g)
    within_rows[[length(within_rows)+1]] <- data.table(
      cohort=coh,
      pc=pc,
      n=nrow(dcoh),
      prev=mean(g),
      mean_PD0=dd$m0,
      mean_PD1=dd$m1,
      d=dd$d,
      p_value=tt$p.value
    )
  }
}
within_tab <- rbindlist(within_rows)

# Multiple-testing correction within this table (BH)
within_tab[, p_adj := p.adjust(p_value, method="BH")]
setorder(within_tab, cohort, pc)

dir.create("results/paper2", recursive=TRUE, showWarnings=FALSE)
out_within <- "results/paper2/table_pd_effect_on_pcs_within_cohort.csv"
fwrite(within_tab, out_within)
cat("Saved:", out_within, "\n")

# 3) Figure: PC1 vs PC2 scatter, colored by cohort, shaped by PD
p_scatter <- ggplot(scores, aes(x=PC1, y=PC2, color=Cohort, shape=factor(PD))) +
  geom_point(alpha=0.75, size=2) +
  theme_bw(base_size=12) +
  labs(shape="PD", x="PC1 score", y="PC2 score",
       title=NULL,
       subtitle=sprintf("Explained variance: PC1=%.1f%%, PC2=%.1f%%", 100*var_exp[1], 100*var_exp[2]))

fig1 <- "results/paper2/fig_pc1_pc2_cohort_pd.png"
ggsave(fig1, p_scatter, width=7.2, height=5.4, dpi=300)
cat("Saved:", fig1, "\n")

# 4) Figure: within-cohort PD shift on PC1 (boxplots)
p_box <- ggplot(scores, aes(x=Cohort, y=PC1, fill=factor(PD))) +
  geom_boxplot(outlier.size=0.8) +
  theme_bw(base_size=12) +
  labs(fill="PD", x=NULL, y="PC1 score", title=NULL)

fig2 <- "results/paper2/fig_pc1_by_cohort_pd_box.png"
ggsave(fig2, p_box, width=7.2, height=5.4, dpi=300)
cat("Saved:", fig2, "\n")

cat("C31 PD ALIGNMENT WITH PC1 DONE\n")