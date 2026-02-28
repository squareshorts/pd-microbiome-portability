cat("C32 STRUCTURAL ASYMMETRY SUMMARY STARTED\n")
cat("WD:", getwd(), "\n")

suppressPackageStartupMessages({
  library(data.table)
})

# ---- Inputs ----
eco_path   <- "results/paper2/table_loco_logreg_summary.csv"  # only to confirm file tree
pc1_by_coh <- "results/paper2/table_pc1_by_cohort.csv"
pd_pc_tab  <- "results/paper2/table_pd_effect_on_pcs_within_cohort.csv"

stopifnot(file.exists(pc1_by_coh),
          file.exists(pd_pc_tab))

pc1_coh <- fread(pc1_by_coh)
pd_tab  <- fread(pd_pc_tab)

# ---- Hard-coded from C22 output (PERMANOVA) ----
R2_cohort <- 0.68918
R2_PD     <- 0.00923

ratio_R2  <- R2_cohort / R2_PD

cat("\n--- PERMANOVA STRUCTURE ---\n")
cat(sprintf("Cohort R2 = %.5f\n", R2_cohort))
cat(sprintf("PD R2     = %.5f\n", R2_PD))
cat(sprintf("Cohort explains %.1f× more variance than PD\n", ratio_R2))

# ---- PC1 variance from C30 ----
# We recompute quickly for robustness

x_path    <- "data/processed/ipd_paper2/ipd_genus_shared2_clr.csv"
meta_path <- "data/processed/ipd_paper2/ipd_metadata_paper2.csv"
Xdt  <- fread(x_path)
meta <- fread(meta_path)
setnames(Xdt, 1, "SampleID")

setkey(meta, SampleID)
setkey(Xdt, SampleID)
dt <- Xdt[meta[, .(SampleID)], on="SampleID", nomatch=0]

feat_cols <- setdiff(names(dt), "SampleID")
X <- as.matrix(dt[, ..feat_cols])

pca <- prcomp(X, center=TRUE, scale.=FALSE)
var_exp <- (pca$sdev^2) / sum(pca$sdev^2)

cat("\n--- PCA STRUCTURE ---\n")
cat(sprintf("PC1 variance explained = %.4f (%.1f%%)\n",
            var_exp[1], 100*var_exp[1]))

# ---- PD alignment with PC1 (from C31 model) ----
scores <- data.table(
  SampleID = dt$SampleID,
  PC1 = as.numeric(pca$x[,1])
)[meta, on="SampleID"]

fit1 <- glm(PD ~ PC1 + Cohort, data=scores, family=binomial())
fit0 <- glm(PD ~ Cohort, data=scores, family=binomial())

s1 <- summary(fit1)$coefficients
beta_pc1 <- s1["PC1","Estimate"]
se_pc1   <- s1["PC1","Std. Error"]
p_pc1    <- s1["PC1","Pr(>|z|)"]
OR_pc1   <- exp(beta_pc1)

ll0 <- as.numeric(logLik(fit0))
ll1 <- as.numeric(logLik(fit1))
lr_stat <- 2*(ll1-ll0)
lr_p <- pchisq(lr_stat, df=1, lower.tail=FALSE)

cat("\n--- PD ~ PC1 + Cohort ---\n")
cat(sprintf("PC1 beta = %.4f (SE=%.4f)\n", beta_pc1, se_pc1))
cat(sprintf("Odds ratio per unit PC1 = %.3f\n", OR_pc1))
cat(sprintf("Wald p = %.4f\n", p_pc1))
cat(sprintf("LR test (add PC1) p = %.4f\n", lr_p))

# ---- Extract strongest within-cohort PD effects ----
pd_tab_sig <- pd_tab[p_adj < 0.01]

cat("\n--- WITHIN-COHORT PD EFFECTS (BH p<0.01) ---\n")
print(pd_tab_sig)

# ---- Save compact structural summary table ----
summary_tab <- data.table(
  metric = c("PERMANOVA_R2_Cohort",
             "PERMANOVA_R2_PD",
             "R2_ratio_Cohort_over_PD",
             "PC1_variance_explained",
             "PC1_beta_PD_adjCohort",
             "PC1_OR_PD_adjCohort",
             "PC1_p_PD_adjCohort",
             "PC1_LR_p_adjCohort"),
  value = c(R2_cohort,
            R2_PD,
            ratio_R2,
            var_exp[1],
            beta_pc1,
            OR_pc1,
            p_pc1,
            lr_p)
)

dir.create("results/paper2", recursive=TRUE, showWarnings=FALSE)
out_path <- "results/paper2/table_structural_asymmetry_summary.csv"
fwrite(summary_tab, out_path)
cat("\nSaved:", out_path, "\n")

cat("\nC32 STRUCTURAL ASYMMETRY SUMMARY DONE\n")