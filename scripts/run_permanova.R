.libPaths('R_libs')
library(vegan)

# Read datasets
clr_df <- read.csv("results/portability_analysis/dataset_clr.csv")

# Extract metadata and features
meta <- clr_df[, c("SampleID", "Cohort", "PD")]
meta$Cohort <- as.factor(meta$Cohort)
meta$PD <- as.factor(meta$PD)

features <- clr_df[, !(names(clr_df) %in% c("SampleID", "Cohort", "PD"))]

# Calculate Euclidean distance on CLR (which is Aitchison distance)
dist_matrix <- vegdist(features, method="euclidean")

# Run PERMANOVA: distance ~ cohort + status (marginal or sequential? the prompt said "report whether R2 values are from sequential terms or marginal tests"). Let's use marginal testing by = "margin".
print("Running PERMANOVA...")
set.seed(42)
perm_res <- adonis2(dist_matrix ~ Cohort + PD, data = meta, by = "margin", permutations = 999)
print(perm_res)

r2_cohort <- perm_res["Cohort", "R2"]
r2_pd <- perm_res["PD", "R2"]

cddr <- r2_cohort / r2_pd

res_df <- data.frame(
  Term = c("Cohort", "PD"),
  R2 = c(r2_cohort, r2_pd),
  P_value = c(perm_res["Cohort", "Pr(>F)"], perm_res["PD", "Pr(>F)"]),
  CDDR = c(cddr, NA)
)

write.csv(res_df, "results/portability_analysis/cohort_disease_dominance_ratio.csv", row.names=FALSE)
print("Saved CDDR.")
