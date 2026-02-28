cat("C22 ECOLOGY STRUCTURE STARTED\n")
cat("WD:", getwd(), "\n")

library(data.table)
library(vegan)

# ---- Paths ----
clr_path  <- "data/processed/ipd_paper2/ipd_genus_shared2_clr.csv"
meta_path <- "data/processed/ipd_paper2/ipd_metadata_paper2.csv"
# ---- Load ----
clr  <- fread(clr_path)
meta <- fread(meta_path)

stopifnot(all(clr$SampleID == meta$SampleID))

X <- as.matrix(clr[, !"SampleID"])
mode(X) <- "numeric"

cat("Matrix dim:", dim(X), "\n")

# ---- PCA ----
pca <- prcomp(X, center=FALSE, scale.=FALSE)

var_expl <- (pca$sdev^2) / sum(pca$sdev^2)

cat("Variance explained PC1:", round(var_expl[1], 4), "\n")
cat("Variance explained PC2:", round(var_expl[2], 4), "\n")
cat("Variance explained PC3:", round(var_expl[3], 4), "\n")

# ---- PERMANOVA ----
set.seed(1)

dist_mat <- dist(X)

perm_cohort <- adonis2(dist_mat ~ Cohort, data=meta, permutations=999)
perm_pd     <- adonis2(dist_mat ~ PD, data=meta, permutations=999)

cat("\nPERMANOVA Cohort:\n")
print(perm_cohort)

cat("\nPERMANOVA PD:\n")
print(perm_pd)

cat("C22 ECOLOGY STRUCTURE DONE\n")