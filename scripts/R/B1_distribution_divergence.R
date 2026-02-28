#!/usr/bin/env Rscript

# ============================================================
# B1_distribution_divergence_ipd.R
# Distribution divergence using IPD union dataset
# ============================================================

library(data.table)
library(compositions)
library(vegan)

cat("Loading IPD union counts and metadata...\n")

counts_path <- "data/processed/ipd_sex/ipd_genus_counts_union.csv"
meta_path   <- "data/processed/ipd_sex/ipd_metadata.csv"

counts <- fread(counts_path)
meta   <- fread(meta_path)

counts <- as.data.frame(counts)
meta   <- as.data.frame(meta)

cat("Rows in counts:", nrow(counts), "\n")
cat("Rows in metadata:", nrow(meta), "\n")

if(nrow(counts) != nrow(meta)) {
  stop("Counts and metadata row mismatch.")
}

# ------------------------------------------------------------
# Remove non-numeric columns (e.g., SampleID)
# ------------------------------------------------------------

numeric_cols <- sapply(counts, is.numeric)
counts <- counts[, numeric_cols]

cat("Numeric genus columns retained:", ncol(counts), "\n")

# ------------------------------------------------------------
# Prepare metadata
# ------------------------------------------------------------

if(!all(c("PD","Cohort") %in% colnames(meta))) {
  stop("Metadata must contain PD and Cohort columns.")
}

meta$cohort <- factor(meta$Cohort)
meta$PD     <- factor(meta$PD)

# ------------------------------------------------------------
# CLR transformation
# ------------------------------------------------------------

cat("Performing CLR transformation...\n")
counts_clr <- clr(as.matrix(counts) + 1)

# ------------------------------------------------------------
# Aitchison distance
# ------------------------------------------------------------

cat("Computing Aitchison distance matrix...\n")
dist_mat <- dist(counts_clr, method = "euclidean")

# ------------------------------------------------------------
# PERMANOVA (Sequential / Type I)
# ------------------------------------------------------------

cat("Running PERMANOVA (sequential effects)...\n")

perm_seq <- adonis2(
  dist_mat ~ cohort + PD,
  data = meta,
  permutations = 999
)

print(perm_seq)

write.csv(as.data.frame(perm_seq),
          "results/paper2_portability/permanova_sequential.csv")

# ------------------------------------------------------------
# PERMANOVA (Marginal / Type III-like)
# ------------------------------------------------------------

cat("Running PERMANOVA (marginal effects)...\n")

perm_margin <- adonis2(
  dist_mat ~ cohort + PD,
  data = meta,
  permutations = 999,
  by = "margin"
)

print(perm_margin)

write.csv(as.data.frame(perm_margin),
          "results/paper2_portability/permanova_marginal.csv")

# ------------------------------------------------------------
# Cohort centroid distances
# ------------------------------------------------------------

cat("Computing cohort centroids...\n")

centroids <- aggregate(counts_clr,
                       by = list(meta$cohort),
                       FUN = mean)

rownames(centroids) <- centroids$Group.1
centroids <- as.matrix(centroids[, -1])

centroid_dist <- as.matrix(dist(centroids, method = "euclidean"))

print(centroid_dist)

write.csv(as.data.frame(centroid_dist),
          "results/paper2_portability/centroid_aitchison_distances.csv")

cat("Distribution divergence analysis complete.\n")
