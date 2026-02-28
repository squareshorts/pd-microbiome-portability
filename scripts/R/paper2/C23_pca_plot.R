cat("C23 PCA PLOTS STARTED\n")
cat("WD:", getwd(), "\n")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# Inputs (produced by C21)
clr_path  <- "data/processed/ipd_paper2/ipd_genus_shared2_clr.csv"
meta_path <- "data/processed/ipd_paper2/ipd_metadata_paper2.csv"
stopifnot(file.exists(clr_path), file.exists(meta_path))

X <- fread(clr_path)
meta <- fread(meta_path)

# Defensive: enforce SampleID name and type
if (!"SampleID" %in% names(X)) stop("CLR file missing SampleID column.")
if (!"SampleID" %in% names(meta)) stop("Metadata missing SampleID column.")

X[, SampleID := as.character(SampleID)]
meta[, SampleID := as.character(SampleID)]

# Align rows strictly to metadata order (same approach as C20b)
setkey(X, SampleID)
setkey(meta, SampleID)

dt <- X[meta[, .(SampleID)], on = "SampleID", nomatch = 0]
stopifnot(nrow(dt) == nrow(meta))

# PCA on CLR features (exclude SampleID)
feat_cols <- setdiff(names(dt), "SampleID")
M <- as.matrix(dt[, ..feat_cols])

# prcomp expects complete finite matrix
if (any(!is.finite(M))) stop("Non-finite values found in CLR matrix.")

pca <- prcomp(M, center = TRUE, scale. = FALSE)

ve <- (pca$sdev^2) / sum(pca$sdev^2)
pc1 <- round(100 * ve[1], 2)
pc2 <- round(100 * ve[2], 2)

scores <- as.data.table(pca$x[, 1:2])
setnames(scores, c("PC1", "PC2"))
scores[, SampleID := dt$SampleID]

plot_dt <- merge(scores, meta, by = "SampleID", all.x = TRUE)
stopifnot(nrow(plot_dt) == nrow(meta))

# Consistent PD labeling
if (!"PD" %in% names(plot_dt)) stop("Metadata missing PD column.")
plot_dt[, PD := as.integer(PD)]
plot_dt[, PD_label := ifelse(PD == 1, "PD", "Control")]

# Output folder
out_dir <- "results/paper2"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 1) PCA colored by Cohort
p_cohort <- ggplot(plot_dt, aes(x = PC1, y = PC2, color = Cohort)) +
  geom_point(alpha = 0.75, size = 1.6) +
  labs(
    x = paste0("PC1 (", pc1, "%)"),
    y = paste0("PC2 (", pc2, "%)")
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_blank()
  )

out1 <- file.path(out_dir, "fig_pca_cohort.png")
ggsave(out1, p_cohort, width = 7.2, height = 5.4, dpi = 300)
cat("Saved:", out1, "\n")

# 2) PCA colored by PD
p_pd <- ggplot(plot_dt, aes(x = PC1, y = PC2, color = PD_label)) +
  geom_point(alpha = 0.75, size = 1.6) +
  labs(
    x = paste0("PC1 (", pc1, "%)"),
    y = paste0("PC2 (", pc2, "%)")
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_blank()
  )

out2 <- file.path(out_dir, "fig_pca_pd.png")
ggsave(out2, p_pd, width = 7.2, height = 5.4, dpi = 300)
cat("Saved:", out2, "\n")

# Optional: quick text dump for manuscript
cat("\nVariance explained:\n")
cat("PC1:", ve[1], "\n")
cat("PC2:", ve[2], "\n")

cat("C23 PCA PLOTS DONE\n")