cat("C30 PC1 LOADINGS STARTED\n")
cat("WD:", getwd(), "\n")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# Inputs (paper2-only)
x_path    <- "data/processed/ipd_paper2/ipd_genus_shared2_clr.csv"
meta_path <- "data/processed/ipd_paper2/ipd_metadata_paper2.csv"
stopifnot(file.exists(x_path), file.exists(meta_path))

Xdt  <- fread(x_path)
meta <- fread(meta_path)
setnames(Xdt, 1, "SampleID")

# Align rows strictly by SampleID
setkey(meta, SampleID)
setkey(Xdt, SampleID)
dt <- Xdt[meta[, .(SampleID)], on="SampleID", nomatch=0]
stopifnot(nrow(dt) == nrow(meta))

feat_cols <- setdiff(names(dt), "SampleID")
X <- as.matrix(dt[, ..feat_cols])

# PCA on CLR (already centered per-sample; still do column centering for PCA)
pca <- prcomp(X, center=TRUE, scale.=FALSE)

var_exp <- (pca$sdev^2) / sum(pca$sdev^2)
cat(sprintf("Variance explained PC1: %.4f\n", var_exp[1]))
cat(sprintf("Variance explained PC2: %.4f\n", var_exp[2]))

# Loadings for PC1
load1 <- pca$rotation[, 1]
load_dt <- data.table(
  taxon = names(load1),
  loading_pc1 = as.numeric(load1),
  abs_loading = abs(as.numeric(load1))
)
setorder(load_dt, -abs_loading)

# Save full loadings table
dir.create("results/paper2", recursive=TRUE, showWarnings=FALSE)
out_tab <- "results/paper2/table_pc1_loadings.csv"
fwrite(load_dt, out_tab)
cat("Saved:", out_tab, "\n")

# Top taxa (both directions)
top_n <- 20
top_pos <- load_dt[order(-loading_pc1)][1:top_n]
top_neg <- load_dt[order(loading_pc1)][1:top_n]

cat("\nTop +PC1 taxa (largest positive loadings):\n")
print(top_pos[, .(taxon, loading_pc1)])

cat("\nTop -PC1 taxa (most negative loadings):\n")
print(top_neg[, .(taxon, loading_pc1)])

# Plot: top absolute loadings (signed)
plot_n <- 30
plot_dt <- load_dt[1:plot_n]
# order for plotting
plot_dt[, taxon := factor(taxon, levels=plot_dt[order(loading_pc1)]$taxon)]

p <- ggplot(plot_dt, aes(x=taxon, y=loading_pc1)) +
  geom_col() +
  coord_flip() +
  labs(x=NULL, y="PC1 loading (signed)") +
  theme_bw(base_size=12)

fig_path <- "results/paper2/fig_pc1_top_loadings.png"
ggsave(fig_path, p, width=7.2, height=7.5, dpi=300)
cat("Saved:", fig_path, "\n")

# Optional: correlate PC1 scores with cohort means (quick numeric sanity)
scores <- data.table(SampleID=dt$SampleID, PC1=pca$x[,1], PC2=pca$x[,2])
scores <- scores[meta, on="SampleID"]
coh_summary <- scores[, .(
  n=.N,
  PC1_mean=mean(PC1),
  PC1_sd=sd(PC1)
), by=Cohort][order(PC1_mean)]
out_coh <- "results/paper2/table_pc1_by_cohort.csv"
fwrite(coh_summary, out_coh)
cat("Saved:", out_coh, "\n")
print(coh_summary)

cat("C30 PC1 LOADINGS DONE\n")