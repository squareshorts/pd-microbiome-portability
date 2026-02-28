cat("C33 PD EFFECT ALIGNMENT ANGLES STARTED\n")
cat("WD:", getwd(), "\n")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ---- Inputs ----
x_path    <- "data/processed/ipd_paper2/ipd_genus_shared2_clr.csv"
meta_path <- "data/processed/ipd_paper2/ipd_metadata_paper2.csv"
stopifnot(file.exists(x_path), file.exists(meta_path))

Xdt  <- fread(x_path)
meta <- fread(meta_path)
setnames(Xdt, 1, "SampleID")

stopifnot(all(c("SampleID","PD","Cohort") %in% names(meta)))

# Align
setkey(meta, SampleID)
setkey(Xdt, SampleID)
dt <- Xdt[meta[, .(SampleID)], on="SampleID", nomatch=0]
stopifnot(nrow(dt) == nrow(meta))

feat_cols <- setdiff(names(dt), "SampleID")
X <- as.matrix(dt[, ..feat_cols])
y <- as.integer(meta$PD)
stopifnot(all(y %in% c(0,1)))

# PCA basis on pooled data (dominant cohort structure lives here)
pca <- prcomp(X, center=TRUE, scale.=FALSE)
V <- pca$rotation  # loadings: features x PCs
var_exp <- (pca$sdev^2) / sum(pca$sdev^2)

cat(sprintf("Variance explained: PC1=%.4f PC2=%.4f PC3=%.4f\n",
            var_exp[1], var_exp[2], var_exp[3]))

# Helper: angle between vectors
cos_sim <- function(a, b) {
  den <- sqrt(sum(a*a)) * sqrt(sum(b*b))
  if (den == 0) return(NA_real_)
  sum(a*b) / den
}
angle_deg <- function(cosv) {
  if (is.na(cosv)) return(NA_real_)
  acos(pmax(-1, pmin(1, cosv))) * 180/pi
}

# Permutation test for alignment with PC1 within a cohort
perm_p_alignment <- function(Xc, yc, pc_vec, B=2000, seed=1) {
  set.seed(seed)
  idx1 <- which(yc==1); idx0 <- which(yc==0)
  if (length(idx1)<5 || length(idx0)<5) return(NA_real_)

  d_obs <- colMeans(Xc[idx1,,drop=FALSE]) - colMeans(Xc[idx0,,drop=FALSE])
  c_obs <- abs(cos_sim(d_obs, pc_vec))

  c_perm <- numeric(B)
  for (b in seq_len(B)) {
    yperm <- sample(yc)
    i1 <- which(yperm==1); i0 <- which(yperm==0)
    d_b <- colMeans(Xc[i1,,drop=FALSE]) - colMeans(Xc[i0,,drop=FALSE])
    c_perm[b] <- abs(cos_sim(d_b, pc_vec))
  }
  # two-sided in terms of absolute alignment
  mean(c_perm >= c_obs)
}

cohorts <- sort(unique(meta$Cohort))

rows <- list()

pc1_vec <- V[,1]
pc2_vec <- V[,2]
pc3_vec <- V[,3]

for (coh in cohorts) {
  idx <- which(meta$Cohort == coh)
  Xc <- X[idx,,drop=FALSE]
  yc <- y[idx]

  n <- length(yc)
  prev <- mean(yc)

  idx1 <- which(yc==1); idx0 <- which(yc==0)
  stopifnot(length(idx1) + length(idx0) == n)

  d_vec <- colMeans(Xc[idx1,,drop=FALSE]) - colMeans(Xc[idx0,,drop=FALSE])

  c1 <- cos_sim(d_vec, pc1_vec); a1 <- angle_deg(abs(c1))
  c2 <- cos_sim(d_vec, pc2_vec); a2 <- angle_deg(abs(c2))
  c3 <- cos_sim(d_vec, pc3_vec); a3 <- angle_deg(abs(c3))

  # alignment “energy” captured by first K PCs (projection length)
  proj1 <- sum(d_vec * pc1_vec)
  proj2 <- sum(d_vec * pc2_vec)
  proj3 <- sum(d_vec * pc3_vec)
  proj_norm <- sqrt(proj1^2 + proj2^2 + proj3^2)
  d_norm <- sqrt(sum(d_vec*d_vec))
  frac_1to3 <- ifelse(d_norm==0, NA_real_, proj_norm / d_norm)

  p_align_pc1 <- perm_p_alignment(Xc, yc, pc1_vec, B=2000, seed=1)

  cat("\n---", coh, "---\n")
  cat(sprintf("n=%d prev=%.3f\n", n, prev))
  cat(sprintf("abs(cos) with PC1=%.3f angle=%.1f deg | perm p=%.4f\n",
              abs(c1), a1, p_align_pc1))
  cat(sprintf("abs(cos) with PC2=%.3f angle=%.1f deg\n", abs(c2), a2))
  cat(sprintf("abs(cos) with PC3=%.3f angle=%.1f deg\n", abs(c3), a3))
  cat(sprintf("Frac(|d|) captured by PC1-3 projection = %.3f\n", frac_1to3))

  rows[[length(rows)+1]] <- data.table(
    cohort = coh,
    n = n,
    prev = prev,
    abs_cos_pc1 = abs(c1),
    angle_pc1_deg = a1,
    abs_cos_pc2 = abs(c2),
    angle_pc2_deg = a2,
    abs_cos_pc3 = abs(c3),
    angle_pc3_deg = a3,
    frac_d_in_pc1to3 = frac_1to3,
    p_perm_align_pc1 = p_align_pc1
  )
}

tab <- rbindlist(rows)
# BH correction for the PC1 permutation p-values
tab[, p_perm_align_pc1_adj := p.adjust(p_perm_align_pc1, method="BH")]

dir.create("results/paper2", recursive=TRUE, showWarnings=FALSE)
out_csv <- "results/paper2/table_pd_effect_alignment_angles.csv"
fwrite(tab, out_csv)
cat("\nSaved:", out_csv, "\n")

# Plot angles (smaller angle = stronger alignment)
plt <- melt(
  tab[, .(cohort, angle_pc1_deg, angle_pc2_deg, angle_pc3_deg)],
  id.vars="cohort",
  variable.name="pc",
  value.name="angle_deg"
)
plt[, pc := fifelse(pc=="angle_pc1_deg","PC1",
             fifelse(pc=="angle_pc2_deg","PC2","PC3"))]

g <- ggplot(plt, aes(x=cohort, y=angle_deg, group=pc, color=pc)) +
  geom_point(size=3) +
  geom_line() +
  labs(x="Cohort", y="Angle between PD effect vector and PC axis (degrees)") +
  theme_bw(base_size=12)

fig_path <- "results/paper2/fig_pd_effect_alignment_angles.png"
ggsave(fig_path, g, width=7.2, height=4.2, dpi=220)
cat("Saved:", fig_path, "\n")

cat("C33 PD EFFECT ALIGNMENT ANGLES DONE\n")