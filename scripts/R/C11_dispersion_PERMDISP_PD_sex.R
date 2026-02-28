cat("C11 STARTED\n")

library(data.table)
library(vegan)
library(ggplot2)

meta_file  <- "data/processed/ipd_sex/ipd_metadata.csv"
count_file <- "data/processed/ipd_sex/ipd_genus_counts_union.csv"

stopifnot(file.exists(meta_file), file.exists(count_file))

meta <- fread(meta_file)
X    <- fread(count_file)

setkey(meta, SampleID)
setkey(X, SampleID)

dt <- X[meta]
stopifnot(nrow(dt) > 0)

required <- c("SampleID","PD","Sex","Cohort")
stopifnot(all(required %in% names(dt)))

genus_cols <- setdiff(names(dt), required)

# ----- preprocessing -----
# relative abundance for Bray
mat_counts <- as.matrix(dt[, ..genus_cols])
storage.mode(mat_counts) <- "numeric"
rs <- rowSums(mat_counts)
if (any(rs == 0)) stop("Found samples with zero total counts.")
mat_rel <- mat_counts / rs

# CLR for Aitchison
clr_transform <- function(mat) {
  mat_pc <- mat + 1
  log_mat <- log(mat_pc)
  log_mat - rowMeans(log_mat)
}
mat_clr <- clr_transform(mat_counts)

# grouping factor
dt[, PD  := factor(PD)]
dt[, Sex := factor(Sex)]
# ensure canonical labels
dt[, PD  := ifelse(PD == "1", "PD", "Control")]
dt[, PD  := factor(PD, levels=c("Control","PD"))]
dt[, Sex := factor(Sex)]
# expected: Female/Male
if (!all(c("Female","Male") %in% levels(dt$Sex))) {
  warning("Sex levels are not exactly Female/Male. Levels are: ", paste(levels(dt$Sex), collapse=", "))
}

dt[, Group := interaction(PD, Sex, sep=":")]
dt[, Group := factor(Group, levels=c("Control:Female","Control:Male","PD:Female","PD:Male"))]

dir.create("results/sex_paper/dispersion", showWarnings=FALSE, recursive=TRUE)

# ----- helper: run permdisp -----
run_permdisp <- function(D, label) {

  # Overall dispersion by PD:Sex (pooled across cohorts)
  bd <- betadisper(D, dt$Group, type="centroid")

  # permutation test; stratify by cohort to avoid “cohort drives dispersion” artifact
  pt <- permutest(bd, permutations = 9999, strata = dt$Cohort)

  # distances to centroid table
  dd <- data.table(
    SampleID = dt$SampleID,
    Cohort   = dt$Cohort,
    PD       = dt$PD,
    Sex      = dt$Sex,
    Group    = dt$Group,
    dist     = bd$distances
  )

  fwrite(dd, sprintf("results/sex_paper/dispersion/dist_to_centroid_%s.csv", label))

  # summary stats
  summ <- dd[, .(
    n = .N,
    mean_dist = mean(dist),
    sd_dist = sd(dist)
  ), by=.(Cohort, Group)]
  fwrite(summ, sprintf("results/sex_paper/dispersion/summary_%s.csv", label))

  # pooled linear model on distances with cohort adjustment
  # interpretation: tests mean distance-to-centroid differences (dispersion) with cohort controlled
  fit <- lm(dist ~ Cohort + PD*Sex, data = dd)
  capture.output(summary(fit),
                 file = sprintf("results/sex_paper/dispersion/lm_dist_%s.txt", label))

  # plot: pooled (faceted by cohort)
  p <- ggplot(dd, aes(x=Group, y=dist)) +
    geom_boxplot(outlier.shape = 16) +
    facet_wrap(~Cohort, scales="free_y") +
    labs(
      title = sprintf("Dispersion (distance-to-centroid) by PD × Sex — %s", label),
      subtitle = sprintf("PERMDISP permutest (strata=Cohort): p=%.4g", pt$tab[1, "Pr(>F)"]),
      x = NULL,
      y = "Distance to centroid"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle=30, hjust=1))

  ggsave(sprintf("results/sex_paper/dispersion/box_dist_%s.png", label),
         p, width=11, height=6, dpi=160)

  # also export the permutest table
  out_tab <- as.data.table(pt$tab, keep.rownames = "term")
  fwrite(out_tab, sprintf("results/sex_paper/dispersion/permutest_%s.csv", label))

  invisible(list(bd=bd, pt=pt, dd=dd))
}

# ----- distances -----
D_bray <- vegdist(mat_rel, method="bray")
D_ait  <- dist(mat_clr)  # Euclidean on CLR = Aitchison distance

res_bray <- run_permdisp(D_bray, "bray")
res_ait  <- run_permdisp(D_ait,  "aitchison")

cat("C11 DONE\n")
