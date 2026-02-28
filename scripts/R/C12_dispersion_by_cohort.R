cat("C12 STARTED\n")

library(data.table)
library(vegan)

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

# matrices
mat_counts <- as.matrix(dt[, ..genus_cols])
storage.mode(mat_counts) <- "numeric"

rs <- rowSums(mat_counts)
if (any(rs == 0)) stop("Found samples with zero total counts.")

mat_rel <- mat_counts / rs

clr_transform <- function(mat) {
  mat_pc <- mat + 1
  log_mat <- log(mat_pc)
  log_mat - rowMeans(log_mat)
}
mat_clr <- clr_transform(mat_counts)

# labels
dt[, PD  := ifelse(PD == 1, "PD", "Control")]
dt[, PD  := factor(PD, levels=c("Control","PD"))]
dt[, Sex := factor(Sex)]
dt[, Group := interaction(PD, Sex, sep=":")]
dt[, Group := factor(Group, levels=c("Control:Female","Control:Male","PD:Female","PD:Male"))]

dir.create("results/sex_paper/dispersion", showWarnings=FALSE, recursive=TRUE)

cohorts <- sort(unique(dt$Cohort))

run_one <- function(idx, cohort_name, label) {

  dsub <- dt[Cohort == cohort_name]
  if (nrow(dsub) < 40) return(NULL)

  # require all 4 groups present; otherwise PERMDISP on PD×Sex is not meaningful
  tab <- dsub[, .N, by=Group]
  if (nrow(tab) < 4) return(NULL)

  if (label == "bray") {
    D <- vegdist(mat_rel[idx, , drop=FALSE], method="bray")
  } else if (label == "aitchison") {
    D <- dist(mat_clr[idx, , drop=FALSE])
  } else stop("unknown label")

  bd <- betadisper(D, dsub$Group, type="centroid")
  pt <- permutest(bd, permutations=9999)

  dd <- data.table(
    SampleID = dsub$SampleID,
    Cohort   = dsub$Cohort,
    Group    = dsub$Group,
    dist     = bd$distances
  )

  # group means to interpret direction
  gmeans <- dd[, .(mean_dist = mean(dist), sd_dist = sd(dist), n=.N), by=Group]
  setorder(gmeans, -mean_dist)

  data.table(
    Cohort = cohort_name,
    metric = label,
    F = unname(pt$tab[1,"F"]),
    p = unname(pt$tab[1,"Pr(>F)"]),
    top_disp_group = as.character(gmeans$Group[1]),
    top_disp_mean  = gmeans$mean_dist[1],
    bottom_disp_group = as.character(gmeans$Group[.N]),
    bottom_disp_mean  = gmeans$mean_dist[.N]
  )
}

out <- rbindlist(lapply(cohorts, function(coh) {
  idx <- which(dt$Cohort == coh)
  rbindlist(list(
    run_one(idx, coh, "bray"),
    run_one(idx, coh, "aitchison")
  ), fill=TRUE)
}), fill=TRUE)

fwrite(out, "results/sex_paper/dispersion/permdisp_by_cohort.csv")

cat("Wrote results/sex_paper/dispersion/permdisp_by_cohort.csv\n")
print(out)

cat("C12 DONE\n")
