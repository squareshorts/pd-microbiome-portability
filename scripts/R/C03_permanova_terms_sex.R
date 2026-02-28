library(data.table)
library(vegan)

meta <- fread("data/processed/ipd_sex/ipd_metadata.csv")
X <- fread("data/processed/ipd_sex/ipd_genus_counts_union.csv")

setkey(meta, SampleID); setkey(X, SampleID)
dt <- X[meta]

genus_cols <- setdiff(names(dt), c("SampleID","PD","Sex","Cohort"))
mat <- as.matrix(dt[, ..genus_cols])
storage.mode(mat) <- "numeric"

# CLR
mat_pc <- mat + 1
log_mat <- log(mat_pc)
clr_mat <- log_mat - rowMeans(log_mat)

bray <- vegdist(mat, method="bray")
aitch <- dist(clr_mat)

md <- dt[, .(PD=factor(PD), Sex=factor(Sex), Cohort=factor(Cohort))]

set.seed(1)
perm_bray_terms  <- adonis2(bray  ~ Cohort + PD*Sex, data=md, permutations=999, by="margin")
perm_aitch_terms <- adonis2(aitch ~ Cohort + PD*Sex, data=md, permutations=999, by="margin")

dir.create("results/sex_paper", showWarnings=FALSE, recursive=TRUE)
capture.output(perm_bray_terms,  file="results/sex_paper/permanova_bray_margin.txt")
capture.output(perm_aitch_terms, file="results/sex_paper/permanova_aitchison_margin.txt")

cat("Wrote margin PERMANOVA results.\n")
