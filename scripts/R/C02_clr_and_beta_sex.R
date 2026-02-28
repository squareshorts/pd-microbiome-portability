library(data.table)
library(vegan)

meta <- fread("data/processed/ipd_sex/ipd_metadata.csv")
X <- fread("data/processed/ipd_sex/ipd_genus_counts_union.csv")

setkey(meta, SampleID)
setkey(X, SampleID)
dt <- X[meta]   # align

# Extract count matrix
genus_cols <- setdiff(names(dt), c("SampleID","PD","Sex","Cohort"))
mat <- as.matrix(dt[, ..genus_cols])
storage.mode(mat) <- "numeric"

# Pseudocount for CLR
mat_pc <- mat + 1

# CLR transform: clr(x) = log(x) - mean(log(x))
log_mat <- log(mat_pc)
clr_mat <- log_mat - rowMeans(log_mat)

# Distances
bray <- vegdist(mat, method="bray")          # on raw counts (you can also use relative)
aitch <- dist(clr_mat)                      # Euclidean on CLR = Aitchison distance

# PERMANOVA with blocking by cohort (strata)
md <- dt[, .(PD=factor(PD), Sex=factor(Sex), Cohort=factor(Cohort))]

set.seed(1)
perm_bray <- adonis2(bray ~ Cohort + PD*Sex, data=md, permutations=999)
perm_aitch <- adonis2(aitch ~ Cohort + PD*Sex, data=md, permutations=999)

dir.create("results/sex_paper", showWarnings=FALSE, recursive=TRUE)

capture.output(perm_bray, file="results/sex_paper/permanova_bray.txt")
capture.output(perm_aitch, file="results/sex_paper/permanova_aitchison.txt")

cat("Wrote PERMANOVA results to results/sex_paper/\n")
