library(data.table)

meta <- fread("data/processed/ipd_sex/ipd_metadata.csv")
X <- fread("data/processed/ipd_sex/ipd_genus_counts_union.csv")

setkey(meta, SampleID); setkey(X, SampleID)
dt <- X[meta]

# identify genera present in all 3 cohorts
genus_cols <- setdiff(names(dt), c("SampleID","PD","Sex","Cohort"))

present_by_cohort <- lapply(split(dt, dt$Cohort), function(d) {
  g <- genus_cols
  # present if any nonzero in cohort
  g[colSums(as.matrix(d[, ..g]) > 0) > 0]
})
shared <- Reduce(intersect, present_by_cohort)

cat("Shared genera across all cohorts:", length(shared), "\n")

# CLR on shared
mat <- as.matrix(dt[, ..shared])
storage.mode(mat) <- "numeric"
mat_pc <- mat + 1
log_mat <- log(mat_pc)
clr_mat <- log_mat - rowMeans(log_mat)

df <- data.table(
  PD = factor(dt$PD),
  Sex = factor(dt$Sex),
  Cohort = factor(dt$Cohort)
)
df[, Sex := relevel(Sex, ref="Female")]
df[, PD := relevel(PD, ref="0")]
df[, Cohort := relevel(Cohort, ref="finland")]

res <- rbindlist(lapply(seq_along(shared), function(j){
  g <- shared[j]
  y <- clr_mat[, j]
  fit <- lm(y ~ Cohort + PD*Sex, data=df)
  co <- summary(fit)$coefficients
  if (!("PD1:SexMale" %in% rownames(co))) return(NULL)
  data.table(
    Genus=g,
    beta_interaction=co["PD1:SexMale",1],
    p_interaction=co["PD1:SexMale",4],
    beta_PD_female=co["PD1",1],
    beta_PD_male=co["PD1",1] + co["PD1:SexMale",1]
  )
}))
res[, q_interaction := p.adjust(p_interaction, method="BH")]
fwrite(res, "results/sex_paper/interaction_shared_genera.csv")

cat("Top shared genera by p_interaction:\n")
print(res[order(p_interaction)][1:min(15,.N)])
